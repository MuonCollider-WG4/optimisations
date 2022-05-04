import json
import sys
import math
import numpy
import scipy.integrate
import scipy.interpolate
import scipy.fft
import matplotlib.pyplot
from platypus import Problem
import platypus
import ROOT
from field_models import FieldSum
from field_models import FlatGaussField
from field_models import LinearInterpolation
from field_models import SineField
from field_models import UniformField
from field_models import CurrentSheet


class BetaFinder2(object):
    """
    Calculate transfer matrix by evolving infinitesimal (decoupled) TM

    Find the periodic solution using the usual alpha/beta relationship (if phase advance real)
    """
    def __init__(self, field, momentum):
        self.field_model = field
        self.period = self.field_model.get_period()
        self.momentum = momentum
        self.hmax = self.period*1e-4
        self.verbose = 0
        self.q = 1

    def field(self, z):
        return self.field_model.get_field(z)

    def matrix_derivatives(self, y, z):
        pz = self.momentum
        bz = self.field(z)
        q = 1
        b0 = 3.0*q*bz
        matrix = numpy.array([[y[0], y[1]], [y[2], y[3]]])
        dmatrix = numpy.array([[0.0, 1/pz], [-b0 **2/4.0/pz, 0.0]])
        mderiv = numpy.dot(dmatrix, matrix)
        delta_matrix = (mderiv[0,0], mderiv[0,1], mderiv[1,0], mderiv[1,1], -b0/2/pz)
        return delta_matrix

    def get_beta_periodic(self):
        zout = [0.0, self.period]
        matrix = scipy.integrate.odeint(self.matrix_derivatives,
                                        [1, 0.0, 0.0, 1.0, 0.0], # m_00 m_01 m_10 m_11 larmor
                                        zout, hmax=self.hmax,
                                        mxstep=int(zout[-1]/self.hmax*2))[-1]
        larmor_angle = matrix[4]
        m2d = numpy.array([[matrix[0], matrix[1]], [matrix[2], matrix[3]]])
        cosmu = (m2d[0,0]+m2d[1,1])/2
        print("Transfer matrix with cosmu", cosmu, "\n", m2d)
        if abs(cosmu) > 1:
            return (0, 0)
        n2d = m2d - numpy.array([[cosmu, 0.0],[0.0,cosmu]])
        sinmu = numpy.linalg.det(n2d)**0.5
        n2d /= sinmu
        v2d = numpy.array([[n2d[0,1], -n2d[0,0]], [-n2d[0,0], -n2d[1, 0]]])
        beta, alpha, phase_advance = v2d[0, 1], -v2d[0,0]/2.0, math.atan2(cosmu, sinmu)
        print("beta alpha phi", beta, alpha, phase_advance)
        return beta, phase_advance

class BetaFinder(object):
    """
    Calculate beta function by evolving equation for beta''

    Find the periodic solution by firing minuit
    """
    def __init__(self, field, momentum):
        self.field_model = field
        self.period = self.field_model.get_period()
        self.momentum = momentum
        self.hmax = self.period*1e-4
        self.minuit = None
        self.verbose = 0

    def field(self, z):
        return self.field_model.get_field(z)

    def beta_derivatives(self, y, z):
        pz = self.momentum
        bz = self.field(z)
        kappa = 0.15*bz/pz
        Ltwiddle = 0
        beta = y[0]
        dbetadz = y[1]
        d2betadz2 = +(dbetadz)**2 \
                    -4*beta**2*kappa**2 \
                    +4*(1+Ltwiddle**2)
        d2betadz2 = d2betadz2/(2*beta)
        dphidz = 1/beta
        dydz = (dbetadz, d2betadz2, dphidz)
        if self.verbose > 0:
            print("beta derivatives z:", z, "bz:", bz, "k:", kappa, "y:", y, "dy/dz:", dydz)
        return dydz

    def evolve(self, beta0, dbetadz0, zout):
        zout = [0.0, zout]
        output = scipy.integrate.odeint(self.beta_derivatives,
                                        [beta0, dbetadz0, 0.],
                                        zout, hmax=self.hmax,
                                        mxstep=int(zout[-1]/self.hmax*2))
        return output[-1]

    def is_not_periodic(self, beta0, beta1):
        test_out = abs(beta0-beta1)*2/(beta0+beta1) > 0.05 or abs(beta0-beta1) > 1
        if test_out:
            print("Not periodic", self.momentum, beta0-beta1, abs(beta0-beta1)*2/(beta0+beta1), test_out)
        return test_out

    def minuit_fitter(self, seed_beta0, err_beta0, max_beta0, n_iterations, tolerance):
        self.minuit = ROOT.TMinuit()
        self.minuit.SetPrintLevel(-1)
        self.minuit.DefineParameter(0, "beta0", seed_beta0, err_beta0, 0.0, max_beta0)
        self.minuit.SetFCN(self.score_function)
        self.minuit.Command("SIMPLEX "+str(n_iterations)+" "+str(tolerance))
        beta0 = self.score_function(0, 0, [0], 0, 0)
        return beta0

    def score_function(self, nvar, parameters, score, jacobian, err):
        beta0 = ROOT.Double()
        err = ROOT.Double()
        self.minuit.GetParameter(0, beta0, err)
        beta1, dbetadz, phi = self.evolve(float(beta0), 0., self.period)
        score[0] = abs(beta1-beta0)+abs(dbetadz*100)**2
        #print (beta0, beta1, dbetadz, score[0])
        return beta0, beta1, dbetadz, phi

    def get_beta_periodic(self, seed_beta0):
        beta0, beta1, dbetadz, phi = self.minuit_fitter(seed_beta0, seed_beta0/10., 100.0, 500, 1e-5)
        print(format(self.momentum, "8.4g"), format(beta0, "12.6g"), format(beta1, "12.6g"), dbetadz, abs(beta0-beta1)*2.0/(beta0+beta1))
        if self.is_not_periodic(beta0, beta1):
            return 0.0, 0.0, 0.0
        return beta1, dbetadz, phi

        beta0, beta1 = seed_beta0, seed_beta0*2
        iterations = 0
        while self.is_not_periodic(beta0, beta1) and iterations < 500:
            beta1, dbetadz, phi = self.evolve(beta0, 0., self.period)
            beta0 = abs(beta1+beta0)/2
            iterations += 1
            if beta0 > 100.:
                return 0., 0., 0.
        print(format(self.momentum, "8.4g"), format(beta0, "12.6g"), format(beta1, "12.6g"), dbetadz, abs(beta0-beta1)*2.0/(beta0+beta1), iterations)
        antibeta, antidbetadz, antiphi = self.evolve(beta0, 0., self.period/2)
        if self.is_not_periodic(beta0, beta1):
            return 0., 0., 0.
        return beta0, antibeta, phi

    def propagate_beta(self, beta0, dbetadz0):
        z_list = [self.period*i/100. for i in range(101)]
        output_list = []
        out, infodict = scipy.integrate.odeint(self.beta_derivatives,
                                     [beta0, dbetadz0, 0.],
                                     z_list, hmax=self.hmax, full_output=True,
                                     mxstep=int(z_list[-1]/self.hmax*2))
        if self.verbose > 0:
            for key, value in infodict.items():
                print(key, value)
        output_list = out
        #print(z, self.momentum, output_list[-1])
        beta_list = [output[0] for output in output_list]
        dbetadz_list = [output[1] for output in output_list]
        phi_list = [output[2] for output in output_list]
        return z_list, beta_list, dbetadz_list, phi_list

class GaussMinimiser(object):
    def __init__(self, period, width_seed, peak, tail, momentum):
        self.gauss_field = FlatGaussField(period, width_seed, peak, tail)
        self.width_seed = width_seed
        self.momentum = momentum
        self.initial_beta = 1./(0.15*peak/self.momentum)
        self.iteration = 0
        self.algorithm = None
        self.score_list = []

    def setup_optimisation_cmaes(self):
        n_parameters = 1
        n_variables = 1
        problem = platypus.Problem(n_parameters, n_variables)
        problem.types[:] = platypus.Real(0., self.width_seed*2)
        problem.directions[:] = Problem.MINIMIZE
        problem.function = self.run_one
        algorithm = platypus.CMAES(problem)
        algorithm.initial_search_point = [self.width_seed] 
        algorithm.sigma = self.width_seed*1e-1
        self.iteration = 0
        self.score_list = []
        self.algorithm = algorithm

    def setup_optimisation_root(self):
        self.minuit = ROOT.TMinuit(1)
        self.minuit.DefineParameter(0,
                                    "fringe_field",
                                    self.width_seed,
                                    abs(self.width_seed/100.),
                                    self.gauss_field.period/20.,
                                    self.gauss_field.period/6.)
        self.minuit.SetFCN(self.minuit_run_one)

    def minuit_run_one(self, nvar, parameters, score, jacobian, err):
        var = ROOT.Double()
        err = ROOT.Double()
        self.minuit.GetParameter(0, var, err)
        output = self.run_one([var])
        score[0] = output[0]

    def run_one(self, arguments):
        self.gauss_field.width = arguments[0]
        beta_finder = BetaFinder(self.gauss_field, self.momentum)
        alpha_half = beta_finder.evolve(self.initial_beta, 0., self.gauss_field.period/2.)[1]
        #alpha_full = beta_finder.evolve(self.initial_beta, 0., self.gauss_field.period)[1]
        print(str(self.iteration).ljust(4), format(arguments[0], "12.8g"), format(alpha_half, "8.4g"))
        self.iteration += 1
        self.score_list.append([alpha_half**2, arguments])
        return [alpha_half**2] #, alpha_full**2]

    def optimise(self, n_trials):
        self.setup_optimisation_root()
        try:
            self.minuit.Command("SIMPLEX"+" "+str(n_trials)+" 1.e-12")
        except Exception:
            sys.excepthook(*sys.exc_info())
            print("Minuit failed")
        self.score_list = sorted(self.score_list)
        self.run_one(self.score_list[0][1])
        return self.gauss_field

    def do_plot(self):
        global fignum
        beta_finder = BetaFinder(self.gauss_field, self.momentum)
        z_list, beta_list, dbetadz_list, phi_list = \
                            beta_finder.propagate_beta(self.initial_beta, 0.)
        figure = matplotlib.pyplot.figure(fignum, figsize=(12, 10))
        axes = figure.add_subplot(1, 2, 1)
        axes.plot(z_list, beta_list, 'b')
        axes.set_xlabel("z [m]")
        axes.set_ylabel("$\\beta$ [m]")

        bz_list = [self.gauss_field.get_field(z) for z in z_list]
        axes = figure.add_subplot(1, 2, 2)
        axes.plot(z_list, bz_list, 'b')
        axes.set_xlabel("z [m]")
        axes.set_ylabel("B$_{z}$ [T]")

        name = "minimiser_L_"+self.gauss_field.get_name()+"_pz_"+str(self.momentum)+".png"

        figure.savefig(name+".png")
        fignum += 1

#    period = 0.5
#    bz0 = 0.0
#    bz1 = 4.0
#    bz2 = 4.0
#    bz3 = 0.0
#    pz0 = 0.07

# period = 0.75
# bz0 = 0.0
# bz1 = 16.0
# bz2 = 0.0
# bz3 = 16.0
# bsquared = 20.0
# pz0 = 0.07

fignum = 1
def do_plots(field, pz0, pz_list):
    global fignum
    #matplotlib.rcParams['text.usetex'] = True
    bsquared = None
    period = field.get_period()
    beta_list = []
    antibeta_list = []
    phi_list = []
    beta_finder = BetaFinder(field, 1.)
    beta_finder2 = BetaFinder2(field, 1.)
    if bsquared != None:
        beta_finder.field_model.normalise_bz_squared(bsquared)
    z_list = [i*period/1000 for i in range(1001)]
    bz_list = [beta_finder.field_model.get_field(z) for z in z_list]
    bz2_list = [bz**2 for bz in bz_list]
    figure = matplotlib.pyplot.figure(fignum, figsize=(12, 10))
    figure.suptitle("$"+field.human_readable()+"$")
    fignum += 1
    if len(figure.get_axes()) == 0:
        axes = [figure.add_subplot(2, 2, 1), figure.add_subplot(2, 2, 2),
                figure.add_subplot(2, 2, 3), figure.add_subplot(2, 2, 4)]
        axes.append(axes[1].twinx())
    else:
        axes = figure.get_axes()
    zero_list = [0. for i in z_list]
    #axes[0].plot(z_list, zero_list, '--g')
    axes[0].plot(z_list, bz_list)
    axes[0].set_xlabel("z [m]")
    axes[0].set_ylabel("B$_{z}$ [T]")
    axes[1].plot(z_list, bz2_list)
    axes[1].set_xlabel("z [m]")
    axes[1].set_ylabel("B$_{z}^2$ [T]")

    beta = 1.0
    print("     pz    beta_0     beta_1      dbetadz       delta      n_iterations")
    for pz in pz_list:
        beta_finder.momentum = pz
        beta_finder2.momentum = pz
        if beta < 1e-9:
            beta = 0.1
        beta, antibeta, phi = beta_finder.get_beta_periodic(beta)
        beta2, phi2 = beta_finder2.get_beta_periodic()
        print ("    ", pz, beta, beta2)
        beta_list.append(beta)
        antibeta_list.append(antibeta)
        phi_list.append(phi)
    axes[2].plot(pz_list, beta_list, label="$\\beta(L)$")
    #axes[1].plot(pz_list, antibeta_list, 'g--', label="$\\beta(L/2)$")
    axes[2].set_xlabel("p$_{z}$ [GeV/c]")
    axes[2].set_ylabel("$\\beta$ [m]")
    axes[2].set_ylim([0.0, 0.5])
    """
    axes[4].set_ylabel("$\\phi$ [rad]", color='r')
    axes[4].tick_params(axis='y', labelcolor='r')
    for i in range(int(min(phi_list)/math.pi)+1, int(max(phi_list)/math.pi)+1):
        pi_list = [i*math.pi, i*math.pi]
        pi_pz_list = [min(pz_list), max(pz_list)]
        label = None
        if i == 0:
            label = "$\\phi$"
        #axes[4].plot(pi_pz_list, pi_list, ':', color='pink', label=label)
    #axes[4].plot(pz_list, phi_list, 'r-.')
    """

    for pz in [pz0]:
        beta_finder.momentum = pz
        beta, antibeta, phi = beta_finder.get_beta_periodic(pz)
        beta_finder.verbose = 0
        z_list, beta_list, dbetadz_list, phi_list = \
                                        beta_finder.propagate_beta(beta, 0.)
        beta_finder.verbose = 0
        axes[3].plot(z_list, beta_list, label="p$_z$ "+format(pz, "6.4g")+" GeV/c")
        axes[3].set_xlabel("z [m]")
        axes[3].set_ylabel("$\\beta$ [m]")
        axes[3].set_ylim([0.0, 0.5])
        axes[3].legend()

    emittance = 0.0003 # metres
    mass = 0.105658 # GeV/c^2
    for emittance in [0.0003]:
        sigma_x_list = [(beta*emittance*mass/beta_finder.momentum)**0.5 for beta in beta_list]
        #axes[3].plot(z_list, sigma_x_list)
    #axes[3].grid(1)
    #axes[3].set_xlabel("z [m]")
    #axes[3].set_ylabel("$\\sigma_x$ [m]")

    name = "optics_L_"+str(period)+field.get_name()+"_pz_"+str(pz)+".png"
    figure.savefig(name)

def make_sine_field(sol_bz2, bz1, bz2, bz3):
    factor = 1
    period = 1.0
    bz0 = 0.0
    #bz1 = 7.206
    #bz2 = 1.0
    #bz3 = 0.0
    f = SineField(factor*bz0, factor*bz1, factor*bz2, factor*bz3, period)#
    if sol_bz2 > 0:
        f_bz2 = f.get_bz2_int()
        factor = (sol_bz2/f_bz2)**0.5
        f = SineField(factor*bz0, factor*bz1, factor*bz2, factor*bz3, period)#
    return f


def make_solenoid_field():
    #b0, zcentre, length, radius, period, nrepeats
    period = 1.0
    bz0 = 7.206
    field_list = [(1.0, 0.1, 0.1, 0.4)]
    solenoid_list = [CurrentSheet(b0, z0, l0, r0, period, 4) for b0, z0, l0, r0 in field_list]
    solenoid_list += [CurrentSheet(-b0, -z0, l0, r0, period, 4) for b0, z0, l0, r0 in field_list]
    field_sum = FieldSum(solenoid_list)
    bz_norm = max([field_sum.get_field(i*period/100.0) for i in range(101)])
    for field in field_sum.field_list:
        field.b0 = field.b0*bz0/bz_norm # watch the sign on the field
    bz_norm = max([field_sum.get_field(i*period/100.0) for i in range(101)])
    return field_sum

def fft_field(truncation):
    source_field = make_solenoid_field()
    period = source_field.get_period()
    field_values = [source_field.get_field(i*period/100.0) for i in range(100)]
    my_fft = scipy.fft.fft(field_values)
    trunc_fft = [x if i < truncation else 0.0 for i, x in enumerate(my_fft)]
    inverse = numpy.real(scipy.fft.ifft(trunc_fft))
    interpolation = LinearInterpolation(inverse, period)
    interpolation.name = "fft_truncated_"+str(interpolation)
    print(trunc_fft, inverse)
    return interpolation

def main():
    global fignum
    pz_list = [pz_i/1000. for pz_i in range(150, 251, 10)]
    sol_field = make_solenoid_field()
    n_points = 2
    for i in range(1, n_points+1):
        bz1 = 7.206 # 10-abs(i)
        bz2 = i
        bz3 = 0
        #sol_field.get_bz2_int()
        sine_field = make_sine_field(-1, bz1, bz2, bz3)
        sine_field.normalise_bz_squared(25)
        print("Bz2", sine_field.get_bz2_int())
        do_plots(sine_field, 0.2, pz_list)
    #fft_field(10000)
    #do_plots(fft_field(10000), 0.2, pz_list)
    fignum -= 1
    matplotlib.pyplot.show(block=False)

if __name__ == "__main__":
    main()
    input("Press <CR> to end")
