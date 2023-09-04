import json
import sys
import shutil
import os
import math
import numpy
import scipy.integrate
import scipy.interpolate
import scipy.fft
import matplotlib.pyplot
try:
    import ROOT
except ImportError:
    raise ImportError("PyROOT import failed - you can still run the code, but minuit fitter won't work")
from models.field_models import FieldSum
from models.field_models import FlatGaussField
from models.field_models import LinearInterpolation
from models.field_models import SineField
from models.field_models import UniformField
from models.field_models import CurrentSheet


"""
Two classes to evolve the beta function
BetaFinder - uses a second order differential equation
BetaFinder2 - integrates a transfer matrix

Some plotting routines as well

Units are T, GeV/c, m
"""

class BetaFinder(object):
    """
    Calculate transfer matrix by evolving infinitesimal (decoupled) TM

    Find the periodic solution using the usual alpha/beta relationship (if phase advance real)
    """
    def __init__(self, field, momentum, use_analytic=True):
        """
        Initialise
        - field is an object of type field_models.Field
        - momentum is a float (momentum in GeV/c)
        """
        self.momentum = momentum
        self.verbose = 0
        self.q = 1
        self.use_analytic = use_analytic
        self.reset_field(field)

    def field(self, z):
        """
        Return the field at position z
        """
        return self.field_model.get_field(z)

    def reset_field(self, field):
        """
        Reset the field, potentially with a different cell period etc
        """
        self.field_model = field
        self.period = self.field_model.get_period()
        self.hmax = self.period*1e-4

    def matrix_derivatives(self, y, z):
        """
        Returns dM/dz at position z given M
        - y is a list of floats like M[0,0], M[0,1], M[1,0], M[1,1], larmor_angle
        - z is a float z position [m]
        Returns a list of floats derivatives with respect to z of 
        M[0,0], M[0,1], M[1,0], M[1,1], larmor_angle 
        (for input into scipy.odeint)
        """
        pz = self.momentum
        bz = self.field(z)
        q = 1
        b0 = 0.3*q*bz
        matrix = numpy.array([[y[0], y[1]], [y[2], y[3]]])
        dmatrix = numpy.array([[0.0, 1/pz], [-b0 **2/4.0/pz, 0.0]])
        mderiv = numpy.dot(dmatrix, matrix)
        delta_matrix = (mderiv[0,0], mderiv[0,1], mderiv[1,0], mderiv[1,1], -b0/2/pz)
        return delta_matrix

    def get_beta_periodic(self):
        if self.use_analytic:
            return self.get_beta_periodic_analytic()
        else:
            return self.get_beta_periodic_minuit(0.5)

    def get_beta_periodic_analytic(self):
        """
        Returns a tuple of (beta, alpha, phase advance) where:
        - beta is in [m]
        - phase_advance is in [rad]
        beta should be periodic, so beta at the start of one period is the same 
        as beta at the end of one period. If there is no solution, returns (0,0,0)
        """
        zout = [0.0, self.period]
        matrix = scipy.integrate.odeint(self.matrix_derivatives,
                                        [1, 0.0, 0.0, 1.0, 0.0], # m_00 m_01 m_10 m_11 larmor
                                        zout, hmax=self.hmax,
                                        mxstep=int(zout[-1]/self.hmax*2))[-1]
        larmor_angle = matrix[4]
        m2d = numpy.array([[matrix[0], matrix[1]], [matrix[2], matrix[3]]])
        cosmu = (m2d[0,0]+m2d[1,1])/2
        #print("Transfer matrix with cosmu", cosmu, "det", numpy.linalg.det(m2d), "\n", m2d)
        if abs(cosmu) > 1:
            return (0.0, 0.0, 0.0)
        n2d = m2d - numpy.array([[cosmu, 0.0],[0.0,cosmu]])
        #print("N2D\n", n2d)
        # absolute value of sin(mu) from square root; we use beta is pos definite to force the sign of sinmu
        sinmu = numpy.linalg.det(n2d)**0.5*numpy.sign(n2d[0,1]) 
        n2d /= sinmu
        #print("N2D over sinmu with sinmu=", sinmu, "\n", n2d)
        v2d = numpy.array([[n2d[0,1], -n2d[0,0]], [-n2d[0,0], -n2d[1, 0]]])
        # v2d[0,0] = beta/p
        beta, alpha, phase_advance = v2d[0, 0]*self.momentum, -v2d[0,1], -math.atan2(cosmu, sinmu)
        #print("beta alpha phi", beta, alpha, phase_advance)
        return beta, alpha, phase_advance

    def beta_derivatives(self, y, z):
        """
        Returns a tuple of optical parameters like (dbeta/dz, d2beta/dz2, dphi/dz)
        - y is a list like [beta, dbeta/dz, phi] with units [m], [], [radians]
        - z is z-position in [m]
        phi is the phase advance. This is used by odeint to evolve beta/etc
        """

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

    def evolve(self, beta0, dbetadz0, zin, zout):
        """
        Calculates beta, dbeta/dz, phi at some position zout
        - beta0: initial optical beta function
        - dbetadz0: initial derivative of optical beta function
        - zout: final z position
        Returns a tuple of beta, dbeta.dz, phi with units [m], [], [rad]
        """
        zpoints = [zin, zout]
        output = scipy.integrate.odeint(self.beta_derivatives,
                                        [beta0, dbetadz0, 0.],
                                        zpoints, hmax=self.hmax,
                                        mxstep=int(zpoints[-1]/self.hmax*2))
        return output[-1]

    def is_not_periodic(self, beta0, beta1):
        """
        Check if beta0 and beta1 are the same within tolerances
        """
        test_out = abs(beta0-beta1)*2/(beta0+beta1) > 0.05 or abs(beta0-beta1) > 1
        if test_out:
            print("Not periodic", self.momentum, beta0-beta1, abs(beta0-beta1)*2/(beta0+beta1), test_out)
        return test_out

    def minuit_fitter(self, seed_beta0, err_beta0, max_beta0, n_iterations, tolerance):
        """
        Drive minuit to find a periodic beta function
        - called by get_beta_periodic_minuit
        """
        self.minuit = ROOT.TMinuit()
        self.minuit.SetPrintLevel(-1)
        self.minuit.DefineParameter(0, "beta0", seed_beta0, err_beta0, 0.0, max_beta0)
        self.minuit.SetFCN(self.score_function)
        self.minuit.Command("SIMPLEX "+str(n_iterations)+" "+str(tolerance))
        beta0 = self.score_function(0, 0, [0], 0, 0)
        return beta0

    def score_function(self, nvar, parameters, score, jacobian, err):
        """
        Minuit score function used by the minuit root finding routines
        """
        beta0 = ROOT.Double()
        err = ROOT.Double()
        self.minuit.GetParameter(0, beta0, err)
        beta1, dbetadz, phi = self.evolve(float(beta0), 0., 0.0, self.period)
        score[0] = abs(beta1-beta0)+abs(dbetadz*100)**2
        #print (beta0, beta1, dbetadz, score[0])
        return beta0, beta1, dbetadz, phi

    def get_beta_periodic_minuit(self, seed_beta0):
        """
        A slow alternative to the transfer matrix approach to get_beta_periodic.
        Use minuit to try to find a periodic beta function numerically.
        - seed_beta0: guess at initial periodic beta function
        Returns a tuple of beta(z=0), beta(z=period/2), phase_advance. If no 
        solution is found, returns (0.0, 0.0, 0.0)
        """
        beta0, beta1, dbetadz, phi = self.minuit_fitter(seed_beta0, seed_beta0/10., 100.0, 500, 1e-5)
        #print(format(self.momentum, "8.4g"), format(beta0, "12.6g"), format(beta1, "12.6g"), dbetadz, abs(beta0-beta1)*2.0/(beta0+beta1))
        if self.is_not_periodic(beta0, beta1):
            return 0.0, 0.0, 0.0
        return beta1, dbetadz, phi

    def propagate_beta(self, beta0, dbetadz0, n_points=101):
        """
        Propagate beta and return a tuple of lists with beta as a function of z
        - beta0: initial beta
        - dbetadz0: initial derivative of beta w.r.t. z
        - n_points: number of z points
        Returns a tuple like 
        - z_list: list of z positions
        - beta_list: list of beta 
        - dbetadz_list: list of first derivatives of beta
        - phi_list: list of phase advance
        """
        z_list = [self.period*i/float(n_points-1) for i in range(n_points)]
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


def clear_dir(a_dir):
    try:
        shutil.rmtree(a_dir)
    except OSError:
        pass
    os.makedirs(a_dir)

def round_delta_sf(a_float_1, a_float_2, n_sf):
    """Round to some number of significant figures in the difference"""
    delta = abs(a_float_1-a_float_2)
    dlog = math.log10(delta)
    delta_rounded = 10**int(dlog)
    a_float_1 = round(a_float_1/delta_rounded, n_sf)*delta_rounded
    a_float_2 = round(a_float_2/delta_rounded, n_sf)*delta_rounded
    return [a_float_1, a_float_2]


def make_kinetic_energy_axis(p_axes):
    mass = 0.105865
    ke_calc = lambda p: (p**2+mass**2)**0.5-mass
    p_calc = lambda ke: ((ke+mass)**2-mass**2)**0.5
    ke_axes = p_axes.twiny()
    p_lim = p_axes.get_xlim()
    ke_lim = (ke_calc(p_lim[0]), ke_calc(p_lim[1]))
    ke_rounded_lim = round_delta_sf(ke_lim[0], ke_lim[1], 1)
    delta = (ke_rounded_lim[1]-ke_rounded_lim[0])/5
    ke_labels = [ke_rounded_lim[0]]
    while ke_labels[-1] < ke_lim[1]:
        ke_labels.append(ke_labels[-1]+delta)
        #print(ke_labels[-1], delta, ke_lim[1])
    tick_position = [p_calc(ke) for ke in ke_labels]
    ke_labels = [format(ke, "2.2g") for ke in ke_labels]
    ke_axes.set_xticks(tick_position)
    ke_axes.set_xticklabels(ke_labels)
    ke_axes.set_xlim(p_lim)
    #ylim = ke_axes.get_ylim()
    #ke_axes.plot([1e-5, 1e-5], ylim)
    #ke_axes.set_ylim(ylim)
    print(ke_axes.get_xlim())
    print("Setting ke axis with ticks at momentum", tick_position, "ke", ke_labels)
    return ke_axes



fignum = 1
def do_plots(field, pz_plot_list, pz_list, plot_dir, fig1, fig2, fig3):
    """
    Plot the beta function for a given field model
    - field: the field model. Should be of type field_model.Field
    - pz0: a reference momentum for plotting vs z
    - pz_list: a set of pz values for finding the periodic solution and plotting
               vs pz
    """
    global fignum
    #matplotlib.rcParams['text.usetex'] = True
    bsquared = None
    period = field.get_period()
    name = "optics_L_"+str(period)+"_"+field.get_name()
    out_data = {"pz_items":[], "field_name":field.get_name(), "field_bi":field.bz_list, "field_b0":field.bz0, "period":field.get_period()}
    beta_list = []
    antibeta_list = []
    phi_list = []
    int_bz2 = field.get_bz2_int()
    inv_kappa_list = [0.15*pz/int_bz2**0.5 for pz in pz_list]
    beta_finder = BetaFinder(field, 1., use_analytic=True)
    if bsquared != None:
        beta_finder.field_model.normalise_bz_squared(bsquared)
    z_list = [i*period/100 for i in range(101)]
    bz_list = [beta_finder.field_model.get_field(z) for z in z_list]
    bz2_list = [bz**2 for bz in bz_list]
    if fig1 == None:
        figure = matplotlib.pyplot.figure(fignum, figsize=(12, 10))
    else:
        figure = fig1
    human = field.human_readable()
    if human:
        human = "$"+human+"$\n"
    figure.suptitle(human+\
           "$\\int B^2(z) dz =$ {0:.4g} T$^2$ m".format(field.get_bz2_int()))
    fignum += 1
    if len(figure.get_axes()) == 0:
        axes = [figure.add_subplot(2, 2, 1), figure.add_subplot(2, 2, 2),
                figure.add_subplot(2, 2, 3), figure.add_subplot(2, 2, 4)]
        axes.append(axes[2].twinx())
    else:
        axes = figure.get_axes()

    zero_list = [0. for i in z_list]
    #axes[0].plot(z_list, zero_list, '--g')
    axes[0].plot(z_list, bz_list)
    axes[0].set_xlabel("z [m]")
    axes[0].set_ylabel("B$_{z}$ [T]")
    axes[1].plot(z_list, bz2_list)
    axes[1].set_xlabel("z [m]")
    axes[1].set_ylabel("B$_{z}^2$ [T$^2$]")

    print("     pz    beta_0     phi      n_iterations")
    for pz in pz_list:
        beta_finder.momentum = pz
        beta, alpha, phi = beta_finder.get_beta_periodic()
        print ("    ", pz, beta, phi)
        beta_list.append(beta)
        antibeta_list.append(0.0)
        phi_list.append(phi)
    out_data["pz_list"] = pz_list
    out_data["beta_list"] = beta_list
    out_data["antibeta_list"] = antibeta_list
    axes[2].plot(pz_list, beta_list, label="$\\beta(L)$")
    axes[1].plot(pz_list, antibeta_list, 'g--', label="$\\beta(L/2)$")
    axes[2].set_xlabel("p$_{z}$ [GeV/c]")
    axes[2].set_ylabel("$\\beta$ [m]")
    axes[2].set_ylim([0.0, 0.5])

    for y in [math.pi/2.0, 0.0, -math.pi/2.0]:
        pass #axes[4].plot([pz_list[0], pz_list[-1]], [y, y], linestyle='dotted', color='pink')
    axes[4].set_ylim(-math.pi, math.pi)
    axes[4].tick_params(axis='y', labelcolor='r')
    axes[4].tick_params(axis='y', labelcolor='r')
    axes[4].set_ylim(-math.pi, math.pi)
    axes.append(make_kinetic_energy_axis(axes[2]))
    this_pz_list, this_phi_list = [], []
    phi_old = None
    plotting = False
    for i, pz in enumerate(pz_list):
        phi = phi_list[i]
        pz = pz_list[i]
        if phi == 0 and i != len(pz_list)-1:
            plotting = True
            continue
        #if i > 0 and abs(phi-phi_list[i-1]) > math.pi:
        #    plotting = True
        if i == len(pz_list)-1:
            plotting = True
        if i > 0 and abs(phi_list[i] - phi_list[i-1]) > math.pi:
            plotting = True
        if plotting:
            print("Plotting", this_pz_list)
            #axes[4].plot(this_pz_list, this_phi_list, 'r-.')
            this_pz_list, this_phi_list = [], []
            plotting = False
        this_pz_list.append(pz)
        this_phi_list.append(phi)


    for pz in pz_plot_list:
        print ("PZ PLOT LIST", pz, pz_plot_list)
        beta_finder.momentum = pz
        beta, alpha, phi = beta_finder.get_beta_periodic()
        if beta < 1e-9:
            continue
        print("PLOTTING", pz, beta)
        beta_finder.verbose = 0
        z_list, beta_list, dbetadz_list, phi_list = \
                                        beta_finder.propagate_beta(beta, 0.)
        print(z_list)
        print(beta_list)
        axes[3].plot(z_list, beta_list, label="p$_z$ "+format(pz, "6.4g")+" GeV/c")
        out_data["pz_items"].append({"z":z_list, "beta_list":beta_list, "pz":pz})
    axes[3].set_xlabel("z [m]")
    axes[3].set_ylabel("$\\beta$ [m]")
    axes[3].set_ylim([0.0, 1.6])
    axes[3].legend()

    emittance = 0.0003 # metres
    mass = 0.105658 # GeV/c^2
    #for emittance in [0.0003]:
        #sigma_x_list = [(beta*emittance*mass/beta_finder.momentum)**0.5 for beta in beta_list]
        #axes[3].plot(z_list, sigma_x_list)
    #axes[3].grid(1)
    #axes[3].set_xlabel("z [m]")
    #axes[3].set_ylabel("$\\sigma_x$ [m]")
    #name = "optics_L_"+str(period)+"_"+field.get_name()+"_pz_"+str(pz)+".png"
    fout = open(os.path.join(plot_dir, name+".json"), "w")
    fout.write(json.dumps(out_data))
    figure.savefig(os.path.join(plot_dir, name+".png"))
    return figure, None, None

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

def main():
    global fignum
    plot_dir = "optics-scan_v22"
    pz_plot_list = [0.19, 0.2, 0.21]
    pz_scan_list = [pz_i/1000. for pz_i in range(130, 251, 1)]
    n_points = 2
    norm = 160
    clear_dir(plot_dir)
    b1 = 7.0
    b2 = 1.0
    b3 = 0.0
    for b1, b2, b3 in [
            (7.0, 1.0, 0.0),
            (7.5, 1.0, 0.0),
            (7.0, 1.25, 0.0),
            (7.0, 1.0, 0.5)
        ]:
    #for bi in range(0, 21, 1):
    #    b3 = (bi-10)*0.05
        solenoid_field = SineField(0.0, b3, b1, 0.0, b2, 0.0, 2.0)
        fig1, fig2, fig3 = do_plots(solenoid_field, pz_plot_list, pz_scan_list, plot_dir, None, None, None)
    matplotlib.pyplot.show(block=False)

if __name__ == "__main__":
    main()
    input("Press <CR> to end")
