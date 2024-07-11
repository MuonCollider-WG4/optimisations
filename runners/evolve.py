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
from models.field_models import CurrentBlock


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
        self.z_start = 0.0
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
        zout = [self.z_start, self.z_start+self.period/4, self.z_start+self.period]
        matrix_list = scipy.integrate.odeint(self.matrix_derivatives,
                                        [1, 0.0, 0.0, 1.0, 0.0], # m_00 m_01 m_10 m_11 larmor
                                        zout, hmax=self.hmax,
                                        mxstep=int(zout[-1]/self.hmax*2))
        assert(len(matrix_list) == len(zout))
        matrix = matrix_list[-1]
        larmor_angle = matrix[4]
        m2d = numpy.array([[matrix[0], matrix[1]], [matrix[2], matrix[3]]])
        cosmu = (m2d[0,0]+m2d[1,1])/2
        #print("Transfer matrix with cosmu", cosmu, "det", numpy.linalg.det(m2d), "\n", m2d)
        if abs(cosmu) > 1:
            return (0.0, 0.0, 0.0, 0.0)
        n2d = m2d - numpy.array([[cosmu, 0.0],[0.0,cosmu]])
        #print("N2D\n", n2d)
        # absolute value of sin(mu) from square root; we use beta is pos definite to force the sign of sinmu
        sinmu = numpy.linalg.det(n2d)**0.5*numpy.sign(n2d[0,1])
        n2d /= sinmu
        #print("N2D over sinmu with sinmu=", sinmu, "\n", n2d)
        v2d = numpy.array([[n2d[0,1], -n2d[0,0]], [-n2d[0,0], -n2d[1, 0]]])
        # v2d[0,0] = beta/p
        beta, alpha, phase_advance = v2d[0, 0]*self.momentum, -v2d[0,1], -math.atan2(cosmu, sinmu)

        #print(numpy.dot( numpy.dot(m2d, v2d), numpy.transpose(m2d)) - v2d)
        matrix = matrix_list[1]
        m2d = numpy.array([[matrix[0], matrix[1]], [matrix[2], matrix[3]]])
        vhalf = numpy.dot( numpy.dot(m2d, v2d), numpy.transpose(m2d))
        beta_half = vhalf[0, 0]*self.momentum
        #print("beta alpha phi", beta, alpha, phase_advance)
        return beta, alpha, phase_advance, beta_half

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
        return beta1, dbetadz, phi, 0.0

    def propagate_beta(self, beta0, dbetadz0, z0=0, z1=None, n_points=101):
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
        if z1 is None:
            z1 = self.period
        z_list = [z0+(z1-z0)*i/float(n_points-1) for i in range(n_points)]
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


def make_energy_axis(p_axes, is_kinetic):
    mass = 0.105865
    if is_kinetic:
        ke_calc = lambda p: (p**2+mass**2)**0.5-mass
    else:
        ke_calc = lambda p: (p**2+mass**2)**0.5
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


class Plotter:
    def __init__(self, field, pz_plot_list, pz_list, plot_dir, beta_max, nominal_emittance = -1):
        self.field = field
        self.pz_plot_list = pz_plot_list
        self.pz_list = pz_list
        self.plot_dir = plot_dir
        self.beta_max = beta_max
        self.nominal_emittance = nominal_emittance
        self.bsquared = None
        self.make_title()
        self.out_data = {}
        self.figure = None

    def make_title(self):
        human = self.field.human_readable()
        if human:
            human = "$"+human+"$\n"
        title = human+\
               "$\\int B^2(z) dz =$ {0:.4g} T$^2$ m".format(self.field.get_bz2_int())
        return title

    def get_name(self):
        return "optics_L_"+str(self.field.period)+"_"+self.field.get_name()


    def prepare_plots(self):
        self.out_data = {"pz_items":[], "field_name":self.field.get_name(),}
        try:
            self.out_data.update({"field_bi":self.field.bz_list, "field_b0":self.field.bz0, "period":self.field.get_period()})
        except:
            pass
        if self.figure == None:
            self.figure = matplotlib.pyplot.figure(figsize=(12, 10))
            self.axes = [
                self.figure.add_subplot(2, 2, 1),
                self.figure.add_subplot(2, 2, 2),
                self.figure.add_subplot(2, 2, 3),
                self.figure.add_subplot(2, 2, 4)]
            self.axes.append(self.axes[2].twinx())
            self.axes.append(self.axes[3].twinx())
            self.figure.suptitle(self.make_title())
        self.beta_finder = BetaFinder(self.field, 1., use_analytic=True)
        if self.bsquared != None:
            self.field.normalise_bz_squared(bsquared)
        self.z_list = [i*self.field.period/100 for i in range(101)]

    def write(self):
        fout = open(os.path.join(self.plot_dir, self.get_name()+".json"), "w")
        fout.write(json.dumps(self.out_data))
        self.figure.savefig(os.path.join(self.plot_dir, self.get_name()+".png"))
        print("Finished plotting lattice")


    def do_plots(self):
        self.prepare_plots()
        self.field_plot()
        self.field_squared_plot()
        self.plot_beta_vs_pz()
        self.plot_beta_vs_z()
        self.write()

    def field_plot(self):
        bz_list = [self.field.get_field(z) for z in self.z_list]
        self.axes[0].plot(self.z_list, bz_list)
        self.axes[0].set_xlabel("z [m]")
        self.axes[0].set_ylabel("B$_{z}$ [T]")

    def field_squared_plot(self):
        bz2_list = [self.field.get_field(z)**2 for z in self.z_list]
        self.axes[1].plot(self.z_list, bz2_list)
        self.axes[1].set_xlabel("z [m]")
        self.axes[1].set_ylabel("B$_{z}^2$ [T$^2$]")

    def plot_beta_vs_pz(self):
        print("     pz    beta_0     antibeta_0 phi      n_iterations")
        beta_list = []
        antibeta_list = []
        for pz in self.pz_list:
            self.beta_finder.momentum = pz
            beta, alpha, phi, beta_half = self.beta_finder.get_beta_periodic()
            print (f"    {pz:12.4g} {beta:12.4g} {beta_half:12.4g} {phi:12.4g}")
            beta_list.append(beta)
            antibeta_list.append(beta_half)
        self.out_data["pz_list"] = self.pz_list
        self.out_data["beta_list"] = beta_list
        self.out_data["antibeta_list"] = antibeta_list
        self.axes[2].plot([0.0, self.pz_list[-1]], [0.0, (beta_list[-1]+antibeta_list[-1])/2.0], c="xkcd:light grey", linestyle="dotted")
        self.axes[2].plot(self.pz_list, beta_list, label="$\\beta(L)$")
        self.axes[2].plot(self.pz_list, antibeta_list, 'g--', label="$\\beta(L/2)$")
        self.axes[2].set_xlabel("p$_{z}$ [GeV/c]")
        self.axes[2].set_ylabel("$\\beta$ [m]")
        self.axes[2].set_ylim([0.0, self.beta_max])

    def plot_beta_vs_z(self):
        for pz in self.pz_plot_list:
            self.beta_finder.momentum = pz
            beta, alpha, phi, beta_mid = self.beta_finder.get_beta_periodic()
            if beta < 1e-9:
                print(f"Failed to plot {pz}; no closed solution found")
                continue
            self.beta_finder.verbose = 0
            z_list, beta_list, dbetadz_list, phi_list = \
                                            self.beta_finder.propagate_beta(beta, 0.)
            self.axes[3].plot(z_list, beta_list, label="$\\beta_\\perp$ "+format(pz, "6.4g")+" GeV/c")
            self.out_data["pz_items"].append({"z":z_list, "beta_list":beta_list, "pz":pz})
        self.axes[3].set_xlabel("z [m]")
        self.axes[3].set_ylabel("$\\beta$ [m]")
        self.axes[3].set_ylim([0.0, self.beta_max])
        mass = 0.105658
        if self.nominal_emittance > 0:
            for pz in self.pz_plot_list:
                sigma_x_list = [(beta*self.nominal_emittance*mass/pz)**0.5 for beta in beta_list]
                self.axes[5].plot(z_list, sigma_x_list, linestyle="dashed", label=f"$\\sigma_x$ {pz} GeV/c")
        self.axes[5].set_ylabel("$\\sigma_x$ [m]")
        self.axes[5].set_ylim([0.0, (self.beta_max*self.nominal_emittance*mass/min(self.pz_plot_list))**0.5])
        try:
            self.axes[3].legend(loc="upper left")
            self.axes[5].legend(loc="upper right")
        except:
            pass

class EqmPlotter:
    def __init__(self):
        self.mass = 0.105658

    def dedx(self):
        pass

    def get_energy(self, pz):
        return (pz**2+self.mass**2)**0.5

    def get_g_l(self, g_t, pz):
        dow = 1-g_t
        gamma_rel = self.get_energy(pz)/self.mass
        beta_rel = pz/gamma_rel/self.mass
        I = 36.5e-9 # GeV
        m_e_c2 = 0.00051099895069 # GeV
        K = 2*m_e_c2/I
        g_l0  = -2/gamma_rel**2 + 2*(1-(beta_rel/gamma_rel)**2)/(log(K*beta_rel**2*gamma_rel**2)-beta_rel**2)
        g_l = g_l0-dow
        return g_l

    def plot_eqm_emittance(self, pz):
        beta_rel = 0.0
        figure = matplotlib.pyplot.figure("emittance")
        gt_list = [0.01*i for i in range(101)]
        gl_list = [self.get_g_l(gt, pz) for gt in gt_list]
        E_mu = self.get_energy(pz)
        beta_t = 0.1
        eps_t = beta_t*E_mu**2/2/beta/g_t/self.mass/L_r
        eps_l = beta_l*m_e*gamma**2*(1-0.0)

def main():
    global fignum
    plot_dir = "output/demonstrator_optics_v3"
    c_x = 0.8
    c_b = 1/0.8
    c_p = 1.0
    pz_plot_list = [c_p*0.19, c_p*0.20, c_p*0.21]
    pz_scan_list = [c_p*pz_i/1000. for pz_i in range(140, 301, 2)]
    n_points = 2
    norm = 100
    #clear_dir(plot_dir)
    b0 = c_b*0.0
    b1 = c_b*7.0
    b2 = c_b*1.0
    b3 = c_b*0.0
    b4 = c_b*0.0
    length = c_x*1.0
    for b2 in [0.0, 0.4, 1.0]:
        b2 = b2*c_b
        solenoid_field = SineField(b0, b1, b2, b3, b4, 0.0, length)
        solenoid_field.normalise_bz_squared(31.25)
        #solenoid_field.normalise_bz_squared(norm)
        #solenoid_field = CurrentSheet(108.02, 0.1655, 0.1, 0.25, 1.0, 1)
        #def __init__(self, b0, zcentre, length, rmin, rmax, period, nrepeats, nsheets):
        plotter = Plotter(solenoid_field, pz_plot_list, pz_scan_list, plot_dir, nominal_emittance = 0.0025, beta_max=0.8)
        plotter.do_plots()

    matplotlib.pyplot.show(block=False)

if __name__ == "__main__":
    main()
    input("Press <CR> to end")
