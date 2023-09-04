import os
import shutil

import ctypes
import matplotlib
import matplotlib.pyplot

import ROOT

import evolve
import field_models
import evolve

class FitCoils(object):
    """
    Fit a lattice for final cooling
    """
    def __init__(self, pz, long_field, peak_field):
        """
        Initialise
        - pz: [GeV/c] the momentum of the cooling cell
        - long_field: [T] the field strength of the long, low field solenoid
        - peak_field: [T] the field strength of the shorter, high field solenoid
        """
        self.long_field = long_field
        self.peak_field = peak_field
        self.m1 = 0
        self.m2 = 0
        self.m3 = 0
        self.z0 = 0
        self.z1 = 0
        self.pz = pz
        self.n_iterations_per_fit = 1000
        self.lattice = None
        self.make_fields()

    def make_fields(self):
        """
        Set up the lattice with some nominal geometry parameters hard coded
        - match_field_1: [T] the field strength of the match coil between the 
                         field flip and the long solenoid
        - match_field_2: [T] the field strength of the first match coil between the
                         long solenoid and the high field solenoid  
        - match_field_3: [T] the field strength of the second match coil between the
                         long solenoid and the high field solenoid  
        """
        period = 40.0
        nrepeats = 20
        nsheets = 1
        field_list = []
        match_field_1, match_field_2, match_field_3 = self.m1, self.m2, self.m3
        length_m1, rmin_m1, rmax_m1, dz_m1 = 0.5, 0.2, 0.3, 1.0
        length_m2, rmin_m2, rmax_m2, dz_m2 = 0.5, 0.2, 0.3, period/4.0-2.0
        length_m3, rmin_m3, rmax_m3, dz_m3 = 0.5, 0.2, 0.3, period/4.0-2.5
        for (b0, zcentre, length, rmin, rmax) in [
                (match_field_1, dz_m1, length_m1, rmin_m1, rmax_m1),
                (match_field_2, dz_m2, length_m2, rmin_m2, rmax_m2),
                (match_field_3, dz_m3, length_m3, rmin_m3, rmax_m3),
                (self.long_field, period/4.0, 19.0, 0.2, 0.3),
                (self.peak_field-self.long_field, period/4.0, 2.0, 0.1, 0.2),
                (match_field_3, period/2-dz_m3, length_m3, rmin_m3, rmax_m3),
                (match_field_2, period/2-dz_m2, length_m2, rmin_m2, rmax_m2),
                (match_field_1, period/2-dz_m1, length_m1, rmin_m1, rmax_m1),
                (-match_field_1, period/2+dz_m1, length_m1, rmin_m1, rmax_m1),
                (-match_field_2, period/2+dz_m2, length_m2, rmin_m2, rmax_m2),
                (-match_field_3, period/2+dz_m3, length_m3, rmin_m3, rmax_m3),
                (-self.long_field, 3.0*period/4.0, 19.0, 0.2, 0.3),
                (-self.peak_field+self.long_field, 3.0*period/4.0, 2.0, 0.1, 0.2),
                (-match_field_3, period-dz_m3, length_m3, rmin_m3, rmax_m3),
                (-match_field_2, period-dz_m2, length_m2, rmin_m2, rmax_m2),
                (-match_field_1, period-dz_m1, length_m1, rmin_m1, rmax_m1),
            ]:
            coil = field_models.CurrentBlock(b0, zcentre, length, rmin, rmax, period, nrepeats, nsheets)
            field_list.append(coil)
        self.lattice = field_models.FieldSum(field_list)
        self.beta_finder = evolve.BetaFinder(self.lattice, self.pz)

    def set_fields_minuit(self):
        """Extract fields from minuit and set them using self.make_fields(...)"""
        match_field_1, match_field_2, match_field_3 = ctypes.double(), ROOT.Double(), ROOT.Double()
        err = ROOT.Double()
        self.minuit.GetParameter(0, match_field_1, err)
        self.minuit.GetParameter(1, match_field_2, err)
        self.minuit.GetParameter(2, match_field_3, err)
        self.m1 = float(match_field_1)
        self.m2 = float(match_field_2)
        self.m3 = float(match_field_3)
        self.make_fields()
        return match_field_1, match_field_2, match_field_3

    def uniform_field_beta(self, z_pos):
        """
        Return the beta having alpha = 0 for a uniform field given by the field
        at z = z_pos
        """
        beta0 = abs(self.beta_finder.momentum/self.lattice.get_field(z_pos)/0.15)
        return beta0


    def fit_2(self):
        """
        Second we match from the constant field region to the high field region
        to find match_field_2 and match_field_3
        """
        global global_self
        # set the z position of the start (z0) and end (z1) of the integration
        self.z0 = self.lattice.period*2.0/8.0 # peak field
        self.z1 = self.lattice.period*3.0/8.0 # long solenoid field
        if self.minuit:
            m1, m2, m3 = self.set_fields_minuit()
        self.minuit = ROOT.TMinuit()
        self.minuit.SetPrintLevel(-1)
        self.minuit.DefineParameter(0, "match_field_1", m1, 0.01, -10, 10)
        self.minuit.DefineParameter(1, "match_field_2", m2, 0.01, -10, 10)
        self.minuit.DefineParameter(2, "match_field_3", m3, 0.01, -10, 10)
        self.minuit.FixParameter(0)
        #self.minuit.FixParameter(1)
        #self.minuit.FixParameter(2)
        global_self = self
        self.minuit.SetFCN(self.score_function_2)
        self.minuit.Command("SIMPLEX "+str(self.n_iterations_per_fit)+" "+str(1e-16))
        beta1 = self.score_function_2(0, 0, [0], 0, 0)
        return beta1
        
    def score_function_2(self, nvar, parameters, score, jacobian, err):
        """
        Score function for second fit - we evolve beta from long solenoid to the
        high field region and require uniform beta function at the centre of the
        high field region (should be set to z1) 
        """
        match_field_1, match_field_2, match_field_3 = self.set_fields_minuit()
        beta0 = self.uniform_field_beta(self.z0)
        beta1_tgt = self.uniform_field_beta(self.z1)
        beta1, dbeta1dz, phi = self.beta_finder.evolve(float(beta0), 0., self.z0, self.z1)
        score[0] = (beta1-beta1_tgt)**2+(10*dbeta1dz)**2
        print("Match1: ", match_field_1, "Match2: ", match_field_2, "Match3: ", match_field_3,
              "beta0", beta0, "beta1", beta1,
              "beta1_tgt", beta1_tgt, "dbeta1dz", dbeta1dz, "score", score[0])
        return beta1, dbeta1dz

    def fit_1(self):
        """
        First we match from the constant field region to the end of the cell, to
        find match_field_1
        """
        # set the z position of the start (z0) and end (z1) of the integration
        self.z0 = self.lattice.period*7.0/8.0 # half way between 0.0 and peak field
        self.z1 = self.lattice.period*8.0/8.0 # first peak field
        self.minuit = ROOT.TMinuit()
        self.minuit.SetPrintLevel(-1)
        self.minuit.DefineParameter(0, "match_field_1", self.m1_seed,  0.01, -10, 10)
        self.minuit.DefineParameter(1, "match_field_2", self.m2_seed, 0.01, -10, 10)
        self.minuit.DefineParameter(2, "match_field_3", self.m3_seed,  0.01, -10, 10)
        #self.minuit.FixParameter(0)
        self.minuit.FixParameter(1)
        self.minuit.FixParameter(2)
        self.minuit.SetFCN(self.score_function_1)
        self.minuit.Command("SIMPLEX "+str(self.n_iterations_per_fit)+" "+str(1e-16))
        beta1 = self.score_function_1(0, 0, [0], 0, 0)
        return beta1
        
    def score_function_1(self, nvar, parameters, score, jacobian, err):
        """
        Score function for first fit - we evolve beta from long solenoid to the
        field flip and require derivative of beta function wrt z at the centre
        of the high field region is 0
        """
        match_field_1, match_field_2, match_field_3 = self.set_fields_minuit()
        beta0 = self.uniform_field_beta(self.z0)
        beta1, dbeta1dz, phi = self.beta_finder.evolve(float(beta0), 0., self.z0, self.z1)
        score[0] = dbeta1dz**2
        print("Match1: ", match_field_1, "Match2: ", match_field_2, "Match3: ", match_field_3,
              "beta0", beta0, "beta1", beta1,
              "dbeta1dz", dbeta1dz, "score", score[0])
        return beta1, dbeta1dz

    def plot_beta(self, beta0, dbetadz0):
        """
        Plot beta as a function of z
        """
        z_list, beta_list, dbetadz_list, phi_list = self.beta_finder.propagate_beta(beta0, dbetadz0, 1001)
        m1, m2, m3 = self.set_fields_minuit()
        figure1 = matplotlib.pyplot.figure()
        axes = figure1.add_subplot(1, 1, 1)
        axes.plot(z_list, beta_list)
        axes.set_xlabel("z [m]")
        axes.set_ylabel("$\\beta$ [m]")
        axes.set_ylim([0.0, axes.get_ylim()[1]])
        axes = axes.twinx()
        axes.plot(z_list, [self.beta_finder.field(z) for z in z_list], color='green', linestyle="--")
        axes.set_xlabel("z [m]")
        axes.set_ylabel("B [T]")
        title = f"M1={m1:.3g} [T]; M2={m2:.3g} [T]; M3={m3:.3g} [T]\n"
        title += f"B$_{{const}}$={self.long_field:.3g} [T]; B$_{{peak}}$={self.peak_field:.3g} [T]"
        title += f"$\\beta_0$={beta_list[0]:.5g} [m]; $\\beta_1$={beta_list[-1]:.5g} [m] p$_z$={self.pz:.5g}"
        figure1.suptitle(title)
        figure2 = matplotlib.pyplot.figure()
        axes = figure2.add_subplot(1, 1, 1)
        i_list = [i for i in range(len(z_list)) if z_list[i] > 9.0 and z_list[i] < 11.0]
        z_zoom_list = [z_list[i] for i in i_list]
        beta_zoom_list = [beta_list[i] for i in i_list]
        axes.plot(z_zoom_list, beta_zoom_list)
        axes.set_xlabel("z [m]")
        axes.set_ylabel("$\\beta$ [m]")
        axes.set_ylim([0.0, axes.get_ylim()[1]])

        return figure1, figure2

def clear_dir(a_dir):
    try:
        shutil.rmtree(a_dir)
    except OSError:
        pass
    os.makedirs(a_dir)

def get_pz(pz, delta_ek):
    mass = 0.105658
    etot = (pz**2+mass**2)**0.5
    pz_new = ((etot+delta_ek)**2-mass**2)**0.5
    return pz_new


def main():
    plot_dir = "final_cooling_matching_v2"
    clear_dir(plot_dir)
    pz = 0.07 # GeV/c
    print("Pz spread", get_pz(pz, 0.0), get_pz(pz, -0.0038), get_pz(pz, 0.0038))
    long_field = 4.0
    peak_field = 30.0
    fitter = FitCoils(pz, long_field, peak_field)
    fitter.n_iterations_per_fit = 1000
    fitter.m1_seed = -7.410838204077985
    fitter.m2_seed = 1.9161292750681458
    fitter.m3_seed = 5.371882228095956
    print("Fit 1")
    beta, dbetadz = fitter.fit_1()
    print("Fit 2")
    fitter.fit_2()
    for fi in range(-3, 4):
        fitter.peak_field = 30.0 # + 10.0*fi
        fitter.pz = pz + fi*0.05
        print("Doing optics with peak field", fitter.peak_field)
        fitter.make_fields()
        figure1, figure2 = fitter.plot_beta(beta, -dbetadz)
        figure1.savefig(os.path.join(plot_dir, f"matched_lattice_{fi}.png"))
        figure2.savefig(os.path.join(plot_dir, f"matched_lattice_{fi}_zoom.png"))
    matplotlib.pyplot.show(block=False)

if __name__ == "__main__":
    main()
    input("Press <CR> to end")
