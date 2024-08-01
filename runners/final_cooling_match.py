import os
import shutil

import ctypes
import matplotlib
import matplotlib.pyplot

import ROOT

import runners.evolve
import models.field_models
import g4bl_interface.g4bl_field_model_wrapper

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
        self.m4 = 0
        self.fit1_z0 = 0
        self.fit1_z1 = 0
        self.pz = pz
        self.n_iterations_per_fit = 1000
        self.lattice = None
        self.minuit = ROOT.TMinuit()
        self.minuit.SetPrintLevel(-1)
        self.make_fields()

    def make_fields(self):
        """
        Set up the lattice with some nominal geometry parameters hard coded
        - match_field_1: [T] the field strength of the first match coil between the
                         field flip and the long solenoid
        - match_field_2: [T] the field strength of the second match coil between the
                         field flip and the long solenoid
        - match_field_3: [T] the field strength of the first match coil between the
                         long solenoid and the high field solenoid  
        - match_field_4: [T] the field strength of the second match coil between the
                         long solenoid and the high field solenoid  
        """
        period = 40.0
        nrepeats = 20
        nsheets = 1
        field_list = []
        match_field_1, match_field_2, match_field_3, match_field_4 = self.m1, self.m2, self.m3, self.m4
        length_m1, rmin_m1, rmax_m1, dz_m1 = 0.5, 0.2, 0.3, 1.0
        length_m2, rmin_m2, rmax_m2, dz_m2 = 0.5, 0.2, 0.3, 1.5
        length_m3, rmin_m3, rmax_m3, dz_m3 = 0.5, 0.2, 0.3, period/4.0-2.0
        length_m4, rmin_m4, rmax_m4, dz_m4 = 0.5, 0.2, 0.3, period/4.0-2.5
        for (b0, zcentre, length, rmin, rmax) in [
                (match_field_1, dz_m1, length_m1, rmin_m1, rmax_m1),
                (match_field_2, dz_m2, length_m2, rmin_m2, rmax_m2),
                (match_field_3, dz_m3, length_m3, rmin_m3, rmax_m3),
                (match_field_4, dz_m4, length_m4, rmin_m4, rmax_m4),
                (self.long_field, period/4.0, 19.0, 0.2, 0.3),
                (self.peak_field-self.long_field, period/4.0, 2.0, 0.1, 0.2),
                (match_field_4, period/2-dz_m4, length_m4, rmin_m4, rmax_m4),
                (match_field_3, period/2-dz_m3, length_m3, rmin_m3, rmax_m3),
                (match_field_2, period/2-dz_m2, length_m2, rmin_m2, rmax_m2),
                (match_field_1, period/2-dz_m1, length_m1, rmin_m1, rmax_m1),
                (-match_field_1, period/2+dz_m1, length_m1, rmin_m1, rmax_m1),
                (-match_field_2, period/2+dz_m2, length_m2, rmin_m2, rmax_m2),
                (-match_field_3, period/2+dz_m3, length_m3, rmin_m3, rmax_m3),
                (-match_field_4, period/2+dz_m4, length_m4, rmin_m4, rmax_m4),
                (-self.long_field, 3.0*period/4.0, 19.0, 0.2, 0.3),
                (-self.peak_field+self.long_field, 3.0*period/4.0, 2.0, 0.1, 0.2),
                (-match_field_4, period-dz_m4, length_m4, rmin_m4, rmax_m4),
                (-match_field_3, period-dz_m3, length_m3, rmin_m3, rmax_m3),
                (-match_field_2, period-dz_m2, length_m2, rmin_m2, rmax_m2),
                (-match_field_1, period-dz_m1, length_m1, rmin_m1, rmax_m1),
            ]:
            coil = models.field_models.CurrentBlock(1.0, zcentre-period*int(nrepeats/2), length, rmin, rmax, period, nrepeats, nsheets)
            coil.set_peak_b0(b0)
            field_list.append(coil)
        self.lattice = models.field_models.FieldSum(field_list)
        self.beta_finder = runners.evolve.BetaFinder(self.lattice, self.pz)

    def set_fields_minuit(self):
        """Extract fields from minuit and set them using self.make_fields(...)"""
        match_field_1, match_field_2 = ctypes.c_double(), ctypes.c_double()
        match_field_3, match_field_4 = ctypes.c_double(), ctypes.c_double()
        err = ctypes.c_double()
        self.minuit.GetParameter(0, match_field_1, err)
        self.minuit.GetParameter(1, match_field_2, err)
        self.minuit.GetParameter(2, match_field_3, err)
        self.minuit.GetParameter(3, match_field_4, err)
        self.m1 = float(match_field_1.value)
        self.m2 = float(match_field_2.value)
        self.m3 = float(match_field_3.value)
        self.m4 = float(match_field_4.value)
        self.make_fields()
        return match_field_1, match_field_2, match_field_3, match_field_4

    def uniform_field_beta(self, z_pos):
        """
        Return the beta having alpha = 0 for a uniform field given by the field
        at z = z_pos
        """
        beta0 = abs(self.beta_finder.momentum/self.lattice.get_field(z_pos)/0.15)
        return beta0

    def setup_minuit(self):
        """
        Setup a new minuit instance with default values
        """
        self.minuit = ROOT.TMinuit()
        self.minuit.SetPrintLevel(-1)
        self.minuit.DefineParameter(0, "match_field_1", self.m1, 0.01, -4, 4)
        self.minuit.DefineParameter(1, "match_field_2", self.m2, 0.01, -1, 1)
        self.minuit.DefineParameter(2, "match_field_3", self.m3, 0.01, -10, 10)
        self.minuit.DefineParameter(3, "match_field_4", self.m4, 0.01, -10, 10)

    def fit_1(self):
        """
        First we match from the constant field region to the end of the cell, to
        find match_field_1
        """
        # set the z position of the start (z0) and end (z1) of the integration
        self.z0 = self.lattice.period*7.0/8.0 # half way between 0.0 and peak field
        self.z1 = self.lattice.period*8.0/8.0 # first peak field
        self.setup_minuit()

        self.minuit.FixParameter(0)
        self.minuit.FixParameter(1)
        self.minuit.FixParameter(2)
        self.minuit.FixParameter(3)

        # ROOT can't take member function as argument to minuit so we need to
        # define a module function
        FitCoils.my_fit_coils = self
        self.minuit.SetFCN(global_score_function_1)
        # Run the minimiser
        self.minuit.Command("SIMPLEX "+str(self.n_iterations_per_fit)+" "+str(1e-16))
        # finally, get the score one more time
        beta1, dbeta1dz = self.score_function_1(0, 0, ctypes.c_double(0.0), 0, 0)
        return beta1

    def score_function_1(self, nvar, parameters, score, jacobian, err):
        """
        Score function for first fit - we evolve beta from long solenoid to the
        field flip and require derivative of beta function wrt z at the centre
        of the high field region is 0
        """
        match_field_1, match_field_2, match_field_3, match_field_4 = self.set_fields_minuit()
        beta0 = self.uniform_field_beta(self.z0)
        beta1, dbeta1dz, phi = self.beta_finder.evolve(float(beta0), 0., self.z0, self.z1)
        score.value = dbeta1dz**2
        print("Match1: ", match_field_1.value, "Match2: ", match_field_2.value,
              "Match3: ", match_field_3.value, "Match4: ", match_field_4.value,
              "beta0", beta0, "beta1", beta1,
              "dbeta1dz", dbeta1dz, "score", score.value)
        return beta1, dbeta1dz

    def fit_2(self):
        """
        Second we match from the constant field region to the high field region
        to find match_field_2 and match_field_3
        """
        global global_self
        # set the z position of the start (z0) and end (z1) of the integration
        self.z0 = self.lattice.period*2.0/8.0 # peak field
        self.z1 = self.lattice.period*3.0/8.0 # long solenoid field
        self.setup_minuit()
        self.minuit.FixParameter(0)
        self.minuit.FixParameter(1)
        self.minuit.FixParameter(2)
        self.minuit.FixParameter(3)

        # ROOT can't take member function as argument to minuit so we need to
        # define a module function
        FitCoils.my_fit_coils = self
        self.minuit.SetFCN(global_score_function_2)
        self.minuit.Command("SIMPLEX "+str(self.n_iterations_per_fit)+" "+str(1e-16))
        beta1 = self.score_function_2(0, 0, ctypes.c_double(), 0, 0)
        return beta1

    def score_function_2(self, nvar, parameters, score, jacobian, err):
        """
        Score function for second fit - we evolve beta from long solenoid to the
        high field region and require uniform beta function at the centre of the
        high field region (should be set to z1)
        """
        match_field_1, match_field_2, match_field_3, match_field_4 = self.set_fields_minuit()
        beta0 = self.uniform_field_beta(self.z0)
        beta1_tgt = self.uniform_field_beta(self.z1)
        beta1, dbeta1dz, phi = self.beta_finder.evolve(float(beta0), 0., self.z0, self.z1)
        score.value = (beta1-beta1_tgt)**2+(10*dbeta1dz)**2
        print("Match1: ", match_field_1.value, "Match2: ", match_field_2.value,
              "Match3: ", match_field_3.value, "Match4: ", match_field_4.value, "  **  "
              "beta0", beta0, "beta1", beta1,
              "beta1_tgt", beta1_tgt, "dbeta1dz", dbeta1dz, "score", score.value)
        return beta1, dbeta1dz

    def plot_beta(self, beta0, dbetadz0, z0, z1):
        """
        Plot beta as a function of z
        """
        z_list, beta_list, dbetadz_list, phi_list = self.beta_finder.propagate_beta(beta0, dbetadz0, z0, z1, 1001) # 1001 data points
        figure1 = matplotlib.pyplot.figure()
        axes = figure1.add_subplot(1, 1, 1)
        axes.plot(z_list, beta_list)
        axes.set_xlabel("z [m]")
        axes.set_ylabel("$\\beta$ [m]")
        axes.set_xlim([0.0, 40.0])
        axes.set_ylim([0.0, 0.5])
        axes = axes.twinx()
        axes.plot(z_list, [self.beta_finder.field(z) for z in z_list], color='green', linestyle="--")
        axes.set_xlabel("z [m]")
        axes.set_ylabel("B [T]")
        title = f"M1={self.m1:.5g} [T]; M2={self.m2:.5g} [T]; M3={self.m3:.5g} [T]; M4={self.m4:.5g} [T]\n"
        title += f"B$_{{const}}$={self.long_field:.3g} [T]; B$_{{peak}}$={self.peak_field:.3g} [T]"
        title += f"$\\beta_0$={beta_list[0]:.5g} [m]; $\\beta_1$={beta_list[-1]:.5g} [m] p$_z$={self.pz:.5g}"
        figure1.suptitle(title)

        return figure1

    def make_g4bl_lattice(self, file_name):
        g4bl_lattice = g4bl_interface.G4BLLinac(filename)
        for current_block in self.lattice:
            g4bl_field_list = g4bl_field_model_wrapper.CurrentBlockG4BL(current_block)
            g4bl_lattice.elements_list += g4bl_field_list
        

    my_fit_coils = None

def global_score_function_1(nvar, parameters, score, jacobian, err):
    return FitCoils.my_fit_coils.score_function_1(nvar, parameters, score, jacobian, err)

def global_score_function_2(nvar, parameters, score, jacobian, err):
    return FitCoils.my_fit_coils.score_function_2(nvar, parameters, score, jacobian, err)


def clear_dir(a_dir):
    try:
        shutil.rmtree(a_dir)
    except OSError:
        pass
    os.makedirs(a_dir)

def main():
    plot_dir = "output/final_cooling_matching_test"
    clear_dir(plot_dir)
    pz = 0.07 # GeV/c
    long_field = 4.0 # T
    peak_field = 30.0 # T
    # set up the fitter
    fitter = FitCoils(pz, long_field, peak_field)
    fitter.n_iterations_per_fit = 100
    # WARNING there is a minimum at +2.8 T which does not provide a match - MINUIT will land here if we start m1 = 0.0
    fitter.m1 = -3.708232350 # T
    fitter.m2 = 0.0 # T
    fitter.m3 = 0.04453909 # T
    fitter.m4 = -0.069320357 # T
    # first do the match from the low field region across the field flip
    print("Doing first fit - through the field flip")
    beta1 = fitter.fit_1()
    # second do the match from the high field region to the low field region
    print("Doing second fit - peak field to low field")
    beta2 = fitter.fit_2()
    # finally plot
    print("Plotting")
    # fitter.plot_beta(initial_beta, initial_dbetadz, z0, z1)
    figure1 = fitter.plot_beta(beta1, 0.0, 0.0, 40.0)
    fname = os.path.join(plot_dir, f"matched_lattice.png")
    figure1.savefig(fname)
    print("Saved to", fname)
    matplotlib.pyplot.show(block=False)

if __name__ == "__main__":
    main()
    input("Press <CR> to end")
