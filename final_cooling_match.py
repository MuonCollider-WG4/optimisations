import os
import shutil

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
        self.z0 = 0
        self.z1 = 0
        self.pz = pz
        self.lattice = None
        self.make_fields(0, 0)

    def make_fields(self, match_field_1, match_field_2):
        """
        Set up the lattice with some nominal geometry parameters hard coded
        - match_field_1: [T] the field strength of the match coil between the 
                         field flip and the long solenoid
        - match_field_2: [T] the field strength of the match coil between the
                         long solenoid and the high field solenoid  
        """
        period = 40.0
        nrepeats = 20
        nsheets = 1
        field_list = []
        length_m1, rmin_m1, rmax_m1, dz_m1 = 0.5, 0.2, 0.3, 0.5
        length_m2, rmin_m2, rmax_m2, dz_m2 = 0.5, 0.2, 0.3, period/4.0-1.0
        for (b0, zcentre, length, rmin, rmax) in [
                (match_field_1, dz_m1, length_m1, rmin_m1, rmax_m1),
                (match_field_2, dz_m2, length_m2, rmin_m2, rmax_m2),
                (self.long_field, period/4.0, 19.0, 0.2, 0.3),
                (self.peak_field, period/4.0, 1.0, 0.1, 0.2),
                (match_field_2, period/2-dz_m2, length_m2, rmin_m2, rmax_m2),
                (match_field_1, period/2-dz_m1, length_m1, rmin_m1, rmax_m1),
                (-match_field_1, period/2+dz_m1, length_m1, rmin_m1, rmax_m1),
                (-match_field_2, period/2+dz_m2, length_m2, rmin_m2, rmax_m2),
                (-self.long_field, 3.0*period/4.0, 19.0, 0.2, 0.3),
                (-self.peak_field, 3.0*period/4.0, 1.0, 0.1, 0.2),
                (-match_field_2, period-dz_m2, length_m2, rmin_m2, rmax_m2),
                (-match_field_1, period-dz_m1, length_m1, rmin_m1, rmax_m1),
            ]:
            coil = field_models.CurrentBlock(b0, zcentre, length, rmin, rmax, period, nrepeats, nsheets)
            field_list.append(coil)
        self.lattice = field_models.FieldSum(field_list)
        self.beta_finder = evolve.BetaFinder(self.lattice, self.pz)

    def set_fields_minuit(self):
        """Extract fields from minuit and set them using self.make_fields(...)"""
        match_field_1, match_field_2 = ROOT.Double(), ROOT.Double()
        err = ROOT.Double()
        self.minuit.GetParameter(0, match_field_1, err)
        self.minuit.GetParameter(1, match_field_2, err)
        self.make_fields(float(match_field_1), float(match_field_2))
        return match_field_1, match_field_2

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
        to find match_field_2
        """
        # set the z position of the start (z0) and end (z1) of the integration
        self.z0 = self.lattice.period*1.0/8.0 # half way between peak field and end field
        self.z1 = self.lattice.period*2.0/8.0 # end of the cell
        if self.minuit:
            m1, m2 = self.set_fields_minuit()
        self.minuit = ROOT.TMinuit()
        self.minuit.SetPrintLevel(-1)
        self.minuit.DefineParameter(0, "match_field_1", m1, 0.5, -10, 10)
        self.minuit.DefineParameter(1, "match_field_2", -1.0, 0.1, -10, 10)
        self.minuit.FixParameter(0)
        self.minuit.SetFCN(self.score_function_2)
        self.minuit.Command("SIMPLEX "+str(1000)+" "+str(1e-12))
        beta1 = self.score_function_2(0, 0, [0], 0, 0)
        return beta1
        
    def score_function_2(self, nvar, parameters, score, jacobian, err):
        """
        Score function for second fit - we evolve beta from long solenoid to the
        high field region and require uniform beta function at the centre of the
        high field region (should be set to z1) 
        """
        match_field_1, match_field_2 = self.set_fields_minuit()
        beta0 = self.uniform_field_beta(self.z0)
        beta1_tgt = self.uniform_field_beta(self.z1)
        beta1, dbeta1dz, phi = self.beta_finder.evolve(float(beta0), 0., self.z0, self.z1)
        score[0] = (beta1-beta1_tgt)**2+(100*dbeta1dz)**2
        print("Match1: ", match_field_1, "Match2: ", match_field_2,
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
        self.minuit.DefineParameter(0, "match_field_1", -6.571660086194663, 0.5, -10, 10)
        self.minuit.DefineParameter(1, "match_field_2", 0, 0.1, -10, 10)
        self.minuit.FixParameter(0)
        self.minuit.FixParameter(1)
        self.minuit.SetFCN(self.score_function_1)
        self.minuit.Command("SIMPLEX "+str(1)+" "+str(1e-12))
        beta1 = self.score_function_1(0, 0, [0], 0, 0)
        return beta1
        
    def score_function_1(self, nvar, parameters, score, jacobian, err):
        """
        Score function for first fit - we evolve beta from long solenoid to the
        field flip and require derivative of beta function wrt z at the centre
        of the high field region is 0
        """
        match_field_1, match_field_2 = self.set_fields_minuit()
        beta0 = self.uniform_field_beta(self.z0)
        beta1, dbeta1dz, phi = self.beta_finder.evolve(float(beta0), 0., self.z0, self.z1)
        score[0] = dbeta1dz**2
        print("Match1: ", match_field_1, "Match2: ", match_field_2,
              "beta0", beta0, "beta1", beta1,
              "dbeta1dz", dbeta1dz, "score", score[0])
        return beta1, dbeta1dz


    def plot_beta(self, beta0, dbetadz0):
        """
        Plot beta as a function of z
        """
        z_list, beta_list, dbetadz_list, phi_list = self.beta_finder.propagate_beta(beta0, dbetadz0, 1001)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 2, 1)
        axes.plot(z_list, beta_list)
        axes.set_xlabel("z [m]")
        axes.set_ylabel("$\\beta$ [m]")
        axes.set_xlim([0.0, axes.get_xlim()[1]])
        axes = figure.add_subplot(1, 2, 2)
        axes.plot(z_list, [self.beta_finder.field(z) for z in z_list])
        axes.set_xlabel("z [m]")
        axes.set_ylabel("B [T]")
        return figure

def match_fields(field, pz):
    pass

def clear_dir(a_dir):
    try:
        shutil.rmtree(a_dir)
    except OSError:
        pass
    os.makedirs(a_dir)

def main():
    plot_dir = "final_cooling_matching_v1"
    clear_dir(plot_dir)
    pz = 0.2 # GeV/c
    fitter = FitCoils(0.2, 4.0, 0.0)
    beta, dbetadz = fitter.fit_1()
    fitter.plot_beta(beta, -dbetadz)
    fitter.peak_field = 30.0
    fitter.fit_2()
    fitter.plot_beta(beta, -dbetadz)
    matplotlib.pyplot.show(block=False)

if __name__ == "__main__":
    main()
    input("Press <CR> to end")
