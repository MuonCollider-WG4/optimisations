import os
import shutil
import json
import matplotlib
import matplotlib.pyplot
import ROOT
from field_models import SineField
from field_models import CurrentSheet
from field_models import CurrentBlock

class CoilFitter(object):
    """
    Routines to fit a bunch of current sheets (aka coils) to some aribtrary on-axis field
    """
    def __init__(self, ri = 0.2):
        """
        Initialisation. 
        """
        self.n_iterations = 10000
        self.n_sheets_per_coil = 1
        self.cell_extent = 3 # n in each direction beyond the first
        self.n_fit_points = 100
        self.period = -1
        self.fit_params = [] # value:seed
        self.limits = []
        self.force_symmetry = -1 # -1, 0 or 1 for antisymmetric, no symmetry, symmetric
        #(b0, zcentre, length, rmin, rmax, period, nrepeats, nsheets)
        self.suffix = "_ri="+str(ri)
        self.coil_list = []
        self.minuit = None
        self.field_to_match = None
        self.iteration = 0

    @classmethod
    def new_from_pixels(self, ri, dr, nr, zi, dz, nz, suffix):
        """
        Initialise a new CoilFitter based on a rectangular set of current pixels
        - ri: inner radius of the grid [m]
        - dr: radial step [m]
        - nr: number of radial steps
        - zi: initial z position of the grid [m]
        - dz: z step [m]
        - nz: number of z steps
        """
        fitter = CoilFitter()
        fitter.coil_list = []
        for rindex in range(nr):
            for zindex in range(nz):           
                fitter.coil_list.append(
                    CurrentBlock(1.0, zi+dz*(zindex+0.5), dz, ri+dr*rindex, ri+dr*(rindex+1), 2.0, 1, 10)
                )
        if suffix[0] != "_":
            suffix = "_"+suffix
        fitter.suffix = suffix
        return fitter

    def fit_coil(self, field_to_match, tolerance):
        self.field_to_match = field_to_match
        self.iteration = 0
        self.period = field_to_match.period
        self.fit_params = ["b0"]#, "zcentre", "length"] # value:seed
        self.limits = [(0, 0)]#, (0.1, self.period/2.0*0.9), (0.1, self.period/2.0*0.9)]
        self.minuit = ROOT.TMinuit()
        param_index = 0
        for coil in self.coil_list:
            for i, key in enumerate(self.fit_params):
                value = coil.__dict__[key]
                print("Setting limits", self.limits[i])
                pmin = min(self.limits[i])
                pmax = max(self.limits[i])
                self.minuit.DefineParameter(param_index, key, value, max(abs(value)/10.0, 0.1), pmin, pmax)
                param_index += 1
        self.minuit.SetPrintLevel(-1)
        self.minuit.SetFCN(self.score_function)
        self.minuit.Command("SIMPLEX "+str(self.n_iterations)+" "+str(tolerance))
        self.print_coil_params()
        self.save_coil_params(os.path.join(self.output_dir, "coil_params"+self.suffix+".out"))

    def print_coil_params(self):
        for i, coil in enumerate(self.coil_list):
            print("Coil parameters for coil "+str(i))
            print("  nominal field [T]:       ", format(coil.b0, "9.4g"))
            print("  current density [A/m^2]:", format(coil.get_current_density(), "9.4g"))
            print("  centre position [m]:    ", format(coil.zcentre, "9.4g"))
            print("  length [m]:             ", format(coil.length, "9.4g"))
            print("  inner radius [m]:       ", format(coil.rmin, "9.4g"))
            print("  outer radius [m]:       ", format(coil.rmax, "9.4g"))
            print("  hoop stress:            ", "NEEDS CALC")

    def save_coil_params(self, file_name):
        # note G4BL units are *mm* hence all the constants
        mm_units = 1e3
        param_dict = {
            "__solenoid_current":lambda coil: coil.get_current_density()/mm_units**2,
            "__solenoid_z_pos":lambda coil: coil.zcentre*mm_units,
            "__solenoid_z_length":lambda coil: coil.length*mm_units,
            "__solenoid_inner_radius":lambda coil: coil.rmin*mm_units,
            "__solenoid_outer_radius":lambda coil: coil.rmax*mm_units,
        }
        coil_dict = {}
        i = 1
        for coil in self.coil_list:
            i_str = "_"+str(i)+"__"
            for key, value_lambda in param_dict.items():
                coil_dict[key+i_str] = value_lambda(coil)
            i += 1
            i_str = "_"+str(i)+"__"
            if True or self.force_symmetry == 0:
                continue
            for key, value_lambda in param_dict.items():
                if "z_pos" in key:
                    value = self.period*1e3-value_lambda(coil)
                elif "current" in key:
                    value = self.force_symmetry*value_lambda(coil)
                else:
                    value = value_lambda(coil)
                coil_dict[key+i_str] = value
            i += 1
        coil_dict["__solenoid_symmetry__"] = self.force_symmetry
        coil_dict["__optimisation_iteration__"] = self.iteration
        coil_dict["__optimisation_score__"] = self.score_function()
        coil_dict["__cell_length__"] = self.period
        json.dump(coil_dict, open(file_name, "w"), indent=2)

    def set_magnets(self):
        param_index = 0
        for i, coil in enumerate(self.coil_list):
            print("Coil"+str(i), end=" ")
            for key in self.fit_params:
                value = ROOT.Double()
                err = ROOT.Double()
                self.minuit.GetParameter(param_index, value, err)
                coil.__dict__[key] = value
                param_index += 1
                print(key, format(value, "6.4g"), end="; ")
            coil.reset()
            print()

    def get_test_field(self, z):
        test_field = 0
        for cell in range(-self.cell_extent, self.cell_extent+1):
            ztest0 = z-cell*self.period
            ztest1 = self.period-z-cell*self.period
            test_field += sum([coil.get_field(ztest0) for coil in self.coil_list])
            if self.force_symmetry == None:
                pass
            else:
                #z1 = self.period-z-cell*self.period
                test_field += self.force_symmetry*sum([coil.get_field(ztest1) for coil in self.coil_list])
        return test_field

    def compare_magnets(self):
        delta = 0.0 # RMS
        for i in range(self.n_fit_points):
            z = self.period*i/float(self.n_fit_points)
            ref_field = self.field_to_match.get_field(z)
            test_field = self.get_test_field(z)
            delta += (ref_field-test_field)**2
        return delta**0.5/self.n_fit_points

    def score_function(self, nvar=None, parameters=None, score=None, jacobian=None, err=None):
        self.iteration += 1
        self.set_magnets()
        if not score:
            score = [0.0]
        score[0] = self.compare_magnets()
        print("Iteration:", self.iteration, "Score:", score[0])
        return score[0]

    def plot_fit(self):
        z_list = [i*self.period/100 for i in range(101)]
        bref_list = [self.field_to_match.get_field(z) for z in z_list]
        bopt_list = [self.get_test_field(z) for z in z_list]
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.plot(z_list, bref_list, label="Reference Field")
        axes.plot(z_list, bopt_list, label="Generated Field")
        axes.legend()
        axes.set_xlabel("z [m]")
        axes.set_ylabel("B$_{z}$ [T]")
        figure.savefig(os.path.join(self.output_dir, "fitted_coils"+self.suffix+".png"))

def clear_dir(a_dir):
    try:
        shutil.rmtree(a_dir)
    except OSError:
        pass
    os.makedirs(a_dir)

def pixel_fit():
    output_dir = "coil_fitter_v4"
    clear_dir(output_dir)
    for b1, b3, l in [(6.0, 0.0, 0.8), (6.0, 3.0, 1.6), (12.0, 6.0, 0.8)]:
        field_to_match = SineField(0.0, b1, 0.0, b3, 0.0, 0.0, l)
        ri, dr, nr = 0.3, 0.1, 1 
        zi, dz, nz = l/16.0, l/16.0, 6
        fitter = CoilFitter.new_from_pixels(ri, dr, nr, zi, dz, nz, "_b1={0:.4g}_b3={1:.4g}_l={2:.4g}".format(b1, b3, l))
        fitter.output_dir = output_dir
        fitter.fit_coil(field_to_match, 1e-8)
        fitter.plot_fit()
    matplotlib.pyplot.show(block=False)

def trim_fit():
    pass

def main():
    pixel_fit()

if __name__ == "__main__":
    main()
    input("Press <CR> to finish")