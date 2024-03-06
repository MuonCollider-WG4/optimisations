import ctypes
import os
import shutil
import json
import matplotlib
import matplotlib.pyplot
import ROOT
from models.field_models import SineField
from models.field_models import CurrentSheet
from models.field_models import CurrentBlock

class CoilFitter(object):
    """
    Routines to fit a bunch of current sheets (aka coils) to some arbitrary
    on-axis field
    """
    def __init__(self, ri = 0.2):
        """
        Initialisation. 
        """
        self.n_iterations = 10000
        self.n_sheets_per_coil = 1
        self.cell_extent = 10 # n in each direction beyond the first
        self.n_fit_points = 100
        self.period = -1
        self.fit_params = ["b0"]#, "zcentre", "length"] # value:seed
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

    @classmethod
    def new_from_coil(self, ri, dr, zc, dz, suffix):
        """
        Initialise a new CoilFitter based on a single coil
        - ri: inner radius of the coil [m]
        - dr: radial thickness [m]
        - zc: position of the coil z centre
        - dz: length of the coil [m]
        - suffix: string
        """
        # b0, zcentre, length, rmin, rmax, period, nrepeats, nsheets
        fitter = CoilFitter()
        fitter.coil_list = []
        fitter.coil_list.append(
            CurrentBlock(1.0, zc, dz, ri, ri+dr, 2.0, 1, 10)
        )
        if len(suffix) and suffix[0] != "_":
            suffix = "_"+suffix
        fitter.suffix = suffix
        return fitter

    def get_plot_value(self, var):
        value = None
        coil = None
        if 'coil' in var:
            coil_number = var['coil']
            coil = self.coil_list[coil_number]
        if var['parameter'] == "sigma":
            value = self.score_function()
        elif var['parameter'] == 'current_density':
            value = abs(coil.get_current_density())
        elif var['parameter'] == 'zmax':
            value = coil.zcentre+coil.length/2
        elif var['parameter'] == 'zmin':
            value = coil.zcentre-coil.length/2
        else:
            value = coil.__dict__[var['parameter']]
        if 'minimum' in var and value < var['minimum']:
            value = var['minimum']
        if 'maximum' in var and value > var['maximum']:
            value = var['maximum']
        return value

    def get_scan_values(self, var):
        if var['n_steps'] == 1:
            value_list = [(var['maximum']+var['minimum'])/2.0]
            bin_list = [var['minimum'], var['maximum']]
            return value_list, bin_list

        delta = (var['maximum']-var['minimum'])/(var['n_steps']-1)

        value_list = [delta*index+var['minimum'] for index in range(var['n_steps'])]
        bin_list = [value-delta*0.5 for value in value_list]
        bin_list.append(bin_list[-1]+delta)
        return value_list, bin_list

    def apply_value(self, var, value):
        param_coil_list = self.coil_list
        if var['coil'] is not None:
            param_coil_list = [self.coil_list[var['coil']]]
        for coil in param_coil_list:
            coil.__dict__[var['parameter']] = value
        return value

    def get_bin_edges(self, var):
        if var['n_steps'] == 1:
            return [var['minimum'], var['maximum']]
        else:
            delta = (var['maximum']-var['minimum'])/(var['n_steps']-1)
            return [var['minimum']+delta*i-0.5 for i in range(var['n_steps']+1)]

    def plot_scan_2d(self, var_x, var_y, var_z_list):
        """
        Scan over variable var_x and var_y, plot the fit function at each point
        - var_x, var_y are dictionaries, which defines the x variable and the y
        variable and have the following string parameters:
            - var['coil']: integer that indexes the coil that var pertains to.
            If None, apply to all coils.
            - var['parameter']: string name that specifies the coil member data
            that var pertains to
            - var['minimum']: minimum value for the variable
            - var['maximum']: maximum value for the variable
            - var['n_steps']: number of steps (>0) that will be made in the scan,
            including the minimum and maximum. If 1, then the variable will be
            stepped only at the midpoint between minimum and maximum
        - var_z: a dictionary that describes the z (i.e. plot) variable and has
        the following string parameters:
            - var['parameter']: string name that specifies either coil member
            data or 'sigma' for the RMS deviation from required field or
            'current_density' for a coil current density
            - var['coil']: integer that indexes the coil that the variable
            pertains to.
            - var['minimum']: if the value is less than minimum, it will be set
            to minimum instead. Ignored if 'minimum' is None or not defined.
            - var['maximum']: if the value is more than maximum, it will be set
            to maximum instead. Ignored if 'maximum' is None or not defined.
        """
        x_list = []
        y_list = []
        z_list = [[] for var_z in var_z_list]
        x_values, x_bins = self.get_scan_values(var_x)
        y_values, y_bins = self.get_scan_values(var_y)
        for x in x_values:
            for y in y_values:
                x_list.append(self.apply_value(var_x, x))
                y_list.append(self.apply_value(var_y, y))
                self.fit_coil(1e-4)
                for iz, var_z in enumerate(var_z_list):
                    z_list[iz].append(self.get_plot_value(var_z))
                self.save_coil_params(os.path.join(self.output_dir, "coil_parameter"), "a", None)
        for iz, var_z in enumerate(var_z_list):
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)
            axes.set_xlabel(self.human_readable(var_x['parameter']))
            axes.set_ylabel(self.human_readable(var_y['parameter']))
            h = axes.hist2d(x_list, y_list, [x_bins, y_bins], weights=z_list[iz])
            figure.colorbar(h[3])
            axes.set_title(self.human_readable(var_z['parameter']))
            fout = f"{var_x['parameter']}_vs_{var_y['parameter']}_vs_{var_z['parameter']}.png"
            figure.savefig(os.path.join(self.output_dir, fout))

    def fit_coil(self, tolerance):
        self.iteration = 0
        self.period = self.field_to_match.period
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
        global my_self
        my_self = self
        self.minuit.SetFCN(score_function)
        self.minuit.Command("SIMPLEX "+str(self.n_iterations)+" "+str(tolerance))
        self.print_coil_params()
        self.save_coil_params(os.path.join(self.output_dir, "coil_params_beta-scale"+self.suffix+".out"))
        my_self = None

    def human_readable(self, key):
        my_var = {
            "sigma":"RMS($\\delta B$) [T]",
            "b0":"nominal field [T]",
            "zcentre":"centre position [m]",
            "length":"length [m]",
            "rmin":"inner radius [m]",
            "rmax":"outer radius [m]",
            "zmax":"maximum coil z [m]",
            "zmin":"minimum coil z [m]",
            "current_density":"current density [A/m^2]",
        }[key]
        return my_var

    def print_coil_params(self):
        for i, coil in enumerate(self.coil_list):
            print("Coil parameters for coil "+str(i))
            print("  nominal field [T]:       ", format(coil.b0, "9.4g"))
            print("  current density [A/mm^2]:", format(coil.get_current_density()*1e-6, "9.4g"))
            print("  centre position [m]:    ", format(coil.zcentre, "9.4g"))
            print("  length [m]:             ", format(coil.length, "9.4g"))
            print("  inner radius [m]:       ", format(coil.rmin, "9.4g"))
            print("  outer radius [m]:       ", format(coil.rmax, "9.4g"))
            print("  hoop stress:            ", "NEEDS CALC")

    def save_coil_params(self, file_name, file_mode="w", indent=2):
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
        with open(file_name, file_mode) as fout:
            json.dump(coil_dict, fout, indent=indent)
            fout.write("\n")


    def set_magnets(self):
        param_index = 0
        for i, coil in enumerate(self.coil_list):
            print("Coil"+str(i), end=" ")
            for key in self.fit_params:
                value = ctypes.c_double()
                err = ctypes.c_double()
                self.minuit.GetParameter(param_index, value, err)
                py_value = float(value.value)
                coil.__dict__[key] = py_value
                param_index += 1
                print(key, format(py_value, "6.4g"), end="; ")
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

    def physicality_penalty(self):
        penalty = 0.0
        for i, coil in enumerate(self.coil_list):
            radial_oops = coil.rmin - coil.rmax
            if radial_oops > 0:
                print(f"coil {i} radial underflow {radial_oops} [m] r0 {coil.rmin} r1 {coil.rmax}")
                penalty += radial_oops
            delta_length = coil.zcentre + coil.length/2 - coil.period/4
            #print("delta length ", delta_length, coil.zcentre, coil.length, coil.period)
            if delta_length > 0:
                print(f"coil {i} longitudinal overflow {delta_length} [m] zc {coil.zcentre} l {coil.length}")
                penalty += delta_length
            delta_length = coil.length/2-coil.zcentre
            if delta_length > 0:
                print(f"coil {i} longitudinal underflow {delta_length} [m] zc {coil.zcentre} l {coil.length}")
                penalty += delta_length
        return penalty

    def score_function(self, nvar=None, parameters=None, score=None, jacobian=None, err=None):
        self.iteration += 1
        self.set_magnets()
        if score is None:
            score = ctypes.c_double(0.0)
        score.value = self.compare_magnets()
        penalty = self.physicality_penalty() # check if the magnet is physical
        score.value += penalty
        score.value *= 1+penalty
        print("Iteration:", self.iteration, "Score:", score.value)
        return score.value

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

my_self = None
def score_function(nvar=None, parameters=None, score=None, jacobian=None, err=None):
    return my_self.score_function(nvar, parameters, score, jacobian, err)

def clear_dir(a_dir):
    try:
        shutil.rmtree(a_dir)
    except OSError:
        pass
    os.makedirs(a_dir)

def pixel_fit():
    output_dir = "output/demo_field_v1"
    clear_dir(output_dir)
    for scale in [1.0]: #0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2.0]:
        b1 = 7.0
        b2 = 1.0
        l = 1.0
        ri, dr, nr = 0.5, 0.1, 1
        #zcentre = zi+dz*(zindex+0.5)
        dz = l/10.0
        zi, dz, nz = l/4-dz/2, l/10.0, 6
        fitter = CoilFitter.new_from_pixels(ri, dr, nr, zi, dz, nz, "_b1={0:.4g}_b2={1:.4g}_l={2:.4g}".format(b1, b2, l))
        fitter.output_dir = output_dir
        fitter.fit_params = ["b0"] #, "length"] # value:seed
        fitter.limits = [[0, 0]]
        fitter.field_to_match = SineField(0.0, b1, b2, 0.0, 0.0, 0.0, l)
        fitter.fit_coil(1e-8)
        fitter.plot_fit()
    matplotlib.pyplot.show(block=False)

def get_demo_old():
    field_to_match = CurrentBlock(1.0, 0.1, 0.1, 0.4, 0.5, 1.0, 20, 10)
    #bnom = 3003.451065/(field_to_match.get_current_density()*1e-6)
    #field_to_match.b0 = bnom
    #field_to_match.reset()
    #print("Current density", field_to_match.get_current_density()*1e-6, "[A/mm^2]")
    #print("n repeats", field_to_match.nrepeats)
    #print("period", field_to_match.period)
    return field_to_match

def scan():
    output_dir = "output/demo_field_v5"
    clear_dir(output_dir)
    fitter = CoilFitter.new_from_coil(0.2, 0.1, 0.1, 0.1, "test")
    fitter.field_to_match = SineField(0.0, 7.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    #fitter.field_to_match = get_demo_old()
    fitter.output_dir = output_dir
    var_x = {
        'coil': 0,
        'parameter':'length',
        'minimum':0.05,
        'maximum':0.4,
        'n_steps':25,
    }
    var_y = {
        'coil': 0,
        'parameter':'rmin',
        'minimum':0.05,
        'maximum':0.5,
        'n_steps':25,
    }
    var_x = {
        'coil': 0,
        'parameter':'length',
        'minimum':0.10,
        'maximum':0.20,
        'n_steps':1,
    }
    var_y = {
        'coil': 0,
        'parameter':'rmin',
        'minimum':0.25,
        'maximum':0.35,
        'n_steps':1,
    }
    fitter.fit_params = ["b0", "zcentre", "rmax"] #, "length"] # value:seed
    fitter.limits = [[0, 0], [0.05, 0.45], [0.0, 0.8]]
    plot_variables = [
        {"parameter":"sigma"},
        {"parameter":"current_density", "coil":0, "minimum":-1000e6, "maximum":1000e6},
        {"parameter":"zcentre", "coil":0},
        {"parameter":"rmax", "coil":0},
        {"parameter":"zmax", "coil":0},
        {"parameter":"zmin", "coil":0},
    ]

    fitter.plot_scan_2d(var_x, var_y, plot_variables)
    fitter.plot_fit()

def main():
    scan()

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")