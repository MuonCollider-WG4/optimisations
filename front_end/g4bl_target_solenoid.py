import matplotlib
import copy
import operator
import os
import json
import g4bl_interface.g4bl_interface
import xboa

class TargetRegion:
    def __init__(self):
        self.element_list = []
        self.target_length = 800.0
        self.target_radius = 15.0
        self.target_z_start = 0.0
        self.solenoid_file = "share/target_solenoid.tex"

    def load_solenoid_from_latex(self):
        self.element_list = []
        fin = open(self.solenoid_file)
        for i in range(3):
            fin.readline()
        for line in fin.readlines():
            words = line.split()
            if len(words) < 3:
                continue
            words = [word for word in words if word != "&"]
            m2mm = 1e3
            coil = {
                "name":words[0],
                "inner_radius":float(words[1])*m2mm,
                "outer_radius":(float(words[1])+float(words[3]))*m2mm,
                "z_position":float(words[2])*m2mm,
                "x_position":0.0,
                "length":float(words[4])*m2mm,
                "current":float(words[7]),
                "type":"solenoid",
            }
            nr = float(words[5])
            nz = float(words[6])
            coil["current"] *= nr*nz/(coil["outer_radius"]-coil["inner_radius"])/coil["length"] #A/mm^2
            self.element_list.append(coil)

    def build_target(self):
        target = {
            "name":"target",
            "outer_radius":self.target_radius,
            "z_position":self.target_z_start+self.target_length/2,
            "x_position":0.0,
            "length":self.target_length,
            "material":"C",
            "type":"tube",
        }
        self.element_list.append(target)

    def build_shielding(self):
        pass


def get_beam(my_linac, n_particles):
    pid = 2212
    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(my_linac.lattice_filename)[0],
        "x_position":0.0,
        "z_position":0.0,
        "pid":pid,
        "beams":[{
            "type":"longitudinal_grid",
            "t_min":0.0, # ns
            "t_max":0.0, # ns
            "n_t_steps":1,
            "e_min":5000,
            "e_max":5000,
            "n_e_steps":n_particles,
            "pid":pid,
            "default_hit":{"x":0.0, "mass":xboa.common.pdg_pid_to_mass[pid], "pid":pid},
        },][0:],
    }
    mass = xboa.common.pdg_pid_to_mass[pid]
    hit = xboa.hit.Hit.new_from_dict({"pid":pid, "energy":120.0+mass, "mass":mass}, "pz")
    reference = {
        "p_start":hit["pz"],
    }

    my_linac.beam = beam_def
    my_linac.reference = reference



def build_target(base_dir, n_particles, cleanup_dir):
    target_region = TargetRegion()
    target_region.load_solenoid_from_latex()
    target_region.build_target()
    #print(json.dumps(target_region.element_list, indent=2))
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(os.path.join(base_dir, "target_region.g4bl"))
    my_linac.cleanup_dir = cleanup_dir
    my_linac.elements = target_region.element_list
    my_linac.do_stochastics = 1 # e.g. decays
    my_linac.max_z = max([sol["z_position"] for sol in my_linac.elements])
    get_beam(my_linac, n_particles)
    return my_linac, target_region

class Analysis:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.z_analysis = 10000.0

    def load_data(self, my_linac, target_region, plot_dir):
        self.out_dir = my_linac.out_dir()
        self.plot_dir = plot_dir
        self.out_filename = os.path.join(self.out_dir, my_linac.output_file)+".txt"
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)
        for bunch in self.bunch_list:
            z_pos = bunch[0]["z"]
            if z_pos >= self.z_analysis:
                break
        self.bunch_analysis = bunch

    def do_plots(self):
        g4bl_interface.g4bl_interface.clean_dir(self.plot_dir, True)
        self.bunch_analysis_plot()
        self.z_plot()

    def z_plot(self):
        z_list = [bunch[0]["z"] for bunch in self.bunch_list]
        bz_list = [bunch[0]["bz"]*1e3 for bunch in self.bunch_list]
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.plot(z_list, bz_list)
        axes.set_xlabel("z [mm]")
        axes.set_ylabel("B$_{z}$ [T]")

    def one_d_bunch_analysis_plot(self, x_label, var):
        bunch = self.bunch_analysis
        bunch_weight = bunch.bunch_weight()
        figure = matplotlib.pyplot.figure()
        x_axes = figure.add_subplot(1, 1, 1)
        x_axes.set_xlabel(x_label)

        n_bins = int(bunch_weight**0.5)
        x_list_of_lists = []
        list_of_labels = []
        for pid in reversed(self.pid_list):
            bunch.clear_weights()
            p_name = "All"
            if pid != 0:
                bunch.cut({"pid":pid}, operator.ne)
                p_name = xboa.common.pdg_pid_to_name[pid]
            list_of_labels.append(p_name)
            x_list_of_lists.append([hit[var] for hit in bunch.hits() if hit["weight"] > 0.5])
        x_axes.hist(x_list_of_lists, bins=n_bins, label=list_of_labels, histtype="step")
        x_axes.legend()
        bunch.clear_weights()

    def bunch_analysis_plot(self):
        self.pid_list = [0, 2212, -13, 13, -11, 11, 211, -211]
        for var, x_label in [("x", "x [mm]"), ("kinetic_energy", "KE [MeV]")]:
            self.one_d_bunch_analysis_plot(x_label, var)

def main():
    do_execute = True
    clear_dirs = False # recalculates field maps
    base_dir = f"output/target_model_v2/"
    my_linac, target_region = build_target(base_dir, 10, clear_dirs)
    my_linac.build_linac()
    if do_execute: # alternatively we just redo the analysis
        my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
        my_execution.execute()
    my_analysis = Analysis(base_dir)
    my_analysis.load_data(my_linac, target_region, os.path.split(my_linac.lattice_filename)[0]+"/plots")
    my_analysis.do_plots()
    print("Finished")

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")
