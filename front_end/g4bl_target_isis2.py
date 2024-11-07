import sys
import glob
import matplotlib
import copy
import operator
import os
import json
import xboa
import g4bl_interface.g4bl_interface
import g4bl_interface.stripper

class TargetRegion:
    def __init__(self):
        self.element_list = []
        self.target_length = 1
        self.target_radius = 15.0
        self.target_z_start = 50.0
        self.target_material = "C"
        self.max_z = 100.0
        self.solenoid_file = "share/target_solenoid.tex"

    def load_solenoid_from_latex(self):
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

    def build_target_dumb(self):
        target = {
            "name":"target",
            "outer_radius":self.target_radius,
            "z_position":self.target_z_start+self.target_length/2,
            "x_position":0.0,
            "length":self.target_length,
            "material":self.target_material,
            "type":"tube",
        }
        self.element_list.append(target)

    def build_target_sensitive(self):
        # WARNING: G4BL throws a segv if we have more than one detector
        elements = [{
            "name":"target",
            "outer_radius":self.target_radius,
            "z_position":self.target_z_start+self.target_length/2.,
            "length":self.target_length,
            "material":self.target_material,
            "type":"tube",
        },{
            "name":"target_outer",
            "inner_radius":self.target_radius+1e-3,
            "outer_radius":self.target_radius+2e-3,
            "z_position":None,
            "length":self.target_length,
            "material":"Vacuum",
            "type":"tube",
        },{
            "name":"target_outer_detector",
            "solid":"target_outer",
            "filename":"target_outer",
            "format":"ascii",
            "z_position":self.target_z_start+self.target_length/2.,
            "coordinates":"global",
            "type":"detector",
        },{
            "name":"target_end",
            "outer_radius":self.target_radius,
            "z_position":None,
            "length":1e-3,
            "material":"Vacuum",
            "type":"tube",
        },{
            "name":"target_end_detector",
            "solid":"target_end",
            "filename":"target_end",
            "format":"ascii",
            "z_position":self.target_z_start+self.target_length+1e-3,
            "coordinates":"global",
            "type":"detector",
        },{
            "name":"target_detector",
            "solid":"target",
            "filename":"target",
            "format":"ascii",
            "z_position":self.target_z_start+self.target_length/2.,
            "coordinates":"global",
            "type":"detector",
        }]
        self.element_list += elements[:3]


    def build_target(self, sensitive):
        if sensitive:
            self.build_target_sensitive()
        else:
            self.build_target_dumb()

    def build_beam_stop(self):
        beam_stop = {
            "name":"beam_stop",
            "outer_radius":1e3,
            "z_position":self.max_z+1,
            "x_position":0.0,
            "length":1.0,
            "material":"kill",
            "type":"tube",
        }

    def build_shielding(self):
        pass


def get_beam(my_linac, start_event, n_particles, energy, pid):
    mass = xboa.common.pdg_pid_to_mass[abs(pid)]
    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(my_linac.lattice_filename)[0],
        "x_position":0.0,
        "z_position":0.0,
        "start_event_number":start_event,
        "pid":pid,
        "beams":[{
            "type":"longitudinal_grid",
            "t_min":0.0, # ns
            "t_max":0.0, # ns
            "n_t_steps":1,
            "e_min":energy,
            "e_max":energy,
            "n_e_steps":n_particles,
            "pid":pid,
            "default_hit":{"x":0.0, "mass":mass, "pid":pid},
        },],
    }
    reference = {
        "p_start":((energy+mass)**2-mass**2)**0.5,
        "particle":"proton"
    }

    my_linac.beam = beam_def
    my_linac.reference = reference

class Analysis:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.z_analysis = 100.0
        self.linac = None

    def load_data(self, my_linac, target_region, beam_file_glob, plot_dir):
        self.linac = my_linac
        self.out_dir = my_linac.out_dir()
        self.plot_dir = plot_dir
        self.load_from_list(sorted(glob.glob(beam_file_glob)))
        for bunch in self.bunch_list:
            z_pos = bunch[0]["z"]
            if z_pos >= self.z_analysis:
                break
        self.bunch_analysis = bunch

    def load_from_list(self, filename_list):
        all_bunch_dict = {}
        self.bunch_list = []
        for filename in filename_list:
            bunch_dict = xboa.bunch.Bunch.new_dict_from_read_builtin("icool_for009", filename)
            for station, bunch in bunch_dict.items():
                if station in all_bunch_dict:
                    all_bunch_dict[station] += bunch.hits()
                else:
                    all_bunch_dict[station] = bunch.hits()
        for station, hit_list in sorted(all_bunch_dict.items()):
            self.bunch_list.append(xboa.bunch.Bunch.new_from_hits(hit_list))


    def do_plots(self, title, will_close):
        g4bl_interface.g4bl_interface.clean_dir(self.plot_dir, True)
        self.will_close = will_close
        #self.bunch_analysis_plot(title)
        #self.z_plot()
        self.z_energy(title, -13)

    def z_start(self, title, pid, e_range):
        event_set = set()
        track_start_list = []
        for bunch in self.bunch_list:
            for hit in bunch:
                if hit["pid"] != pid:
                    continue
                event_id = (hit["event_number"], hit["particle_number"])
                if event_id not in event_set:
                    event_set.add(event_id)
                    if hit["kinetic_energy"] < e_range[1] and hit["kinetic_energy"] > e_range[0]:
                        track_start_list.append(hit["z"])
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.hist(track_start_list, bins=100)
        axes.set_xlabel("z [mm]")
        axes.set_title(title+f"\nInitial {xboa.common.pdg_pid_to_name[pid]} position")
        if self.will_close:
            figure.savefig(f"{self.plot_dir}/z_start.png")
        json_data = {
            "energy":self.linac.beam["beams"][0]["e_min"],
            "n_particles":self.linac.beam["beams"][0]["n_e_steps"],
            "data":track_start_list
        }
        json.dump(json_data, open(f"{self.plot_dir}/z_start.json", "w"))

    def z_energy(self, title, pid):
        x_list = []
        y_list = []
        for bunch in self.bunch_list:
            for hit in bunch:
                if hit["pid"] != pid:
                    continue
                x_list.append(hit["z"])
                y_list.append(hit["kinetic_energy"])
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.hist2d(x_list, y_list, bins=[10, 100])
        axes.set_xlabel("z [mm]")
        axes.set_ylabel("kinetic energy [MeV]")
        axes.set_title(title+f"\nInitial {xboa.common.pdg_pid_to_name[pid]} position")
        figure.savefig(f"{self.plot_dir}/z_energy.png")


    def z_plot(self, title):
        z_list = [bunch[0]["z"] for bunch in self.bunch_list]
        bz_list = [bunch[0]["bz"]*1e3 for bunch in self.bunch_list]
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.plot(z_list, bz_list)
        axes.set_xlabel("z [mm]")
        axes.set_ylabel("B$_{z}$ [T]")
        axes.set_title(title)
        figure.savefig("z_plot.png")
        if self.will_close:
            matplotlib.pyplot.close(figure)

    def one_d_bunch_analysis_plot(self, x_label, var, bins, title):
        bunch = self.bunch_analysis
        bunch_in = self.bunch_list[0]
        bunch_weight = bunch.bunch_weight()
        figure = matplotlib.pyplot.figure()
        x_axes = figure.add_subplot(1, 1, 1)
        x_axes.set_xlabel(x_label)
        x_axes.set_xlim([min(bins), max(bins)])

        x_list_of_lists = []
        list_of_labels = []
        for pid in reversed(self.pid_list):
            bunch.clear_weights()
            p_name = "All"
            if pid != 0:
                bunch.cut({"pid":pid}, operator.ne)
                p_name = xboa.common.pdg_pid_to_name[pid]
            x_list_of_lists.append([hit[var] for hit in bunch.hits() if hit["weight"] > 0.5])
            list_of_labels.append(p_name+f" ({len(x_list_of_lists[-1])})")
        x_axes.hist(x_list_of_lists, bins=bins, label=list_of_labels, histtype="step")
        x_axes.legend()
        x_axes.set_title(title)
        bunch.clear_weights()
        figure.savefig(f"{self.plot_dir}/one_d_bunch_analysis_plot_{var}.png")
        if self.will_close:
            matplotlib.pyplot.close(figure)

    def two_d_bunch_analysis_plot(self, label_x, var_x, bins_x, label_y, var_y, bins_y, z_range, title):
        bunch = self.bunch_analysis
        bunch_in = self.bunch_list[0]
        bunch_weight = bunch.bunch_weight()
        for i, pid in enumerate(self.pid_list):
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)
            axes.set_xlabel(label_x)
            axes.set_xlim([min(bins_x), max(bins_x)])
            axes.set_ylabel(label_y)
            axes.set_ylim([min(bins_y), max(bins_y)])

            bunch.clear_weights()
            p_name = "All"
            if pid != 0:
                bunch.cut({"pid":pid}, operator.ne)
                p_name = xboa.common.pdg_pid_to_name[pid]
            x_list = [hit[var_x] for hit in bunch.hits() if hit["weight"] > 0.5]
            y_list = [hit[var_y] for hit in bunch.hits() if hit["weight"] > 0.5]
            axes.hist2d(x_list, y_list, bins=[bins_x, bins_y], vmin=z_range[0], vmax=z_range[1])
            axes.set_title(p_name+" "+title)
            bunch.clear_weights()
            figure.savefig(f"{self.plot_dir}/two_d_bunch_analysis_plot_{p_name}_{var_x}_vs_{var_y}.png")
            if self.will_close or i != 0:
                matplotlib.pyplot.close(figure)

    def bunch_analysis_plot(self, title):
        self.pid_list = self.linac.track_cuts.keep
        self.pid_list.remove(2212)
        for var, x_label, x_range in [("x'", "x' []", [-3+i*0.1 for i in range(61)]), ("kinetic_energy", "KE [MeV]", [i*1 for i in range(81)])]:
            self.one_d_bunch_analysis_plot(x_label, var, x_range, title)
        self.two_d_bunch_analysis_plot("x' []", "x'", [-3+i*0.2 for i in range(31)], "KE [MeV]", "kinetic_energy", [i*2 for i in range(41)], [0, 10], title)

def build_target(base_dir, start_event, n_particles, energy, cleanup_dir, will_gen_beam, pid, target_radius):
    target_region = TargetRegion()
    target_region.max_z = 1000.0
    target_region.target_z_start = 0.0001
    target_region.target_length = 1000.0 # 1.0*2/19.30 #
    target_region.target_radius = target_radius # 1.0*2/19.30 #
    target_region.target_material ="C"
    target_region.build_target(True)

    #print(json.dumps(target_region.element_list, indent=2))
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(os.path.join(base_dir, "target_region.g4bl"))
    my_linac.cleanup_dir = cleanup_dir
    my_linac.elements = target_region.element_list
    my_linac.physics.do_stochastics = 1 # e.g. decays
    #my_linac.physics.disable = "Decay" # e.g. decays
    #my_linac.physics.physics_list = "QGSP_BIC"
    my_linac.max_z = target_region.max_z+1
    my_linac.track_cuts = g4bl_interface.g4bl_interface.TrackCuts()
    my_linac.track_cuts.keep = [2212, 211, -211, 13, -13]
    my_linac.z_spacing = 1000.0
    if will_gen_beam:
        get_beam(my_linac, start_event, n_particles, energy, pid)
    return my_linac, target_region

def strip_file(run_dir, for009_file_list, g4bl_file_list):
    for009_stripper = g4bl_interface.stripper.Stripper.For009Stripper()
    g4bl_stripper = g4bl_interface.stripper.Stripper.G4BLStripper()
    for my_file in for009_file_list:
        my_file = f"{run_dir}/{my_file}"
        for009_stripper.strip(my_file)
    for my_file in g4bl_file_list:
        my_file = f"{run_dir}/{my_file}"
        g4bl_stripper.strip(my_file)

def dir_name(config, strip_glob):
    name = ""
    for key, value in config.items():
        if strip_glob and "?" in str(value) or "*" in str(value):
            continue
        name += f"{key}={value};_"
    return name[:-2]

def emu():
    mpi = xboa.common.pdg_pid_to_mass[211]
    mmu = xboa.common.pdg_pid_to_mass[13]
    pnu = (mpi**2-mmu**2)/2/mpi
    emu = (pnu**2+mmu**2)**0.5-mmu
    return emu

def main():
    do_execute = True
    do_analysis = True
    clear_dirs = False # recalculates field maps
    base_dir = "output/isis_target_model_v9/"
    for energy in [i*100 for i in range(3, 21)]:
        for target_radius in [10, 20]:
            n_versions = 10
            n_events = 1000000
            beam_pid = 2212
            for version in range(0, n_versions):
                param = {"energy":energy, "version":version, "radius":target_radius}
                run_dir = base_dir+dir_name(param, False)
                my_linac, target_region = build_target(run_dir, version*n_events, n_events, energy, clear_dirs, do_execute, beam_pid, target_radius)
                my_linac.build_linac()
                if do_execute: # alternatively we just redo the analysis
                    my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
                    my_execution.execute()
                    strip_file(run_dir, ["output_data.txt"], ["beam.txt", "target_outer.txt"])
            beam_file_glob = base_dir+dir_name(param, False)+"/output_data.txt"
            if do_analysis:
                title = dir_name(param, True).replace(";_", " ")
                beam_name = xboa.common.pdg_pid_to_name[beam_pid]
                my_analysis = Analysis(base_dir)
                plot_dir =  base_dir+dir_name(param, True)+"_plots"
                will_close = True #energy not in [400, 800, 1800]
                try:
                    my_analysis.load_data(my_linac, target_region, beam_file_glob, plot_dir)
                    my_analysis.do_plots(title, will_close)
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    print("Failed to plot")
        print("Finished")

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")
