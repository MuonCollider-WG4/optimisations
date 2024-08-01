import math
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
import numpy

class TargetRegion:
    def __init__(self):
        self.element_list = []
        self.target_length = 800.0
        self.target_radius = 30.0
        self.target_z_start = 50.0
        self.beam_sigma = 5.0
        self.target_material = "C"
        self.beam_pipe_radius = 500.0
        self.max_z = 50000.0
        self.solenoid_file = "share/target_solenoid.tex"
        # tanh solenoid:
        self.end_length = 1000.0 # mm
        self.peak_field = 20 # T
        self.low_field = 1.5 # T
        self.b_target_end = 19.0 # field at target end (defines the z position of the tanh field)
        self.n_centre = 10 # centre_length = end_length*self.n_centre


    def build_solenoid_from_latex(self):
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


    def build_tanh_solenoid(self):
        centre_length = self.end_length*self.n_centre # make sure it is body-dominated
        field_scale = (self.b_target_end-self.low_field)/(self.peak_field-self.low_field)
        z_target_end = self.target_length+self.target_z_start
        dx = get_tanh_offset(self.end_length, centre_length, z_target_end, field_scale)
        derivatives_solenoid = {
            "name":"target_solenoid",
            "type":"derivatives_solenoid",
            "z_position":-dx,
            "length":20*self.end_length,
            "nominal_field":self.peak_field-self.low_field,
            "end_length":self.end_length,
            "centre_length":centre_length,
            "field_model":"tanh",
        }
        uniform_field = {
            "name":"uniform_field",
            "type":"uniform_field",
            "length":1000000.0, # 1000 m
            "radius":self.beam_pipe_radius, # 1 m
            "bz":self.low_field,
        }

        self.element_list.append(derivatives_solenoid)
        self.element_list.append(uniform_field)

    def build_target(self):
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

    def build_beam_stop(self):
        beam_stop = {
            "name":"beam_stop",
            "outer_radius":self.beam_pipe_radius,
            "z_position":self.max_z+1,
            "x_position":0.0,
            "length":1.0,
            "material":"kill",
            "type":"tube",
        }
        beam_start = copy.deepcopy(beam_stop)
        beam_start["z_position"] = -1.0
        beam_start["name"] = "beam_start"
        self.element_list.append(beam_start)
        self.element_list.append(beam_stop)


    def build_beam_pipe(self):
        beam_pipe = {
            "name":"beam_pipe",
            "inner_radius":self.beam_pipe_radius+0.1,
            "outer_radius":2*self.beam_pipe_radius+0.1,
            "z_position":self.max_z/2,
            "x_position":0.0,
            "length":self.max_z*2,
            "material":"kill",
            "type":"tube",
        }
        self.element_list.append(beam_pipe)

def get_beam(my_linac, start_event, n_particles, kinetic_energy):
    pid = 2212
    mass = xboa.common.pdg_pid_to_mass[pid]
    p = ((kinetic_energy+mass)**2-mass**2)**0.5
    px = 1e-3*p
    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(my_linac.lattice_filename)[0],
        "x_position":0.0,
        "z_position":0.0,
        "start_event_number":start_event,
        "pid":pid,
        "beams":[{
            "type":"beam_ellipse",
            "variables":["x", "px", "y", "py"], # ns
            "ellipse":[[25.0, 0.0, 0.0, 0.0], [0.0, px**2, 0.0, 0.0], [0.0, 0.0, 25.0, 0.0], [0.0, 0.0, 0.0, px**2]], # ns
            "mean":[0.0, 0.0, 0.0, 0.0],
            "n_particles":n_particles,
            "mass_shell_condition":"pz",
            "default_hit":{"energy":kinetic_energy+mass, "mass":xboa.common.pdg_pid_to_mass[pid], "pid":pid},
        },],
    }
    reference = {
        "p_start":((kinetic_energy+mass)**2-mass**2)**0.5,
        "particle":"proton"
    }

    my_linac.beam = beam_def
    my_linac.reference = reference

class Analysis:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.z_analysis = 100.0
        self.energy_range = [10, 500] # kinetic
        self.virtual_pipe_radius = 500
        self.virtual_bz = 1.5
        self.linac = None
        self.out_dir = None
        self.plot_dir = None
        self.bunch_analysis = None
        self.pid_list = []
        self.output_plot_data = {
            "cut_scan":[],
        }

    def load_data(self, my_linac, target_region, beam_file_glob, plot_dir):
        self.linac = my_linac
        self.virtual_bz = target_region.low_field
        self.out_dir = my_linac.out_dir()
        self.plot_dir = plot_dir
        self.load_from_list(sorted(glob.glob(beam_file_glob)))
        for bunch in self.bunch_list:
            z_pos = bunch[0]["z"]
            if z_pos >= self.z_analysis:
                break
        self.bunch_analysis = bunch
        self.pid_list = self.linac.track_cuts.keep
        self.pid_list.remove(2212)

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

    def do_cut(self, pid):
        self.bunch_analysis.clear_weights()
        self.bunch_analysis.cut({"kinetic_energy":self.energy_range[0]}, operator.lt)
        self.bunch_analysis.cut({"kinetic_energy":self.energy_range[1]}, operator.gt)
        if pid != 0:
            self.bunch_analysis.cut({"pid":pid}, operator.ne)
        means = {"x":0.0, "px":0.0, "y":0.0, "py":0.0}
        charge = xboa.common.pdg_pid_to_charge[pid] # particle charge sign??
        for hit in self.bunch_analysis:
            if hit["local_weight"] == 0.0:
                continue
            max_r2 = self.get_max_r2(hit, self.virtual_bz)
            print_color = ""
            if max_r2 > self.virtual_pipe_radius**2:
                hit["local_weight"] = 0.0
                #print_color = "\033[92m"
                #print(print_color+f"r cut, {hit['x']:8.4g} {hit['y']:8.4g} {hit['px']:8.4g} {hit['py']:8.4g} r: {hit['r']:8.4g} {max_r2**0.5:8.4g} {self.virtual_pipe_radius:8.4g}")#+"\033[0m")

    def get_max_r2(self, hit, bz):
        kilotesla = 1e3
        c_light = xboa.common.constants['c_light']
        q = hit["charge"]
        # coordinates of helix centre
        phi0 = math.atan2(hit["px"], hit["py"])
        pt = hit["pt"]
        if pt == 0.0: # probably the reference particle
            return 0.0
        x0 = hit["x"] + pt*math.cos(phi0)/q/bz/c_light*kilotesla
        y0 = hit["y"] - pt*math.sin(phi0)/q/bz/c_light*kilotesla
        r0 = abs(pt/q/bz/c_light*kilotesla) # wavenumber is qBm/pz dz
        max_r2 = ((x0**2 + y0**2)**0.5 + r0)**2 # vector to centre + radius of track
        return max_r2

    def do_plots(self, title, will_close):
        g4bl_interface.g4bl_interface.clean_dir(self.plot_dir, True)
        self.will_close = will_close
        self.bunch_analysis_plot(title)
        self.cut_scan_plots(title)
        self.print_beam(f"{self.plot_dir}/beam_out.txt", "g4beamline_bl_track_file", title)
        self.write_output_data(f"{self.plot_dir}/data_out.json")

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
                self.do_cut(pid)
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

    def one_cut_scan_plot(self, x_list, x_label, scan_lambda, pid, title, series_label, figure = None):
        """
        Apply scan_lambda(x), then self.do_cut(pid) and plot the resultant bunch_weight
        """
        y_list = []
        for x in x_list:
            scan_lambda(x)
            self.do_cut(pid)
            y_list.append(self.bunch_analysis.bunch_weight())
        if figure == None:
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)
            axes.set_title(title)
            axes.set_xlabel(x_label)
            axes.set_ylabel("Transmission")
        else:
            axes = figure.axes[0]
        axes.plot(x_list, y_list, label=series_label)

        self.output_plot_data["cut_scan"].append({
                "energy_range":self.energy_range,
                "pid":pid,
                "x_label":x_label,
                "x_data":x_list,
                "y_data":y_list
            })

        return figure

    def print_beam(self, file_name, file_format, title):
        bunch = self.bunch_analysis.deepcopy()
        bunch.clear_weights()
        bunch.clear_global_weights()
        title += f" z: {bunch[0]['z']} [m]"
        for i, hit in enumerate(bunch):
            hit["z"] = 0.0
        bunch = xboa.bunch.Bunch.new_from_hits([hit for hit in bunch if hit["pid"] in self.pid_list])
        bunch.hit_write_builtin(file_format, file_name, user_comment=title)

    def write_output_data(self, file_name):
        fin = open(file_name, "w", encoding="utf-8")
        output_str = json.dumps(self.output_plot_data, indent=2)
        fin.write(output_str)
        fin.close()

    def set_pipe_radius(self, r):
        self.virtual_pipe_radius = r

    def cut_scan_plots(self, title):
        print("Doing cut scan plots")
        energy_range = copy.deepcopy(self.energy_range)
        pipe_radius = copy.deepcopy(self.virtual_pipe_radius)
        label = f"{self.energy_range[0]} < E < {self.energy_range[1]} [MeV]"
        figure = self.one_cut_scan_plot([50+i*10 for i in range(51)], "Virtual pipe radius [mm]", self.set_pipe_radius, -13, title, label, figure = None)
        self.energy_range = [40, 180]
        label = f"{self.energy_range[0]} < E < {self.energy_range[1]} [MeV]"
        figure = self.one_cut_scan_plot([50+i*10 for i in range(51)], "Virtual pipe radius [mm]", self.set_pipe_radius, -13, title, label, figure)
        self.energy_range = [40, 500]
        label = f"{self.energy_range[0]} < E < {self.energy_range[1]} [MeV]"
        figure = self.one_cut_scan_plot([50+i*10 for i in range(51)], "Virtual pipe radius [mm]", self.set_pipe_radius, -13, title, label, figure)
        figure.axes[0].legend()
        figure.savefig(f"{self.plot_dir}/cut_scan_plot.png")
        if self.will_close:
            matplotlib.pyplot.close(figure)
        self.energy_range = energy_range
        self.virtual_pipe_radius = pipe_radius
        print("     ... Done")

    def two_d_bunch_analysis_plot(self, label_x, var_x, bins_x, label_y, var_y, bins_y, z_range, title):
        bunch = self.bunch_analysis
        bunch_in = self.bunch_list[0]
        bunch_weight = bunch.bunch_weight()
        lambda_dict = {
            "max_r2":(lambda hit: self.get_max_r2(hit, self.virtual_bz)),
            "max_r":(lambda hit: self.get_max_r2(hit, self.virtual_bz)**0.5),
        }
        for i, pid in enumerate(self.pid_list):
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)
            axes.set_xlabel(label_x)
            axes.set_xlim([min(bins_x), max(bins_x)])
            axes.set_ylabel(label_y)
            axes.set_ylim([min(bins_y), max(bins_y)])

            self.do_cut(pid)
            p_name = "All"
            if pid != 0:
                p_name = xboa.common.pdg_pid_to_name[pid]
            if var_x in lambda_dict:
                a_lambda = lambda_dict[var_x]
                x_list = [a_lambda(hit) for hit in bunch.hits() if hit["weight"] > 0.5]
            else:
                x_list = [hit[var_x] for hit in bunch.hits() if hit["weight"] > 0.5]
            if var_y in lambda_dict:
                a_lambda = lambda_dict[var_y]
                y_list = [a_lambda(hit) for hit in bunch.hits() if hit["weight"] > 0.5]
            else:
                y_list = [hit[var_y] for hit in bunch.hits() if hit["weight"] > 0.5]

            h2d = axes.hist2d(x_list, y_list, bins=[bins_x, bins_y], vmin=z_range[0], vmax=z_range[1])
            figure.colorbar(h2d[3])
            axes.set_title(title+" "+p_name)
            figure.savefig(f"{self.plot_dir}/two_d_bunch_analysis_plot_{p_name}_{var_x}_vs_{var_y}.png")
            if self.will_close or i != 0:
                matplotlib.pyplot.close(figure)

    def bunch_analysis_plot(self, title):
        for var, x_label, x_range in [("x'", "x' []", [-0.5+i*0.01 for i in range(101)]), ("kinetic_energy", "KE [MeV]", [i*1 for i in range(81)]), ("kinetic_energy", "KE [MeV]", [i*10 for i in range(81)])]:
            self.one_d_bunch_analysis_plot(x_label, var, x_range, title)
        self.two_d_bunch_analysis_plot("x' []", "x'", [-0.5+i*0.01 for i in range(101)], "KE [MeV]", "kinetic_energy", [i*10 for i in range(61)], [None, None], title)
        self.two_d_bunch_analysis_plot("x [mm]", "x", [i*10 for i in range(-40, 41)], "x' []", "x'", [-0.5+i*0.01 for i in range(101)], [None, None], title)
        self.two_d_bunch_analysis_plot("x [mm]", "x", [i*10 for i in range(-40, 41)], "y [mm]", "y", [i*10 for i in range(-40, 41)], [0, 1], title)
        self.two_d_bunch_analysis_plot("r max [mm]", "max_r", [i*10 for i in range(50)], "kinetic energy [MeV]", "kinetic_energy", [i*10 for i in range(61)], [None, None], title)

def bz_plot(file_list, plot_file):
    figure = matplotlib.pyplot.figure()
    axes = figure.add_subplot(1, 1, 1)
    for a_file in sorted(file_list):
        end_length = a_file.split("/")[-2]
        end_length = end_length.split("=")[-1]
        fin = open(a_file)
        for i in range(2):
            fin.readline() # header
        z_list = []
        bz_list = []
        for line in fin.readlines():
            [x, y, z, t, bx, by, bz, ex, ey, ez] = line.split()
            if float(x) != 0.0 or float(y) != 0.0:
                continue
            z_list.append(float(z))
            bz_list.append(float(bz))
        axes.plot(z_list, bz_list, label="End Length: "+end_length+" mm")
    axes.legend()
    axes.set_xlabel("z [mm]")
    axes.set_ylabel("B$_{z}$ [T]")
    axes.set_xlim([0, 10000])
    figure.savefig(plot_file)

def build_target(base_dir, start_event, n_particles, energy, end_field_length, target_magnitude, cw_magnitude, cleanup_dir, will_gen_beam):
    target_region = TargetRegion()
    target_region.max_z = 50000.0
    target_region.z_start = 50.0
    target_region.target_length = 800.0 # 1.0*2/19.30 #
    target_region.target_material ="C" # "W" #
    target_region.build_target()

    target_region.end_length = end_field_length # mm
    target_region.peak_field = target_magnitude # T
    target_region.low_field = cw_magnitude # T
    target_region.build_tanh_solenoid()
    #target_region.build_solenoid_from_latex()

    #print(json.dumps(target_region.element_list, indent=2))
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(os.path.join(base_dir, "target_region.g4bl"))
    my_linac.cleanup_dir = cleanup_dir
    my_linac.elements = target_region.element_list
    my_linac.do_stochastics = 1 # e.g. decays
    my_linac.max_z = target_region.max_z+1
    my_linac.track_cuts = g4bl_interface.g4bl_interface.TrackCuts()
    my_linac.track_cuts.keep = [2212, -13, 13, 211, -211]
    my_linac.z_spacing = 10000.0
    if will_gen_beam:
        get_beam(my_linac, start_event, n_particles, energy)
    return my_linac, target_region

def strip_file(run_dir):
    my_file = f"{run_dir}/output_data.txt"
    stripper = g4bl_interface.stripper.Stripper()
    stripper.strip(my_file)
    beam_file = f"{run_dir}/beam.txt"
    os.unlink(beam_file)

def dir_name(config):
    name = ""
    for key, value in config.items():
        name += f"{key}={value};_"
    return name[:-2]

def main():
    do_execute = False
    do_analysis = True
    clear_dirs = False # recalculates field maps
    base_dir = "output/mucol_target_v2/"
    target_field = 20
    uniform_field = 1.5
    for energy in [5000]:#, 10000]:
        for end_length in [1000*i for i in range(1, 3)]:
            n_versions = 100
            n_events = 1000
            for version in range(0, n_versions):
                run_dir = base_dir+dir_name({"energy":energy, "version":version, "end_length":end_length})
                my_linac, target_region = build_target(run_dir, version*n_events, n_events, energy, end_length, target_field, uniform_field, clear_dirs, do_execute)
                my_linac.build_linac()
                if do_execute: # else we just redo the analysis
                    my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
                    my_execution.execute()
                    #strip_file(run_dir)
            beam_file_glob = base_dir+dir_name({"energy":energy, "version":"*", "end_length":end_length})+"/output_data.txt"
            if do_analysis:
                target = f"{target_region.target_length:6.3g} mm {target_region.target_material}"
                title = f"energy: {energy} [MeV]   n$_{{protons}}$: {float(n_versions*n_events):4.6g} end length: {end_length/1000} [m]"
                my_analysis = Analysis(base_dir)
                my_analysis.z_analysis = target_region.max_z
                plot_dir =  base_dir+dir_name({"energy":energy, "end_length":end_length})+"_plots"
                will_close = end_length != 1000 and end_length != 9000
                try:
                    my_analysis.load_data(my_linac, target_region, beam_file_glob, plot_dir)
                    my_analysis.do_plots(title, will_close)
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    print("Failed to plot")
    print("Finished")
    if do_analysis:
        bz_plot(glob.glob(base_dir+"/energy=5000*version=0;*/fieldmap.dat.txt"), base_dir+"bz.png")
        matplotlib.pyplot.show(block=False)
        input("Press <CR> to finish")

def get_tanh_offset(end_length, centre_length, z_pos, field_scale):
    """
    Given end length, what is the off set of the derivatives solenoid that
    yields a field B = field_scale*B_max at position (z_pos)
    """
    x0 = centre_length/2.0
    dx = x0-z_pos-numpy.arctanh(2*field_scale-1)*end_length
    b = (numpy.tanh((x0-dx-z_pos)/end_length)+1)/2
    print(f"Tanh offset calculation: end_length {end_length} centre_length {centre_length} dx {dx:6.4g} b(dx): {b:6.4g}")
    return dx

def get_tanh(end_length, centre_length, b0):
    for i in range(41):
        z = 100*i
        print(z, b0*(numpy.tanh((centre_length/2-z)/end_length)+1)/2)

if __name__ == "__main__":
    main()
