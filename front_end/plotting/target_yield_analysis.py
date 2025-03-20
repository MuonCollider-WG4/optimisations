import json
import copy
import math
import operator
import matplotlib
import glob
import xboa
import g4bl_interface.g4bl_interface

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
