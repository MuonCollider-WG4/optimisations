import copy
import os
import matplotlib
import g4bl_interface.g4bl_interface
import xboa.hit
import xboa.bunch

def multiharmonic_rf(n_cells, principle_harmonic, v1, v2, v3, v4, phi1, phi2, phi3, phi4):
    cavities = [{
            "name":f"pillbox_{i}.1",
            "inner_length":200.0,
            "frequency":principle_harmonic*1,
            "max_gradient":v1,
            "phase":phi1,
            "rgb":[0.5,1,0],
            "z_position":150+1000.0*i,
            "type":"cavity",
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.2",
            "inner_length":200.0,
            "frequency":principle_harmonic*2,
            "max_gradient":v2,
            "phase":phi2,
            "rgb":[1,1,0],
            "z_position":400+1000.0*i,
            "type":"cavity",
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.3",
            "inner_length":200.0,
            "frequency":principle_harmonic*3,
            "max_gradient":v3,
            "phase":phi3,
            "rgb":[1,0.5,0],
            "z_position":650+1000.0*i,
            "type":"cavity",
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.4",
            "inner_length":200.0,
            "frequency":principle_harmonic*4,
            "max_gradient":v4,
            "phase":phi4,
            "rgb":[0.5,0.5,0],
            "z_position":900.0+1000.0*i,
            "type":"cavity",
        } for i in range(n_cells)
    ]
    return cavities


def build_final_cooling_lattice(cleanup = True):
    lattice_filename = "output/multiharmonic_v6/linac.g4bl"
    # f=0.02, v1=0.9, v2=0.0
    frequency = 0.020
    voltage_1 = 1.7
    voltage_2 = -0.0
    phi_1 = 25.0
    tune_coeff = 1*int(6/(frequency*voltage_1)**0.5)
    print("Simulating", tune_coeff, "cells")
    if tune_coeff < 1:
        tune_coeff = 1
    cavities = multiharmonic_rf(tune_coeff, frequency,
                                voltage_1, voltage_2, 0.0, 0.0,
                                phi_1, 0.0, 0.0, 0.0)
    mass = xboa.common.pdg_pid_to_mass[13]
    pz_in = 70.2 # table from elena
    sigma_e = 5.5 # is this full width or half width (or something else)
    sigma_z = 1400
    energy_in = (pz_in**2+mass**2)**0.5-mass # 21.2 -> 25.7 MeV
    speed_in = pz_in/(energy_in+mass)*xboa.common.constants["c_light"]
    sigma_t = sigma_z/speed_in

    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(lattice_filename)[0],
        "beams":[{
            "type":"longitudinal_grid",
            "t_min":-sigma_t, # ns
            "t_max":+sigma_t, # ns
            "n_t_steps":1,
            "e_min":energy_in-sigma_e,
            "e_max":energy_in+sigma_e,
            "n_e_steps":11,
            "pid":-13,
            "default_hit":{"mass":xboa.common.pdg_pid_to_mass[13]}
        }, {
            "type":"longitudinal_grid",
            "t_min":-sigma_t, # ns
            "t_max":+sigma_t, # ns
            "n_t_steps":11,
            "e_min":energy_in-sigma_e,
            "e_max":energy_in+sigma_e,
            "n_e_steps":1,
            "pid":-13,
            "default_hit":{"mass":xboa.common.pdg_pid_to_mass[13]}
        }, {
            "type":"longitudinal_grid",
            "t_min":-1/frequency/2, # ns
            "t_max":+1/frequency/2, # ns
            "n_t_steps":3,
            "e_min":energy_in,
            "e_max":energy_in,
            "n_e_steps":1,
            "pid":-13,
            "default_hit":{"mass":xboa.common.pdg_pid_to_mass[13]}
        }, {
            "type":"longitudinal_ellipse",
            "delta_t":sigma_t, # ns
            "delta_e":sigma_e, # ns
            "t_centre":0.0,
            "e_centre":energy_in,
            "n_per_dimension":8,
            "pid":-13,
        }][0:],
    }
    reference = {
        "p_start":pz_in,
    }
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(lattice_filename)
    my_linac.cleanup_dir = cleanup
    my_linac.max_step = 1.0
    my_linac.elements = cavities
    my_linac.beam = beam_def
    my_linac.reference = reference
    my_linac.do_stochastics = 0 # e.g. decays
    my_linac.max_z = max([cav["z_position"] for cav in cavities])
    my_linac.build_linac()
    return my_linac

class Analysis():
    def __init__(self, linac, plot_dir):
        self.out_dir = linac.out_dir()
        self.plot_dir = plot_dir
        self.clean_dir = True
        self.out_filename = os.path.join(self.out_dir, linac.output_file)+".txt"
        self.t_targets = []
        self.e_targets = []

    def load_data(self):
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)

    def do_plots(self):
        g4bl_interface.g4bl_interface.clean_dir(self.plot_dir, self.clean_dir)
        self.plot_time_energy_event()
        self.plot_time_energy_station()
        self.plot_s_energy_station()

    def get_time_energy(self, bunch):
        t_list = bunch.list_get_hit_variable(["t"], ["ns"])[0]
        e_list = bunch.list_get_hit_variable(["kinetic_energy"], ["ns"])[0]
        t_list = [t-t_list[0] for t in t_list]
        return t_list, e_list

    def get_s_energy(self, bunch):
        clight = xboa.common.constants["c_light"]
        s_list = [(hit["t"]-bunch[0]["t"])*hit["pz"]/hit["energy"]*clight for hit in bunch]
        e_list = bunch.list_get_hit_variable(["kinetic_energy"], ["ns"])[0]
        s_list = [s-s_list[0] for s in s_list]
        return s_list, e_list

    def plot_time_energy_station(self):
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)

        bunch_start = self.bunch_list[0]
        t_list, e_list = self.get_time_energy(bunch_start)
        axes.scatter(t_list, e_list, label=f"z {bunch_start[0]['z']} mm")

        bunch_end = self.bunch_list[-1]
        t_list, e_list = self.get_time_energy(bunch_end)
        axes.scatter(t_list, e_list, label=f"z {bunch_end[0]['z']} mm")

        axes.set_xlabel("$\\Delta$t [ns]")
        axes.set_ylabel("Kinetic energy [MeV]")
        axes.legend()
        self.plot_lines(axes)

        figure.savefig(os.path.join(self.plot_dir, "time_energy.png"))

    def plot_s_energy_station(self):
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)

        bunch_start = self.bunch_list[0]
        s_list, e_list = self.get_s_energy(bunch_start)
        axes.scatter(s_list, e_list, label=f"z {bunch_start[0]['z']} mm")

        bunch_end = self.bunch_list[-1]
        s_list, e_list = self.get_s_energy(bunch_end)
        axes.scatter(s_list, e_list, label=f"z {bunch_end[0]['z']} mm")

        axes.set_xlabel("$\\Delta$z [ns]")
        axes.set_ylabel("Kinetic energy [MeV]")
        axes.legend()
        vz = self.bunch_list[0][0]["pz"]/self.bunch_list[0][0]["energy"]*xboa.common.constants["c_light"]
        self.t_targets = [t*vz for t in self.t_targets]
        self.plot_lines(axes)
        self.t_targets = [t/vz for t in self.t_targets]

        figure.savefig(os.path.join(self.plot_dir, "s_energy.png"))

    def plot_time_energy_event(self):
        track_dict = {}
        bunch_list = copy.deepcopy(self.bunch_list)
        e_0, t_0 = [], []
        for bunch in bunch_list:
            hit_list = bunch.hits()
            for hit in reversed(hit_list):
                hit['t'] -= hit_list[0]['t']
                hit['energy'] -= hit_list[0]['energy'] # no longer on shell!!!
                ev = hit['event_number']
                if ev not in track_dict:
                    track_dict[ev] = []
                track_dict[ev].append(hit)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        for ev, track in track_dict.items():
            e_list = [hit['energy'] for hit in track]
            t_list = [hit['t'] for hit in track]
            t_0.append(t_list[0])
            e_0.append(e_list[0])
            axes.scatter(t_list, e_list, s=1)
        axes.scatter(t_0, e_0, color="orange")
        figure.savefig(os.path.join(self.plot_dir, "contours.png"))

    def plot_lines(self, axes):
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        for e in self.e_targets:
            axes.plot(xlim, [e,e], color="grey", linestyle="dashed")
        for t in self.t_targets:
            axes.plot([t,t], ylim, color="grey", linestyle="dashed")
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)


def main():
    do_execute = True
    my_linac = build_final_cooling_lattice(do_execute)
    if do_execute:
        my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
        my_execution.execute()
    my_analysis = Analysis(my_linac, os.path.split(my_linac.lattice_filename)[0]+"/plots")
    e_centre = 25.7 # 21.2 #
    my_analysis.e_targets = [e_centre-3.8, e_centre, e_centre+3.8]
    dt = 1700.0/0.55/xboa.common.constants["c_light"]
    my_analysis.t_targets = [-dt, 0.0, dt]
    my_analysis.load_data()
    my_analysis.do_plots()



if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to end")
