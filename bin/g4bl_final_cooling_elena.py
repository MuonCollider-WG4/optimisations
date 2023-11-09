import numpy
import math
import json
import copy
import os
import matplotlib
import g4bl_interface.g4bl_interface
import xboa.hit
import xboa.bunch

def absorber(absorber_length):
    absorber = [{
        "name":"absorber",
        "length":absorber_length,
        "outer_radius":1000.0,
        "material":"lH2",
        "z_position":absorber_length/2+0.001,
        "type":"tube",
    }]
    return absorber

def multiharmonic_rf(n_cells, principle_harmonic, v1, v2, v3, v4, phi1, phi2, phi3, phi4, z_offset=0.0):
    cavities = [{
            "name":f"pillbox_{i}.1",
            "inner_length":200.0,
            "frequency":principle_harmonic*1,
            "max_gradient":v1,
            "phase":phi1,
            "rgb":[0.5,1,0],
            "z_position":150+1000.0*i+z_offset,
            "type":"cavity",
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.2",
            "inner_length":200.0,
            "frequency":principle_harmonic*2,
            "max_gradient":v2,
            "phase":phi2,
            "rgb":[1,1,0],
            "z_position":400+1000.0*i+z_offset,
            "type":"cavity",
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.3",
            "inner_length":200.0,
            "frequency":principle_harmonic*3,
            "max_gradient":v3,
            "phase":phi3,
            "rgb":[1,0.5,0],
            "z_position":650+1000.0*i+z_offset,
            "type":"cavity",
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.4",
            "inner_length":200.0,
            "frequency":principle_harmonic*4,
            "max_gradient":v4,
            "phase":phi4,
            "rgb":[0.5,0.5,0],
            "z_position":900.0+1000.0*i+z_offset,
            "type":"cavity",
        } for i in range(n_cells)
    ]
    return cavities

def build_final_cooling_lattice(cleanup = True, target_cell=1):
    lattice_filename = "output/elena_v2/linac.g4bl"
    src_path = "../elena_rftrack_to_g4bl/"

    lattice_def = open(os.path.join(src_path, "channel_params_solenoid_RF.json")).read()
    lattice_def = json.loads(lattice_def)

    # absorber
    absorber_length = lattice_def["cell"][target_cell-1]["abs_len"]*1000.0
    my_absorber = absorber(absorber_length)
    ref_t = 0.0 #20.5 # I haven't got the z position right

    # beam
    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(lattice_filename)[0],
        "z_position":0.0,
        "beams":[{
            "type":"longitudinal_g4bl_trackfile",
            "file_name":os.path.join(src_path, f"g4bl_beam_init"), #g4bl_beam_cooled_cell_{target_cell}"),
            "n_particles":1000,
        }],
    }
    reference = {
        "t_start":ref_t,
        "p_start":lattice_def["cell"][target_cell-1]["pz"],
    }

    # RF
    frequency = 0.07115
    voltage_1 = 0.0 #15.8
    voltage_2 = 0.00001
    phi_1 = 0.0
    tune_coeff = 2 #*int(6/(frequency*voltage_1)**0.5)
    if tune_coeff < 1:
        tune_coeff = 1
    if tune_coeff > 100:
        tune_coeff = 100
    print("Simulating", tune_coeff, "cells")
    cavities = multiharmonic_rf(tune_coeff, frequency,
                                voltage_1, voltage_2, 0.0, 0.0,
                                phi_1, 0.0, 0.0, 0.0,
                                z_offset=absorber_length)

    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(lattice_filename)
    my_linac.cleanup_dir = cleanup
    my_linac.max_step = 1.0
    my_linac.elements = my_absorber+cavities
    my_linac.beam = beam_def
    my_linac.reference = reference
    my_linac.do_stochastics = 0 # e.g. decays
    my_linac.min_z = 0 #(math.floor(absorber_length / 100)+1)*100
    my_linac.max_z = max([cav["z_position"] for cav in cavities])
    my_linac.build_linac()
    return my_linac

class Analysis():
    def __init__(self, linac, plot_dir, source_bunch_filename, target_bunch_filename):
        self.out_dir = linac.out_dir()
        self.plot_dir = plot_dir
        self.clean_dir = True
        self.out_filename = os.path.join(self.out_dir, linac.output_file)+".txt"
        self.source_bunch_filename = source_bunch_filename
        self.target_bunch_filename = target_bunch_filename
        self.t_targets = []
        self.e_targets = []

    def load_data(self):
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)
        self.target_bunch = xboa.bunch.Bunch.new_from_read_builtin("g4beamline_bl_track_file", self.target_bunch_filename)
        self.source_bunch = xboa.bunch.Bunch.new_from_read_builtin("g4beamline_bl_track_file", self.source_bunch_filename)

    def do_plots(self):
        g4bl_interface.g4bl_interface.clean_dir(self.plot_dir, self.clean_dir)
        self.plot_time_energy_event()
        self.plot_time_energy_station()
        self.plot_s_energy_station()
        self.plot_time_amplitude()

    def get_time_energy(self, bunch, first_is_reference = True):
        t_list = bunch.list_get_hit_variable(["t"], ["ns"])[0]
        e_list = bunch.list_get_hit_variable(["kinetic_energy"], ["ns"])[0]
        if first_is_reference:
            t_list = [t-t_list[0] for t in t_list]
        else:
            mean_t = numpy.mean(t_list)
            t_list = [t-mean_t for t in t_list]
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
        psize = 2

        t_list, e_list = self.get_time_energy(self.source_bunch, False)
        axes.scatter(t_list, e_list, s=psize, label=os.path.split(self.source_bunch_filename)[1])

        t_list, e_list = self.get_time_energy(self.target_bunch, False)
        axes.scatter(t_list, e_list, s=psize, label=os.path.split(self.target_bunch_filename)[1])

        bunch_start = self.bunch_list[0]
        t_list, e_list = self.get_time_energy(bunch_start)
        axes.scatter(t_list, e_list, s=psize, label=f"z {bunch_start[0]['z']} mm")

        bunch_end = self.bunch_list[-1]
        t_list, e_list = self.get_time_energy(bunch_end)
        axes.scatter(t_list, e_list, s=psize, label=f"z {bunch_end[0]['z']} mm")

        axes.set_xlabel("$\\Delta$t [ns]")
        axes.set_ylabel("Kinetic energy [MeV]")
        axes.legend()
        self.plot_lines(axes)

        figure.savefig(os.path.join(self.plot_dir, "time_energy.png"))


    def plot_time_amplitude(self):
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)

        t_list, e_list = self.get_time_energy(self.source_bunch, False)
        amplitude_list = [self.source_bunch.get_hit_variable(hit, "amplitude x y") for hit in self.source_bunch]

        axes.scatter(t_list, amplitude_list)
        axes.set_xlabel("$\\Delta$t [ns]")
        axes.set_ylabel("Amplitude [mm]")
        axes.legend()
        self.plot_lines(axes)

        figure.savefig(os.path.join(self.plot_dir, "time_amplitude.png"))


    def plot_s_energy_station(self):
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)

        bunch_start = self.bunch_list[0]
        s_list, e_list = self.get_s_energy(bunch_start)
        axes.scatter(s_list, e_list, label=f"z {bunch_start[0]['z']} mm")

        bunch_end = self.bunch_list[-1]
        s_list, e_list = self.get_s_energy(bunch_end)
        axes.scatter(s_list, e_list, label=f"z {bunch_end[0]['z']} mm")

        axes.set_xlabel("$\\Delta$z [mm]")
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
        axes.set_xlabel("$\\Delta$t [ns]")
        axes.set_ylabel("$\\Delta$ energy [MeV]")
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
    source_file = "../elena_rftrack_to_g4bl/g4bl_beam_init"
    #source_file = "../elena_rftrack_to_g4bl/g4bl_beam_cooled_cell_1"
    target_file = "../elena_rftrack_to_g4bl/g4bl_beam_cooled_cell_1"
    my_analysis = Analysis(my_linac, os.path.split(my_linac.lattice_filename)[0]+"/plots", source_file, target_file)
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
