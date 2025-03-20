import math
import sys
import glob
import matplotlib
import copy
import operator
import os
import json
import xboa
import numpy
import g4bl_interface.g4bl_interface
import g4bl_interface.stripper

import front_end.plotting.target_yield_analysis
import front_end.plotting.target_outer_detector_analysis

class TargetRegion:
    def __init__(self):
        self.element_list = []
        self.target_length = 800.0
        self.target_radius = 30.0
        self.target_tilt = 0.0 # degrees
        self.target_z_start = 50.0
        self.beam_sigma = 5.0
        self.target_material = "C"
        self.beam_pipe_radius = 605.0
        self.detector_radius = 600.0
        self.detector_length = 10000.0
        self.max_z = 50000.0
        self.solenoid_file = "share/target_solenoid.tex"
        self.solenoid_max_derivative = 3
        self.target_rgb = [0.8, 0.8, 0.8]
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
            "max_derivative":self.solenoid_max_derivative,
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
            "y_rotation":self.target_tilt,
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

    def build_outer_radius_detector(self):
        detector_tube = {
            "type":"tube",
            "name":"outer_radius_detector_tube",
            "inner_radius":self.detector_radius,
            "outer_radius":self.detector_radius+1,
            "z_position":None,
            "x_position":0.0,
            "length":self.detector_length,
            "material":"Vacuum",
            "type":"tube",
        }
        detector = {
            "type":"detector",
            "name":"outer_radius_detector",
            "solid":"outer_radius_detector_tube",
            "coordinates":"global",
            "format":"for009.dat",
            "filename":"outer_radius_detector.txt",
            "z_position":self.detector_length/2.0,
        }
        self.element_list.append(detector_tube)
        self.element_list.append(detector)

    def build_beam_pipe(self):
        beam_pipe = {
            "name":"beam_pipe",
            "inner_radius":self.beam_pipe_radius+0.1,
            "outer_radius":2*self.beam_pipe_radius+0.1,
            "z_position":self.max_z/2,
            "x_position":0.0,
            "length":self.max_z*2,
            "material":"Vacuum",
            "will_kill":"True",
            "type":"tube",
        }
        self.element_list.append(beam_pipe)

def get_beam(my_linac, target_region, start_event, n_particles, kinetic_energy):
    pid = 2212
    mass = xboa.common.pdg_pid_to_mass[pid]
    p = ((kinetic_energy+mass)**2-mass**2)**0.5
    print("GET BEAM p", p, "m",  mass, "ke", kinetic_energy)
    px = 1e-3*p
    target_tilt = target_region.target_tilt
    target_length = target_region.target_length/2+1.0
    target_z = target_region.target_length/2+target_region.target_z_start
    beam_x = -target_length*math.sin(math.radians(target_tilt))
    beam_z = target_z-target_length*math.cos(math.radians(target_tilt))
    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(my_linac.lattice_filename)[0],
        "x_position":beam_x,
        "z_position":beam_z,
        "y_rotation":target_tilt,
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

def build_target(base_dir, start_event, n_particles, energy, end_field_length, target_magnitude, cw_magnitude, target_tilt, cleanup_dir, will_gen_beam):
    target_region = TargetRegion()
    target_region.max_z = 50000.0
    target_region.z_start = 50.0
    target_region.target_tilt = target_tilt # degree
    target_region.target_length = 800.0 # 1.0*2/19.30 #
    target_region.target_material ="C" # "W" #
    target_region.build_target()

    target_region.end_length = end_field_length # mm
    target_region.peak_field = target_magnitude # T
    target_region.low_field = cw_magnitude # T
    target_region.build_outer_radius_detector()
    target_region.build_tanh_solenoid()
    target_region.build_beam_pipe()
    #target_region.build_solenoid_from_latex()

    #print(json.dumps(target_region.element_list, indent=2))
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(os.path.join(base_dir, "target_region.g4bl"))
    my_linac.cleanup_dir = cleanup_dir
    my_linac.elements = target_region.element_list
    my_linac.do_stochastics = 1 # e.g. decays
    my_linac.max_z = target_region.max_z+1
    my_linac.physics.physics_list = "QGSP_BERT"
    my_linac.track_cuts = g4bl_interface.g4bl_interface.TrackCuts()
    my_linac.track_cuts.keep = [2212, -13, 13, 211, -211]
    my_linac.z_spacing = 10000.0
    if will_gen_beam:
        get_beam(my_linac, target_region, start_event, n_particles, energy)
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
    clear_dirs = False # recalculates field maps
    do_execute = False
    do_analysis = True
    base_dir = "output/mucol_target_v6/"
    target_field = 20
    uniform_field = 1.5
    end_length = 1000
    energy = 10000
    n_versions = 10
    n_events = 1000
    for target_tilt, end_length in [(10, 1.0)]:#, (0, 1.0), (0, 2.0), (8, 2.0), (6, 2.0)]:
        end_length *= 1000
        for energy in [10000, 5000]:
            for version in range(0, n_versions):
                run_dir = base_dir+dir_name({"energy":energy, "version":version, "end_length":end_length, "target_tilt":target_tilt})
                my_linac, target_region = build_target(run_dir, version*n_events, n_events, energy, end_length, target_field, uniform_field, target_tilt, clear_dirs, do_execute)
                my_linac.build_linac()
                if do_execute: # else we just redo the analysis
                    my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
                    my_execution.execute()
                    #strip_file(run_dir)
            beam_file_glob = base_dir+dir_name({"energy":energy, "version":"*", "end_length":end_length, "target_tilt":target_tilt})+"/output_data.txt"
            if do_analysis:
                target = f"{target_region.target_length:6.3g} mm {target_region.target_material}"
                title = f"energy: {energy} [MeV]   n$_{{protons}}$: {float(n_versions*n_events):4.6g} end length: {end_length/1000} [m] tilt: {target_tilt} [deg]"
                my_analysis = front_end.plotting.target_yield_analysis.Analysis(base_dir)
                my_analysis.z_analysis = target_region.max_z
                plot_dir =  base_dir+dir_name({"energy":energy, "end_length":end_length, "target_tilt":target_tilt})+"_plots"
                will_close = end_length != 1000 and end_length != 9000
                try:
                    my_analysis.load_data(my_linac, target_region, beam_file_glob, plot_dir)
                    my_analysis.do_plots(title, will_close)
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    print("Failed to plot yield")
                analysis_outer = front_end.plotting.target_outer_detector_analysis.Analysis()
                analysis_outer.filename = base_dir+dir_name({"energy":energy, "version":"*", "end_length":end_length, "target_tilt":target_tilt})+"/outer_radius_detector.txt"
                analysis_outer.plot_dir = plot_dir
                analysis_outer.title = title
                try:
                    analysis_outer.load_data()
                    analysis_outer.do_plots()
                except Exception:
                    sys.excepthook(*sys.exc_info())
                    print("Failed to plot detector")

    print("Finished")
    if False: #do_analysis:
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
