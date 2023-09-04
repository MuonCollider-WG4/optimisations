import math
import copy
import os
import matplotlib
import xboa.hit
import xboa.bunch

import g4bl_interface.g4bl_interface
import models.field_models

class Chicane:
    def __init__(self):
        # chicane parameters
        self.bend_angle = 30 # degrees
        self.radius_of_curvature = 2000.0 # mm
        self.b_field = 2 # T
        self.pre_length = 1000.0 # mm, length before the beam is injected
        self.length_before = 1000.0 # mm
        self.length_after = 5000.0 # mm
        self.stroke_length = 100.0 # mm
        self.name = "chicane"
        self.coil_list = []

        # coil parameters
        self.coil_length = 95.0
        self.inner_radius = 200.0
        self.outer_radius = 250.0
        self.current = 10.0

    def calculate_current(self):
        block = models.field_models.CurrentBlock(
            self.b_field, self.stroke_length*1e-3/2.0, self.coil_length*1e-3,
            self.inner_radius*1e-3, self.outer_radius*1e-3, self.stroke_length*1e-3, 20, 10)
        bz_avg = 0.0
        for iz in range(21):
            z_pos = self.stroke_length*1e-3/10*iz
            bz = block.get_field(z_pos)
            bz_avg += bz/21
        block.b0 *= self.b_field/bz_avg
        block.reset()
        bz_avg = 0.0
        for iz in range(21):
            z_pos = self.stroke_length*1e-3/10*iz
            bz = block.get_field(z_pos)
            bz_avg += bz/21
        self.current = block.get_current_density()*1e-6
        print(f"Adjusted current to {self.current} A mm^-2 to get bz_avg {bz_avg} T")



    def build(self):
        self.calculate_current()
        self.build_line(self.length_before+self.pre_length, 0.0, -self.pre_length)
        self.build_bend(self.coil_list[-1], +1.0, {"rgb":[0.5, 1.0, 0.0]})
        self.build_bend(self.coil_list[-1], -1.0, {"rgb":[0.0, 1.0, 0.5]})
        self.build_line(self.length_after, self.coil_list[-1]["x_position"], self.coil_list[-1]["z_position"]+self.stroke_length)
        for coil in self.coil_list:
            coil["x_position"] -= self.coil_list[-1]["x_position"]
        return self.coil_list

    def build_line(self, length, x_start, z_start):
        n_coils = int(length/self.stroke_length)
        for i in range(n_coils+1):
            coil = self.default_coil()
            coil["z_position"] += z_start+self.stroke_length*i
            coil["x_position"] = x_start
            coil["name"] = f"chicane_straight_{i}"
            self.coil_list.append(coil)

    def build_bend(self, last_coil, bend_scale, modifiers={}):
        rot_start = math.radians(last_coil["y_rotation"])
        x_start = last_coil["x_position"]+math.sin(rot_start)*self.stroke_length
        z_start = last_coil["z_position"]+math.cos(rot_start)*self.stroke_length
        dtheta = self.stroke_length/self.radius_of_curvature
        n_coils = int(math.radians(self.bend_angle)/dtheta)
        theta_list = [dtheta*i for i in range(1, n_coils+1)]
        x_list = [self.radius_of_curvature*math.cos(theta)*bend_scale for theta in theta_list]
        z_list = [self.radius_of_curvature*math.sin(theta) for theta in theta_list]
        for i in range(n_coils):
            x_tmp = x_list[i]*math.cos(rot_start)+z_list[i]*math.sin(rot_start)
            z_tmp = -x_list[i]*math.sin(rot_start)+z_list[i]*math.cos(rot_start)
            x_list[i] = x_tmp
            z_list[i] = z_tmp
        theta_list = [(theta-theta_list[0])*bend_scale-rot_start for theta in theta_list]
        x_list = [x - x_list[0]+x_start for x in x_list]
        z_list = [z - z_list[0]+z_start for z in z_list]

        for i in range(n_coils):
            coil = self.default_coil()
            coil.update(modifiers)
            coil["x_position"] = x_list[i]
            coil["z_position"] = z_list[i]
            coil["y_rotation"] = -math.degrees(theta_list[i])
            coil["name"] = f"chicane_bend_{i}"
            self.coil_list.append(coil)

    def default_coil(self):
        coil = {
            "x_position":0.0,
            "y_rotation":0.0, # rotation about the y axis i.e. in x z plane
            "z_position":0.0,
            "inner_radius":self.inner_radius,
            "outer_radius":self.outer_radius,
            "length":self.coil_length,
            "current":self.current,
            "name":self.name,
            "type":"solenoid",
        }
        return coil

def get_beam(my_linac):
    x_start = my_linac.elements[0]["x_position"]
    print(f"Offsetting beam to x={x_start} mm")
    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(my_linac.lattice_filename)[0],

        "beams":[{
            "type":"longitudinal_grid",
            "t_min":0.0, # ns
            "t_max":0.0, # ns
            "n_t_steps":1,
            "e_min":50,
            "e_max":500,
            "n_e_steps":10,
            "pid":-13,
            "default_hit":{"x":x_start, "mass":xboa.common.pdg_pid_to_mass[13], "pid":-13},
        }, {
            "type":"longitudinal_grid",
            "t_min":0.0, # ns
            "t_max":0.0, # ns
            "n_t_steps":1,
            "e_min":4000,
            "e_max":10000,
            "n_e_steps":7,
            "pid":2212,
            "default_hit":{"x":x_start, "mass":xboa.common.pdg_pid_to_mass[2212], "pid":2212},
        }],
    }
    mass = xboa.common.pdg_pid_to_mass[13]
    hit = xboa.hit.Hit.new_from_dict({"pid":-13, "energy":120.0+mass, "mass":mass}, "pz")
    reference = {
        "p_start":hit["pz"],
    }

    my_linac.beam = beam_def
    my_linac.reference = reference

def build_chicane(cleanup = True):
    lattice_filename = "output/chicane_1/chicane.g4bl"
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(lattice_filename)
    my_linac.cleanup_dir = cleanup
    my_linac.max_step = 1.0
    my_linac.elements = Chicane().build()
    get_beam(my_linac)
    my_linac.do_stochastics = 0 # e.g. decays
    my_linac.max_z = max([sol["z_position"] for sol in my_linac.elements])
    my_linac.build_linac()
    return my_linac

class Analysis():
    def __init__(self, linac, plot_dir):
        self.linac = linac
        self.out_dir = linac.out_dir()
        self.plot_dir = plot_dir
        self.clean_dir = True
        self.out_filename = os.path.join(self.out_dir, linac.output_file)+".txt"

    def load_data(self):
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)

    def do_plots(self):
        g4bl_interface.g4bl_interface.clean_dir(self.plot_dir, self.clean_dir)
        self.plot_y_z()

    @classmethod
    def get_b_tot(cls, hit):
        return 1e3*(hit["bx"]**2+hit["by"]**2+hit["bz"]**2)**0.5

    def get_amplitude(cls, hit):
        p = hit["p"]
        btot = cls.get_b_tot(hit)
        beta_t = p/0.15/btot
        mean = {"x":0.0, "y":0.0, "px":0.0, "py":0.0}
        cov = xboa.bunch.Bunch.build_penn_ellipse(1.0, hit["mass"], beta_t, 0, p, 0, btot*1e-3, 1)
        if False and p < 220 and p > 160:
            print(f"Get Amplitude {p}, {btot}, {beta_t}, {hit['mass']}\n{cov}")
        bunchish = xboa.bunch.Bunch.new_from_hits([hit])
        amplitude = xboa.bunch.Bunch.get_amplitude(bunchish, hit, ["x", "y"], covariance_matrix=cov, mean_dict=mean)
        return amplitude

    def plot_y_z(self):
        plot_reference = False
        track_dict = {}
        bunch_list = copy.deepcopy(self.bunch_list)
        for bunch in bunch_list:
            hit_list = bunch.hits()
            for hit in reversed(hit_list):
                ev = hit['event_number']
                if ev not in track_dict:
                    track_dict[ev] = []
                track_dict[ev].append(hit)
        if not plot_reference:
            del track_dict[0]
        get_dict = {
            "x":lambda hit: hit["x"], "y":lambda hit: hit["y"], "z":lambda hit: hit["x"],
            "btot":self.get_b_tot, "amp":self.get_amplitude,
        }
        label_dict = {"x":"x [mm]", "y":"y [mm]", "btot":"B$_{tot}$ [T]", "amp":"A [mm]"}
        for yvar in ["x", "y", "btot", "amp"]:
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)
            linestyle_dict = {-13:"-", 2212:"dotted"}
            for ev, track in sorted(track_dict.items()):
                ls = linestyle_dict[track[0]["pid"]]
                x_list = [hit['z'] for hit in track]
                y_list = [get_dict[yvar](hit) for hit in track]
                axes.plot(x_list, y_list, linestyle=ls, label=f"KE {round(track[0]['kinetic_energy'])} MeV")
            if yvar == "x" and len(track_dict):
                axes.plot(x_list, [self.linac.elements[0]["x_position"] for point in x_list], linestyle="dashed", color="grey")
                axes.plot(x_list, [self.linac.elements[-1]["x_position"] for point in x_list], linestyle="dashed", color="grey")
            axes.set_xlabel("z [mm]")
            axes.set_ylabel(label_dict[yvar])
            axes.legend()
            figure.savefig(os.path.join(self.plot_dir, f"{yvar}-z.png"))


def main():
    do_execute = True
    my_linac = build_chicane(do_execute)
    if do_execute: # alternatively we just redo the analysis
        my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
        my_execution.execute()
    my_analysis = Analysis(my_linac, os.path.split(my_linac.lattice_filename)[0]+"/plots")
    my_analysis.load_data()
    my_analysis.do_plots()



if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to end")
