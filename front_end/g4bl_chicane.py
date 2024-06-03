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
        self.bend_angle = 12.5 # degrees
        self.radius_of_curvature_1 = 20000.0 # mm
        self.radius_of_curvature_2 = 20000.0 # mm
        self.b_field = 2 # T
        self.middle_straight = 0.0
        self.pre_length = 1000.0 # mm, length before the beam is injected
        self.length_before = 1000.0 # mm
        self.length_after = 1000.0 # mm
        self.post_length_1 = 100.0 # mm
        self.post_length_2 = 1000.0 # mm, analysis will exclude the last bit where field uniformity is poor
        self.stroke_length = 100.0 # mm
        self.name = "chicane"
        self.coil_list = []

        # coil parameters
        self.coil_length = 95.0
        self.inner_radius = 500.0
        self.outer_radius = 600.0
        self.current = 10.0

    def calculate_current(self, coil = None):
        if coil == None:
            coil = self.default_coil()
        block = models.field_models.CurrentBlock(
            self.b_field, self.stroke_length*1e-3/2.0, coil["length"]*1e-3,
            coil["inner_radius"]*1e-3, coil["outer_radius"]*1e-3, self.stroke_length*1e-3, 20, 10)
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
        return block.get_current_density()*1e-6

    def bend_1_modifier(self, i_scale, coil):
        return {"rgb":[0.5, 1.0, 0.0]}

    def bend_2_modifier(self, i_scale, coil):
        return {"rgb":[0.0, 1.0, 0.5]}

    def build(self):
        self.current = self.calculate_current()
        self.build_line(self.length_before+self.pre_length, 0.0, -self.pre_length)
        self.build_bend(self.coil_list[-1], +1.0, self.bend_1_modifier)
        self.build_middle_straight()
        self.build_bend(self.coil_list[-1], -1.0, self.bend_2_modifier)
        self.build_line(self.length_after, self.coil_list[-1]["x_position"], self.coil_list[-1]["z_position"]+self.stroke_length)
        self.build_line(self.post_length_1, self.coil_list[-1]["x_position"],
                        self.coil_list[-1]["z_position"]+self.stroke_length)
        self.build_line(self.post_length_2, self.coil_list[-1]["x_position"],
                        self.coil_list[-1]["z_position"]+self.stroke_length)
        for coil in self.coil_list:
            coil["x_position"] -= self.coil_list[-1]["x_position"]
        return self.coil_list

    def build_middle_straight(self):
        y_rot = math.radians(self.coil_list[-1]["y_rotation"])
        x_start = self.coil_list[-1]["x_position"]+math.sin(y_rot)*self.stroke_length
        z_start = self.coil_list[-1]["z_position"]+math.cos(y_rot)*self.stroke_length
        if self.middle_straight <= 0:
            return
        self.build_line(self.middle_straight, x_start, z_start, y_rot)

    def apply_mod(self, i_scale, coil, modifier):
        if modifier is None or modifier == {}:
            return
        if type(modifier) != type({}):
            modifier = modifier(i_scale, coil) # should be a dict
        coil.update(modifier)

    def build_line(self, length, x_start, z_start, y_rotation=0.0, modifier={}):
        n_coils = int(length/self.stroke_length)
        for i in range(n_coils+1):
            coil = self.default_coil()
            coil["x_position"] = x_start+math.sin(y_rotation)*self.stroke_length*i
            coil["z_position"] += z_start+math.cos(y_rotation)*self.stroke_length*i
            coil["y_rotation"] = -y_rotation
            coil["name"] = f"chicane_straight_{i}"
            self.apply_mod(i/n_coils, coil, modifier)
            self.coil_list.append(coil)

    def build_bend(self, last_coil, bend_scale, modifier={}):
        rot_start = math.radians(last_coil["y_rotation"])
        x_start = last_coil["x_position"]+math.sin(rot_start)*self.stroke_length
        z_start = last_coil["z_position"]+math.cos(rot_start)*self.stroke_length
        dtheta = self.stroke_length/self.radius_of_curvature_1
        n_coils = int(math.radians(self.bend_angle)/dtheta)
        theta_list = [dtheta*i for i in range(1, n_coils+1)]
        x_list = [self.radius_of_curvature_1*math.cos(theta)*bend_scale for theta in theta_list]
        z_list = [self.radius_of_curvature_1*math.sin(theta) for theta in theta_list]
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
            coil["x_position"] = x_list[i]
            coil["z_position"] = z_list[i]
            coil["y_rotation"] = -math.degrees(theta_list[i])
            coil["name"] = f"chicane_bend_{i}"
            self.apply_mod(i/(n_coils-1), coil, modifier)
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
            "e_min":1,
            "e_max":1,
            "n_e_steps":1,
            "pid":-13,
            "default_hit":{"x":x_start, "mass":xboa.common.pdg_pid_to_mass[13], "pid":-13},
        }, {
            "type":"longitudinal_grid",
            "t_min":0.0, # ns
            "t_max":0.0, # ns
            "n_t_steps":1,
            "e_min":100,
            "e_max":1000,
            "n_e_steps":10,
            "pid":-13,
            "default_hit":{"x":x_start, "mass":xboa.common.pdg_pid_to_mass[13], "pid":-13},
        }, {
            "type":"longitudinal_grid",
            "t_min":0.0, # ns
            "t_max":0.0, # ns
            "n_t_steps":1,
            "e_min":4999,
            "e_max":5001,
            "n_e_steps":1,
            "pid":2212,
            "default_hit":{"x":x_start, "mass":xboa.common.pdg_pid_to_mass[2212], "pid":2212},
        }][0:],
    }
    mass = xboa.common.pdg_pid_to_mass[13]
    hit = xboa.hit.Hit.new_from_dict({"pid":-13, "energy":120.0+mass, "mass":mass}, "pz")
    reference = {
        "p_start":hit["pz"],
    }

    my_linac.beam = beam_def
    my_linac.reference = reference

chicane = None
def coil_mod(i_scale, coil):
    global chicane
    coil_variant = {"rgb":[0.0, 0.0, 1.0], "inner_radius":200.0, "outer_radius":300.0}
    modded_coil = chicane.default_coil()
    modded_coil.update(coil_variant)
    current = chicane.calculate_current(modded_coil)
    coil_variant["current"] = current
    if abs(coil["x_position"]) > 400.0:
        coil.update(coil_variant)
    else:
        pass
    ##import json
    #print(json.dumps(coil))
    return coil

def build_chicane(root_dir, b_field, bend_angle, r_curv, cleanup = True):
    global chicane
    lattice_filename = root_dir+"/chicane.g4bl"
    g4bl_interface.g4bl_interface.Solenoid.clear()
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(lattice_filename)
    my_linac.cleanup_dir = cleanup
    my_linac.max_step = 10.0
    chicane = Chicane()
    chicane.bend_angle = bend_angle # degrees
    chicane.radius_of_curvature_1 = r_curv # mm
    chicane.radius_of_curvature_2 = r_curv # mm
    chicane.b_field = b_field # T
    chicane.middle_straight = 0.0
    # get bend 2
    chicane.bend_1_modifier = coil_mod
    chicane.bend_2_modifier = coil_mod
    # build_linac
    my_linac.elements = chicane.build()
    get_beam(my_linac)
    my_linac.do_stochastics = 0 # e.g. decays
    my_linac.max_z = max([sol["z_position"] for sol in my_linac.elements])-chicane.post_length_2
    return my_linac, chicane

def build_charge_separation(root_dir, b_field, bend_angle, r_curv_1, r_curv_2, cleanup = True):
    global chicane
    lattice_filename = root_dir+"/chicane.g4bl"
    g4bl_interface.g4bl_interface.Solenoid.clear()
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(lattice_filename)
    my_linac.cleanup_dir = cleanup
    my_linac.max_step = 10.0
    chicane = Chicane()
    chicane.bend_angle = bend_angle # degrees
    chicane.radius_of_curvature_1 = r_curv_1 # mm
    chicane.radius_of_curvature_2 = r_curv_2 # mm
    chicane.b_field = b_field # T
    chicane.middle_straight = 0.0
    # get bend 2
    chicane.bend_1_modifier = coil_mod
    chicane.bend_2_modifier = coil_mod
    # build_linac
    my_linac.elements = chicane.build()
    get_beam(my_linac)
    my_linac.do_stochastics = 0 # e.g. decays
    my_linac.max_z = max([sol["z_position"] for sol in my_linac.elements])-chicane.post_length_2
    return my_linac, chicane

class Analysis():
    def __init__(self, base_dir):
        self.clean_dir = True
        self.base_dir = base_dir
        self.amp_axes = None
        self.z_axes = None
        self.r_axes = None
        self.plot_energy_range = [0, 500]

    def load_data(self, linac, chicane, plot_dir):
        self.chicane = chicane
        self.linac = linac
        self.out_dir = linac.out_dir()
        self.plot_dir = plot_dir
        self.out_filename = os.path.join(self.out_dir, linac.output_file)+".txt"
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)

    def do_plots(self):
        g4bl_interface.g4bl_interface.clean_dir(self.plot_dir, self.clean_dir)
        self.get_track_dict()
        self.plot_y_z()
        self.plot_amp_max()
        self.plot_z_max()
        self.plot_r_max()

    @classmethod
    def get_b_tot(cls, hit):
        return 1e3*(hit["bx"]**2+hit["by"]**2+hit["bz"]**2)**0.5

    def post_linac_start(self):
        return self.linac.max_z - self.chicane.length_after

    def post_linac_end(self):
        return self.linac.max_z

    def pre_linac_end(self):
        return self.chicane.length_before

    def get_amplitude(self, hit):
        p = hit["p"]
        if hit["z"] < self.pre_linac_end(): # if we are still in the pre linac
            hit["x"] -= self.linac.elements[0]["x_position"]
        elif hit["z"] < self.post_linac_start():
            return 0.0
        btot = self.get_b_tot(hit)
        beta_t = p/0.15/btot
        mean = {"x":0.0, "y":0.0, "px":0.0, "py":0.0}
        cov = xboa.bunch.Bunch.build_penn_ellipse(1.0, hit["mass"], beta_t, 0, p, 0, btot*1e-3, 1)
        if False and p < 220 and p > 160:
            print(f"Get Amplitude {p}, {btot}, {beta_t}, {hit['mass']}\n{cov}")
        bunchish = xboa.bunch.Bunch.new_from_hits([hit])
        amplitude = xboa.bunch.Bunch.get_amplitude(bunchish, hit, ["x", "y"], covariance_matrix=cov, mean_dict=mean)
        return amplitude

    def get_track_dict(self):
        plot_reference = False
        self.track_dict = {}
        bunch_list = copy.deepcopy(self.bunch_list)
        for bunch in bunch_list:
            hit_list = bunch.hits()
            for hit in reversed(hit_list):
                ev = hit['event_number']
                if ev not in self.track_dict:
                    self.track_dict[ev] = []
                self.track_dict[ev].append(hit)
        if not plot_reference:
            del self.track_dict[0]

    def run_label(self, chicane, linac):
        return f"b$_z$ {chicane.b_field:.3g} T; r$_{{curv}}$ {chicane.radius_of_curvature_1*1e-3:.3g}, {chicane.radius_of_curvature_2*1e-3:.3g} m; $\\theta$ {chicane.bend_angle}$^\\circ$"

    def plot_amp_max(self):
        energy_list = [track[0]["kinetic_energy"] for ev, track in sorted(self.track_dict.items())]
        amp_list = []
        for ev, track in sorted(self.track_dict.items()):
            if track[-1]["z"] > self.post_linac_start():
                amp_list.append(self.get_amplitude(track[-1]))
            else:
                amp_list.append(-1.0)
        x_range = [min(energy_list), max(energy_list)]
        if self.plot_energy_range:
            x_range = self.plot_energy_range
        red_energy_list, red_amp_list = [], []
        for i, amp in enumerate(amp_list):
            if amp >= 0:
                red_amp_list.append(amp)
                red_energy_list.append(energy_list[i])
        if self.amp_axes is None:
            figure = matplotlib.pyplot.figure()
            self.amp_axes = figure.add_subplot(1, 1, 1)
        self.amp_axes.plot(red_energy_list, red_amp_list, label=self.run_label(self.chicane, self.linac))
        self.amp_axes.set_xlabel("Kinetic Energy [MeV]")
        self.amp_axes.set_ylabel("A$_{\\perp}$ [mm]")
        self.amp_axes.set_xlim(x_range)
        self.amp_axes.legend()
        self.amp_axes.figure.savefig(os.path.join(self.base_dir, f"energy-amp.png"))
        self.amp_axes.set_xlim([0.0, 500.0])
        self.amp_axes.set_ylim([0.0, 75.0])


    def plot_r_max(self):
        energy_list = [track[0]["kinetic_energy"] for ev, track in sorted(self.track_dict.items())]
        ref_track = [track for track in sorted(self.track_dict.values(), key = lambda x: x[0]["energy"])][0]
        print("plot r max energy", ref_track[0]["energy"])
        x_ref = [hit["x"] for hit in ref_track]
        y_ref = [hit["y"] for hit in ref_track]
        r_list = []
        proton_r = None
        for ev, track in sorted(self.track_dict.items()):
            if track[-1]["z"] > self.post_linac_start():
                r_hit = []
                for i, hit in enumerate(track):
                    r_hit.append( ((hit["x"]-ref_track[i]["x"])**2+\
                                   (hit["y"]-ref_track[i]["y"])**2)**0.5 )
                r_list.append(max(r_hit))
                if track[0]["pid"] == 2212:
                    proton_r = max(r_hit)
            else:
                r_list.append(-1.0)
        x_range = [min(energy_list), max(energy_list)]
        if self.plot_energy_range:
            x_range = self.plot_energy_range
        if self.r_axes is None:
            figure = matplotlib.pyplot.figure()
            self.r_axes = figure.add_subplot(1, 1, 1)
        line_2d = self.r_axes.plot(energy_list, r_list, label=self.run_label(self.chicane, self.linac))[0]
        if proton_r != None:
            myc = line_2d.get_color()
            self.r_axes.plot([min(energy_list), max(energy_list)], [proton_r, proton_r], color=myc, linestyle="--")
        self.r_axes.set_xlabel("Kinetic Energy [MeV]")
        self.r_axes.set_ylabel("Max. Radial Deviation [mm]")
        self.r_axes.set_xlim()
        self.r_axes.legend()
        self.r_axes.figure.savefig(os.path.join(self.base_dir, f"energy-rmax.png"))
        self.r_axes.set_xlim([0.0, 500.0])


    def plot_z_max(self):
        chicane_end = self.post_linac_start()
        chicane_length = self.post_linac_start()-self.pre_linac_end()
        energy_list = [track[0]["kinetic_energy"] for ev, track in sorted(self.track_dict.items())]
        z_list = []
        for ev, track in sorted(self.track_dict.items()):
            my_z = [hit["z"] for hit in track]
            z_list.append((max(my_z)-chicane_end)/chicane_length)
            #print(z_list[-1], my_z)
        if self.z_axes is None:
            figure = matplotlib.pyplot.figure()
            self.z_axes = figure.add_subplot(1, 1, 1)
        # cut out particles that survived
        x_range = [min(energy_list), max(energy_list)]
        if self.plot_energy_range:
            x_range = self.plot_energy_range
        z_cut = (self.post_linac_end()-chicane_end)/chicane_length-0.1
        energy_list = [e for i,e in enumerate(energy_list) if z_list[i] < z_cut]
        z_list = [z for i,z in enumerate(z_list) if z_list[i] < z_cut]

        self.z_axes.plot(energy_list, z_list, label=self.run_label(self.chicane, self.linac))
        self.z_axes.set_xlabel("Kinetic Energy [MeV]")
        self.z_axes.set_ylabel("$\\frac{\\mathrm{Max(z)} - \\mathrm{Chicane\\ End} }{\\mathrm{Chicane\\ Length} }$")
        a_range = self.z_axes.get_ylim()
        if a_range[0] < -2:
            self.z_axes.set_ylim([-2, a_range[1]])
        self.z_axes.set_xlim(x_range)
        self.z_axes.legend()
        self.z_axes.figure.savefig(os.path.join(self.base_dir, f"energy-z.png"))

    def plot_y_z(self):
        track_dict = self.track_dict
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
                axes.plot([0.0, self.post_linac_end()], [self.linac.elements[0]["x_position"]]*2, linestyle="dashed", color="grey")
                axes.plot([0.0, self.post_linac_end()], [self.linac.elements[-1]["x_position"]]*2, linestyle="dashed", color="grey")
            axes.set_xlabel("z [mm]")
            axes.set_ylabel(label_dict[yvar])
            axes_lim = axes.get_xlim()
            axes.set_xlim(axes_lim[0], axes_lim[1]+(axes_lim[1]-axes_lim[0])*0.5)
            axes.legend()
            figure.savefig(os.path.join(self.plot_dir, f"{yvar}-z.png"))
            if figure.number > 10:
                matplotlib.pyplot.close(figure)

def proton_extraction():
    do_execute = True
    b_field = 1.5
    bend_angle = 15
    r_curv = 20 # metres
    middle_straight = 0.0
    base_dir = f"output/chicane_scan_v19/"
    my_analysis = Analysis(base_dir)
    scale = 1.8
    for r_curv, bend_angle, b_field, in [(40, 10.0, 1.0)]:
        run_dir = os.path.join(base_dir, f"bz={b_field:.3g}_angle={bend_angle:.3g}_rcurv={r_curv:.3g}")
        my_linac, my_chicane = build_chicane(run_dir, b_field, bend_angle, r_curv*1e3, do_execute)
        my_linac.build_linac()
        if do_execute: # alternatively we just redo the analysis
            my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
            my_execution.execute()
        my_analysis.load_data(my_linac, my_chicane, os.path.split(my_linac.lattice_filename)[0]+"/plots")
        my_analysis.do_plots()
        print()

def main():
    do_execute = True
    b_field = 1.5
    bend_angle = 15
    r_curv = 20 # metres
    middle_straight = 0.0
    base_dir = f"output/charge_sep_v1/"
    my_analysis = Analysis(base_dir)
    scale = 1.8
    for r_curv, bend_angle, b_field, in [(41, 8.33, 1.5)]:
        run_dir = os.path.join(base_dir, f"bz={b_field:.3g}_angle={bend_angle:.3g}_rcurv={r_curv:.3g}")
        my_linac, my_chicane = build_charge_separation(run_dir, b_field, bend_angle, r_curv*1e3, r_curv*1e3, do_execute)
        my_linac.build_linac()
        if do_execute: # alternatively we just redo the analysis
            my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
            my_execution.execute()
        my_analysis.load_data(my_linac, my_chicane, os.path.split(my_linac.lattice_filename)[0]+"/plots")
        my_analysis.do_plots()
        print()


if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to end")
