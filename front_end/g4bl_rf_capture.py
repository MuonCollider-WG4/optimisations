import sys
import os
import math
import matplotlib

import xboa.common
import g4bl_interface.g4bl_interface


class RFCapture:
    def __init__(self):
        # lengths
        self.ref_p0 = 200.0 # MeV/c
        self.ref_p1 = 400.0 # MeV/c
        self.ref_t0 = 0.0 # ns
        self.ref_t1 = 0.0 # ns
        self.z_pos = 0.0
        self.line_param = []

        self.rf_capture = [{
            "name":"default",
            "n_rf_repeats":10, # integer
            "rf_cell_length":200.0, # mm
            "rf_repeat_length":250.0, # mm
            "rf_phi0":0.0, # degree
            "rf_phi1":0.0, # degree
            "rf_e_z0":0.0, # MV/m
            "rf_e_z1":1.0, # MV/m
            "n_bunches":1, # integer
            "bz":1.5, # kT
        }]
        self.default_coil = {
            "x_position":0.0,
            "y_rotation":0.0, # rotation about the y axis i.e. in x z plane
            "z_position":0.0,
            "inner_radius":200.0,
            "outer_radius":300.0,
            "length":100.0,
            "current":1,
            "name":"coil",
            "type":"solenoid",
        }

    def build(self):
        element_list = []
        for line in self.rf_capture:
            element_list += self.build_line(line)
        return element_list

    def kick_momentum(self, momentum, phi, length, gradient, frequency):
        energy = (momentum**2+self.mass**2)**0.5
        beta = momentum/energy
        ttf_arg = math.pi*length*frequency/beta/xboa.common.constants["c_light"]
        ttf = math.sin(ttf_arg)/ttf_arg
        print("TTF:", ttf, "L:", length, "f:", frequency, "beta:", beta, "P:", momentum)
        energy += ttf*gradient*length*1e-3*math.sin(math.radians(phi))
        if type(momentum) != type(0.0) or energy < 0:
            raise RuntimeError(f"negative energy p {momentum} phi {phi} l {length} v {gradient}")
        momentum = (energy**2-self.mass**2)**0.5
        return momentum

    def get_delta_time(self, momentum, distance):
        energy = (momentum**2+self.mass**2)**0.5
        beta = momentum/energy
        dt = distance/(beta*xboa.common.constants["c_light"])
        return dt

    def get_gradient(self, i, line):
        if line["rf_e_z0"] is None or line["rf_e_z1"] is None:
            return None
        if line["n_rf_repeats"] == 1:
            return (line["rf_e_z0"]+line["rf_e_z1"])/2
        return line["rf_e_z0"] + (line["rf_e_z1"]-line["rf_e_z0"])*i/(line["n_rf_repeats"]-1)

    def get_frequency(self, t0, t1, phi0, phi1, n_bunches):
        old_frequency = 0.0
        if abs(t1 - t0) < 1e-3:
            return 1e3
        frequency = abs(n_bunches/(t1 - t0))
        print("Setting frequency")
        while (frequency - old_frequency)/frequency > 1e-6:
            old_frequency = frequency
            delta0 = phi0/(360*frequency)
            delta1 = phi1/(360*frequency)
            #BUG - WIP
            frequency = abs(n_bunches/(t1 - delta1 - t0 + delta0))
            print("    ", frequency, t1, delta1, t0, delta0)
        return frequency

    def build_line(self, line):
        element_list = []
        sol_length = line['n_rf_repeats']*line['rf_repeat_length']
        element_list.append({
            "name":f"{line['name']}_solenoid",
            "radius":line['radius'],
            "bz":line['bz'],
            "z_position":self.z_pos+sol_length/2,
            "length":sol_length,
            "x_position":0.0,
            "type":"uniform_field",
        })
        for i in range(line["n_rf_repeats"]):
            rf_cavity = g4bl_interface.g4bl_interface.Cavity()
            new_z_pos = self.z_pos+line["rf_repeat_length"]
            self.ref_t0 += self.get_delta_time(self.ref_p0, new_z_pos - self.z_pos)
            self.ref_t1 += self.get_delta_time(self.ref_p1, new_z_pos - self.z_pos)
            virtual_n_bunches = line["n_bunches"]+(line["rf_phi1"]-line["rf_phi0"])/360.0
            frequency = 1/abs((self.ref_t1-self.ref_t0)/virtual_n_bunches)
            print("Old frequency", frequency)
            frequency = self.get_frequency(self.ref_t0, self.ref_t1, line["rf_phi0"], line["rf_phi1"], line["n_bunches"])
            print("New frequency", frequency)
            gradient = self.get_gradient(i, line)
            rf_cavity = {
                "name":f"{line['name']}_rf_{i}",
                "inner_length":line["rf_cell_length"],
                "frequency":frequency,
                "max_gradient":gradient,
                "x_position":0.0,
                "z_position":new_z_pos,
                "rgb":"0.75,0.5,0.0",
                "phase":line["rf_phi0"],
                "type":"cavity",
            }
            if gradient is not None:
                element_list.append(rf_cavity)
            if "static_energy_change" in line:
                dE = line["static_energy_change"]
                length = line["rf_repeat_length"]-line["rf_cell_length"]
                stopper = {
                    "name":f"{line['name']}_stopper_{i}",
                    "inner_length":length,
                    "frequency":1e-12,
                    "max_gradient":dE/length*1e3,
                    "x_position":0.0,
                    "z_position":new_z_pos+line["rf_repeat_length"]/2,
                    "rgb":"0.75,0.5,0.0",
                    "phase":90.0,
                    "type":"cavity",
                }
                element_list.append(stopper)
            else:
                stopper = None

            self.z_pos = new_z_pos
            self.line_param.append({
                "name":line['name'],
                "rf_cavity":rf_cavity,
                "stopper":stopper,
                "z_pos":new_z_pos,
                "ref_t0":1.0*self.ref_t0,
                "ref_t1":1.0*self.ref_t1,
                "ref_p0":1.0*self.ref_p0,
                "ref_p1":1.0*self.ref_p1,
                "length":line["rf_repeat_length"],
            })
            if gradient is not None:
                self.ref_p0 = self.kick_momentum(self.ref_p0, line["rf_phi0"], rf_cavity["inner_length"], rf_cavity["max_gradient"], rf_cavity["frequency"])
                self.ref_p1 = self.kick_momentum(self.ref_p1, line["rf_phi1"], rf_cavity["inner_length"], rf_cavity["max_gradient"], rf_cavity["frequency"])
            print("z", new_z_pos, "t0, p0", self.ref_t0, self.ref_p0, "t1, p1", self.ref_t1, self.ref_p1, "freq", rf_cavity["frequency"])
        return element_list

    mass = xboa.common.pdg_pid_to_mass[13]

def get_beam(my_linac, rf_capture, n_e_steps, n_t_steps):
    x_start = my_linac.elements[0]["x_position"]
    print(f"Offsetting beam to x={x_start} mm")
    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(my_linac.lattice_filename)[0],

        "beams":[{
            "type":"longitudinal_grid",
            "t_min":-5.0, # ns
            "t_max":5.0, # ns
            "n_t_steps":n_t_steps,
            "e_min":50,
            "e_max":400,
            "n_e_steps":n_e_steps,
            "pid":-13,
            "default_hit":{"x":x_start, "mass":xboa.common.pdg_pid_to_mass[13], "pid":-13},
        }],
    }
    reference = [{
        "p_start":rf_capture.line_param[0]["ref_p0"],
    },{
        "p_start":rf_capture.line_param[0]["ref_p1"],
    }]

    my_linac.beam = beam_def
    my_linac.reference = reference

class Analysis:
    def __init__(self, base_dir):
        self.base_dir = base_dir
        self.target_ke = (300**2+105.658**2)**0.5-105.658

    def load_data(self, my_linac, rf_capture, plot_dir):
        self.out_dir = my_linac.out_dir()
        self.plot_dir = plot_dir
        self.out_filename = os.path.join(self.out_dir, my_linac.output_file)+".txt"
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)
        self.rf_capture = rf_capture

    def t_periodicity(self, t_list, frequency):
        dt_list = [None]*len(t_list)
        for i, t in enumerate(t_list):
            nu = t*frequency
            dt_list[i] = (nu - math.floor(nu+0.5))/frequency # domain between -0.5/f, 0.5/f
        return dt_list

    def get_trajectory(self, t, p, cell):
        t_list, p_list = [t], [p]
        cavity = cell["rf_cavity"]
        omega = 360*cavity["frequency"]
        for i in range(100):
            delta_z = cell["length"]
            t += self.rf_capture.get_delta_time(p, delta_z)
            phi = cavity["phase"]+omega*t
            print(phi, cavity["max_gradient"], cavity["inner_length"], end = " ")
            p = self.rf_capture.kick_momentum(p, phi, cavity["inner_length"], cavity["max_gradient"], cavity["frequency"])
            t_list.append(t)
            p_list.append(p)
        print("\n    p", p_list)
        mass = self.rf_capture.mass
        ke_list = [(p**2+mass**2)**0.5-mass for p in p_list]
        return t_list, ke_list

    def plot_contours(self, cell, axes):
        ref_p0 = cell["ref_p0"]
        freq = cell["rf_cavity"]["frequency"]
        ref_t_list, ref_e_list = self.get_trajectory(0.0, ref_p0, cell)
        for ti in [-0.5, 0.5]:
            print("ti", ti)
            t = ti/freq
            t_list, e_list = self.get_trajectory(t, ref_p0, cell)
            dt_list = [t-ref_t_list[i] for i, t in enumerate(t_list)]
            de_list = [e-ref_e_list[i]+ref_e_list[0] for i, e in enumerate(e_list)]
            axes.plot(dt_list, de_list)
            print("    dt", dt_list)
            print("    de", de_list)

    def do_plots(self):
        t_range = []
        e_range = [0]
        for bunch in self.bunch_list:
            t_list = [hit["t"] - bunch[0]["t"] for hit in bunch]+t_range
            e_list = [hit["energy"]-hit["mass"] for hit in bunch]+e_range
            t_range = [min(t_list), max(t_list)]
            e_range = [min(e_list), max(e_list)]
        t_range = [t_range[0]-(t_range[1]-t_range[0])/10.0, t_range[1]+(t_range[1]-t_range[0])/10.0]
        e_range = [e_range[0]-(e_range[1]-e_range[0])/10.0, e_range[1]+(e_range[1]-e_range[0])/10.0]
        max_period = max([1/item['rf_cavity']['frequency'] for item in self.rf_capture.line_param])
        dt_range = [-max_period, max_period]
        mass = self.rf_capture.mass

        figure_dt = matplotlib.pyplot.figure("dt vs E")
        axes_dt = figure_dt.add_subplot(1, 1, 1)

        figure = matplotlib.pyplot.figure("t vs E")
        axes = figure.add_subplot(1, 1, 1)
        for bunch in self.bunch_list:
            z = round(bunch[0]["z"])
            cell = self.rf_capture.line_param[0]
            for test_cell in self.rf_capture.line_param:
                if test_cell["z_pos"] > z:
                    break
                cell = test_cell
            print(f"\rplotting {cell['z_pos']:10.6g}  {cell['ref_t0']:8.4g}  {cell['ref_t1']:8.4g} p0 {bunch[0]['p']:8.4g} p1  {bunch[1]['p']:8.4g}", end="            ")
            axes.clear()
            axes_dt.clear()
            t_list = [hit["t"] - bunch[0]["t"] for hit in bunch]
            e_list = [hit["energy"]-hit["mass"] for hit in bunch]
            axes.scatter(t_list[2:], e_list[2:], s=1)
            axes.scatter(t_list[:2], e_list[:2], s=2, color="red")

            dt_list = self.t_periodicity(t_list, cell['rf_cavity']['frequency'])
            axes_dt.scatter(dt_list[2:], e_list[2:], s=1)
            axes_dt.scatter(dt_list[:2], e_list[:2], s=2, color="red")


            t_list = [cell["ref_t0"] - bunch[0]["t"], cell["ref_t1"] - bunch[0]["t"]]
            dt_list = self.t_periodicity(t_list, cell['rf_cavity']['frequency'])
            e_list = [cell["ref_p0"], cell["ref_p1"]]
            e_list = [(e**2+mass**2)**0.5-mass for e in e_list]
            axes.scatter(t_list, e_list, s=4, edgecolors="blue", facecolors="none")
            axes_dt.scatter(dt_list, e_list, s=4, edgecolors="blue", facecolors="none")
            axes.set_xlabel("t-t$_0$ [ns]")
            axes.set_ylabel("Kinetic Energy [MeV]")
            axes.plot(t_range, [self.target_ke, self.target_ke], linestyle="--")
            axes.set_xlim(t_range)
            axes.set_ylim(e_range)
            axes_dt.set_xlim(dt_range)
            axes_dt.set_ylim(e_range)
            freq = cell['rf_cavity']['frequency']
            volts = cell['rf_cavity']['max_gradient']
            text = f"{cell['name']}\nz={z*1e-3:8.4g} m\nt$_{{ref0}}$={cell['ref_t0']:8.4g}\n"+\
                    f"t$_{{ref1}}$={cell['ref_t1']:8.4g}\n$t_{{0}}$={bunch[0]['t']:8.4g} ns\n"+\
                    f"$t_{{1}}$={bunch[1]['t']:8.4g} ns\nf={freq:8.4g} GHz\nv={volts:8.4g} MV/m"
            axes.text(0.76, 0.55, text, transform=axes.transAxes)
            axes_dt.text(0.76, 0.55, text, transform=axes.transAxes)
            figure.savefig(f"{self.plot_dir}/t_vs_e_z={int(z):0>6}.png")
            figure_dt.savefig(f"{self.plot_dir}/dt_vs_e_z={int(z):0>6}.png")
        print("\n\nDone")

def get_test_capture():
    n_bunches = 1
    bz = 1.5
    rf_capture = [{
        "name":"drift",
        "n_rf_repeats":20, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":0.0, # degree
        "rf_phi1":0.0, # degree
        "rf_e_z0":0.0, # MV/m; None to ignore the RF cavity (just update the position)
        "rf_e_z1":0.0, # MV/m
        "n_bunches":n_bunches, # integer
        "bz":bz, # kT?
        "radius":500.0,
    },{
        "name":"rotator",
        "n_rf_repeats":5, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":12.0, # degree
        "rf_phi1":-12.0, # degree
        "rf_e_z0":12.0, # MV/m
        "rf_e_z1":12.0, # MV/m
        "n_bunches":n_bunches, # integer
        "bz":bz, # kT?
        "radius":500.0,
        "static_energy_change":-1.0,
    }]
    return rf_capture


def get_longitudinal_capture():
    n_bunches = 7
    bz = 1.5
    rf_capture = [{
        "name":"drift",
        "n_rf_repeats":202, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":0.0, # degree
        "rf_phi1":0.0, # degree
        "rf_e_z0":0.0, # MV/m; None to ignore the RF cavity (just update the position)
        "rf_e_z1":0.0, # MV/m
        "n_bunches":n_bunches, # integer
        "bz":bz, # kT?
        "radius":500.0,
    },{
        "name":"buncher",
        "n_rf_repeats":200, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":0.0, # degree
        "rf_phi1":0.0, # degree
        "rf_e_z0":0.0, # MV/m
        "rf_e_z1":12.0, # MV/m
        "n_bunches":n_bunches, # integer
        "bz":bz, # kT?
        "radius":500.0,
    },{
        "name":"rotator",
        "n_rf_repeats":200, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":11.4, # degree
        "rf_phi1":-11.8, # degree
        "rf_e_z0":12.0, # MV/m
        "rf_e_z1":12.0, # MV/m
        "n_bunches":n_bunches, # integer
        "bz":bz, # kT?
        "radius":500.0,
    },{
        "name":"continuation",
        "n_rf_repeats":200, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":0.0, # degree
        "rf_phi1":0.0, # degree
        "rf_e_z0":12.0, # MV/m
        "rf_e_z1":12.0, # MV/m
        "n_bunches":n_bunches, # integer
        "bz":bz, # kT?
        "radius":500.0,
        "static_energy_change":-0.5,
    }]
    return rf_capture

def build_longitudinal_capture(run_dir, do_execute):
    cleanup = do_execute
    lattice_filename = run_dir+"/longitudinal_capture.g4bl"
    g4bl_interface.g4bl_interface.Solenoid.clear()
    my_linac = g4bl_interface.g4bl_interface.G4BLLinac(lattice_filename)
    my_linac.cleanup_dir = cleanup
    my_linac.max_step = 10.0
    my_linac.z_spacing = 250.0
    rf_capture = RFCapture()
    rf_capture.rf_capture = get_longitudinal_capture()
    my_linac.elements = rf_capture.build()
    get_beam(my_linac, rf_capture, 21, 5)# 21, 5)
    my_linac.do_stochastics = 0 # e.g. decays
    my_linac.max_z = max([sol["z_position"] for sol in my_linac.elements])
    my_linac.build_linac()
    return my_linac, rf_capture

def main(do_execute):
    base_dir = f"output/rf_capture_v7/"
    my_analysis = Analysis(base_dir)
    for r_curv, bend_angle, b_field, in [(40, 10.0, 1.0)]:
        run_dir = os.path.join(base_dir)
        my_linac, rf_capture = build_longitudinal_capture(base_dir, do_execute)
        if do_execute: # alternatively we just redo the analysis
            my_execution = g4bl_interface.g4bl_interface.G4BLExecution(my_linac)
            my_execution.execute()
        my_analysis.load_data(my_linac, rf_capture, base_dir)
        my_analysis.do_plots()
        print()

if __name__ == "__main__":
    do_execute = "--only-plot" not in sys.argv
    main(do_execute)
    if "--no-pause" not in sys.argv:
        matplotlib.pyplot.show(block=False)
        input("Press <CR> to end")

