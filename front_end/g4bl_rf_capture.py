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

    def calculate_current(self, coil = None):
        if coil == None:
            coil = self.default_coil
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
        #print("TTF:", ttf, "L:", length, "P:", momentum)
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

    def build_line(self, line):
        element_list = []
        for i in range(line["n_rf_repeats"]):
            rf_cavity = g4bl_interface.g4bl_interface.Cavity()
            new_z_pos = self.z_pos+line["rf_repeat_length"]
            self.ref_t0 += self.get_delta_time(self.ref_p0, new_z_pos - self.z_pos)
            self.ref_t1 += self.get_delta_time(self.ref_p1, new_z_pos - self.z_pos)
            virtual_n_bunches = line["n_bunches"]+(line["rf_phi1"]-line["rf_phi0"])/360.0
            frequency = 1/abs((self.ref_t1-self.ref_t0)/virtual_n_bunches)
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
            self.z_pos = new_z_pos
            self.line_param.append({
                "rf_cavity":rf_cavity,
                "z_pos":new_z_pos,
                "ref_t0":1.0*self.ref_t0,
                "ref_t1":1.0*self.ref_t1,
                "ref_p0":1.0*self.ref_p0,
                "ref_p1":1.0*self.ref_p1,

            })
            if gradient is not None:
                self.ref_p0 = self.kick_momentum(self.ref_p0, line["rf_phi0"], rf_cavity["inner_length"], rf_cavity["max_gradient"], rf_cavity["frequency"])
                self.ref_p1 = self.kick_momentum(self.ref_p1, line["rf_phi1"], rf_cavity["inner_length"], rf_cavity["max_gradient"], rf_cavity["frequency"])
            print("z", new_z_pos, "t0, p0", self.ref_t0, self.ref_p0, "t1, p1", self.ref_t1, self.ref_p1)
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

    def load_data(self, my_linac, rf_capture, plot_dir):
        self.out_dir = my_linac.out_dir()
        self.plot_dir = plot_dir
        self.out_filename = os.path.join(self.out_dir, my_linac.output_file)+".txt"
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)
        self.rf_capture = rf_capture

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
        mass = self.rf_capture.mass

        figure = matplotlib.pyplot.figure("t vs E")
        axes = figure.add_subplot(1, 1, 1)
        for bunch in self.bunch_list:
            z = round(bunch[0]["z"])
            cell = self.rf_capture.line_param[0]
            for test_cell in self.rf_capture.line_param:
                print("checking", len(self.rf_capture.line_param), test_cell["z_pos"], z)
                if test_cell["z_pos"] > z:
                    break
                cell = test_cell
            print(cell["z_pos"], cell["ref_t0"], cell["ref_t1"])
            axes.clear()
            t_list = [hit["t"] - bunch[0]["t"] for hit in bunch[2:]]
            e_list = [hit["energy"]-hit["mass"] for hit in bunch[2:]]
            axes.scatter(t_list, e_list, s=1)
            t_list = [hit["t"] - bunch[0]["t"] for hit in bunch[:2]]
            e_list = [hit["energy"]-hit["mass"] for hit in bunch[:2]]
            axes.scatter(t_list, e_list, s=2, color="red")
            t_list = [cell["ref_t0"] - bunch[0]["t"], cell["ref_t1"] - bunch[0]["t"]]
            e_list = [cell["ref_p0"], cell["ref_p1"]]
            e_list = [(e**2+mass**2)**0.5-mass for e in e_list]
            axes.scatter(t_list, e_list, s=4, edgecolors="blue", facecolors="none")
            axes.set_xlabel("t-t$_0$ [ns]")
            axes.set_ylabel("Kinetic Energy [MeV]")
            axes.set_xlim(t_range)
            axes.set_ylim(e_range)
            freq = cell['rf_cavity']['frequency']
            volts = cell['rf_cavity']['max_gradient']
            axes.text(0.78, 0.70, f"z={z*1e-3} m\nt$_0$={bunch[0]['t']:8.4g} ns\nref$_{{t0}}$={cell['ref_t0']:8.4g}\nref$_{{t1}}$={cell['ref_t1']:8.4g}\nf={freq:8.4g} GHz\nv={volts:8.4g} MV/m", transform=axes.transAxes)
            figure.savefig(f"{self.plot_dir}/t_vs_e_z={int(z):0>6}.png")
        print("\n\nDone")

def get_longitudinal_capture():
    rf_capture = [{
        "name":"drift",
        "n_rf_repeats":400, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":0.0, # degree
        "rf_phi1":0.0, # degree
        "rf_e_z0":0.0, # MV/m; None to ignore the RF cavity (just update the position)
        "rf_e_z1":0.0, # MV/m
        "n_bunches":21, # integer
        "bz":1.5, # kT?
    },{
        "name":"buncher",
        "n_rf_repeats":200, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":0.0, # degree
        "rf_phi1":0.0, # degree
        "rf_e_z0":0.0, # MV/m
        "rf_e_z1":10.0, # MV/m
        "n_bunches":21, # integer
        "bz":1.5, # kT?
    },{
        "name":"rotator",
        "n_rf_repeats":200, # integer
        "rf_cell_length":200.0, # mm
        "rf_repeat_length":250.0, # mm
        "rf_phi0":8.0, # degree
        "rf_phi1":-20.0, # degree
        "rf_e_z0":10.0, # MV/m
        "rf_e_z1":10.0, # MV/m
        "n_bunches":21, # integer
        "bz":1.5, # kT?
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
    get_beam(my_linac, rf_capture, 2, 2)
    my_linac.do_stochastics = 0 # e.g. decays
    my_linac.max_z = max([sol["z_position"] for sol in my_linac.elements])
    my_linac.build_linac()
    return my_linac, rf_capture

def main(do_execute):
    base_dir = f"output/rf_capture_v2/"
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

















def sine_hack():
    sine_list = [
        (1.0, 1, 0),
        (1.05, 0.9, 0),
        (1.10, 0.2, 0),
        (1.15, 0.01, 0),
    ]
    t_list = [t/1000 for t in range(20001)]
    v_list = [0 for t in t_list]
    zero_crossing_t_list = []
    zero_crossing_v_list = []
    zero_crossing_v_max = [0]
    for i, t in enumerate(t_list):
        for freq, amp, phase in sine_list:
            v_list[i] += amp*math.sin(freq*2*math.pi*t + phase)
        zero_crossing_v_max[-1] = max(zero_crossing_v_max[-1], abs(v_list[i]))

        if i > 0 and (v_list[i] > 0 and v_list[i-1] < 0):
            t = (t-t_list[i-1])/(v_list[i]-v_list[i-1])*abs(v_list[i-1])+t_list[i-1]
            v = 0
            for freq, amp, phase in sine_list:
                v += amp*math.sin(freq*2*math.pi*t + phase)
            zero_crossing_t_list.append(t)
            zero_crossing_v_list.append(v)
            zero_crossing_v_max.append(0)


    figure = matplotlib.pyplot.figure()
    axes = figure.add_subplot(1, 1, 1)
    axes.plot(t_list, v_list)
    axes.scatter(zero_crossing_t_list, zero_crossing_v_list)
    axes.set_xlabel("t")
    axes.set_ylabel("V")

    dt_list = [t1-zero_crossing_t_list[i] for i, t1 in enumerate(zero_crossing_t_list[1:])]
    figure = matplotlib.pyplot.figure()
    axes = figure.add_subplot(1, 1, 1)
    axes.plot(zero_crossing_t_list[:-1], dt_list)
    axes2 = axes.twinx()
    axes2.plot(zero_crossing_t_list, zero_crossing_v_max[1:], color="red")
    axes.set_xlabel("t")
    axes.set_ylabel("dt", color="blue")
    axes2.set_ylabel("v$_{peak}$", color="red")
