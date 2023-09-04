import copy
import os
import matplotlib
import g4bl_longitudinal
import xboa.hit
import xboa.bunch

def single_harmonic_rf():
    n_cells = 100
    cavities = [
        {
            "name":f"pillbox_{i}.1",
            "inner_length":250.0,
            "frequency":0.010,
            "max_gradient":3.0,
            "phase":10.0,
            "z_position":-150.0+1000.0*i,
        } for i in range(1, n_cells)
    ]+[{
            "name":f"pillbox_{i}.2",
            "inner_length":250.0,
            "frequency":0.010,
            "max_gradient":3.0,
            "phase":10.0,
            "z_position":150.0+1000.0*i,
        } for i in range(1, n_cells)
    ]
    return cavities

def dual_harmonic_rf():
    n_cells = 100
    cavities = [{
            "name":f"pillbox_{i}.1",
            "inner_length":250.0,
            "frequency":0.010,
            "max_gradient":3.0,
            "phase":10.0,
            "z_position":-300.0+1000.0*i,
        } for i in range(1, n_cells)
    ]+[{
            "name":f"pillbox_{i}.2",
            "inner_length":250.0,
            "frequency":0.010,
            "max_gradient":3.0,
            "phase":10.0,
            "z_position":300.0+1000.0*i,
        } for i in range(1, n_cells)
    ]+[{
            "name":f"pillbox_{i}.1",
            "inner_length":125.0,
            "frequency":0.020,
            "max_gradient":-1.8,
            "phase":0.0,
            "z_position":1000.0*i-75,
        } for i in range(1, n_cells)
    ]
    return cavities

def multiharmonic_rf(n_cells, principle_harmonic, v1, v2, v3, v4, phi1, phi2, phi3, phi4):
    cavities = [{
            "name":f"pillbox_{i}.1",
            "inner_length":200.0,
            "frequency":principle_harmonic*1,
            "max_gradient":v1,
            "phase":phi1,
            "z_position":150+1000.0*i,
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.2",
            "inner_length":200.0,
            "frequency":principle_harmonic*2,
            "max_gradient":v2,
            "phase":phi2,
            "z_position":400+1000.0*i,
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.3",
            "inner_length":200.0,
            "frequency":principle_harmonic*3,
            "max_gradient":v3,
            "phase":phi3,
            "z_position":650+1000.0*i,
        } for i in range(n_cells)
    ]+[{
            "name":f"pillbox_{i}.4",
            "inner_length":200.0,
            "frequency":principle_harmonic*4,
            "max_gradient":v4,
            "phase":phi4,
            "z_position":900.0+1000.0*i,
        } for i in range(n_cells)
    ]
    return cavities


def build_final_cooling_lattice(cleanup = True):
    lattice_filename = "output/multiharmonic_3/linac.g4bl"
    #cavities = multiharmonic_rf(130, 0.010,
    #                            3.0,   0.0, 0.0, 0.0,
    #                            10.0,  0.0, 0.0, 0.0)
    #cavities = multiharmonic_rf(130, 0.010,
    #                            3.0,  -0.4, 0.0, 0.0,
    #                            10.0,  0.0, 0.0, 0.0)
    cavities = multiharmonic_rf(130, 0.010,
                                3.0,  -0.4, 0.0, 0.12,
                                10.0,  0.0, 0.0, 0.0)

    beam_def = {
        "filename":"beam.txt",
        "out_dir":os.path.split(lattice_filename)[0],
        "beams":[{
            "type":"longitudinal_grid",
            "t_min":-50.0, # ns
            "t_max":50.0, # ns
            "n_t_steps":1,
            "e_min":26,
            "e_max":70,
            "n_e_steps":21,
        }],
    }
    mass = xboa.common.pdg_pid_to_mass[13]
    hit = xboa.hit.Hit.new_from_dict({"pid":-13, "energy":26+mass, "mass":mass}, "pz")
    reference = {
        "p_start":hit["pz"],
    }
    my_linac = g4bl_longitudinal.G4BLLinac(lattice_filename)
    my_linac.cleanup_dir = cleanup
    my_linac.max_step = 1.0
    my_linac.rf_cavities = cavities
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

    def load_data(self):
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)

    def do_plots(self):
        g4bl_longitudinal.clean_dir(self.plot_dir, self.clean_dir)
        self.plot_time_energy_station()
        self.plot_time_energy_event()

    def get_time_energy(self, bunch):
        t_list = bunch.list_get_hit_variable(["t"], ["ns"])[0]
        e_list = bunch.list_get_hit_variable(["kinetic_energy"], ["ns"])[0]
        t_list = [t-t_list[0] for t in t_list]
        return t_list, e_list

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
        figure.savefig(os.path.join(self.plot_dir, "time_energy.png"))

    def plot_time_energy_event(self):
        track_dict = {}
        bunch_list = copy.deepcopy(self.bunch_list)
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
            axes.scatter(t_list, e_list, s=1)
        #dedt = (track_dict[2][-1]['energy']/track_dict[2][-1]['t'])
        #t_lim = axes.get_xlim()
        #axes.plot([0, t_lim[1]], [0, dedt*t_lim[1]], color="grey", linestyle="dashed")
        #axes.set_xlim(t_lim)
        figure.savefig(os.path.join(self.plot_dir, "contours.png"))


def main():
    my_linac = build_final_cooling_lattice(True)
    my_execution = g4bl_longitudinal.G4BLExecution(my_linac)
    my_execution.execute()
    my_analysis = Analysis(my_linac, os.path.split(my_linac.lattice_filename)[0]+"/plots")
    my_analysis.load_data()
    my_analysis.do_plots()



if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to end")
