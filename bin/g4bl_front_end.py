import os

import g4bl_interface.g4bl_interface
import front_end.g4bl_target_solenoid
import front_end.g4bl_chicane
import xboa
import matplotlib
import numpy

class FrontEnd:
    def __init__(self):
        self.chicane_bfield = 1.5
        self.chicane_bend_angle = 8.3
        self.chicane_r_curv = 41.2
        self.clear_dir = False
        self.n_protons = 100000
        self.run_dir = "output/front_end_v2/1"

    def build_front_end(self):
        self.target_linac, target = front_end.g4bl_target_solenoid.build_target("/tmp/", self.n_protons, self.clear_dir)
        self.chicane_linac, chicane = front_end.g4bl_chicane.build_chicane("/tmp/", self.chicane_bfield, self.chicane_bend_angle, self.chicane_r_curv*1e3, False)
        x_offset = self.chicane_linac.elements[0]["x_position"]
        z_offset = max([ele["z_position"]+ele["length"]/2.0 for ele in self.target_linac.elements if ele["type"] == "solenoid"])
        print(f"Offseting target by x: {x_offset} and chicane by z: {z_offset}")
        self.target_linac.offset_linac(x_offset, 0.0)
        self.chicane_linac.offset_linac(0.0, z_offset)
        self.front_end_linac = g4bl_interface.g4bl_interface.G4BLLinac(os.path.join(self.run_dir, "front_end.g4bl"))
        self.front_end_linac.cleanup_dir = self.clear_dir
        self.front_end_linac.elements = self.target_linac.elements + self.chicane_linac.elements
        self.front_end_linac.max_z = max([ele["z_position"] for ele in self.target_linac.elements])
        self.front_end_linac.beam = self.target_linac.beam
        self.front_end_linac.beam["out_dir"] = self.run_dir
        self.front_end_linac.beam["beams"][0]["n_e_steps"] = self.n_protons
        self.front_end_linac.beam["x_position"] = x_offset
        self.front_end_linac.build_linac()

class Analysis:
    def __init__(self):
        pass

    def load_data(self, linac, plot_dir):
        self.linac = linac
        self.out_dir = linac.out_dir()
        self.plot_dir = os.path.join(linac.out_dir(), plot_dir)
        self.out_filename = os.path.join(self.out_dir, linac.output_file)+".txt"
        self.bunch_list = xboa.bunch.Bunch.new_list_from_read_builtin("icool_for009", self.out_filename)

    def plots(self):
        self.plot_2d("x", "weight", [0.0, 25000.], [600.0, 1000.0])

    def plot_2d(self, var_y, var_z, range_x, range_y):
        bunch_list = [bunch for bunch in self.bunch_list if bunch[0]["z"] > range_x[0] and bunch[0]["z"] < range_x[1]]
        bins_x = [bunch[0]["z"] for bunch in bunch_list]
        bins_x = [(bins_x[i]+z1)/2.0 for i, z1 in enumerate(bins_x[1:])]
        bins_x = [bins_x[0]-(bins_x[1]-bins_x[0])]+bins_x+[bins_x[-1]+(bins_x[-1]-bins_x[-2])]
        x_list, y_list, z_list = [], [], []
        for bunch in bunch_list:
            for hit in bunch:
                x_list.append(hit["z"])
                y_list.append(hit[var_y])
                z_list.append(hit[var_z])
        ny_bins = (len(bunch_list[0])/10)**0.5
        bins_y = numpy.linspace(range_y[0], range_y[1], int(ny_bins))
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        hist2d = axes.hist2d(x_list, y_list, weights=z_list, bins=[bins_x, bins_y])

def main():
    do_execute = True
    front_end = FrontEnd()
    front_end.clear_dir = do_execute
    front_end.build_front_end()
    if do_execute: # alternatively we just redo the analysis
        my_execution = g4bl_interface.g4bl_interface.G4BLExecution(front_end.front_end_linac)
        my_execution.execute()
    analysis = Analysis()
    analysis.load_data(front_end.front_end_linac, "plots")
    analysis.plots()

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")

