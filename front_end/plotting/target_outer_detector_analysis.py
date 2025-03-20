import math
import operator

import matplotlib
import numpy

import xboa

class Analysis:
    def __init__(self):
        self.momentum_cutoff = 1000.0
        self.linac = None
        self.filename = None
        self.plot_dir = None
        self.bunch = None
        self.title = ""

    def load_data(self):
        self.bunch = xboa.bunch.Bunch.new_from_read_builtin("g4beamline_bl_track_file_2", self.filename)
        print("Loaded bunch of length", len(self.bunch))
        if self.momentum_cutoff:
            self.bunch.cut({"p":self.momentum_cutoff}, operator.lt)

    def do_plots(self):
        a_h2d = self.plot_energy_deposit(0.0, 10, 36, 0.0, 100.0, 100, "outer_detector_energy_deposit.png")

    def plot_energy_deposit(self, phi0, dphi, nphi, z0, dz, nz, fname):
        phi_list = [math.degrees(hit["phi"])%360 for hit in self.bunch if hit["p"] > self.momentum_cutoff]
        z_list = [hit["z"] for hit in self.bunch if hit["p"] > self.momentum_cutoff]
        weight_list = [hit["energy"] for hit in self.bunch if hit["p"]  > self.momentum_cutoff]

        figure = matplotlib.pyplot.figure()
        z_bins = [z0 + dz*i for i in range(0, nz+1)]
        phi_bins = [phi0 + dphi*i for i in range(0, nphi+1)]
        axes = figure.add_subplot(1, 1, 1)
        a_h2d = axes.hist2d(z_list, phi_list, weights=weight_list, bins=[z_bins, phi_bins])
        axes.text(0.01, 1.01, self.title, transform=axes.transAxes)

        axes.set_xlabel("z [mm]")
        axes.set_ylabel("$\\Phi$ [degree]")
        figure.savefig(f"{self.plot_dir}/{fname}")
        return a_h2d

def main():
    analysis = Analysis()

if __name__ == "__main__":
    main()

