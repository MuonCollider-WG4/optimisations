import json
import os
import glob
import matplotlib
import matplotlib.pyplot
import xboa.bunch
import numpy

class Analysis:
    def __init__(self):
        self.bunch = xboa.bunch.Bunch()
        self.filename = None
        self.plot_dir = None
        self.z_list = None
        self.e_list = None

    def load_data(self, filename_glob):
        filename_list = sorted(glob.glob(filename_glob))
        for filename in filename_list:
            self.filename = filename
            bunch = xboa.bunch.Bunch.new_from_read_builtin("g4beamline_bl_track_file_2", filename)
            self.bunch = xboa.bunch.Bunch.new_from_hits(self.bunch.hits() + bunch.hits())
            print(f"Loaded {len(bunch)} hits from {filename}. Now have {len(self.bunch)}")
        self.parse_data()

    def parse_data(self):
        self.z_list = []
        self.e_list = []
        for hit in self.bunch:
            if hit["pid"] == -13:
                if hit["kinetic_energy"] < 4.12:
                    self.z_list.append(hit["z"])
                    self.e_list.append(hit["kinetic_energy"])

    def plot_data(self):
        if not os.path.exists(self.plot_dir):
            os.makedirs(self.plot_dir)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.hist(self.z_list, bins=[i*10 for i in range(51)])
        axes.set_xlabel("z [mm]")
        figure.savefig(self.plot_dir+"/z.png")

        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.set_xlim([0.0, 5.0])
        axes.hist(self.e_list, bins=[0.2*i for i in range(26)])
        axes.set_xlabel("Energy [MeV]")
        figure.savefig(self.plot_dir+"/e.png")

    def write_data(self, a_data, filehandle):
        a_data["e_list"] = self.e_list
        a_data["z_list"] = self.z_list
        filehandle.write(json.dumps(a_data)+"\n")

def plot_bulk_data(filename):
    plot_dir = os.path.split(filename)[0]
    fin = open(filename)
    my_data = []
    for line in fin:
        if line == "":
            break
        my_data.append(json.loads(line))

    x_bins = [0.2*i for i in range(0, 23)]
    z_bins = [50*i for i in range(0, 21)]
    dy = (my_data[1]["energy"]-my_data[0]["energy"])
    y_bins = [item["energy"] for item in my_data]
    y_bins = sorted(list(set(y_bins)))
    print(y_bins)
    dy = y_bins[1]-y_bins[0]
    y_bins = [y-dy/2 for y in y_bins]
    y_bins.append(y_bins[-1] + dy)
    print(y_bins)
    #def plot_hist2d_ratio(axes, h1, h2, xedges, yedges, on_nan=0.0, hist_args={}):

    for radius in [10, 20]:
        title = f"target radius = {radius} mm"
        z_figure = matplotlib.pyplot.figure()
        z_axes = z_figure.add_subplot(1, 1, 1)
        e_figure = matplotlib.pyplot.figure()
        e_axes = e_figure.add_subplot(1, 1, 1)
        e_list_ratio = []
        e_list = []
        z_list = []
        proton_energy = []
        for item in my_data:
            if item["radius"] != radius:
                continue
            print("Loaded item energy", item["energy"], item["radius"])
            e_list += item["e_list"]
            z_list += item["z_list"]
            proton_energy += [item["energy"] for e in item["e_list"]]
        e_h = e_axes.hist2d(e_list, proton_energy, bins=[x_bins, y_bins])
        e_axes.set_xlabel("Muon energy [MeV]")
        e_axes.set_ylabel("Proton energy [MeV]")
        e_axes.set_title(title)
        e_figure.colorbar(e_h[3])
        e_figure.savefig(f"{plot_dir}/yield_vs_energy_r={radius}.png")
        z_h = z_axes.hist2d(z_list, proton_energy, bins=[z_bins, y_bins])
        z_axes.set_xlabel("Target position [mm]")
        z_axes.set_ylabel("Proton energy [MeV]")
        z_axes.set_title(title)
        z_figure.colorbar(z_h[3])
        z_figure.savefig(f"{plot_dir}/yield_vs_z_r={radius}.png")


def single_plots(run_dir):
    fout = open(f"{run_dir}/data_out.txt", "w")
    for energy in [i*100 for i in range(3, 21)]:
        for radius in [10, 20]:
            my_analysis = Analysis()
            my_analysis.plot_dir = f"{run_dir}/plot_energy={energy}_radius={radius}/"
            file_glob = f"{run_dir}/energy={energy};_version=*;_radius={radius}/target_outer.txt"
            for filename in sorted(glob.glob(file_glob)):
                my_analysis.load_data(filename)
            my_analysis.plot_data()
            a_data = {"energy":energy, "radius":radius, "file_glob":file_glob}
            my_analysis.write_data(a_data, fout)
    fout.close()

def main():
    run_dir = "output/isis_target_model_v9/"
    #single_plots(run_dir)
    plot_bulk_data(f"{run_dir}/data_out.txt")


if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")