import json
import glob
import matplotlib
import matplotlib.pyplot

class Plot():
    def __init__(self):
        self.my_data = []
        self.data_glob = ""
        self.filename_keys_depth = 2
        self.plot_dir = ""

    def load_data(self):
        self.my_data = []
        for filename in sorted(glob.glob(self.data_glob)):
            fin = open(filename)
            my_data = json.loads(fin.read())
            my_data["filename"] = filename
            self.parse_filename(my_data, self.filename_keys_depth)
            self.my_data.append(my_data)

    def parse_filename(self, data, target_element):
        filename = data["filename"]
        my_name = filename.split("/")[target_element]
        my_kvs = my_name.split(";_")
        for item in my_kvs:
            key = item.split("=")[0]
            value = float(item.split("=")[1])
            data[key] = value


    def plot(self):
        sub_data = [item for item in self.my_data if item["t_max"] - item["t_min"] > 2.5]
        sub_data = sorted(sub_data, key = lambda x: x["t_min"]+x["t_max"])
        x_list = [(item["t_min"]+item["t_max"])/2.0 for item in sub_data]
        y_list = [item["fractional_yield_of_mu+_per_proton_on_target"] for item in sub_data]
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        axes.plot(x_list, y_list)
        axes.set_xlabel("$\\delta$t [ns]")
        axes.set_ylabel("Number of mu+ per proton on target")
        axes.set_ylim(0.0, 0.10)
        figure.savefig(f"{self.plot_dir}/mu+_yield_per_pot.png")


def main():
    plotter = Plot()
    plotter.plot_dir = "output/rf_capture_v27/"
    plotter.data_glob = plotter.plot_dir+"/t_min=*;_t_max=*/performance.json"
    plotter.load_data()
    plotter.plot()

if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to finish")