import glob
import bisect
import json
import matplotlib
import matplotlib.pyplot

class PostPlot():
	def __init__(self):
		self.data = []
		self.pz_target = 0.2 # GeV/c
		self.harmonic_target = 1
		self.plot_dir = "./"

	def get_beta_at_pz_target(self, item):
		index = bisect.bisect_left(item["pz_list"], self.pz_target)
		x_item = item["beta_list"][index]
		return x_item

	def get_bi_harmonic(self, item):
		return item["field_bi"][self.harmonic_target]

	def load_data(self, file_glob):
		self.data = []
		print(f"Loading... {file_glob}")
		for filename in sorted(glob.glob(file_glob)):
			fin = open(filename)
			self.data.append(json.loads(fin.read()))
			print(f"    {filename}")

	def plot(self, x_lambda, y_lambda, x_label, y_label):
		x_list = []
		y_list = []
		for item in self.data:
			x_list.append(x_lambda(item))
			y_list.append(y_lambda(item))
		print(x_list)
		print(y_list)
		zipped = sorted(zip(x_list, y_list))
		zipped = zip(*zipped)
		zipped = [x for x in zipped]
		x_list, y_list = zipped
		figure = matplotlib.pyplot.figure()
		axes = figure.add_subplot(1, 1, 1)
		axes.plot(x_list, y_list)
		axes.set_xlabel(x_label)
		axes.set_ylabel(y_label)
		figure.savefig(f"{self.plot_dir}/figure_1.png")

def main():
	plotter = PostPlot()
	plotter.plot_dir = "optics-scan_v21"
	plotter.harmonic_target = 0
	plotter.pz_target = 0.2
	plotter.load_data(f"{plotter.plot_dir}/*.json")
	plotter.plot(plotter.get_bi_harmonic, plotter.get_beta_at_pz_target, f"$b_{{{plotter.harmonic_target}}}$ [T]", f"$\\beta_{{\\perp}}$ [mm] at {plotter.pz_target} GeV/c")

#NEXT JOB DO b3 and then go back and do b-1

if __name__ == "__main__":
	main()
