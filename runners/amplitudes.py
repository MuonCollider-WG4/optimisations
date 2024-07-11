import json
import sys
import shutil
import os
import math
import numpy
import matplotlib.pyplot
import scipy.stats
import runners.evolve
try:
    import ROOT
except ImportError:
    raise ImportError("PyROOT import failed - you can still run the code, but minuit fitter won't work")
from models.field_models import FieldSum
from models.field_models import FlatGaussField
from models.field_models import LinearInterpolation
from models.field_models import SineField
from models.field_models import UniformField
from models.field_models import CurrentSheet
from models.field_models import CurrentBlock


"""
Plot the transverse amplitude distribution and compare with physical aperture
"""

class AmplitudeDistribution:
    def __init__(self):
        self.field = None
        self.rms_emittance = 0.0
        self.cavity_frequency = -1.0
        self.iris_factor = 0.5
        self.p_list = []
        figure = matplotlib.pyplot.figure()
        self.axes_pdf = figure.add_subplot(1, 1, 1)
        self.axes_aperture = self.axes_pdf.twinx()

    def get_beta_max(self, pz):
        self.finder = runners.evolve.BetaFinder(self.field, pz)
        beta, alpha, phase_advance, beta_mid = self.finder.get_beta_periodic()
        if beta == 0.0:
            return -1
        print("Found periodic pz, beta", pz, beta)
        z_list, beta_list, dbetadz_list, phi_list = \
                                    self.finder.propagate_beta(beta, 0.)
        beta_max = max(beta_list)
        return beta_max

    def get_amplitude_max(self, pz):
        physical_aperture = self.ideal_cavity_aperture()
        beta_max = self.get_beta_max(pz)
        if beta_max == 0.0:
            return 0.0
        amplitude = physical_aperture**2/beta_max
        return amplitude

    def ideal_cavity_aperture(self):
        # frequency in GHz
        c_light = 0.299 # m/ns
        cavity_radius = 2.40483*c_light/2/math.pi/self.cavity_frequency
        iris_radius = cavity_radius*self.iris_factor
        return iris_radius

    def plot_amplitude_distribution(self):
        a_list = [i/1000 for i in range(201)]
        degrees_of_freedom = 4
        n_list = [scipy.stats.chi2.pdf(a/self.rms_emittance, degrees_of_freedom) for a in a_list]
        c_list = [scipy.stats.chi2.cdf(a/self.rms_emittance, degrees_of_freedom) for a in a_list]
        self.axes_pdf.plot(a_list, n_list)
        self.axes_pdf.set_xlabel("Amplitude [m]")
        self.axes_pdf.set_ylabel("pdf")

    def plot_aperture(self, label):
        aperture_list = [self.get_amplitude_max(pz) for pz in self.p_list]
        self.axes_aperture.plot(aperture_list, self.p_list, "--", label=label)
        self.axes_aperture.set_ylabel("momentum [GeV/c]")

def main():
    #clear_dir(plot_dir)
    amplitude_distribution = AmplitudeDistribution()
    amplitude_distribution.rms_emittance = 0.0176
    amplitude_distribution.p_list = [0.001*i for i in range(150, 401, 10)]
    amplitude_distribution.plot_amplitude_distribution()

    amplitude_distribution.cavity_frequency = 0.352/2


    for b0, b1, length, label in [
                        (2.0, 0.0, 1.0, "2 T solenoid"),
                        (0.0, 2.5, 1.8, "ZR A1"),
                        (0.0, 3.0, 1.8, "A-type"),]:
        solenoid_field = SineField(b0, b1, 0.0, 0.0, 0.0, 0.0, length)
        amplitude_distribution.field = solenoid_field
        amplitude_distribution.plot_aperture(label+f" L: {length} b0: {b0} b1: {b1}")
    amplitude_distribution.axes_aperture.legend()


if __name__ == "__main__":
    main()
    matplotlib.pyplot.show(block=False)
    input("Press <CR> to end")
