import math
import unittest

import matplotlib
import matplotlib.pyplot

import models.rf_model.pulsed_rf

class TestPulsedRF(unittest.TestCase):
    def setUp(self):
        self.structure = models.rf_model.pulsed_rf.PulsedRF()
        self.structure.set_tanh_field_model(2.0, 100.0, 10.0)
        self.structure.set_max_derivative(4)
        self.structure.set_length(1000.0)
        self.structure.set_max_r(100.0)
        self.structure.set_v0(0.01)
        self.structure.set_z_offset(-400.0)
        self.var_list = ["bx", "by", "bz", "ex", "ey", "ez"]
        self.delta = [0.01, 0.01, 0.01, 0.01] # dx dy dz dt

    def tearDown(self):
        pass

    def _test_construction(self):
        self.assertEqual(self.structure.get_max_derivative(), 4)
        self.assertEqual(self.structure.get_length(), 1000.0)
        self.assertEqual(self.structure.get_max_r(), 100.0)
        self.assertEqual(self.structure.get_v0(), 0.01)
        self.assertEqual(self.structure.get_z_offset(), -400.0)

    def _test_field_on_axis(self):
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        for ti in range(2):
            t = ti * 100
            field_lists = dict([(var, []) for var in self.var_list])
            z_list = []
            for zi in range(-50, 50):
                z = zi * 10
                z_list.append(z)
                field = self.structure.get_field_value(0, 0, z, t)
                for i, value in enumerate(field):
                    field_lists[self.var_list[i]].append(value)
                    if i != 5:
                        self.assertEqual(value, 0.0)
            for var in self.var_list[5:]:
                axes.plot(z_list, field_lists[var], label=f"{var} @ t={t}")
        axes.legend()
        figure.savefig("field_on_axis.png")

    def _test_field_off_axis_zeroes(self):
        for ti in range(2):
            for xi in range(3):
                for yi in range(3):
                    for zi in range(-10, 10):
                        point = (xi*10, yi*10, zi*10-400, ti*10)
                        field = self.structure.get_field_value(*point)
                        for i in [1, 2, 3]:
                            self.assertEqual(field[i], 0.0)

    def _test_field_off_axis_x_indie(self):
        point = [0.0, 10.0, -350, 5.0]
        refField = None
        for xi in range(-9, 10):
            point[0] = xi*10
            field = self.structure.get_field_value(*point)
            if refField != None:
                for i in range(6):
                    self.assertEqual(field[i], refField[i])
            refField = field

    def _test_bb(self):
        self.structure.set_tanh_field_model(2.0, 10000.0, 10.0) # very long, fill the space
        for point, in_bounds in[
                ([0.0, 99.0, -350, 5.0], True),
                ([0.0, 101.0, -350, 5.0], False),
                ([99.0, 0.0, -350, 5.0], True),
                ([101.0, 0.0, -350, 5.0], False),
                ([0.0, 0.0, -499, 5.0], True),
                ([0.0, 0.0, -501, 5.0], False),
                ([0.0, 0.0, 499, 5.0], True),
                ([0.0, 0.0, 501, 5.0], False),
            ]:
            field = self.structure.get_field_value(*point)
            if in_bounds:
                self.assertNotEqual(field[5], 0.0, msg=point)
            else:
                self.assertEqual(field[5], 0.0, msg=point)

    def _test_field_map(self):
        y_step = 1
        z_step = 1
        n_t_steps = 10
        for ti in range(n_t_steps+1):
            figure = matplotlib.pyplot.figure()
            bx_axes = figure.add_subplot(2, 2, 1)
            ey_axes = figure.add_subplot(2, 2, 2)
            ez_axes = figure.add_subplot(2, 2, 3)
            y_bins = []
            y_list = []
            z_list = []
            ez_field_list = []
            ey_field_list = []
            bx_field_list = []
            t = ti/n_t_steps*1000/300/self.structure.get_v0()
            for yi in range(-20, 21):
                z_bins = []
                y = yi * y_step
                y_bins.append(y-y_step/2)
                for zi in range(-600, 600):
                    z = zi * z_step
                    z_bins.append(z-z_step/2)
                    field = self.structure.get_field_value(0, y, z, t)
                    y_list.append(y)
                    z_list.append(z)
                    bx_field_list.append(field[0])
                    ey_field_list.append(field[4])
                    ez_field_list.append(field[5])
            y_bins.append(y_bins[-1]+y_step)
            z_bins.append(z_bins[-1]+z_step)
            h2d_bx = bx_axes.hist2d(z_list, y_list, weights = bx_field_list, bins = [z_bins, y_bins])
            bx_axes.set_title(f"bx t={round(t)}")
            h2d_ey = ey_axes.hist2d(z_list, y_list, weights = ey_field_list, bins = [z_bins, y_bins])
            ey_axes.set_title(f"ey t={round(t)}")
            h2d_ez = ez_axes.hist2d(z_list, y_list, weights = ez_field_list, bins = [z_bins, y_bins])
            ez_axes.set_title(f"ez t={round(t)}")
            figure.colorbar(h2d_bx[3], ax=bx_axes)
            figure.colorbar(h2d_ey[3], ax=ey_axes)
            figure.colorbar(h2d_ez[3], ax=ez_axes)
            figure.savefig(f"field_{ti}.png")

    def get_dfdx(self, point, xi):
        point = list(point)
        point[xi] += self.delta[xi]
        field1 = self.structure.get_field_value(*point)
        point[xi] -= 2*self.delta[xi]
        field0 = self.structure.get_field_value(*point)
        point[xi] += self.delta[xi]
        deriv = [(field1[i] - field0[i])/2/self.delta[xi] for i in range(6)]
        return deriv

    def get_div_e(self, point):
        deriv = [self.get_dfdx(point, xi)[xi+3] for xi in range(3)]
        return sum(deriv)

    def get_curl_e(self, point):
        deriv = [self.get_dfdx(point, xi) for xi in range(3)] # d(b,e)/dxi
        curl = [
            deriv[1][5] - deriv[2][4], # dy Ez - dz Ey
            deriv[2][3] - deriv[0][5], # dz Ex - dx Ez
            deriv[0][4] - deriv[1][3], # dx Ey - dy Ex
        ]
        return curl


    def test_div_curl_e(self):
        self.structure.set_v0(0.5)
        figure = matplotlib.pyplot.figure()
        axes = figure.add_subplot(1, 1, 1)
        ez_list, ey_list, curl_list, div_list, max_deriv_list = [], [], [], [], []
        y0 = 1.0
        z0 = -355.0
        print(f" i f         Bx         Ey         Ez        Div(E)     Curl(E)")
        for max_deriv in range(10):
            self.structure.set_max_derivative(max_deriv)
            max_deriv_list.append(max_deriv)
            point = (0.0, y0, z0, 0.0)
            field = self.structure.get_field_value(*point)
            div_list.append(self.get_div_e(point))
            deriv = [self.get_dfdx(point, xi) for xi in range(3)]
            curl = self.get_curl_e(point)
            curl_list.append((curl[0]**2+curl[1]**2+curl[2]**2)**0.5)
            fn = self.structure.get_on_axis_field(z0-self.structure.get_z_offset(), max_deriv)
            dBdt = self.get_dfdx(point, 3)
            print(f"{max_deriv: 2d} {fn:+.06f} {field[0]:+.06f} {field[4]:+.06f} {field[5]:+.06f} Div: {div_list[-1]:+.06f} Curl: x {curl[0]:+.06f} {dBdt[0]:+.06f} {curl[0]+dBdt[0]:+.06f} y {curl[1]:+.06f}  {dBdt[1]:+.06f} z {curl[2]:+.06f} {dBdt[2]:+.06f}")
            #print(f"                                                dez/dy: {deriv[1][5]:+.06f}  dey/dz: {deriv[2][4]:+.06f}")

        print("  order", max_deriv_list)
        print("  div  ", div_list)
        print("  curl ", curl_list)
        axes.plot(max_deriv_list, div_list, label=f"$\\nabla . E$")
        axes.plot(max_deriv_list, curl_list, label=f"$\\nabla \\times E$")
        axes.legend()
        figure.savefig("div_e.png")
        print("Fringe field function - derivatives")
        for deriv in range(11):
            print(deriv, self.structure.get_on_axis_field(z0-self.structure.get_z_offset(), deriv))





if __name__ == "__main__":
    unittest.main()
