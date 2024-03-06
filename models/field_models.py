import json
import sys
import math
import numpy
import scipy.integrate
import scipy.interpolate
import scipy.fft

mu0=1.25663706212e-6 # [kg m s^-2 A^-2]; for reference field is in [T] = [kg s^-2 A^-1]

class Field(object):
    def __init__(self):
        pass

    def reset(self):
        pass # if needed, reset the fields based on member parameters

    def get_field(self, z):
        raise NotImplementedError("Should overload this function")

    def get_period(self):
        raise NotImplementedError("Should overload this function")

    def bz2_int(self, y, z):
        dbzsquareddz = self.get_field(z)**2
        return dbzsquareddz

    def get_bz2_int(self):
        intbz2 = scipy.integrate.odeint(self.bz2_int, 0., [0., self.period])[1][0]
        return intbz2

    def get_name(self):
        return "field"

    def human_readable(self):
        return ""

class FieldSum(Field):
    def __init__(self, field_list):
        self.field_list = field_list
        self.period = field_list[0].get_period()
        for field in self.field_list:
            if abs(field.period - self.period) > self.period/1e-9:
                raise RuntimeError("Field sum but periods dont match")

    def get_field(self, z):
        bz = 0.0
        for field in self.field_list:
            bz += field.get_field(z)
        return bz

    def get_period(self):
        return self.period


class FlatGaussField(Field):
    def __init__(self, period, width, peak, tail):
        self.period = period
        self.width = width
        self.peak = peak
        self.tail = tail

    def get_field(self, z0):
        z1 = z0-self.period
        b0 = (self.peak-self.tail)
        bz = self.tail + \
             b0*math.exp(-z0**2/2./self.width**2) + \
             b0*math.exp(-z1**2/2./self.width**2)
        return bz

    def get_period(self):
        return self.period

    def normalise_bz_squared(self, intbz2dz):
        intbz2 = self.get_bz2_int()
        print("Integral", intbz2)
        self.peak *= (intbz2dz/intbz2)**0.5
        self.tail *= (intbz2dz/intbz2)**0.5
        intbz2_new = scipy.integrate.odeint(self.bz2_int, [0.], [self.period])[0]
        print(intbz2, intbz2_new, intbz2dz)

    def get_name(self):
        return "gauss"

class LinearInterpolation(Field):
    def __init__(self, points, period):
        self.points = points
        self.period = period
        self.name = None
        x = numpy.linspace(0, period, len(points))
        y = numpy.array(self.points)
        print("Linear interpolation in x", x)
        self.interpolation = scipy.interpolate.interp1d(x, y, kind="cubic")

    def get_field(self, z0):
        while z0 > self.period:
            z0 -= self.period 
        while z0 < 0.0:
            z0 += self.period 
        return numpy.real(self.interpolation(z0))

    def get_period(self):
        return self.period

    def get_name(self):
        if self.name != None:
            return self.name
        return "interpolation"

class SineField(Field):
    def __init__(self, bz0, bz1, bz2, bz3, bz4, bz5, period):
        self.bz0 = bz0
        self.bz_list = [bz1, bz2, bz3, bz4, bz5]
        self.period = period

    def get_field(self, z):
        bz = self.bz0
        for i, b0 in enumerate(self.bz_list):
            bz += b0*math.sin(2*(i+1)*math.pi*z/self.period)
        return bz

    def get_period(self):   
        return self.period

    def normalise_bz_squared(self, intbz2dz):
        intbz2 = self.get_bz2_int()
        print("Integral", intbz2)
        self.bz0 *= (intbz2dz/intbz2)**0.5
        for i, bz in enumerate(self.bz_list):
            self.bz_list[i] *= (intbz2dz/intbz2)**0.5
        intbz2_new = scipy.integrate.odeint(self.bz2_int, [0.], [self.period])[0]
        print(intbz2, intbz2_new, intbz2dz)

    def get_name(self):
        name = "sine_"+str(self.bz0)
        for b0 in self.bz_list:
            name += "_"+"{:4.4f}".format(b0)
        return name

    def human_readable(self):
        name = "L={1:.4g}; b_0={0:.4g}".format(self.bz0, self.period)
        for i, b0 in enumerate(self.bz_list):
            name += "; b_{0}={1:.4g}".format(i+1, b0)
        return name


class UniformField(Field):
    def __init__(self, bz0):
        self.bz0 = bz0

    def get_field(self, z):
        return self.bz0

    def get_period(self):
        return 1.0
  
    def get_name(self):
        return "uniform_field"
 
    def normalise_bz_squared(self, intbz2dz):
        pass

class CurrentSheet(Field):
    def __init__(self, b0, zcentre, length, radius, period, nrepeats):
        self.zcentre = zcentre
        self.length = length
        self.radius = radius
        self.period = period
        self.nrepeats = nrepeats
        self.b0 = b0

    def get_field(self, z):
        bz = self.get_one_field(z, self.zcentre)
        for i in range(1, self.nrepeats):
            zcentre = self.zcentre+self.period*i
            bz += self.get_one_field(z, zcentre)
            zcentre = self.zcentre-self.period*i
            bz += self.get_one_field(z, zcentre)
        return bz

    def get_period(self):
        return self.period

    def get_one_field(self, z, zcentre):
        delta_z = z-zcentre
        z1 = delta_z-self.length/2.0
        z2 = delta_z+self.length/2.0
        b1 = z1/(z1**2+self.radius**2)**0.5
        b2 = z2/(z2**2+self.radius**2)**0.5
        field = self.b0*(b2-b1)/2.0
        return field

    def get_current_per_length(self):
        global mu0
        # b0 = mu0 I / 2 -> I/l = 2 b0/mu0/l
        return self.b0/mu0 #[A/m]

    def get_off_axis_field(self, z, r):
        """
        Return the field off-axis
        - z position along the current sheet axis
        - r position radially (perpendicular to the current sheet axis)
        Returns a two vector bz, br
        """
        return sum([coil.get_field(z) for coil in self.coil_list])

    def get_name(self):
        return "sheet"

class CurrentBlock(Field):
    def __init__(self, b0, zcentre, length, rmin, rmax, period, nrepeats, nsheets):
        self.zcentre = zcentre
        self.length = length
        self.rmin = rmin
        self.rmax = rmax
        self.period = period
        self.nrepeats = nrepeats
        self.b0 = b0
        self.nsheets = nsheets
        self.coil_list = []
        self.reset()

    def reset(self):
        rstep = (self.rmax-self.rmin)/self.nsheets
        self.coil_list = [
            CurrentSheet(self.b0/self.nsheets, self.zcentre, self.length, self.rmin+rstep*(i+0.5), self.period, self.nrepeats)
            for i in range(self.nsheets)
        ]

    def get_field(self, z): # not yet
        return sum([coil.get_field(z) for coil in self.coil_list])

    def get_off_axis_field(self, z, r): # not yet
        raise NotImplementedError()
        return sum([coil.get_field(z) for coil in self.coil_list])


    def get_current_density(self):
        # assumes nothing weird like some current sheets longer than others
        current_per_length = 0.0
        for coil in self.coil_list:
            current_per_length += coil.get_current_per_length()
        current_density = current_per_length / (self.rmax-self.rmin)
        return current_density # A m^-2

    def get_period(self):
        return self.coil_list[0].get_period()
