import copy
import numpy
import shutil
import os
import subprocess
import xboa.hit
import xboa.bunch

class Setup():
    def __init__(self):
        self.type = None

    def setup(self, config):
        for key, value in config.items():
            if key not in self.__dict__:
                raise KeyError("Did not recognise element configuration", key)
            if self.__dict__[key] is not None:
                value = type(self.__dict__[key])(value)
            self.__dict__[key] = value

class Solenoid(Setup):
    def __init__(self):
        super().__init__()
        self.x_position = 0.0
        self.y_rotation = 0.0 # rotation about the y axis i.e. in x z plane
        self.z_position = 0.0
        self.inner_radius = 0.0
        self.outer_radius = 0.0
        self.length = 0.0
        self.current = 0.0
        self.name = "my_coil"
        self.rgb = [0, 1, 0]
        self.set_coil_name()

    @classmethod
    def clear(cls):
        cls.coil_list = []

    def set_coil_name(self):
        self.coil_name = f"coil_{self.inner_radius:.5g}_{self.outer_radius}_{self.length}"

    def build(self):
        my_solenoid = ""
        self.set_coil_name()
        if self.coil_name not in self.coil_list:
            self.coil_list.append(self.coil_name)
            my_solenoid = \
f"coil {self.coil_name} innerRadius={self.inner_radius} outerRadius={self.outer_radius} length={self.length}\n"
        my_solenoid += \
f"solenoid {self.name} coilName={self.coil_name} current=1 kill=1 color={self.rgb[0]},{self.rgb[1]},{self.rgb[2]}\n"
        my_solenoid += \
f"place {self.name} x={self.x_position} z={self.z_position} current={self.current} rotation=Y{self.y_rotation}\n\n"
        return my_solenoid

    coil_list = []


class Cavity(Setup):
    def __init__(self):
        super().__init__()
        self.name = ""
        self.inner_length = 0.0
        self.frequency = 0.0
        self.max_gradient = 0.0
        self.x_position = 0.0
        self.z_position = 0.0
        self.rgb = None
        self.phase = None
        self.time_offset = None

    def build(self):
        my_cavity = """
pillbox {0} innerLength={1} frequency={2} \\
    maxGradient={3} irisRadius=100.0 \\
    win1Thick=0.0 win2Thick=0.0 wallThick=0.0 collarThick=0.0 \\
    kill=1 maxStep=0.1 innerRadius=500.0""".format(self.name, self.inner_length, self.frequency, self.max_gradient, self.time_offset)
        if self.phase is not None:
            my_cavity += f" phaseAcc={self.phase}"
        if self.time_offset is not None:
            my_cavity += f" timeOffset={self.phase}"
        my_cavity += f"\nplace {self.name} z={self.z_position} x={self.x_position}"
        if self.rgb is not None:
            my_cavity += f" color={self.rgb[0]},{self.rgb[1]},{self.rgb[2]}"
        my_cavity += "\n"
        return my_cavity

class Tube(Setup):
    def __init__(self):
        super().__init__()
        self.name = ""
        self.length = 0.0
        self.outer_radius = 0.0
        self.material = ""
        self.x_position = 0.0
        self.z_position = 0.0

    def build(self):
        my_tube = \
            f"tube {self.name} outerRadius={self.outer_radius} length={self.length} material={self.material}\n"
        my_tube += f"\nplace {self.name} x={self.x_position} z={self.z_position} color=1,0,0\n"
        return my_tube


class Reference(Setup):
    def __init__(self):
        self.particle = "mu+"
        self.p_start = 100.0
        self.z_start = 0.0
        self.t_start = 0.0
        self.no_e_field = 0
        self.no_e_loss = 0

    def build(self):
        my_reference = \
            f"reference particle={self.particle} referenceMomentum={self.p_start} "+\
            f" beamZ={self.z_start} beamX=0.0 beamT={self.t_start} "+\
            f"noEfield={self.no_e_field} noEloss={self.no_e_loss}"
        return my_reference

class Beam(Setup):
    def __init__(self):
        self.filename = ""
        self.out_dir = ""
        self.pid = -13
        self.x_position = 0.0
        self.z_position = 0.0
        self.beams = []
        self.particles = []
        self.default_hit = {"pid":self.pid, "mass":xboa.common.pdg_pid_to_mass[abs(self.pid)]}

    def build(self):
        self.build_beam_file()
        my_beam = """
beam ascii particle={0} nEvents={1} filename={2} format=BLTrackFile beamZ={3} beamX={4}
""".format(self.pid, len(self.particles), self.filename, self.z_position, self.x_position)
        return my_beam

    def build_beam_file(self):
        self.particles = []
        for a_beam in self.beams:
            self.build_a_beam(a_beam)
        bunch = xboa.bunch.Bunch.new_from_hits(self.particles)
        file_name = os.path.join(self.out_dir, self.filename)
        bunch.hit_write_builtin("g4beamline_bl_track_file", file_name)
        print("Wrote file to", file_name)

    def build_a_beam(self, a_beam):
        beam_type = a_beam["type"]
        beam_builder = {
            "longitudinal_grid":self.longitudinal_grid,
            "longitudinal_ellipse":self.longitudinal_ellipse,
            "longitudinal_g4bl_trackfile":self.longitudinal_g4bl,
        }[beam_type]
        beam_builder(a_beam)

    def my_linspace(self, start, stop, num):
        if num == 1:
            return [(start+stop)/2]
        else:
            return numpy.linspace(start, stop, num).tolist()

    def longitudinal_grid(self, a_beam):
        t_list = self.my_linspace(a_beam["t_min"], a_beam["t_max"], a_beam["n_t_steps"])
        e_list = self.my_linspace(a_beam["e_min"], a_beam["e_max"], a_beam["n_e_steps"])
        mass = xboa.common.pdg_pid_to_mass[abs(a_beam["pid"])]
        for t in t_list:
            for e in e_list:
                hit_dict = copy.deepcopy(a_beam["default_hit"])
                hit_dict.update({"t":t, "energy":e+mass, "mass":mass, "pid":a_beam["pid"], "event_number":len(self.particles)+1})
                self.particles.append(xboa.hit.Hit.new_from_dict(hit_dict, "pz"))


    def longitudinal_g4bl(self, a_beam):
        bunch = xboa.bunch.Bunch.new_from_read_builtin("g4beamline_bl_track_file", a_beam["file_name"])
        for i, hit in enumerate(bunch):
            hit["x"] = 0.0
            hit["y"] = 0.0
            hit["z"] = 0.0
            hit["px"] = 0.0
            hit["py"] = 0.0
            if i < 5:
                print(i, "pz", hit["pz"])
            hit.mass_shell_condition("pz")
            if i < 5:
                print(i, "pz", hit["pz"])
            self.particles.append(hit)
        if a_beam["n_particles"]:
            ev_numbers = numpy.linspace(0, len(self.particles)-1, a_beam["n_particles"])
            self.particles = [self.particles[int(ev)] for ev in ev_numbers]

    def longitudinal_ellipse(self, a_beam):
        delta_t = a_beam["delta_t"]
        delta_e = a_beam["delta_e"]
        e_centre = a_beam["e_centre"]
        t_centre = a_beam["t_centre"]
        n_per_dimension = a_beam["n_per_dimension"]
        pid = a_beam["pid"]
        ellipse = [[delta_t**2, 0.0], [0.0, delta_e**2]]
        mean = [t_centre, e_centre]
        shell = xboa.common.make_shell(n_per_dimension, ellipse, mean)
        mass = xboa.common.pdg_pid_to_mass[abs(a_beam["pid"])]
        for item in shell:
            hit_dict = {}#copy.deepcopy(a_beam["default_hit"])
            t = item[0, 0]
            e = item[0, 1]
            hit_dict.update({"t":t, "energy":e+mass, "mass":mass, "pid":pid, "event_number":len(self.particles)+1})
            self.particles.append(xboa.hit.Hit.new_from_dict(hit_dict, "pz"))


class G4BLExecution:
    def __init__(self, linac):
        self.g4bl_path = os.path.expandvars("${HOME}/Software/install/bin/g4bl")
        self.lattice_filename = linac.lattice_filename
        self.linac = linac
        self.guess_logfile()
        self.command_line = []

    def guess_logfile(self):
        out_dir = self.linac.out_dir()
        self.log_filename = os.path.join(out_dir, "log")

    def execute(self):
        command = [self.g4bl_path, os.path.split(self.lattice_filename)[1]]
        cwd = os.getcwd()
        os.chdir(self.linac.out_dir())
        print("Running", command, "in", os.getcwd())
        with open("log", "w") as logfile:
            proc = subprocess.run(
                    command,
                    stdout=logfile, stderr=subprocess.STDOUT)
        print("   ... completed with return code", proc.returncode)
        os.chdir(cwd)
        if proc.returncode:
            raise RuntimeError("G4BL did not execute successfully")

class G4BLLinac:
    def __init__(self, lattice_filename):
        self.lattice_file = None
        self.lattice_filename = lattice_filename
        self.elements = []
        self.reference = {}
        self.beam = {"filename":os.path.join(self.out_dir(), "beam.txt")}
        self.do_stochastics=1
        self.z_spacing = 100.0 # mm
        self.min_z = 0.0 # mm
        self.max_z = 10000.0 # mm
        self.max_step = 100.0 # mm
        self.eps_max = 0.01
        self.output_file = "output_data" # g4bl puts this in the run directory and adds ".txt" as suffix
        self.cleanup_dir = True
        self.linac_name = ""

    def out_dir(self):
        return os.path.split(self.lattice_filename)[0]

    def build_topmatter(self):
        topmatter = ""
        topmatter += f"physics default doStochastics={self.do_stochastics}\n"
        topmatter += f"zntuple cooling_monitor zloop={self.min_z}:{self.max_z}:{self.z_spacing} format=for009 file=output_data coordinates=c\n"
        topmatter += f"param epsMax={self.eps_max}\n" # g4bl bug
        topmatter += f"param maxStep={self.max_step}\n" # g4bl bug
        topmatter += f'#g4ui "/run/beamOn 100" when=4 # visualisation tracking\n'
        self.lattice_file.write(topmatter)

    def build_reference(self):
        my_reference = Reference()
        my_reference.setup(self.reference)
        ref_string = my_reference.build()
        self.lattice_file.write(ref_string)

    def build_beam(self):
        my_beam = Beam()
        my_beam.setup(self.beam)
        beam_string = my_beam.build()
        self.lattice_file.write(beam_string)

    def build_elements(self):
        for element_json in self.elements:
            Element = self.element_dict[element_json["type"]]
            my_element = Element()
            my_element.setup(element_json)
            element_string = my_element.build()
            self.lattice_file.write(element_string)

    def build_solenoids(self):
        for solenoid in self.solenoids:
            my_solenoid = Solenoid()
            my_solenoid.setup(solenoid)
            solenoid_string = my_solenoid.build()
            self.lattice_file.write(solenoid_string)

    def build_linac(self):
        clean_dir(self.out_dir(), self.cleanup_dir)
        with open(self.lattice_filename, "w") as self.lattice_file:
            self.build_topmatter()
            self.build_reference()
            self.build_beam()
            self.build_elements()

    def offset_linac(self, x_offset, z_offset):
        for item in self.elements:
            if "x_position" in item:
                item["x_position"] += x_offset
            if "z_position" in item:
                item["z_position"] += z_offset

    element_dict = {
        "solenoid":Solenoid,
        "cavity":Cavity,
        "tube":Tube,
    }


def clean_dir(my_dir, cleanup):
    if os.path.exists(my_dir):
        if cleanup:
            shutil.rmtree(my_dir)
        else:
            return
    os.makedirs(my_dir)
