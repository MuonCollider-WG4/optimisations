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
                raise KeyError(f"Did not recognise element configuration {key} for {config['type']} element {config['name']}")
            if self.__dict__[key] is not None and config[key] is not None:
                value = type(self.__dict__[key])(value)
            self.__dict__[key] = value

class GenericElement(Setup):
    def __init__(self):
        super().__init__()
        self.element_type = ""
        self.name = ""
        self.x_position = 0.0
        self.y_position = 0.0
        self.z_position = 0.0
        self.rgb = [1, 0, 0]
        self.element_parameters = {}
        self.place_parameters = {}

    def build(self):
        raise NotImplementedError("Not sure this is the right way")
        my_element = f"{self.element_type} {self.name}"
        for key, value in self.element_parameters:
            my_element += f" {key}={value}"
        my_element += "\n"
        if self.place_parameters != None:
            for key, value in self.place_parameters:
                my_element += f" {key}={value}"

    def setup(self, config):
        """Setup overloading base class; this type does not throw an exception"""
        for key, value in config.items():
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
        self.nsheets = 10
        self.name = "my_coil"
        self.rgb = [0, 1, 0]
        self.set_coil_name()

    @classmethod
    def clear(cls):
        cls.coil_list = []

    def set_coil_name(self):
        self.coil_name = f"coil_{self.inner_radius:.5g}_{self.outer_radius}_{self.length}_{self.nsheets}"

    def build(self):
        my_solenoid = ""
        self.set_coil_name()
        if self.coil_name not in self.coil_list:
            self.coil_list.append(self.coil_name)
            my_solenoid = f"coil {self.coil_name} innerRadius={self.inner_radius} outerRadius={self.outer_radius} length={self.length} nSheets={self.nsheets}\n"
        colour = " ".join(self.rgb)
        my_solenoid += f"solenoid {self.name} coilName={self.coil_name} current=1 kill=1 color={colour}\n"
        my_solenoid +=  f"place {self.name} x={self.x_position} z={self.z_position} current={self.current} rotation=Y{self.y_rotation}\n\n"
        return my_solenoid

    coil_list = []

class DerivativesSolenoid(Setup):
    def __init__(self):
        super().__init__()
        self.x_position = 0.0
        self.y_rotation = 0.0 # rotation about the y axis i.e. in x z plane
        self.z_position = 0.0
        self.length = None
        self.max_derivative = 9
        self.max_r = 500
        self.nominal_field = 0.0
        self.end_length = 100.0
        self.centre_length = 1000.0
        self.name = ""
        self.field_model = "tanh"
        self.harmonics = [] # field model = "fourier"

    def build(self):
        my_ds = f"""
derivativessolenoid {self.name} fieldModel={self.field_model} length={self.length} maxDerivative={self.max_derivative} maxR={self.max_r} \
                    """
        if self.field_model == "tanh":
            my_ds += f"nominalField={self.nominal_field} endLength={self.end_length} centreLength={self.centre_length}"
        elif self.field_model == "fourier":
            for i, value in enumerate(self.harmonics):
                my_ds += f" harmonic{i}={value}"
        my_ds += "\nplace {self.name} x={self.x_position} z={self.z_position} rotation=Y{self.y_rotation}\n\n"
        return my_ds

class Cavity(Setup):
    def __init__(self):
        super().__init__()
        self.name = ""
        self.inner_length = 0.0
        self.frequency = 0.0
        self.max_gradient = 0.0
        self.x_position = 0.0
        self.z_position = 0.0
        self.collarRadialThick = 0.0001
        self.iris_radius = 100.0
        self.rgb = None
        self.phase = None
        self.time_offset = None

    def build(self):
        my_cavity = """
pillbox {0} innerLength={1} frequency={2} \\
    maxGradient={3} irisRadius={5} \\
    win1Thick=0.0 win2Thick=0.0 wallThick=1.0 collarThick=1.0 collarRadialThick={6} \\
    kill=1 maxStep=0.1 innerRadius=500.0""".format(self.name, self.inner_length, self.frequency, self.max_gradient, self.time_offset, self.iris_radius, self.collarRadialThick)
        if self.phase is not None:
            my_cavity += f" phaseAcc={self.phase}"
        if self.time_offset is not None:
            my_cavity += f" timeOffset={self.phase}"
        my_cavity += f"\nplace {self.name} z={self.z_position} x={self.x_position}"
        if self.rgb is not None:
            my_cavity += f" color={self.rgb[0]},{self.rgb[1]},{self.rgb[2]}"
        my_cavity += "\n"
        return my_cavity

class UniformField(Setup):
    def __init__(self):
        super().__init__()
        self.name = ""
        self.radius = 0.0
        self.bz = 0.0
        self.z_start = 0.0
        self.x_position = 0.0
        self.z_position = 0.0
        self.length = 0.0

    def build(self):
        if self.z_start != 0.0:
            if self.z_centre != 0.0:
                raise ValueError(f"Can only define either z_centre or z_start, not both in UniformField {self.name}")
            self.z_position = self.z_start + self.length/2.0

        my_field = \
            f"fieldexpr {self.name} Bz={self.bz} length={self.length} radius={self.radius}\n"
        my_field += f"place {self.name} z={self.z_position}\n"
        return my_field

class Tube(Setup):
    def __init__(self):
        super().__init__()
        self.name = ""
        self.length = 0.0
        self.inner_radius = 0.0
        self.outer_radius = 0.0
        self.material = ""
        self.x_position = 0.0
        self.z_position = 0.0 # if None will not place

    def build(self):
        my_tube = \
            f"tube {self.name} innerRadius={self.inner_radius} outerRadius={self.outer_radius} length={self.length} material={self.material}\n"
        if self.z_position != None:
            my_tube += f"\nplace {self.name} x={self.x_position} z={self.z_position} color=1,0,0\n"
        else:
            my_tube += "\n"
        return my_tube




class Detector(Setup):
    def __init__(self):
        super().__init__()
        self.name = ""
        self.solid = ""
        self.format = ""
        self.filename = ""
        self.coordinates = "local"
        self.x_position = 0.0
        self.z_position = 0.0

    def build(self):
        my_detector = \
            f"detector {self.name} solid={self.solid} format={self.format} filename={self.filename} coordinates={self.coordinates}\n"
        my_detector += f"place {self.name} x={self.x_position} z={self.z_position} color=1,0,0\n\n"
        return my_detector


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
        self.start_event_number = 0
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
        print("Wrote beam file to", file_name)

    def build_a_beam(self, a_beam):
        beam_type = a_beam["type"]
        beam_builder = {
            "longitudinal_grid":self.longitudinal_grid,
            "longitudinal_ellipse":self.longitudinal_ellipse,
            "longitudinal_g4bl_trackfile":self.longitudinal_g4bl,
            "full_g4bl_trackfile":self.full_g4bl,
            "beam_ellipse":self.beam_ellipse,
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
                hit_dict.update({"t":t, "energy":e+mass, "mass":mass, "pid":a_beam["pid"], "event_number":len(self.particles)+self.start_event_number+1})
                self.particles.append(xboa.hit.Hit.new_from_dict(hit_dict, "pz"))

    def beam_ellipse(self, a_beam):
        var_list = a_beam["variables"]
        ellipse = a_beam["ellipse"]
        mean = a_beam["mean"]
        n_particles = a_beam["n_particles"]
        mass_shell = a_beam["mass_shell_condition"]
        psv_list = numpy.random.multivariate_normal(mean, ellipse, n_particles)
        for psv in psv_list:
            hit_dict = copy.deepcopy(a_beam["default_hit"])
            for i, var in enumerate(var_list):
                hit_dict[var] = psv[i]
            hit_dict["event_number"] = len(self.particles)+self.start_event_number+1
            self.particles.append(xboa.hit.Hit.new_from_dict(hit_dict, mass_shell))

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
            ev_numbers = numpy.linspace(self.start_event_number, self.start_event_number+len(self.particles)-1, a_beam["n_particles"])
            self.particles = [self.particles[int(ev)] for ev in ev_numbers]

    def full_g4bl(self, a_beam):
        bunch = xboa.bunch.Bunch.new_from_read_builtin("g4beamline_bl_track_file", a_beam["file_in"])
        if "time_distribution" in a_beam:
            if a_beam["time_distribution"] == "uniform":
                t_list = numpy.random.uniform(a_beam["t_min"], a_beam["t_max"], len(bunch))
            elif a_beam["time_distribution"] == "gaussian":
                t_list = numpy.random.norm(a_beam["mean"], a_beam["scale"], len(bunch))
            if a_beam["mode"] == "add":
                for i, hit in enumerate(bunch):
                    hit["t"] += t_list[i]
            elif a_beam["mode"] == "replace":
                for i, hit in enumerate(bunch):
                    hit["t"] = t_list[i]
            else:
                raise KeyError("Did not recognise bunch time distribuction mode")
        for i, hit in enumerate(bunch):
            self.particles.append(hit)
        if a_beam["n_particles"]:
            ev_numbers = numpy.linspace(self.start_event_number, self.start_event_number+len(self.particles)-1, a_beam["n_particles"])
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
            hit_dict.update({"t":t, "energy":e+mass, "mass":mass, "pid":pid, "event_number":len(self.particles)+1+self.start_event_number})
            self.particles.append(xboa.hit.Hit.new_from_dict(hit_dict, "pz"))


class G4BLExecution:
    def __init__(self, linac):
        self.g4bl_path = os.path.expandvars("${HOME}/Software/install/bin/g4bl")
        self.lattice_filename = linac.lattice_filename
        self.linac = linac
        self.np = 1
        self.guess_logfile()
        self.command_line = []

    def guess_logfile(self):
        out_dir = self.linac.out_dir()
        self.log_filename = os.path.join(out_dir, "log")

    def execute(self):
        command = [self.g4bl_path, os.path.split(self.lattice_filename)[1]]
        if self.np != 1:
            command += ["-np", "4"]
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

class TrackCuts(Setup):
    def __init__(self):
        self.keep = []
        self.kill = []

    def build(self):
        track_cuts = "trackcuts my_track_cuts "
        if len(self.keep):
            track_cuts += "keep="
            for item in self.keep:
                name = self.string_lookup(item)
                track_cuts += f"{name},"
            track_cuts = track_cuts[:-1]
        if len(self.kill):
            track_cuts += "kill="
            for item in self.kill:
                track_cuts += f"{item},"
            track_cuts = track_cuts[:-1]
        track_cuts += "\n"
        return track_cuts

    @classmethod
    def string_lookup(cls, pid):
        return {
            2212:"proton",
            211:"pi+",
            -211:"pi-",
            13:"mu-",
            -13:"mu+",
            11:"e-",
            -11:"e+",
        }[pid]

class Physics(Setup):
    def __init__(self):
        self.physics_list = "default"
        self.do_stochastics = None
        self.disable = None

    def build(self):
        physics = "physics"
        if not self.physics_list is None:
            physics += f" {self.physics_list}"
        if not self.do_stochastics is None:
            physics += f" doStochastics={self.do_stochastics}"
        if not self.disable is None:
            physics += f" disable={self.disable}"
        physics += "\n"
        return physics


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
        self.track_cuts = TrackCuts()
        self.fieldntuple_step = 10 # 0/negative to disable field nturple
        self.physics = Physics()

    def out_dir(self):
        return os.path.split(self.lattice_filename)[0]

    def build_topmatter(self):
        topmatter = ""
        topmatter += self.physics.build()
        topmatter += f"zntuple cooling_monitor zloop={self.min_z}:{self.max_z}:{self.z_spacing} format=for009 file=output_data coordinates=c\n"
        topmatter += f"param epsMax={self.eps_max}\n" # g4bl bug
        topmatter += f"param maxStep={self.max_step}\n" # g4bl bug
        if self.fieldntuple_step > 0:
            topmatter += f"fieldntuple field_out format=ascii filename=fieldmap.dat x=0,1,100 z=0,{self.max_z},{self.fieldntuple_step}\n"
        topmatter += f'#g4ui "/run/beamOn 100" when=4 # visualisation tracking\n'
        if self.track_cuts:
            track_cuts = self.track_cuts.build()
            topmatter += track_cuts
        self.lattice_file.write(topmatter)

    def build_reference(self):
        # note I had to hack G4BL for009 output to cope with multi ref particles
        if type(self.reference) != type([]):
            self.reference = [self.reference]
        for a_ref_dict in self.reference:
            my_reference = Reference()
            my_reference.setup(a_ref_dict)
            ref_string = my_reference.build()
            self.lattice_file.write(ref_string+"\n")

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
        "derivatives_solenoid":DerivativesSolenoid,
        "solenoid":Solenoid,
        "uniform_field":UniformField,
        "cavity":Cavity,
        "tube":Tube,
        "detector":Detector,
    }


def clean_dir(my_dir, cleanup):
    if os.path.exists(my_dir):
        if cleanup:
            shutil.rmtree(my_dir)
        else:
            return
    os.makedirs(my_dir)
