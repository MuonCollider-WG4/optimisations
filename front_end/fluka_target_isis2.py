import multiprocessing
import sys
import os
import shutil
import subprocess
import xboa.common
import time
import datetime

import g4bl_interface.stripper

class RunFluka:
    def __init__(self, config):
        self.config = config

    def setup_run_dir(self, run_index):
        a_config = self.config["config_list"][run_index]
        run_dir = a_config["run_dir"]
        if self.config["clean_dir"]:
            if os.path.exists(run_dir):
                shutil.rmtree(run_dir)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        print("    Setup dir", run_dir)

    def get_target_file_name(self, run_index):
        a_config = self.config["config_list"][run_index]
        lattice_file = os.path.split(self.config["lattice_source"])[1]
        run_dir = a_config["run_dir"]
        target_file_name = os.path.join(run_dir, lattice_file)
        return target_file_name

    def setup_lattice(self, run_index):
        a_config = self.config["config_list"][run_index]
        self.setup_run_dir(run_index)
        target_file_name = self.get_target_file_name(run_index)
        for key, value in a_config["subs"].items():
            try:
                value = float(value)
                a_config["subs"][key] = f"{value:4.2e}"
            except ValueError:
                pass
        xboa.common.substitute(self.config["lattice_source"], target_file_name, a_config["subs"])

    def run_lattice(self, run_index):
        exe = self.config["executable"]
        target_dir, target_file_name = os.path.split(self.get_target_file_name(run_index))
        here = os.getcwd()
        os.chdir(target_dir)
        log = open(self.config["logfile"], "w")
        command = [exe, target_file_name]
        print("    Running ", " ".join(command))
        print()
        try:
            proc = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            if len(self.config["g4bl_strip"]):
                for my_file in self.config["g4bl_strip"]:
                    g4bl_interface.stripper.Stripper.G4BLStripper().strip(my_file)
            print(f"    Completed job in dir '{target_dir}' with return {proc.returncode}")
        except Exception:
            sys.excepthook(*sys.exc_info())
            print(f"    Caught exception")
        os.chdir(here)

    def one_function(self, run_index):
        print("Run config", run_index+1)
        print("    Date time", datetime.datetime.now())
        self.setup_lattice(run_index)
        self.run_lattice(run_index)

    def run_single_threaded(self):
        for run_index, a_config in enumerate(self.config["config_list"]):
            self.one_function(run_index)

    def run_pool(self):
        with multiprocessing.Pool(5) as my_pool:
            result_list = []
            print("Building pool")
            for i in range(len(self.config["config_list"])):
                result_list.append(my_pool.apply_async(self.one_function, (i,)))
                time.sleep(1)
            start_time = time.time()
            running_count, done_count = 1, 0
            while running_count > 0:
                running_count, done_count = 0, 0
                for result in result_list:
                    if result.ready():
                        done_count += 1
                    else:
                        running_count += 1
                now_time = round(time.time() - start_time)
                print(f"\r{now_time:5d}   Waiting for {running_count} processes with {done_count} finished", end="")
                time.sleep(5)
            print("Finished")

def energy_scan_fluka_config():
    config = {
        "executable":"/home/cr67/Software/fluka/fluka4-4.1/bin/rfluka",
        "lattice_source":"lattice/plate_target/plate_target.inp",
        "clean_dir":True,
        "logfile":"fluka.log",
        "g4bl_strip":[],
        "config_list":[{
            "subs":{
                "__energ__":energy, # GeV
                "__n_p__":2e7, # number of primaries
            },
            "run_dir":f"output/fluka_model_v8/energy_{-energy:2.1f}/"
        } for energy in [-0.1*i for i in range(22, 31, 2)]],
    }
    return config

def energy_scan_geant4_config(do_left):
    config_list = []
    for energy in [0.1*i for i in range(4, 31, 2)]:
        for run_id in range(5):
            config_list.append({
                "subs":{
                    "__momentum__":1e3*((energy+0.938272)**2-0.938272**2)**0.5, # GeV
                    "__n_p__":2e8, # number of primaries
                    "__run_id__":run_id,
                    "__do_left__":do_left,
                },
                "run_dir":f"output/g4_plate_model_v3/energy_{energy:2.1f}_run-id_{run_id}/"
            })

    config = {
        "executable":"/home/cr67/Software/G4Beamline/G4beamline-3.08-source/build/g4bl/g4bl",
        "lattice_source":"lattice/plate_target_geant4/target_region.g4bl",
        "clean_dir":True,
        "logfile":"sim.log",
        "config_list":config_list,
        "g4bl_strip":["target_right.txt", "target_left.txt"][do_left:do_left+1],
    }
    return config

def main():
    run = RunFluka(energy_scan_geant4_config(do_left=0))
    try:
        run.run_pool()
    except Exception:
        sys.excepthook(*sys.exc_info())
        print("Terminated abnormally...")

if __name__ == "__main__":
    main()

