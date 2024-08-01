import json
import math

import xboa.hit

import Configuration
import maus_cpp.globals
import maus_cpp.global_error_tracking
import maus_cpp.field

def get_phi_start_end(tp_0, tp_1, bz):
    """
    Get azimuthal angles
    - tp_0: track point with lower z
    - tp_1: track point with higher z
    - bz: magnetic field
    Returns azimuthal angles, sorted for direction of rotation. If more than one
    rotation, add the delta phi on (so that we get extra turns).
    """
    c_light = xboa.common.constants['c_light']
    q = tp_0["hit"]["charge"]
    if q*bz > 0:
        phi_start = math.atan2(tp_0["hit"]["px"], tp_0["hit"]["py"])
        phi_end = math.atan2(tp_1["hit"]["px"], tp_1["hit"]["py"])
    elif q*bz < 0:
        phi_start = math.atan2(tp_1["hit"]["px"], tp_1["hit"]["py"])
        phi_end = math.atan2(tp_0["hit"]["px"], tp_0["hit"]["py"])
    k0 = q*bz*c_light/tp_1["hit"]["pz"]
    while phi_start < 0.: # phi_start is always in domain 0 < 2 pi
        phi_start += 2*math.pi
    while phi_end < phi_start: # phi_end is always > phi_start
        phi_end += 2*math.pi
    # for larger steps, we can make more than one rotation
    delta_phi = k0*(tp_0["hit"]["z"] - tp_1["hit"]["z"])
    n_rot = abs(delta_phi)/math.pi/2.
    #print "    delta_phi:", delta_phi, "n_rot:", n_rot, "z:", tp_0["hit"]["z"], tp_1["hit"]["z"]
    phi_end += math.floor(n_rot)*2.*math.pi
    return phi_start, phi_end

def get_max_phi_r2(hit, bz, phi_start):
    """
    Get the azimuthal angle at which the trajectory has maximum excursion from
    the axis
    - tp_0: track point with lower z
    - tp_1: track point with higher z
    - bz: magnetic field
    - phi_start: phi at which we start
    Returns azimuthal angles, sorted for direction of rotation. We add on 2pi so
    that it is always greater than phi_start.
    """
    c_light = xboa.common.constants['c_light']
    q = tp_0["hit"]["charge"]
    if q*bz > 0:
        tp_start = tp_0
    elif q*bz < 0:
        tp_start = tp_1
    # coordinates of helix centre
    pt = tp_start["hit"]["pt"]
    x0 = tp_start["hit"]["x"] + pt*math.cos(phi_start)/q/bz/c_light
    y0 = tp_start["hit"]["y"] - pt*math.sin(phi_start)/q/bz/c_light
    r0 = abs(pt/q/bz/c_light) # wavenumber is qBm/pz dz
    #print 'x0:', x0, 'y0:', y0, 'r0:', r0, 'k0:', k0
    # azimuthal angle of maximum extent (when r is in line with x0, y0)
    direction = abs(q*bz)/q/bz
    phi_max = math.atan2(direction*y0, -direction*x0)
    while phi_max < phi_start: # phi_max is always > phi_start
        phi_max += 2*math.pi
    max_r2 = ((x0**2 + y0**2)**0.5 + r0)**2 # vector to centre + radius of track
    return phi_max, max_r2


def get_r_max(track_point_list, bz):
    """
    Set the maximum r2 on each track point between tp and the next tp that
    the track excurts from the axis. Estimate the region between adjacent track
    points by extrapolating assuming uniform field and calculate the excursion
    analytically. Sets the max_r2 parameter on the track.
    """
    try:
        q = track_point_list[0]["hit"]["charge"]
    except KeyError:
        print track_point_list[0]["hit"], track_point_list[0]["detector"]
        raise
    if abs(q) < 1e-9:
        raise ValueError("Attempt to find r max on track with no charge")
    track_point_list = sorted(track_point_list, key = lambda tp: tp["hit"]["z"])
    c_light = xboa.common.constants['c_light']
    for i, tp_0 in enumerate(track_point_list[:-1]):
        tp_1 = track_point_list[i+1]
        # azimuthal angle at start and end
        phi_start, phi_end = get_phi_start_end(tp_0, tp_1, bz)
        # azimuthal angle and radius at maximum excursion from z-axis
        phi_max, max_r2 = get_max_phi_r2(tp_0, tp_1, bz, phi_start)
        r2_start = tp_0["hit"]["x"]**2+tp_0["hit"]["y"]**2
        r2_end = tp_1["hit"]["x"]**2+tp_1["hit"]["y"]**2
        if phi_max > phi_start and phi_end > phi_max: # phi_max is between phi_end and phi_start
            tp_0["max_r2"] = max_r2
        else:
            tp_0["max_r2"] = max(r2_start, r2_end)
        #print "   tp0 p:", tp_0["hit"]["px"], tp_0["hit"]["py"], "tp1 p:", tp_1["hit"]["px"], tp_1["hit"]["py"]
        #print "   ", str(round(phi_start, 1)).rjust(8), str(round(phi_end, 1)).rjust(8), str(round(phi_max, 1)).rjust(8), str(round(r2_start, 1)).rjust(8), str(round(max_r2, 1)).rjust(8), str(round(r2_end, 1)).rjust(8)

def hit_from_psv(psv):
    hit_dict = {
        "t":psv[0],
        "x":psv[1],
        "y":psv[2],
        "z":psv[3],
        "energy":psv[4],
        "px":psv[5],
        "py":psv[6],
        "pz":psv[7],
        "charge":1.,
        "mass":xboa.common.pdg_pid_to_mass[13],
    }
    hit = xboa.hit.Hit.new_from_dict(hit_dict, "")
    detector_hit = {
        "hit":hit,
        "max_r2":0.,
    }
    return detector_hit

def test_get_r_max(bz, x0, y0, px0, py0, pz0, z_list):
    setup(bz)
    energy = (px0**2+py0**2+pz0**2+105.658**2)**0.5
    var = [0., x0, y0, float(z_list[0]), energy, px0, py0, pz0]
    ellipse = [[0. for i in range(6)] for j in range(6)]
    tracking = maus_cpp.global_error_tracking.GlobalErrorTracking()
    track_list = []
    z_0 = z_list[0]
    for z_1 in z_list:
        max_r2 = var[1]**2+var[2]**2
        for i in range(z_0, z_1):
            var, ellipse = tracking.propagate_errors(var, ellipse, i)
            r2 = var[1]**2+var[2]**2
            max_r2 = max(r2, max_r2)
        track_list.append((hit_from_psv(var), max_r2))
        z_0 = z_1

    print "Tracking output:"
    for tp, max_r2 in track_list:
        hit = tp["hit"]
        r2 = hit["y"]**2+hit["x"]**2
        print "   ", str(round(hit["x"], 1)).rjust(8), str(round(hit["y"], 1)).rjust(8), str(round(hit["z"], 1)).rjust(8), str(round(r2, 1)).rjust(8), str(round(max_r2, 1)).rjust(8)
    max_r2 = max([r2 for tp, r2 in track_list])
    track_list = [tp for tp, r2 in track_list] 
    print "Helix extrapolation output:"
    get_r_max(track_list, bz)
    extrapolation_max = max([tp["max_r2"] for tp in track_list])
    print "Tracking:            ", max_r2
    print "Helix Extrapolation: ", extrapolation_max
    if abs(extrapolation_max - max_r2) > 1:
        raise RuntimeError("Test Fail!")


def setup(bz):
    if not maus_cpp.globals.has_instance():
        config = Configuration.Configuration().getConfigJSON()
        maus_cpp.globals.birth(config)
    module = maus_cpp.mice_module.MiceModule("scripts/dedx/test_track.dat")
    children = module.get_children()
    for child in children:
        if child.get_name() == "Field":
            child.set_property("ConstantField", "Hep3Vector", {"x":0., "y":0., "z":bz})
            break
    module.set_children(children)
    maus_cpp.globals.set_monte_carlo_mice_modules(module)
    print "Set bz to ", maus_cpp.field.get_field_value(0., 0., 0., 0.)[2]

def main():
    """
    test_get_r_max(3e-3, 50., 100., 20., 30., 140., [100, 200, 300, 400, 500])
    print "\n\n*****************************************************************"
    test_get_r_max(3e-3, 50., 100., 20., 30., 140., [0, 5000])
    print "\n\n*****************************************************************"
    test_get_r_max(3e-3, 50., 100., 20., 30., 140., [0, 50])
    print "\n\n*****************************************************************"
    test_get_r_max(3e-3, 50., 100., 20., 30., 140., [0, -50])
    """
    for px_quad in [1, -1]:
        for py_quad in [1, -1]:
            for bz in [3e-3, -3e-3]:
                print "\n\n*****************************************************************"
                test_get_r_max(bz, 50., 100., px_quad*20., py_quad*30., 140., [100, 300]) # not quite through phi_max (bz)
                print "\n\n*****************************************************************"
                test_get_r_max(bz, 50., 100., px_quad*20., py_quad*30., 140., [100, 400]) # not quite through phi_max (-bz); just pass phi_max (bz)
                print "\n\n*****************************************************************"
                test_get_r_max(bz, 50., 100., px_quad*20., py_quad*30., 140., [100, 500]) # just past phi_max (-bz)
                print "\n\n*****************************************************************"
                test_get_r_max(bz, 50., 100., px_quad*20., py_quad*30., 140., [100, 1000]) # almost 2pi
                print "\n\n*****************************************************************"
                test_get_r_max(bz, 50., 100., px_quad*20., py_quad*30., 140., [100, 1100]) # just past 2pi

if __name__ == "__main__":
    main()
  
