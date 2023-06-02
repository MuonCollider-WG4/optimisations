# solenoid_optics

You will need:
* ROOT
* Python
* numpy/scipy/matplotlib
I am currently (June 2023) using Python 3.10.6 and ROOT 6.28/02

evolve: routines to calculate matched beta function through a lattice. This is the python script I use most often.

coil_fitter: uses ROOT tminuit to find coils that can be used to "mimic" fields from field_models

final_cooling_match: uses ROOT tminuit to do matching for the final cooling lattice

field_models: set of python classes that contain a bunch of field models

movie: small script to drive mencoder (and ffmpeg) to make animations. This has an additional dependency on mencoder and ffmpeg.

gauss_minimiser: I was experimenting with some genetic algorithm stuff. It is not used.


