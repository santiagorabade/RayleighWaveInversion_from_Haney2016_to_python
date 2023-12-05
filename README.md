# Haney_2016_RayleighInversion
Haney 2016 Rayleigh wave Inversion in Python

This is the Python version for the code of the paper 
Haney, M. M., & Tsai, V. C. (2017). Perturbational 
and nonperturbational inversion of Rayleigh-wave velocities. 
In GEOPHYSICS (Vol. 82, Issue 3, pp. F15–F28). 
Society of Exploration Geophysicists. 
https://doi.org/10.1190/geo2016-0397.1 

You can find the original version here: 
https://github.com/matt-haney/raylee_codes

Code requirements:
Python 3.8
-Numpy
-Scipy

##Instructions

All the inputs are generated using the code get_initial_model.py  
from the dsp_file_1.dat file
The dsp_file.dat contains 5 columns. Period (s), Velocity (m/s), Uncertainty (%),
Mode number (Fundamental=1, First overtone=2, etc), Vector of velocity type, either phase (0) or group (1)

Inputs:
vp_init.txt               Initial P-wave velocity model in solid
vs_init.txt               Initial S-wave velocity model in solid
rho_init.txt              Initial density model in solid
vpf.txt                   P-wave velocity in fluid layer (if no water layer, this is not used)
rhof.txt                  Density in fluid layer (if no water layer, this is not used)
grid_values_solid.txt     Finite element grid in solid layer
grid_values_fluid.txt     Finite element grid in fluid layer (if no water layer, this is not used)
input_params.txt          See description below

The file "input_params.txt" contains

% flag for fixed poisson's ratio (0=no,1=yes)
% smoothness scale (m)
% a priori model standard deviation factor
% maximum number of updates (iterations)
% number of measurements
% number of elements in solid part of model
% number of elements in fluid part of model
% lower chi squared window
% higher chi squared window

The outputs are
final_vs.txt              Depth in (km) and final Vs (m/s)
fw_group.txt              Period in (s) surface wave velocity (m/s)

To run:
python3 driver.py 1 &
