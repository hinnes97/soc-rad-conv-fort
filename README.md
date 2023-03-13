# soc-rad-conv-fort
Fortran radiative-convective equilibrium code using met office SOCRATES radiation code

## Getting started
The main runscript for the code is `run.py`. The code should run on calling `python run.py`. Under the hood, this runscript calls `./build`, executing the compilation of the fortran code. The makefile is created using the GFDL tool `mkmf`, which should automatically handle the fortran dependencies of the main code

SOCRATES source code is proprietry and not under version control here. The `src` directory should be put in the `src/socrates/` directory next to the `interface` directory.

The parameters of the scheme are set in `input.nml`. 

Currently, the code is only supported using the SOCRATES radiation scheme. There are other radiation schemes available (e.g. semi-grey) that can be activated by changing the `-r SOC` option at the beginning of `run.py` to another compiler option. This prevents the SOCRATES source from being compiled.

This code was created with modelling sub-Neptunes with water in mind. The radiation scheme only includes CH4, H2, He and H2O in ratios specified in `flux.f90` (hard-coded). To add more species, this will need to be hacked manually. Development on a generalised tracer treatment is on a long to do list (see issues).

The code is setup to treat water vapour as a special component in the atmopshere. Depending on the value of the namelist parameter `moisture_scheme`, water will be treated as follows:

- 'none' - Specific humidity is left alone and set to $q = q_0$
- 'surface' - There is a surface water ocean, and the atmosphere is assumed to be saturated everywhere (except above the cold trap). Special care is taken when $q_sat(T,p)>1$. At this point, the atmosphere is set to $q=1$ and the temperature is limited to $q(T^*, p) = 1$. Once this point is reached at the bottom of the atmosphere, the model rarely converges and special treatment is required.
- 'deep' - The specific humidity is set to $\text{min}(q_sat, q_0)$. 
