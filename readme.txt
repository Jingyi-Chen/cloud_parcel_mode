maincpm.f90: 
- temporary wrapper
- call chem to set the aerosol chemistry properties
- call aerospec to set aerosol size distribution
- call vertvelo to set the vertical velocity (default 0)
- call entrain to set the entrianment mixing properties (defaul entrainment rate is 0)
- call dropspec to set the droplet size distribution
- call cpm to calculate the deliquesce/efflorescence and condensation/evaporation


cpm.f90: 
- the most primary code to calculate the deliquesce/efflorescence and condensation/evaporation
- call ODE solver to solve the equations in cond.f90
- units are g/cm/s
- aerostate:
     0: dry aerosol
     1: intistitial aerosol
     2: larger than r_c but smaller than r_ki
     4: larger than both r_c and r_ki
     r_c:  critical radius defined in Reutter et al., 2009
     r_ki: Kohler critical radius for each bin

cond.f90:
- time tendency equations

DVODE_F90_M.f90
- ODE solver
