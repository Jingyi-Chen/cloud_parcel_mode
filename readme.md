# Cloud Parcel Model
<a href="https://zenodo.org/badge/latestdoi/460960645"><img src="https://zenodo.org/badge/460960645.svg" alt="DOI"></a>

## Description
This cloud parcel model simulates the updraft cooling of a cloud parcel with entrainment-mixing and entrained aerosols in Fortran 95. The initial development of this model was supported by the U.S. Department of Energy’s Earth System Modeling program via the Fast-physics System Testbed and Research project and Atmospheric System Research program through Stony Brook University and Brookhaven National Laboratory. The adiabatic version was developed by comparing with two other similar models by Graham Feingold/Huiwen Xue and Yiran Peng in 2014. The author would thank them for the kind sharing. The entrainment-mixing version was developed between 2015-2017.

## Author
  Jingyi Chen, PhD
  
  Research Scientist, Pacific Northwest National Laboratory
  
  Email: jingyi.jenna.chen@outlook.com
                
## Terms of Use
1. The parcel model code may be used for educational or non-profit purposes only. Any other usage must be first approved by the author.
2. The code cannot be modified in any way or form or distributed without the author's prior consent.
3. No portion of the source code can be used in other codes without the author's prior consent.
4. The code is provided on an as-is basis, and the author bears no liability from its usage.
5. Publications resulting from the usage of the code must cite at least one of the references below, as appropriate, for proper acknowledgment.
         
## References
1. Chen, J., Y. Liu, M. Zhang, and Y. Peng (2016), New understanding and quantification of the regime dependence of aerosol-cloud interaction for studying aerosol indirect effects, Geophysical Research Letter, 43, 1780–1787, doi:10.1002/2016GL067683.
2. Chen, J., Y. Liu, M. Zhang, and Y. Peng (2018), Height dependency of aerosol-cloud interaction regimes, Journal of Geophysical Research: Atmospheres, doi: 10.1002/2017JD027431.
3. Chen, J., Y. Liu and M. Zhang (2020), Effects of Lateral Entrainment-Mixing with Entrained Aerosols on Cloud Microphysics, Geophysical Research Letter, doi: 10.1029/2020GL087667.


## Code Structure
**maincpm.f90**
- temporary wrapper
- call chem to set the aerosol chemistry properties
- call aerospec to set aerosol size distribution
- call vertvelo to set the vertical velocity
- call entrain to set the entrianment mixing properties (defaul entrainment rate is 0)
- call dropspec to set the droplet size distribution
- call cpm to calculate the deliquesce/efflorescence and condensation/evaporation


**cpm.f90**
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

**cond.f90**
- time tendency equations

**DVODE_F90_M.f90**
- VODE(variable-coefficient ODE solver) is a general purpose solver written by P.N. Brown, G.D. Byrne and A.C. Hindmarsh from Lawrence Livermore National Laboratory ODEPACK suite. 
- version updated on Aug. 2013
- P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE, A Variable- Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10 (1989), pp. 1038-1051

