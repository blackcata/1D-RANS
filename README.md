# 1D-RANS
1D RANS model simulation at fully developed turbulent channel flow.

This code is programmed for 1D RANS simulation using k-epsilon model to fully developed channel flows. This code is parallelized by OpenMP that you have to set the number of threads.

Various damping function is applied for 1D RANS.
  - Van Driest (1954)
  - Launder and Sharma (1974)
  - Lam and Bremhorst (1981)
  - Park et al (1997)
  
These simulation datas are validated with DNS results (available in http://turbulence.ices.utexas.edu/)

### 0. Related papers & data 
  - Related theory : https://www.dropbox.com/s/fxk5e604x3v1yua/4.RANS_model_wallmodel.pdf?dl=0
  - Related paper : https://www.dropbox.com/s/7ik7s9a8zik6u4o/Turbulence%20Modeling%20-%20Project1.pdf?dl=0
  - Related report : https://www.dropbox.com/s/h8er705m7pk2lsb/Report.pdf?dl=0

### 1. Setting for 1D RANS code
  - Make 'RESULT' folder 
  - Set the channel-half height
  - Set Re_tau and kinetic viscosity
  - Set damping function type
  - Adjust relaxation factors
  
### 2. Output files
  - Velocity profile ( U.plt )
  - Turbulent kinetic energy ( k.plt )
  - Dissipation ( dissipation.plt )
  - Damping function ( fm.plt )
  
### 3. Code composition
  - RANS_main.f90
  - RANS_module.f90
  - RANS_setup.f90
  - RANS_poiseuille.f90
  - RANS_getfm.f90
  - RANS_getnut.f90
  - RANS_getu.f90
  - RANS_getprod.f90
  - RANS_getk.f90
  - RANS_getdis.f90
  - RANS_output.f90
