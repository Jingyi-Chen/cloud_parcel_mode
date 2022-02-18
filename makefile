# This is makefile 
#

# debug flag
 FFLAGS = 
 MyComp = ifort

# FFLAGS = -fopenmp
# MyComp = gfortran
	
exemain: constant.o empirical.o DVODE_F90_M.o chem.o  \
	 aerospec.o dropspec.o dynam.o deli.o cloudspec.o cond.o cpm.o maincpm.o
	$(MyComp) $(FFLAGS) -o exemain constant.o \
	empirical.o DVODE_F90_M.o \
	chem.o aerospec.o dropspec.o dynam.o deli.o cloudspec.o cond.o \
	cpm.o maincpm.o 


# compile the module without using other module 

constant.o: constant.f90
	$(MyComp) $(FFLAGS) -c constant.f90

empirical.o: empirical.f90
	$(MyComp) $(FFLAGS)  -c empirical.f90

chem.o: chem.f90
	$(MyComp)  $(FFLAGS) -c chem.f90

constant.mod: constant.o constant.f90
	$(MyComp)  $(FFLAGS) -c constant.f90

empirical.mod: empirical.o empirical.f90
	$(MyComp)  $(FFLAGS) -c empirical.f90

chem.mod: chem.o chem.f90
	$(MyComp)  $(FFLAGS) -c chem.f90

# 1st grade

aerospec.o: constant.mod aerospec.f90
	$(MyComp)  $(FFLAGS) -c aerospec.f90

dropspec.o: constant.mod chem.mod dropspec.f90
	$(MyComp)  $(FFLAGS) -c dropspec.f90

dynam.o: constant.mod empirical.mod dynam.f90
	$(MyComp)  $(FFLAGS) -c dynam.f90

deli.o: constant.mod empirical.mod chem.mod deli.f90
	$(MyComp)  $(FFLAGS) -c deli.f90

cloudspec.o: constant.mod cloudspec.f90
	$(MyComp)  $(FFLAGS) -c cloudspec.f90

aerospec.mod: aerospec.o aerospec.f90
	$(MyComp)  $(FFLAGS) -c aerospec.f90

dropspec.mod: dropspec.o chem.o dropspec.f90
	$(MyComp)  $(FFLAGS) -c dropspec.f90

dynam.mod: dynam.o dynam.f90
	$(MyComp)  $(FFLAGS) -c dynam.f90

deli.mod: deli.o deli.f90
	$(MyComp)  $(FFLAGS) -c deli.f90

cloudspec.mod: cloudspec.o cloudspec.f90
	$(MyComp)  $(FFLAGS) -c cloudspec.f90

# 2rd grade

cond.o: constant.mod empirical.mod chem.mod aerospec.mod dynam.mod cond.f90
	$(MyComp)  $(FFLAGS) -c cond.f90

cond.mod: cond.o cond.f90
	$(MyComp)  $(FFLAGS) -c cond.f90

# 3rd grade

DVODE_F90_M.o: DVODE_F90_M.f90
	$(MyComp)  $(FFLAGS) -c DVODE_F90_M.f90

DVODE_F90_M.mod: DVODE_F90_M.o DVODE_F90_M.f90
	$(MyComp)  $(FFLAGS) -c DVODE_F90_M.f90

cpm.o: constant.mod empirical.mod DVODE_F90_M.mod \
	deli.mod cond.mod cloudspec.mod cpm.f90
	$(MyComp)  $(FFLAGS) -c cpm.f90

cpm.mod: cpm.o cpm.f90
	$(MyComp)  $(FFLAGS) -c cpm.f90

# 4th grade

maincpm.o: cpm.mod aerospec.mod dropspec.mod chem.mod dynam.mod maincpm.f90
	$(MyComp)  $(FFLAGS) -c maincpm.f90


clean : 
	rm *.o *.mod
#END of the makefile

