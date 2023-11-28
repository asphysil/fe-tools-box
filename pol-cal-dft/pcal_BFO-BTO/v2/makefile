#!/usr/bin/bash

FC=gfortran
FCFLAG= -Wall -g
LIBPATH=
LIB= 
#${LIBPATH} -llapack -lblas
objects=polarizationdftv3.o

exe_frac.x:$(objects) 
	$(FC) ${FCFLAG} -o exe_frac.x $(objects) $(LIB) 

polarizationdftv3.o:polarizationdftv3.f90
	$(FC) ${FCFLAG} -c polarizationdftv3.f90 $(LIB) 
clean:
	rm $(objects)  exe_frac.x
