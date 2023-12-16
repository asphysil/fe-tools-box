#!/usr/bin/bash

FC=gfortran
FCFLAG=
FFTLIB= # -L /usr/local/lib/libfftw3.a  # directory
# directory is the directory where this library is stored
# a file with the name libfftw3.a
FFTFLAG=  #-lfftw3
objects=ABO3_dtype.o ABO3_angle.o ABO3_domain.o ABO3_main.o

domain.x:$(objects) 
	$(FC) -o domain.x $(objects) $(FFTLIB) $(FFTFLAG)

ABO3_dtype.o: ABO3_dtype.f90
	$(FC) -c ABO3_dtype.f90 $(FFTLIB) $(FFTFLAG)

ABO3_angle.o: ABO3_angle.f90
	$(FC) -c ABO3_angle.f90 $(FFTLIB) $(FFTFLAG)
ABO3_domain.o:ABO3_domain.f90
	$(FC) -c ABO3_domain.f90 $(FFTLIB) $(FFTFLAG)
ABO3_main.o:ABO3_main.f90
	$(FC) -c ABO3_main.f90 $(FFTLIB) $(FFTFLAG)

clean:
	rm $(objects)  domain.x *.mod*
