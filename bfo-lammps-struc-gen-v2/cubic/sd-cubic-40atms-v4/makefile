#!/usr/bin/bash

FC=gfortran
#mpif90
# for debug
FCFLAG = -Wall -Wextra -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
#-Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace
# for production run
#FCFLAG = -Wall -Wextra -Wimplicit-interface -fPIC -fmax-errors=1 -O3 -march=native -ffast-math -funroll-loops

FFTLIB=-I /opt/cray/pe/fftw/3.3.10.5/x86_milan/
#/global/share/mpi/mpich-ofi-gnu9-ohpc/3.4.2/lib/ 
# -L /usr/local/lib/libfftw3.a  # directory
# directory is the directory where this library is stored
# a file with the name libfftw3.a
FFTFLAG=  #-lfftw3
objects=ABO3_dtype.o ABO3_struc.o\
       ABO3_AFM.o ABO3_axis_dm.o\
       ABO3_domain.o\
      ABO3_output_vasp.o\
      ABO3_output_lmp_jiahao.o ABO3_output_lmp_shi.o\
       ABO3_main.o
       

#ABO3_angle.o ABO3_domain.o ABO3_main.o

domain.x:$(objects) 
	$(FC) ${FCFLAG} -o domain.x $(objects) $(FFTLIB) $(FFTFLAG)

ABO3_dtype.o: ABO3_dtype.f90
	$(FC)  ${FCFLAG} -c ABO3_dtype.f90 $(FFTLIB) $(FFTFLAG)

ABO3_axis_dm.o: ABO3_axis_dm.f90
	${FC}  ${FCFLAG} -c ABO3_axis_dm.f90 $(FFTLIB) $(FFTFLAG)

#ABO3_angle.o: ABO3_angle.f90
#	$(FC) -c ABO3_angle.f90 $(FFTLIB) $(FFTFLAG)
ABO3_domain.o:ABO3_domain.f90
	$(FC)  ${FCFLAG} -c ABO3_domain.f90 $(FFTLIB) $(FFTFLAG)

ABO3_main.o:ABO3_main.f90
	$(FC)  ${FCFLAG} -c ABO3_main.f90 $(FFTLIB) $(FFTFLAG)

ABO3_struc.o: ABO3_struc.f90
	${FC}  ${FCFLAG} -c  ABO3_struc.f90 $(FFTLIB) $(FFTFLAG)

ABO3_AFM.o: ABO3_AFM.f90
	${FC}  ${FCFLAG} -c ABO3_AFM.f90 $(FFTLIB) $(FFTFLAG)

ABO3_output_vasp.o: ABO3_output_vasp.f90
	${FC}  ${FCFLAG} -c ABO3_output_vasp.f90 $(FFTLIB) $(FFTFLAG)
#
ABO3_output_lmp_jiahao.o: ABO3_output_lmp_jiahao.f90
	${FC}  ${FCFLAG} -c ABO3_output_lmp_jiahao.f90 $(FFTLIB) $(FFTFLAG)
ABO3_output_lmp_shi.o: ABO3_output_lmp_shi.f90
	${FC}  ${FCFLAG} -c ABO3_output_lmp_shi.f90 $(FFTLIB) $(FFTFLAG)
clean:
	rm $(objects)  domain.x *.mod*
