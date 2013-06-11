include ./make.inc
#
# Make runsim_rbc
#

UTIL_FILES = rouwenhorst.o linspace.o

HELPER_FILES = loaddata.o writemat.o writemat_int.o 

PROG_FILE = runsim_rbc.o 

.SUFFIXES : .f90

all : runsim_rbc 

runsim_rbc : 	$(UTIL_FILES) $(HELPER_FILES) $(PROG_FILE)
		$(FC) -o runsim_rbc $(PROG_FILE) $(UTIL_FILES) $(HELPER_FILES) -m64 -openmp -openmp-report


.f90.o:
	$(FC) $(OPTS) -c $<

.f.o:
	$(FC) $(OPTS) -c $<

clean : 
	rm -f *.o *~ runsim_skander *.mod *.a *__genmod.f90
