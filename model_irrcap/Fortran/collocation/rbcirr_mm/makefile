include ./make.inc
#
# Make runsim_rbc and minpack library
#

MODULE = model_parameters.o wrapper_vars.o

RBC_FILES = rfuncirr.o cfunc.o invcfunc.o rfuncirr_wrapper.o rfunc.o rfunc_wrapper.o rfunc_wrapperls.o rfuncirr_wrapperls.o rfuncirr_mm_wrapperls.o rfuncirr_mm.o

HELPER_FILES = chebypol.o chebyroots.o rouwenhorst.o markovappr.o linspace.o normcdf.o sort.o constructpol.o dotproduct.o writemat_screen.o writevec_screen.o writevec.o

PROG_FILE = runsim_rbcirr.o

MINPACK_FILE = minpack.o

.SUFFIXES : .f90

all : minpack runsim_rbcirr 

runsim_rbcirr : 	$(MODULE) $(RBC_FILES) $(HELPER_FILES) $(PROG_FILE) $(DERIVS_FILES) 
		$(FC) -o runsim_rbcirr $(PROG_FILE) $(MODULE) $(RBC_FILES) $(HELPER_FILES) $(MINPACKLIB) -parallel

minpack :       $(MINPACK_FILE)
		ar rcs libminpack.a $(MINPACK_FILE)
		ranlib libminpack.a 

.f90.o:
	$(FC) $(OPTS) -c $<

.f.o:
	$(FC) $(OPTS) -c $<

clean : 
	rm -f *.o *~ runsim_rbcirr *.mod *.a
