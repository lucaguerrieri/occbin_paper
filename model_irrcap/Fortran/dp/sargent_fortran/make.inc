SHELL = /bin/sh

#
# Fortran Compilers
#

# to link to NAG use f95
FC = ifort
#
# Fortran Module Compiler
#
M2C = $(FC)

#
# Options
#
OPTS	   = -m64 -fast -openmp

# Debugging options for ifort compiler
#OPTS = -debug -check -traceback -gen-interfaces
#

#
# Path to MKL libraries
#
MKLPATH    = /Developer/opt/intel/Compiler/11.1/080/Frameworks/mkl/lib/em64t
MKLPATH2   = /Developer/opt/intel/Compiler/11.1/080/lib

#MINPACKPATH = /ofs/research2/Guerrieri/projection_skander
MINPACKPATH = /Users/Jason/Programs/projection_skander

NAGPATH = /opt/NAG/fll6a21dfl/lib
#
# Loader
#
LOADER   = ifort
LOADOPTS = 

#
# AR, RANLIB
#
ARCH		= ar
ARCHFLAGS	= -cr
RANLIB		= ranlib

#
# Path to Libraries
#
#BLAS		= -L/users/Jason/programs/libraries/lapack-3.2.1/ -lblas_gfortran
#LAPACK		= -L/users/Jason/programs/libraries/lapack-3.2.1/ -llapack_gfortran
#AIM 		= -L/users/Jason/programs/test2_gfortran/LIBAIM2/ -laim
