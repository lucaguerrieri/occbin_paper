SHELL = /bin/sh

#
# Fortran Compilers
#
FC = ifort
#
# Fortran Module Compiler
#
M2C = $(FC)

#
# Options
#
OPTS	   = -m64 -fast -parallel
#OPTS = -debug -check -traceback -gen-interfaces

#
# Path to MKL libraries
#
MKLPATH    = /Developer/opt/intel/Compiler/11.1/080/Frameworks/mkl/lib/em64t
MKLPATH2   = /Developer/opt/intel/Compiler/11.1/080/lib

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

MINPACKLIB 	=-L/Users/Jason/Programs/projection_collection/rbc_collocation/ -lminpack

