# AMREX_HOME defines the directory in which we will find all the AMReX code
AMREX_HOME ?= ../../../amrex

CVODE_LIB_DIR ?= $(CVODE_LIB)

# TOP defines the directory in which we will find Source, Exec, etc
TOP = ../..

# compilation options
COMP    = pgi
USE_MPI = TRUE
USE_OMP = FALSE
USE_CUDA = TRUE

PROFILE       = FALSE
TRACE_PROFILE = FALSE
COMM_PROFILE  = FALSE
TINY_PROFILE  = TRUE

PRECISION = DOUBLE
USE_SINGLE_PRECISION_PARTICLES = TRUE
DEBUG     = FALSE

GIMLET = FALSE
REEBER = FALSE

# physics
DIM      = 3
USE_GRAV = TRUE
USE_HEATCOOL = TRUE
USE_AGN = FALSE
USE_CVODE = FALSE

USE_CUVODE = TRUE
CUVODE_HOME ?= $(AMREX_HOME)/../Microphysics/integration/VODE90/cuVODE/
#include ../../Source/HeatCool/Make.cuVODE

Bpack := ./Make.package
Blocs := .

Bdirs := source linear_algebra_modules
Bpack += $(foreach dir, $(Bdirs), $(CUVODE_HOME)/$(dir)/Make.package)
Blocs += $(foreach dir, $(Bdirs), $(CUVODE_HOME)/$(dir))


TOP = ../..

include $(TOP)/Exec/Make.Nyx

