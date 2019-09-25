# AMREX_HOME defines the directory in which we will find all the BoxLib code
# If you set AMREX_HOME as an environment variable, this line will be ignored
AMREX_HOME ?= ../../../amrex

DEBUG            = FALSE
USE_MPI          = FALSE #True should also work here
USE_OMP          = FALSE
USE_CUDA         = TRUE
COMP             = pgi # gnu should also work here
DIM              = 3
PRECISION        = DOUBLE
EBASE            = main

CVODE_LIB_DIR ?= ${CVODE_LIB}
USE_SUNDIALS_3x4x = TRUE

#For consistency testing:
USE_CUDA=FALSE 
USE_HEATCOOL=TRUE
USE_ARKODE_LIBS=TRUE

ifeq ($(USE_HEATCOOL), TRUE)
  DEFINES += -DHEATCOOL
endif

include $(AMREX_HOME)/Tools/GNUMake/Make.defs

include ./Make.package
include $(AMREX_HOME)/Src/Base/Make.package

#ifeq (($USE_CUDA),FALSE)
include ../VODE/Make.package
include ../BLAS/Make.package
#endif

CUVODE_HOME=$(AMREX_HOME)/../Microphysics/integration/VODE90/cuVODE/

#Blocs += $(AMREX_HOME)/../Microphysics/integration/VODE90/cuVODE/source

#include $(AMREX_HOME)/../Microphysics/integration/VODE90/cuVODE/Make.cuVODE
include Make.cuVODE

#ifeq (($USE_CUDA),FALSE)
INCLUDE_LOCATIONS += $(Blocs) ../BLAS ../VODE
VPATH_LOCATIONS   += $(Blocs) ../BLAS ../VODE
#endif

include $(AMREX_HOME)/Tools/GNUMake/Make.rules



