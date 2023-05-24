# Copyright (c) 2018-2023, Bremen Center for Computational Materials Science (BCCMS), M. Farzalipour Tabriz
# Copyrights licensed under the 2-Clause BSD License.
# See the accompanying LICENSE.txt file for terms.

# todo: debug option
# test mkl

# these variable are read from the environment:
# $CC: C compiler (default: gcc)
# $CXX: C++ compiler (default: g++)
# $FFTW_HOME: path to FFTW library home
# $FFTW_LIB: FFTW library flag (default: -lfftw3)
# $BLAS_HOME path to BLAS library home
# $BLAS_LIB: BLAS library flags (default: -lblas -llapack -lpthread)
# $EXTRA_FLAGS: extra compiler flags for CC and CXX
# $LD_EXTRA_FLAGS: extra linker flags

# Please read the “Validation” section in the manual/readme file for more information 
# on the compiler/libraries in our test environment.

##============ DEFAULT PARAMETERS =============##
## 
## COMPILER
## 

CC ?= gcc
CXX ?= g++

OPT := -O2
CPP_STD := c++14

# compilation objects are put here
OBJ_PATH := slabcc_obj
## 
## EXTERNAL LIBRARIES
## 
NLOPT_PREFIX := nlopt


ifeq ($(VERBOSE),1)
  LD_EXTRA_FLAGS := $(LD_EXTRA_FLAGS) -Wl,--verbose
endif

ifeq ($(MKL),1)
  DMKL := -DMKL

  FFTW_INC_PATH := ${MKLROOT}/include/fftw
  BLAS_INC_PATH := ${MKLROOT}/include
  BLAS_LIB_PATH := ${MKLROOT}/lib/intel64


  ifneq (,$(findstring icpc,$(CXX)))
    BLAS_LIB ?= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
  endif
  ifneq (,$(findstring icpx,$(CXX)))
    BLAS_LIB ?= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
  endif
  ifneq (,$(findstring g++,$(CXX)))
    BLAS_LIB ?=  -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
  endif

else

  #OpenBLAS
  BLAS_LIB ?= -lblas -llapack -lpthread


  ifdef $(BLAS_HOME)
    BLAS_INC_PATH := ${BLAS_HOME}/include/
    BLAS_LIB_PATH := ${BLAS_HOME}/lib/
  endif

  # FFTW library
  FFTW_LIB ?= -lfftw3


  ifdef $(FFTW_HOME)
    FFTW_INC_PATH := ${FFTW_HOME}/include/
    FFTW_LIB_PATH := ${FFTW_HOME}/lib/
  endif

endif

##================================================##
# these flags are also used by the NLOPT configure/makefile
COMMON_FLAGS = $(OMP) $(OPT) $(EXTRA_FLAGS) 

ifeq ($(CXX),icpc)
  OMP := -qopenmp
  LD_EXTRA_FLAGS := $(LD_EXTRA_FLAGS) -lifcore $(OMP)
endif

ifeq ($(CXX),g++)
  OMP := -fopenmp
  LD_EXTRA_FLAGS := $(LD_EXTRA_FLAGS) -lgfortran -lquadmath $(OMP)
endif

CFLAGS   = $(COMMON_FLAGS) 
CXXFLAGS = $(COMMON_FLAGS) -std=$(CPP_STD)# -Wall

##included libraries
NLOPT_SRC_PATH := $(shell realpath ../src/nlopt/nlopt-2.4.2)

NLOPT_PATH = $(shell realpath $(NLOPT_PREFIX))
NLOPT_INC_PATH = $(NLOPT_PREFIX)/include
NLOPT_LIB_PATH = $(NLOPT_PREFIX)/lib
NLOPT_LIB_FILE = $(NLOPT_LIB_PATH)/libnlopt.a
NLOPT_LIB := -lnlopt

CPP_DEFS = -DARMA_NO_DEBUG -DARMA_DONT_USE_WRAPPER $(DMKL)

LIB_PATH_FLAGS = $(addprefix -L,$(strip $(NLOPT_LIB_PATH) $(FFTW_LIB_PATH) $(BLAS_LIB_PATH)))
LD_LIBS := $(NLOPT_LIB) $(BLAS_LIB) $(FFTW_LIB) $(EXTRA_LIBS)

SOURCE_INC_PATHS := $(addprefix ../src/,$(strip . armadillo/include inih/cpp clara/single_include spline spdlog))
INC_FLAGS = $(addprefix -I,$(strip $(SOURCE_INC_PATHS) $(NLOPT_INC_PATH) $(FFTW_INC_PATH) $(BLAS_INC_PATH)))

CPPFLAGS = $(CPP_DEFS) $(INC_FLAGS)

# generate source file paths from these directories
vpath %.cpp ../src:../src/inih/cpp
vpath %.c ../src/inih

SOURCES := general_io.cpp slabcc_math.cpp vasp.cpp slabcc.cpp stdafx.cpp slabcc_model.cpp slabcc_input.cpp ini.c INIReader.cpp madelung.cpp isolated.cpp
OBJECTS := $(addprefix $(OBJ_PATH)/,$(patsubst %.c,%.o,$(SOURCES:.cpp=.o)))
LIBRARIES = $(NLOPT_LIB_FILE)
EXECUTABLE := slabcc


all: $(SOURCES) $(LIBRARIES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIB_PATH_FLAGS) $(LD_LIBS) $(LD_EXTRA_FLAGS) -o $@

$(OBJ_PATH)/%.o: %.cpp
	@mkdir -p $(shell dirname $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

$(OBJ_PATH)/%.o: %.c
	@mkdir -p $(shell dirname $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

##compile nlopt
$(NLOPT_LIB_FILE):
	mkdir -p $(NLOPT_PREFIX)
	cd $(NLOPT_PREFIX) && $(NLOPT_SRC_PATH)/configure CC=$(CC) CXX=$(CXX) CFLAGS='$(COMMON_FLAGS)' CXXFLAGS='$(COMMON_FLAGS)' --prefix=$(NLOPT_PATH) || cat $(NLOPT_PREFIX)/config.log;
	make --directory=$(NLOPT_PATH)
	make --directory=$(NLOPT_PATH) install

.PHONY: clean distclean

clean:
	rm -fr $(OBJECTS) $(EXECUTABLE)

distclean: clean
	rm -fr $(NLOPT_PREFIX)
