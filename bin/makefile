# Copyright (c) 2018-2023, Bremen Center for Computational Materials Science (BCCMS), M. Farzalipour Tabriz
# Copyrights licensed under the 2-Clause BSD License.
# See the accompanying LICENSE.txt file for terms.

# todo: debug option

# these variable are read from the environment:
# $CC: C compiler (default: gcc)
# $CXX: C++ compiler (default: g++)
# $FFTW_HOME: path to FFTW library home
# $FFTW_LIB_FLAGS: FFTW library flag (default: -lfftw3)
# $BLAS_HOME path to BLAS library home
# $BLAS_LIB_FLAGS: BLAS library flags (default: -lblas -llapack -lpthread)
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
TESTS_PATH := ../tests
BUILD_INFO_PATH = ../src
UTILS_PATH := ../utils

## EXTERNAL LIBRARIES
## 
NLOPT_PREFIX := nlopt


##=================COMPILER FLAGS==============##

ifeq ($(VERBOSE),1)
  LD_EXTRA_FLAGS := $(LD_EXTRA_FLAGS) -Wl,--verbose
endif

ifeq ($(MKL),1)
  DMKL := -DMKL

  FFTW_INC_PATH := ${MKLROOT}/include/fftw
  BLAS_INC_PATH := ${MKLROOT}/include
  BLAS_LIB_PATH := ${MKLROOT}/lib/intel64

  # Note: compiler executables may have names such as gcc-5
  ifneq (,$(findstring icpc,$(CXX)))
    BLAS_LIB_FLAGS ?= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
  endif
  ifneq (,$(findstring icpx,$(CXX)))
    BLAS_LIB_FLAGS ?= -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
  endif
  ifneq (,$(findstring g++,$(CXX)))
    BLAS_LIB_FLAGS ?=  -Wl,--no-as-needed -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl
  endif

else

  #OpenBLAS
  BLAS_LIB_FLAGS ?= -lblas -llapack -lpthread


  ifdef $(BLAS_HOME)
    BLAS_INC_PATH := ${BLAS_HOME}/include/
    BLAS_LIB_PATH := ${BLAS_HOME}/lib/
  endif

  # FFTW library
  FFTW_LIB_FLAGS ?= -lfftw3


  ifdef $(FFTW_HOME)
    FFTW_INC_PATH := ${FFTW_HOME}/include/
    FFTW_LIB_PATH := ${FFTW_HOME}/lib/
  endif

endif


# these flags are also used by the NLOPT configure/makefile
COMMON_FLAGS = $(strip $(OMP) $(OPT) $(EXTRA_FLAGS)) -g

ifneq (,$(findstring clang,$(CXX)))
  WARN_FLAGS := -Wall -Wextra -Wno-unused-parameter -Wno-tautological-constant-compare -Werror
else ifneq (,$(findstring icpx,$(CXX)))
  WARN_FLAGS := -Wall -Wextra -Wno-unused-parameter -Wno-tautological-constant-compare -Werror
else
  WARN_FLAGS := -Wall -Wextra -Wno-unused-parameter -Werror
endif

ifneq (,$(findstring ic,$(CXX)))
  OMP := -qopenmp
endif

ifneq (,$(findstring g++,$(CXX)))
  OMP := -fopenmp
endif

LD_EXTRA_FLAGS :=  $(strip $(LD_EXTRA_FLAGS) $(OMP))

CFLAGS   = $(COMMON_FLAGS) 
CXXFLAGS = $(COMMON_FLAGS) -std=$(CPP_STD) $(WARN_FLAGS)

##included libraries
NLOPT_SRC_PATH := $(shell realpath ../src/nlopt)

NLOPT_PATH = $(shell realpath $(NLOPT_PREFIX))
NLOPT_INC_PATH = $(NLOPT_PREFIX)/include
NLOPT_LIB_PATH = $(NLOPT_PREFIX)/lib
NLOPT_LIB_FILE = $(NLOPT_LIB_PATH)/libnlopt.a
NLOPT_HEADER_FILES = $(NLOPT_INC_PATH)/nlopt.hpp
NLOPT_LIB_FLAGS := -lnlopt

# preprocessor flags
CPP_DEFS = $(strip -DARMA_NO_DEBUG -DARMA_DONT_USE_WRAPPER $(DMKL))

LIB_PATH_FLAGS = $(addprefix -L,$(strip $(NLOPT_LIB_PATH) $(FFTW_LIB_PATH) $(BLAS_LIB_PATH)))
LD_LIB_FLAGS := $(strip $(NLOPT_LIB_FLAGS) $(BLAS_LIB_FLAGS) $(FFTW_LIB_FLAGS))

SOURCE_INC_PATHS := $(addprefix ../src/,$(strip . armadillo/include inih/cpp clara/single_include spline spdlog))
INC_FLAGS = $(addprefix -I,$(strip $(SOURCE_INC_PATHS) $(NLOPT_INC_PATH) $(FFTW_INC_PATH) $(BLAS_INC_PATH)))

CPPFLAGS = $(CPP_DEFS) $(INC_FLAGS)

##====================SOURCES=====================##
# generate source file paths from these directories
vpath %.cpp ../src:../src/inih/cpp
vpath %.c ../src/inih

SOURCES := general_io.cpp slabcc_math.cpp vasp.cpp slabcc.cpp slabcc_model.cpp slabcc_input.cpp ini.c INIReader.cpp madelung.cpp isolated.cpp
BUILD_INFO := $(addprefix $(BUILD_INFO_PATH)/,compilerinfo.inc gitinfo.inc)
OBJECTS := $(addprefix $(OBJ_PATH)/,$(patsubst %.c,%.o,$(SOURCES:.cpp=.o)))
LIBRARIES = $(NLOPT_LIB_FILE)
LIBRARY_HEADERS = $(NLOPT_HEADER_FILES)
EXECUTABLE := slabcc

##====================BUILD INFO=====================##
COMPILER_VER := $(shell $(CXX) --version | head -n 1)
COMPILER_INFO := "const char *cxx_compiler=\"$(CXX)\", *cxx_compilerver=\"$(COMPILER_VER)\", *cxx_compilerflags=\"$(CXXFLAGS)\";"
COMMIT_SHA := $(shell $(UTILS_PATH)/git_commit.sh)
COMMIT_INFO := "const char * git_info=\"Git commit $(COMMIT_SHA)\";"

##===================================================##

.PHONY: all clean distclean test man force

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) $(LIBRARIES)
	$(CXX) $(OBJECTS) $(LIB_PATH_FLAGS) $(LD_LIB_FLAGS) $(LD_EXTRA_FLAGS) -o $@

$(OBJECTS): $(LIBRARY_HEADERS) $(BUILD_INFO) $(SOURCES)

$(OBJ_PATH)/%.o: %.cpp
	@mkdir -p $(shell dirname $@)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

$(OBJ_PATH)/%.o: %.c
	@mkdir -p $(shell dirname $@)
	$(CC) $(CFLAGS) $(CPPFLAGS) -c -o $@ $<

##compile nlopt
$(NLOPT_LIB_FILE) $(NLOPT_HEADER_FILES):
	mkdir -p $(NLOPT_PREFIX)
	cd $(NLOPT_PREFIX) && $(NLOPT_SRC_PATH)/configure CC=$(CC) CXX=$(CXX) CFLAGS='$(COMMON_FLAGS)' CXXFLAGS='$(COMMON_FLAGS)' --prefix=$(NLOPT_PATH) || cat config.log;
	make --directory=$(NLOPT_PATH)
	make --directory=$(NLOPT_PATH) install

$(BUILD_INFO_PATH)/compilerinfo.inc: force
	echo $(COMPILER_INFO) | cmp -s - $@ || echo $(COMPILER_INFO) > $@

$(BUILD_INFO_PATH)/gitinfo.inc: force
	echo $(COMMIT_INFO) | cmp -s - $@ || echo $(COMMIT_INFO) > $@

test: $(EXECUTABLE)
	$(MAKE) -C $(TESTS_PATH)

man:
	cd .. && utils/generate_docs.sh

clean:
	rm -fr $(OBJECTS) $(EXECUTABLE)
	rm -f $(BUILD_INFO_PATH)/*.inc
	$(MAKE) -C $(TESTS_PATH) clean

distclean: clean
	rm -fr $(NLOPT_PREFIX)
