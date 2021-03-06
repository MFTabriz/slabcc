# Copyright (c) 2018, Bremen Center for Computational Materials Science (BCCMS), M. Farzalipour Tabriz
# Copyrights licensed under the 2-Clause BSD License.
# See the accompanying LICENSE.txt file for terms.


CXX =icpc
#CXX =g++

# Please edit/choose the path to one of the available
# libraries on your system. 
# Please read the “Validation” section in the manual/readme file for more information 
# about the compiler/libraries in our test environment.
##===========EXTERNAL LIBRARIES============##

#FFTW_INC_PATH = -I${FFTW3_HOME}/include/
#FFTW_LIB_PATH = -L${FFTW3_HOME}/lib/

#FFTW_INC_PATH = -I/cluster/fftw/3.3.6p2/intel2016/include/
#FFTW_LIB_PATH = -L/cluster/fftw/3.3.6p2/intel2016/lib/

#FFTW_LIB = -lfftw3

#BLAS_INC_PATH = -I/cluster/OpenBLAS/0.2.19/gcc62/include/
#BLAS_LIB_PATH = -L/cluster/OpenBLAS/0.2.19/gcc62/lib/

#BLAS_INC_PATH = -I/cluster/OpenBLAS/1f217/intel2016/include/
#BLAS_LIB_PATH = -L/cluster/OpenBLAS/1f217/intel2016/lib/

FFTW_INC_PATH = -I${MKLROOT}/include/fftw
BLAS_INC_PATH = -I${MKLROOT}/include
BLAS_LIB_PATH = -L${MKLROOT}/lib/intel64

#BLAS_LIB = -lblas -llapack #OpenBLAS
BLAS_LIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -ldl #MKL_INTEL
#BLAS_LIB =-lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lpthread -ldl #MKL_GCC

##============COMPILER FLAGS=============##
OPT      = -O2
MARCH    = -march=native #change this if the architecture of the machine compiling and running the program is not the same
EXTRAOPT = #add yours own compiler optimization flags here!
OPTFLAGS = $(OPT) $(MARCH) $(EXTRAOPT)


ifeq ($(CXX),icpc)
  OMP = -qopenmp
  CC = icc
  EXTRA_LIBS = -lifcore $(OMP)
endif

ifeq ($(CXX),g++)
  OMP = -fopenmp
  CC  = gcc
  EXTRA_LIBS = -lgfortran -lquadmath $(OMP)
endif

ALLFLAGS = $(OMP) $(OPTFLAGS) # these flags and the compiler are also used by the NLOPT configure/makefile

CFLAGS   = $(ALLFLAGS) 
CXXFLAGS = $(ALLFLAGS) -std=c++14 #-Wall

##included libraries
NLOPT_PATH = ../src/nlopt
NLOPT_INC_PATH = -I$(NLOPT_PATH)/include/
NLOPT_LIB_PATH = -L$(NLOPT_PATH)/lib/
NLOPT_LIB_FILE = $(NLOPT_PATH)/libnlopt.a
NLOPT_LIB = -lnlopt

#Add the -DMKL if you are linking against it!
CPP_DEFS = -DARMA_NO_DEBUG -DARMA_DONT_USE_WRAPPER #-DMKL 

LDLIBS = $(NLOPT_LIB) $(BLAS_LIB) $(FFTW_LIB) $(EXTRA_LIBS)
LIB_PATHS = $(NLOPT_LIB_PATH) $(FFTW_LIB_PATH) $(BLAS_LIB_PATH)
#LD_EXTRA_FLAGS = -Wl,--verbose

SOURCE_INC_PATHS = -I../src/ -I../src/armadillo/include/ -I../src/inih/cpp/ -I../src/clara/single_include/ -I../src/spline/ -I../src/spdlog/
CPPFLAGS = $(CPP_DEFS) $(SOURCE_INC_PATHS) $(NLOPT_INC_PATH) $(FFTW_INC_PATH) $(BLAS_INC_PATH)

SOURCES = general_io.cpp slabcc_math.cpp vasp.cpp slabcc.cpp stdafx.cpp slabcc_model.cpp slabcc_input.cpp ini.c INIReader.cpp madelung.cpp isolated.cpp
OBJECTS = $(patsubst %.c,%.o,$(SOURCES:.cpp=.o))
EXECUTABLE = slabcc

vpath %.cpp ../src:../src/inih/cpp
vpath %.c ../src/inih

all: $(SOURCES) $(NLOPT_LIB_FILE) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(LIB_PATHS) $(LD_EXTRA_FLAGS) $(OBJECTS) $(LDLIBS) -o $@
	rm -f $(OBJECTS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -c

##compile the nlopt
#.ONESHELL: #not supported in GNU Make 3.81 and earlier!
$(NLOPT_LIB_FILE):
	cd ../src/nlopt/nlopt-2.4.2;\
	chmod +x configure;\
	./configure CC=$(CC) CXX=$(CXX) CFLAGS='$(ALLFLAGS)' CXXFLAGS='$(ALLFLAGS)' --prefix=$$(cd ../ && pwd);\
	make;\
	make install

.PHONY : clean distclean

clean :
	rm -f $(OBJECTS) $(NLOPT_LIB_FILE)

distclean: clean
	rm -fr $(NLOPT_PATH)/include $(NLOPT_PATH)/lib $(NLOPT_PATH)/share
	make distclean -C ../src/nlopt/nlopt-2.4.2/
