# Copyright (c) 2018, Bremen Center for Computational Materials Science (BCCMS), M. Farzalipour Tabriz
# Copyrights licensed under the 2-Clause BSD License.
# See the accompanying LICENSE.txt file for terms.


# Please edit/choose the path to one of the available 
# libraries on your system. The code have been tested 
# using g++/icpc and OpenBLAS/MKL
##===========Path to the external libraries============##

#FFTW_INC_PATH = -I/cluster/fftw/3.3.6p2/gcc63/include
#FFTW_LIB_PATH = -L/cluster/fftw/3.3.6p2/gcc63/lib/

FFTW_INC_PATH = -I/cluster/fftw/3.3.6p2/intel2016/include/
FFTW_LIB_PATH = -L/cluster/fftw/3.3.6p2/intel2016/lib/

FFTW_LIB = -lfftw3

#BLAS_INC_PATH = -I/cluster/OpenBLAS/0.2.19/gcc62/include/
#BLAS_LIB_PATH = -L/cluster/OpenBLAS/0.2.19/gcc62/lib/

BLAS_INC_PATH = -I/cluster/OpenBLAS/1f217/intel2016/include/
BLAS_LIB_PATH = -L/cluster/OpenBLAS/1f217/intel2016/lib/

#BLAS_INC_PATH =  -I${MKLROOT}/include
#BLAS_LIB_PATH =  -L${MKLROOT}/lib/intel64

BLAS_LIB = -lopenblas_omp
#BLAS_LIB = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm -ldl


##Compiler
#CXX = g++
CXX = icpc

CXXFLAGS = -std=c++14 -O2
##=====================================================##



# To compile or make static linking with other compilers, 
# you must include your own switches here!
# These are pre-set for icpc/g++ with Openblas
##===========Compiler specific options=================##
ifeq ($(CXX), icpc)
	EXTRA_LIBS = -lifcore -qopenmp
	ifdef STATIC
		LDSTATIC = -static-intel -Wl,-Bstatic
	endif
endif

ifeq ($(CXX), g++)
	EXTRA_LIBS = -lgfortran -lquadmath -fopenmp
	ifdef STATIC
		LDSTATIC = -static-libgcc -static-libstdc++ -Wl,-Bstatic
	endif
endif
##=====================================================##


##included libraries
NLOPT_PATH = ../src/nlopt
NLOPT_INC_PATH = -I../src/nlopt/include/
NLOPT_LIB_PATH = -L../src/nlopt/lib/
NLOPT_LIB_FILE = ../src/nlopt/lib/libnlopt.a
NLOPT_LIB = -lnlopt

CPP_DEFS = -DARMA_NO_DEBUG -DARMA_DONT_USE_WRAPPER

LDLIBS = $(NLOPT_LIB) $(BLAS_LIB) $(FFTW_LIB) $(EXTRA_LIBS)
LIB_PATHS = $(NLOPT_LIB_PATH) $(FFTW_LIB_PATH) $(BLAS_LIB_PATH)
#LD_EXTRA_FLAGS = -Wl,--verbose

SOURCE_INC_PATHS = -I../src/ -I../src/armadillo/include/ -I../src/inih/cpp/ -I../src/clara/single_include/ -I../src/spline/
CPPFLAGS = $(CPP_DEFS) $(SOURCE_INC_PATHS) $(NLOPT_INC_PATH) $(FFTW_INC_PATH) $(BLAS_INC_PATH)

SOURCES = general_io.cpp slabcc_math.cpp vasp.cpp slabcc.cpp stdafx.cpp slabcc_core.cpp ini.c INIReader.cpp madelung.cpp
OBJECTS = $(patsubst %.c,%.o,$(SOURCES:.cpp=.o))
EXECUTABLE = slabcc


vpath %.cpp ../src:../src/inih/cpp
vpath %.c ../src/inih

all: $(SOURCES) $(NLOPT_LIB_FILE) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LIB_PATHS) $(LD_EXTRA_FLAGS) $(OBJECTS) $(LDSTATIC) $(LDLIBS) -o $@
	rm -f $(OBJECTS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< -c

##unzip and compile the nlopt
##.ONESHELL: #not supported in GNU Make 3.81 and earlier!
$(NLOPT_LIB_FILE):
	cd ../src/nlopt/nlopt-2.4.2;\
	chmod +x configure;\
	./configure --prefix=$$(cd ../ && pwd);\
	make;\
	make install
	
.PHONY : clean distclean

clean :
	rm -f $(OBJECTS) $(EXECUTABLE) $(NLOPT_LIB_FILE)
	
distclean: clean
	rm -fr $(NLOPT_PATH)/include $(NLOPT_PATH)/lib $(NLOPT_PATH)/share	
	make distclean -C ../src/nlopt/nlopt-2.4.2/ 


