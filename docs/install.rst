============
Installation
============
1. **Prerequisites:**

 #. **Source:** Download the `latest stable release <https://codeberg.org/meisam/slabcc/releases>`__ and extract the files. You can also clone the `git repository <https://codeberg.org/meisam/slabcc.git>`__ and use the latest commit on the master branch to get the latest changes.
 #. **Compiler:** You need a C++ compiler with `C++14 standard support <https://en.cppreference.com/w/cpp/compiler_support#C.2B.2B14_features>`_ (e.g. `g++ <https://gcc.gnu.org/>`_ 5.0 or later)
 #. **BLAS/OpenBLAS/MKL:** You can use BLAS+LAPACK for the matrix operations inside the slabcc but it is highly recommended to use one of the high performance replacements, e.g., the `OpenBLAS <https://github.com/xianyi/OpenBLAS/releases>`_/`MKL <https://software.intel.com/en-us/mkl>`_ instead. If you don't have OpenBLAS installed on your system, follow the guide on the `OpenBLAS website <https://www.openblas.net>`_. Please refer to the `Armadillo documentation <https://gitlab.com/conradsnicta/armadillo-code/-/blob/9.900.x/README.md>`_ for linking to other BLAS replacements.
 #. **FFTW:** If you don't have FFTW installed on your system, follow the guide on the `FFTW website <https://www.fftw.org/download.html>`_. Alternatively, you can use the FFTW interface of the MKL.

2. **Configuration:** Set compilation parameters through environment variables.

 #. **$CC:** C compiler (default: gcc)
 #. **$CXX:** C++ compiler (default: g++)
 #. **$FFTW_HOME:** path to FFTW library home
 #. **$FFTW_LIB_FLAG:** FFTW library flag (default: -lfftw3)
 #. **$BLAS_HOME:** path to BLAS library home
 #. **$BLAS_LIB_FLAG:** BLAS library flags (default: -lblas -llapack -lpthread)
 #. **$EXTRA_FLAGS:** extra compiler flags for CC and CXX
 #. **$LD_EXTRA_FLAGS:** extra linker flags

3. **Compilation:** Run the command ``make`` in the ``bin/`` to compile the slabcc.
4. **Cleanup:** You can run ``make clean`` to remove the compiled objects. ``make distclean`` additionally removes all the compiled objects of the bundled external libraries.
