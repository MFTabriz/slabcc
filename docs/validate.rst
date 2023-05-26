==========
Validation
==========
We are trying to keep the slabcc compatible with as many compilers as possible by using only the standard features of the C++ language. But it is not possible to guarantee this due to the dependency on third-party components. 
The current version of the slabcc has been `build and validated <https://ci.codeberg.org/meisam/slabcc/branches/master>`_ on:

- Ubuntu Linux 16.04

 - with GNU C/C++ compilers (5), OpenBLAS, FFTW

- Ubuntu Linux 22.04

 - with GNU C/C++ compilers (9,11), OpenBLAS, FFTW
 - with GNU C/C++ compilers (11), MKL (2023)
 - with Intel oneAPI DPC++/C++ Compiler (2023), MKL (2023)
 - with LLVM Clang (14), OpenBLAS, FFTW

- AlmaLinux 8.7

 - with GNU C/C++ compilers (8), BLAS, FFTW

- openSUSE Leap 15.4

 - with GNU C/C++ compilers (10), BLAS, FFTW

Test set
--------

You can download a complete test set including input files, input parameters, and expected output `here <https://doi.org/10.5281/zenodo.1323558>`__!
You can also run the regression tests and verify their results with ``make test``. You'll need `numdiff <https://www.nongnu.org/numdiff/>`__ for these tests.
