==========
Validation
==========
We are trying to keep the slabcc compatible with as many compilers as possible by using only the standard features of the C++ language. However, it is not possible to guarantee this due to the dependency on third-party components.
The current version of the slabcc has been `built and validated <https://ci.codeberg.org/meisam/slabcc/branches/master>`_ inside containers based on:

- Ubuntu Linux 16.04 (`dockerfile <https://codeberg.org/meisam/slabcc/src/branch/master/ci/dockerfile.ubuntu.16.04>`__)

 - with GNU C/C++ compilers (5), OpenBLAS, FFTW

- Ubuntu Linux 22.04 (dockerfiles `1 <https://codeberg.org/meisam/slabcc/src/branch/master/ci/dockerfile.ubuntu.22.04>`__ `2 <https://raw.githubusercontent.com/intel/oneapi-containers/master/images/docker/basekit/Dockerfile.ubuntu-22.04>`__)

 - with GNU C/C++ compilers (9,11), OpenBLAS, FFTW
 - with GNU C/C++ compilers (11), MKL (2023)
 - with Intel oneAPI DPC++/C++ Compiler (2023), MKL (2023)
 - with LLVM Clang (14), OpenBLAS, FFTW

- AlmaLinux 8.7 (`dockerfile <https://codeberg.org/meisam/slabcc/src/branch/master/ci/dockerfile.almalinux.8.7>`__)

 - with GNU C/C++ compilers (8), BLAS, FFTW

- openSUSE Leap 15.4 (`dockerfile <https://codeberg.org/meisam/slabcc/src/branch/master/ci/dockerfile.opensuse.leap.15.4>`__)

 - with GNU C/C++ compilers (10), BLAS, FFTW

Older versions of the code were also being tested on MS Windows 10 (latest toolchains: Intel C/C++ compilers 19.0.4 and Microsoft C/C++ compilers 19.20.27508 linked to MKL 19.0.4 and FFTW 3.3.5), support for which is currently dropped.

Test set
--------

You can download a complete test set, including input files, input parameters, and expected output, `here <https://doi.org/10.5281/zenodo.1323558>`__!
You can also run regression tests and verify their results with ``make test``. You'll need `numdiff <https://www.nongnu.org/numdiff/>`__ for these tests.
