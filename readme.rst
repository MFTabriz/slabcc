.. |Version| image:: https://img.shields.io/badge/version-1.0.0-blue.svg
   :alt: Version 1.0.0
   :target: https://codeberg.org/meisam/slabcc/releases
.. |Manual| image:: https://img.shields.io/badge/Manual-HTML-informational
   :alt: HTML Manual
   :target: https://meisam.codeberg.page/slabcc
.. |Standard| image:: https://img.shields.io/badge/c%2B%2B-%3E%3D14-informational
   :alt: C++ Standard
   :target: https://en.cppreference.com/w/cpp/compiler_support#cpp14
.. |Style| image:: https://img.shields.io/badge/code%20style-LLVM-black
   :alt: Code style
   :target: https://llvm.org/docs/CodingStandards.html
.. |Woodpecker| image:: https://ci.codeberg.org/api/badges/meisam/slabcc/status.svg
   :alt: Build Status
   :target: https://ci.codeberg.org/meisam/slabcc/branches/master
.. |Zenodo| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.1323558.svg
   :alt: Test Set
   :target: https://doi.org/10.5281/zenodo.1323558
.. |License| image:: https://img.shields.io/badge/License-BSD--2--Clause-blue
   :alt: License: BSD-2-Clause
   :target: https://opensource.org/licenses/BSD-2-Clause

- For users: |Version| |Manual| |License|
- For developers and maintainers: |Standard| |Style| |Woodpecker|
- Large test set: |Zenodo|

.. sectnum::

.. contents::

=============
Introduction
=============
SLABCC calculates *a posteriori* energy correction for charged slab models under 3D periodic boundary conditions (PBC) based on the method proposed in:

 Hannu-Pekka Komsa and Alfredo Pasquarello, Finite-Size Supercell Correction for Charged Defects at Surfaces and Interfaces, Physical Review Letters 110, 095505 (2013) DOI: `10.1103/PhysRevLett.110.095505 <https://doi.org/10.1103/PhysRevLett.110.095505>`_ `(Supplements) <https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.110.095505/supplR1.pdf>`_

This method estimates the error in the total energy of the charged models under 3D PBC due to the excess charge in the real system using Gaussian models.
The model charge is assumed to be embedded in a medium with a dielectric-tensor profile ε(k) depending only on a single Cartesian space axis (k) that is orthogonal to the slab.
The energy correction is calculated as:

    E\ :sub:`corr` \  = E\ :sub:`isolated` \ - E\ :sub:`periodic` \ - qΔV

where E\ :sub:`corr` \ is the total energy correction for the model;
E\ :sub:`periodic` \ is the energy of the model charge, calculated by solving the periodic Poisson equation. E\ :sub:`isolated` \ is the energy of the model charge embedded in the dielectric medium and can be determined by extrapolation.
q is the total extra charge, and ΔV is the difference between the potential of the Gaussian model charge system and the DFT calculations.

The code can also calculate the charge correction for the 2D models under PBC. The isolated energies for the 2D models are calculated by extrapolation based on the method proposed in:

 Hannu-Pekka Komsa, Natalia Berseneva, Arkady V. Krasheninnikov, and Risto M. Nieminen, Charged Point Defects in the Flatland: Accurate Formation Energy Calculations in Two-Dimensional Materials, Physical Review X 4, 031044 (2014) DOI: `10.1103/PhysRevX.4.031044 <https://doi.org/10.1103/PhysRevX.4.031044>`_ `(Erratum) <https://doi.org/10.1103/PhysRevX.8.039902>`_

And by the cylindrical Bessel expansion of the Poisson equation as proposed in:

 Ravishankar Sundararaman, and Yuan Ping, First-principles electrostatic potentials for reliable alignment at interfaces and defects, The Journal of Chemical Physics 146, 104109 (2017) DOI: `10.1063/1.4978238 <https://doi.org/10.1063/1.4978238>`_

| SLABCC have been initially developed for the `Bremen Center for Computational Materials Science (BCCMS) <https://www.uni-bremen.de/bccms>`_. The manual is accessible `online <https://meisam.codeberg.page/slabcc/>`__.

=================
Quick start guide
=================
To calculate the charge correction slabcc needs the following files:

- Input parameters file (default: ``slabcc.in``)
- CHGCAR of the neutral system from the VASP calculation (default: ``CHGCAR.N``)
- CHGCAR of the charged system from the VASP calculation (default: ``CHGCAR.C``)
- LOCPOT of the neutral system from the VASP calculation (default: ``LOCPOT.N``)
- LOCPOT of the charged system from the VASP calculation (default: ``LOCPOT.C``)

Input parameters file for a slab should minimally include (all in relative scale [0 1]):

- ``charge_position``: position of the localized charge
- ``diel_in``: dielectric tensor of the slab
- ``normal_direction``: direction normal to the surface
- ``interfaces``: position of the surfaces of the slab (in the normal direction)


Example
--------
The following examples list the `input parameters`_ to be defined in ``slabcc.in`` file, assuming the VASP outputs (LOCPOT and CHGCAR files) are in the same directory.

1. **Minimum input**: The program models the extra charge with a Gaussian charge distribution localized around the position (``charge_position= 0.24  0.56  0.65``) in a slab model with a normal direction of (``normal_direction = y``) and surfaces at (``interfaces = 0.25  0.75``). The dielectric tensor inside the slab is assumed to be isotropic (``diel_in = 4.8``)::

    charge_position = 0.24  0.56  0.65
    diel_in = 4.8
    normal_direction = y
    interfaces = 0.25 0.75

 The program will use the default values for the other parameters. Afterwards, slabcc will:

- Calculate the total extra charge from the difference between the charged and neutralized CHGCARs.
- Optimize the ``charge_position``, ``interfaces`` and ``charge_sigma``.
- Calculate the total energy correction for the charged system.
- Write all the input parameters used for calculation, the optimized parameters, and the results to the output file.

2. **Correction with multiple localized Gaussian charges:** If a single charge cannot represent your localized charge properly, you can use multiple Gaussian charges in your model. You have to define the positions of each Gaussian charge, as shown in the example below::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    charge_position = 0.24  0.56  0.65; 0.20  0.50  0.65
    diel_in = 4.8
    normal_direction = a
    interfaces = 0.25 0.75

3. **Correction for the uniform dielectric medium, e.g., bulk models:** You must have the same dielectric tensor inside and outside::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    charge_position = 0.24  0.56  0.65
    diel_in = 4.8
    diel_out = 4.8

4. **Correction for the monolayers, i.e., 2D models (without extrapolation):** To use the Bessel expansion of the Poisson equation for calculating the isolated energy of the 2D models, the in-plane dielectric constants must be equal and the model must be surrounded by a vacuum. Use the extrapolation method (``extrapolate=yes``) for more general cases::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    2D_model = yes
    charge_position = 0.5 0.4 0.56
    interfaces = 0.66 0.46
    normal_direction = z
    diel_in = 6.28 6.28 1.83
    diel_out = 1

5. **Correction for the monolayers, i.e., 2D models (with extrapolation):** To calculate the isolated energy by fitting the extrapolation results with the non-linear formula, extrapolation to relatively large cell sizes (1/α < 0.2) is necessary. To avoid large discretization errors, the size of the extrapolation grid will be automatically increased::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    2D_model = yes
    extrapolate = yes
    charge_position = 0.5 0.4 0.56
    interfaces = 0.66 0.46
    normal_direction = z
    diel_in = 6.28 6.28 1.83

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

==================================
Known issues and limitations
==================================
- Only orthogonal cells are supported.
- The maximum line length of the input file (slabcc.in) is 4000 bytes.
- Bessel expansion of the Poisson equation cannot be used for the calculation of isolated energies for the 2D models with anisotropic in-plane screening, trivariate Gaussian model change, or the models that are not surrounded by the vacuum (diel_out > 1). The extrapolation method must be used in these cases.

==========================
Release history highlights
==========================
* 2023-06-05: version 1.0 - Improved build, error mitigation, and correctness checks
* 2019-06-13: version 0.8 - OO redesign
* 2019-05-14: version 0.7 - Added discretization error mitigation
* 2019-04-04: version 0.6 - Added trivariate Gaussian model charge and selective charge optimization support
* 2019-03-18: version 0.5 - Added 2D model support
* 2018-10-10: version 0.4 - Added spdlog and several user interface and performance improvements
* 2018-07-29: version 0.3 - First public release

===========================
Copyright and attributions
===========================
Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz

The source code and all documentation are available under the 2-Clause BSD License. For more information, see license_.

| This code uses several open-source components, each of which is located under a separate sub-directory of ``src/``. The copyrights of these libraries belong to their respective owners. Any modification made to those codes is also published under the same license. We acknowledge and are grateful to these developers and maintainers for their valuable contributions to this software and, more importantly, to the free software society.
| Attributions are also present in the binary file and can be accessed by using ``--copyright`` flag.

Included third-party components
-------------------------------

- `Armadillo C++ Linear Algebra Library <https://arma.sourceforge.net>`_ licensed under the Apache License 2.0

 - Copyright 2008-2023 Conrad Sanderson
 - Copyright 2008-2016 National ICT Australia (NICTA)
 - Copyright 2017-2023 Data61 / CSIRO
 - This product includes software developed by Conrad Sanderson
 - This product includes software developed at National ICT Australia (NICTA)
 - This product includes software developed at Data61 / CSIRO

- `inih <https://github.com/benhoyt/inih>`_ (INI Not Invented Here) licensed under the 3-clause BSD license

 - © 2009, Ben Hoyt, `et al. <https://github.com/benhoyt/inih/contributors>`__

- `clara <https://github.com/catchorg/Clara>`_ licensed under the Boost Software License 1.0

 - © 2014, Phil Nash, Martin Hořeňovský, `et al. <https://github.com/catchorg/Clara/contributors>`__

- `spline <https://shiftedbits.org/2011/01/30/cubic-spline-interpolation/>`_ (Cubic Spline Interpolation) licensed under the Beer-Ware License 42

 - © 2011, Devin Lane

- `NLOPT <https://nlopt.readthedocs.io>`_ licensed under the GNU LGPL

 - © 2007-2014, Massachusetts Institute of Technology
 - © 2007-2014, Steven G. Johnson `et al. <https://github.com/stevengj/nlopt/contributors>`__

- `spdlog <https://github.com/gabime/spdlog>`_ licensed under the MIT License

 - © 2016, Gabi Melman, `et al. <https://github.com/gabime/spdlog/contributors>`__

- `Boost.Predef <https://github.com/boostorg/predef>`_ licensed under the Boost Software License 1.0

 - © 2005-2018 Rene Rivera
 - © 2015 Charly Chevalier
 - © 2015 Joel Falcou, `et al. <https://github.com/boostorg/predef/contributors>`__

License
-------
Copyright (c) 2018-2023, University of Bremen, M. Farzalipour Tabriz

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

