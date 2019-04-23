:Last updated: 22 April 2019
:version: 0.6.4

.. sectnum::

.. contents::

=============
Introduction
=============
SLABCC calculates an *a posteriori* energy correction for charged slab models under 3D periodic boundary condition (PBC) based on the method proposed in:

 Hannu-Pekka Komsa and Alfredo Pasquarello, Finite-Size Supercell Correction for Charged Defects at Surfaces and Interfaces, Physical Review Letters 110, 095505 (2013) DOI: `10.1103/PhysRevLett.110.095505 <https://doi.org/10.1103/PhysRevLett.110.095505>`_ `(Supplements) <https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.110.095505/supplR1.pdf>`_
 
This method estimates the error in the total energy of the charged models under 3D PBC, due to the excess charge in the real system using a simple Gaussian model.
The model charge is assumed to be embedded in a medium with dielectric-tensor profile depending only on one Cartesian space axis orthogonal to the slab.
The energy correction is calculated as:

    E\ :sub:`corr` \  = E\ :sub:`isolated` \ - E\ :sub:`periodic` \ - qΔV

where, E\ :sub:`corr` \ is the total energy correction for the model, 
E\ :sub:`periodic` \ is the energy of the model charge, calculated by solving the periodic Poisson equation. E\ :sub:`isolated` \ is the energy of the model charge embedded in the dielectric medium and can be determined by extrapolation.
q is the total extra charge and ΔV is the difference between the potential of the Gaussian model charge system and the DFT calculations.

The code can also calculate the charge correction for the 2D models under PBC. The isolated energies for the 2D models are calculated by extrapolation based on the method proposed in:

 Hannu-Pekka Komsa, Natalia Berseneva, Arkady V. Krasheninnikov, and Risto M. Nieminen, Charged Point Defects in the Flatland: Accurate Formation Energy Calculations in Two-Dimensional Materials, Physical Review X 4, 031044 (2014) DOI: `10.1103/PhysRevX.4.031044 <https://doi.org/10.1103/PhysRevX.4.031044>`_ `(Erratum) <https://doi.org/10.1103/PhysRevX.8.039902>`_ 
 
And by the cylindrical Bessel expansion of the Poisson equation as proposed in:

 Ravishankar Sundararaman, and Yuan Ping, First-principles electrostatic potentials for reliable alignment at interfaces and defects, The Journal of Chemical Physics 146, 104109 (2017) DOI: `10.1063/1.4978238 <https://doi.org/10.1063/1.4978238>`_

| SLABCC have been initially developed for the `Bremen Center for Computational Materials Science (BCCMS) <http://www.bccms.uni-bremen.de>`_

Please read the `manual.html <http://htmlpreview.github.io/?https://github.com/MFTabriz/slabcc/blob/master/manual.html>`_ for more information on the algorithm, input and outputs.

=================
Quick start guide
=================
To calculate the charge correction slabcc needs the following files:

- Input parameters file (default: `slabcc.in`)
- CHGCAR of the neutral system from the VASP calculation (default: `CHGCAR.N`)
- CHGCAR of the charged system from the VASP calculation (default: `CHGCAR.C`)
- LOCPOT of the neutral system from the VASP calculation (default: `LOCPOT.N`)
- LOCPOT of the charged system from the VASP calculation (default: `LOCPOT.C`)

Input parameters file for a slab should minimally include (all in relative scale [0 1]):

- ``charge_position``: position of the localized charge
- ``diel_in``: dielectric tensor of the slab
- ``normal_direction``: direction normal to the surface
- ``interfaces``: position of the surfaces of the slab (in the normal direction)


Example
--------
The following examples list the input parameters to be defined in `slabcc.in` file, assuming the VASP outputs (LOCPOT and CHGCAR files) to be in the same directory. Please read the `manual.html <http://htmlpreview.github.io/?https://github.com/MFTabriz/slabcc/blob/master/manual.html>`_ for complete list of the input parameters.

1. **Minimum input**: The program models the extra charge with a Gaussian charge distribution localized around the position (``charge_position= 0.24  0.56  0.65``) in a slab model with normal direction of (``normal_direction = y``) and surfaces at (``interfaces = 0.25  0.75``). The dielectric tensor inside of the slab is assumed to be isotropic (``diel_in = 4.8``)::

    charge_position = 0.24  0.56  0.65
    diel_in = 4.8
    normal_direction = y
    interfaces = 0.25 0.75

 By default the program will also:

 - Load the CHGCAR of charged and neutralized systems. 
 - Load the LOCPOT of charged and neutralized systems.  
 - Calculate the total extra charge from the difference between the charged and neutralized CHGCARs.
 - Optimize the ``charge_position``, ``interfaces`` and ``charge_sigma``.
 - Calculate the total energy correction for the charged system.
 - Write all the input parameters used for calculation, optimized parameters and the results to output file.

2. **Correction with multiple localized Gaussian charges:** If a single charge cannot represent your localized charge properly, you can use multiple Gaussian charges in your model. You have to define the positions of each Gaussian charge as shown in example below::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    charge_position = 0.24  0.56  0.65; 0.20  0.50  0.65
    diel_in = 4.8
    normal_direction = a
    interfaces = 0.25 0.75

3. **Correction for the uniform dielectric medium e.g. bulk models:** You must have the same dielectric tensor inside and outside::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    charge_position = 0.24  0.56  0.65
    diel_in = 4.8
    diel_out = 4.8

4. **Correction for the monolayers i.e. 2D models (without extrapolation):** To use the Bessel expansion of the Poisson equation for calculating the isolated energy of the 2D models, in-plane dielectric constants must be equal and the model must be surrounded by the vacuum. Use the extrapolation method (``extrapolate=yes``) for more general cases::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    2D_model = yes
    charge_position = 0.5 0.4 0.56
    interfaces =  0.66 0.46
    normal_direction = z
    diel_in = 6.28 6.28 1.83
    diel_out = 1

5. **Correction for the monolayers i.e. 2D models (with extrapolation):** To calculate the isolated energy by fitting the extrapolation results with the non-linear formula, extrapolation to relatively large cell sizes (α < 0.2) is necessary. A large extrapolation grid size is also required to calculate the energies accurately::

    LOCPOT_charged = CHARGED_LOCPOT
    LOCPOT_neutral = UNCHARGED_LOCPOT
    CHGCAR_charged = CHARGED_CHGCAR
    CHGCAR_neutral = UNCHARGED_CHGCAR
    2D_model = yes
    extrapolate = yes
    extrapolate_steps_number = 20
    extrapolate_grid_x = 4.5
    charge_position = 0.5 0.4 0.56
    interfaces =  0.66 0.46
    normal_direction = z
    diel_in = 6.28 6.28 1.83

Test set
--------

You can download a complete test set including input files, input parameters and expected output `here <https://doi.org/10.5281/zenodo.1323558>`_!

============
Installation 
============
1. **Prerequisites:**

 #. **Compiler:** You need a C++ compiler with C++14 standard support (e.g. `g++ <https://gcc.gnu.org/>`_ 5.0 or later, `icpc <https://software.intel.com/en-us/c-compilers>`_ 15.0 or later, etc.) 
 #. **FFTW:** If you don't have FFTW installed on your system follow the guide on the `FFTW website <http://www.fftw.org/download.html>`_
 #. **BLAS/OpenBLAS/MKL:** You can use BLAS for the matrix operations inside the slabcc but it is highly recommended to use the `OpenBLAS <https://github.com/xianyi/OpenBLAS/releases>`_/`MKL <https://software.intel.com/en-us/mkl>`_ instead. If you don't have OpenBLAS installed on your system, follow the guide on the `OpenBLAS website <http://www.openblas.net>`_

2. **Configuration:** You must edit the `src/makefile` to choose your compiler and add the paths to FFTW and OpenBLAS libraries. 
3. **Compilation:** Run the command `make` in the `src/` to compile the slabcc.

**Note**: By default, the code will be compiled for the specific architecture of your compilation machine. If you are compiling and running the slabcc on different machines, you must edit the makefile and change the ``-march`` flag.

==================================
Known issues and limitations
==================================
- Shape of the VASP files cell is limited to orthogonal cells.
- Maximum line length of the input file (slabcc.in) is 4000 bytes.
- Bessel expansion of the Poisson equation cannot be used for the calculation of isolated energies for the 2D models with anisotropic in-plane screening, trivariate Gaussian model change, or the models which are not surrounded by the vacuum (diel_out > 1). Extrapolation method must be used in these cases.

===============
Release history
===============
* 2019-04-04: version 0.6 - Added trivariate Gaussian model charge and selective charge optimization support
* 2019-03-18: version 0.5 - Added 2D model support
* 2018-10-10: version 0.4 - Added spdlog. General interface and performance improvements
* 2018-07-29: version 0.3 - First public release

===========================
Copyright and attributions
===========================
Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz

The source codes and all the documentations are available under The 2-Clause BSD License. For more information see license_.

| This code uses several open source components each of which are located under a separate sub-directory in the `src/`. The copyright of these libraries belong to their respective owners. Any modification made to those codes is also published under the same license. We acknowledge and are grateful to these developers and maintainers for their valuable contributions to this software and more importantly to the free software society.
| The attributions are also present in the binary file and can be accessed using the command-line parameters.

Included libraries
------------------

- `Armadillo C++ Linear Algebra Library <http://arma.sourceforge.net>`_ licensed under the Apache License 2.0
 
 - Copyright 2008-2018, Conrad Sanderson
 - Copyright 2008-2016, National ICT Australia (NICTA)
 - Copyright 2017-2018, Arroyo Consortium
 - Copyright 2017-2018, Data61, CSIRO
 - This product includes software developed by Conrad Sanderson
 - This product includes software developed at National ICT Australia (NICTA)
 - This product includes software developed at Arroyo Consortium
 - This product includes software developed at Data61, CSIRO

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
 
Linked libraries
---------------------

- `FFTW <http://www.fftw.org>`_ licensed under GNU General Public License
- `OpenBLAS <http://www.openblas.net>`_ licensed under the 3-clause BSD license 

License
-------
Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


.. image:: https://api.codacy.com/project/badge/Grade/85d640604a454870b04081617ed538c9
   :alt: Codacy Badge
   :target: https://app.codacy.com/app/MFTabriz/slabcc?utm_source=github.com&utm_medium=referral&utm_content=MFTabriz/slabcc&utm_campaign=Badge_Grade_Dashboard