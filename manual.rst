**Note**: github `does not support <https://github.com/github/markup/issues/274#issuecomment-77102262>`_ math equations the in reStructuredText format. Please check the `manual.html <http://htmlpreview.github.io/?https://github.com/MFTabriz/slabcc/blob/master/manual.html>`_ for the proper rendering!

:Last updated: 13 June 2019
:version: 0.8.0

.. sectnum::

.. contents::

=============
Introduction
=============
SLABCC calculates an *a posteriori* energy correction for charged slab models under 3D periodic boundary condition (PBC) based on the method proposed in:

 Hannu-Pekka Komsa and Alfredo Pasquarello, Finite-Size Supercell Correction for Charged Defects at Surfaces and Interfaces, Physical Review Letters 110, 095505 (2013) DOI: `10.1103/PhysRevLett.110.095505 <https://doi.org/10.1103/PhysRevLett.110.095505>`_ `(Supplements) <https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.110.095505/supplR1.pdf>`_
 
This method estimates the error in the total energy of the charged models under 3D PBC, due to the excess charge in the real system using a simple Gaussian model.
The model charge is assumed to be embedded in a medium with dielectric-tensor profile ε(k) depending only on a single Cartesian space axis (k) which is orthogonal to the slab.
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

Algorithm
----------
* slabcc reads the `VASP <https://www.vasp.at>`_ output files (CHGCAR and LOCPOT), and calculates the extra charge distribution and the potential due to the extra charge. All the input files should correspond to the same geometry.
* The extra charge is modeled by sum of Gaussian charge distributions as:

.. math::

 \rho(r) = \sum_{i}\frac{q_i}{\sigma_{i}^{3}(2\pi)^{3/2}} \exp \left ({- \frac{r_{i}^{2}}{2\sigma_{i}^{2}} } \right ) 
  
or `trivariate Gaussian <http://mathworld.wolfram.com/TrivariateNormalDistribution.html>`_ charge distributions as:

.. math::
 
 \rho(r) = \sum_{i}\frac{q_i}{\sigma_{i,x}\sigma_{i,y}\sigma_{i,z}(2\pi)^{3/2}} \exp \left ({- \frac{r_{i,x}^{2}}{2\sigma_{i,x}^{2}} - \frac{r_{i,y}^{2}}{2\sigma_{i,y}^{2}}- \frac{r_{i,z}^{2}}{2\sigma_{i,z}^{2}} } \right )

which gives a charge distribution normalized to q\ :sub:`i` \ with standard deviation of σ\ :sub:`i` \ (``charge_sigma``) for each Gaussian distribution i centered at r\ :sub:`i` \ (``charge_position``).

* The generated Gaussian model charge will be embedded in a dielectric medium with profile of the form:

.. math::
  \epsilon (k) =  \frac{\epsilon_2-\epsilon_1}{2} \text{erf}\left(\frac{k-k_0 }{\beta}\right)+\frac{\epsilon_2+\epsilon_1}{2}

where k\ :sub:`0` \ is the interface position in Cartesian k-direction, ε\ :sub:`1` \ and ε\ :sub:`2` \ are dielectric tensors on either side of interface (``diel_in`` & ``diel_out``) and β (``diel_taper``) defines the smoothness of transition assuming anisotropic dielectric tensor as:

.. math::
 \epsilon = 
 \left| \begin{bmatrix}
    \epsilon_{11} & 0 & 0 \\
    0 & \epsilon_{22} & 0 \\
    0 & 0&  \epsilon_{33}
 \end{bmatrix}\right|

* The potential due to the charge distribution ρ under 3D PBC embedded in the dielectric medium ε(k) is calculated by solving the Poisson equation in Fourier space:

.. math::
	 \epsilon(k) \nabla^2 V(r)+\frac{\partial}{\partial k} \epsilon(k)\frac{\partial}{\partial k}V(r) = -\rho(r)

* A non-linear optimization routine minimizes the difference of our calculated V(r) for the model charge and the V resulted from the VASP calculation by changing the position of the model Gaussian charge, its width, and the position of the slab interfaces.

* The E\ :sub:`periodic` is calculated as:

.. math::
	E = \frac{1}{2} \int V(r) \rho(r) \, dr

* E\ :sub:`isolated` is calculated the same way as E\ :sub:`periodic` but with extrapolation of the fixed model charge embedded in an infinitely large dielectric medium. For the bulk and the slab models, the extrapolation is done linearly. For the monolayer models (2D systems) the following equation is used for the extrapolation [`10.1103/PhysRevX.8.039902 <https://doi.org/10.1103/PhysRevX.8.039902>`_]:

.. math::
	E = c_0 + c_1 x + c_2 x^2 + d e^{-c_3 x}

 where c\ :sub:`i` are the fitting parameters and 

.. math::
	d =  \frac{c_1 - \frac{\partial E_M}{\partial x}}{c_3}

 guarantees the correct energy gradient at x(=1/α)→0. E\ :sub:`M` being the Madelung energy.

* ΔV is calculated at the position least affected by the model charge.

More information about the algorithms and the implementation details can be found `here`__.

__ cite_
	 
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
The following examples list the `input parameters`_ to be defined in `slabcc.in` file, assuming the VASP outputs (LOCPOT and CHGCAR files) to be in the same directory.

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

You can download a complete test set including input files, input parameters and expected output `here <https://doi.org/10.5281/zenodo.1323558>`_! Bitwise reproducibility of the results is not guaranteed across the different versions.

============
Installation
============
1. **Prerequisites:**

 #. **Compiler:** You need a C++ compiler with C++14 standard support (e.g. `g++ <https://gcc.gnu.org/>`_ 5.0 or later, `icpc <https://software.intel.com/en-us/c-compilers>`_ 15.0 or later, etc.) 
 #. **BLAS/OpenBLAS/MKL:** You can use BLAS+LAPACK for the matrix operations inside the slabcc but it is highly recommended to use one of the high performance replacements e.g. the `OpenBLAS <https://github.com/xianyi/OpenBLAS/releases>`_/`MKL <https://software.intel.com/en-us/mkl>`_ instead. If you don't have OpenBLAS installed on your system, follow the guide on the `OpenBLAS website <http://www.openblas.net>`_. Please refer to the `Armadillo documentations <https://gitlab.com/conradsnicta/armadillo-code/blob/9.100.x/README.md>`_ for linking to the other BLAS replacements.
 #. **FFTW:** If you don't have FFTW installed on your system follow the guide on the `FFTW website <http://www.fftw.org/download.html>`_. Alternatively, you can use the FFTW interface of the MKL.

2. **Configuration:** You must edit the `src/makefile` to choose your compiler and add the paths to FFTW and BLAS libraries. 
3. **Compilation:** Run the command `make` in the `src/` to compile the slabcc.
4. **Cleanup:** You can run `make clean` to remove the compiled objects, and static library files. `make distclean` additionally removes all the compiled objects in the external libraries.

**Note**: By default, the code will be compiled for the specific architecture of your compilation machine. If you are compiling and running the slabcc on different machines, you must edit the makefile and change the ``-march`` flag.

==========
Validation
==========
We are trying to keep the slabcc compatible with as many compilers as possible by using only the standard features of the C++ language. But it is not possible to guarantee this due to the dependency on the third-party components. 
The current version of the slabcc has been validated on:

- CentOS Linux release 7.6.1810
 - with Intel C++ compiler 18.0.3, MKL 18.0.3, FFTW (from MKL)
- Ubuntu Linux release 16.04.6 (`Travis <https://travis-ci.org/MFTabriz/slabcc>`_)
 - with GNU C++ compilers (5.5.0/6.5.0/8.1.0/9.1.0), OpenBLAS 0.2.18, FFTW 3.3.4
- Microsoft Windows version 10.0.17134
 - with Intel C++ compiler 19.0.4, MKL 19.0.4, FFTW 3.3.5
 
=======================
Command-line parameters
=======================
You can run slabcc without any additional options. Alternatively, you can use the following options to modify its behavior:

-h, --help						Display the usage information (this list)
-i, --input <input_file>			slabcc input file name
-o, --output <input_file>			slabcc output file name
-l, --log <log_file>			slabcc log file name
-d, --diff						Calculate the charge and the potential differences only
-m, --manual					Show the quick start guide
-v, --version					Show the slabcc version and its compilation date
-c, --copyright					Show the copyright information and the attributions

======================
Input parameters
======================
slabcc reads all its parameters from the input file (by default: `slabcc.in`) You can change the input file's name using the `command-line parameters`_.
The input file is processed as follows:

- Lines starting with # will be treated as comments. Inline comments are also allowed
- Double quotation marks will be removed from the strings
- A warning will be issued for any unidentified parameter
- A warning will be issued for the use of the deprecated parameters
- All the coordinates must be in fractional form [0-1]
- Boolean (True/False) parameters can be also declared as 0/1, on/off, yes/no, .true./.false.
- Parameter names can be written in the small or the CAPITAL letters
- For vectors and matrices, columns are separated by a “ ”(space), while the rows are separated by a “;” (semicolon)
- Lines starting with a space “ ” will be treated as the continuation of the last parameter's value
- Subsequent definitions for any parameter will be concatenated to the existing definition

 
+------------------------------+-------------------------------------------------------+---------------+
| Parameter                    | Description and the options / ``example``             | Default value |
+==============================+=======================================================+===============+
| ``2d_model``                 | Calculate the charge correction for a 2D model        |  false        |
|                              |                                                       |               |
|                              |                                                       |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Fraction of the total extra charge in each localized   |*The extra     |
|                              |Gaussian model charge (in the case of multiple Gaussian|charge will be |
| ``charge_fraction``          |charges)                                               |equally divided|
|                              |                                                       |among all      |
|                              |``charge_fraction = 0.4 0.6``                          |positions*     |
|                              |                                                       |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Center of the model Gaussian charges                   |               |
| ``charge_position``          |                                                       |               |
|                              |``charge_position = 0.2 0.5 0.3``                      |               |
|                              |                                                       |               |
|                              |``charge_position = 0.2 0.2 0.2; 0.3 0.4 0.3``         |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Rotation angles around each axis for the trivariate    |               |
|                              |Gaussian charges in arc degree (-90, 90)               |       0       |
| ``charge_rotation``          |                                                       |               |
|                              |``charge_rotation = 0 45 0``                           |               |
|                              |                                                       |               |
|                              |``charge_rotation = 45 0 0; 0 -10 0``                  |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |width of each localized Gaussian charge. It can be 1   |               |
|                              |or in case of trivariate models, 3 parameters per      |               |
|                              |localized Gaussian charge. For trivariate Gaussian     |               |
|                              |models, defining a single parameter per charge, sets   |               |
|                              |the sigma values to be equal in all directions.        |               |
|                              |                                                       |               |
|                              |for a single Gaussian charge                           |               |
| ``charge_sigma``             |``charge_sigma = 1``                                   |1 (for each    |
|                              |                                                       |charge in each |
|                              |for multiple Gaussian charges                          |direction)     |
|                              |``charge_sigma = 1 1.5``                               |               |
|                              |                                                       |               |
|                              |for two trivariate Gaussian charges                    |               |
|                              |``charge_sigma = 1 2 3; 1.5 2.5 3.5;``                 |               |
|                              |                                                       |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``charge_trivariate``        |Use trivariate Gaussian model along the main axis      |   false       |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Charge density file (CHGCAR) of the charged system     |               |
| ``CHGCAR_charged``           |                                                       | CHGCAR.C      |
|                              |``CHGCAR_charged = CHGCAR1``                           |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Charge density file (CHGCAR) of the neutral system     |               |
| ``CHGCAR_neutral``           |                                                       | CHGCAR.N      |
|                              |``CHGCAR_neutral = CHGCAR2``                           |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``diel_in``                  |Diagonal elements of the static dielectric tensor      |       1       |
|                              |inside of the slab. If only a single value is given,   |               |
|                              |all of them will be assumed to be equal.               |               |
|                              |                                                       |               |
|                              |``diel_in = 3``                                        |               |
|                              |                                                       |               |
|                              |``diel_in = 3 4 5``                                    |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``diel_out``                 |Diagonal elements of the static dielectric tensor      |       1       |
|                              |outside of the slab                                    |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``diel_taper``               |The steepness of the transition between diel_in and    |       1       |
|                              |diel_out (β in the dielectric profile formula)         |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``extrapolate``              |Calculate the isolated energy using the extrapolation  |opposite of the|
|                              |method                                                 |``2d_model``   |
|                              |                                                       |parameter      |
|                              |                                                       |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Extrapolation grid size multiplier. The number of the  |               |
|                              |grid points in each direction will be multiplied by    |               |
|                              |this value.                                            |               |
|                              |                                                       |               |
|                              |extrapolate_grid_x > 1 will use a larger grid in the   |               |
|``extrapolate_grid_x``        |extrapolations which will increase its accuracy but    |       1       |
|                              |will requires more memory and the computational power. |               |
|                              |                                                       |               |
|                              |extrapolate_grid_x = 1 will use the same grid size as  |               |
|                              |the VASP input files in the extrapolation.             |               |
|                              |                                                       |               |
|                              |extrapolate_grid_x < 1 will use a smaller grid for the |               |
|                              |extrapolations which increases the speed and decreases |               |
|                              |the memory usage but the energies for the higher orders|               |
|                              |of the extrapolation may not be accurate!              |               |
|                              |                                                       |               |
|                              |``extrapolate_grid_x = 1.8``                           |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Number of the extrapolation steps in calculation of    |               |
| ``extrapolate_steps_number`` |E\ :sub:`isolated` \ [#]_                              |       4       |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Size of extrapolation steps with respect to the initial|               |
| ``extrapolate_steps_size``   |supercell size                                         |       0.5     |
+------------------------------+-------------------------------------------------------+---------------+
| ``interfaces``               |Interfaces of the slab in normal direction             |   0.25 0.75   |
|                              |                                                       |               |
|                              |``interfaces = 0.11 0.40``                             |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Local potential file (LOCPOT) of the charged system    |               |
| ``LOCPOT_charged``           |                                                       |   LOCPOT.C    |
|                              |``LOCPOT_charged = LOCPOT1``                           |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Local potential file (LOCPOT) of the neutral system    |               |
| ``LOCPOT_neutral``           |                                                       |   LOCPOT.N    |
|                              |``LOCPOT_neutral = LOCPOT2``                           |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``normal_direction``         |Normal direction of the slab: one of x/y/z or a/b/c    |      z        |
|                              |corresponding to the 1st, 2nd and 3rd vectors in the   |               |
|                              |input file's cell vectors                              |               |
|                              |                                                       |               |
|                              |``normal_direction = b``                               |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_algorithm``       |Optimization algorithm in the NLOPT library            |    BOBYQA     |
|                              |                                                       |               |
|                              |`BOBYQA <https://en.wikipedia.org/wiki/BOBYQA>`_ :     |               |
|                              |Bound Optimization BY Quadratic Approximation [#]_     |               |
|                              |                                                       |               |
|                              |`COBYLA <https://en.wikipedia.org/wiki/COBYLA>`_:      |               |
|                              |Constrained Optimization BY Linear Approximation [#]_  |               |
|                              |                                                       |               |
|                              |SBPLX: S.G. Johnson's implementation of the            |               |
|                              |Subplex (subspace-searching simplex) algorithm [#]_    |               |
|                              |                                                       |               |
|                              |``optimize_algorithm = SBPLX``                         |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_charge_fraction`` |**true**: find the optimal values for the model's      |     true      |
|                              |charge_fraction parameter to construct the best model  |               |
|                              |charge which mimics the potential obtained from the    |               |
|                              |VASP calculation                                       |               |
|                              |                                                       |               |
|                              |**false**: do not change the charge_fraction parameter |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_charge_position`` |**true**: find the optimal values for the model's      |     true      |
|                              |charge_position parameter to construct the best model  |               |
|                              |charge which mimics the potential obtained from the    |               |
|                              |VASP calculation                                       |               |
|                              |                                                       |               |
|                              |**false**: do not change the charge_position parameter |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_charge_rotation`` |**true**: find the optimal values for the model's      |     false     |
|                              |charge_rotation parameter to construct the best model  |               |
|                              |charge which mimics the potential obtained from the    |               |
|                              |VASP calculation. This can only be used for the        |               |
|                              |trivariate Gaussian models.                            |               |
|                              |                                                       |               |
|                              |**false**: do not change the charge_rotation parameter |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_charge_sigma``    |**true**: find the optimal values for the model's      |     true      |
|                              |charge_sigma parameter to construct the best model     |               |
|                              |charge which mimics the potential obtained from the    |               |
|                              |VASP calculation                                       |               |
|                              |                                                       |               |
|                              |**false**: do not change the charge_sigma parameter    |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Optimization grid size multiplier. The number of the   |               |
|                              |grid points in each direction will be multiplied by    |               |
|                              |this value.                                            |               |
|                              |                                                       |               |
|                              |optimize_grid_x > 1 will use a larger grid in the      |               |
| ``optimize_grid_x``          |optimization which will increase its accuracy but will |       0.8     |
|                              |requires more memory and the computational power.      |               |
|                              |[usually this is not necessary]                        |               |
|                              |                                                       |               |
|                              |optimize_grid_x = 1 will use the same grid as the      |               |
|                              |VASP input files in the optimization                   |               |
|                              |                                                       |               |
|                              |optimize_grid_x < 1 will use a smaller grid for the    |               |
|                              |optimization which increases the speed and decreases   |               |
|                              |the memory usage but the parameters obtained using very|               |
|                              |small grid sizes may be inaccurate!                    |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_interfaces``      |**true**: find the optimal values for the model's      |               |
|                              |interfaces to construct the best model which mimics    |     true      |
|                              |the potential obtained from the VASP calculation       |               |
|                              |                                                       |               |
|                              |**false**: do not change the position of interfaces in |               |
|                              |the model charge                                       |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_maxsteps``        |Maximum number of optimization steps                   |               |
|                              |                                                       |               |
|                              |``optimize_maxsteps = 2000``                           |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_maxtime``         |Maximum time for optimization in minutes               |               |
|                              |                                                       |               |
|                              |``optimize_maxtime = 1440``                            |               |
+------------------------------+-------------------------------------------------------+---------------+
| ``optimize_tolerance``       |Relative optimization tolerance (convergence criteria) |    0.01       |
|                              |for root mean square error of the model potential      |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Center of the slab. During the calculations, everything|               |
| ``slab_center``              |will be shifted to keep this point at the center. This |  0.5 0.5 0.5  |
|                              |point must be inside of the slab.                      |               |
|                              |                                                       |               |
|                              |``slab_center = 0.2 0.7 0.5``                          |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Verbosity of the program [#]_                          |               |
| ``verbosity``                |                                                       |       1       |
|                              |**0**: No extra info. Only write the output file.      |               |
|                              |Logging is disabled.                                   |               |
|                              |                                                       |               |
|                              |**1**: Display calculated energy correction terms.     |               |
|                              |Write the planar averaged potential and charge for the |               |
|                              |Gaussian model charge and the extra-charge of QM       |               |
|                              |calculations in the direction normal to the slab       |               |
|                              |surface.                                               |               |
|                              |                                                       |               |
|                              |**2**: Write extra-charge density, extra-charge        |               |
|                              |potential and dielectric profiles. Display debug info  |               |
|                              |including the compilation machine info and a few       |               |
|                              |important enviroment variables.                        |               |
|                              |                                                       |               |
|                              |**3**: Write the planar averaged files in all          |               |
|                              |directions.                                            |               |
|                              |                                                       |               |
|                              |**4**: Display the time passed since the start of      |               |
|                              |slabcc (in seconds) and a description of each          |               |
|                              |calculation step (trace mode)                          |               |
+------------------------------+-------------------------------------------------------+---------------+

.. [#] Extrapolating the model to very large order will accumulate errors due to energy calculations for large systems over a coarse grid size.
.. [#] M.J.D. Powell, `The BOBYQA algorithm for bound constrained optimization without derivatives <http://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf>`_, Department of Applied Mathematics and Theoretical Physics, Cambridge England, technical report NA2009/06 (2009).
.. [#] M.J.D. Powell, `Direct search algorithms for optimization calculations <https://doi.org/10.1017/S0962492900002841>`_, Acta Numerica, Vol. 7(1998) pp. 287-336
.. [#] T.H. Rowan, `Functional Stability Analysis of Numerical Algorithms <https://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.31.5708>`_, Ph.D. thesis, Department of Computer Sciences, University of Texas at Austin, 1990.
.. [#] Each verbosity level includes all the outputs from the lower verbosity options. Check `the files table`_ for complete list of the output files.

===============================
Results and the generated files
===============================
slabcc writes its calculated energy correction values to the standard output as well as the output file. All reported energy values are in eV.

Depending on the verbosity level of your choice, you may get additional reports from each part of the calculation in the standard output and/or extra output files. 


Output files
------------------
The parsed input variables or their default values and the calculation results will be written to the output file (by default: slabcc.out) You can change this file’s name using the `command-line parameters`_. A typical output file is shown below::

	# Parameters read from the file or their default values:
	2d_model = no
	charge_fraction = 1
	charge_position = 0.5 0.5 0.37; 
	charge_rotation = 0 0 0;
	charge_sigma = 1;
	charge_trivariate = no
	CHGCAR_charged = ../03-V_Cl_pos/CHGCAR
	CHGCAR_neutral = ../02-V_Cl/CHGCAR
	diel_in = 2.45
	diel_out = 1
	diel_taper = 1
	extrapolate = yes
	extrapolate_grid_x = 1
	extrapolate_steps_number = 4
	extrapolate_steps_size = 0.5
	interfaces = 0 0.375
	LOCPOT_charged = ../03-V_Cl_pos/LOCPOT
	LOCPOT_neutral = ../02-V_Cl/LOCPOT
	normal_direction = z
	optimize_algorithm = COBYLA
	optimize_charge_fraction = yes
	optimize_charge_position = yes
	optimize_charge_rotation = no
	optimize_charge_sigma = yes
	optimize_grid_x = 0.8
	optimize_interfaces = yes
	optimize_maxsteps = 0
	optimize_maxtime = 0
	optimize_tolerance = 0.01
	slab_center = 0.5 0.5 0.25
	verbosity = 5

	[Optimized_model_parameters]
	interfaces_optimized =  0.942000748357 0.455672787711
	charge_sigma_optimized = 1.4132676877
	charge_position_optimized = 0.501460639345 0.50145532106 0.385476689493;

	[Results]
	dV = -0.00291385176718
	E_periodic of the model charge = 2.0404453156
	E_isolated of the model charge = 2.59716677886
	Energy correction for the model charge (E_iso-E_per-q*dV) = 0.559635314929

Planar average files are written as a single column in plain text format and named as: "slabcc_{1}{2}{XXX}.dat" where:

- {1}: **N**: Neutral system, **C**: Charged system, **D**: Difference
- {2}: **X**/**Y**/**Z**: Corresponds to the 1st, 2nd, and the 3rd axis in the input files
- {XXX}: **CHG**: CHGCAR, **POT**: LOCPOT

.. _`the files table`:

All the possible output files and the minimum value of the verbosity parameter for activation of each are listed in the table below:

+------------------------+-------------------------------------------------------+---------------+
| file name              | Description                                           |   verbosity   |
+========================+=======================================================+===============+
|`slabcc_CXCHG.dat`      |Planar average of charged CHGCAR file in X direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_CXPOT.dat`      |Planar average of charged LOCPOT file in X direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_CYCHG.dat`      |Planar average of charged CHGCAR file in Y direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_CYPOT.dat`      |Planar average of charged LOCPOT file in Y direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_CZCHG.dat`      |Planar average of charged CHGCAR file in Z direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_CZPOT.dat`      |Planar average of charged LOCPOT file in Z direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_D.CHGCAR`       |Difference in the neutral and charged CHGCAR files     |2              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_D.LOCPOT`       |Difference in the neutral and charged LOCPOT files     |2              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_DIEL.dat`       |Generated dielectric profile (ε\ :sub:`11` ε\ :sub:`22`|3              |
|                        |ε\ :sub:`33`) along the normal axis to the surface     |               |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_DXCHG.dat`      |Planar average of extra charge (neutral and charged    |`3`__          |
|                        |difference) CHGCAR file in X direction                 |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_DXPOT.dat`      |Planar average of extra charge (neutral and charged    |`3`__          |
|                        |difference) LOCPOT file in X direction                 |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_DYCHG.dat`      |Planar average of extra charge (neutral and charged    |`3`__          |
|                        |difference) CHGCAR file in Y direction                 |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_DYPOT.dat`      |Planar average of extra charge (neutral and charged    |`3`__          |
|                        |difference) LOCPOT file in Y direction                 |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_DZCHG.dat`      |Planar average of extra charge (neutral and charged    |`3`__          |
|                        |difference) CHGCAR file in Z direction                 |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_DZPOT.dat`      |Planar average of extra charge (neutral and charged    |`3`__          |
|                        |difference) LOCPOT file in Z direction                 |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_M.CHGCAR`       |CHGCAR of the Gaussian model                           |2              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_M.LOCPOT`       |LOCPOT of the Gaussian model                           |2              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_MXCHG.dat`      |Planar average of model charge in X direction          |`3`__          |
|                        |                                                       |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_MXPOT.dat`      |Planar average of model potential in X direction       |`3`__          |
|                        |                                                       |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_MYCHG.dat`      |Planar average of model charge in Y direction          |`3`__          |
|                        |                                                       |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_MYPOT.dat`      |Planar average of model potential in Y direction       |`3`__          |
|                        |                                                       |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_MZCHG.dat`      |Planar average of model charge in Z direction          |`3`__          |
|                        |                                                       |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_MZPOT.dat`      |Planar average of model potential in Z direction       |`3`__          |
|                        |                                                       |               |
|                        |                                                       |__ avg_note_   |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_NXCHG.dat`      |Planar average of neutral CHGCAR file in X direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_NXPOT.dat`      |Planar average of neutral LOCPOT file in X direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_NYCHG.dat`      |Planar average of neutral CHGCAR file in Y direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_NYPOT.dat`      |Planar average of neutral LOCPOT file in Y direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_NZCHG.dat`      |Planar average of neutral CHGCAR file in Z direction   |3              |
+------------------------+-------------------------------------------------------+---------------+
|`slabcc_NZPOT.dat`      |Planar average of neutral LOCPOT file in Z direction   |3              |
+------------------------+-------------------------------------------------------+---------------+

.. _avg_note:

**Note:** The planar averaged potential and charge files corresponding to the normal direction will be written in verbosity = 1

===
FAQ
===

1. **How to obtain the CHGCAR and LOCPOT files from VASP calculations?** You can add the following tags to your INCAR file to get the LOCPOT and CHGCAR files::

    LVTOT = .TRUE.
    LVHAR = .TRUE.
    LCHARG = .TRUE.

 After obtaining the files for your charged system, do the calculation again *without relaxing (changing) the geometry* to get the necessary files for the neutralized system.

2. **Do I need to perform spin polarized calculation in VASP?**  Although, the slabcc only reads the sum of both spins, but for proper description of the charge distribution in your system you may need to perform spin polarized calculation.

3. **How can I speed-up the model parameters optimization process?** You can try using a different optimization algorithm or improve the initial guess for the model parameters to speed-up the optimization. As a last resort, you can also use a smaller computation grid for the optimization (``optimize_grid_x < 1``), or increase the optimization convergence criteria (``optimize_tolerance``) to speed up the process but the accuracy of the obtained results in these cases must be always checked.

4. **Why do I need to provide an initial guess for the parameters which will be optimized?** The optimization algorithms used in slabcc are local error minimization algorithms. Their success and performance highly depend on the initial guess for the provided parameters.

5. **How should I decide on the initial guess for the parameters which will be optimized?** As a rule of thumb, start by a single Gaussian charge as your model. Set its position to your expected position of the charge localization. Use the location of the surface atoms as the interface position. You can use the “-d” switch in the command line (./slabcc -d) to just generate the CHGCAR and the LOCPOT file for the extra charge and their planar averages without shifting the input files to the `slab_center`. These files will guide you on how to provide the initial guess for the input parameters.

6. **Can I turn off the optimization for the input parameters?** Yes. But optimization ensures the model charge mimics the original localized charge in large distances as close as possible. If you turn off the optimization, you must be aware of the possible side-effects and definitely `check your results`__.

__ check_

7. **Can I run the slabcc on a computational cluster?** Yes. BUT… Although slabcc hugely benefits from the multicore architecture of the computation nodes using OpenMP, it has not yet been parallelized using MPI. Therefore, It won’t use more than one machine at a time.

8. **Is the slabcc free? Can I use its source code in my own software?** slabcc is released under the 2-Clause BSD license_ which permits this software to be modified, redistributed and/or used for commercial purposes provided that the source retains the original copyright owner's name (University of Bremen, M. Farzalipour Tabriz) and full text of the license (LICENSE.txt)

9. **How accurate are the slabcc results?** The accuracy of the final results depends on various factors including the accuracy/grid-size of the input files and provided input parameters. The optimization algorithm used for parameters estimation is a non-linear local optimizer which means that the result will highly depend on its initial conditions. Models with different number of Gaussian charges have different accuracy and may be compared with caution. In case of the models with multiple charges, the results must be vigorously checked. You must always do your own testing before using the results. There are a few `known issues and limitations`_ to the slabcc code and its algorithm. Also keep in mind that this is a free software and as the license_ explicitly mentions: there is absolutely no warranty for its fitness for any particular purpose.

.. _check:

10. **How can I check the slabcc results?** slabcc can calculate the planar averaged potential and charge files for the extra charge in the input files and the model Gaussian charge. You should compare the model charge distribution and potential specially in the direction normal to the surface and compare them to the original VASP results. For example, if z is the normal direction in your slab model (``normal_direction = z``), then you should compare `slabcc_MZCHG.dat` and `slabcc_MZPOT.dat`, with `slabcc_DZCHG.dat` and `slabcc_DZPOT.dat`, respectively. Check `the files table`_ for complete list of the output files.

 Another method to test the effectiveness of the charge correction is to increase the thickness of the vacuum in your slab model and check the energies. If the charge correction is done properly, the energy values must be independent of the (adequately large) vacuum thickness.

.. _cite:

11. **How should I cite the slabcc?** Please cite the slabcc as: (You can `download the citation in the RIS format from here <https://www.sciencedirect.com/sdfe/arp/cite?pii=S0010465519300700&format=application%2Fx-research-info-systems&withabstract=true>`_!)

 Meisam Farzalipour Tabriz, Bálint Aradi, Thomas Frauenheim, Peter Deák, *SLABCC: Total energy correction code for charged periodic slab models*, Computer Physics Communications, Vol. 240C (2019), pp. 101-105, DOI: `10.1016/j.cpc.2019.02.018 <https://doi.org/10.1016/j.cpc.2019.02.018>`_
  
==================================
Known issues and limitations
==================================
- Shape of the VASP files cell is limited to orthogonal cells.
- Maximum line length of the input file (slabcc.in) is 4000 bytes.
- Bessel expansion of the Poisson equation cannot be used for the calculation of isolated energies for the 2D models with anisotropic in-plane screening, trivariate Gaussian model change, or the models which are not surrounded by the vacuum (diel_out > 1). Extrapolation method must be used in these cases.

==========================
Release history highlights
==========================
* 2019-06-13: version 0.8 - OO redesign
* 2019-05-14: version 0.7 - Added discretization error mitigation
* 2019-04-04: version 0.6 - Added trivariate Gaussian model charge and selective charge optimization support
* 2019-03-18: version 0.5 - Added 2D model support
* 2018-10-10: version 0.4 - Added spdlog and several user interface and performance improvements
* 2018-07-29: version 0.3 - First public release

===========================
Copyright and attributions
===========================
Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz

The source codes and all the documentations are available under The 2-Clause BSD License. For more information see license_.

| This code uses several open source components each of which are located under a separate sub-directory in the `src/`. The copyright of these libraries belong to their respective owners. Any modification made to those codes is also published under the same license. We acknowledge and are grateful to these developers and maintainers for their valuable contributions to this software and more importantly to the free software society.
| The attributions are also present in the binary file and can be accessed using the command-line parameters.

Included third-party components
-------------------------------

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

License
-------
Copyright (c) 2018-2019, University of Bremen, M. Farzalipour Tabriz

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
