======================
Input parameters
======================
slabcc reads all its parameters from the input file (by default, ``slabcc.in``). You can change the input file's name using the `command-line parameters`_.
The input file is processed as follows:

- Lines starting with # will be treated as comments. Inline comments are also allowed.
- Double quotation marks will be removed from the strings.
- A warning will be issued for any unidentified parameter.
- A warning will be issued for the use of the deprecated parameters.
- All the coordinates must be in fractional form [0-1].
- Boolean (True/False) parameters can also be declared as 0/1, on/off, yes/no, .true./.false.
- Parameter names can be written in small or capital letters.
- For vectors and matrices, columns are separated by a " "(space), while rows are separated by a ";" (semicolon).
- Lines starting with a space " " will be treated as the continuation of the last parameter's value.
- Subsequent definitions for any parameter will be concatenated with the existing definition.

 
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
|                              |Gaussian charges in arc degrees (-90, 90)              |       0       |
| ``charge_rotation``          |                                                       |               |
|                              |``charge_rotation = 0 45 0``                           |               |
|                              |                                                       |               |
|                              |``charge_rotation = 45 0 0; 0 -10 0``                  |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |width of each localized Gaussian charge. It can be 1   |               |
|                              |or, in the case of trivariate models, 3 parameters per |               |
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
|                              |inside the slab. If only a single value is given,      |               |
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
|                              |Extrapolation grid size multiplier. The number of      |               |
|                              |grid points in each direction will be multiplied by    |               |
|                              |this value.                                            |               |
|                              |                                                       |               |
|                              |extrapolate_grid_x > 1 will use a larger grid in the   |               |
|``extrapolate_grid_x``        |extrapolations, which will increase its accuracy but   |       1       |
|                              |will requires more memory and computational power.     |               |
|                              |                                                       |               |
|                              |extrapolate_grid_x = 1 will use the same grid size as  |               |
|                              |the VASP input files in the extrapolation.             |               |
|                              |                                                       |               |
|                              |extrapolate_grid_x < 1 will use a smaller grid for the |               |
|                              |extrapolations, which increases the speed and decreases|               |
|                              |the memory usage, but the energies for the higher      |               |
|                              |orders of the extrapolation may not be accurate!       |               |
|                              |                                                       |               |
|                              |``extrapolate_grid_x = 1.8``                           |               |
+------------------------------+-------------------------------------------------------+---------------+
|                              |Number of extrapolation steps in the calculation of    |10: for 2D     |
| ``extrapolate_steps_number`` |E\ :sub:`isolated` \                                   |models         |
|                              |                                                       |               |
|                              |                                                       |4: for the rest|
+------------------------------+-------------------------------------------------------+---------------+
|                              |Size of extrapolation steps with respect to the initial|1: for 2D      |
| ``extrapolate_steps_size``   |supercell size                                         |models         |
|                              |                                                       |               |
|                              |                                                       |0.5: for the   |
|                              |                                                       |rest           |
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
| ``optimize``                 |Optimizer master switch. Upon deactivation, it takes   |               |
|                              |precedence over all the other optimization options:    |    true       |
|                              |``optimize_charge_fraction``,                          |               |
|                              |``optimize_charge_position``,                          |               |
|                              |``optimize_charge_rotation``,                          |               |
|                              |``optimize_charge_sigma``, and ``optimize_interfaces`` |               |
|                              |                                                       |               |
|                              |**true**: evaluate each of the optimization switches   |               |
|                              |individually                                           |               |
|                              |                                                       |               |
|                              |**false**: deactivate all optimization switches        |               |
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

.. [#] M.J.D. Powell, `The BOBYQA algorithm for bound constrained optimization without derivatives <https://www.damtp.cam.ac.uk/user/na/NA_papers/NA2009_06.pdf>`_, Department of Applied Mathematics and Theoretical Physics, Cambridge England, technical report NA2009/06 (2009).
.. [#] M.J.D. Powell, `Direct search algorithms for optimization calculations <https://doi.org/10.1017/S0962492900002841>`_, Acta Numerica, Vol. 7(1998) pp. 287-336
.. [#] T.H. Rowan, `Functional Stability Analysis of Numerical Algorithms <https://dl.acm.org/doi/book/10.5555/100816>`_, Ph.D. thesis, Department of Computer Sciences, University of Texas at Austin, 1990.
.. [#] Each verbosity level includes all the outputs from the lower verbosity options. Check `the files table`_ for complete list of the output files.
