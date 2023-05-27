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
  interfaces_optimized = 0.942000748357 0.455672787711
  charge_sigma_optimized = 1.4132676877
  charge_position_optimized = 0.501460639345 0.50145532106 0.385476689493;

  [Results]
  dV = -0.00291385176718
  E_periodic of the model charge = 2.0404453156f
  E_isolated of the model charge = 2.59716677886
  Energy correction for the model charge (E_iso-E_per-q*dV) = 0.559635314929

Planar average files are written as the double column in plain text format. The first column represents the coordinates along the axis (in Angstrom) and the second column is the planar average value. The files are named as: "slabcc_{1}{2}{XXX}.dat" where:

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
