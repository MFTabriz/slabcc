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
