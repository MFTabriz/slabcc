=============
Introduction
=============
SLABCC calculates an *a posteriori* energy correction for charged slab models under 3D periodic boundary conditions (PBC) based on the method proposed in:

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

| SLABCC have been initially developed for the `Bremen Center for Computational Materials Science (BCCMS) <https://www.uni-bremen.de/bccms>`_
