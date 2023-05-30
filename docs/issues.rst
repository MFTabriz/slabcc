==================================
Known issues and limitations
==================================
- Only orthogonal cells are supported.
- The maximum line length of the input file (slabcc.in) is 4000 bytes.
- Bessel expansion of the Poisson equation cannot be used for the calculation of isolated energies for the 2D models with anisotropic in-plane screening, trivariate Gaussian model change, or the models that are not surrounded by the vacuum (diel_out > 1). The extrapolation method must be used in these cases.
