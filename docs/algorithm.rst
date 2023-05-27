Algorithm
----------
* slabcc reads the `VASP <https://www.vasp.at>`_ output files (CHGCAR and LOCPOT), and calculates the extra charge distribution and the potential due to the extra charge. All the input files should correspond to the same geometry.
* The extra charge is modeled by the sum of Gaussian charge distributions as follows:

.. math::

 \rho(r) = \sum_{i}\frac{q_i}{\sigma_{i}^{3}(2\pi)^{3/2}} \exp \left ({- \frac{r_{i}^{2}}{2\sigma_{i}^{2}} } \right )

or `trivariate Gaussian <https://mathworld.wolfram.com/TrivariateNormalDistribution.html>`_ charge distributions as:

.. math::

 \rho(r) = \sum_{i}\frac{q_i}{\sigma_{i,x}\sigma_{i,y}\sigma_{i,z}(2\pi)^{3/2}} \exp \left ({- \frac{r_{i,x}^{2}}{2\sigma_{i,x}^{2}} - \frac{r_{i,y}^{2}}{2\sigma_{i,y}^{2}}- \frac{r_{i,z}^{2}}{2\sigma_{i,z}^{2}} } \right )

which gives a charge distribution normalized to q\ :sub:`i` \ with a standard deviation of σ\ :sub:`i` \ (``charge_sigma``) for each Gaussian distribution i centered at r\ :sub:`i` \ (``charge_position``).

* The generated Gaussian model charge will be embedded in a dielectric medium with a profile of the form:

.. math::
  \epsilon (k) =  \frac{\epsilon_2-\epsilon_1}{2} \text{erf}\left(\frac{k-k_0 }{\beta}\right)+\frac{\epsilon_2+\epsilon_1}{2}

where k\ :sub:`0` \ is the interface position in Cartesian k-direction, ε\ :sub:`1` \ and ε\ :sub:`2` \ are dielectric tensors on either side of the interface (``diel_in`` & ``diel_out``) and β (``diel_taper``) defines the smoothness of transition assuming anisotropic dielectric tensor as:

.. math::
 \epsilon =
  \left | \begin{matrix}
    \epsilon_{11} & 0 & 0 \\
    0 & \epsilon_{22} & 0 \\
    0 & 0&  \epsilon_{33}
 \end{matrix}  \right |

* The potential due to the charge distribution ρ under 3D PBC embedded in the dielectric medium ε(k) is calculated by solving the Poisson equation in Fourier space:

.. math::
  \epsilon(k) \nabla^2 V(r)+\frac{\partial}{\partial k} \epsilon(k)\frac{\partial}{\partial k}V(r) = -\rho(r)

* A non-linear optimization routine minimizes the difference between our calculated V(r) for the model charge and the V resulted from the VASP calculation by changing the position of the model Gaussian charge, its width, and the position of the slab interfaces.

* The E\ :sub:`periodic` is calculated as:

.. math::
  E = \frac{1}{2} \int V(r) \rho(r) \, dr

* E\ :sub:`isolated` is calculated the same way as E\ :sub:`periodic` but with extrapolation of the fixed model charge embedded in an infinitely large dielectric medium. For the bulk and slab models, the extrapolation is done linearly. For the monolayer models (2D systems), the following equation is used for the extrapolation [`10.1103/PhysRevX.8.039902 <https://doi.org/10.1103/PhysRevX.8.039902>`_]:

.. math::
  E = c_0 + c_1 x + c_2 x^2 + d e^{-c_3 x}

where c\ :sub:`i` are the fitting parameters and

.. math::
  d =  \frac{c_1 - \frac{\partial E_M}{\partial x}}{c_3}

guarantees the correct energy gradient at x(=1/α)→0. E\ :sub:`M` being the Madelung energy.

* ΔV is calculated at the position least affected by the model charge.

More information about the algorithms and the implementation details can be found `here`__.

__ cite_
