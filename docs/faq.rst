===
FAQ
===

1. **How do I obtain the CHGCAR and LOCPOT files from VASP calculations?** You can add the following tags to your INCAR file to get the LOCPOT and CHGCAR files::

    LVTOT = .TRUE.
    LVHAR = .TRUE.
    LCHARG = .TRUE.

 After obtaining the files for your charged system, do the calculation again *without relaxing (changing) the geometry* to get the necessary files for the neutralized system.

2. **Do I need to perform spin-polarized calculations in VASP?**  Although the slabcc only reads the sum of both spins, for a proper description of the charge distribution in your system, you may need to perform a spin-polarized calculation.

3. **How can I speed up the model parameter optimization process?** You can try using a different optimization algorithm or improving the initial guess for the model parameters to speed up the optimization. As a last resort, you can also use a smaller computation grid for the optimization (``optimize_grid_x < 1``), or increase the optimization convergence criteria (``optimize_tolerance``) to speed up the process, but the accuracy of the obtained results may not be adequate. Ideally, you should use the optimzied parameters on a coarser grid as the initial guess for optimization on a finer (initial) grid.

4. **Why do I need to provide an initial guess for the parameters that will be optimized?** The optimization algorithms used in slabcc are local error minimization algorithms. Their success and performance highly depend on the initial guess for the provided parameters.

5. **How should I decide on the initial guess for the parameters that will be optimized?** As a rule of thumb, start with a single Gaussian charge as your model. Set its position to the expected position of the charge localization. Use the location of the surface atoms as the interface position. You can use the "-d" switch in the command line (./slabcc -d) to just generate the CHGCAR and the LOCPOT files for the extra charge and their planar averages without shifting the input files to the `slab_center`. These files will guide you on how to provide the initial guess for the input parameters.

6. **Can I turn off the optimization for the input parameters?** Yes. But optimization ensures the model charge mimics the original localized charge at large distances as closely as possible. If you turn off the optimization, you must be aware of the possible side-effects and definitely `check your results`__.

__ check_

7. **Can I run the slabcc on a computational cluster?** Yes. BUT… Although slabcc hugely benefits from the multicore architecture of the computation nodes using OpenMP, it has not yet been parallelized using MPI. Therefore, it won’t use more than one machine at a time.

8. **Is the slabcc free? Can I use its source code in my own software?** slabcc is released under the 2-Clause BSD license_ which permits this software to be modified, redistributed, and/or used for commercial purposes provided that the source retains the original copyright owner's name (University of Bremen, M. Farzalipour Tabriz) and full text of the license (LICENSE.txt)

9. **How accurate are the slabcc results?** The accuracy of the final results depends on various factors, including the accuracy and grid size of the input files and the provided input parameters. The optimization algorithm used for parameter estimation is a non-linear local optimizer, which means that the result will highly depend on its initial conditions. Models with different numbers of Gaussian charges have different accuracy and may be compared with caution. In the case of models with multiple charges, the results must be vigorously checked. You must always do your own testing before using the results. There are a few `known issues and limitations`_ to the slabcc code and its algorithm. Also keep in mind that this is free software, and as the license_ explicitly mentions, there is absolutely no warranty for its fitness for any particular purpose.

.. _check:

10. **How can I check the slabcc results?** slabcc can calculate the planar averaged potential and charge files for the extra charge in the input files and the model Gaussian charge. You should compare the model charge distribution and potential, especially in the direction normal to the surface, and compare them to the original VASP results. For example, if z is the normal direction in your slab model (``normal_direction = z``), then you should compare ``slabcc_MZCHG.dat`` and ``slabcc_MZPOT.dat``, with ``slabcc_DZCHG.dat`` and ``slabcc_DZPOT.dat``, respectively. Check `the files table`_ for a complete list of the output files.

 Another method to test the effectiveness of the charge correction is to increase the thickness of the vacuum in your slab model and check the charge-corrected total energies. If the charge correction is done properly, the energy values must be independent of the (adequately large) vacuum thickness.

.. _cite:

11. **How should I cite the slabcc?** Please cite the slabcc as: (You can `download the citation in the RIS format from here <https://www.sciencedirect.com/sdfe/arp/cite?pii=S0010465519300700&format=application%2Fx-research-info-systems&withabstract=true>`_!)

 Meisam Farzalipour Tabriz, Bálint Aradi, Thomas Frauenheim, Peter Deák, *SLABCC: Total energy correction code for charged periodic slab models*, Computer Physics Communications, Vol. 240C (2019), pp. 101-105, DOI: `10.1016/j.cpc.2019.02.018 <https://doi.org/10.1016/j.cpc.2019.02.018>`_

12. **How can I extract the files in slabcc_data.tar.xz?** You can use `Tar <https://www.gnu.org/software/tar/>`_ + `XZ Utils <https://tukaani.org/xz/>`_ as:

    tar -xvf slabcc_data.tar.xz

 Alternatively, you can use `WinRAR <https://www.rarlab.com>`_ or `7zip <https://www.7-zip.org>`_.

13. **Something is not working! What should I do?**

 * If you need help with compiling the code or running it on a cluster, please contact your `system administrator <https://en.wikipedia.org/wiki/System_administrator>`_.
 * If you have found a bug in the code, please report it `here <https://codeberg.org/meisam/slabcc/issues/new>`__.

