It solves the Lippmann-Schwinger equation via volume integral equation using GMRES, an iterative method. A preconditioner is available that is based on HODLR direct solver. The iterative solver that uses preconditioner is termed Hybrid solver, and that which does not use preconditioner is simply termed GMRES solver. The matrix-vector products encountered in each of GMRES' iteration have been computed using DAFMM (Directional Algebraic Fast Multipole Method). All the low rank factorizations encountered were formed using Nested Cross Approximation (NCA). For details on HODLR library please refer [[1]](#1). This code can be used to validate the methods used, to time-profile the code, to obtain the error and scattered field plots.

To run GMRES solver set the INPUT_FILE variable in CMakeLists.txt file as LippmannSchwinger/main/hybridTrueValidation_gmres.cpp

To run Hybrid solver set the INPUT_FILE variable in CMakeLists.txt file as LippmannSchwinger/main/hybridTrueValidation_hybrid.cpp

To run the project make sure boost, eigen, openmp libraries are linked.

Both of them takes the following inputs at run-time:

nCones_LFR: number of cones a box at the highest level (in number) of high frequency regime is to be divided into.

nChebNodes: number of gridPoints in leaf box in 1D

treeAdaptivity: tolerance requested from adaptive discretization of the tree.

L: half side length of the square domain.

kappa: wavenumber

yes2DFMM: 1 for Directional AFMM; and 0 for AFMM.

degreeOfBases: number of polynomials used in a box for discretization of Lippmann-Schwinger equation

TOL_POW: tolerance used for ACA - to be inputed in negative powers of 10.

m: number of iterations of GMRES to be run - currently not used (the current version of GMRES stops when the tolerance is reached)

restart: number of restarts of Restarted GMRES - currently not used

pre_conditioner_tolerance: tolerance used for constructing low rank factorizations of the low rank sub-blocks of preconditioner matrix. (currently not used)

preconditioner_target_rank: rank used for constructing low rank factorisations of the low rank sub-blocks of preconditioner matrix.

Qchoice: Choice for type of contrast functions <br />
0 - Gaussian <br />
1 - Multiple Gaussians <br />
2 - Plasma <br />
3 - Flat Bump <br />
4 - Cavity <br />
5 - Lens <br />

nLevelsUniform: used in construction of error function

A sample run of the code in terminal is here:

cmake.. <br />
make <br />

./exec1 16 8 8 0.5 40.0 1 6 10 10 1 0.0001 10 0 4

After a successful run of Hybrid solver, the results get stored in the following files:

build/result/result_0_8_8_40_10/errU_scattered_Hybrid_imag/ : error in imaginary part of scattered field <br />
build/result/result_0_8_8_40_10/errU_scattered_Hybrid_real : error in real part of scattered field <br />
build/result/result_0_8_8_40_10/solutionI : imaginary part of scattered field <br />
build/result/result_0_8_8_40_10/solutionR : real part of scattered field <br />
build/result/result_0_8_8_40_10/gridPointsY : target gridpoints <br />
build/result/result_0_8_8_40_10/gridPointsX : charge gridpoints (both target and charge gridpoints are same in the current version) <br />


## References
<a id="1">[1]</a>
Ambikasaran, S., Singh, K. R., & Sankaran, S. S. (2019). Hodlrlib: A library for hierarchical matrices. Journal of Open Source Software, 4(34), 1167.
