It solves the Lippmann-Schwinger equation via volume integral equation using HODLR, a direct solver. For details on HODLR library please refer [[1]](#1).

To run the project make sure boost, eigen, openmp libraries are linked.

It takes the following inputs at run-time:

nCones_LFR: number of cones a box at the highest level (in number) of high frequency regime is to be divided into.

nChebNodes: number of gridPoints in leaf box in 1D

treeAdaptivity: tolerance requested from adaptive discretization of the tree.

L: half side length of the square domain.

kappa: wavenumber

yes2DFMM: 1 for Directional AFMM; and 0 for AFMM.

degreeOfBases: number of polynomials used in a box for discretization of Lippmann-Schwinger equation

TOL_POW: tolerance used for ACA - to be inputed in negative powers of 10.

nLevelsUniform: In case a uniform tree is to be wanted, input the number of levels to be present in the uniform tree. By default the code makes the tree adaptive. (To make it uniform appropriate changes in the code are to be made.)

Qchoice: Choice for type of contrast functions <br />
0 - Gaussian <br />
1 - Multiple Gaussians <br />
2 - Plasma <br />
3 - Flat Bump <br />
4 - Cavity <br />
5 - Lens <br />

To run the code input in terminal:

cmake.. <br />
make <br />

./exec1 16 6 5 1.0 40.0 1 6 8 4 0

The output looks like:

Wavenumber:		40 <br />
Wavelength:		0.15708 <br />
no. of full cycles:	62.8319 <br />
nLevels: 6 <br />
level_LFR: 5 <br />
Number of particles: 7056 <br />
Time taken by kernel initialisation : 1.36059 <br />
========================= Problem Parameters ========================= <br />
Matrix Size                        :7056 <br />
Leaf Size                          :36 <br />
Dimensionality                     :2 <br />
Tolerance                          :1e-08 <br />

========================= Assembly Time ========================= <br />
assembleTime: 46.9169 <br />
========================= Factorization ========================= <br />
Time to factorize HODLR form       :1.2943 <br />
========================= Solving ========================= <br />
Time to solve HODLR form           : 0.023343 <br />
Err in HODLR Phi : 7.05145e-15 <br />
Time taken by HODLR solver: 1.31765 <br />
========================= Matrix-Vector Multiplication ========================= <br />
nLevels: 6 <br />
level_LFR: 5 <br />
Number of particles: 7056 <br />
Time taken for field computation      :1.35282


## References
<a id="1">[1]</a>
Ambikasaran, S., Singh, K. R., & Sankaran, S. S. (2019). Hodlrlib: A library for hierarchical matrices. Journal of Open Source Software, 4(34), 1167.
