It solves the Lippmann-Schwinger equation via volume integral equation using GMRES, an iterative method. The matrix-vector products encountered in each of its iteration have been computed using DAFMM (Directional Algebraic Fast Multipole Method). All the low rank factorisations encountered were formed using Nested Cross Approximation (NCA).

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

m: number of iterations of GMRES to be run

restart: number of restarts of Restarted GMRES

nLevelsUniform: In case a uniform tree is to be wanted, input the number of levels to be present in the uniform tree. By default the code makes the tree adaptive. (To make it uniform appropriate changes in the code are to be made.)

Qchoice: Choice for type of contrast functions <br />
0 - Gaussian <br />
1 - Multiple Gaussians <br />
2 - Plasma <br />
3 - Flat Bump <br />
4 - Cavity <br />
5 - Lens <br />

To run the code input in terminal:

make -f Makefile2D_GMRES.mk clean

make -f Makefile2D_GMRES.mk

./testFMM2D_GMRES 16 6 5 0.5 40.0 1 6 8 10 1 5 0

The output looks like:

Wavenumber:		40 <br />
Wavelength:		0.15708 <br />
no. of full cycles:	31.4159 <br />
nLevels: 5 <br />
level_LFR: 4 <br />
Number of particles: 5328 <br />
O;	j: 5	Nboxes: 48	rows,cols: 576,36	Crank: 18 <br />
I;	j: 5	Nboxes: 48	rows,cols: 36,540	Crank: 21 <br />
O;	j: 4	Nboxes: 64	rows,cols: 1665,70	Crank: 21 <br />
I;	j: 4	Nboxes: 64	rows,cols: 81,1644	Crank: 23 <br />
O;	j: 3	Nboxes: 64	k: 30	rows,cols: 2360,36	Crank: 13 <br />
I;	j: 3	Nboxes: 64	k: 30	rows,cols: 36,2300	Crank: 12 <br />
O;	j: 2	Nboxes: 16	k: 0	rows,cols: 0,35	Crank: 4 <br />
I;	j: 2	Nboxes: 16	k: 1	rows,cols: 16,0	Crank: 4 <br />

Time taken to assemble: 285.211

Time taken by GMRES Solver 1,10: 0.277782 <br />
Err in GMRES: 4.59412e-06


Time taken by GMRES Solver 1,50: 1.47135 <br />
Err in GMRES: 1.40294e-15


Time taken by GMRES Solver 1,100: 2.90013 <br />
Err in GMRES: 1.30732e-15


Time taken by GMRES Solver 1,200: 6.22192 <br />
Err in GMRES: 1.39836e-15


Time taken by GMRES Solver 1,500: 19.1974 <br />
Err in GMRES: 1.40164e-15


Time taken by GMRES Solver 2,10: 0.751761 <br />
Err in GMRES: 1.28443e-10


Time taken by GMRES Solver 2,50: 2.63477 <br />
Err in GMRES: 8.96364e-16


Time taken by GMRES Solver 2,100: 5.29393 <br />
Err in GMRES: 7.94185e-16


Time taken by GMRES Solver 2,200: 11.7579 <br />
Err in GMRES: 8.37373e-16


Time taken by GMRES Solver 2,500: 38.9819 <br />
Err in GMRES: 7.82761e-16

nLevels: 5 <br />
level_LFR: 4 <br />
Number of particles: 5328 <br />
O;	j: 5	Nboxes: 48	rows,cols: 540,36	Crank: 12 <br />
I;	j: 5	Nboxes: 48	rows,cols: 36,540	Crank: 14 <br />
O;	j: 4	Nboxes: 64	rows,cols: 1635,36	Crank: 15 <br />
I;	j: 4	Nboxes: 64	rows,cols: 36,906	Crank: 18 <br />
O;	j: 3	Nboxes: 64	k: 22	rows,cols: 1726,36	Crank: 18 <br />
I;	j: 3	Nboxes: 64	k: 22	rows,cols: 36,1648	Crank: 20 <br />
O;	j: 2	Nboxes: 16	k: 1	rows,cols: 0,69	Crank: 7 <br />
I;	j: 2	Nboxes: 16	k: 0	rows,cols: 60,0	Crank: 8 <br />

Time taken for MatVec product: 0.011503

Time taken for field computation: 38.9934
