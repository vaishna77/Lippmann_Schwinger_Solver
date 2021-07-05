It solves the Lippmann-Schwinger equation via volume integral equation using GMRES, an iterative method. HODLR, a direct solver with low precision is used as pre-conditioner. The matrix-vector products encountered in each of GMRES' iteration have been computed using DAFMM (Directional Algebraic Fast Multipole Method). All the low rank factorisations encountered were formed using Nested Cross Approximation (NCA). For details on HODLR library please refer [[1]](#1).

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

pre_conditioner_tolerance: tolerance used for constructing low rank factorisations of the low rank sub-blocks of preconditioner matrix.

preconditioner_target_rank: rank used for constructing low rank factorisations of the low rank sub-blocks of preconditioner matrix. The basis of the low rank sub-blocks are constructed until the pre_conditioner_tolerance or preconditioner_target_rank is reached.

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

./exec1 16 6 5 0.5 40.0 1 6 8 10 1 5 25 0

Wavenumber:		40 <br />
Wavelength:		0.15708 <br />
no. of full cycles:	31.4159 <br />
Begining GMRES initialisation <br />
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
I;	j: 2	Nboxes: 16	k: 1	rows,cols: 16,0	Crank: 4

Time taken to assemble: 271.235 <br />
GMRES initialisation done <br />
nLevels: 5 <br />
level_LFR: 4 <br />
Number of particles: 5328 <br />
Time taken to preconditioner factorization: 0.026538 <br />
HODLR pre-conditioner initialisation done <br />
Err in GMRES: 3.53921e-07 <br />
Time taken by GMRES Solver 1,10: 0.313262

Err in GMRES: 2.17642e-13 <br />
Time taken by GMRES Solver 2,10: 0.515279

Err in GMRES: 2.62095e-14 <br />
Time taken by GMRES Solver 1,50: 1.2031

Err in GMRES: 1.48022e-13 <br />
Time taken by GMRES Solver 1,100: 2.51867

------------------------------------------------------ <br />
Err in GMRES without precond: 4.59412e-06 <br />
Time taken by GMRES Solver without precond 1,10: 0.230691

Err in GMRES without precond: 1.28443e-10 <br />
Time taken by GMRES Solver without precond 2,10: 0.437019

Err in GMRES without precond: 1.40294e-15 <br />
Time taken by GMRES Solver without precond 1,50: 1.0451

Err in GMRES without precond: 1.30732e-15 <br />
Time taken by GMRES Solver without precond 1,100: 2.26183

Err in GMRES without precond: 7.94185e-16 <br />
Time taken by GMRES Solver without precond 2,100: 4.72124

Err in GMRES without precond: 1.39836e-15 <br />
Time taken by GMRES Solver without precond 1,200: 5.07559

Err in GMRES without precond: 8.37373e-16 <br />
Time taken by GMRES Solver without precond 2,200: 14.5496

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
I;	j: 2	Nboxes: 16	k: 0	rows,cols: 60,0	Crank: 8

Time taken for MatVec product: 0.009641

Time taken for field computation: 10.6825

## References
<a id="1">[1]</a>
Ambikasaran, S., Singh, K. R., & Sankaran, S. S. (2019). Hodlrlib: A library for hierarchical matrices. Journal of Open Source Software, 4(34), 1167.
