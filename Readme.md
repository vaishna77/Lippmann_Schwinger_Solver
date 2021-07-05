This repository has 3 projects to solve the Lippmann-Schwinger equation via volume integral equation using 3 techniques.

1. HODLR: HODLR based direct solver.

2. GMRES: GMRES based iterative solver. The matrix-vector products encountered in each of its iteration have been computed using DAFMM (Directional Algebraic Fast Multipole Method). All the low rank factorisations encountered were formed using Nested Cross Approximation (NCA).

3. Hybrid: GMRES based iterative solver with HODLR as pre-conditioner.

For details on matrix assembly please refer the article [[1]](#1).

The inputs to the codes are to be given at run-time, the details of which have been given in the Readme.md files of the respective projects.

The value of Phi and the real part of the field get stored in a directory called result, which gets created at run-time.

## References
<a id="1">[1]</a>
Ambikasaran, S., Borges, C., Imbert-Gerard, L. M., & Greengard, L. (2016). Fast, adaptive, high-order accurate discretization of the Lippmann--Schwinger equation in two dimensions. SIAM Journal on Scientific Computing, 38(3), A1770-A1787.
