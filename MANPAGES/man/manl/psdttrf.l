.TH PSDTTRF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSDTTRF - compute a LU factorization of an N-by-N real tridiagonal diagonally dominant-like distributed matrix A(1:N, JA:JA+N-1)
.SH SYNOPSIS
.TP 20
SUBROUTINE PSDTTRF(
N, DL, D, DU, JA, DESCA, AF, LAF, WORK, LWORK,
INFO )
.TP 20
.ti +4
INTEGER
INFO, JA, LAF, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
AF( * ), D( * ), DL( * ), DU( * ), WORK( * )
.SH PURPOSE
PSDTTRF computes a LU factorization of an N-by-N real tridiagonal diagonally dominant-like distributed matrix A(1:N, JA:JA+N-1). Reordering is used to increase parallelism in the factorization.
This reordering results in factors that are DIFFERENT from those
produced by equivalent sequential codes. These factors cannot
be used directly by users; however, they can be used in
.br
subsequent calls to PSDTTRS to solve linear systems.
.br

The factorization has the form
.br

        P A(1:N, JA:JA+N-1) P^T = L U
.br

where U is a tridiagonal upper triangular matrix and L is tridiagonal
lower triangular, and P is a permutation matrix.
.br

