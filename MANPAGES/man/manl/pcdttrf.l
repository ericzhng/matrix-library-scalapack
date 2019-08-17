.TH PCDTTRF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCDTTRF - compute a LU factorization of an N-by-N complex tridiagonal diagonally dominant-like distributed matrix A(1:N, JA:JA+N-1)
.SH SYNOPSIS
.TP 20
SUBROUTINE PCDTTRF(
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
COMPLEX
AF( * ), D( * ), DL( * ), DU( * ), WORK( * )
.SH PURPOSE
PCDTTRF computes a LU factorization of an N-by-N complex tridiagonal diagonally dominant-like distributed matrix A(1:N, JA:JA+N-1). Reordering is used to increase parallelism in the factorization.
This reordering results in factors that are DIFFERENT from those
produced by equivalent sequential codes. These factors cannot
be used directly by users; however, they can be used in
.br
subsequent calls to PCDTTRS to solve linear systems.
.br

The factorization has the form
.br

        P A(1:N, JA:JA+N-1) P^T = L U
.br

where U is a tridiagonal upper triangular matrix and L is tridiagonal
lower triangular, and P is a permutation matrix.
.br

