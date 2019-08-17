.TH PZPTTRF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZPTTRF - compute a Cholesky factorization of an N-by-N complex tridiagonal symmetric positive definite distributed matrix A(1:N, JA:JA+N-1)
.SH SYNOPSIS
.TP 20
SUBROUTINE PZPTTRF(
N, D, E, JA, DESCA, AF, LAF, WORK, LWORK,
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
COMPLEX*16
AF( * ), E( * ), WORK( * )
.TP 20
.ti +4
DOUBLE
PRECISION D( * )
.SH PURPOSE
PZPTTRF computes a Cholesky factorization of an N-by-N complex tridiagonal symmetric positive definite distributed matrix A(1:N, JA:JA+N-1). Reordering is used to increase parallelism in the factorization.
This reordering results in factors that are DIFFERENT from those
produced by equivalent sequential codes. These factors cannot
be used directly by users; however, they can be used in
.br
subsequent calls to PZPTTRS to solve linear systems.
.br

The factorization has the form
.br

        P A(1:N, JA:JA+N-1) P^T = U' D U  or
.br

        P A(1:N, JA:JA+N-1) P^T = L D L',
.br

where U is a tridiagonal upper triangular matrix and L is tridiagonal
lower triangular, and P is a permutation matrix.
.br

