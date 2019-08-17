.TH PZPBTRF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZPBTRF - compute a Cholesky factorization of an N-by-N complex banded symmetric positive definite distributed matrix with bandwidth BW
.SH SYNOPSIS
.TP 20
SUBROUTINE PZPBTRF(
UPLO, N, BW, A, JA, DESCA, AF, LAF, WORK,
LWORK, INFO )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
BW, INFO, JA, LAF, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
COMPLEX*16
A( * ), AF( * ), WORK( * )
.SH PURPOSE
PZPBTRF computes a Cholesky factorization of an N-by-N complex banded symmetric positive definite distributed matrix with bandwidth BW: A(1:N, JA:JA+N-1). Reordering is used to increase parallelism in the factorization.
This reordering results in factors that are DIFFERENT from those
produced by equivalent sequential codes. These factors cannot
be used directly by users; however, they can be used in
.br
subsequent calls to PZPBTRS to solve linear systems.
.br

The factorization has the form
.br

        P A(1:N, JA:JA+N-1) P^T = U' U ,  if UPLO = 'U', or

        P A(1:N, JA:JA+N-1) P^T = L L', if UPLO = 'L'
.br

where U is a banded upper triangular matrix and L is banded
lower triangular, and P is a permutation matrix.
.br

