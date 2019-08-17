.TH PCDBTRF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCDBTRF - compute a LU factorization of an N-by-N complex banded diagonally dominant-like distributed matrix with bandwidth BWL, BWU
.SH SYNOPSIS
.TP 20
SUBROUTINE PCDBTRF(
N, BWL, BWU, A, JA, DESCA, AF, LAF, WORK,
LWORK, INFO )
.TP 20
.ti +4
INTEGER
BWL, BWU, INFO, JA, LAF, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
COMPLEX
A( * ), AF( * ), WORK( * )
.SH PURPOSE
PCDBTRF computes a LU factorization of an N-by-N complex banded diagonally dominant-like distributed matrix with bandwidth BWL, BWU: A(1:N, JA:JA+N-1). Reordering is used to increase parallelism in the factorization.
This reordering results in factors that are DIFFERENT from those
produced by equivalent sequential codes. These factors cannot
be used directly by users; however, they can be used in
.br
subsequent calls to PCDBTRS to solve linear systems.
.br

The factorization has the form
.br

        P A(1:N, JA:JA+N-1) P^T = L U
.br

where U is a banded upper triangular matrix and L is banded
lower triangular, and P is a permutation matrix.
.br

