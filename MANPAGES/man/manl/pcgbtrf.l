.TH PCGBTRF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCGBTRF - compute a LU factorization of an N-by-N complex banded distributed matrix with bandwidth BWL, BWU
.SH SYNOPSIS
.TP 20
SUBROUTINE PCGBTRF(
N, BWL, BWU, A, JA, DESCA, IPIV, AF, LAF,
WORK, LWORK, INFO )
.TP 20
.ti +4
INTEGER
BWL, BWU, INFO, JA, LAF, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * ), IPIV( * )
.TP 20
.ti +4
COMPLEX
A( * ), AF( * ), WORK( * )
.SH PURPOSE
PCGBTRF computes a LU factorization of an N-by-N complex banded distributed matrix with bandwidth BWL, BWU: A(1:N, JA:JA+N-1). Reordering is used to increase parallelism in the factorization.
This reordering results in factors that are DIFFERENT from those
produced by equivalent sequential codes. These factors cannot
be used directly by users; however, they can be used in
.br
subsequent calls to PCGBTRS to solve linear systems.
.br

The factorization has the form
.br

        P A(1:N, JA:JA+N-1) Q = L U
.br

where U is a banded upper triangular matrix and L is banded
lower triangular, and P and Q are permutation matrices.
.br
The matrix Q represents reordering of columns
.br
for parallelism's sake, while P represents
.br
reordering of rows for numerical stability using
.br
classic partial pivoting.
.br

