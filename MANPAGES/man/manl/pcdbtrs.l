.TH PCDBTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCDBTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PCDBTRS(
TRANS, N, BWL, BWU, NRHS, A, JA, DESCA, B, IB,
DESCB, AF, LAF, WORK, LWORK, INFO )
.TP 20
.ti +4
CHARACTER
TRANS
.TP 20
.ti +4
INTEGER
BWL, BWU, IB, INFO, JA, LAF, LWORK, N, NRHS
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 20
.ti +4
COMPLEX
A( * ), AF( * ), B( * ), WORK( * )
.SH PURPOSE
PCDBTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)                                   or
.br
          A(1:N, JA:JA+N-1)' * X = B(IB:IB+N-1, 1:NRHS)
.br

where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PCDBTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N complex
.br
banded diagonally dominant-like distributed
.br
matrix with bandwidth BWL, BWU.
.br

Routine PCDBTRF MUST be called first.
.br

