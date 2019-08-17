.TH PDDBTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDDBTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PDDBTRS(
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
DOUBLE
PRECISION A( * ), AF( * ), B( * ), WORK( * )
.SH PURPOSE
PDDBTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)                                   or
.br
          A(1:N, JA:JA+N-1)' * X = B(IB:IB+N-1, 1:NRHS)
.br

where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PDDBTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N real
.br
banded diagonally dominant-like distributed
.br
matrix with bandwidth BWL, BWU.
.br

Routine PDDBTRF MUST be called first.
.br

