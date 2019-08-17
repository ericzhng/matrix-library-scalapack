.TH PDDTTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDDTTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PDDTTRS(
TRANS, N, NRHS, DL, D, DU, JA, DESCA, B, IB,
DESCB, AF, LAF, WORK, LWORK, INFO )
.TP 20
.ti +4
CHARACTER
TRANS
.TP 20
.ti +4
INTEGER
IB, INFO, JA, LAF, LWORK, N, NRHS
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 20
.ti +4
DOUBLE
PRECISION AF( * ), B( * ), D( * ), DL( * ), DU( * ),
WORK( * )
.SH PURPOSE
PDDTTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)                                   or
.br
          A(1:N, JA:JA+N-1)' * X = B(IB:IB+N-1, 1:NRHS)
.br

where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PDDTTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N real
.br
tridiagonal diagonally dominant-like distributed
.br
matrix.
.br

Routine PDDTTRF MUST be called first.
.br

