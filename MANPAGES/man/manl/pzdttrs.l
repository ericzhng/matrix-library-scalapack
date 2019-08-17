.TH PZDTTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZDTTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PZDTTRS(
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
COMPLEX*16
AF( * ), B( * ), D( * ), DL( * ), DU( * ),
WORK( * )
.SH PURPOSE
PZDTTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)                                   or
.br
          A(1:N, JA:JA+N-1)' * X = B(IB:IB+N-1, 1:NRHS)
.br

where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PZDTTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N complex
.br
tridiagonal diagonally dominant-like distributed
.br
matrix.
.br

Routine PZDTTRF MUST be called first.
.br

