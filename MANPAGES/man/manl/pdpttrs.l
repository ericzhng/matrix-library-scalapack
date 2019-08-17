.TH PDPTTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDPTTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PDPTTRS(
N, NRHS, D, E, JA, DESCA, B, IB, DESCB, AF,
LAF, WORK, LWORK, INFO )
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
PRECISION AF( * ), B( * ), D( * ), E( * ), WORK( * )
.SH PURPOSE
PDPTTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PDPTTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N real
.br
tridiagonal symmetric positive definite distributed
.br
matrix.
.br

Routine PDPTTRF MUST be called first.
.br

