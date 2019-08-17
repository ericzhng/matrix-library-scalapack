.TH PCDTSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCDTSV - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 19
SUBROUTINE PCDTSV(
N, NRHS, DL, D, DU, JA, DESCA, B, IB, DESCB,
WORK, LWORK, INFO )
.TP 19
.ti +4
INTEGER
IB, INFO, JA, LWORK, N, NRHS
.TP 19
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 19
.ti +4
COMPLEX
B( * ), D( * ), DL( * ), DU( * ), WORK( * )
.SH PURPOSE
PCDTSV solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is an N-by-N complex
.br
tridiagonal diagonally dominant-like distributed
.br
matrix.
.br

Gaussian elimination without pivoting
.br
is used to factor a reordering
.br
of the matrix into L U.
.br

See PCDTTRF and PCDTTRS for details.
.br

