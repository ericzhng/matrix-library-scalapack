.TH PDDBSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDDBSV - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 19
SUBROUTINE PDDBSV(
N, BWL, BWU, NRHS, A, JA, DESCA, B, IB, DESCB,
WORK, LWORK, INFO )
.TP 19
.ti +4
INTEGER
BWL, BWU, IB, INFO, JA, LWORK, N, NRHS
.TP 19
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 19
.ti +4
DOUBLE
PRECISION A( * ), B( * ), WORK( * )
.SH PURPOSE
PDDBSV solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is an N-by-N real
.br
banded diagonally dominant-like distributed
.br
matrix with bandwidth BWL, BWU.
.br

Gaussian elimination without pivoting
.br
is used to factor a reordering
.br
of the matrix into L U.
.br

See PDDBTRF and PDDBTRS for details.
.br

