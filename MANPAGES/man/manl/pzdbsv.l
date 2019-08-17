.TH PZDBSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZDBSV - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 19
SUBROUTINE PZDBSV(
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
COMPLEX*16
A( * ), B( * ), WORK( * )
.SH PURPOSE
PZDBSV solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is an N-by-N complex
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

See PZDBTRF and PZDBTRS for details.
.br

