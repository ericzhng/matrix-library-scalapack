.TH PDGBSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDGBSV - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 19
SUBROUTINE PDGBSV(
N, BWL, BWU, NRHS, A, JA, DESCA, IPIV, B, IB,
DESCB, WORK, LWORK, INFO )
.TP 19
.ti +4
INTEGER
BWL, BWU, IB, INFO, JA, LWORK, N, NRHS
.TP 19
.ti +4
INTEGER
DESCA( * ), DESCB( * ), IPIV( * )
.TP 19
.ti +4
DOUBLE
PRECISION A( * ), B( * ), WORK( * )
.SH PURPOSE
PDGBSV solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is an N-by-N real
.br
banded distributed
.br
matrix with bandwidth BWL, BWU.
.br

Gaussian elimination with pivoting
.br
is used to factor a reordering
.br
of the matrix into P L U.
.br

See PDGBTRF and PDGBTRS for details.
.br

