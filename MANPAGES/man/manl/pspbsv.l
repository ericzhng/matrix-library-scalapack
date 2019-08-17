.TH PSPBSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSPBSV - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 19
SUBROUTINE PSPBSV(
UPLO, N, BW, NRHS, A, JA, DESCA, B, IB, DESCB,
WORK, LWORK, INFO )
.TP 19
.ti +4
CHARACTER
UPLO
.TP 19
.ti +4
INTEGER
BW, IB, INFO, JA, LWORK, N, NRHS
.TP 19
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 19
.ti +4
REAL
A( * ), B( * ), WORK( * )
.SH PURPOSE
PSPBSV solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is an N-by-N real
.br
banded symmetric positive definite distributed
.br
matrix with bandwidth BW.
.br

Cholesky factorization is used to factor a reordering of
.br
the matrix into L L'.
.br

See PSPBTRF and PSPBTRS for details.
.br

