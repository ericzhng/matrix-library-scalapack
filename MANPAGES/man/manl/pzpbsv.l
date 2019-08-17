.TH PZPBSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZPBSV - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 19
SUBROUTINE PZPBSV(
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
COMPLEX*16
A( * ), B( * ), WORK( * )
.SH PURPOSE
PZPBSV solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is an N-by-N complex
.br
banded symmetric positive definite distributed
.br
matrix with bandwidth BW.
.br

Cholesky factorization is used to factor a reordering of
.br
the matrix into L L'.
.br

See PZPBTRF and PZPBTRS for details.
.br

