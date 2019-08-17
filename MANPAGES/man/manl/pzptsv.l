.TH PZPTSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZPTSV - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 19
SUBROUTINE PZPTSV(
UPLO, N, NRHS, D, E, JA, DESCA, B, IB, DESCB,
WORK, LWORK, INFO )
.TP 19
.ti +4
CHARACTER
UPLO
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
COMPLEX*16
B( * ), E( * ), WORK( * )
.TP 19
.ti +4
DOUBLE
PRECISION D( * )
.SH PURPOSE
PZPTSV solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is an N-by-N complex
.br
tridiagonal symmetric positive definite distributed
.br
matrix.
.br

Cholesky factorization is used to factor a reordering of
.br
the matrix into L L'.
.br

See PZPTTRF and PZPTTRS for details.
.br

