.TH PZPBTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZPBTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PZPBTRS(
UPLO, N, BW, NRHS, A, JA, DESCA, B, IB, DESCB,
AF, LAF, WORK, LWORK, INFO )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
BW, IB, INFO, JA, LAF, LWORK, N, NRHS
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 20
.ti +4
COMPLEX*16
A( * ), AF( * ), B( * ), WORK( * )
.SH PURPOSE
PZPBTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PZPBTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N complex
.br
banded symmetric positive definite distributed
.br
matrix with bandwidth BW.
.br
Depending on the value of UPLO, A stores either U or L in the equn
A(1:N, JA:JA+N-1) = U'*U or L*L' as computed by PZPBTRF.
.br

Routine PZPBTRF MUST be called first.
.br

