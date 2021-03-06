.TH PCPTTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCPTTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PCPTTRS(
UPLO, N, NRHS, D, E, JA, DESCA, B, IB, DESCB,
AF, LAF, WORK, LWORK, INFO )
.TP 20
.ti +4
CHARACTER
UPLO
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
COMPLEX
AF( * ), B( * ), E( * ), WORK( * )
.TP 20
.ti +4
REAL
D( * )
.SH PURPOSE
PCPTTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS) 
where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PCPTTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N complex
.br
tridiagonal symmetric positive definite distributed
.br
matrix.
.br
Depending on the value of UPLO, A stores either U or L in the equn
A(1:N, JA:JA+N-1) = U'D *U or L*D L' as computed by PCPTTRF.

Routine PCPTTRF MUST be called first.
.br

