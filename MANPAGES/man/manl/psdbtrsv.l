.TH PSDBTRSV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSDBTRSV - solve a banded triangular system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 21
SUBROUTINE PSDBTRSV(
UPLO, TRANS, N, BWL, BWU, NRHS, A, JA, DESCA,
B, IB, DESCB, AF, LAF, WORK, LWORK, INFO )
.TP 21
.ti +4
CHARACTER
TRANS, UPLO
.TP 21
.ti +4
INTEGER
BWL, BWU, IB, INFO, JA, LAF, LWORK, N, NRHS
.TP 21
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 21
.ti +4
REAL
A( * ), AF( * ), B( * ), WORK( * )
.SH PURPOSE
PSDBTRSV solves a banded triangular system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)                                or
.br
                    A(1:N, JA:JA+N-1)^T * X = B(IB:IB+N-1, 1:NRHS)

where A(1:N, JA:JA+N-1) is a banded
.br
triangular matrix factor produced by the
.br
Gaussian elimination code PD@(dom_pre)BTRF
.br
and is stored in A(1:N,JA:JA+N-1) and AF.
.br
The matrix stored in A(1:N, JA:JA+N-1) is either
.br
upper or lower triangular according to UPLO,
.br
and the choice of solving A(1:N, JA:JA+N-1) or A(1:N, JA:JA+N-1)^T
is dictated by the user by the parameter TRANS.
.br

Routine PSDBTRF MUST be called first.
.br

