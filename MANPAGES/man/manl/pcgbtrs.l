.TH PCGBTRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCGBTRS - solve a system of linear equations  A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)
.SH SYNOPSIS
.TP 20
SUBROUTINE PCGBTRS(
TRANS, N, BWL, BWU, NRHS, A, JA, DESCA, IPIV,
B, IB, DESCB, AF, LAF, WORK, LWORK, INFO )
.TP 20
.ti +4
CHARACTER
TRANS
.TP 20
.ti +4
INTEGER
BWU, BWL, IB, INFO, JA, LAF, LWORK, N, NRHS
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCB( * ), IPIV(*)
.TP 20
.ti +4
COMPLEX
A( * ), AF( * ), B( * ), WORK( * )
.SH PURPOSE
PCGBTRS solves a system of linear equations A(1:N, JA:JA+N-1) * X = B(IB:IB+N-1, 1:NRHS)                                   or
.br
          A(1:N, JA:JA+N-1)' * X = B(IB:IB+N-1, 1:NRHS)
.br

where A(1:N, JA:JA+N-1) is the matrix used to produce the factors
stored in A(1:N,JA:JA+N-1) and AF by PCGBTRF.
.br
A(1:N, JA:JA+N-1) is an N-by-N complex
.br
banded distributed
.br
matrix with bandwidth BWL, BWU.
.br

Routine PCGBTRF MUST be called first.
.br

