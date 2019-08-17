.TH PSLAEDZ l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAEDZ - Form the z-vector which consists of the last row of Q_1 and the first row of Q_2
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLAEDZ(
N, N1, ID, Q, IQ, JQ, LDQ, DESCQ, Z, WORK )
.TP 20
.ti +4
INTEGER
ID, IQ, JQ, LDQ, N, N1
.TP 20
.ti +4
INTEGER
DESCQ( * )
.TP 20
.ti +4
REAL
Q( LDQ, * ), WORK( * ), Z( * )
.SH PURPOSE
PSLAEDZ Form the z-vector which consists of the last row of Q_1 and the first row of Q_2. 