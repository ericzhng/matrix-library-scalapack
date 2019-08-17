.TH PCLATRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCLATRS - solve a triangular system
.SH SYNOPSIS
.TP 20
SUBROUTINE PCLATRS(
UPLO, TRANS, DIAG, NORMIN, N, A, IA,
JA, DESCA, X, IX, JX, DESCX, SCALE, CNORM,
WORK )
.TP 20
.ti +4
CHARACTER
DIAG, NORMIN, TRANS, UPLO
.TP 20
.ti +4
INTEGER
IA, IX, JA, JX, N
.TP 20
.ti +4
REAL
SCALE
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCX( * )
.TP 20
.ti +4
REAL
CNORM( * )
.TP 20
.ti +4
COMPLEX
A( * ), X( * ), WORK( * )
.SH PURPOSE
PCLATRS solves a triangular system. This routine in unfinished at this time, but will be part of the next release.
.br

