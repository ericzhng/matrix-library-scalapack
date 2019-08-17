.TH PSLATRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLATRS - solve a triangular system
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLATRS(
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
A( * ), CNORM( * ),
X( * ), WORK( * )
.SH PURPOSE
PSLATRS solves a triangular system. This routine in unfinished at this time, but will be part of the next release.
.br

