.TH PZLATRS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZLATRS - solve a triangular system
.SH SYNOPSIS
.TP 20
SUBROUTINE PZLATRS(
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
DOUBLE
PRECISION SCALE
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCX( * )
.TP 20
.ti +4
DOUBLE
PRECISION CNORM( * )
.TP 20
.ti +4
COMPLEX*16
A( * ), X( * ), WORK( * )
.SH PURPOSE
PZLATRS solves a triangular system. This routine in unfinished at this time, but will be part of the next release.
.br

