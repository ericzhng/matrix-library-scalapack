.TH SLAREF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
SLAREF - applie one or several Householder reflectors of size 3 to one or two matrices (if column is specified) on either their rows or columns
.SH SYNOPSIS
.TP 19
SUBROUTINE SLAREF(
TYPE, A, LDA, WANTZ, Z, LDZ, BLOCK, IROW1,
ICOL1, ISTART, ISTOP, ITMP1, ITMP2, LILOZ,
LIHIZ, VECS, V2, V3, T1, T2, T3 )
.TP 19
.ti +4
LOGICAL
BLOCK, WANTZ
.TP 19
.ti +4
CHARACTER
TYPE
.TP 19
.ti +4
INTEGER
ICOL1, IROW1, ISTART, ISTOP, ITMP1, ITMP2, LDA,
LDZ, LIHIZ, LILOZ
.TP 19
.ti +4
REAL
T1, T2, T3, V2, V3
.TP 19
.ti +4
REAL
A( LDA, * ), VECS( * ), Z( LDZ, * )
.SH PURPOSE
SLAREF applies one or several Householder reflectors of size 3 to one or two matrices (if column is specified) on either their rows or columns. 
.SH ARGUMENTS
.TP 8
TYPE    (global input) CHARACTER*1
If 'R': Apply reflectors to the rows of the matrix
(apply from left)
Otherwise: Apply reflectors to the columns of the matrix
Unchanged on exit.
.TP 8
A       (global input/output) REAL array, (LDA,*)
On entry, the matrix to receive the reflections.
The updated matrix on exit.
.TP 8
LDA     (local input) INTEGER
On entry, the leading dimension of A.  Unchanged on exit.
.TP 8
WANTZ   (global input) LOGICAL
If .TRUE., then apply any column reflections to Z as well.
If .FALSE., then do no additional work on Z.
.TP 8
Z       (global input/output) REAL array, (LDZ,*)
On entry, the second matrix to receive column reflections.
This is changed only if WANTZ is set.
.TP 8
LDZ     (local input) INTEGER
On entry, the leading dimension of Z.  Unchanged on exit.
.TP 8
BLOCK   (global input) LOGICAL
If .TRUE., then apply several reflectors at once and read
their data from the VECS array.
If .FALSE., apply the single reflector given by V2, V3,
T1, T2, and T3.
.TP 8
IROW1   (local input/output) INTEGER
On entry, the local row element of A.
Undefined on output.
.TP 8
ICOL1   (local input/output) INTEGER
On entry, the local column element of A.
Undefined on output.
.TP 8
ISTART  (global input) INTEGER
Specifies the "number" of the first reflector.  This is
used as an index into VECS if BLOCK is set.
ISTART is ignored if BLOCK is .FALSE..
.TP 8
ISTOP   (global input) INTEGER
Specifies the "number" of the last reflector.  This is
used as an index into VECS if BLOCK is set.
ISTOP is ignored if BLOCK is .FALSE..
.TP 8
ITMP1   (local input) INTEGER
Starting range into A.  For rows, this is the local
first column.  For columns, this is the local first row.
.TP 8
ITMP2   (local input) INTEGER
Ending range into A.  For rows, this is the local last
column.  For columns, this is the local last row.

LILOZ
LIHIZ   (local input) INTEGER
These serve the same purpose as ITMP1,ITMP2 but for Z
when WANTZ is set.
.TP 8
VECS    (global input) REAL array of size 3*N (matrix
size)
This holds the size 3 reflectors one after another and this
is only accessed when BLOCK is .TRUE.

V2
V3
T1
T2
T3      (global input/output) REAL
This holds information on a single size 3 Householder
reflector and is read when BLOCK is .FALSE., and
overwritten when BLOCK is .TRUE.

Implemented by:  G. Henry, November 17, 1996
