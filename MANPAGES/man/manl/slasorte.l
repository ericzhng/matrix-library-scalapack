.TH SLASORTE l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
SLASORTE - sort eigenpairs so that real eigenpairs are together and complex are together
.SH SYNOPSIS
.TP 21
SUBROUTINE SLASORTE(
S, LDS, J, OUT, INFO )
.TP 21
.ti +4
INTEGER
INFO, J, LDS
.TP 21
.ti +4
REAL
OUT( J, * ), S( LDS, * )
.SH PURPOSE
SLASORTE sorts eigenpairs so that real eigenpairs are together and complex are together. This way one can employ 2x2 shifts easily   since every 2nd subdiagonal is guaranteed to be zero.
.br
This routine does no parallel work.
.br

.SH ARGUMENTS
.TP 8
S       (local input/output) REAL array, dimension LDS
On entry, a matrix already in Schur form.
On exit, the diagonal blocks of S have been rewritten to pair
the eigenvalues.  The resulting matrix is no longer
similar to the input.
.TP 8
LDS     (local input) INTEGER
On entry, the leading dimension of the local array S.
Unchanged on exit.
.TP 8
J       (local input) INTEGER
On entry, the order of the matrix S.
Unchanged on exit.
.TP 8
OUT     (local input/output) REAL array, dimension Jx2
This is the work buffer required by this routine.
.TP 8
INFO    (local input) INTEGER
This is set if the input matrix had an odd number of real
eigenvalues and things couldn't be paired or if the input
matrix S was not originally in Schur form.
0 indicates successful completion.

Implemented by:  G. Henry, November 17, 1996
