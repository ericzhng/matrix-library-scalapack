.TH PSLAED0 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAED0 - compute all eigenvalues and corresponding eigenvectors of a symmetric tridiagonal matrix using the divide and conquer method
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLAED0(
N, D, E, Q, IQ, JQ, DESCQ, WORK, IWORK, INFO )
.TP 20
.ti +4
INTEGER
INFO, IQ, JQ, N
.TP 20
.ti +4
INTEGER
DESCQ( * ), IWORK( * )
.TP 20
.ti +4
REAL
D( * ), E( * ), Q( * ), WORK( * )
.SH PURPOSE
PSLAED0 computes all eigenvalues and corresponding eigenvectors of a symmetric tridiagonal matrix using the divide and conquer method. 

.SH ARGUMENTS
.TP 8
N       (global input) INTEGER
The order of the tridiagonal matrix T.  N >= 0.
.TP 8
D       (global input/output) REAL array, dimension (N)
On entry, the diagonal elements of the tridiagonal matrix.
On exit, if INFO = 0, the eigenvalues in descending order.
.TP 8
E       (global input/output) REAL array, dimension (N-1)
On entry, the subdiagonal elements of the tridiagonal matrix.
On exit, E has been destroyed.
.TP 8
Q       (local output) REAL array,
global dimension (N, N),
local dimension ( LLD_Q, LOCc(JQ+N-1))
Q  contains the orthonormal eigenvectors of the symmetric
tridiagonal matrix.
On output, Q is distributed across the P processes in block
cyclic format.
.TP 8
IQ      (global input) INTEGER
Q's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JQ      (global input) INTEGER
Q's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCQ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Z.
.TP 8
WORK    (local workspace ) REAL array, dimension (LWORK)
LWORK = 6*N + 2*NP*NQ, with
NP = NUMROC( N, MB_Q, MYROW, IQROW, NPROW )
NQ = NUMROC( N, NB_Q, MYCOL, IQCOL, NPCOL )
IQROW = INDXG2P( IQ, NB_Q, MYROW, RSRC_Q, NPROW )
IQCOL = INDXG2P( JQ, MB_Q, MYCOL, CSRC_Q, NPCOL )
.TP 8
IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
LIWORK = 2 + 7*N + 8*NPCOL
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  The algorithm failed to compute the INFO/(N+1) th
eigenvalue while working on the submatrix lying in
global rows and columns mod(INFO,N+1).
