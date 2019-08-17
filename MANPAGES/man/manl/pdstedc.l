.TH PDSTEDC l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDSTEDC - tridiagonal matrix in parallel, using the divide and conquer algorithm
.SH SYNOPSIS
.TP 20
SUBROUTINE PDSTEDC(
COMPZ, N, D, E, Q, IQ, JQ, DESCQ, WORK, LWORK,
IWORK, LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
COMPZ
.TP 20
.ti +4
INTEGER
INFO, IQ, JQ, LIWORK, LWORK, N
.TP 20
.ti +4
INTEGER
DESCQ( * ), IWORK( * )
.TP 20
.ti +4
DOUBLE
PRECISION D( * ), E( * ), Q( * ), WORK( * )
.SH PURPOSE
symmetric tridiagonal matrix in parallel, using the divide and conquer algorithm. 
This code makes very mild assumptions about floating point
arithmetic. It will work on machines with a guard digit in
add/subtract, or on those binary machines without guard digits
which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
It could conceivably fail on hexadecimal or decimal machines
without guard digits, but we know of none.  See DLAED3 for details.

.SH ARGUMENTS
.TP 8
COMPZ   (input) CHARACTER*1
= 'N':  Compute eigenvalues only.    (NOT IMPLEMENTED YET)
.br
= 'I':  Compute eigenvectors of tridiagonal matrix also.
.br
= 'V':  Compute eigenvectors of original dense symmetric
matrix also.  On entry, Z contains the orthogonal
matrix used to reduce the original matrix to
tridiagonal form.            (NOT IMPLEMENTED YET)
.TP 8
N       (global input) INTEGER
The order of the tridiagonal matrix T.  N >= 0.
.TP 8
D       (global input/output) DOUBLE PRECISION array, dimension (N)
On entry, the diagonal elements of the tridiagonal matrix.
On exit, if INFO = 0, the eigenvalues in descending order.
.TP 8
E       (global input/output) DOUBLE PRECISION array, dimension (N-1)
On entry, the subdiagonal elements of the tridiagonal matrix.
On exit, E has been destroyed.
.TP 8
Q       (local output) DOUBLE PRECISION array,
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
WORK    (local workspace/output) DOUBLE PRECISION array,
dimension (LWORK)
On output, WORK(1) returns the workspace needed.
.TP 8
LWORK   (local input/output) INTEGER,
the dimension of the array WORK.
LWORK = 6*N + 2*NP*NQ
NP = NUMROC( N, NB, MYROW, DESCQ( RSRC_ ), NPROW )
NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )

If LWORK = -1, the LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
size for the WORK array.  The required workspace is returned
as the first element of WORK and no error message is issued
by PXERBLA.
.TP 8
IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
.TP 8
LIWORK  (input) INTEGER
The dimension of the array IWORK.
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
.SH FURTHER DETAILS
Contributed by Francoise Tisseur, University of Manchester.

Reference:  F. Tisseur and J. Dongarra, "A Parallel Divide and
            Conquer Algorithm for the Symmetric Eigenvalue Problem
            on Distributed Memory Architectures",
.br
            SIAM J. Sci. Comput., 6:20 (1999), pp. 2223--2236.
            (see also LAPACK Working Note 132)
.br
              http://www.netlib.org/lapack/lawns/lawn132.ps

