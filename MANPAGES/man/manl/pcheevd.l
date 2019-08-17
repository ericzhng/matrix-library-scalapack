.TH PCHEEVD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCHEEVD - compute all the eigenvalues and eigenvectors of a Hermitian matrix A by using a divide and conquer algorithm
.SH SYNOPSIS
.TP 20
SUBROUTINE PCHEEVD(
JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
DESCZ, WORK, LWORK, RWORK, LRWORK, IWORK,
LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
JOBZ, UPLO
.TP 20
.ti +4
INTEGER
IA, INFO, IZ, JA, JZ, LIWORK, LRWORK, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCZ( * ), IWORK( * )
.TP 20
.ti +4
REAL
RWORK( * ), W( * )
.TP 20
.ti +4
COMPLEX
A( * ), WORK( * ), Z( * )
.SH PURPOSE
PCHEEVD computes all the eigenvalues and eigenvectors of a Hermitian matrix A by using a divide and conquer algorithm. 
.SH ARGUMENTS

NP = the number of rows local to a given process.
NQ = the number of columns local to a given process.
.TP 8
JOBZ    (input) CHARACTER*1
= 'N':  Compute eigenvalues only;    (NOT IMPLEMENTED YET)
.br
= 'V':  Compute eigenvalues and eigenvectors.
.TP 8
UPLO    (global input) CHARACTER*1
Specifies whether the upper or lower triangular part of the
symmetric matrix A is stored:
.br
= 'U':  Upper triangular
.br
= 'L':  Lower triangular
.TP 8
N       (global input) INTEGER
The number of rows and columns of the matrix A.  N >= 0.
.TP 8
A       (local input/workspace) block cyclic COMPLEX array,
global dimension (N, N), local dimension ( LLD_A,
LOCc(JA+N-1) )

On entry, the symmetric matrix A.  If UPLO = 'U', only the
upper triangular part of A is used to define the elements of
the symmetric matrix.  If UPLO = 'L', only the lower
triangular part of A is used to define the elements of the
symmetric matrix.

On exit, the lower triangle (if UPLO='L') or the upper
triangle (if UPLO='U') of A, including the diagonal, is
destroyed.
.TP 8
IA      (global input) INTEGER
A's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JA      (global input) INTEGER
A's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
If DESCA( CTXT_ ) is incorrect, PCHEEVD cannot guarantee
correct error reporting.
.TP 8
W       (global output) REAL array, dimension (N)
If INFO=0, the eigenvalues in ascending order.
.TP 8
Z       (local output) COMPLEX array,
global dimension (N, N),
local dimension ( LLD_Z, LOCc(JZ+N-1) )
Z contains the orthonormal eigenvectors of the matrix A.
.TP 8
IZ      (global input) INTEGER
Z's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JZ      (global input) INTEGER
Z's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCZ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Z.
DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
.TP 8
WORK    (local workspace/output) COMPLEX array,
dimension (LWORK)
On output, WORK(1) returns the workspace needed for the
computation.
.TP 8
LWORK   (local input) INTEGER
If eigenvectors are requested:
LWORK = N + ( NP0 + MQ0 + NB ) * NB,
with  NP0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPROW )
MQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine calculates the size for all
work arrays. Each of these values is returned in the first
entry of the corresponding work array, and no error message
is issued by PXERBLA.
.TP 8
RWORK   (local workspace/output) REAL array,
dimension (LRWORK)
On output RWORK(1) returns the real workspace needed to
guarantee completion.  If the input parameters are incorrect,
RWORK(1) may also be incorrect.
.TP 8
LRWORK  (local input) INTEGER
Size of RWORK array.
RWORK >= 1 + 8*N + 2*NP*NQ,
NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )
.TP 8
IWORK   (local workspace/output) INTEGER array, dimension (LIWORK)
On output IWORK(1) returns the integer workspace needed.
.TP 8
LIWORK  (input) INTEGER
The dimension of the array IWORK.
LIWORK = 7*N + 8*NPCOL + 2
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  If INFO = 1 through N, the i(th) eigenvalue did not
converge in PSLAED3.

Alignment requirements
======================

The distributed submatrices sub( A ), sub( Z ) must verify
some alignment properties, namely the following expression
should be true:
( MB_A.EQ.NB_A.EQ.MB_Z.EQ.NB_Z .AND. IROFFA.EQ.ICOFFA .AND.
IROFFA.EQ.0 .AND.IROFFA.EQ.IROFFZ. AND. IAROW.EQ.IZROW)
with IROFFA = MOD( IA-1, MB_A )
and ICOFFA = MOD( JA-1, NB_A ).
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

