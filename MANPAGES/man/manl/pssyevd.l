.TH PSSYEVD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSSYEVD - compute all the eigenvalues and eigenvectors of a real symmetric matrix A by calling the recommended sequence of ScaLAPACK routines
.SH SYNOPSIS
.TP 20
SUBROUTINE PSSYEVD(
JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
DESCZ, WORK, LWORK, IWORK, LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
JOBZ, UPLO
.TP 20
.ti +4
INTEGER
IA, INFO, IZ, JA, JZ, LIWORK, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCZ( * ), IWORK( * )
.TP 20
.ti +4
REAL
A( * ), W( * ), WORK( * ), Z( * )
.SH PURPOSE
PSSYEVD computes all the eigenvalues and eigenvectors of a real symmetric matrix A by calling the recommended sequence of ScaLAPACK routines. 
In its present form, PSSYEVD assumes a homogeneous system and makes
no checks for consistency of the eigenvalues or eigenvectors across
the different processes.  Because of this, it is possible that a
heterogeneous system may return incorrect results without any error
messages.
.br

.SH ARGUMENTS

NP = the number of rows local to a given process.
NQ = the number of columns local to a given process.
.TP 8
JOBZ    (input) CHARACTER*1
= 'N':  Compute eigenvalues only;     (NOT IMPLEMENTED YET)
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
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/workspace) block cyclic REAL array,
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
.TP 8
W       (global output) REAL array, dimension (N)
If INFO=0, the eigenvalues in ascending order.
.TP 8
Z       (local output) REAL array,
global dimension (N, N),
local dimension ( LLD_Z, LOCc(JZ+N-1) )
Z contains the orthonormal eigenvectors
of the symmetric matrix A.
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
WORK    (local workspace/output) REAL array,
dimension (LWORK)
On output, WORK(1) returns the workspace required.
.TP 8
LWORK   (local input) INTEGER
LWORK >= MAX( 1+6*N+2*NP*NQ, TRILWMIN ) + 2*N
TRILWMIN = 3*N + MAX( NB*( NP+1 ), 3*NB )
NP = NUMROC( N, NB, MYROW, IAROW, NPROW )
NQ = NUMROC( N, NB, MYCOL, IACOL, NPCOL )

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
LIWORK = 7*N + 8*NPCOL + 2
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

