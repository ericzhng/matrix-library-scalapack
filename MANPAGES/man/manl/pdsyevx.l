.TH PDSYEVX l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDSYEVX - compute selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A by calling the recommended sequence of ScaLAPACK routines
.SH SYNOPSIS
.TP 20
SUBROUTINE PDSYEVX(
JOBZ, RANGE, UPLO, N, A, IA, JA, DESCA, VL,
VU, IL, IU, ABSTOL, M, NZ, W, ORFAC, Z, IZ,
JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL,
ICLUSTR, GAP, INFO )
.TP 20
.ti +4
CHARACTER
JOBZ, RANGE, UPLO
.TP 20
.ti +4
INTEGER
IA, IL, INFO, IU, IZ, JA, JZ, LIWORK, LWORK, M,
N, NZ
.TP 20
.ti +4
DOUBLE
PRECISION ABSTOL, ORFAC, VL, VU
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCZ( * ), ICLUSTR( * ),
IFAIL( * ), IWORK( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), GAP( * ), W( * ), WORK( * ), Z( * )
.SH PURPOSE
PDSYEVX computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A by calling the recommended sequence of ScaLAPACK routines. Eigenvalues/vectors can be selected by specifying a range of values or a range of indices for the desired
eigenvalues.
.br

Notes
.br
=====
.br

Each global data object is described by an associated description
vector.  This vector stores the information required to establish
the mapping between an object element and its corresponding process
and memory location.
.br

Let A be a generic term for any 2D block cyclicly distributed array.
Such a global array has an associated description vector DESCA.
In the following comments, the character _ should be read as
"of the global array".
.br

NOTATION        STORED IN      EXPLANATION
.br
--------------- -------------- --------------------------------------
DTYPE_A(global) DESCA( DTYPE_ )The descriptor type.  In this case,
                               DTYPE_A = 1.
.br
CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
                               the BLACS process grid A is distribu-
                               ted over. The context itself is glo-
                               bal, but the handle (the integer
                               value) may vary.
.br
M_A    (global) DESCA( M_ )    The number of rows in the global
                               array A.
.br
N_A    (global) DESCA( N_ )    The number of columns in the global
                               array A.
.br
MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
                               the rows of the array.
.br
NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
                               the columns of the array.
.br
RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
                               row of the array A is distributed.
CSRC_A (global) DESCA( CSRC_ ) The process column over which the
                               first column of the array A is
                               distributed.
.br
LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
                               array.  LLD_A >= MAX(1,LOCr(M_A)).

Let K be the number of rows or columns of a distributed matrix,
and assume that its process grid has dimension p x q.
.br
LOCr( K ) denotes the number of elements of K that a process
would receive if K were distributed over the p processes of its
process column.
.br
Similarly, LOCc( K ) denotes the number of elements of K that a
process would receive if K were distributed over the q processes of
its process row.
.br
The values of LOCr() and LOCc() may be determined via a call to the
ScaLAPACK tool function, NUMROC:
.br
        LOCr( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
        LOCc( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
An upper bound for these quantities may be computed by:
.br
        LOCr( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
.br
        LOCc( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
.br

PDSYEVX assumes IEEE 754 standard compliant arithmetic.  To port
to a system which does not have IEEE 754 arithmetic, modify
the appropriate SLmake.inc file to include the compiler switch
-DNO_IEEE.  This switch only affects the compilation of pdlaiect.c.

.SH ARGUMENTS

NP = the number of rows local to a given process.
NQ = the number of columns local to a given process.
.TP 8
JOBZ    (global input) CHARACTER*1
Specifies whether or not to compute the eigenvectors:
.br
= 'N':  Compute eigenvalues only.
.br
= 'V':  Compute eigenvalues and eigenvectors.
.TP 8
RANGE   (global input) CHARACTER*1
.br
= 'A': all eigenvalues will be found.
.br
= 'V': all eigenvalues in the interval [VL,VU] will be found.
.br
= 'I': the IL-th through IU-th eigenvalues will be found.
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
A       (local input/workspace) block cyclic DOUBLE PRECISION array,
global dimension (N, N),
local dimension ( LLD_A, LOCc(JA+N-1) )

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
If DESCA( CTXT_ ) is incorrect, PDSYEVX cannot guarantee
correct error reporting.
.TP 8
VL      (global input) DOUBLE PRECISION
If RANGE='V', the lower bound of the interval to be searched
for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
.TP 8
VU      (global input) DOUBLE PRECISION
If RANGE='V', the upper bound of the interval to be searched
for eigenvalues.  Not referenced if RANGE = 'A' or 'I'.
.TP 8
IL      (global input) INTEGER
If RANGE='I', the index (from smallest to largest) of the
smallest eigenvalue to be returned.  IL >= 1.
Not referenced if RANGE = 'A' or 'V'.
.TP 8
IU      (global input) INTEGER
If RANGE='I', the index (from smallest to largest) of the
largest eigenvalue to be returned.  min(IL,N) <= IU <= N.
Not referenced if RANGE = 'A' or 'V'.
.TP 8
ABSTOL  (global input) DOUBLE PRECISION
If JOBZ='V', setting ABSTOL to PDLAMCH( CONTEXT, 'U') yields
the most orthogonal eigenvectors.

The absolute error tolerance for the eigenvalues.
An approximate eigenvalue is accepted as converged
when it is determined to lie in an interval [a,b]
of width less than or equal to

ABSTOL + EPS *   max( |a|,|b| ) ,

where EPS is the machine precision.  If ABSTOL is less than
or equal to zero, then EPS*norm(T) will be used in its place,
where norm(T) is the 1-norm of the tridiagonal matrix
obtained by reducing A to tridiagonal form.

Eigenvalues will be computed most accurately when ABSTOL is
set to twice the underflow threshold 2*PDLAMCH('S') not zero.
If this routine returns with ((MOD(INFO,2).NE.0) .OR.
(MOD(INFO/8,2).NE.0)), indicating that some eigenvalues or
eigenvectors did not converge, try setting ABSTOL to
2*PDLAMCH('S').

See "Computing Small Singular Values of Bidiagonal Matrices
with Guaranteed High Relative Accuracy," by Demmel and
Kahan, LAPACK Working Note #3.

See "On the correctness of Parallel Bisection in Floating
Point" by Demmel, Dhillon and Ren, LAPACK Working Note #70
.TP 8
M       (global output) INTEGER
Total number of eigenvalues found.  0 <= M <= N.
.TP 8
NZ      (global output) INTEGER
Total number of eigenvectors computed.  0 <= NZ <= M.
The number of columns of Z that are filled.
If JOBZ .NE. 'V', NZ is not referenced.
If JOBZ .EQ. 'V', NZ = M unless the user supplies
insufficient space and PDSYEVX is not able to detect this
before beginning computation.  To get all the eigenvectors
requested, the user must supply both sufficient
space to hold the eigenvectors in Z (M .LE. DESCZ(N_))
and sufficient workspace to compute them.  (See LWORK below.)
PDSYEVX is always able to detect insufficient space without
computation unless RANGE .EQ. 'V'.
.TP 8
W       (global output) DOUBLE PRECISION array, dimension (N)
On normal exit, the first M entries contain the selected
eigenvalues in ascending order.
.TP 8
ORFAC   (global input) DOUBLE PRECISION
Specifies which eigenvectors should be reorthogonalized.
Eigenvectors that correspond to eigenvalues which are within
tol=ORFAC*norm(A) of each other are to be reorthogonalized.
However, if the workspace is insufficient (see LWORK),
tol may be decreased until all eigenvectors to be
reorthogonalized can be stored in one process.
No reorthogonalization will be done if ORFAC equals zero.
A default value of 10^-3 is used if ORFAC is negative.
ORFAC should be identical on all processes.
.TP 8
Z       (local output) DOUBLE PRECISION array,
global dimension (N, N),
local dimension ( LLD_Z, LOCc(JZ+N-1) )
If JOBZ = 'V', then on normal exit the first M columns of Z
contain the orthonormal eigenvectors of the matrix
corresponding to the selected eigenvalues.  If an eigenvector
fails to converge, then that column of Z contains the latest
approximation to the eigenvector, and the index of the
eigenvector is returned in IFAIL.
If JOBZ = 'N', then Z is not referenced.
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
WORK    (local workspace/output) DOUBLE PRECISION array,
dimension (LWORK)
On return, WROK(1) contains the optimal amount of
workspace required for efficient execution.
if JOBZ='N' WORK(1) = optimal amount of workspace
required to compute eigenvalues efficiently
if JOBZ='V' WORK(1) = optimal amount of workspace
required to compute eigenvalues and eigenvectors
efficiently with no guarantee on orthogonality.
If RANGE='V', it is assumed that all eigenvectors
may be required.
.TP 8
LWORK   (local input) INTEGER
Size of WORK
See below for definitions of variables used to define LWORK.
If no eigenvectors are requested (JOBZ = 'N') then
LWORK >= 5 * N + MAX( 5 * NN, NB * ( NP0 + 1 ) )
If eigenvectors are requested (JOBZ = 'V' ) then
the amount of workspace required to guarantee that all
eigenvectors are computed is:
LWORK >= 5*N + MAX( 5*NN, NP0 * MQ0 + 2 * NB * NB ) +
ICEIL( NEIG, NPROW*NPCOL)*NN

The computed eigenvectors may not be orthogonal if the
minimal workspace is supplied and ORFAC is too small.
If you want to guarantee orthogonality (at the cost
of potentially poor performance) you should add
the following to LWORK:
(CLUSTERSIZE-1)*N
where CLUSTERSIZE is the number of eigenvalues in the
largest cluster, where a cluster is defined as a set of
close eigenvalues: { W(K),...,W(K+CLUSTERSIZE-1) |
W(J+1) <= W(J) + ORFAC*2*norm(A) }
Variable definitions:
NEIG = number of eigenvectors requested
NB = DESCA( MB_ ) = DESCA( NB_ ) =
DESCZ( MB_ ) = DESCZ( NB_ )
NN = MAX( N, NB, 2 )
DESCA( RSRC_ ) = DESCA( NB_ ) = DESCZ( RSRC_ ) =
DESCZ( CSRC_ ) = 0
NP0 = NUMROC( NN, NB, 0, 0, NPROW )
MQ0 = NUMROC( MAX( NEIG, NB, 2 ), NB, 0, 0, NPCOL )
ICEIL( X, Y ) is a ScaLAPACK function returning
ceiling(X/Y)

When LWORK is too small:
If LWORK is too small to guarantee orthogonality,
PDSYEVX attempts to maintain orthogonality in
the clusters with the smallest
spacing between the eigenvalues.
If LWORK is too small to compute all the eigenvectors
requested, no computation is performed and INFO=-23
is returned.  Note that when RANGE='V', PDSYEVX does
not know how many eigenvectors are requested until
the eigenvalues are computed.  Therefore, when RANGE='V'
and as long as LWORK is large enough to allow PDSYEVX to
compute the eigenvalues, PDSYEVX will compute the
eigenvalues and as many eigenvectors as it can.

Relationship between workspace, orthogonality & performance:
Greater performance can be achieved if adequate workspace
is provided.  On the other hand, in some situations,
performance can decrease as the workspace provided
increases above the workspace amount shown below:

For optimal performance, greater workspace may be
needed, i.e.
LWORK >=  MAX( LWORK, 5*N + NSYTRD_LWOPT )
Where:
LWORK, as defined previously, depends upon the number
of eigenvectors requested, and
NSYTRD_LWOPT = N + 2*( ANB+1 )*( 4*NPS+2 ) +
( NPS + 3 ) *  NPS

ANB = PJLAENV( DESCA( CTXT_), 3, 'PDSYTTRD', 'L',
0, 0, 0, 0)
SQNPC = INT( SQRT( DBLE( NPROW * NPCOL ) ) )
NPS = MAX( NUMROC( N, 1, 0, 0, SQNPC ), 2*ANB )

NUMROC is a ScaLAPACK tool functions;
PJLAENV is a ScaLAPACK envionmental inquiry function
MYROW, MYCOL, NPROW and NPCOL can be determined by
calling the subroutine BLACS_GRIDINFO.

For large N, no extra workspace is needed, however the
biggest boost in performance comes for small N, so it
is wise to provide the extra workspace (typically less
than a Megabyte per process).

If CLUSTERSIZE >= N/SQRT(NPROW*NPCOL), then providing
enough space to compute all the eigenvectors
orthogonally will cause serious degradation in
performance. In the limit (i.e. CLUSTERSIZE = N-1)
PDSTEIN will perform no better than DSTEIN on 1
processor.
For CLUSTERSIZE = N/SQRT(NPROW*NPCOL) reorthogonalizing
all eigenvectors will increase the total execution time
by a factor of 2 or more.
For CLUSTERSIZE > N/SQRT(NPROW*NPCOL) execution time will
grow as the square of the cluster size, all other factors
remaining equal and assuming enough workspace.  Less
workspace means less reorthogonalization but faster
execution.

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the size
required for optimal performance for all work arrays. Each of
these values is returned in the first entry of the
corresponding work arrays, and no error message is issued by
PXERBLA.
.TP 8
IWORK   (local workspace) INTEGER array
On return, IWORK(1) contains the amount of integer workspace
required.
.TP 8
LIWORK  (local input) INTEGER
size of IWORK
LIWORK >= 6 * NNP
Where:
NNP = MAX( N, NPROW*NPCOL + 1, 4 )
If LIWORK = -1, then LIWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
IFAIL   (global output) INTEGER array, dimension (N)
If JOBZ = 'V', then on normal exit, the first M elements of
IFAIL are zero.  If (MOD(INFO,2).NE.0) on exit, then
IFAIL contains the
indices of the eigenvectors that failed to converge.
If JOBZ = 'N', then IFAIL is not referenced.

ICLUSTR (global output) integer array, dimension (2*NPROW*NPCOL)
This array contains indices of eigenvectors corresponding to
a cluster of eigenvalues that could not be reorthogonalized
due to insufficient workspace (see LWORK, ORFAC and INFO).
Eigenvectors corresponding to clusters of eigenvalues indexed
ICLUSTR(2*I-1) to ICLUSTR(2*I), could not be
reorthogonalized due to lack of workspace. Hence the
eigenvectors corresponding to these clusters may not be
orthogonal.  ICLUSTR() is a zero terminated array.
(ICLUSTR(2*K).NE.0 .AND. ICLUSTR(2*K+1).EQ.0) if and only if
K is the number of clusters
ICLUSTR is not referenced if JOBZ = 'N'
.TP 8
GAP     (global output) DOUBLE PRECISION array,
dimension (NPROW*NPCOL)
This array contains the gap between eigenvalues whose
eigenvectors could not be reorthogonalized. The output
values in this array correspond to the clusters indicated
by the array ICLUSTR. As a result, the dot product between
eigenvectors correspoding to the I^th cluster may be as high
as ( C * n ) / GAP(I) where C is a small constant.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  if (MOD(INFO,2).NE.0), then one or more eigenvectors
failed to converge.  Their indices are stored
in IFAIL.  Ensure ABSTOL=2.0*PDLAMCH( 'U' )
Send e-mail to scalapack@cs.utk.edu
if (MOD(INFO/2,2).NE.0),then eigenvectors corresponding
to one or more clusters of eigenvalues could not be
reorthogonalized because of insufficient workspace.
The indices of the clusters are stored in the array
ICLUSTR.
if (MOD(INFO/4,2).NE.0), then space limit prevented
PDSYEVX from computing all of the eigenvectors
between VL and VU.  The number of eigenvectors
computed is returned in NZ.
if (MOD(INFO/8,2).NE.0), then PDSTEBZ failed to compute
eigenvalues.  Ensure ABSTOL=2.0*PDLAMCH( 'U' )
Send e-mail to scalapack@cs.utk.edu

Alignment requirements
======================

The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
must verify some alignment properties, namely the following
expressions should be true:

( MB_A.EQ.NB_A.EQ.MB_Z .AND. IROFFA.EQ.IROFFZ .AND. IROFFA.EQ.0 .AND.
IAROW.EQ.IZROW )
where
IROFFA = MOD( IA-1, MB_A ) and ICOFFA = MOD( JA-1, NB_A ).
