.TH PCSTEIN l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCSTEIN - compute the eigenvectors of a symmetric tridiagonal matrix in parallel, using inverse iteration
.SH SYNOPSIS
.TP 20
SUBROUTINE PCSTEIN(
N, D, E, M, W, IBLOCK, ISPLIT, ORFAC, Z, IZ,
JZ, DESCZ, WORK, LWORK, IWORK, LIWORK, IFAIL,
ICLUSTR, GAP, INFO )
.TP 20
.ti +4
INTEGER
INFO, IZ, JZ, LIWORK, LWORK, M, N
.TP 20
.ti +4
REAL
ORFAC
.TP 20
.ti +4
INTEGER
DESCZ( * ), IBLOCK( * ), ICLUSTR( * ),
IFAIL( * ), ISPLIT( * ), IWORK( * )
.TP 20
.ti +4
REAL
D( * ), E( * ), GAP( * ), W( * ), WORK( * )
.TP 20
.ti +4
COMPLEX
Z( * )
.SH PURPOSE
PCSTEIN computes the eigenvectors of a symmetric tridiagonal matrix in parallel, using inverse iteration. The eigenvectors found correspond to user specified eigenvalues. PCSTEIN does not
orthogonalize vectors that are on different processes. The extent
of orthogonalization is controlled by the input parameter LWORK.
Eigenvectors that are to be orthogonalized are computed by the same
process. PCSTEIN decides on the allocation of work among the
processes and then calls SSTEIN2 (modified LAPACK routine) on each
individual process. If insufficient workspace is allocated, the
expected orthogonalization may not be done.
.br

Note : If the eigenvectors obtained are not orthogonal, increase
       LWORK and run the code again.
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
and assume that its process grid has dimension r x c.
.br
LOCr( K ) denotes the number of elements of K that a process
would receive if K were distributed over the r processes of its
process column.
.br
Similarly, LOCc( K ) denotes the number of elements of K that a
process would receive if K were distributed over the c processes of
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


.SH ARGUMENTS

P = NPROW * NPCOL is the total number of processes
.TP 8
N       (global input) INTEGER
The order of the tridiagonal matrix T.  N >= 0.
.TP 8
D       (global input) REAL array, dimension (N)
The n diagonal elements of the tridiagonal matrix T.
.TP 8
E       (global input) REAL array, dimension (N-1)
The (n-1) off-diagonal elements of the tridiagonal matrix T.
.TP 8
M       (global input) INTEGER
The total number of eigenvectors to be found. 0 <= M <= N.
.TP 8
W       (global input/global output) REAL array, dim (M)
On input, the first M elements of W contain all the
eigenvalues for which eigenvectors are to be computed. The
eigenvalues should be grouped by split-off block and ordered
from smallest to largest within the block (The output array
W from PSSTEBZ with ORDER='b' is expected here). This
array should be replicated on all processes.
On output, the first M elements contain the input
eigenvalues in ascending order.

Note : To obtain orthogonal vectors, it is best if
eigenvalues are computed to highest accuracy ( this can be
done by setting ABSTOL to the underflow threshold =
SLAMCH('U') --- ABSTOL is an input parameter
to PSSTEBZ )
.TP 8
IBLOCK  (global input) INTEGER array, dimension (N)
The submatrix indices associated with the corresponding
eigenvalues in W -- 1 for eigenvalues belonging to the
first submatrix from the top, 2 for those belonging to
the second submatrix, etc. (The output array IBLOCK
from PSSTEBZ is expected here).
.TP 8
ISPLIT  (global input) INTEGER array, dimension (N)
The splitting points, at which T breaks up into submatrices.
The first submatrix consists of rows/columns 1 to ISPLIT(1),
the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
etc., and the NSPLIT-th consists of rows/columns
ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N (The output array
ISPLIT from PSSTEBZ is expected here.)
.TP 8
ORFAC   (global input) REAL
ORFAC specifies which eigenvectors should be orthogonalized.
Eigenvectors that correspond to eigenvalues which are within
ORFAC*||T|| of each other are to be orthogonalized.
However, if the workspace is insufficient (see LWORK), this
tolerance may be decreased until all eigenvectors to be
orthogonalized can be stored in one process.
No orthogonalization will be done if ORFAC equals zero.
A default value of 10^-3 is used if ORFAC is negative.
ORFAC should be identical on all processes.
.TP 8
Z       (local output) COMPLEX array,
dimension (DESCZ(DLEN_), N/npcol + NB)
Z contains the computed eigenvectors associated with the
specified eigenvalues. Any vector which fails to converge is
set to its current iterate after MAXITS iterations ( See
SSTEIN2 ).
On output, Z is distributed across the P processes in block
cyclic format.
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
.TP 8
WORK    (local workspace/global output) REAL array,
dimension ( LWORK )
On output, WORK(1) gives a lower bound on the
workspace ( LWORK ) that guarantees the user desired
orthogonalization (see ORFAC).
Note that this may overestimate the minimum workspace needed.
.TP 8
LWORK   (local input) integer
LWORK  controls the extent of orthogonalization which can be
done. The number of eigenvectors for which storage is
allocated on each process is
NVEC = floor(( LWORK- max(5*N,NP00*MQ00) )/N).
Eigenvectors corresponding to eigenvalue clusters of size
NVEC - ceil(M/P) + 1 are guaranteed to be orthogonal ( the
orthogonality is similar to that obtained from CSTEIN2).
Note : LWORK  must be no smaller than:
max(5*N,NP00*MQ00) + ceil(M/P)*N,
and should have the same input value on all processes.
It is the minimum value of LWORK input on different processes
that is significant.

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
IWORK   (local workspace/global output) INTEGER array,
dimension ( 3*N+P+1 )
On return, IWORK(1) contains the amount of integer workspace
required.
On return, the IWORK(2) through IWORK(P+2) indicate
the eigenvectors computed by each process. Process I computes
eigenvectors indexed IWORK(I+2)+1 thru' IWORK(I+3).
.TP 8
LIWORK  (local input) INTEGER
Size of array IWORK.  Must be >= 3*N + P + 1

If LIWORK = -1, then LIWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
IFAIL   (global output) integer array, dimension (M)
On normal exit, all elements of IFAIL are zero.
If one or more eigenvectors fail to converge after MAXITS
iterations (as in CSTEIN), then INFO > 0 is returned.
If mod(INFO,M+1)>0, then
for I=1 to mod(INFO,M+1), the eigenvector
corresponding to the eigenvalue W(IFAIL(I)) failed to
converge ( W refers to the array of eigenvalues on output ).

ICLUSTR (global output) integer array, dimension (2*P)
This output array contains indices of eigenvectors
corresponding to a cluster of eigenvalues that could not be
orthogonalized due to insufficient workspace (see LWORK,
ORFAC and INFO). Eigenvectors corresponding to clusters of
eigenvalues indexed ICLUSTR(2*I-1) to ICLUSTR(2*I), I = 1 to
INFO/(M+1), could not be orthogonalized due to lack of
workspace. Hence the eigenvectors corresponding to these
clusters may not be orthogonal. ICLUSTR is a zero terminated
array --- ( ICLUSTR(2*K).NE.0 .AND. ICLUSTR(2*K+1).EQ.0 )
if and only if K is the number of clusters.
.TP 8
GAP     (global output) REAL array, dimension (P)
This output array contains the gap between eigenvalues whose
eigenvectors could not be orthogonalized. The INFO/M output
values in this array correspond to the INFO/(M+1) clusters
indicated by the array ICLUSTR. As a result, the dot product
between eigenvectors corresponding to the I^th cluster may be
as high as ( O(n)*macheps ) / GAP(I).
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
< 0 :  if INFO = -I, the I-th argument had an illegal value
.br
> 0 :  if mod(INFO,M+1) = I, then I eigenvectors failed to
converge in MAXITS iterations. Their indices are
stored in the array IFAIL.
if INFO/(M+1) = I, then eigenvectors corresponding to
I clusters of eigenvalues could not be orthogonalized
due to insufficient workspace. The indices of the
clusters are stored in the array ICLUSTR.
