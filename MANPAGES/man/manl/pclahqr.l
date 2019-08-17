.TH PCLAHQR l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCLAHQR - i an auxiliary routine used to find the Schur decomposition and or eigenvalues of a matrix already in Hessenberg form from cols ILO to IHI
.SH SYNOPSIS
.TP 20
SUBROUTINE PCLAHQR(
WANTT, WANTZ, N, ILO, IHI, A, DESCA, W, ILOZ,
IHIZ, Z, DESCZ, WORK, LWORK, IWORK, ILWORK,
INFO )
.TP 20
.ti +4
LOGICAL
WANTT, WANTZ
.TP 20
.ti +4
INTEGER
IHI, IHIZ, ILO, ILOZ, ILWORK, INFO, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCZ( * ), IWORK( * )
.TP 20
.ti +4
COMPLEX
A( * ), W( * ), WORK( * ), Z( * )
.SH PURPOSE
PCLAHQR is an auxiliary routine used to find the Schur decomposition and or eigenvalues of a matrix already in Hessenberg form from cols ILO to IHI. If Z = I, and WANTT=WANTZ=.TRUE., H gets replaced with Z'HZ,
  with Z'Z=I, and H in Schur form.
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
                               array.  LLD_A >= MAX(1,LOCp(M_A)).

Let K be the number of rows or columns of a distributed matrix,
and assume that its process grid has dimension p x q.
.br
LOCp( K ) denotes the number of elements of K that a process
would receive if K were distributed over the p processes of its
process column.
.br
Similarly, LOCq( K ) denotes the number of elements of K that a
process would receive if K were distributed over the q processes of
its process row.
.br
The values of LOCp() and LOCq() may be determined via a call to the
ScaLAPACK tool function, NUMROC:
.br
        LOCp( M ) = NUMROC( M, MB_A, MYROW, RSRC_A, NPROW ),
        LOCq( N ) = NUMROC( N, NB_A, MYCOL, CSRC_A, NPCOL ).
An upper bound for these quantities may be computed by:
.br
        LOCp( M ) <= ceil( ceil(M/MB_A)/NPROW )*MB_A
.br
        LOCq( N ) <= ceil( ceil(N/NB_A)/NPCOL )*NB_A
.br

.SH ARGUMENTS
.TP 8
WANTT   (global input) LOGICAL
= .TRUE. : the full Schur form T is required;
.br
= .FALSE.: only eigenvalues are required.
.TP 8
WANTZ   (global input) LOGICAL
.br
= .TRUE. : the matrix of Schur vectors Z is required;
.br
= .FALSE.: Schur vectors are not required.
.TP 8
N       (global input) INTEGER
The order of the Hessenberg matrix A (and Z if WANTZ).
N >= 0.
.TP 8
ILO     (global input) INTEGER
IHI     (global input) INTEGER
It is assumed that A is already upper quasi-triangular in
rows and columns IHI+1:N, and that A(ILO,ILO-1) = 0 (unless
ILO = 1). PCLAHQR works primarily with the Hessenberg
submatrix in rows and columns ILO to IHI, but applies
transformations to all of H if WANTT is .TRUE..
1 <= ILO <= max(1,IHI); IHI <= N.
.TP 8
A       (global input/output) COMPLEX array, dimension
(DESCA(LLD_),*)
On entry, the upper Hessenberg matrix A.
On exit, if WANTT is .TRUE., A is upper triangular in rows
and columns ILO:IHI.  If WANTT is .FALSE., the contents of
A are unspecified on exit.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 7
W      (global replicated output) COMPLEX array, dimension (N)
The computed eigenvalues ILO to IHI are stored in the
corresponding elements of W.  If WANTT is .TRUE., the
eigenvalues are stored in the same order as on the diagonal
of the Schur form returned in A.  A may be returned with
larger diagonal blocks until the next release.
.TP 8
ILOZ    (global input) INTEGER
IHIZ    (global input) INTEGER
Specify the rows of Z to which transformations must be
applied if WANTZ is .TRUE..
1 <= ILOZ <= ILO; IHI <= IHIZ <= N.
.TP 8
Z       (global input/output) COMPLEX array.
If WANTZ is .TRUE., on entry Z must contain the current
matrix Z of transformations accumulated by PCHSEQR, and on
exit Z has been updated; transformations are applied only to
the submatrix Z(ILOZ:IHIZ,ILO:IHI).
If WANTZ is .FALSE., Z is not referenced.
.TP 8
DESCZ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Z.
.TP 8
WORK    (local output) COMPLEX array of size LWORK
(Unless LWORK=-1, in which case WORK must be at least size 1)
.TP 8
LWORK   (local input) INTEGER
WORK(LWORK) is a local array and LWORK is assumed big enough
so that LWORK >= 3*N +
MAX( 2*MAX(DESCZ(LLD_),DESCA(LLD_)) + 2*LOCq(N),
7*Ceil(N/HBL)/LCM(NPROW,NPCOL)) +
MAX( 2*N, (8*LCM(NPROW,NPCOL)+2)**2 )
If LWORK=-1, then WORK(1) gets set to the above number and
the code returns immediately.
.TP 8
IWORK   (global and local input) INTEGER array of size ILWORK
This will hold some of the IBLK integer arrays.
This is held as a place holder for a future release.
Currently unreferenced.
.TP 8
ILWORK  (local input) INTEGER
This will hold the size of the IWORK array.
This is held as a place holder for a future release.
Currently unreferenced.
.TP 8
INFO    (global output) INTEGER
< 0: parameter number -INFO incorrect or inconsistent
.br
= 0: successful exit
.br
> 0: PCLAHQR failed to compute all the eigenvalues ILO to IHI
in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
elements i+1:ihi of W contains those eigenvalues
which have been successfully computed.

Logic:
This algorithm is very similar to SLAHQR.  Unlike SLAHQR,
instead of sending one double shift through the largest
unreduced submatrix, this algorithm sends multiple double shifts
and spaces them apart so that there can be parallelism across
several processor row/columns.  Another critical difference is
that this algorithm aggregrates multiple transforms together in
order to apply them in a block fashion.

Important Local Variables:
IBLK = The maximum number of bulges that can be computed.
Currently fixed.  Future releases this won't be fixed.
HBL  = The square block size (HBL=DESCA(MB_)=DESCA(NB_))
ROTN = The number of transforms to block together
NBULGE = The number of bulges that will be attempted on the
current submatrix.
IBULGE = The current number of bulges started.
K1(*),K2(*) = The current bulge loops from K1(*) to K2(*).

Subroutines:
From LAPACK, this routine calls:
CLAHQR     -> Serial QR used to determine shifts and
eigenvalues
CLARFG     -> Determine the Householder transforms

This ScaLAPACK, this routine calls:
PCLACONSB  -> To determine where to start each iteration
CLAMSH     -> Sends multiple shifts through a small
submatrix to see how the consecutive
subdiagonals change (if PCLACONSB indicates
we can start a run in the middle)
PCLAWIL    -> Given the shift, get the transformation
PCLACP3    -> Parallel array to local replicated array copy
& back.
CLAREF     -> Row/column reflector applier.  Core routine
here.
PCLASMSUB  -> Finds negligible subdiagonal elements.

Current Notes and/or Restrictions:
1.) This code requires the distributed block size to be square
and at least six (6); unlike simpler codes like LU, this
algorithm is extremely sensitive to block size.  Unwise
choices of too small a block size can lead to bad
performance.
2.) This code requires A and Z to be distributed identically
and have identical contxts.  A future version may allow Z to
have a different contxt to 1D row map it to all nodes (so no
communication on Z is necessary.)
3.) This code does not currently block the initial transforms
so that none of the rows or columns for any bulge are
completed until all are started.  To offset pipeline
start-up it is recommended that at least 2*LCM(NPROW,NPCOL)
bulges are used (if possible)
4.) The maximum number of bulges currently supported is fixed at
32.  In future versions this will be limited only by the
incoming WORK and IWORK array.
5.) The matrix A must be in upper Hessenberg form.  If elements
below the subdiagonal are nonzero, the resulting transforms
may be nonsimilar.  This is also true with the LAPACK
routine CLAHQR.
6.) For this release, this code has only been tested for
RSRC_=CSRC_=0, but it has been written for the general case.
7.) Currently, all the eigenvalues are distributed to all the
nodes.  Future releases will probably distribute the
eigenvalues by the column partitioning.
8.) The internals of this routine are subject to change.
9.) To optimize this for your architecture, try tuning CLAREF.
10.) This code has only been tested for WANTZ = .TRUE. and may
behave unpredictably for WANTZ set to .FALSE.
.SH FURTHER DETAILS
Contributed by Mark Fahey, June, 2000.
.br

