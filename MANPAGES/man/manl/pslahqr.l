.TH PSLAHQR l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAHQR - i an auxiliary routine used to find the Schur decomposition and or eigenvalues of a matrix already in Hessenberg form from cols ILO to IHI
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLAHQR(
WANTT, WANTZ, N, ILO, IHI, A, DESCA, WR, WI,
ILOZ, IHIZ, Z, DESCZ, WORK, LWORK, IWORK,
ILWORK, INFO )
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
REAL
A( * ), WI( * ), WORK( * ), WR( * ), Z( * )
.SH PURPOSE
PSLAHQR is an auxiliary routine used to find the Schur decomposition and or eigenvalues of a matrix already in Hessenberg form from cols ILO to IHI. 
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
ILO = 1). PSLAHQR works primarily with the Hessenberg
submatrix in rows and columns ILO to IHI, but applies
transformations to all of H if WANTT is .TRUE..
1 <= ILO <= max(1,IHI); IHI <= N.
.TP 8
A       (global input/output) REAL array, dimension
(DESCA(LLD_),*)
On entry, the upper Hessenberg matrix A.
On exit, if WANTT is .TRUE., A is upper quasi-triangular in
rows and columns ILO:IHI, with any 2-by-2 or larger diagonal
blocks not yet in standard form. If WANTT is .FALSE., the
contents of A are unspecified on exit.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
WR      (global replicated output) REAL array,
dimension (N)
WI      (global replicated output) REAL array,
dimension (N)
The real and imaginary parts, respectively, of the computed
eigenvalues ILO to IHI are stored in the corresponding
elements of WR and WI. If two eigenvalues are computed as a
complex conjugate pair, they are stored in consecutive
elements of WR and WI, say the i-th and (i+1)th, with
WI(i) > 0 and WI(i+1) < 0. If WANTT is .TRUE., the
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
Z       (global input/output) REAL array.
If WANTZ is .TRUE., on entry Z must contain the current
matrix Z of transformations accumulated by PDHSEQR, and on
exit Z has been updated; transformations are applied only to
the submatrix Z(ILOZ:IHIZ,ILO:IHI).
If WANTZ is .FALSE., Z is not referenced.
.TP 8
DESCZ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Z.
.TP 8
WORK    (local output) REAL array of size LWORK
.TP 8
LWORK   (local input) INTEGER
WORK(LWORK) is a local array and LWORK is assumed big enough
so that LWORK >= 3*N +
MAX( 2*MAX(DESCZ(LLD_),DESCA(LLD_)) + 2*LOCc(N),
7*Ceil(N/HBL)/LCM(NPROW,NPCOL)) )
.TP 8
IWORK   (global and local input) INTEGER array of size ILWORK
.TP 8
ILWORK  (local input) INTEGER
This holds the some of the IBLK integer arrays.  This is held
as a place holder for the next release.
.TP 8
INFO    (global output) INTEGER
< 0: parameter number -INFO incorrect or inconsistent
.br
= 0: successful exit
.br
> 0: PSLAHQR failed to compute all the eigenvalues ILO to IHI
in a total of 30*(IHI-ILO+1) iterations; if INFO = i,
elements i+1:ihi of WR and WI contain those eigenvalues
which have been successfully computed.

Logic:
This algorithm is very similar to _LAHQR.  Unlike _LAHQR,
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
This routine calls:
PSLACONSB   -> To determine where to start each iteration
PSLAWIL   -> Given the shift, get the transformation
SLASORTE   -> Pair up eigenvalues so that reals are paired.
PSLACP3   -> Parallel array to local replicated array copy &
back.
SLAREF   -> Row/column reflector applier.  Core routine
here.
PSLASMSUB   -> Finds negligible subdiagonal elements.

Current Notes and/or Restrictions:
1.) This code requires the distributed block size to be square
and at least six (6); unlike simpler codes like LU, this
algorithm is extremely sensitive to block size.  Unwise
choices of too small a block size can lead to bad
performance.
2.) This code requires A and Z to be distributed identically
and have identical contxts.
3.) This release currently does not have a routine for
resolving the Schur blocks into regular 2x2 form after
this code is completed.  Because of this, a significant
performance impact is required while the deflation is done
by sometimes a single column of processors.
4.) This code does not currently block the initial transforms
so that none of the rows or columns for any bulge are
completed until all are started.  To offset pipeline
start-up it is recommended that at least 2*LCM(NPROW,NPCOL)
bulges are used (if possible)
5.) The maximum number of bulges currently supported is fixed at
32.  In future versions this will be limited only by the
incoming WORK array.
6.) The matrix A must be in upper Hessenberg form.  If elements
below the subdiagonal are nonzero, the resulting transforms
may be nonsimilar.  This is also true with the LAPACK
routine.
7.) For this release, it is assumed RSRC_=CSRC_=0
8.) Currently, all the eigenvalues are distributed to all the
nodes.  Future releases will probably distribute the
eigenvalues by the column partitioning.
9.) The internals of this routine are subject to change.

Implemented by:  G. Henry, November 17, 1996
