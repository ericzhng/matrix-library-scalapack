.TH PDLASMSUB l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLASMSUB - look for a small subdiagonal element from the bottom of the matrix that it can safely set to zero
.SH SYNOPSIS
.TP 22
SUBROUTINE PDLASMSUB(
A, DESCA, I, L, K, SMLNUM, BUF, LWORK )
.TP 22
.ti +4
INTEGER
I, K, L, LWORK
.TP 22
.ti +4
DOUBLE
PRECISION SMLNUM
.TP 22
.ti +4
INTEGER
DESCA( * )
.TP 22
.ti +4
DOUBLE
PRECISION A( * ), BUF( * )
.SH PURPOSE
PDLASMSUB looks for a small subdiagonal element from the bottom of the matrix that it can safely set to zero. 
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
A       (global input) DOUBLE PRECISION array, dimension
(DESCA(LLD_),*)
On entry, the Hessenberg matrix whose tridiagonal part is
being scanned.
Unchanged on exit.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
I       (global input) INTEGER
The global location of the bottom of the unreduced
submatrix of A.
Unchanged on exit.
.TP 8
L       (global input) INTEGER
The global location of the top of the unreduced submatrix
of A.
Unchanged on exit.
.TP 8
K       (global output) INTEGER
On exit, this yields the bottom portion of the unreduced
submatrix.  This will satisfy: L <= M  <= I-1.
.TP 8
SMLNUM  (global input) DOUBLE PRECISION
On entry, a "small number" for the given matrix.
Unchanged on exit.
.TP 8
BUF     (local output) DOUBLE PRECISION array of size LWORK.
.TP 8
LWORK   (global input) INTEGER
On exit, LWORK is the size of the work buffer.
This must be at least 2*Ceil( Ceil( (I-L)/HBL ) /
LCM(NPROW,NPCOL) )
Here LCM is least common multiple, and NPROWxNPCOL is the
logical grid size.

Notes:

This routine does a global maximum and must be called by all
processes.

This code is basically a parallelization of the following snip
of LAPACK code from DLAHQR:

Look for a single small subdiagonal element.

DO 20 K = I, L + 1, -1
TST1 = ABS( H( K-1, K-1 ) ) + ABS( H( K, K ) )
IF( TST1.EQ.ZERO )
$         TST1 = DLANHS( '1', I-L+1, H( L, L ), LDH, WORK )
IF( ABS( H( K, K-1 ) ).LE.MAX( ULP*TST1, SMLNUM ) )
$         GO TO 30
20    CONTINUE
30    CONTINUE

Implemented by:  G. Henry, November 17, 1996
