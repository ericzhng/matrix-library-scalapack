.TH PDGEHRD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDGEHRD - reduce a real general distributed matrix sub( A ) to upper Hessenberg form H by an orthogonal similarity transforma- tion
.SH SYNOPSIS
.TP 20
SUBROUTINE PDGEHRD(
N, ILO, IHI, A, IA, JA, DESCA, TAU, WORK,
LWORK, INFO )
.TP 20
.ti +4
INTEGER
IA, IHI, ILO, INFO, JA, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), TAU( * ), WORK( * )
.SH PURPOSE
PDGEHRD reduces a real general distributed matrix sub( A ) to upper Hessenberg form H by an orthogonal similarity transforma- tion: Q' * sub( A ) * Q = H, where sub( A ) = A(IA+N-1:IA+N-1,JA+N-1:JA+N-1).
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

.SH ARGUMENTS
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
ILO     (global input) INTEGER
IHI     (global input) INTEGER
It is assumed that sub( A ) is already upper triangular in
rows IA:IA+ILO-2 and IA+IHI:IA+N-1 and columns JA:JA+ILO-2
and JA+IHI:JA+N-1. See Further Details. If N > 0,
.TP 8
A       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
On entry, this array contains the local pieces of the N-by-N
general distributed matrix sub( A ) to be reduced. On exit,
the upper triangle and the first subdiagonal of sub( A ) are
overwritten with the upper Hessenberg matrix H, and the ele-
ments below the first subdiagonal, with the array TAU, repre-
sent the orthogonal matrix Q as a product of elementary
reflectors. See Further Details.
IA      (global input) INTEGER
The row index in the global array A indicating the first
row of sub( A ).
.TP 8
JA      (global input) INTEGER
The column index in the global array A indicating the
first column of sub( A ).
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
TAU     (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-2)
The scalar factors of the elementary reflectors (see Further
Details). Elements JA:JA+ILO-2 and JA+IHI:JA+N-2 of TAU are
set to zero. TAU is tied to the distributed matrix A.
.TP 8
WORK    (local workspace/local output) DOUBLE PRECISION array,
dimension (LWORK)
On exit, WORK( 1 ) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= NB*NB + NB*MAX( IHIP+1, IHLP+INLQ )

where NB = MB_A = NB_A, IROFFA = MOD( IA-1, NB ),
ICOFFA = MOD( JA-1, NB ), IOFF = MOD( IA+ILO-2, NB ),
IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
IHIP = NUMROC( IHI+IROFFA, NB, MYROW, IAROW, NPROW ),
ILROW = INDXG2P( IA+ILO-1, NB, MYROW, RSRC_A, NPROW ),
IHLP = NUMROC( IHI-ILO+IOFF+1, NB, MYROW, ILROW, NPROW ),
ILCOL = INDXG2P( JA+ILO-1, NB, MYCOL, CSRC_A, NPCOL ),
INLQ = NUMROC( N-ILO+IOFF+1, NB, MYCOL, ILCOL, NPCOL ),

INDXG2P and NUMROC are ScaLAPACK tool functions;
MYROW, MYCOL, NPROW and NPCOL can be determined by calling
the subroutine BLACS_GRIDINFO.

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
.SH FURTHER DETAILS
The matrix Q is represented as a product of (ihi-ilo) elementary
reflectors
.br

   Q = H(ilo) H(ilo+1) . . . H(ihi-1).
.br

Each H(i) has the form
.br

   H(i) = I - tau * v * v'
.br

where tau is a real scalar, and v is a real vector with
.br
v(1:I) = 0, v(I+1) = 1 and v(IHI+1:N) = 0; v(I+2:IHI) is stored on
exit in A(IA+ILO+I:IA+IHI-1,JA+ILO+I-2), and tau in TAU(JA+ILO+I-2).

The contents of A(IA:IA+N-1,JA:JA+N-1) are illustrated by the follow-
ing example, with N = 7, ILO = 2 and IHI = 6:
.br

on entry                         on exit
.br

( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
(     a   a   a   a   a   a )    (      a   h   h   h   h   a )
(     a   a   a   a   a   a )    (      h   h   h   h   h   h )
(     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
(     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
(     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
(                         a )    (                          a )

where a denotes an element of the original matrix sub( A ), H denotes
a modified element of the upper Hessenberg matrix H, and vi denotes
an element of the vector defining H(JA+ILO+I-2).
.br

Alignment requirements
.br
======================
.br

The distributed submatrix sub( A ) must verify some alignment proper-
ties, namely the following expression should be true:
.br
( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
.br

