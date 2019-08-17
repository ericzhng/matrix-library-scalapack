.TH PDORMQL l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDORMQL - overwrite the general real M-by-N distributed matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with  SIDE = 'L' SIDE = 'R' TRANS = 'N'
.SH SYNOPSIS
.TP 20
SUBROUTINE PDORMQL(
SIDE, TRANS, M, N, K, A, IA, JA, DESCA, TAU,
C, IC, JC, DESCC, WORK, LWORK, INFO )
.TP 20
.ti +4
CHARACTER
SIDE, TRANS
.TP 20
.ti +4
INTEGER
IA, IC, INFO, JA, JC, K, LWORK, M, N
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCC( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), C( * ), TAU( * ), WORK( * )
.SH PURPOSE
PDORMQL overwrites the general real M-by-N distributed matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with SIDE = 'L' SIDE = 'R' TRANS = 'N': Q * sub( C ) sub( C ) * Q TRANS = 'T':      Q**T * sub( C )       sub( C ) * Q**T
.br

where Q is a real orthogonal distributed matrix defined as the
product of K elementary reflectors
.br

      Q = H(k) . . . H(2) H(1)
.br

as returned by PDGEQLF. Q is of order M if SIDE = 'L' and of order N
if SIDE = 'R'.
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
SIDE    (global input) CHARACTER
= 'L': apply Q or Q**T from the Left;
.br
= 'R': apply Q or Q**T from the Right.
.TP 8
TRANS   (global input) CHARACTER
.br
= 'N':  No transpose, apply Q;
.br
= 'T':  Transpose, apply Q**T.
.TP 8
M       (global input) INTEGER
The number of rows to be operated on i.e the number of rows
of the distributed submatrix sub( C ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( C ). N >= 0.
.TP 8
K       (global input) INTEGER
The number of elementary reflectors whose product defines the
matrix Q.  If SIDE = 'L', M >= K >= 0, if SIDE = 'R',
N >= K >= 0.
.TP 8
A       (local input) DOUBLE PRECISION pointer into the local memory
to an array of dimension (LLD_A,LOCc(JA+K-1)). On entry, the
j-th column must contain the vector which defines the elemen-
tary reflector H(j), JA <= j <= JA+K-1, as returned by
PDGEQLF in the K columns of its distributed matrix
argument A(IA:*,JA:JA+K-1). A(IA:*,JA:JA+K-1) is modified by
the routine but restored on exit.
If SIDE = 'L', LLD_A >= MAX( 1, LOCr(IA+M-1) ),
if SIDE = 'R', LLD_A >= MAX( 1, LOCr(IA+N-1) ).
.TP 8
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
TAU     (local input) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
This array contains the scalar factors TAU(j) of the
elementary reflectors H(j) as returned by PDGEQLF.
TAU is tied to the distributed matrix A.
.TP 8
C       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
On entry, the local pieces of the distributed matrix sub(C).
On exit, sub( C ) is overwritten by Q*sub( C ) or Q'*sub( C )
or sub( C )*Q' or sub( C )*Q.
.TP 8
IC      (global input) INTEGER
The row index in the global array C indicating the first
row of sub( C ).
.TP 8
JC      (global input) INTEGER
The column index in the global array C indicating the
first column of sub( C ).
.TP 8
DESCC   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix C.
.TP 8
WORK    (local workspace/local output) DOUBLE PRECISION array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
If SIDE = 'L',
LWORK >= MAX( (NB_A*(NB_A-1))/2, (NqC0 + MpC0)*NB_A ) +
NB_A * NB_A
else if SIDE = 'R',
LWORK >= MAX( (NB_A*(NB_A-1))/2, ( NqC0 + MAX( NpA0 +
NUMROC( NUMROC( N+ICOFFC, NB_A, 0, 0, NPCOL ),
NB_A, 0, 0, LCMQ ), MpC0 ) )*NB_A ) +
NB_A * NB_A
end if

where LCMQ = LCM / NPCOL with LCM = ICLM( NPROW, NPCOL ),

IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
NpA0 = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),

IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),

ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
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

Alignment requirements
======================

The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
must verify some alignment properties, namely the following
expressions should be true:

If SIDE = 'L',
( MB_A.EQ.MB_C .AND. IROFFA.EQ.IROFFC .AND. IAROW.EQ.ICROW )
If SIDE = 'R',
( MB_A.EQ.NB_C .AND. IROFFA.EQ.ICOFFC )
