.TH PDLARFB l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLARFB - applie a real block reflector Q or its transpose Q**T to a real distributed M-by-N matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1)
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLARFB(
SIDE, TRANS, DIRECT, STOREV, M, N, K, V, IV,
JV, DESCV, T, C, IC, JC, DESCC, WORK )
.TP 20
.ti +4
CHARACTER
SIDE, TRANS, DIRECT, STOREV
.TP 20
.ti +4
INTEGER
IC, IV, JC, JV, K, M, N
.TP 20
.ti +4
INTEGER
DESCC( * ), DESCV( * )
.TP 20
.ti +4
DOUBLE
PRECISION C( * ), T( * ), V( * ), WORK( * )
.SH PURPOSE
PDLARFB applies a real block reflector Q or its transpose Q**T to a real distributed M-by-N matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) from the left or the right.
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
DIRECT  (global input) CHARACTER
Indicates how Q is formed from a product of elementary
reflectors
= 'F': Q = H(1) H(2) . . . H(k) (Forward)
.br
= 'B': Q = H(k) . . . H(2) H(1) (Backward)
.TP 8
STOREV  (global input) CHARACTER
Indicates how the vectors which define the elementary
reflectors are stored:
.br
= 'C': Columnwise
.br
= 'R': Rowwise
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
The order of the matrix T (= the number of elementary
reflectors whose product defines the block reflector).
.TP 8
V       (local input) DOUBLE PRECISION pointer into the local memory
to an array of dimension ( LLD_V, LOCc(JV+K-1) ) if
STOREV = 'C', ( LLD_V, LOCc(JV+M-1)) if STOREV = 'R' and
SIDE = 'L', ( LLD_V, LOCc(JV+N-1) ) if STOREV = 'R' and
SIDE = 'R'. It contains the local pieces of the distributed
vectors V representing the Householder transformation.
See further details.
If STOREV = 'C' and SIDE = 'L', LLD_V >= MAX(1,LOCr(IV+M-1));
if STOREV = 'C' and SIDE = 'R', LLD_V >= MAX(1,LOCr(IV+N-1));
if STOREV = 'R', LLD_V >= LOCr(IV+K-1).
.TP 8
IV      (global input) INTEGER
The row index in the global array V indicating the first
row of sub( V ).
.TP 8
JV      (global input) INTEGER
The column index in the global array V indicating the
first column of sub( V ).
.TP 8
DESCV   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix V.
.TP 8
T       (local input) DOUBLE PRECISION array, dimension MB_V by MB_V
if STOREV = 'R' and NB_V by NB_V if STOREV = 'C'. The trian-
gular matrix T in the representation of the block reflector.
.TP 8
C       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
On entry, the M-by-N distributed matrix sub( C ). On exit,
sub( C ) is overwritten by Q*sub( C ) or Q'*sub( C ) or
sub( C )*Q or sub( C )*Q'.
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
WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
If STOREV = 'C',
if SIDE = 'L',
LWORK >= ( NqC0 + MpC0 ) * K
else if SIDE = 'R',
LWORK >= ( NqC0 + MAX( NpV0 + NUMROC( NUMROC( N+ICOFFC,
NB_V, 0, 0, NPCOL ), NB_V, 0, 0, LCMQ ),
MpC0 ) ) * K
end if
else if STOREV = 'R',
if SIDE = 'L',
LWORK >= ( MpC0 + MAX( MqV0 + NUMROC( NUMROC( M+IROFFC,
MB_V, 0, 0, NPROW ), MB_V, 0, 0, LCMP ),
NqC0 ) ) * K
else if SIDE = 'R',
LWORK >= ( MpC0 + NqC0 ) * K
end if
end if

where LCMQ = LCM / NPCOL with LCM = ICLM( NPROW, NPCOL ),

IROFFV = MOD( IV-1, MB_V ), ICOFFV = MOD( JV-1, NB_V ),
IVROW = INDXG2P( IV, MB_V, MYROW, RSRC_V, NPROW ),
IVCOL = INDXG2P( JV, NB_V, MYCOL, CSRC_V, NPCOL ),
MqV0 = NUMROC( M+ICOFFV, NB_V, MYCOL, IVCOL, NPCOL ),
NpV0 = NUMROC( N+IROFFV, MB_V, MYROW, IVROW, NPROW ),

IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
NpC0 = NUMROC( N+ICOFFC, MB_C, MYROW, ICROW, NPROW ),
NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),

ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
MYROW, MYCOL, NPROW and NPCOL can be determined by calling
the subroutine BLACS_GRIDINFO.

Alignment requirements
======================

The distributed submatrices V(IV:*, JV:*) and C(IC:IC+M-1,JC:JC+N-1)
must verify some alignment properties, namely the following
expressions should be true:

If STOREV = 'Columnwise'
If SIDE = 'Left',
( MB_V.EQ.MB_C .AND. IROFFV.EQ.IROFFC .AND. IVROW.EQ.ICROW )
If SIDE = 'Right',
( MB_V.EQ.NB_C .AND. IROFFV.EQ.ICOFFC )
else if STOREV = 'Rowwise'
If SIDE = 'Left',
( NB_V.EQ.MB_C .AND. ICOFFV.EQ.IROFFC )
If SIDE = 'Right',
( NB_V.EQ.NB_C .AND. ICOFFV.EQ.ICOFFC .AND. IVCOL.EQ.ICCOL )
end if
