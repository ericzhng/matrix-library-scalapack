.TH PCLARZ l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCLARZ - applie a complex elementary reflector Q to a complex M-by-N distributed matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1), from either the left or the right
.SH SYNOPSIS
.TP 19
SUBROUTINE PCLARZ(
SIDE, M, N, L, V, IV, JV, DESCV, INCV, TAU, C,
IC, JC, DESCC, WORK )
.TP 19
.ti +4
CHARACTER
SIDE
.TP 19
.ti +4
INTEGER
IC, INCV, IV, JC, JV, L, M, N
.TP 19
.ti +4
INTEGER
DESCC( * ), DESCV( * )
.TP 19
.ti +4
COMPLEX
C( * ), TAU( * ), V( * ), WORK( * )
.SH PURPOSE
PCLARZ applies a complex elementary reflector Q to a complex M-by-N distributed matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1), from either the left or the right. Q is represented in the form 
      Q = I - tau * v * v'
.br

where tau is a complex scalar and v is a complex vector.
.br

If tau = 0, then Q is taken to be the unit matrix.
.br

Q is a product of k elementary reflectors as returned by PCTZRZF.

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

Because vectors may be viewed as a subclass of matrices, a
distributed vector is considered to be a distributed matrix.

Restrictions
.br
============
.br

If SIDE = 'Left' and INCV = 1, then the row process having the first
entry V(IV,JV) must also own C(IC+M-L,JC:JC+N-1). Moreover,
MOD(IV-1,MB_V) must be equal to MOD(IC+N-L-1,MB_C), if INCV=M_V, only
the last equality must be satisfied.
.br

If SIDE = 'Right' and INCV = M_V then the column process having the
first entry V(IV,JV) must also own C(IC:IC+M-1,JC+N-L) and
MOD(JV-1,NB_V) must be equal to MOD(JC+N-L-1,NB_C), if INCV = 1 only
the last equality must be satisfied.
.br

.SH ARGUMENTS
.TP 8
SIDE    (global input) CHARACTER
= 'L': form  Q * sub( C ),
.br
= 'R': form  sub( C ) * Q.
.TP 8
M       (global input) INTEGER
The number of rows to be operated on i.e the number of rows
of the distributed submatrix sub( C ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( C ). N >= 0.
.TP 8
L       (global input) INTEGER
The columns of the distributed submatrix sub( A ) containing
the meaningful part of the Householder reflectors.
If SIDE = 'L', M >= L >= 0, if SIDE = 'R', N >= L >= 0.
.TP 8
V       (local input) COMPLEX pointer into the local memory
to an array of dimension (LLD_V,*) containing the local
pieces of the distributed vectors V representing the
Householder transformation Q,
V(IV:IV+L-1,JV) if SIDE = 'L' and INCV = 1,
.br
V(IV,JV:JV+L-1) if SIDE = 'L' and INCV = M_V,
.br
V(IV:IV+L-1,JV) if SIDE = 'R' and INCV = 1,
.br
V(IV,JV:JV+L-1) if SIDE = 'R' and INCV = M_V,

The vector v in the representation of Q. V is not used if
TAU = 0.
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
INCV    (global input) INTEGER
The global increment for the elements of V. Only two values
of INCV are supported in this version, namely 1 and M_V.
INCV must not be zero.
.TP 8
TAU     (local input) COMPLEX, array, dimension  LOCc(JV) if
INCV = 1, and LOCr(IV) otherwise. This array contains the
Householder scalars related to the Householder vectors.
TAU is tied to the distributed matrix V.
.TP 8
C       (local input/local output) COMPLEX pointer into the
local memory to an array of dimension (LLD_C, LOCc(JC+N-1) ),
containing the local pieces of sub( C ). On exit, sub( C )
is overwritten by the Q * sub( C ) if SIDE = 'L', or
sub( C ) * Q if SIDE = 'R'.
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
WORK    (local workspace) COMPLEX array, dimension (LWORK)
If INCV = 1,
if SIDE = 'L',
if IVCOL = ICCOL,
LWORK >= NqC0
else
LWORK >= MpC0 + MAX( 1, NqC0 )
end if
else if SIDE = 'R',
LWORK >= NqC0 + MAX( MAX( 1, MpC0 ), NUMROC( NUMROC(
N+ICOFFC,NB_V,0,0,NPCOL ),NB_V,0,0,LCMQ ) )
end if
else if INCV = M_V,
if SIDE = 'L',
LWORK >= MpC0 + MAX( MAX( 1, NqC0 ), NUMROC( NUMROC(
M+IROFFC,MB_V,0,0,NPROW ),MB_V,0,0,LCMP ) )
else if SIDE = 'R',
if IVROW = ICROW,
LWORK >= MpC0
else
LWORK >= NqC0 + MAX( 1, MpC0 )
end if
end if
end if

where LCM is the least common multiple of NPROW and NPCOL and
LCM = ILCM( NPROW, NPCOL ), LCMP = LCM / NPROW,
LCMQ = LCM / NPCOL,

IROFFC = MOD( IC-1, MB_C ), ICOFFC = MOD( JC-1, NB_C ),
ICROW = INDXG2P( IC, MB_C, MYROW, RSRC_C, NPROW ),
ICCOL = INDXG2P( JC, NB_C, MYCOL, CSRC_C, NPCOL ),
MpC0 = NUMROC( M+IROFFC, MB_C, MYROW, ICROW, NPROW ),
NqC0 = NUMROC( N+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),

ILCM, INDXG2P and NUMROC are ScaLAPACK tool functions;
MYROW, MYCOL, NPROW and NPCOL can be determined by calling
the subroutine BLACS_GRIDINFO.

Alignment requirements
======================

The distributed submatrices V(IV:*, JV:*) and C(IC:IC+M-1,JC:JC+N-1)
must verify some alignment properties, namely the following
expressions should be true:

MB_V = NB_V,

If INCV = 1,
If SIDE = 'Left',
( MB_V.EQ.MB_C .AND. IROFFV.EQ.IROFFC .AND. IVROW.EQ.ICROW )
If SIDE = 'Right',
( MB_V.EQ.NB_A .AND. MB_V.EQ.NB_C .AND. IROFFV.EQ.ICOFFC )
else if INCV = M_V,
If SIDE = 'Left',
( MB_V.EQ.NB_V .AND. MB_V.EQ.MB_C .AND. ICOFFV.EQ.IROFFC )
If SIDE = 'Right',
( NB_V.EQ.NB_C .AND. ICOFFV.EQ.ICOFFC .AND. IVCOL.EQ.ICCOL )
end if
