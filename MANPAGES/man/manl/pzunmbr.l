.TH PZUNMBR l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZUNMBR - VECT = 'Q', PZUNMBR overwrites the general complex distributed M-by-N matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with  SIDE = 'L' SIDE = 'R' TRANS = 'N'
.SH SYNOPSIS
.TP 20
SUBROUTINE PZUNMBR(
VECT, SIDE, TRANS, M, N, K, A, IA, JA, DESCA,
TAU, C, IC, JC, DESCC, WORK, LWORK, INFO )
.TP 20
.ti +4
CHARACTER
SIDE, TRANS, VECT
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
COMPLEX*16
A( * ), C( * ), TAU( * ), WORK( * )
.SH PURPOSE
If VECT = 'Q', PZUNMBR overwrites the general complex distributed M-by-N matrix sub( C ) = C(IC:IC+M-1,JC:JC+N-1) with SIDE = 'L' SIDE = 'R' TRANS = 'N': Q * sub( C ) sub( C ) * Q TRANS = 'C':      Q**H * sub( C )       sub( C ) * Q**H
.br

If VECT = 'P', PZUNMBR overwrites sub( C ) with
.br

                     SIDE = 'L'           SIDE = 'R'
.br
TRANS = 'N':      P * sub( C )          sub( C ) * P
.br
TRANS = 'C':      P**H * sub( C )       sub( C ) * P**H
.br

Here Q and P**H are the unitary distributed matrices determined by
PZGEBRD when reducing a complex distributed matrix A(IA:*,JA:*) to
bidiagonal form: A(IA:*,JA:*) = Q * B * P**H. Q and P**H are defined
as products of elementary reflectors H(i) and G(i) respectively.

Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
order of the unitary matrix Q or P**H that is applied.
.br

If VECT = 'Q', A(IA:*,JA:*) is assumed to have been an NQ-by-K
matrix:
.br
if nq >= k, Q = H(1) H(2) . . . H(k);
.br
if nq < k, Q = H(1) H(2) . . . H(nq-1).
.br

If VECT = 'P', A(IA:*,JA:*) is assumed to have been a K-by-NQ
matrix:
.br
if k < nq, P = G(1) G(2) . . . G(k);
.br
if k >= nq, P = G(1) G(2) . . . G(nq-1).
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
VECT    (global input) CHARACTER
= 'Q': apply Q or Q**H;
.br
= 'P': apply P or P**H.
.TP 8
SIDE    (global input) CHARACTER
.br
= 'L': apply Q, Q**H, P or P**H from the Left;
.br
= 'R': apply Q, Q**H, P or P**H from the Right.
.TP 8
TRANS   (global input) CHARACTER
.br
= 'N':  No transpose, apply Q or P;
.br
= 'C':  Conjugate transpose, apply Q**H or P**H.
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
If VECT = 'Q', the number of columns in the original
distributed matrix reduced by PZGEBRD.
If VECT = 'P', the number of rows in the original
distributed matrix reduced by PZGEBRD.
K >= 0.
.TP 8
A       (local input) COMPLEX*16 pointer into the local memory
to an array of dimension (LLD_A,LOCc(JA+MIN(NQ,K)-1)) if
VECT='Q', and (LLD_A,LOCc(JA+NQ-1)) if VECT = 'P'. NQ = M
if SIDE = 'L', and NQ = N otherwise. The vectors which
define the elementary reflectors H(i) and G(i), whose
products determine the matrices Q and P, as returned by
PZGEBRD.
If VECT = 'Q', LLD_A >= max(1,LOCr(IA+NQ-1));
if VECT = 'P', LLD_A >= max(1,LOCr(IA+MIN(NQ,K)-1)).
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
TAU     (local input) COMPLEX*16 array, dimension
LOCc(JA+MIN(NQ,K)-1) if VECT = 'Q', LOCr(IA+MIN(NQ,K)-1) if
VECT = 'P', TAU(i) must contain the scalar factor of the
elementary  reflector H(i) or G(i), which determines Q or P,
as returned by PDGEBRD in its array argument TAUQ or TAUP.
TAU is tied to the distributed matrix A.
.TP 8
C       (local input/local output) COMPLEX*16 pointer into the
local memory to an array of dimension (LLD_C,LOCc(JC+N-1)).
On entry, the local pieces of the distributed matrix sub(C).
On exit, if VECT='Q', sub( C ) is overwritten by Q*sub( C )
or Q'*sub( C ) or sub( C )*Q' or sub( C )*Q; if VECT='P,
sub( C ) is overwritten by P*sub( C ) or P'*sub( C ) or
sub( C )*P or sub( C )*P'.
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
WORK    (local workspace/local output) COMPLEX*16 array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
If SIDE = 'L',
NQ = M;
if( (VECT = 'Q' and NQ >= K) or (VECT <> 'Q' and NQ > K) ),
IAA=IA; JAA=JA; MI=M; NI=N; ICC=IC; JCC=JC;
else
IAA=IA+1; JAA=JA; MI=M-1; NI=N; ICC=IC+1; JCC=JC;
end if
else if SIDE = 'R',
NQ = N;
if( (VECT = 'Q' and NQ >= K) or (VECT <> 'Q' and NQ > K) ),
IAA=IA; JAA=JA; MI=M; NI=N; ICC=IC; JCC=JC;
else
IAA=IA; JAA=JA+1; MI=M; NI=N-1; ICC=IC; JCC=JC+1;
end if
end if

If VECT = 'Q',
If SIDE = 'L',
LWORK >= MAX( (NB_A*(NB_A-1))/2, (NqC0 + MpC0)*NB_A ) +
NB_A * NB_A
else if SIDE = 'R',
LWORK >= MAX( (NB_A*(NB_A-1))/2, ( NqC0 + MAX( NpA0 +
NUMROC( NUMROC( NI+ICOFFC, NB_A, 0, 0, NPCOL ),
NB_A, 0, 0, LCMQ ), MpC0 ) )*NB_A ) +
NB_A * NB_A
end if
else if VECT <> 'Q',
if SIDE = 'L',
LWORK >= MAX( (MB_A*(MB_A-1))/2, ( MpC0 + MAX( MqA0 +
NUMROC( NUMROC( MI+IROFFC, MB_A, 0, 0, NPROW ),
MB_A, 0, 0, LCMP ), NqC0 ) )*MB_A ) +
MB_A * MB_A
else if SIDE = 'R',
LWORK >= MAX( (MB_A*(MB_A-1))/2, (MpC0 + NqC0)*MB_A ) +
MB_A * MB_A
end if
end if

where LCMP = LCM / NPROW, LCMQ = LCM / NPCOL, with
LCM = ICLM( NPROW, NPCOL ),

IROFFA = MOD( IAA-1, MB_A ), ICOFFA = MOD( JAA-1, NB_A ),
IAROW = INDXG2P( IAA, MB_A, MYROW, RSRC_A, NPROW ),
IACOL = INDXG2P( JAA, NB_A, MYCOL, CSRC_A, NPCOL ),
MqA0 = NUMROC( MI+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),
NpA0 = NUMROC( NI+IROFFA, MB_A, MYROW, IAROW, NPROW ),

IROFFC = MOD( ICC-1, MB_C ), ICOFFC = MOD( JCC-1, NB_C ),
ICROW = INDXG2P( ICC, MB_C, MYROW, RSRC_C, NPROW ),
ICCOL = INDXG2P( JCC, NB_C, MYCOL, CSRC_C, NPCOL ),
MpC0 = NUMROC( MI+IROFFC, MB_C, MYROW, ICROW, NPROW ),
NqC0 = NUMROC( NI+ICOFFC, NB_C, MYCOL, ICCOL, NPCOL ),

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

Alignment requirements
======================

The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
must verify some alignment properties, namely the following
expressions should be true:

If VECT = 'Q',
If SIDE = 'L',
( MB_A.EQ.MB_C .AND. IROFFA.EQ.IROFFC .AND. IAROW.EQ.ICROW )
If SIDE = 'R',
( MB_A.EQ.NB_C .AND. IROFFA.EQ.ICOFFC )
else
If SIDE = 'L',
( MB_A.EQ.MB_C .AND. ICOFFA.EQ.IROFFC )
If SIDE = 'R',
( NB_A.EQ.NB_C .AND. ICOFFA.EQ.ICOFFC .AND. IACOL.EQ.ICCOL )
end if
