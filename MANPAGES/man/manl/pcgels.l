.TH PCGELS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCGELS - solve overdetermined or underdetermined complex linear systems involving an M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1),
.SH SYNOPSIS
.TP 19
SUBROUTINE PCGELS(
TRANS, M, N, NRHS, A, IA, JA, DESCA, B, IB, JB,
DESCB, WORK, LWORK, INFO )
.TP 19
.ti +4
CHARACTER
TRANS
.TP 19
.ti +4
INTEGER
IA, IB, INFO, JA, JB, LWORK, M, N, NRHS
.TP 19
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 19
.ti +4
COMPLEX
A( * ), B( * ), WORK( * )
.SH PURPOSE
PCGELS solves overdetermined or underdetermined complex linear systems involving an M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1), or its conjugate-transpose, using a QR or LQ factorization of
sub( A ).  It is assumed that sub( A ) has full rank.
.br

The following options are provided:
.br

1. If TRANS = 'N' and m >= n:  find the least squares solution of
   an overdetermined system, i.e., solve the least squares problem
                minimize || sub( B ) - sub( A )*X ||.
.br

2. If TRANS = 'N' and m < n:  find the minimum norm solution of
   an underdetermined system sub( A ) * X = sub( B ).
.br

3. If TRANS = 'C' and m >= n:  find the minimum norm solution of
   an undetermined system sub( A )**H * X = sub( B ).
.br

4. If TRANS = 'C' and m < n:  find the least squares solution of
   an overdetermined system, i.e., solve the least squares problem
                minimize || sub( B ) - sub( A )**H * X ||.

where sub( B ) denotes B( IB:IB+M-1, JB:JB+NRHS-1 ) when TRANS = 'N'
and B( IB:IB+N-1, JB:JB+NRHS-1 ) otherwise. Several right hand side
vectors b and solution vectors x can be handled in a single call;
When TRANS = 'N', the solution vectors are stored as the columns of
the N-by-NRHS right hand side matrix sub( B ) and the M-by-NRHS
right hand side matrix sub( B ) otherwise.
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
TRANS   (global input) CHARACTER
= 'N': the linear system involves sub( A );
.br
= 'C': the linear system involves sub( A )**H.
.TP 8
M       (global input) INTEGER
The number of rows to be operated on, i.e. the number of
rows of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on, i.e. the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
NRHS    (global input) INTEGER
The number of right hand sides, i.e. the number of columns
of the distributed submatrices sub( B ) and X.  NRHS >= 0.
.TP 8
A       (local input/local output) COMPLEX pointer into the
local memory to an array of local dimension
( LLD_A, LOCc(JA+N-1) ).  On entry, the M-by-N matrix A.
if M >= N, sub( A ) is overwritten by details of its QR
factorization as returned by PCGEQRF;
if M <  N, sub( A ) is overwritten by details of its LQ
factorization as returned by PCGELQF.
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
B       (local input/local output) COMPLEX pointer into the
local memory to an array of local dimension
(LLD_B, LOCc(JB+NRHS-1)).  On entry, this array contains the
local pieces of the distributed matrix B of right hand side
vectors, stored columnwise;
sub( B ) is M-by-NRHS if TRANS='N', and N-by-NRHS otherwise.
On exit, sub( B ) is overwritten by the solution vectors,
stored columnwise:  if TRANS = 'N' and M >= N, rows 1 to N
of sub( B ) contain the least squares solution vectors; the
residual sum of squares for the solution in each column is
given by the sum of squares of elements N+1 to M in that
column; if TRANS = 'N' and M < N, rows 1 to N of sub( B )
contain the minimum norm solution vectors; if TRANS = 'C'
and M >= N, rows 1 to M of sub( B ) contain the minimum norm
solution vectors; if TRANS = 'C' and M < N, rows 1 to M of
sub( B ) contain the least squares solution vectors; the
residual sum of squares for the solution in each column is
given by the sum of squares of elements M+1 to N in that
column.
.TP 8
IB      (global input) INTEGER
The row index in the global array B indicating the first
row of sub( B ).
.TP 8
JB      (global input) INTEGER
The column index in the global array B indicating the
first column of sub( B ).
.TP 8
DESCB   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix B.
.TP 8
WORK    (local workspace/local output) COMPLEX array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= LTAU + MAX( LWF, LWS ) where
If M >= N, then
LTAU = NUMROC( JA+MIN(M,N)-1, NB_A, MYCOL, CSRC_A, NPCOL ),
LWF  = NB_A * ( MpA0 + NqA0 + NB_A )
LWS  = MAX( (NB_A*(NB_A-1))/2, (NRHSqB0 + MpB0)*NB_A ) +
NB_A * NB_A
Else
LTAU = NUMROC( IA+MIN(M,N)-1, MB_A, MYROW, RSRC_A, NPROW ),
LWF  = MB_A * ( MpA0 + NqA0 + MB_A )
LWS  = MAX( (MB_A*(MB_A-1))/2, ( NpB0 + MAX( NqA0 +
NUMROC( NUMROC( N+IROFFB, MB_A, 0, 0, NPROW ),
MB_A, 0, 0, LCMP ), NRHSqB0 ) )*MB_A ) +
MB_A * MB_A
End if

where LCMP = LCM / NPROW with LCM = ILCM( NPROW, NPCOL ),

IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
MpA0 = NUMROC( M+IROFFA, MB_A, MYROW, IAROW, NPROW ),
NqA0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),

IROFFB = MOD( IB-1, MB_B ), ICOFFB = MOD( JB-1, NB_B ),
IBROW = INDXG2P( IB, MB_B, MYROW, RSRC_B, NPROW ),
IBCOL = INDXG2P( JB, NB_B, MYCOL, CSRC_B, NPCOL ),
MpB0 = NUMROC( M+IROFFB, MB_B, MYROW, IBROW, NPROW ),
NpB0 = NUMROC( N+IROFFB, MB_B, MYROW, IBROW, NPROW ),
NRHSqB0 = NUMROC( NRHS+ICOFFB, NB_B, MYCOL, IBCOL, NPCOL ),

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
