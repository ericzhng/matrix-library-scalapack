.TH PDGGQRF l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDGGQRF - compute a generalized QR factorization of an N-by-M matrix sub( A ) = A(IA:IA+N-1,JA:JA+M-1) and an N-by-P matrix sub( B ) = B(IB:IB+N-1,JB:JB+P-1)
.SH SYNOPSIS
.TP 20
SUBROUTINE PDGGQRF(
N, M, P, A, IA, JA, DESCA, TAUA, B, IB, JB,
DESCB, TAUB, WORK, LWORK, INFO )
.TP 20
.ti +4
INTEGER
IA, IB, INFO, JA, JB, LWORK, M, N, P
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), B( * ), TAUA( * ), TAUB( * ), WORK( * )
.SH PURPOSE
PDGGQRF computes a generalized QR factorization of an N-by-M matrix sub( A ) = A(IA:IA+N-1,JA:JA+M-1) and an N-by-P matrix sub( B ) = B(IB:IB+N-1,JB:JB+P-1): 
            sub( A ) = Q*R,        sub( B ) = Q*T*Z,
.br

where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal
matrix, and R and T assume one of the forms:
.br

if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N,
                (  0  ) N-M                         N   M-N
                   M
.br

where R11 is upper triangular, and
.br

if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P,
                 P-N  N                           ( T21 ) P
                                                     P
.br

where T12 or T21 is upper triangular.
.br

In particular, if sub( B ) is square and nonsingular, the GQR
factorization of sub( A ) and sub( B ) implicitly gives the QR
factorization of inv( sub( B ) )* sub( A ):
.br

             inv( sub( B ) )*sub( A )= Z'*(inv(T)*R)
.br

where inv( sub( B ) ) denotes the inverse of the matrix sub( B ),
and Z' denotes the transpose of matrix Z.
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
The number of rows to be operated on i.e the number of rows
of the distributed submatrices sub( A ) and sub( B ). N >= 0.
.TP 8
M       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( A ).  M >= 0.
.TP 8
P       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( B ).  P >= 0.
.TP 8
A       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+M-1)).
On entry, the local pieces of the N-by-M distributed matrix
sub( A ) which is to be factored.  On exit, the elements on
and above the diagonal of sub( A ) contain the min(N,M) by M
upper trapezoidal matrix R (R is upper triangular if N >= M);
the elements below the diagonal, with the array TAUA,
represent the orthogonal matrix Q as a product of min(N,M)
elementary reflectors (see Further Details).
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
TAUA    (local output) DOUBLE PRECISION array, dimension
LOCc(JA+MIN(N,M)-1). This array contains the scalar factors
TAUA of the elementary reflectors which represent the
orthogonal matrix Q. TAUA is tied to the distributed matrix
A. (see Further Details).
B       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_B, LOCc(JB+P-1)).
On entry, the local pieces of the N-by-P distributed matrix
sub( B ) which is to be factored. On exit, if N <= P, the
upper triangle of B(IB:IB+N-1,JB+P-N:JB+P-1) contains the
N by N upper triangular matrix T; if N > P, the elements on
and above the (N-P)-th subdiagonal contain the N by P upper
trapezoidal matrix T; the remaining elements, with the array
TAUB, represent the orthogonal matrix Z as a product of
elementary reflectors (see Further Details).
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
TAUB    (local output) DOUBLE PRECISION array, dimension LOCr(IB+N-1)
This array contains the scalar factors of the elementary
reflectors which represent the orthogonal unitary matrix Z.
TAUB is tied to the distributed matrix B (see Further
Details).
.TP 8
WORK    (local workspace/local output) DOUBLE PRECISION array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= MAX( NB_A * ( NpA0 + MqA0 + NB_A ),
MAX( (NB_A*(NB_A-1))/2, (PqB0 + NpB0)*NB_A ) +
NB_A * NB_A,
MB_B * ( NpB0 + PqB0 + MB_B ) ), where

IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
IAROW  = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
IACOL  = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
NpA0   = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),
MqA0   = NUMROC( M+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),

IROFFB = MOD( IB-1, MB_B ), ICOFFB = MOD( JB-1, NB_B ),
IBROW  = INDXG2P( IB, MB_B, MYROW, RSRC_B, NPROW ),
IBCOL  = INDXG2P( JB, NB_B, MYCOL, CSRC_B, NPCOL ),
NpB0   = NUMROC( N+IROFFB, MB_B, MYROW, IBROW, NPROW ),
PqB0   = NUMROC( P+ICOFFB, NB_B, MYCOL, IBCOL, NPCOL ),

and NUMROC, INDXG2P are ScaLAPACK tool functions;
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
The matrix Q is represented as a product of elementary reflectors

   Q = H(ja) H(ja+1) . . . H(ja+k-1), where k = min(n,m).

Each H(i) has the form
.br

   H(i) = I - taua * v * v'
.br

where taua is a real scalar, and v is a real vector with
.br
v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in
.br
A(ia+i:ia+n-1,ja+i-1), and taua in TAUA(ja+i-1).
.br
To form Q explicitly, use ScaLAPACK subroutine PDORGQR.
.br
To use Q to update another matrix, use ScaLAPACK subroutine PDORMQR.

The matrix Z is represented as a product of elementary reflectors

   Z = H(ib) H(ib+1) . . . H(ib+k-1), where k = min(n,p).

Each H(i) has the form
.br

   H(i) = I - taub * v * v'
.br

where taub is a real scalar, and v is a real vector with
.br
v(p-k+i+1:p) = 0 and v(p-k+i) = 1; v(1:p-k+i-1) is stored on exit in
B(ib+n-k+i-1,jb:jb+p-k+i-2), and taub in TAUB(ib+n-k+i-1).
To form Z explicitly, use ScaLAPACK subroutine PDORGRQ.
.br
To use Z to update another matrix, use ScaLAPACK subroutine PDORMRQ.

Alignment requirements
.br
======================
.br

The distributed submatrices sub( A ) and sub( B ) must verify some
alignment properties, namely the following expression should be true:

( MB_A.EQ.MB_B .AND. IROFFA.EQ.IROFFB .AND. IAROW.EQ.IBROW )

