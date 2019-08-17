.TH PSGELQ2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSGELQ2 - compute a LQ factorization of a real distributed M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = L * Q
.SH SYNOPSIS
.TP 20
SUBROUTINE PSGELQ2(
M, N, A, IA, JA, DESCA, TAU, WORK, LWORK,
INFO )
.TP 20
.ti +4
INTEGER
IA, INFO, JA, LWORK, M, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
A( * ), TAU( * ), WORK( * )
.SH PURPOSE
PSGELQ2 computes a LQ factorization of a real distributed M-by-N matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) = L * Q. 
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
M       (global input) INTEGER
The number of rows to be operated on, i.e. the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on, i.e. the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) REAL pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
On entry, the local pieces of the M-by-N distributed matrix
sub( A ) which is to be factored. On exit, the elements on
and below the diagonal of sub( A ) contain the M by min(M,N)
lower trapezoidal matrix L (L is lower triangular if M <= N);
the elements above the diagonal, with the array TAU, repre-
sent the orthogonal matrix Q as a product of elementary
reflectors (see Further Details).
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
TAU     (local output) REAL, array, dimension
LOCr(IA+MIN(M,N)-1).  This array contains the scalar factors
of the elementary reflectors. TAU is tied to the distributed
matrix A.
.TP 8
WORK    (local workspace/local output) REAL array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= Nq0 + MAX( 1, Mp0 ), where

IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ),
Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ),

and NUMROC, INDXG2P are ScaLAPACK tool functions;
MYROW, MYCOL, NPROW and NPCOL can be determined by calling
the subroutine BLACS_GRIDINFO.

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
INFO    (local output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
.SH FURTHER DETAILS
The matrix Q is represented as a product of elementary reflectors

   Q = H(ia+k-1) H(ia+k-2) . . . H(ia), where k = min(m,n).

Each H(i) has the form
.br

   H(i) = I - tau * v * v'
.br

where tau is a real scalar, and v is a real vector with v(1:i-1)=0
and v(i) = 1; v(i+1:n) is stored on exit in A(ia+i-1,ja+i:ja+n-1),
and tau in TAU(ia+i-1).
.br

