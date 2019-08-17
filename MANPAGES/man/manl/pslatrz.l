.TH PSLATRZ l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLATRZ - reduce the M-by-N ( M<=N ) real upper trapezoidal matrix sub( A ) = [ A(IA:IA+M-1,JA:JA+M-1) A(IA:IA+M-1,JA+N-L:JA+N-1) ] to upper triangular form by means of orthogonal transformations
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLATRZ(
M, N, L, A, IA, JA, DESCA, TAU, WORK )
.TP 20
.ti +4
INTEGER
IA, JA, L, M, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
A( * ), TAU( * ), WORK( * )
.SH PURPOSE
PSLATRZ reduces the M-by-N ( M<=N ) real upper trapezoidal matrix sub( A ) = [ A(IA:IA+M-1,JA:JA+M-1) A(IA:IA+M-1,JA+N-L:JA+N-1) ] to upper triangular form by means of orthogonal transformations. 
The upper trapezoidal matrix sub( A ) is factored as
.br

   sub( A ) = ( R  0 ) * Z,
.br

where Z is an N-by-N orthogonal matrix and R is an M-by-M upper
triangular matrix.
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
M       (global input) INTEGER
The number of rows to be operated on, i.e. the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on, i.e. the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
L       (global input) INTEGER
The columns of the distributed submatrix sub( A ) containing
the meaningful part of the Householder reflectors. L > 0.
.TP 8
A       (local input/local output) REAL pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
On entry, the local pieces of the M-by-N distributed matrix
sub( A ) which is to be factored. On exit, the leading M-by-M
upper triangular part of sub( A ) contains the upper trian-
gular matrix R, and elements N-L+1 to N of the first M rows
of sub( A ), with the array TAU, represent the orthogonal
matrix Z as a product of M elementary reflectors.
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
TAU     (local output) REAL, array, dimension LOCr(IA+M-1)
This array contains the scalar factors of the elementary
reflectors. TAU is tied to the distributed matrix A.
.TP 8
WORK    (local workspace) REAL array, dimension (LWORK)
LWORK >= Nq0 + MAX( 1, Mp0 ), where

IROFF = MOD( IA-1, MB_A ), ICOFF = MOD( JA-1, NB_A ),
IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
Mp0   = NUMROC( M+IROFF, MB_A, MYROW, IAROW, NPROW ),
Nq0   = NUMROC( N+ICOFF, NB_A, MYCOL, IACOL, NPCOL ),

and NUMROC, INDXG2P are ScaLAPACK tool functions;
MYROW, MYCOL, NPROW and NPCOL can be determined by calling
the subroutine BLACS_GRIDINFO.
.SH FURTHER DETAILS
The  factorization is obtained by Householder's method.  The kth
transformation matrix, Z( k ), which is used to introduce zeros into
the (m - k + 1)th row of sub( A ), is given in the form
.br

   Z( k ) = ( I     0   ),
.br
            ( 0  T( k ) )
.br

where
.br

   T( k ) = I - tau*u( k )*u( k )',   u( k ) = (   1    ),
                                               (   0    )
                                               ( z( k ) )

tau is a scalar and z( k ) is an ( n - m ) element vector.
tau and z( k ) are chosen to annihilate the elements of the kth row
of sub( A ).
.br

The scalar tau is returned in the kth element of TAU and the vector
u( k ) in the kth row of sub( A ), such that the elements of z( k )
are in  a( k, m + 1 ), ..., a( k, n ). The elements of R are returned
in the upper triangular part of sub( A ).
.br

Z is given by
.br

   Z =  Z( 1 ) * Z( 2 ) * ... * Z( m ).
.br

