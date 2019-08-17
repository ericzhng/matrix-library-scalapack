.TH PDSYTD2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDSYTD2 - reduce a real symmetric matrix sub( A ) to symmetric tridiagonal form T by an orthogonal similarity transformation
.SH SYNOPSIS
.TP 20
SUBROUTINE PDSYTD2(
UPLO, N, A, IA, JA, DESCA, D, E, TAU, WORK,
LWORK, INFO )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
IA, INFO, JA, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), D( * ), E( * ), TAU( * ), WORK( * )
.SH PURPOSE
PDSYTD2 reduces a real symmetric matrix sub( A ) to symmetric tridiagonal form T by an orthogonal similarity transformation: Q' * sub( A ) * Q = T, where sub( A ) = A(IA:IA+N-1,JA:JA+N-1).

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
UPLO    (global input) CHARACTER
Specifies whether the upper or lower triangular part of the
symmetric matrix sub( A ) is stored:
.br
= 'U':  Upper triangular
.br
= 'L':  Lower triangular
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
On entry, this array contains the local pieces of the
symmetric distributed matrix sub( A ).  If UPLO = 'U', the
leading N-by-N upper triangular part of sub( A ) contains
the upper triangular part of the matrix, and its strictly
lower triangular part is not referenced. If UPLO = 'L', the
leading N-by-N lower triangular part of sub( A ) contains the
lower triangular part of the matrix, and its strictly upper
triangular part is not referenced. On exit, if UPLO = 'U',
the diagonal and first superdiagonal of sub( A ) are over-
written by the corresponding elements of the tridiagonal
matrix T, and the elements above the first superdiagonal,
with the array TAU, represent the orthogonal matrix Q as a
product of elementary reflectors; if UPLO = 'L', the diagonal
and first subdiagonal of sub( A ) are overwritten by the
corresponding elements of the tridiagonal matrix T, and the
elements below the first subdiagonal, with the array TAU,
represent the orthogonal matrix Q as a product of elementary
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
D       (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
The diagonal elements of the tridiagonal matrix T:
D(i) = A(i,i). D is tied to the distributed matrix A.
.TP 8
E       (local output) DOUBLE PRECISION array, dimension LOCc(JA+N-1)
if UPLO = 'U', LOCc(JA+N-2) otherwise. The off-diagonal
elements of the tridiagonal matrix T: E(i) = A(i,i+1) if
UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. E is tied to the
distributed matrix A.
.TP 8
TAU     (local output) DOUBLE PRECISION array, dimension
LOCc(JA+N-1). This array contains the scalar factors TAU of
the elementary reflectors. TAU is tied to the distributed
matrix A.
.TP 8
WORK    (local workspace/local output) DOUBLE PRECISION array,
dimension (LWORK)
On exit, WORK( 1 ) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= 3*N.

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
If UPLO = 'U', the matrix Q is represented as a product of elementary
reflectors
.br

   Q = H(n-1) . . . H(2) H(1).
.br

Each H(i) has the form
.br

   H(i) = I - tau * v * v'
.br

where tau is a real scalar, and v is a real vector with
.br
v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
.br
A(ia:ia+i-2,ja+i), and tau in TAU(ja+i-1).
.br

If UPLO = 'L', the matrix Q is represented as a product of elementary
reflectors
.br

   Q = H(1) H(2) . . . H(n-1).
.br

Each H(i) has the form
.br

   H(i) = I - tau * v * v'
.br

where tau is a real scalar, and v is a real vector with
.br
v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
.br
A(ia+i+1:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
.br

The contents of sub( A ) on exit are illustrated by the following
examples with n = 5:
.br

if UPLO = 'U':                       if UPLO = 'L':
.br

  (  d   e   v2  v3  v4 )              (  d                  )
  (      d   e   v3  v4 )              (  e   d              )
  (          d   e   v4 )              (  v1  e   d          )
  (              d   e  )              (  v1  v2  e   d      )
  (                  d  )              (  v1  v2  v3  e   d  )

where d and e denote diagonal and off-diagonal elements of T, and vi
denotes an element of the vector defining H(i).
.br

Alignment requirements
.br
======================
.br

The distributed submatrix sub( A ) must verify some alignment proper-
ties, namely the following expression should be true:
.br
( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA ) with
.br
IROFFA = MOD( IA-1, MB_A ) and ICOFFA = MOD( JA-1, NB_A ).

