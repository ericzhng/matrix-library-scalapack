.TH PZLATRD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZLATRD - reduce NB rows and columns of a complex Hermitian distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) to complex tridiagonal form by an unitary similarity transformation Q' * sub( A ) * Q, and returns the matrices V and W which are needed to apply the transformation to the unreduced part of sub( A )
.SH SYNOPSIS
.TP 20
SUBROUTINE PZLATRD(
UPLO, N, NB, A, IA, JA, DESCA, D, E, TAU, W,
IW, JW, DESCW, WORK )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
IA, IW, JA, JW, N, NB
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCW( * )
.TP 20
.ti +4
DOUBLE
PRECISION D( * ), E( * )
.TP 20
.ti +4
COMPLEX*16
A( * ), TAU( * ), W( * ), WORK( * )
.SH PURPOSE
PZLATRD reduces NB rows and columns of a complex Hermitian distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) to complex tridiagonal form by an unitary similarity transformation Q' * sub( A ) * Q, and returns the matrices V and W which are needed to apply the transformation to the unreduced part of sub( A ). 
If UPLO = 'U', PZLATRD reduces the last NB rows and columns of a
matrix, of which the upper triangle is supplied;
.br
if UPLO = 'L', PZLATRD reduces the first NB rows and columns of a
matrix, of which the lower triangle is supplied.
.br

This is an auxiliary routine called by PZHETRD.
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
UPLO    (global input) CHARACTER
Specifies whether the upper or lower triangular part of the
Hermitian matrix sub( A ) is stored:
.br
= 'U': Upper triangular
.br
= 'L': Lower triangular
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
NB      (global input) INTEGER
The number of rows and columns to be reduced.
.TP 8
A       (local input/local output) COMPLEX*16 pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
On entry, this array contains the local pieces of the
Hermitian distributed matrix sub( A ).  If UPLO = 'U', the
leading N-by-N upper triangular part of sub( A ) contains
the upper triangular part of the matrix, and its strictly
lower triangular part is not referenced. If UPLO = 'L', the
leading N-by-N lower triangular part of sub( A ) contains the
lower triangular part of the matrix, and its strictly upper
triangular part is not referenced.
On exit, if UPLO = 'U', the last NB columns have been reduced
to tridiagonal form, with the diagonal elements overwriting
the diagonal elements of sub( A ); the elements above the
diagonal with the array TAU, represent the unitary matrix Q
as a product of elementary reflectors. If UPLO = 'L', the
first NB columns have been reduced to tridiagonal form, with
the diagonal elements overwriting the diagonal elements of
sub( A ); the elements below the diagonal with the array TAU,
represent the unitary matrix Q as a product of elementary
reflectors; See Further Details.
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
TAU     (local output) COMPLEX*16, array, dimension
LOCc(JA+N-1). This array contains the scalar factors TAU of
the elementary reflectors. TAU is tied to the distributed
matrix A.
.TP 8
W       (local output) COMPLEX*16 pointer into the local memory
to an array of dimension (LLD_W,NB_W), This array contains
the local pieces of the N-by-NB_W matrix W required to
update the unreduced part of sub( A ).
.TP 8
IW      (global input) INTEGER
The row index in the global array W indicating the first
row of sub( W ).
.TP 8
JW      (global input) INTEGER
The column index in the global array W indicating the
first column of sub( W ).
.TP 8
DESCW   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix W.
.TP 8
WORK    (local workspace) COMPLEX*16 array, dimension (NB_A)
.SH FURTHER DETAILS
If UPLO = 'U', the matrix Q is represented as a product of elementary
reflectors
.br

   Q = H(n) H(n-1) . . . H(n-nb+1).
.br

Each H(i) has the form
.br

   H(i) = I - tau * v * v'
.br

where tau is a complex scalar, and v is a complex vector with
v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in
.br
A(ia:ia+i-2,ja+i), and tau in TAU(ja+i-1).
.br

If UPLO = 'L', the matrix Q is represented as a product of elementary
reflectors
.br

   Q = H(1) H(2) . . . H(nb).
.br

Each H(i) has the form
.br

   H(i) = I - tau * v * v'
.br

where tau is a complex scalar, and v is a complex vector with
v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in
.br
A(ia+i+1:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
.br

The elements of the vectors v together form the N-by-NB matrix V
which is needed, with W, to apply the transformation to the unreduced
part of the matrix, using a Hermitian rank-2k update of the form:
sub( A ) := sub( A ) - V*W' - W*V'.
.br

The contents of A on exit are illustrated by the following examples
with n = 5 and nb = 2:
.br

if UPLO = 'U':                       if UPLO = 'L':
.br

  (  a   a   a   v4  v5 )              (  d                  )
  (      a   a   v4  v5 )              (  1   d              )
  (          a   1   v5 )              (  v1  1   a          )
  (              d   1  )              (  v1  v2  a   a      )
  (                  d  )              (  v1  v2  a   a   a  )

where d denotes a diagonal element of the reduced matrix, a denotes
an element of the original matrix that is unchanged, and vi denotes
an element of the vector defining H(i).
.br

