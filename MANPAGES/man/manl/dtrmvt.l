.TH DTRMVT l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
DTRMVT - perform the matrix-vector operations  x := T' *y, and w := T *z,
.SH SYNOPSIS
.TP 19
SUBROUTINE DTRMVT(
UPLO, N, T, LDT, X, INCX, Y, INCY, W, INCW, Z,
INCZ )
.TP 19
.ti +4
CHARACTER
UPLO
.TP 19
.ti +4
INTEGER
INCW, INCX, INCY, INCZ, LDT, N
.TP 19
.ti +4
DOUBLE
PRECISION T( LDT, * ), W( * ), X( * ), Y( * ), Z( * )
.SH PURPOSE
DTRMVT performs the matrix-vector operations x := T' *y, and w := T *z, 
where x is an n element vector and  T is an n by n
.br
upper or lower triangular matrix.
.br

.SH ARGUMENTS
.TP 7
UPLO   - CHARACTER*1.
On entry, UPLO specifies whether the matrix is an upper or
lower triangular matrix as follows:

UPLO = 'U' or 'u'   A is an upper triangular matrix.

UPLO = 'L' or 'l'   A is a lower triangular matrix.

Unchanged on exit.
.TP 7
N      - INTEGER.
On entry, N specifies the order of the matrix A.
N must be at least zero.
Unchanged on exit.
.TP 7
T      - DOUBLE PRECISION array of DIMENSION ( LDT, n ).
Before entry with  UPLO = 'U' or 'u', the leading n by n
upper triangular part of the array T must contain the upper
triangular matrix and the strictly lower triangular part of
T is not referenced.
Before entry with UPLO = 'L' or 'l', the leading n by n
lower triangular part of the array T must contain the lower
triangular matrix and the strictly upper triangular part of
T is not referenced.
.TP 7
LDT    - INTEGER.
On entry, LDA specifies the first dimension of A as declared
in the calling (sub) program. LDA must be at least
max( 1, n ).
Unchanged on exit.
.TP 7
X      - DOUBLE PRECISION array of dimension at least
( 1 + ( n - 1 )*abs( INCX ) ).
On exit, X = T' * y
.TP 7
INCX   - INTEGER.
On entry, INCX specifies the increment for the elements of
X. INCX must not be zero.
Unchanged on exit.
.TP 7
Y      - DOUBLE PRECISION array of dimension at least
( 1 + ( n - 1 )*abs( INCY ) ).
Before entry, the incremented array Y must contain the n
element vector y.  Unchanged on exit.
.TP 7
INCY   - INTEGER.
On entry, INCY specifies the increment for the elements of
Y. INCY must not be zero.
Unchanged on exit.
.TP 7
W      - DOUBLE PRECISION array of dimension at least
( 1 + ( n - 1 )*abs( INCW ) ).
On exit, W = T * z
.TP 7
INCW   - INTEGER.
On entry, INCW specifies the increment for the elements of
W. INCW must not be zero.
Unchanged on exit.
.TP 7
Z      - DOUBLE PRECISION array of dimension at least
( 1 + ( n - 1 )*abs( INCZ ) ).
Before entry, the incremented array Z must contain the n
element vector z.  Unchanged on exit.
.TP 7
INCY   - INTEGER.
On entrz, INCY specifies the increment for the elements of
Y. INCY must not be zero.
Unchanged on exit.

Level 2 Blas routine.

.. Local Scalars ..

..
.. External Functions ..

..
.. External Subroutines ..

..
.. Intrinsic Functions ..

..
.. Executable Statements ..

Test the input parameters.

Quick return if possible.

End of DTRMVT.
