.SH NAME
DPTTRSV - solve one of the triangular systems L**T* X = B, or L * X = B,
.SH SYNOPSIS
.TP 20
SUBROUTINE DPTTRSV(
TRANS, N, NRHS, D, E, B, LDB,
INFO )
.TP 20
.ti +4
CHARACTER
TRANS
.TP 20
.ti +4
INTEGER
INFO, LDB, N, NRHS
.TP 20
.ti +4
DOUBLE
PRECISION D( * )
.TP 20
.ti +4
DOUBLE
PRECISION B( LDB, * ), E( * )
.SH PURPOSE
DPTTRSV solves one of the triangular systems L**T* X = B, or L * X = B, where L is the Cholesky factor of a Hermitian positive
.br
definite tridiagonal matrix A such that
.br
A = L*D*L**H (computed by DPTTRF).
.br

.SH ARGUMENTS
.TP 8
TRANS   (input) CHARACTER
Specifies the form of the system of equations:
.br
= 'N':  L * X = B     (No transpose)
.br
= 'T':  L**T * X = B  (Transpose)
.TP 8
N       (input) INTEGER
The order of the tridiagonal matrix A.  N >= 0.
.TP 8
NRHS    (input) INTEGER
The number of right hand sides, i.e., the number of columns
of the matrix B.  NRHS >= 0.
.TP 8
D       (input) REAL array, dimension (N)
The n diagonal elements of the diagonal matrix D from the
factorization computed by DPTTRF.
.TP 8
E       (input) COMPLEX array, dimension (N-1)
The (n-1) off-diagonal elements of the unit bidiagonal
factor U or L from the factorization computed by DPTTRF
(see UPLO).
.TP 8
B       (input/output) COMPLEX array, dimension (LDB,NRHS)
On entry, the right hand side matrix B.
On exit, the solution matrix X.
.TP 8
LDB     (input) INTEGER
The leading dimension of the array B.  LDB >= max(1,N).
.TP 8
INFO    (output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
