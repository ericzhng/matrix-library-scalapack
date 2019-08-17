.SH NAME
ZPTTRSV - solve one of the triangular systems L * X = B, or L**H * X = B,
.SH SYNOPSIS
.TP 20
SUBROUTINE ZPTTRSV(
UPLO, TRANS, N, NRHS, D, E, B, LDB,
INFO )
.TP 20
.ti +4
CHARACTER
UPLO, TRANS
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
COMPLEX*16
B( LDB, * ), E( * )
.SH PURPOSE
ZPTTRSV solves one of the triangular systems L * X = B, or L**H * X = B,    U * X = B, or  U**H * X = B,
.br
where L or U is the Cholesky factor of a Hermitian positive
definite tridiagonal matrix A such that
.br
A = U**H*D*U or A = L*D*L**H (computed by ZPTTRF).
.br

.SH ARGUMENTS
.TP 8
UPLO    (input) CHARACTER*1
Specifies whether the superdiagonal or the subdiagonal
of the tridiagonal matrix A is stored and the form of the
factorization:
.br
= 'U':  E is the superdiagonal of U, and A = U'*D*U;
.br
= 'L':  E is the subdiagonal of L, and A = L*D*L'.
(The two forms are equivalent if A is real.)
.TP 8
TRANS   (input) CHARACTER
Specifies the form of the system of equations:
.br
= 'N':  L * X = B     (No transpose)
.br
= 'N':  L * X = B     (No transpose)
.br
= 'C':  U**H * X = B  (Conjugate transpose)
.br
= 'C':  L**H * X = B  (Conjugate transpose)
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
factorization computed by ZPTTRF.
.TP 8
E       (input) COMPLEX array, dimension (N-1)
The (n-1) off-diagonal elements of the unit bidiagonal
factor U or L from the factorization computed by ZPTTRF
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
