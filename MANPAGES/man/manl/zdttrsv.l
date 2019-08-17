.SH NAME
ZDTTRSV - solve one of the systems of equations L * X = B, L**T * X = B, or L**H * X = B,
.SH SYNOPSIS
.TP 20
SUBROUTINE ZDTTRSV(
UPLO, TRANS, N, NRHS, DL, D, DU,
B, LDB, INFO )
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
COMPLEX*16
B( LDB, * ), D( * ), DL( * ), DU( * )
.SH PURPOSE
ZDTTRSV solves one of the systems of equations L * X = B, L**T * X = B, or L**H * X = B,    U * X = B,  U**T * X = B,  or  U**H * X = B,
.br
with factors of the tridiagonal matrix A from the LU factorization
computed by ZDTTRF.
.br

.SH ARGUMENTS
.TP 8
UPLO    (input) CHARACTER*1
Specifies whether to solve with L or U.
.TP 8
TRANS   (input) CHARACTER
Specifies the form of the system of equations:
.br
= 'N':  A * X = B     (No transpose)
.br
= 'T':  A**T * X = B  (Transpose)
.br
= 'C':  A**H * X = B  (Conjugate transpose)
.TP 8
N       (input) INTEGER
The order of the matrix A.  N >= 0.
.TP 8
NRHS    (input) INTEGER
The number of right hand sides, i.e., the number of columns
of the matrix B.  NRHS >= 0.
.TP 8
DL      (input) COMPLEX array, dimension (N-1)
The (n-1) multipliers that define the matrix L from the
LU factorization of A.
.TP 8
D       (input) COMPLEX array, dimension (N)
The n diagonal elements of the upper triangular matrix U from
the LU factorization of A.
.TP 8
DU      (input) COMPLEX array, dimension (N-1)
The (n-1) elements of the first superdiagonal of U.
.TP 8
B       (input/output) COMPLEX array, dimension (LDB,NRHS)
On entry, the right hand side matrix B.
On exit, B is overwritten by the solution matrix X.
.TP 8
LDB     (input) INTEGER
The leading dimension of the array B.  LDB >= max(1,N).
.TP 8
INFO    (output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
