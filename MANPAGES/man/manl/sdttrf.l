.SH NAME
SDTTRF - compute an LU factorization of a complex tridiagonal matrix A using elimination without partial pivoting
.SH SYNOPSIS
.TP 19
SUBROUTINE SDTTRF(
N, DL, D, DU, INFO )
.TP 19
.ti +4
INTEGER
INFO, N
.TP 19
.ti +4
REAL
D( * ), DL( * ), DU( * )
.SH PURPOSE
SDTTRF computes an LU factorization of a complex tridiagonal matrix A using elimination without partial pivoting. 
The factorization has the form
.br
   A = L * U
.br
where L is a product of unit lower bidiagonal
.br
matrices and U is upper triangular with nonzeros in only the main
diagonal and first superdiagonal.
.br

.SH ARGUMENTS
.TP 8
N       (input) INTEGER
The order of the matrix A.  N >= 0.
.TP 8
DL      (input/output) COMPLEX array, dimension (N-1)
On entry, DL must contain the (n-1) subdiagonal elements of
A.
On exit, DL is overwritten by the (n-1) multipliers that
define the matrix L from the LU factorization of A.
.TP 8
D       (input/output) COMPLEX array, dimension (N)
On entry, D must contain the diagonal elements of A.
On exit, D is overwritten by the n diagonal elements of the
upper triangular matrix U from the LU factorization of A.
.TP 8
DU      (input/output) COMPLEX array, dimension (N-1)
On entry, DU must contain the (n-1) superdiagonal elements
of A.
On exit, DU is overwritten by the (n-1) elements of the first
superdiagonal of U.
.TP 8
INFO    (output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
.br
> 0:  if INFO = i, U(i,i) is exactly zero. The factorization
has been completed, but the factor U is exactly
singular, and division by zero will occur if it is used
to solve a system of equations.
