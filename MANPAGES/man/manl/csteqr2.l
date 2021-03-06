.TH CSTEQR2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
CSTEQR2 - i a modified version of LAPACK routine CSTEQR
.SH SYNOPSIS
.TP 20
SUBROUTINE CSTEQR2(
COMPZ, N, D, E, Z, LDZ, NR, WORK, INFO )
.TP 20
.ti +4
CHARACTER
COMPZ
.TP 20
.ti +4
INTEGER
INFO, LDZ, N, NR
.TP 20
.ti +4
REAL
D( * ), E( * ), WORK( * )
.TP 20
.ti +4
COMPLEX
Z( LDZ, * )
.SH PURPOSE
CSTEQR2 is a modified version of LAPACK routine CSTEQR. CSTEQR2 computes all eigenvalues and, optionally, eigenvectors of a
symmetric tridiagonal matrix using the implicit QL or QR method.
CSTEQR2 is modified from CSTEQR to allow each ScaLAPACK process
running CSTEQR2 to perform updates on a distributed matrix Q.
Proper usage of CSTEQR2 can be gleaned from
.br
examination of ScaLAPACK's *  PCHEEV.
.br
CSTEQR2 incorporates changes attributed to Greg Henry.
.br

.SH ARGUMENTS
.TP 8
COMPZ   (input) CHARACTER*1
= 'N':  Compute eigenvalues only.
.br
= 'I':  Compute eigenvalues and eigenvectors of the
tridiagonal matrix.  Z must be initialized to the
identity matrix by PCLASET or CLASET prior
to entering this subroutine.
.TP 8
N       (input) INTEGER
The order of the matrix.  N >= 0.
.TP 8
D       (input/output) REAL array, dimension (N)
On entry, the diagonal elements of the tridiagonal matrix.
On exit, if INFO = 0, the eigenvalues in ascending order.
.TP 8
E       (input/output) REAL array, dimension (N-1)
On entry, the (n-1) subdiagonal elements of the tridiagonal
matrix.
On exit, E has been destroyed.
.TP 8
Z       (local input/local output) COMPLEX array, global
dimension (N, N), local dimension (LDZ, NR).
On entry, if  COMPZ = 'V', then Z contains the orthogonal
matrix used in the reduction to tridiagonal form.
On exit, if INFO = 0, then if  COMPZ = 'V', Z contains the
orthonormal eigenvectors of the original symmetric matrix,
and if COMPZ = 'I', Z contains the orthonormal eigenvectors
of the symmetric tridiagonal matrix.
If COMPZ = 'N', then Z is not referenced.
.TP 8
LDZ     (input) INTEGER
The leading dimension of the array Z.  LDZ >= 1, and if
eigenvectors are desired, then  LDZ >= max(1,N).
.TP 8
NR      (input) INTEGER
NR = MAX(1, NUMROC( N, NB, MYPROW, 0, NPROCS ) ).
If COMPZ = 'N', then NR is not referenced.
.TP 8
WORK    (workspace) REAL array, dimension (max(1,2*N-2))
If COMPZ = 'N', then WORK is not referenced.
.TP 8
INFO    (output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
.br
> 0:  the algorithm has failed to find all the eigenvalues in
a total of 30*N iterations; if INFO = i, then i
elements of E have not converged to zero; on exit, D
and E contain the elements of a symmetric tridiagonal
matrix which is orthogonally similar to the original
matrix.
