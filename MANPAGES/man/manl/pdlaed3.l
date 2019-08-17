.TH PDLAED3 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLAED3 - find the roots of the secular equation, as defined by the values in D, W, and RHO, between 1 and K
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLAED3(
ICTXT, K, N, NB, D, DROW, DCOL, RHO, DLAMDA,
W, Z, U, LDU, BUF, INDX, INDCOL, INDROW,
INDXR, INDXC, CTOT, NPCOL, INFO )
.TP 20
.ti +4
INTEGER
DCOL, DROW, ICTXT, INFO, K, LDU, N, NB, NPCOL
.TP 20
.ti +4
DOUBLE
PRECISION RHO
.TP 20
.ti +4
INTEGER
CTOT( 0: NPCOL-1, 4 ), INDCOL( * ),
INDROW( * ), INDX( * ), INDXC( * ), INDXR( * )
.TP 20
.ti +4
DOUBLE
PRECISION BUF( * ), D( * ), DLAMDA( * ), U( LDU, * ),
W( * ), Z( * )
.SH PURPOSE
PDLAED3 finds the roots of the secular equation, as defined by the values in D, W, and RHO, between 1 and K. It makes the appropriate calls to SLAED4
.br

This code makes very mild assumptions about floating point
arithmetic. It will work on machines with a guard digit in
add/subtract, or on those binary machines without guard digits
which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or Cray-2.
It could conceivably fail on hexadecimal or decimal machines
without guard digits, but we know of none.
.br

.SH ARGUMENTS
.TP 7
ICTXT  (global input) INTEGER
The BLACS context handle, indicating the global context of
the operation on the matrix. The context itself is global.
.TP 7
K      (output) INTEGER
The number of non-deflated eigenvalues, and the order of the
related secular equation. 0 <= K <=N.
.TP 7
N      (input) INTEGER
The dimension of the symmetric tridiagonal matrix.  N >= 0.
.TP 8
NB      (global input) INTEGER
The blocking factor used to distribute the columns of the
matrix. NB >= 1.
.TP 7
D      (input/output) DOUBLE PRECISION array, dimension (N)
On entry, D contains the eigenvalues of the two submatrices to
be combined.
On exit, D contains the trailing (N-K) updated eigenvalues
(those which were deflated) sorted into increasing order.
.TP 7
DROW   (global input) INTEGER
The process row over which the first row of the matrix D is
distributed. 0 <= DROW < NPROW.
.TP 7
DCOL   (global input) INTEGER
The process column over which the first column of the
matrix D is distributed. 0 <= DCOL < NPCOL.
.TP 7
Q      (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
On entry, Q contains the eigenvectors of two submatrices in
the two square blocks with corners at (1,1), (N1,N1)
and (N1+1, N1+1), (N,N).
On exit, Q contains the trailing (N-K) updated eigenvectors
(those which were deflated) in its last N-K columns.
.TP 7
LDQ    (input) INTEGER
The leading dimension of the array Q.  LDQ >= max(1,NQ).
.TP 7
RHO    (global input/output) DOUBLE PRECISION
On entry, the off-diagonal element associated with the rank-1
cut which originally split the two submatrices which are now
being recombined.
On exit, RHO has been modified to the value required by
PDLAED3.

DLAMDA (global output) DOUBLE PRECISION array, dimension (N)
A copy of the first K eigenvalues which will be used by
SLAED3 to form the secular equation.
.TP 7
W      (global output) DOUBLE PRECISION array, dimension (N)
The first k values of the final deflation-altered z-vector
which will be passed to SLAED3.
.TP 7
Z      (global input) DOUBLE PRECISION array, dimension (N)
On entry, Z contains the updating vector (the last
row of the first sub-eigenvector matrix and the first row of
the second sub-eigenvector matrix).
On exit, the contents of Z have been destroyed by the updating
process.
.TP 6
U     (global output) DOUBLE PRECISION array
global dimension (N, N), local dimension (LDU, NQ).
Q  contains the orthonormal eigenvectors of the symmetric
tridiagonal matrix.
.TP 7
LDU    (input) INTEGER
The leading dimension of the array U.
.TP 7
QBUF   (workspace) DOUBLE PRECISION array, dimension 3*N
.TP 7
INDX   (workspace) INTEGER array, dimension (N)
The permutation used to sort the contents of DLAMDA into
ascending order.

INDCOL (workspace) INTEGER array, dimension (N)

INDROW (workspace) INTEGER array, dimension (N)

INDXR (workspace) INTEGER array, dimension (N)

INDXC (workspace) INTEGER array, dimension (N)
.TP 7
CTOT   (workspace) INTEGER array, dimension( NPCOL, 4)
.TP 8
NPCOL   (global input) INTEGER
The total number of columns over which the distributed
submatrix is distributed.
.TP 7
INFO   (output) INTEGER
= 0:  successful exit.
.br
< 0:  if INFO = -i, the i-th argument had an illegal value.
.br
> 0:  The algorithm failed to compute the ith eigenvalue.
