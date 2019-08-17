.TH PSLAED2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAED2 - sort the two sets of eigenvalues together into a single sorted set
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLAED2(
ICTXT, K, N, N1, NB, D, DROW, DCOL, Q, LDQ,
RHO, Z, W, DLAMDA, Q2, LDQ2, QBUF, CTOT, PSM,
NPCOL, INDX, INDXC, INDXP, INDCOL, COLTYP, NN,
NN1, NN2, IB1, IB2 )
.TP 20
.ti +4
INTEGER
DCOL, DROW, IB1, IB2, ICTXT, K, LDQ, LDQ2, N,
N1, NB, NN, NN1, NN2, NPCOL
.TP 20
.ti +4
REAL
RHO
.TP 20
.ti +4
INTEGER
COLTYP( * ), CTOT( 0: NPCOL-1, 4 ),
INDCOL( N ), INDX( * ), INDXC( * ), INDXP( * ),
PSM( 0: NPCOL-1, 4 )
.TP 20
.ti +4
REAL
D( * ), DLAMDA( * ), Q( LDQ, * ),
Q2( LDQ2, * ), QBUF( * ), W( * ), Z( * )
.SH PURPOSE
PSLAED2 sorts the two sets of eigenvalues together into a single sorted set. Then it tries to deflate the size of the problem. There are two ways in which deflation can occur:  when two or more
eigenvalues are close together or if there is a tiny entry in the
Z vector.  For each such occurrence the order of the related secular
equation problem is reduced by one.
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
.TP 7
N1     (input) INTEGER
The location of the last eigenvalue in the leading sub-matrix.
min(1,N) < N1 < N.
.TP 8
NB      (global input) INTEGER
The blocking factor used to distribute the columns of the
matrix. NB >= 1.
.TP 7
D      (input/output) REAL array, dimension (N)
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
Q      (input/output) REAL array, dimension (LDQ, N)
On entry, Q contains the eigenvectors of two submatrices in
the two square blocks with corners at (1,1), (N1,N1)
and (N1+1, N1+1), (N,N).
On exit, Q contains the trailing (N-K) updated eigenvectors
(those which were deflated) in its last N-K columns.
.TP 7
LDQ    (input) INTEGER
The leading dimension of the array Q.  LDQ >= max(1,NQ).
.TP 7
RHO    (global input/output) REAL
On entry, the off-diagonal element associated with the rank-1
cut which originally split the two submatrices which are now
being recombined.
On exit, RHO has been modified to the value required by
PSLAED3.
.TP 7
Z      (global input) REAL array, dimension (N)
On entry, Z contains the updating vector (the last
row of the first sub-eigenvector matrix and the first row of
the second sub-eigenvector matrix).
On exit, the contents of Z have been destroyed by the updating
process.

DLAMDA (global output) REAL array, dimension (N)
A copy of the first K eigenvalues which will be used by
SLAED3 to form the secular equation.
.TP 7
W      (global output) REAL array, dimension (N)
The first k values of the final deflation-altered z-vector
which will be passed to SLAED3.
.TP 7
Q2     (output) REAL array, dimension (LDQ2, NQ)
A copy of the first K eigenvectors which will be used by
.TP 8
LDQ2    (input) INTEGER
The leading dimension of the array Q2.
.TP 7
QBUF   (workspace) REAL array, dimension 3*N
.TP 7
CTOT   (workspace) INTEGER array, dimension( NPCOL, 4)
.TP 7
PSM    (workspace) INTEGER array, dimension( NPCOL, 4)
.TP 8
NPCOL   (global input) INTEGER
The total number of columns over which the distributed
submatrix is distributed.
.TP 7
INDX   (workspace) INTEGER array, dimension (N)
The permutation used to sort the contents of DLAMDA into
ascending order.
.TP 7
INDXC  (output) INTEGER array, dimension (N)
The permutation used to arrange the columns of the deflated
Q matrix into three groups:  the first group contains non-zero
elements only at and above N1, the second contains
non-zero elements only below N1, and the third is dense.
.TP 7
INDXP  (workspace) INTEGER array, dimension (N)
The permutation used to place deflated values of D at the end
of the array.  INDXP(1:K) points to the nondeflated D-values
.br
and INDXP(K+1:N) points to the deflated eigenvalues.

INDCOL (workspace) INTEGER array, dimension (N)

COLTYP (workspace/output) INTEGER array, dimension (N)
During execution, a label which will indicate which of the
following types a column in the Q2 matrix is:
.br
1 : non-zero in the upper half only;
.br
2 : dense;
.br
3 : non-zero in the lower half only;
.br
4 : deflated.
.TP 7
NN     (global output) INTEGER, the order of matrix U, (PSLAED1).
NN1    (global output) INTEGER, the order of matrix Q1, (PSLAED1).
NN2    (global output) INTEGER, the order of matrix Q2, (PSLAED1).
IB1    (global output) INTEGER, pointeur on Q1, (PSLAED1).
IB2    (global output) INTEGER, pointeur on Q2, (PSLAED1).
