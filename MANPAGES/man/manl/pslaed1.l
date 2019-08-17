.TH PSLAED1 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAED1 - compute the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix,
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLAED1(
N, N1, D, ID, Q, IQ, JQ, DESCQ, RHO, WORK,
IWORK, INFO )
.TP 20
.ti +4
INTEGER
ID, INFO, IQ, JQ, N, N1
.TP 20
.ti +4
REAL
RHO
.TP 20
.ti +4
INTEGER
DESCQ( * ), IWORK( * )
.TP 20
.ti +4
REAL
D( * ), Q( * ), WORK( * )
.SH PURPOSE
PSLAED1 computes the updated eigensystem of a diagonal matrix after modification by a rank-one symmetric matrix, in parallel.
.br

  T = Q(in) ( D(in) + RHO * Z*Z' ) Q'(in) = Q(out) * D(out) * Q'(out)

   where Z = Q'u, u is a vector of length N with ones in the
   N1 and N1 + 1 th elements and zeros elsewhere.
.br

   The eigenvectors of the original matrix are stored in Q, and the
   eigenvalues are in D.  The algorithm consists of three stages:

      The first stage consists of deflating the size of the problem
      when there are multiple eigenvalues or if there is a zero in
      the Z vector.  For each such occurence the dimension of the
      secular equation problem is reduced by one.  This stage is
      performed by the routine PSLAED2.
.br

      The second stage consists of calculating the updated
      eigenvalues. This is done by finding the roots of the secular
      equation via the routine SLAED4 (as called by PSLAED3).
      This routine also calculates the eigenvectors of the current
      problem.
.br

      The final stage consists of computing the updated eigenvectors
      directly using the updated eigenvalues.  The eigenvectors for
      the current problem are multiplied with the eigenvectors from
      the overall problem.
.br

.SH ARGUMENTS
.TP 8
N       (global input) INTEGER
The order of the tridiagonal matrix T.  N >= 0.
.TP 8
N1      (input) INTEGER
The location of the last eigenvalue in the leading
sub-matrix.
min(1,N) <= N1 <= N.
.TP 8
D       (global input/output) REAL array, dimension (N)
On entry,the eigenvalues of the rank-1-perturbed matrix.
On exit, the eigenvalues of the repaired matrix.
.TP 8
ID      (global input) INTEGER
Q's global row/col index, which points to the beginning
of the submatrix which is to be operated on.
.TP 8
Q       (local output) REAL array,
global dimension (N, N),
local dimension ( LLD_Q, LOCc(JQ+N-1))
Q  contains the orthonormal eigenvectors of the symmetric
tridiagonal matrix.
.TP 8
IQ      (global input) INTEGER
Q's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JQ      (global input) INTEGER
Q's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCQ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Z.
.TP 7
RHO    (input) REAL
The subdiagonal entry used to create the rank-1 modification.
.TP 8
WORK    (local workspace/output) REAL array,
dimension 6*N + 2*NP*NQ
.TP 8
IWORK   (local workspace/output) INTEGER array,
dimension 7*N + 8*NPCOL + 2
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  The algorithm failed to compute the ith eigenvalue.
