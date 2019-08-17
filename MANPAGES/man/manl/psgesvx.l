.TH PSGESVX l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSGESVX - use the LU factorization to compute the solution to a real system of linear equations  A(IA:IA+N-1,JA:JA+N-1) * X = B(IB:IB+N-1,JB:JB+NRHS-1),
.SH SYNOPSIS
.TP 20
SUBROUTINE PSGESVX(
FACT, TRANS, N, NRHS, A, IA, JA, DESCA, AF,
IAF, JAF, DESCAF, IPIV, EQUED, R, C, B, IB,
JB, DESCB, X, IX, JX, DESCX, RCOND, FERR,
BERR, WORK, LWORK, IWORK, LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
EQUED, FACT, TRANS
.TP 20
.ti +4
INTEGER
IA, IAF, IB, INFO, IX, JA, JAF, JB, JX, LIWORK,
LWORK, N, NRHS
.TP 20
.ti +4
REAL
RCOND
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCAF( * ), DESCB( * ),
DESCX( * ), IPIV( * ), IWORK( * )
.TP 20
.ti +4
REAL
A( * ), AF( * ), B( * ), BERR( * ), C( * ),
FERR( * ), R( * ), WORK( * ), X( * )
.SH PURPOSE
PSGESVX uses the LU factorization to compute the solution to a real system of linear equations A(IA:IA+N-1,JA:JA+N-1) * X = B(IB:IB+N-1,JB:JB+NRHS-1), 
where A(IA:IA+N-1,JA:JA+N-1) is an N-by-N matrix and X and
B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS matrices.
.br

Error bounds on the solution and a condition estimate are also
provided.
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

.SH DESCRIPTION
In the following description, A denotes A(IA:IA+N-1,JA:JA+N-1),
B denotes B(IB:IB+N-1,JB:JB+NRHS-1) and X denotes
.br
X(IX:IX+N-1,JX:JX+NRHS-1).
.br

The following steps are performed:
.br

1. If FACT = 'E', real scaling factors are computed to equilibrate
   the system:
.br
      TRANS = 'N':  diag(R)*A*diag(C)     *inv(diag(C))*X = diag(R)*B
      TRANS = 'T': (diag(R)*A*diag(C))**T *inv(diag(R))*X = diag(C)*B
      TRANS = 'C': (diag(R)*A*diag(C))**H *inv(diag(R))*X = diag(C)*B
   Whether or not the system will be equilibrated depends on the
   scaling of the matrix A, but if equilibration is used, A is
   overwritten by diag(R)*A*diag(C) and B by diag(R)*B (if TRANS='N')
   or diag(C)*B (if TRANS = 'T' or 'C').
.br

2. If FACT = 'N' or 'E', the LU decomposition is used to factor the
   matrix A (after equilibration if FACT = 'E') as
.br
      A = P * L * U,
.br
   where P is a permutation matrix, L is a unit lower triangular
   matrix, and U is upper triangular.
.br

3. The factored form of A is used to estimate the condition number
   of the matrix A.  If the reciprocal of the condition number is
   less than machine precision, steps 4-6 are skipped.
.br

4. The system of equations is solved for X using the factored form
   of A.
.br

5. Iterative refinement is applied to improve the computed solution
   matrix and calculate error bounds and backward error estimates
   for it.
.br

6. If FACT = 'E' and equilibration was used, the matrix X is
   premultiplied by diag(C) (if TRANS = 'N') or diag(R) (if
   TRANS = 'T' or 'C') so that it solves the original system
   before equilibration.
.br

.SH ARGUMENTS
.TP 8
FACT    (global input) CHARACTER
Specifies whether or not the factored form of the matrix
A(IA:IA+N-1,JA:JA+N-1) is supplied on entry, and if not,
.br
whether the matrix A(IA:IA+N-1,JA:JA+N-1) should be
equilibrated before it is factored.
= 'F':  On entry, AF(IAF:IAF+N-1,JAF:JAF+N-1) and IPIV con-
.br
tain the factored form of A(IA:IA+N-1,JA:JA+N-1).
If EQUED is not 'N', the matrix
A(IA:IA+N-1,JA:JA+N-1) has been equilibrated with
scaling factors given by R and C.
A(IA:IA+N-1,JA:JA+N-1), AF(IAF:IAF+N-1,JAF:JAF+N-1),
and IPIV are not modified.
= 'N':  The matrix A(IA:IA+N-1,JA:JA+N-1) will be copied to
.br
AF(IAF:IAF+N-1,JAF:JAF+N-1) and factored.
.br
= 'E':  The matrix A(IA:IA+N-1,JA:JA+N-1) will be equili-
brated if necessary, then copied to
AF(IAF:IAF+N-1,JAF:JAF+N-1) and factored.
.TP 8
TRANS   (global input) CHARACTER
.br
Specifies the form of the system of equations:
.br
= 'N':  A(IA:IA+N-1,JA:JA+N-1) * X(IX:IX+N-1,JX:JX+NRHS-1)
.br
= B(IB:IB+N-1,JB:JB+NRHS-1)     (No transpose)
.br
= 'T':  A(IA:IA+N-1,JA:JA+N-1)**T * X(IX:IX+N-1,JX:JX+NRHS-1)
.br
= B(IB:IB+N-1,JB:JB+NRHS-1)  (Transpose)
.br
= 'C':  A(IA:IA+N-1,JA:JA+N-1)**H * X(IX:IX+N-1,JX:JX+NRHS-1)
.br
= B(IB:IB+N-1,JB:JB+NRHS-1)  (Transpose)
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix A(IA:IA+N-1,JA:JA+N-1).
N >= 0.
.TP 8
NRHS    (global input) INTEGER
The number of right-hand sides, i.e., the number of columns
of the distributed submatrices B(IB:IB+N-1,JB:JB+NRHS-1) and
.br
X(IX:IX+N-1,JX:JX+NRHS-1).  NRHS >= 0.
.TP 8
A       (local input/local output) REAL pointer into
the local memory to an array of local dimension
(LLD_A,LOCc(JA+N-1)).  On entry, the N-by-N matrix
A(IA:IA+N-1,JA:JA+N-1).  If FACT = 'F' and EQUED is not 'N',
.br
then A(IA:IA+N-1,JA:JA+N-1) must have been equilibrated by
.br
the scaling factors in R and/or C.  A(IA:IA+N-1,JA:JA+N-1) is
not modified if FACT = 'F' or  'N', or if FACT = 'E' and
EQUED = 'N' on exit.

On exit, if EQUED .ne. 'N', A(IA:IA+N-1,JA:JA+N-1) is scaled
as follows:
.br
EQUED = 'R':  A(IA:IA+N-1,JA:JA+N-1) :=
.br
diag(R) * A(IA:IA+N-1,JA:JA+N-1)
.br
EQUED = 'C':  A(IA:IA+N-1,JA:JA+N-1) :=
.br
A(IA:IA+N-1,JA:JA+N-1) * diag(C)
.br
EQUED = 'B':  A(IA:IA+N-1,JA:JA+N-1) :=
.br
diag(R) * A(IA:IA+N-1,JA:JA+N-1) * diag(C).
.TP 8
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
AF      (local input or local output) REAL pointer
into the local memory to an array of local dimension
(LLD_AF,LOCc(JA+N-1)).  If FACT = 'F', then
AF(IAF:IAF+N-1,JAF:JAF+N-1) is an input argument and on
entry contains the factors L and U from the factorization
A(IA:IA+N-1,JA:JA+N-1) = P*L*U as computed by PSGETRF.
If EQUED .ne. 'N', then AF is the factored form of the
equilibrated matrix A(IA:IA+N-1,JA:JA+N-1).

If FACT = 'N', then AF(IAF:IAF+N-1,JAF:JAF+N-1) is an output
argument and on exit returns the factors L and U from the
factorization A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the original
.br
matrix A(IA:IA+N-1,JA:JA+N-1).

If FACT = 'E', then AF(IAF:IAF+N-1,JAF:JAF+N-1) is an output
argument and on exit returns the factors L and U from the
factorization A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the equili-
.br
brated matrix A(IA:IA+N-1,JA:JA+N-1) (see the description of
.br
A(IA:IA+N-1,JA:JA+N-1) for the form of the equilibrated
matrix).
.TP 8
IAF     (global input) INTEGER
The row index in the global array AF indicating the first
row of sub( AF ).
.TP 8
JAF     (global input) INTEGER
The column index in the global array AF indicating the
first column of sub( AF ).
.TP 8
DESCAF  (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix AF.
.TP 8
IPIV    (local input or local output) INTEGER array, dimension
LOCr(M_A)+MB_A. If FACT = 'F', then IPIV is an input argu-
ment and on entry contains the pivot indices from the fac-
torization A(IA:IA+N-1,JA:JA+N-1) = P*L*U as computed by
PSGETRF; IPIV(i) -> The global row local row i was
swapped with.  This array must be aligned with
A( IA:IA+N-1, * ).

If FACT = 'N', then IPIV is an output argument and on exit
contains the pivot indices from the factorization
A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the original matrix
.br
A(IA:IA+N-1,JA:JA+N-1).

If FACT = 'E', then IPIV is an output argument and on exit
contains the pivot indices from the factorization
A(IA:IA+N-1,JA:JA+N-1) = P*L*U of the equilibrated matrix
.br
A(IA:IA+N-1,JA:JA+N-1).
.TP 8
EQUED   (global input or global output) CHARACTER
Specifies the form of equilibration that was done.
= 'N':  No equilibration (always true if FACT = 'N').
.br
= 'R':  Row equilibration, i.e., A(IA:IA+N-1,JA:JA+N-1) has
been premultiplied by diag(R).
= 'C':  Column equilibration, i.e., A(IA:IA+N-1,JA:JA+N-1)
has been postmultiplied by diag(C).
= 'B':  Both row and column equilibration, i.e.,
.br
A(IA:IA+N-1,JA:JA+N-1) has been replaced by
.br
diag(R) * A(IA:IA+N-1,JA:JA+N-1) * diag(C).
EQUED is an input variable if FACT = 'F'; otherwise, it is an
output variable.
.TP 8
R       (local input or local output) REAL array,
dimension LOCr(M_A).
The row scale factors for A(IA:IA+N-1,JA:JA+N-1).
.br
If EQUED = 'R' or 'B', A(IA:IA+N-1,JA:JA+N-1) is multiplied
on the left by diag(R); if EQUED='N' or 'C', R is not acces-
sed.  R is an input variable if FACT = 'F'; otherwise, R is
an output variable.
If FACT = 'F' and EQUED = 'R' or 'B', each element of R must
be positive.
R is replicated in every process column, and is aligned
with the distributed matrix A.
.TP 8
C       (local input or local output) REAL array,
dimension LOCc(N_A).
The column scale factors for A(IA:IA+N-1,JA:JA+N-1).
.br
If EQUED = 'C' or 'B', A(IA:IA+N-1,JA:JA+N-1) is multiplied
on the right by diag(C); if EQUED = 'N' or 'R', C is not
accessed.  C is an input variable if FACT = 'F'; otherwise,
C is an output variable.  If FACT = 'F' and EQUED = 'C' or
'B', each element of C must be positive.
C is replicated in every process row, and is aligned with
the distributed matrix A.
.TP 8
B       (local input/local output) REAL pointer
into the local memory to an array of local dimension
(LLD_B,LOCc(JB+NRHS-1) ).  On entry, the N-by-NRHS right-hand
side matrix B(IB:IB+N-1,JB:JB+NRHS-1). On exit, if
.br
EQUED = 'N', B(IB:IB+N-1,JB:JB+NRHS-1) is not modified; if
TRANS = 'N' and EQUED = 'R' or 'B', B is overwritten by
diag(R)*B(IB:IB+N-1,JB:JB+NRHS-1); if TRANS = 'T' or 'C'
.br
and EQUED = 'C' or 'B', B(IB:IB+N-1,JB:JB+NRHS-1) is over-
.br
written by diag(C)*B(IB:IB+N-1,JB:JB+NRHS-1).
.TP 8
IB      (global input) INTEGER
The row index in the global array B indicating the first
row of sub( B ).
.TP 8
JB      (global input) INTEGER
The column index in the global array B indicating the
first column of sub( B ).
.TP 8
DESCB   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix B.
.TP 8
X       (local input/local output) REAL pointer
into the local memory to an array of local dimension
(LLD_X, LOCc(JX+NRHS-1)).  If INFO = 0, the N-by-NRHS
solution matrix X(IX:IX+N-1,JX:JX+NRHS-1) to the original
.br
system of equations.  Note that A(IA:IA+N-1,JA:JA+N-1) and
.br
B(IB:IB+N-1,JB:JB+NRHS-1) are modified on exit if
EQUED .ne. 'N', and the solution to the equilibrated system
is inv(diag(C))*X(IX:IX+N-1,JX:JX+NRHS-1) if TRANS = 'N'
and EQUED = 'C' or 'B', or
inv(diag(R))*X(IX:IX+N-1,JX:JX+NRHS-1) if TRANS = 'T' or 'C'
and EQUED = 'R' or 'B'.
.TP 8
IX      (global input) INTEGER
The row index in the global array X indicating the first
row of sub( X ).
.TP 8
JX      (global input) INTEGER
The column index in the global array X indicating the
first column of sub( X ).
.TP 8
DESCX   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix X.
.TP 8
RCOND   (global output) REAL
The estimate of the reciprocal condition number of the matrix
A(IA:IA+N-1,JA:JA+N-1) after equilibration (if done).  If
RCOND is less than the machine precision (in particular, if
RCOND = 0), the matrix is singular to working precision.
This condition is indicated by a return code of INFO > 0.
.TP 8
FERR    (local output) REAL array, dimension LOCc(N_B)
The estimated forward error bounds for each solution vector
X(j) (the j-th column of the solution matrix
X(IX:IX+N-1,JX:JX+NRHS-1). If XTRUE is the true solution,
FERR(j) bounds the magnitude of the largest entry in
(X(j) - XTRUE) divided by the magnitude of the largest entry
in X(j).  The estimate is as reliable as the estimate for
RCOND, and is almost always a slight overestimate of the
true error.  FERR is replicated in every process row, and is
aligned with the matrices B and X.
.TP 8
BERR    (local output) REAL array, dimension LOCc(N_B).
The componentwise relative backward error of each solution
vector X(j) (i.e., the smallest relative change in
any entry of A(IA:IA+N-1,JA:JA+N-1) or
.br
B(IB:IB+N-1,JB:JB+NRHS-1) that makes X(j) an exact solution).
BERR is replicated in every process row, and is aligned
with the matrices B and X.
.TP 8
WORK    (local workspace/local output) REAL array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK = MAX( PSGECON( LWORK ), PSGERFS( LWORK ) )
+ LOCr( N_A ).

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
IWORK   (local workspace/local output) INTEGER array,
dimension (LIWORK)
On exit, IWORK(1) returns the minimal and optimal LIWORK.
.TP 8
LIWORK  (local or global input) INTEGER
The dimension of the array IWORK.
LIWORK is local input and must be at least
LIWORK = LOCr(N_A).

If LIWORK = -1, then LIWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
.br
> 0:  if INFO = i, and i is
.br
<= N:  U(IA+I-1,IA+I-1) is exactly zero.  The
factorization has been completed, but the
factor U is exactly singular, so the solution
and error bounds could not be computed.
= N+1: RCOND is less than machine precision.  The
factorization has been completed, but the
matrix is singular to working precision, and
the solution and error bounds have not been
computed.
