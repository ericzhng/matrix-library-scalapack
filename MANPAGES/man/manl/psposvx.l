.TH PSPOSVX l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSPOSVX - use the Cholesky factorization A = U**T*U or A = L*L**T to compute the solution to a real system of linear equations  A(IA:IA+N-1,JA:JA+N-1) * X = B(IB:IB+N-1,JB:JB+NRHS-1),
.SH SYNOPSIS
.TP 20
SUBROUTINE PSPOSVX(
FACT, UPLO, N, NRHS, A, IA, JA, DESCA, AF,
IAF, JAF, DESCAF, EQUED, SR, SC, B, IB, JB,
DESCB, X, IX, JX, DESCX, RCOND, FERR, BERR,
WORK, LWORK, IWORK, LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
EQUED, FACT, UPLO
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
DESCX( * ), IWORK( * )
.TP 20
.ti +4
REAL
A( * ), AF( * ),
B( * ), BERR( * ), FERR( * ),
SC( * ), SR( * ), WORK( * ), X( * )
.SH PURPOSE
PSPOSVX uses the Cholesky factorization A = U**T*U or A = L*L**T to compute the solution to a real system of linear equations A(IA:IA+N-1,JA:JA+N-1) * X = B(IB:IB+N-1,JB:JB+NRHS-1), 
where A(IA:IA+N-1,JA:JA+N-1) is an N-by-N matrix and X and
B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS matrices.
.br

Error bounds on the solution and a condition estimate are also
provided.  In the following comments Y denotes Y(IY:IY+M-1,JY:JY+K-1)
a M-by-K matrix where Y can be A, AF, B and X.
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
The following steps are performed:
.br

1. If FACT = 'E', real scaling factors are computed to equilibrate
   the system:
.br
      diag(SR) * A * diag(SC) * inv(diag(SC)) * X = diag(SR) * B
   Whether or not the system will be equilibrated depends on the
   scaling of the matrix A, but if equilibration is used, A is
   overwritten by diag(SR)*A*diag(SC) and B by diag(SR)*B.

2. If FACT = 'N' or 'E', the Cholesky decomposition is used to
   factor the matrix A (after equilibration if FACT = 'E') as
      A = U**T* U,  if UPLO = 'U', or
.br
      A = L * L**T,  if UPLO = 'L',
.br
   where U is an upper triangular matrix and L is a lower triangular
   matrix.
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

6. If equilibration was used, the matrix X is premultiplied by
   diag(SR) so that it solves the original system before
.br
   equilibration.
.br

.SH ARGUMENTS
.TP 8
FACT    (global input) CHARACTER
Specifies whether or not the factored form of the matrix A is
supplied on entry, and if not, whether the matrix A should be
equilibrated before it is factored.
= 'F':  On entry, AF contains the factored form of A.
If EQUED = 'Y', the matrix A has been equilibrated
with scaling factors given by S.  A and AF will not
be modified.
= 'N':  The matrix A will be copied to AF and factored.
.br
= 'E':  The matrix A will be equilibrated if necessary, then
copied to AF and factored.
.TP 8
UPLO    (global input) CHARACTER
= 'U':  Upper triangle of A is stored;
.br
= 'L':  Lower triangle of A is stored.
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix A(IA:IA+N-1,JA:JA+N-1).
N >= 0.
.TP 8
NRHS    (global input) INTEGER
The number of right hand sides, i.e., the number of columns
of the distributed submatrices B and X.  NRHS >= 0.
.TP 8
A       (local input/local output) REAL pointer into
the local memory to an array of local dimension
( LLD_A, LOCc(JA+N-1) ).
On entry, the symmetric matrix A, except if FACT = 'F' and
EQUED = 'Y', then A must contain the equilibrated matrix
diag(SR)*A*diag(SC).  If UPLO = 'U', the leading
N-by-N upper triangular part of A contains the upper
triangular part of the matrix A, and the strictly lower
triangular part of A is not referenced.  If UPLO = 'L', the
leading N-by-N lower triangular part of A contains the lower
triangular part of the matrix A, and the strictly upper
triangular part of A is not referenced.  A is not modified if
FACT = 'F' or 'N', or if FACT = 'E' and EQUED = 'N' on exit.

On exit, if FACT = 'E' and EQUED = 'Y', A is overwritten by
diag(SR)*A*diag(SC).
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
( LLD_AF, LOCc(JA+N-1)).
If FACT = 'F', then AF is an input argument and on entry
contains the triangular factor U or L from the Cholesky
factorization A = U**T*U or A = L*L**T, in the same storage
format as A.  If EQUED .ne. 'N', then AF is the factored form
of the equilibrated matrix diag(SR)*A*diag(SC).

If FACT = 'N', then AF is an output argument and on exit
returns the triangular factor U or L from the Cholesky
factorization A = U**T*U or A = L*L**T of the original
matrix A.

If FACT = 'E', then AF is an output argument and on exit
returns the triangular factor U or L from the Cholesky
factorization A = U**T*U or A = L*L**T of the equilibrated
matrix A (see the description of A for the form of the
equilibrated matrix).
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
EQUED   (global input/global output) CHARACTER
Specifies the form of equilibration that was done.
= 'N':  No equilibration (always true if FACT = 'N').
.br
= 'Y':  Equilibration was done, i.e., A has been replaced by
diag(SR) * A * diag(SC).
EQUED is an input variable if FACT = 'F'; otherwise, it is an
output variable.
.TP 8
SR      (local input/local output) REAL array,
dimension (LLD_A)
The scale factors for A distributed across process rows;
not accessed if EQUED = 'N'.  SR is an input variable if
FACT = 'F'; otherwise, SR is an output variable.
If FACT = 'F' and EQUED = 'Y', each element of SR must be
positive.
.TP 8
SC      (local input/local output) REAL array,
dimension (LOC(N_A))
The scale factors for A distributed across
process columns; not accessed if EQUED = 'N'. SC is an input
variable if FACT = 'F'; otherwise, SC is an output variable.
If FACT = 'F' and EQUED = 'Y', each element of SC must be
positive.
.TP 8
B       (local input/local output) REAL pointer into
the local memory to an array of local dimension
( LLD_B, LOCc(JB+NRHS-1) ).
On entry, the N-by-NRHS right-hand side matrix B.
On exit, if EQUED = 'N', B is not modified; if TRANS = 'N'
and EQUED = 'R' or 'B', B is overwritten by diag(R)*B; if
TRANS = 'T' or 'C' and EQUED = 'C' or 'B', B is overwritten
by diag(C)*B.
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
X       (local input/local output) REAL pointer into
the local memory to an array of local dimension
( LLD_X, LOCc(JX+NRHS-1) ).
If INFO = 0, the N-by-NRHS solution matrix X to the original
system of equations.  Note that A and B are modified on exit
if EQUED .ne. 'N', and the solution to the equilibrated
system is inv(diag(SC))*X if TRANS = 'N' and EQUED = 'C' or
'B', or inv(diag(SR))*X if TRANS = 'T' or 'C' and EQUED = 'R'
or 'B'.
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
A after equilibration (if done).  If RCOND is less than the
machine precision (in particular, if RCOND = 0), the matrix
is singular to working precision.  This condition is
indicated by a return code of INFO > 0, and the solution and
error bounds are not computed.
.TP 8
FERR    (local output) REAL array, dimension (LOC(N_B))
The estimated forward error bounds for each solution vector
X(j) (the j-th column of the solution matrix X).
If XTRUE is the true solution, FERR(j) bounds the magnitude
of the largest entry in (X(j) - XTRUE) divided by
the magnitude of the largest entry in X(j).  The quality of
the error bound depends on the quality of the estimate of
norm(inv(A)) computed in the code; if the estimate of
norm(inv(A)) is accurate, the error bound is guaranteed.
.TP 8
BERR    (local output) REAL array, dimension (LOC(N_B))
The componentwise relative backward error of each solution
vector X(j) (i.e., the smallest relative change in
any entry of A or B that makes X(j) an exact solution).
.TP 8
WORK    (local workspace/local output) REAL array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK = MAX( PSPOCON( LWORK ), PSPORFS( LWORK ) )
+ LOCr( N_A ).
LWORK = 3*DESCA( LLD_ )

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
LIWORK = DESCA( LLD_ )
LIWORK = LOCr(N_A).

If LIWORK = -1, then LIWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
INFO    (global output) INTEGER
= 0: successful exit
.br
< 0: if INFO = -i, the i-th argument had an illegal value
.br
> 0: if INFO = i, and i is
.br
<= N: if INFO = i, the leading minor of order i of A
is not positive definite, so the factorization
could not be completed, and the solution and error
bounds could not be computed.
= N+1: RCOND is less than machine precision.  The
factorization has been completed, but the matrix
is singular to working precision, and the solution
and error bounds have not been computed.
