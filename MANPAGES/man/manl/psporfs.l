.TH PSPORFS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSPORFS - improve the computed solution to a system of linear equations when the coefficient matrix is symmetric positive definite and provides error bounds and backward error estimates for the solutions
.SH SYNOPSIS
.TP 20
SUBROUTINE PSPORFS(
UPLO, N, NRHS, A, IA, JA, DESCA, AF, IAF, JAF,
DESCAF, B, IB, JB, DESCB, X, IX, JX, DESCX,
FERR, BERR, WORK, LWORK, IWORK, LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
IA, IAF, IB, INFO, IX, JA, JAF, JB, JX,
LIWORK, LWORK, N, NRHS
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCAF( * ), DESCB( * ),
DESCX( * ), IWORK( * )
.TP 20
.ti +4
REAL
A( * ), AF( * ), B( * ),
BERR( * ), FERR( * ), WORK( * ), X( * )
.SH PURPOSE
PSPORFS improves the computed solution to a system of linear equations when the coefficient matrix is symmetric positive definite and provides error bounds and backward error estimates for the solutions. 
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

In the following comments, sub( A ), sub( X ) and sub( B ) denote
respectively A(IA:IA+N-1,JA:JA+N-1), X(IX:IX+N-1,JX:JX+NRHS-1) and
B(IB:IB+N-1,JB:JB+NRHS-1).
.br

.SH ARGUMENTS
.TP 8
UPLO    (global input) CHARACTER*1
Specifies whether the upper or lower triangular part of the
symmetric matrix sub( A ) is stored.
= 'U':  Upper triangular
.br
= 'L':  Lower triangular
.TP 8
N       (global input) INTEGER
The order of the matrix sub( A ).  N >= 0.
.TP 8
NRHS    (global input) INTEGER
The number of right hand sides, i.e., the number of columns
of the matrices sub( B ) and sub( X ).  NRHS >= 0.
.TP 8
A       (local input) REAL pointer into the local
memory to an array of local dimension (LLD_A,LOCc(JA+N-1) ).
This array contains the local pieces of the N-by-N symmetric
distributed matrix sub( A ) to be factored.
If UPLO = 'U', the leading N-by-N upper triangular part of
sub( A ) contains the upper triangular part of the matrix,
and its strictly lower triangular part is not referenced.
If UPLO = 'L', the leading N-by-N lower triangular part of
sub( A ) contains the lower triangular part of the distribu-
ted matrix, and its strictly upper triangular part is not
referenced.
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
AF      (local input) REAL pointer into the local memory
to an array of local dimension (LLD_AF,LOCc(JA+N-1)).
On entry, this array contains the factors L or U from the
Cholesky factorization sub( A ) = L*L**T or U**T*U, as
computed by PSPOTRF.
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
B       (local input) REAL pointer into the local memory
to an array of local dimension (LLD_B, LOCc(JB+NRHS-1) ).
On entry, this array contains the the local pieces of the
right hand sides sub( B ).
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
X       (local input) REAL pointer into the local memory
to an array of local dimension (LLD_X, LOCc(JX+NRHS-1) ).
On entry, this array contains the the local pieces of the
solution vectors sub( X ). On exit, it contains the
improved solution vectors.
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
FERR    (local output) REAL array of local dimension
LOCc(JB+NRHS-1).
The estimated forward error bound for each solution vector
of sub( X ).  If XTRUE is the true solution corresponding
to sub( X ), FERR is an estimated upper bound for the
magnitude of the largest element in (sub( X ) - XTRUE)
divided by the magnitude of the largest element in sub( X ).
The estimate is as reliable as the estimate for RCOND, and
is almost always a slight overestimate of the true error.
This array is tied to the distributed matrix X.
.TP 8
BERR    (local output) REAL array of local dimension
LOCc(JB+NRHS-1). The componentwise relative backward
error of each solution vector (i.e., the smallest re-
lative change in any entry of sub( A ) or sub( B )
that makes sub( X ) an exact solution).
This array is tied to the distributed matrix X.
.TP 8
WORK    (local workspace/local output) REAL array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= 3*LOCr( N + MOD( IA-1, MB_A ) )

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
LIWORK >= LOCr( N + MOD( IB-1, MB_B ) ).

If LIWORK = -1, then LIWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
.SH PARAMETERS

ITMAX is the maximum number of steps of iterative refinement.

Notes
=====

This routine temporarily returns when N <= 1.

The distributed submatrices op( A ) and op( AF ) (respectively
sub( X ) and sub( B ) ) should be distributed the same way on the
same processes. These conditions ensure that sub( A ) and sub( AF )
(resp. sub( X ) and sub( B ) ) are "perfectly" aligned.

Moreover, this routine requires the distributed submatrices sub( A ),
sub( AF ), sub( X ), and sub( B ) to be aligned on a block boundary,
i.e., if f(x,y) = MOD( x-1, y ):
f( IA, DESCA( MB_ ) ) = f( JA, DESCA( NB_ ) ) = 0,
f( IAF, DESCAF( MB_ ) ) = f( JAF, DESCAF( NB_ ) ) = 0,
f( IB, DESCB( MB_ ) ) = f( JB, DESCB( NB_ ) ) = 0, and
f( IX, DESCX( MB_ ) ) = f( JX, DESCX( NB_ ) ) = 0.
