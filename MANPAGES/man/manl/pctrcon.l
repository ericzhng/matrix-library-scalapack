.TH PCTRCON l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCTRCON - estimate the reciprocal of the condition number of a triangular distributed matrix A(IA:IA+N-1,JA:JA+N-1), in either the 1-norm or the infinity-norm
.SH SYNOPSIS
.TP 20
SUBROUTINE PCTRCON(
NORM, UPLO, DIAG, N, A, IA, JA, DESCA, RCOND,
WORK, LWORK, RWORK, LRWORK, INFO )
.TP 20
.ti +4
CHARACTER
DIAG, NORM, UPLO
.TP 20
.ti +4
INTEGER
IA, JA, INFO, LRWORK, LWORK, N
.TP 20
.ti +4
REAL
RCOND
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
RWORK( * )
.TP 20
.ti +4
COMPLEX
A( * ), WORK( * )
.SH PURPOSE
PCTRCON estimates the reciprocal of the condition number of a triangular distributed matrix A(IA:IA+N-1,JA:JA+N-1), in either the 1-norm or the infinity-norm. 
The norm of A(IA:IA+N-1,JA:JA+N-1) is computed and an estimate is
obtained for norm(inv(A(IA:IA+N-1,JA:JA+N-1))), then the reciprocal
of the condition number is computed as
.br
           RCOND = 1 / ( norm( A(IA:IA+N-1,JA:JA+N-1)      ) *
                         norm( inv(A(IA:IA+N-1,JA:JA+N-1)) ) ).

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

.SH ARGUMENTS
.TP 8
NORM    (global input) CHARACTER
Specifies whether the 1-norm condition number or the
infinity-norm condition number is required:
.br
= '1' or 'O':  1-norm;
.br
= 'I':         Infinity-norm.
.TP 8
UPLO    (global input) CHARACTER
.br
= 'U':  A(IA:IA+N-1,JA:JA+N-1) is upper triangular;
.br
= 'L':  A(IA:IA+N-1,JA:JA+N-1) is lower triangular.
.TP 8
DIAG    (global input) CHARACTER
.br
= 'N':  A(IA:IA+N-1,JA:JA+N-1) is non-unit triangular;
.br
= 'U':  A(IA:IA+N-1,JA:JA+N-1) is unit triangular.
.TP 8
N       (global input) INTEGER
.br
The order of the distributed matrix A(IA:IA+N-1,JA:JA+N-1).
N >= 0.
.TP 8
A       (local input) COMPLEX pointer into the local memory
to an array of dimension ( LLD_A, LOCc(JA+N-1) ). This array
contains the local pieces of the triangular distributed
matrix A(IA:IA+N-1,JA:JA+N-1). If UPLO = 'U', the leading
N-by-N upper triangular part of this distributed matrix con-
tains the upper triangular matrix, and its strictly lower
triangular part is not referenced.  If UPLO = 'L', the
leading N-by-N lower triangular part of this ditributed
matrix contains the lower triangular matrix, and the strictly
upper triangular part is not referenced. If DIAG = 'U', the
diagonal elements of A(IA:IA+N-1,JA:JA+N-1) are also not
referenced and are assumed to be 1.
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
RCOND   (global output) REAL
The reciprocal of the condition number of the distributed
matrix A(IA:IA+N-1,JA:JA+N-1), computed as
.br
RCOND = 1 / ( norm( A(IA:IA+N-1,JA:JA+N-1)      ) *
.br
norm( inv(A(IA:IA+N-1,JA:JA+N-1)) ) ).
.TP 8
WORK    (local workspace/local output) COMPLEX array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= 2*LOCr(N+MOD(IA-1,MB_A)) +
MAX( 2, MAX(NB_A*CEIL(P-1,Q),LOCc(N+MOD(JA-1,NB_A)) +
NB_A*CEIL(Q-1,P)) ).

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
RWORK   (local workspace/local output) REAL array,
dimension (LRWORK)
On exit, RWORK(1) returns the minimal and optimal LRWORK.
.TP 8
LRWORK  (local or global input) INTEGER
The dimension of the array RWORK.
LRWORK is local input and must be at least
LRWORK >= LOCc(N+MOD(JA-1,NB_A)).

If LRWORK = -1, then LRWORK is global input and a workspace
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
