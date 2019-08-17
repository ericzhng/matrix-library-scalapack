.TH PSSYGST l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSSYGST - reduce a real symmetric-definite generalized eigenproblem to standard form
.SH SYNOPSIS
.TP 20
SUBROUTINE PSSYGST(
IBTYPE, UPLO, N, A, IA, JA, DESCA, B, IB, JB,
DESCB, SCALE, INFO )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
IA, IB, IBTYPE, INFO, JA, JB, N
.TP 20
.ti +4
REAL
SCALE
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCB( * )
.TP 20
.ti +4
REAL
A( * ), B( * )
.SH PURPOSE
PSSYGST reduces a real symmetric-definite generalized eigenproblem to standard form. 
In the following sub( A ) denotes A( IA:IA+N-1, JA:JA+N-1 ) and
sub( B ) denotes B( IB:IB+N-1, JB:JB+N-1 ).
.br

If IBTYPE = 1, the problem is sub( A )*x = lambda*sub( B )*x,
and sub( A ) is overwritten by inv(U**T)*sub( A )*inv(U) or
inv(L)*sub( A )*inv(L**T)
.br

If IBTYPE = 2 or 3, the problem is sub( A )*sub( B )*x = lambda*x or
sub( B )*sub( A )*x = lambda*x, and sub( A ) is overwritten by
U*sub( A )*U**T or L**T*sub( A )*L.
.br

sub( B ) must have been previously factorized as U**T*U or L*L**T by
PSPOTRF.
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

.SH ARGUMENTS
.TP 9
IBTYPE   (global input) INTEGER
= 1: compute inv(U**T)*sub( A )*inv(U) or
inv(L)*sub( A )*inv(L**T);
= 2 or 3: compute U*sub( A )*U**T or L**T*sub( A )*L.
.TP 8
UPLO    (global input) CHARACTER
.br
= 'U':  Upper triangle of sub( A ) is stored and sub( B ) is
factored as U**T*U;
= 'L':  Lower triangle of sub( A ) is stored and sub( B ) is
factored as L*L**T.
.TP 8
N       (global input) INTEGER
The order of the matrices sub( A ) and sub( B ).  N >= 0.
.TP 8
A       (local input/local output) REAL pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
On entry, this array contains the local pieces of the
N-by-N symmetric distributed matrix sub( A ). If UPLO = 'U',
the leading N-by-N upper triangular part of sub( A ) contains
the upper triangular part of the matrix, and its strictly
lower triangular part is not referenced.  If UPLO = 'L', the
leading N-by-N lower triangular part of sub( A ) contains
the lower triangular part of the matrix, and its strictly
upper triangular part is not referenced.

On exit, if INFO = 0, the transformed matrix, stored in the
same format as sub( A ).
.TP 8
IA      (global input) INTEGER
A's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JA      (global input) INTEGER
A's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
B       (local input) REAL pointer into the local memory
to an array of dimension (LLD_B, LOCc(JB+N-1)). On entry,
this array contains the local pieces of the triangular factor
from the Cholesky factorization of sub( B ), as returned by
PSPOTRF.
.TP 8
IB      (global input) INTEGER
B's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JB      (global input) INTEGER
B's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCB   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix B.
.TP 8
SCALE   (global output) REAL
Amount by which the eigenvalues should be scaled to
compensate for the scaling performed in this routine.
At present, SCALE is always returned as 1.0, it is
returned here to allow for future enhancement.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
