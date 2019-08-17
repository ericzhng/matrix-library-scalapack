.TH PDPOTRI l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDPOTRI - compute the inverse of a real symmetric positive definite distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the Cholesky factorization sub( A ) = U**T*U or L*L**T computed by PDPOTRF
.SH SYNOPSIS
.TP 20
SUBROUTINE PDPOTRI(
UPLO, N, A, IA, JA, DESCA, INFO )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
IA, INFO, JA, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * )
.SH PURPOSE
PDPOTRI computes the inverse of a real symmetric positive definite distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the Cholesky factorization sub( A ) = U**T*U or L*L**T computed by PDPOTRF. 
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
UPLO    (global input) CHARACTER*1
= 'U':  Upper triangle of sub( A ) is stored;
.br
= 'L':  Lower triangle of sub( A ) is stored.
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
On entry, the local pieces of the triangular factor U or L
from the Cholesky factorization of the distributed matrix
sub( A ) = U**T*U or  L*L**T, as computed by PDPOTRF.
On exit, the local pieces of the upper or lower triangle of
the (symmetric) inverse of sub( A ), overwriting the input
factor U or L.
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
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  If INFO = i, the (i,i) element of the factor U or L is
zero, and the inverse could not be computed.
