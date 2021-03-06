.TH PZLAUU2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZLAUU2 - compute the product U * U' or L' * L, where the triangular factor U or L is stored in the upper or lower triangular part of the matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1)
.SH SYNOPSIS
.TP 20
SUBROUTINE PZLAUU2(
UPLO, N, A, IA, JA, DESCA )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
IA, JA, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
COMPLEX*16
A( * )
.SH PURPOSE
PZLAUU2 computes the product U * U' or L' * L, where the triangular factor U or L is stored in the upper or lower triangular part of the matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1). 
If UPLO = 'U' or 'u' then the upper triangle of the result is stored,
overwriting the factor U in sub( A ).
.br
If UPLO = 'L' or 'l' then the lower triangle of the result is stored,
overwriting the factor L in sub( A ).
.br

This is the unblocked form of the algorithm, calling Level 2 BLAS.
No communication is performed by this routine, the matrix to operate
on should be strictly local to one process.
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
.TP 8
UPLO    (global input) CHARACTER*1
Specifies whether the triangular factor stored in the matrix
sub( A ) is upper or lower triangular:
.br
= 'U':  Upper triangular,
.br
= 'L':  Lower triangular.
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the order of the triangular factor U or L.  N >= 0.
.TP 8
A       (local input/local output) COMPLEX*16 pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
On entry, the local pieces of the triangular factor L or U.
On exit, if UPLO = 'U', the upper triangle of the distributed
matrix sub( A ) is overwritten with the upper triangle of the
product U * U'; if UPLO = 'L', the lower triangle of sub( A )
is overwritten with the lower triangle of the product L' * L.
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
