.TH PZGETRI l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZGETRI - compute the inverse of a distributed matrix using the LU factorization computed by PZGETRF
.SH SYNOPSIS
.TP 20
SUBROUTINE PZGETRI(
N, A, IA, JA, DESCA, IPIV, WORK, LWORK,
IWORK, LIWORK, INFO )
.TP 20
.ti +4
INTEGER
IA, INFO, JA, LIWORK, LWORK, N
.TP 20
.ti +4
INTEGER
DESCA( * ), IPIV( * ), IWORK( * )
.TP 20
.ti +4
COMPLEX*16
A( * ), WORK( * )
.SH PURPOSE
PZGETRI computes the inverse of a distributed matrix using the LU factorization computed by PZGETRF. This method inverts U and then computes the inverse of sub( A ) = A(IA:IA+N-1,JA:JA+N-1) denoted
InvA by solving the system InvA*L = inv(U) for InvA.
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
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) COMPLEX*16 pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
On entry, the local pieces of the L and U obtained by the
factorization sub( A ) = P*L*U computed by PZGETRF. On
exit, if INFO = 0, sub( A ) contains the inverse of the
original distributed matrix sub( A ).
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
IPIV    (local input) INTEGER array, dimension LOCr(M_A)+MB_A
keeps track of the pivoting information. IPIV(i) is the
global row index the local row i was swapped with.  This
array is tied to the distributed matrix A.
.TP 8
WORK    (local workspace/local output) COMPLEX*16 array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK = LOCr(N+MOD(IA-1,MB_A))*NB_A. WORK is used to keep a
copy of at most an entire column block of sub( A ).

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
The dimension of the array IWORK used as workspace for
physically transposing the pivots.
LIWORK is local input and must be at least
if NPROW == NPCOL then
LIWORK = LOCc( N_A + MOD(JA-1, NB_A) ) + NB_A,
else
LIWORK =  LOCc( N_A + MOD(JA-1, NB_A) ) +
MAX( CEIL(CEIL(LOCr(M_A)/MB_A)/(LCM/NPROW)),
NB_A )
where LCM is the least common multiple of process
rows and columns (NPROW and NPCOL).
end if

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
> 0:  If INFO = K, U(IA+K-1,IA+K-1) is exactly zero; the
matrix is singular and its inverse could not be
computed.
