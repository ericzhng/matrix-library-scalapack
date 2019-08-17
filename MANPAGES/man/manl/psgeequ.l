.TH PSGEEQU l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSGEEQU - compute row and column scalings intended to equilibrate an M-by-N distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA:JA+N-1) and reduce its condition number
.SH SYNOPSIS
.TP 20
SUBROUTINE PSGEEQU(
M, N, A, IA, JA, DESCA, R, C, ROWCND, COLCND,
AMAX, INFO )
.TP 20
.ti +4
INTEGER
IA, INFO, JA, M, N
.TP 20
.ti +4
REAL
AMAX, COLCND, ROWCND
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
A( * ), C( * ), R( * )
.SH PURPOSE
PSGEEQU computes row and column scalings intended to equilibrate an M-by-N distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA:JA+N-1) and reduce its condition number. R returns the row scale factors and C the column scale factors, chosen to try to make the largest entry in
each row and column of the distributed matrix B with elements
B(i,j) = R(i) * A(i,j) * C(j) have absolute value 1.
.br

R(i) and C(j) are restricted to be between SMLNUM = smallest safe
number and BIGNUM = largest safe number.  Use of these scaling
factors is not guaranteed to reduce the condition number of
sub( A ) but works well in practice.
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
M       (global input) INTEGER
The number of rows to be operated on i.e the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input) REAL pointer into the local memory
to an array of dimension ( LLD_A, LOCc(JA+N-1) ), the
local pieces of the M-by-N distributed matrix whose
equilibration factors are to be computed.
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
R       (local output) REAL array, dimension LOCr(M_A)
If INFO = 0 or INFO > IA+M-1, R(IA:IA+M-1) contains the row
scale factors for sub( A ). R is aligned with the distributed
matrix A, and replicated across every process column. R is
tied to the distributed matrix A.
.TP 8
C       (local output) REAL array, dimension LOCc(N_A)
If INFO = 0,  C(JA:JA+N-1) contains the column scale factors
for sub( A ). C is aligned with the distributed matrix A, and
replicated down every process row. C is tied to the distri-
buted matrix A.
.TP 8
ROWCND  (global output) REAL
If INFO = 0 or INFO > IA+M-1, ROWCND contains the ratio of
the smallest R(i) to the largest R(i) (IA <= i <= IA+M-1).
If ROWCND >= 0.1 and AMAX is neither too large nor too small,
it is not worth scaling by R(IA:IA+M-1).
.TP 8
COLCND  (global output) REAL
If INFO = 0, COLCND contains the ratio of the smallest C(j)
to the largest C(j) (JA <= j <= JA+N-1). If COLCND >= 0.1, it
is not worth scaling by C(JA:JA+N-1).
.TP 8
AMAX    (global output) REAL
Absolute value of largest distributed matrix element.  If
AMAX is very close to overflow or very close to underflow,
the matrix should be scaled.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  If INFO = i,  and i is
.br
<= M:  the i-th row of the distributed matrix sub( A )
is exactly zero,
>  M:  the (i-M)-th column of the distributed
matrix sub( A ) is exactly zero.
