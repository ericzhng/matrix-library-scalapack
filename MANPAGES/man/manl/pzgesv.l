.TH PZGESV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZGESV - compute the solution to a complex system of linear equations  sub( A ) * X = sub( B ),
.SH SYNOPSIS
.TP 19
SUBROUTINE PZGESV(
N, NRHS, A, IA, JA, DESCA, IPIV, B, IB, JB,
DESCB, INFO )
.TP 19
.ti +4
INTEGER
IA, IB, INFO, JA, JB, N, NRHS
.TP 19
.ti +4
INTEGER
DESCA( * ), DESCB( * ), IPIV( * )
.TP 19
.ti +4
COMPLEX*16
A( * ), B( * )
.SH PURPOSE
PZGESV computes the solution to a complex system of linear equations sub( A ) * X = sub( B ), 
where sub( A ) = A(IA:IA+N-1,JA:JA+N-1) is an N-by-N distributed
matrix and X and sub( B ) = B(IB:IB+N-1,JB:JB+NRHS-1) are N-by-NRHS
distributed matrices.
.br

The LU decomposition with partial pivoting and row interchanges is
used to factor sub( A ) as sub( A ) = P * L * U, where P is a permu-
tation matrix, L is unit lower triangular, and U is upper triangular.
L and U are stored in sub( A ). The factored form of sub( A ) is then
used to solve the system of equations sub( A ) * X = sub( B ).

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

This routine requires square block decomposition ( MB_A = NB_A ).

.SH ARGUMENTS
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
NRHS    (global input) INTEGER
The number of right hand sides, i.e., the number of columns
of the distributed submatrix sub( A ). NRHS >= 0.
.TP 8
A       (local input/local output) COMPLEX*16 pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
On entry, the local pieces of the N-by-N distributed matrix
sub( A ) to be factored. On exit, this array contains the
local pieces of the factors L and U from the factorization
sub( A ) = P*L*U; the unit diagonal elements of L are not
stored.
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
IPIV    (local output) INTEGER array, dimension ( LOCr(M_A)+MB_A )
This array contains the pivoting information.
IPIV(i) -> The global row local row i was swapped with.
This array is tied to the distributed matrix A.
.TP 8
B       (local input/local output) COMPLEX*16 pointer into the
local memory to an array of dimension
(LLD_B,LOCc(JB+NRHS-1)).  On entry, the right hand side
distributed matrix sub( B ). On exit, if INFO = 0, sub( B )
is overwritten by the solution distributed matrix X.
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
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  If INFO = K, U(IA+K-1,JA+K-1) is exactly zero.
The factorization has been completed, but the factor U
is exactly singular, so the solution could not be
computed.
