.TH PCPOEQU l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCPOEQU - compute row and column scalings intended to equilibrate a distributed Hermitian positive definite matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) and reduce its condition number (with respect to the two-norm)
.SH SYNOPSIS
.TP 20
SUBROUTINE PCPOEQU(
N, A, IA, JA, DESCA, SR, SC, SCOND, AMAX,
INFO )
.TP 20
.ti +4
INTEGER
IA, INFO, JA, N
.TP 20
.ti +4
REAL
AMAX, SCOND
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
SC( * ), SR( * )
.TP 20
.ti +4
COMPLEX
A( * )
.SH PURPOSE
PCPOEQU computes row and column scalings intended to equilibrate a distributed Hermitian positive definite matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) and reduce its condition number (with respect to the two-norm). SR and SC contain the scale factors, S(i) = 1/sqrt(A(i,i)), chosen so that the scaled distri-
buted matrix B with elements B(i,j) = S(i)*A(i,j)*S(j) has ones on
the  diagonal.  This choice of SR and SC puts the condition number
of B within a factor N of the smallest possible condition number
over all possible diagonal scalings.
.br

The scaling factor are stored along process rows in SR and along
process columns in SC. The duplication of information simplifies
greatly the application of the factors.
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
The number of rows and columns to be operated on i.e the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input) COMPLEX pointer into the local memory to an
array of local dimension ( LLD_A, LOCc(JA+N-1) ), the
N-by-N Hermitian positive definite distributed matrix
sub( A ) whose scaling factors are to be computed.  Only the
diagonal elements of sub( A ) are referenced.
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
SR      (local output) REAL array, dimension LOCr(M_A)
If INFO = 0, SR(IA:IA+N-1) contains the row scale factors
for sub( A ). SR is aligned with the distributed matrix A,
and replicated across every process column. SR is tied to the
distributed matrix A.
.TP 8
SC      (local output) REAL array, dimension LOCc(N_A)
If INFO = 0, SC(JA:JA+N-1) contains the column scale factors
.br
for A(IA:IA+M-1,JA:JA+N-1). SC is aligned with the distribu-
ted matrix A, and replicated down every process row. SC is
tied to the distributed matrix A.
.TP 8
SCOND   (global output) REAL
If INFO = 0, SCOND contains the ratio of the smallest SR(i)
(or SC(j)) to the largest SR(i) (or SC(j)), with
IA <= i <= IA+N-1 and JA <= j <= JA+N-1. If SCOND >= 0.1
and AMAX is neither too large nor too small, it is not worth
scaling by SR (or SC).
.TP 8
AMAX    (global output) REAL
Absolute value of largest matrix element.  If AMAX is very
close to overflow or very close to underflow, the matrix
should be scaled.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  If INFO = K, the K-th diagonal entry of sub( A ) is
nonpositive.
