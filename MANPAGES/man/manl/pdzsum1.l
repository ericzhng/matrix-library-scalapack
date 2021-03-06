.TH PDZSUM1 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDZSUM1 - return the sum of absolute values of a complex distributed vector sub( X ) in ASUM,
.SH SYNOPSIS
.TP 20
SUBROUTINE PDZSUM1(
N, ASUM, X, IX, JX, DESCX, INCX )
.TP 20
.ti +4
INTEGER
IX, INCX, JX, N
.TP 20
.ti +4
DOUBLE
PRECISION ASUM
.TP 20
.ti +4
INTEGER
DESCX( * )
.TP 20
.ti +4
COMPLEX*16
X( * )
.SH PURPOSE
PDZSUM1 returns the sum of absolute values of a complex distributed vector sub( X ) in ASUM, 
where sub( X ) denotes X(IX:IX+N-1,JX:JX), if INCX = 1,
.br
                       X(IX:IX,JX:JX+N-1), if INCX = M_X.

Based on PDZASUM from the Level 1 PBLAS. The change is
.br
to use the 'genuine' absolute value.
.br

The serial version of this routine was originally contributed by
Nick Higham for use with ZLACON.
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

Because vectors may be viewed as a subclass of matrices, a
distributed vector is considered to be a distributed matrix.

When the result of a vector-oriented PBLAS call is a scalar, it will
be made available only within the scope which owns the vector(s)
being operated on.  Let X be a generic term for the input vector(s).
Then, the processes which receive the answer will be (note that if
an operation involves more than one vector, the processes which re-
ceive the result will be the union of the following calculation for
each vector):
.br

If N = 1, M_X = 1 and INCX = 1, then one can't determine if a process
row or process column owns the vector operand, therefore only the
process of coordinate {RSRC_X, CSRC_X} receives the result;

If INCX = M_X, then sub( X ) is a vector distributed over a process
row. Each process part of this row receives the result;
.br

If INCX = 1, then sub( X ) is a vector distributed over a process
column. Each process part of this column receives the result;

.SH PARAMETERS
.TP 8
N       (global input) pointer to INTEGER
The number of components of the distributed vector sub( X ).
N >= 0.
.TP 8
ASUM    (local output) pointer to DOUBLE PRECISION
The sum of absolute values of the distributed vector sub( X )
only in its scope.
.TP 8
X       (local input) COMPLEX*16 array containing the local
pieces of a distributed matrix of dimension of at least
( (JX-1)*M_X + IX + ( N - 1 )*abs( INCX ) )
This array contains the entries of the distributed vector
sub( X ).
.TP 8
IX      (global input) pointer to INTEGER
The global row index of the submatrix of the distributed
matrix X to operate on.
.TP 8
JX      (global input) pointer to INTEGER
The global column index of the submatrix of the distributed
matrix X to operate on.
.TP 8
DESCX   (global and local input) INTEGER array of dimension 8.
The array descriptor of the distributed matrix X.
.TP 8
INCX    (global input) pointer to INTEGER
The global increment for the elements of X. Only two values
of INCX are supported in this version, namely 1 and M_X.
