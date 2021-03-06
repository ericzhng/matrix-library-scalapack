.TH PZDRSCL l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZDRSCL - multiplie an N-element complex distributed vector sub( X ) by the real scalar 1/a
.SH SYNOPSIS
.TP 20
SUBROUTINE PZDRSCL(
N, SA, SX, IX, JX, DESCX, INCX )
.TP 20
.ti +4
INTEGER
IX, INCX, JX, N
.TP 20
.ti +4
DOUBLE
PRECISION SA
.TP 20
.ti +4
INTEGER
DESCX( * )
.TP 20
.ti +4
COMPLEX*16
SX( * )
.SH PURPOSE
PZDRSCL multiplies an N-element complex distributed vector sub( X ) by the real scalar 1/a. This is done without overflow or underflow as long as the final sub( X )/a does not overflow or
underflow.
.br

where sub( X ) denotes X(IX:IX+N-1,JX:JX), if INCX = 1,
.br
                       X(IX:IX,JX:JX+N-1), if INCX = M_X.

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
Such a global array has an associated description vector descA.
In the following comments, the character _ should be read as
"of the global array".
.br

NOTATION        STORED IN      EXPLANATION
.br
--------------- -------------- --------------------------------------
DT_A   (global) descA[ DT_ ]   The descriptor type.  In this case,
                               DT_A = 1.
.br
CTXT_A (global) descA[ CTXT_ ] The BLACS context handle, indicating
                               the BLACS process grid A is distribu-
                               ted over. The context itself is glo-
                               bal, but the handle (the integer
                               value) may vary.
.br
M_A    (global) descA[ M_ ]    The number of rows in the global
                               array A.
.br
N_A    (global) descA[ N_ ]    The number of columns in the global
                               array A.
.br
MB_A   (global) descA[ MB_ ]   The blocking factor used to distribu-
                               te the rows of the array.
.br
NB_A   (global) descA[ NB_ ]   The blocking factor used to distribu-
                               te the columns of the array.
RSRC_A (global) descA[ RSRC_ ] The process row over which the first
                               row of the array A is distributed.
CSRC_A (global) descA[ CSRC_ ] The process column over which the
                               first column of the array A is
                               distributed.
.br
LLD_A  (local)  descA[ LLD_ ]  The leading dimension of the local
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

Because vectors may be seen as particular matrices, a distributed
vector is considered to be a distributed matrix.
.br

.SH ARGUMENTS
.TP 8
N       (global input) pointer to INTEGER
The number of components of the distributed vector sub( X ).
N >= 0.
.TP 8
SA      (global input) DOUBLE PRECISION
The scalar a which is used to divide each component of
sub( X ).  SA must be >= 0, or the subroutine will divide by
zero.
.TP 8
SX      (local input/local output) COMPLEX*16 array
containing the local pieces of a distributed matrix of
dimension of at least
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
