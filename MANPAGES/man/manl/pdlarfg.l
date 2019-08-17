.TH PDLARFG l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLARFG - generate a real elementary reflector H of order n, such that  H * sub( X ) = H * ( x(iax,jax) ) = ( alpha ), H' * H = I
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLARFG(
N, ALPHA, IAX, JAX, X, IX, JX, DESCX, INCX,
TAU )
.TP 20
.ti +4
INTEGER
IAX, INCX, IX, JAX, JX, N
.TP 20
.ti +4
DOUBLE
PRECISION ALPHA
.TP 20
.ti +4
INTEGER
DESCX( * )
.TP 20
.ti +4
DOUBLE
PRECISION TAU( * ), X( * )
.SH PURPOSE
PDLARFG generates a real elementary reflector H of order n, such that H * sub( X ) = H * ( x(iax,jax) ) = ( alpha ), H' * H = I.                       (      x     )   (   0   )
.br

where alpha is a scalar, and sub( X ) is an (N-1)-element real
distributed vector X(IX:IX+N-2,JX) if INCX = 1 and X(IX,JX:JX+N-2) if
INCX = DESCX(M_).  H is represented in the form
.br

      H = I - tau * ( 1 ) * ( 1 v' ) ,
.br
                    ( v )
.br

where tau is a real scalar and v is a real (N-1)-element
.br
vector.
.br

If the elements of sub( X ) are all zero, then tau = 0 and H is
taken to be the unit matrix.
.br

Otherwise  1 <= tau <= 2.
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

.SH ARGUMENTS
.TP 8
N       (global input) INTEGER
The global order of the elementary reflector. N >= 0.
.TP 8
ALPHA   (local output) DOUBLE PRECISION
On exit, alpha is computed in the process scope having the
vector sub( X ).
.TP 8
IAX     (global input) INTEGER
The global row index in X of X(IAX,JAX).
.TP 8
JAX     (global input) INTEGER
The global column index in X of X(IAX,JAX).
.TP 8
X       (local input/local output) DOUBLE PRECISION, pointer into the
local memory to an array of dimension (LLD_X,*). This array
contains the local pieces of the distributed vector sub( X ).
Before entry, the incremented array sub( X ) must contain
the vector x. On exit, it is overwritten with the vector v.
.TP 8
IX      (global input) INTEGER
The row index in the global array X indicating the first
row of sub( X ).
.TP 8
JX      (global input) INTEGER
The column index in the global array X indicating the
first column of sub( X ).
.TP 8
DESCX   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix X.
.TP 8
INCX    (global input) INTEGER
The global increment for the elements of X. Only two values
of INCX are supported in this version, namely 1 and M_X.
INCX must not be zero.
.TP 8
TAU     (local output) DOUBLE PRECISION array, dimension  LOCc(JX)
if INCX = 1, and LOCr(IX) otherwise. This array contains the
Householder scalars related to the Householder vectors.
TAU is tied to the distributed matrix X.
