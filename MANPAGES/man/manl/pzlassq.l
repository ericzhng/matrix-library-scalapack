.TH PZLASSQ l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZLASSQ - return the values scl and smsq such that  ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
.SH SYNOPSIS
.TP 20
SUBROUTINE PZLASSQ(
N, X, IX, JX, DESCX, INCX, SCALE, SUMSQ )
.TP 20
.ti +4
INTEGER
IX, INCX, JX, N
.TP 20
.ti +4
DOUBLE
PRECISION SCALE, SUMSQ
.TP 20
.ti +4
INTEGER
DESCX( * )
.TP 20
.ti +4
COMPLEX*16
X( * )
.SH PURPOSE
PZLASSQ returns the values scl and smsq such that ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq, 
where x( i ) = sub( X ) = abs( X( IX+(JX-1)*DESCX(M_)+(i-1)*INCX ) ).
The value of sumsq is assumed to be at least unity and the value of
ssq will then satisfy
.br

   1.0 .le. ssq .le. ( sumsq + 2*n ).
.br

scale is assumed to be non-negative and scl returns the value

   scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
          i
.br

scale and sumsq must be supplied in SCALE and SUMSQ respectively.
SCALE and SUMSQ are overwritten by scl and ssq respectively.

The routine makes only one pass through the vector sub( X ).

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

The result are only available in the scope of sub( X ), i.e if
sub( X ) is distributed along a process row, the correct results are
only available in this process row of the grid. Similarly if sub( X )
is distributed along a process column, the correct results are only
available in this process column of the grid.
.br

.SH ARGUMENTS
.TP 8
N       (global input) INTEGER
The length of the distributed vector sub( X ).
.TP 8
X       (input) COMPLEX*16
The vector for which a scaled sum of squares is computed.
x( i )  = X(IX+(JX-1)*M_X +(i-1)*INCX ), 1 <= i <= n.
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
SCALE   (local input/local output) DOUBLE PRECISION
On entry, the value  scale  in the equation above.
On exit, SCALE is overwritten with  scl , the scaling factor
for the sum of squares.
.TP 8
SUMSQ   (local input/local output) DOUBLE PRECISION
On entry, the value  sumsq  in the equation above.
On exit, SUMSQ is overwritten with  smsq , the basic sum of
squares from which  scl  has been factored out.
