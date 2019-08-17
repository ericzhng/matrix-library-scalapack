.TH PZLACON l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZLACON - estimate the 1-norm of a square, complex distributed matrix A
.SH SYNOPSIS
.TP 20
SUBROUTINE PZLACON(
N, V, IV, JV, DESCV, X, IX, JX, DESCX, EST,
KASE )
.TP 20
.ti +4
INTEGER
IV, IX, JV, JX, KASE, N
.TP 20
.ti +4
DOUBLE
PRECISION EST
.TP 20
.ti +4
INTEGER
DESCV( * ), DESCX( * )
.TP 20
.ti +4
COMPLEX*16
V( * ), X( * )
.SH PURPOSE
PZLACON estimates the 1-norm of a square, complex distributed matrix A. Reverse communication is used for evaluating matrix-vector products. X and V are aligned with the distributed matrix A, this
information is implicitly contained within IV, IX, DESCV, and DESCX.

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
The length of the distributed vectors V and X.  N >= 0.
.TP 8
V       (local workspace) COMPLEX*16 pointer into the local
memory to an array of dimension LOCr(N+MOD(IV-1,MB_V)). On
the final return, V = A*W, where EST = norm(V)/norm(W)
(W is not returned).
.TP 8
IV      (global input) INTEGER
The row index in the global array V indicating the first
row of sub( V ).
.TP 8
JV      (global input) INTEGER
The column index in the global array V indicating the
first column of sub( V ).
.TP 8
DESCV   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix V.
.TP 8
X       (local input/local output) COMPLEX*16 pointer into the
local memory to an array of dimension
LOCr(N+MOD(IX-1,MB_X)). On an intermediate return, X
should be overwritten by
A * X,   if KASE=1,
A' * X,  if KASE=2,
where A' is the conjugate transpose of A, and PZLACON must
be re-called with all the other parameters unchanged.
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
EST     (global output) DOUBLE PRECISION
An estimate (a lower bound) for norm(A).
.TP 8
KASE    (local input/local output) INTEGER
On the initial call to PZLACON, KASE should be 0.
On an intermediate return, KASE will be 1 or 2, indicating
whether X should be overwritten by A * X  or A' * X.
On the final return from PZLACON, KASE will again be 0.
.SH FURTHER DETAILS
The serial version ZLACON has been contributed by Nick Higham,
University of Manchester. It was originally named SONEST, dated
March 16, 1988.
.br

Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
a real or complex matrix, with applications to condition estimation",
ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.

