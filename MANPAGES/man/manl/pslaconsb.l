.TH PSLACONSB l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLACONSB - look for two consecutive small subdiagonal elements by seeing the effect of starting a double shift QR iteration given by H44, H33, & H43H34 and see if this would make a subdiagonal negligible
.SH SYNOPSIS
.TP 22
SUBROUTINE PSLACONSB(
A, DESCA, I, L, M, H44, H33, H43H34, BUF,
LWORK )
.TP 22
.ti +4
INTEGER
I, L, LWORK, M
.TP 22
.ti +4
REAL
H33, H43H34, H44
.TP 22
.ti +4
INTEGER
DESCA( * )
.TP 22
.ti +4
REAL
A( * ), BUF( * )
.SH PURPOSE
PSLACONSB looks for two consecutive small subdiagonal elements by seeing the effect of starting a double shift QR iteration given by H44, H33, & H43H34 and see if this would make a subdiagonal negligible. 
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
A       (global input) REAL array, dimension
(DESCA(LLD_),*)
On entry, the Hessenberg matrix whose tridiagonal part is
being scanned.
Unchanged on exit.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
I       (global input) INTEGER
The global location of the bottom of the unreduced
submatrix of A.
Unchanged on exit.
.TP 8
L       (global input) INTEGER
The global location of the top of the unreduced submatrix
of A.
Unchanged on exit.
.TP 8
M       (global output) INTEGER
On exit, this yields the starting location of the QR double
shift.  This will satisfy: L <= M  <= I-2.

H44
H33
H43H34  (global input) REAL
These three values are for the double shift QR iteration.
.TP 8
BUF     (local output) REAL array of size LWORK.
.TP 8
LWORK   (global input) INTEGER
On exit, LWORK is the size of the work buffer.
This must be at least 7*Ceil( Ceil( (I-L)/HBL ) /
LCM(NPROW,NPCOL) )
Here LCM is least common multiple, and NPROWxNPCOL is the
logical grid size.

Logic:
======

Two consecutive small subdiagonal elements will stall
convergence of a double shift if their product is small
relatively even if each is not very small.  Thus it is
necessary to scan the "tridiagonal portion of the matrix."  In
the LAPACK algorithm DLAHQR, a loop of M goes from I-2 down to
L and examines
H(m,m),H(m+1,m+1),H(m+1,m),H(m,m+1),H(m-1,m-1),H(m,m-1), and
H(m+2,m-1).  Since these elements may be on separate
processors, the first major loop (10) goes over the tridiagonal
and has each node store whatever values of the 7 it has that
the node owning H(m,m) does not.  This will occur on a border
and can happen in no more than 3 locations per block assuming
square blocks.  There are 5 buffers that each node stores these
values:  a buffer to send diagonally down and right, a buffer
to send up, a buffer to send left, a buffer to send diagonally
up and left and a buffer to send right.  Each of these buffers
is actually stored in one buffer BUF where BUF(ISTR1+1) starts
the first buffer, BUF(ISTR2+1) starts the second, etc..  After
the values are stored, if there are any values that a node
needs, they will be sent and received.  Then the next major
loop passes over the data and searches for two consecutive
small subdiagonals.

Notes:

This routine does a global maximum and must be called by all
processes.

Implemented by:  G. Henry, November 17, 1996
