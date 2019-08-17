.TH PSLACP3 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLACP3 - i an auxiliary routine that copies from a global parallel array into a local replicated array or vise versa
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLACP3(
M, I, A, DESCA, B, LDB, II, JJ, REV )
.TP 20
.ti +4
INTEGER
I, II, JJ, LDB, M, REV
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
A( * ), B( LDB, * )
.SH PURPOSE
PSLACP3 is an auxiliary routine that copies from a global parallel array into a local replicated array or vise versa. Notice that   the entire submatrix that is copied gets placed on one node or
  more.  The receiving node can be specified precisely, or all nodes
  can receive, or just one row or column of nodes.
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
M is the order of the square submatrix that is copied.
M >= 0.
Unchanged on exit
.TP 8
I       (global input) INTEGER
A(I,I) is the global location that the copying starts from.
Unchanged on exit.
.TP 8
A       (global input/output) REAL array, dimension
(DESCA(LLD_),*)
On entry, the parallel matrix to be copied into or from.
On exit, if REV=1, the copied data.
Unchanged on exit if REV=0.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
B       (local input/output) REAL array of size (LDB,M)
If REV=0, this is the global portion of the array
A(I:I+M-1,I:I+M-1).
If REV=1, this is the unchanged on exit.
.TP 8
LDB     (local input) INTEGER
The leading dimension of B.
.TP 8
II      (global input) INTEGER
By using REV 0 & 1, data can be sent out and returned again.
If REV=0, then II is destination row index for the node(s)
receiving the replicated B.  
If II>=0,JJ>=0, then node (II,JJ) receives the data
If II=-1,JJ>=0, then all rows in column JJ receive the
data
If II>=0,JJ=-1, then all cols in row II receive the data
If II=-1,JJ=-1, then all nodes receive the data
If REV<>0, then II is the source row index for the node(s)
sending the replicated B.
.TP 8
JJ      (global input) INTEGER
Similar description as II above
.TP 8
REV     (global input) INTEGER
Use REV = 0 to send global A into locally replicated B 
(on node (II,JJ)).
Use REV <> 0 to send locally replicated B from node (II,JJ)
to its owner (which changes depending on its location in
A) into the global A.

Implemented by:  G. Henry, May 1, 1997
