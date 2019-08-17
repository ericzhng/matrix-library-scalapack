.TH PDLAWIL l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLAWIL - get the transform given by H44,H33, & H43H34 into V starting at row M
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLAWIL(
II, JJ, M, A, DESCA, H44, H33, H43H34, V )
.TP 20
.ti +4
INTEGER
II, JJ, M
.TP 20
.ti +4
DOUBLE
PRECISION H33, H43H34, H44
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), V( * )
.SH PURPOSE
PDLAWIL gets the transform given by H44,H33, & H43H34 into V starting at row M. 
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
II      (global input) INTEGER
Row owner of H(M+2,M+2)
.TP 8
JJ      (global input) INTEGER
Column owner of H(M+2,M+2)
.TP 8
M       (global input) INTEGER
On entry, this is where the transform starts (row M.)
Unchanged on exit.
.TP 8
A       (global input) DOUBLE PRECISION array, dimension
(DESCA(LLD_),*)
On entry, the Hessenberg matrix.
Unchanged on exit.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
Unchanged on exit.

H44
H33
H43H34  (global input) DOUBLE PRECISION
These three values are for the double shift QR iteration.
Unchanged on exit.
.TP 8
V       (global output) DOUBLE PRECISION array of size 3.
Contains the transform on ouput.

Implemented by:  G. Henry, November 17, 1996
