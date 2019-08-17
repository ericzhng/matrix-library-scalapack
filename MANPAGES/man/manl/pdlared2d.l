.TH PDLARED2D l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLARED2D - redistribute a 1D array  It assumes that the input array, BYROW, is distributed across columns and that all process rows contain the same copy of BYROW
.SH SYNOPSIS
.TP 22
SUBROUTINE PDLARED2D(
N, IA, JA, DESC, BYROW, BYALL, WORK, LWORK )
.TP 22
.ti +4
INTEGER
IA, JA, LWORK, N
.TP 22
.ti +4
INTEGER
DESC( * )
.TP 22
.ti +4
DOUBLE
PRECISION BYALL( * ), BYROW( * ), WORK( LWORK )
.SH PURPOSE
PDLARED2D redistributes a 1D array It assumes that the input array, BYROW, is distributed across columns and that all process rows contain the same copy of BYROW. The output array, BYALL, will be identical on all processes and will contain the entire array.
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

NP = Number of local rows in BYROW()
.TP 8
N       (global input) INTEGER
The number of elements to be redistributed.  N >= 0.
.TP 8
IA      (global input) INTEGER
IA must be equal to 1
.TP 8
JA      (global input) INTEGER
JA must be equal to 1
.TP 8
DESC    (global/local input) INTEGER Array of dimension DLEN_
A 2D array descriptor, which describes BYROW
.TP 8
BYROW   (local input) distributed block cyclic DOUBLE PRECISION array
global dimension (N), local dimension NP
BYCOL is distributed across the process columns
All process rows are assumed to contain the same value
.TP 8
BYALL   (global output) DOUBLE PRECISION global dimension( N )
local dimension (N)
BYALL is exactly duplicated on all processes
It contains the same values as BYCOL, but it is replicated
across all processes rather than being distributed

BYALL(i) = BYCOL( NUMROC(i,NB,MYROW,0,NPROW ) on the procs
whose MYROW == mod((i-1)/NB,NPROW)
.TP 8
WORK    (local workspace) DOUBLE PRECISION dimension (LWORK)
Used to hold the buffers sent from one process to another
.TP 8
LWORK   (local input) INTEGER size of WORK array LWORK >=
LWORK >= NUMROC(N, DESC( NB_ ), 0, 0, NPCOL)
