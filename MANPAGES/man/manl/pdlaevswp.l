.TH PDLAEVSWP l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLAEVSWP - move the eigenvectors (potentially unsorted) from where they are computed, to a ScaLAPACK standard block cyclic array, sorted so that the corresponding eigenvalues are sorted
.SH SYNOPSIS
.TP 22
SUBROUTINE PDLAEVSWP(
N, ZIN, LDZI, Z, IZ, JZ, DESCZ, NVS, KEY,
WORK, LWORK )
.TP 22
.ti +4
INTEGER
IZ, JZ, LDZI, LWORK, N
.TP 22
.ti +4
INTEGER
DESCZ( * ), KEY( * ), NVS( * )
.TP 22
.ti +4
DOUBLE
PRECISION WORK( * ), Z( * ), ZIN( LDZI, * )
.SH PURPOSE
PDLAEVSWP moves the eigenvectors (potentially unsorted) from where they are computed, to a ScaLAPACK standard block cyclic array, sorted so that the corresponding eigenvalues are sorted. 
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

NP = the number of rows local to a given process.
NQ = the number of columns local to a given process.
.TP 8
N       (global input) INTEGER
The order of the matrix A.  N >= 0.
.TP 8
ZIN     (local input) DOUBLE PRECISION array,
dimension ( LDZI, NVS(iam) )
The eigenvectors on input.  Each eigenvector resides entirely
in one process.  Each process holds a contiguous set of
NVS(iam) eigenvectors.  The first eigenvector which the
process holds is:  sum for i=[0,iam-1) of NVS(i)
.TP 8
LDZI    (locl input) INTEGER
leading dimension of the ZIN array
.TP 8
Z       (local output) DOUBLE PRECISION array
global dimension (N, N), local dimension (DESCZ(DLEN_), NQ)
The eigenvectors on output.  The eigenvectors are distributed
in a block cyclic manner in both dimensions, with a
block size of NB.
.TP 8
IZ      (global input) INTEGER
Z's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JZ      (global input) INTEGER
Z's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCZ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Z.
.TP 8
NVS     (global input) INTEGER array, dimension( nprocs+1 )
nvs(i) = number of processes
number of eigenvectors held by processes [0,i-1)
nvs(1) = number of eigen vectors held by [0,1-1) == 0
nvs(nprocs+1) = number of eigen vectors held by [0,nprocs) ==
total number of eigenvectors
.TP 8
KEY     (global input) INTEGER array, dimension( N )
Indicates the actual index (after sorting) for each of the
eigenvectors.
.TP 8
WORK    (local workspace) DOUBLE PRECISION array, dimension (LWORK)
.TP 8
LWORK   (local input) INTEGER dimension of WORK
