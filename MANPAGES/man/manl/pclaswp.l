.TH PCLASWP l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCLASWP - perform a series of row or column interchanges on the distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1)
.SH SYNOPSIS
.TP 20
SUBROUTINE PCLASWP(
DIREC, ROWCOL, N, A, IA, JA, DESCA, K1, K2,
IPIV )
.TP 20
.ti +4
CHARACTER
DIREC, ROWCOL
.TP 20
.ti +4
INTEGER
IA, JA, K1, K2, N
.TP 20
.ti +4
INTEGER
DESCA( * ), IPIV( * )
.TP 20
.ti +4
COMPLEX
A( * )
.SH PURPOSE
PCLASWP performs a series of row or column interchanges on the distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1). One interchange is initiated for each of rows or columns K1 trough K2 of
sub( A ). This routine assumes that the pivoting information has
already been broadcast along the process row or column.
.br
Also note that this routine will only work for K1-K2 being in the
same MB (or NB) block.  If you want to pivot a full matrix, use
PCLAPIV.
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
DIREC   (global input) CHARACTER
Specifies in which order the permutation is applied:
= 'F' (Forward)
= 'B' (Backward)
.TP 8
ROWCOL  (global input) CHARACTER
Specifies if the rows or columns are permuted:
= 'R' (Rows)
= 'C' (Columns)
.TP 8
N       (global input) INTEGER
If ROWCOL = 'R', the length of the rows of the distributed
matrix A(*,JA:JA+N-1) to be permuted;
If ROWCOL = 'C', the length of the columns of the distributed
matrix A(IA:IA+N-1,*) to be permuted.
.TP 8
A       (local input/local output) COMPLEX pointer into the
local memory to an array of dimension (LLD_A, * ).
On entry, this array contains the local pieces of the distri-
buted matrix to which the row/columns interchanges will be
applied. On exit the permuted distributed matrix.
.TP 8
IA      (global input) INTEGER
The row index in the global array A indicating the first
row of sub( A ).
.TP 8
JA      (global input) INTEGER
The column index in the global array A indicating the
first column of sub( A ).
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
K1      (global input) INTEGER
The first element of IPIV for which a row or column inter-
change will be done.
.TP 8
K2      (global input) INTEGER
The last element of IPIV for which a row or column inter-
change will be done.
.TP 8
IPIV    (local input) INTEGER array, dimension LOCr(M_A)+MB_A for
row pivoting and LOCc(N_A)+NB_A for column pivoting.  This
array is tied to the matrix A, IPIV(K) = L implies rows
(or columns) K and L are to be interchanged.
