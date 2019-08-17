.TH PSLAPV2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAPV2 - applie either P (permutation matrix indicated by IPIV) or inv( P ) to a M-by-N distributed matrix sub( A ) denoting A(IA:IA+M-1,JA:JA+N-1), resulting in row or column pivoting
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLAPV2(
DIREC, ROWCOL, M, N, A, IA, JA, DESCA, IPIV,
IP, JP, DESCIP )
.TP 20
.ti +4
CHARACTER
DIREC, ROWCOL
.TP 20
.ti +4
INTEGER
IA, IP, JA, JP, M, N
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCIP( * ), IPIV( * )
.TP 20
.ti +4
REAL
A( * )
.SH PURPOSE
PSLAPV2 applies either P (permutation matrix indicated by IPIV) or inv( P ) to a M-by-N distributed matrix sub( A ) denoting A(IA:IA+M-1,JA:JA+N-1), resulting in row or column pivoting. The pivot vector should be aligned with the distributed matrix A.  For
pivoting the rows of sub( A ), IPIV should be distributed along a
process column and replicated over all process rows.  Similarly,
IPIV should be distributed along a process row and replicated over
all process columns for column pivoting.
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
= 'F' (Forward) Applies pivots Forward from top of matrix.
Computes P * sub( A );
= 'B' (Backward) Applies pivots Backward from bottom of
matrix. Computes inv( P ) * sub( A ).
.TP 8
ROWCOL  (global input) CHARACTER
Specifies if the rows or columns are to be permuted:
= 'R' Rows will be permuted,
= 'C' Columns will be permuted.
.TP 8
M       (global input) INTEGER
The number of rows to be operated on, i.e. the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on, i.e. the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) REAL pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
On entry, this local array contains the local pieces of the
distributed matrix sub( A ) to which the row or columns
interchanges will be applied. On exit, this array contains
the local pieces of the permuted distributed matrix.
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
IPIV    (input) INTEGER array, dimension >= LOCr(M_A)+MB_A if
ROWCOL = 'R', LOCc(N_A)+NB_A otherwise. It contains
the pivoting information. IPIV(i) is the global row (column),
local row (column) i was swapped with.  The last piece of the
array of size MB_A (resp. NB_A) is used as workspace. IPIV is
tied to the distributed matrix A.
.TP 8
IP      (global input) INTEGER
IPIV's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JP      (global input) INTEGER
IPIV's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCIP  (global and local input) INTEGER array of dimension 8
The array descriptor for the distributed matrix IPIV.
