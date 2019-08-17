.TH PDLAPIV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLAPIV - applie either P (permutation matrix indicated by IPIV) or inv( P ) to a general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1), resulting in row or column pivoting
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLAPIV(
DIREC, ROWCOL, PIVROC, M, N, A, IA, JA,
DESCA, IPIV, IP, JP, DESCIP, IWORK )
.TP 20
.ti +4
CHARACTER*1
DIREC, PIVROC, ROWCOL
.TP 20
.ti +4
INTEGER
IA, IP, JA, JP, M, N
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCIP( * ), IPIV( * ), IWORK( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * )
.SH PURPOSE
PDLAPIV applies either P (permutation matrix indicated by IPIV) or inv( P ) to a general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1), resulting in row or column pivoting. The pivot vector may be distributed across a process row or a column. The pivot vector should be aligned with the distributed
matrix A. This routine will transpose the pivot vector if necessary.
For example if the row pivots should be applied to the columns of
sub( A ), pass ROWCOL='C' and PIVROC='C'.
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

Restrictions
.br
============
.br

IPIV must always be a distributed vector (not a matrix).  Thus:
IF( ROWPIV .EQ. 'C' ) THEN
.br
   JP must be 1
.br
ELSE
.br
   IP must be 1
.br
END IF
.br

The following restrictions apply when IPIV must be transposed:
IF( ROWPIV.EQ.'C' .AND. PIVROC.EQ.'C') THEN
.br
    DESCIP(MB_) must equal DESCA(NB_)
.br
ELSE IF( ROWPIV.EQ.'R" .AND. PIVROC.EQ.'R') THEN
.br
    DESCIP(NB_) must equal DESCA(MB_)
.br
END IF
.br

.SH ARGUMENTS
.TP 8
DIREC   (global input) CHARACTER*1
Specifies in which order the permutation is applied:
= 'F' (Forward) Applies pivots Forward from top of matrix.
Computes P*sub( A ).
= 'B' (Backward) Applies pivots Backward from bottom of
matrix. Computes inv( P )*sub( A ).
.TP 8
ROWCOL  (global input) CHARACTER*1
Specifies if the rows or columns are to be permuted:
= 'R' Rows will be permuted,
= 'C' Columns will be permuted.
.TP 8
PIVROC  (global input) CHARACTER*1
Specifies whether IPIV is distributed over a process row
or column:
= 'R' IPIV distributed over a process row
= 'C' IPIV distributed over a process column
.TP 8
M       (global input) INTEGER
The number of rows to be operated on, i.e. the number of
rows of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on, i.e. the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_A, LOCc(JA+N-1)).
On entry, this array contains the local pieces of the
distributed submatrix sub( A ) to which the row or column
interchanges will be applied. On exit, the local pieces
of the permuted distributed submatrix.
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
IPIV    (local input) INTEGER array, dimension (LIPIV) where LIPIV is
when ROWCOL='R' or 'r':
>= LOCr( IA+M-1 ) + MB_A      if PIVROC='C' or 'c',
>= LOCc( M + MOD(JP-1,NB_P) ) if PIVROC='R' or 'r', and,
when ROWCOL='C' or 'c':
>= LOCr( N + MOD(IP-1,MB_P) ) if PIVROC='C' or 'c',
>= LOCc( JA+N-1 ) + NB_A      if PIVROC='R' or 'r'.
This array contains the pivoting information. IPIV(i) is the
global row (column), local row (column) i was swapped with.
When ROWCOL='R' or 'r' and PIVROC='C' or 'c', or ROWCOL='C'
or 'c' and PIVROC='R' or 'r', the last piece of this array of
size MB_A (resp. NB_A) is used as workspace. In those cases,
this array is tied to the distributed matrix A.
.TP 8
IP      (global input) INTEGER
The row index in the global array P indicating the first
row of sub( P ).
.TP 8
JP      (global input) INTEGER
The column index in the global array P indicating the
first column of sub( P ).
.TP 8
DESCIP  (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed vector IPIV.
.TP 8
IWORK   (local workspace) INTEGER array, dimension (LDW)
where LDW is equal to the workspace necessary for
transposition, and the storage of the tranposed IPIV:

Let LCM be the least common multiple of NPROW and NPCOL.
IF( ROWCOL.EQ.'R' .AND. PIVROC.EQ.'R' ) THEN
IF( NPROW.EQ.NPCOL ) THEN
LDW = LOCr( N_P + MOD(JP-1, NB_P) ) + NB_P
ELSE
LDW = LOCr( N_P + MOD(JP-1, NB_P) ) +
NB_P * CEIL( CEIL(LOCc(N_P)/NB_P) / (LCM/NPCOL) )
END IF
ELSE IF( ROWCOL.EQ.'C' .AND. PIVROC.EQ.'C' ) THEN
IF( NPROW.EQ.NPCOL ) THEN
LDW = LOCc( M_P + MOD(IP-1, MB_P) ) + MB_P
ELSE
LDW = LOCc( M_P + MOD(IP-1, MB_P) ) +
MB_P * CEIL( CEIL(LOCr(M_P)/MB_P) / (LCM/NPROW) )
END IF
ELSE
IWORK is not referenced.
END IF
