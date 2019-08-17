.TH PDLASE2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLASE2 - initialize an M-by-N distributed matrix sub( A ) denoting A(IA:IA+M-1,JA:JA+N-1) to BETA on the diagonal and ALPHA on the offdiagonals
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLASE2(
UPLO, M, N, ALPHA, BETA, A, IA, JA, DESCA )
.TP 20
.ti +4
CHARACTER
UPLO
.TP 20
.ti +4
INTEGER
IA, JA, M, N
.TP 20
.ti +4
DOUBLE
PRECISION ALPHA, BETA
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * )
.SH PURPOSE
PDLASE2 initializes an M-by-N distributed matrix sub( A ) denoting A(IA:IA+M-1,JA:JA+N-1) to BETA on the diagonal and ALPHA on the offdiagonals. PDLASE2 requires that only dimension of the matrix operand is distributed.
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
UPLO    (global input) CHARACTER
Specifies the part of the distributed matrix sub( A ) to be
set:
.br
= 'U':      Upper triangular part is set; the strictly lower
triangular part of sub( A ) is not changed;
= 'L':      Lower triangular part is set; the strictly upper
triangular part of sub( A ) is not changed;
Otherwise:  All of the matrix sub( A ) is set.
.TP 8
M       (global input) INTEGER
The number of rows to be operated on i.e the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
ALPHA   (global input) DOUBLE PRECISION
The constant to which the offdiagonal elements are to be
set.
.TP 8
BETA    (global input) DOUBLE PRECISION
The constant to which the diagonal elements are to be set.
.TP 8
A       (local output) DOUBLE PRECISION pointer into the local memory
to an array of dimension (LLD_A,LOCc(JA+N-1)).  This array
contains the local pieces of the distributed matrix sub( A )
to be set.  On exit, the leading M-by-N submatrix sub( A )
is set as follows:

if UPLO = 'U', A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=j-1, 1<=j<=N,
if UPLO = 'L', A(IA+i-1,JA+j-1) = ALPHA, j+1<=i<=M, 1<=j<=N,
otherwise,     A(IA+i-1,JA+j-1) = ALPHA, 1<=i<=M, 1<=j<=N,
IA+i.NE.JA+j,
and, for all UPLO, A(IA+i-1,JA+i-1) = BETA, 1<=i<=min(M,N).
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
