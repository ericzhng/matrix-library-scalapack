.TH PCLAQGE l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCLAQGE - equilibrate a general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) using the row and scaling factors in the vectors R and C
.SH SYNOPSIS
.TP 20
SUBROUTINE PCLAQGE(
M, N, A, IA, JA, DESCA, R, C, ROWCND, COLCND,
AMAX, EQUED )
.TP 20
.ti +4
CHARACTER
EQUED
.TP 20
.ti +4
INTEGER
IA, JA, M, N
.TP 20
.ti +4
REAL
AMAX, COLCND, ROWCND
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
C( * ), R( * )
.TP 20
.ti +4
COMPLEX
A( * )
.SH PURPOSE
PCLAQGE equilibrates a general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) using the row and scaling factors in the vectors R and C. 
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
The number of rows to be operated on i.e the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) COMPLEX pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1))
containing on entry the M-by-N matrix sub( A ). On exit,
the equilibrated distributed matrix.  See EQUED for the
form of the equilibrated distributed submatrix.
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
R       (local input) REAL array, dimension LOCr(M_A)
The row scale factors for sub( A ). R is aligned with the
distributed matrix A, and replicated across every process
column. R is tied to the distributed matrix A.
.TP 8
C       (local input) REAL array, dimension LOCc(N_A)
The column scale factors of sub( A ). C is aligned with the
distributed matrix A, and replicated down every process
row. C is tied to the distributed matrix A.
.TP 8
ROWCND  (global input) REAL
The global ratio of the smallest R(i) to the largest R(i),
IA <= i <= IA+M-1.
.TP 8
COLCND  (global input) REAL
The global ratio of the smallest C(i) to the largest C(i),
JA <= j <= JA+N-1.
.TP 8
AMAX    (global input) REAL
Absolute value of largest distributed submatrix entry.
.TP 8
EQUED   (global output) CHARACTER
Specifies the form of equilibration that was done.
= 'N':  No equilibration
.br
= 'R':  Row equilibration, i.e., sub( A ) has been pre-
.br
multiplied by diag(R(IA:IA+M-1)),
.br
= 'C':  Column equilibration, i.e., sub( A ) has been post-
.br
multiplied by diag(C(JA:JA+N-1)),
.br
= 'B':  Both row and column equilibration, i.e., sub( A )
has been replaced by
diag(R(IA:IA+M-1)) * sub( A ) * diag(C(JA:JA+N-1)).
.SH PARAMETERS

THRESH is a threshold value used to decide if row or column scaling
should be done based on the ratio of the row or column scaling
factors.  If ROWCND < THRESH, row scaling is done, and if
COLCND < THRESH, column scaling is done.

LARGE and SMALL are threshold values used to decide if row scaling
should be done based on the absolute size of the largest matrix
element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.
