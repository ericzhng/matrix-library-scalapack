.TH PCLAQSY l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCLAQSY - equilibrate a symmetric distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the scaling factors in the vectors SR and SC
.SH SYNOPSIS
.TP 20
SUBROUTINE PCLAQSY(
UPLO, N, A, IA, JA, DESCA, SR, SC, SCOND,
AMAX, EQUED )
.TP 20
.ti +4
CHARACTER
EQUED, UPLO
.TP 20
.ti +4
INTEGER
IA, JA, N
.TP 20
.ti +4
REAL
AMAX, SCOND
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
REAL
SC( * ), SR( * )
.TP 20
.ti +4
COMPLEX
A( * )
.SH PURPOSE
PCLAQSY equilibrates a symmetric distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1) using the scaling factors in the vectors SR and SC. 
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
Specifies whether the upper or lower triangular part of the
symmetric distributed matrix sub( A ) is to be referenced:
.br
= 'U':  Upper triangular
.br
= 'L':  Lower triangular
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (input/output) COMPLEX pointer into the local
memory to an array of local dimension (LLD_A,LOCc(JA+N-1)).
On entry, the local pieces of the distributed symmetric
matrix sub( A ). If UPLO = 'U', the leading N-by-N upper
triangular part of sub( A ) contains the upper triangular
part of the matrix, and the strictly lower triangular part
of sub( A ) is not referenced.  If UPLO = 'L', the leading
N-by-N lower triangular part of sub( A ) contains the lower
triangular part of the matrix, and the strictly upper trian-
gular part of sub( A ) is not referenced.
On exit, if EQUED = 'Y', the equilibrated matrix:
.br
diag(SR(IA:IA+N-1)) * sub( A ) * diag(SC(JA:JA+N-1)).
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
SR      (local input) REAL array, dimension LOCr(M_A)
The scale factors for A(IA:IA+M-1,JA:JA+N-1). SR is aligned
with the distributed matrix A, and replicated across every
process column. SR is tied to the distributed matrix A.
.TP 8
SC      (local input) REAL array, dimension LOCc(N_A)
The scale factors for sub( A ). SC is aligned with the dis-
tributed matrix A, and replicated down every process row.
SC is tied to the distributed matrix A.
.TP 8
SCOND   (global input) REAL
Ratio of the smallest SR(i) (respectively SC(j)) to the
largest SR(i) (respectively SC(j)), with IA <= i <= IA+N-1
and JA <= j <= JA+N-1.
.TP 8
AMAX    (global input) REAL
Absolute value of the largest distributed submatrix entry.
.TP 8
EQUED   (output) CHARACTER*1
Specifies whether or not equilibration was done.
= 'N':  No equilibration.
.br
= 'Y':  Equilibration was done, i.e., sub( A ) has been re-
.br
placed by:
.br
diag(SR(IA:IA+N-1)) * sub( A ) * diag(SC(JA:JA+N-1)).
.SH PARAMETERS

THRESH is a threshold value used to decide if scaling should be done
based on the ratio of the scaling factors.  If SCOND < THRESH,
scaling is done.

LARGE and SMALL are threshold values used to decide if scaling should
be done based on the absolute size of the largest matrix element.
If AMAX > LARGE or AMAX < SMALL, scaling is done.
