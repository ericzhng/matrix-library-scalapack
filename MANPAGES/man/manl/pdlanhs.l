.TH PDLANHS l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLANHS - return the value of the one norm, or the Frobenius norm,
.SH SYNOPSIS
.TP 17
DOUBLE PRECISION
FUNCTION PDLANHS( NORM, N, A, IA, JA, DESCA,
WORK )
.TP 17
.ti +4
CHARACTER
NORM
.TP 17
.ti +4
INTEGER
IA, JA, N
.TP 17
.ti +4
INTEGER
DESCA( * )
.TP 17
.ti +4
DOUBLE
PRECISION A( * ), WORK( * )
.SH PURPOSE
PDLANHS returns the value of the one norm, or the Frobenius norm, or the infinity norm, or the element of largest absolute value of a
Hessenberg distributed matrix sub( A ) = A(IA:IA+N-1,JA:JA+N-1).

PDLANHS returns the value
.br

   ( max(abs(A(i,j))),  NORM = 'M' or 'm' with IA <= i <= IA+N-1,
   (                                      and  JA <= j <= JA+N-1,
   (
.br
   ( norm1( sub( A ) ), NORM = '1', 'O' or 'o'
.br
   (
.br
   ( normI( sub( A ) ), NORM = 'I' or 'i'
.br
   (
.br
   ( normF( sub( A ) ), NORM = 'F', 'f', 'E' or 'e'
.br

where norm1 denotes the  one norm of a matrix (maximum column sum),
normI denotes the  infinity norm  of a matrix  (maximum row sum) and
normF denotes the  Frobenius norm of a matrix (square root of sum of
squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.

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
NORM    (global input) CHARACTER
Specifies the value to be returned in PDLANHS as described
above.
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on i.e the
number of rows and columns of the distributed submatrix
sub( A ). When N = 0, PDLANHS is set to zero. N >= 0.
.TP 8
A       (local input) DOUBLE PRECISION pointer into the local memory
to an array of dimension (LLD_A, LOCc(JA+N-1) ) containing
the local pieces of sub( A ).
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
WORK    (local workspace) DOUBLE PRECISION array dimension (LWORK)
LWORK >=   0 if NORM = 'M' or 'm' (not referenced),
Nq0 if NORM = '1', 'O' or 'o',
Mp0 if NORM = 'I' or 'i',
0 if NORM = 'F', 'f', 'E' or 'e' (not referenced),
where

IROFFA = MOD( IA-1, MB_A ), ICOFFA = MOD( JA-1, NB_A ),
IAROW = INDXG2P( IA, MB_A, MYROW, RSRC_A, NPROW ),
IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL ),
Np0 = NUMROC( N+IROFFA, MB_A, MYROW, IAROW, NPROW ),
Nq0 = NUMROC( N+ICOFFA, NB_A, MYCOL, IACOL, NPCOL ),

INDXG2P and NUMROC are ScaLAPACK tool functions; MYROW,
MYCOL, NPROW and NPCOL can be determined by calling the
subroutine BLACS_GRIDINFO.
