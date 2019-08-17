.TH PDGEBD2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDGEBD2 - reduce a real general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper or lower bidiagonal form B by an orthogonal transformation
.SH SYNOPSIS
.TP 20
SUBROUTINE PDGEBD2(
M, N, A, IA, JA, DESCA, D, E, TAUQ, TAUP,
WORK, LWORK, INFO )
.TP 20
.ti +4
INTEGER
IA, INFO, JA, LWORK, M, N
.TP 20
.ti +4
INTEGER
DESCA( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), D( * ), E( * ), TAUP( * ), TAUQ( * ),
WORK( * )
.SH PURPOSE
PDGEBD2 reduces a real general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper or lower bidiagonal form B by an orthogonal transformation: Q' * sub( A ) * P = B. 
If M >= N, B is upper bidiagonal; if M < N, B is lower bidiagonal.

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
The number of rows to be operated on, i.e. the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on, i.e. the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
A       (local input/local output) DOUBLE PRECISION pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
On entry, this array contains the local pieces of the
general distributed matrix sub( A ). On exit, if M >= N,
the diagonal and the first superdiagonal of sub( A ) are
overwritten with the upper bidiagonal matrix B; the elements
below the diagonal, with the array TAUQ, represent the
orthogonal matrix Q as a product of elementary reflectors,
and the elements above the first superdiagonal, with the
array TAUP, represent the orthogonal matrix P as a product
of elementary reflectors. If M < N, the diagonal and the
first subdiagonal are overwritten with the lower bidiagonal
matrix B; the elements below the first subdiagonal, with the
array TAUQ, represent the orthogonal matrix Q as a product of
elementary reflectors, and the elements above the diagonal,
with the array TAUP, represent the orthogonal matrix P as a
product of elementary reflectors. See Further Details.
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
D       (local output) DOUBLE PRECISION array, dimension
LOCc(JA+MIN(M,N)-1) if M >= N; LOCr(IA+MIN(M,N)-1) otherwise.
The distributed diagonal elements of the bidiagonal matrix
B: D(i) = A(i,i). D is tied to the distributed matrix A.
.TP 8
E       (local output) DOUBLE PRECISION array, dimension
LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-2) otherwise.
The distributed off-diagonal elements of the bidiagonal
distributed matrix B:
if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
E is tied to the distributed matrix A.
.TP 8
TAUQ    (local output) DOUBLE PRECISION array dimension
LOCc(JA+MIN(M,N)-1). The scalar factors of the elementary
reflectors which represent the orthogonal matrix Q. TAUQ
is tied to the distributed matrix A. See Further Details.
TAUP    (local output) DOUBLE PRECISION array, dimension
LOCr(IA+MIN(M,N)-1). The scalar factors of the elementary
reflectors which represent the orthogonal matrix P. TAUP
is tied to the distributed matrix A. See Further Details.
WORK    (local workspace/local output) DOUBLE PRECISION array,
dimension (LWORK)
On exit, WORK(1) returns the minimal and optimal LWORK.
.TP 8
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK is local input and must be at least
LWORK >= MAX( MpA0, NqA0 )

where NB = MB_A = NB_A, IROFFA = MOD( IA-1, NB )
IAROW = INDXG2P( IA, NB, MYROW, RSRC_A, NPROW ),
IACOL = INDXG2P( JA, NB, MYCOL, CSRC_A, NPCOL ),
MpA0 = NUMROC( M+IROFFA, NB, MYROW, IAROW, NPROW ),
NqA0 = NUMROC( N+IROFFA, NB, MYCOL, IACOL, NPCOL ).

INDXG2P and NUMROC are ScaLAPACK tool functions;
MYROW, MYCOL, NPROW and NPCOL can be determined by calling
the subroutine BLACS_GRIDINFO.

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
INFO    (local output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
.SH FURTHER DETAILS
The matrices Q and P are represented as products of elementary
reflectors:
.br

If m >= n,
.br

   Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)

Each H(i) and G(i) has the form:
.br

   H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'

where tauq and taup are real scalars, and v and u are real vectors;
v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in
A(ia+i:ia+m-1,ja+i-1);
.br
u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in
A(ia+i-1,ja+i+1:ja+n-1);
.br
tauq is stored in TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
.br

If m < n,
.br

   Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)

Each H(i) and G(i) has the form:
.br

   H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'

where tauq and taup are real scalars, and v and u are real vectors;
v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in
A(ia+i+1:ia+m-1,ja+i-1);
.br
u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in
A(ia+i-1,ja+i:ja+n-1);
.br
tauq is stored in TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
.br

The contents of sub( A ) on exit are illustrated by the following
examples:
.br

m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):

  (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
  (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
  (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
  (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
  (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
  (  v1  v2  v3  v4  v5 )
.br

where d and e denote diagonal and off-diagonal elements of B, vi
denotes an element of the vector defining H(i), and ui an element of
the vector defining G(i).
.br

Alignment requirements
.br
======================
.br

The distributed submatrix sub( A ) must verify some alignment proper-
ties, namely the following expressions should be true:
.br
                ( MB_A.EQ.NB_A .AND. IROFFA.EQ.ICOFFA )
.br

