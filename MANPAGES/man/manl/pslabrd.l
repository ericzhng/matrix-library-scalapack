.TH PSLABRD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLABRD - reduce the first NB rows and columns of a real general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper or lower bidiagonal form by an orthogonal transformation Q' * A * P,
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLABRD(
M, N, NB, A, IA, JA, DESCA, D, E, TAUQ, TAUP,
X, IX, JX, DESCX, Y, IY, JY, DESCY, WORK )
.TP 20
.ti +4
INTEGER
IA, IX, IY, JA, JX, JY, M, N, NB
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCX( * ), DESCY( * )
.TP 20
.ti +4
REAL
A( * ), D( * ), E( * ), TAUP( * ),
TAUQ( * ), X( * ), Y( * ), WORK( * )
.SH PURPOSE
PSLABRD reduces the first NB rows and columns of a real general M-by-N distributed matrix sub( A ) = A(IA:IA+M-1,JA:JA+N-1) to upper or lower bidiagonal form by an orthogonal transformation Q' * A * P, and returns the matrices X and Y which are needed to apply the
transformation to the unreduced part of sub( A ).
.br

If M >= N, sub( A ) is reduced to upper bidiagonal form; if M < N, to
lower bidiagonal form.
.br

This is an auxiliary routine called by PSGEBRD.
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
M       (global input) INTEGER
The number of rows to be operated on, i.e. the number of rows
of the distributed submatrix sub( A ). M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns to be operated on, i.e. the number of
columns of the distributed submatrix sub( A ). N >= 0.
.TP 8
NB      (global input) INTEGER
The number of leading rows and columns of sub( A ) to be
reduced.
.TP 8
A       (local input/local output) REAL pointer into the
local memory to an array of dimension (LLD_A,LOCc(JA+N-1)).
On entry, this array contains the local pieces of the
general distributed matrix sub( A ) to be reduced. On exit,
the first NB rows and columns of the matrix are overwritten;
the rest of the distributed matrix sub( A ) is unchanged.
If m >= n, elements on and below the diagonal in the first NB
columns, with the array TAUQ, represent the orthogonal
matrix Q as a product of elementary reflectors; and
elements above the diagonal in the first NB rows, with the
array TAUP, represent the orthogonal matrix P as a product
of elementary reflectors.
If m < n, elements below the diagonal in the first NB
columns, with the array TAUQ, represent the orthogonal
matrix Q as a product of elementary reflectors, and
elements on and above the diagonal in the first NB rows,
with the array TAUP, represent the orthogonal matrix P as
a product of elementary reflectors.
See Further Details.
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
D       (local output) REAL array, dimension
LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-1) otherwise.
The distributed diagonal elements of the bidiagonal matrix
B: D(i) = A(ia+i-1,ja+i-1). D is tied to the distributed
matrix A.
.TP 8
E       (local output) REAL array, dimension
LOCr(IA+MIN(M,N)-1) if M >= N; LOCc(JA+MIN(M,N)-2) otherwise.
The distributed off-diagonal elements of the bidiagonal
distributed matrix B:
if m >= n, E(i) = A(ia+i-1,ja+i) for i = 1,2,...,n-1;
if m < n, E(i) = A(ia+i,ja+i-1) for i = 1,2,...,m-1.
E is tied to the distributed matrix A.
.TP 8
TAUQ    (local output) REAL array dimension
LOCc(JA+MIN(M,N)-1). The scalar factors of the elementary
reflectors which represent the orthogonal matrix Q. TAUQ
is tied to the distributed matrix A. See Further Details.
TAUP    (local output) REAL array, dimension
LOCr(IA+MIN(M,N)-1). The scalar factors of the elementary
reflectors which represent the orthogonal matrix P. TAUP
is tied to the distributed matrix A. See Further Details.
X       (local output) REAL pointer into the local memory
to an array of dimension (LLD_X,NB). On exit, the local
pieces of the distributed M-by-NB matrix
X(IX:IX+M-1,JX:JX+NB-1) required to update the unreduced
part of sub( A ).
.TP 8
IX      (global input) INTEGER
The row index in the global array X indicating the first
row of sub( X ).
.TP 8
JX      (global input) INTEGER
The column index in the global array X indicating the
first column of sub( X ).
.TP 8
DESCX   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix X.
.TP 8
Y       (local output) REAL pointer into the local memory
to an array of dimension (LLD_Y,NB).  On exit, the local
pieces of the distributed N-by-NB matrix
Y(IY:IY+N-1,JY:JY+NB-1) required to update the unreduced
part of sub( A ).
.TP 8
IY      (global input) INTEGER
The row index in the global array Y indicating the first
row of sub( Y ).
.TP 8
JY      (global input) INTEGER
The column index in the global array Y indicating the
first column of sub( Y ).
.TP 8
DESCY   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Y.
.TP 8
WORK    (local workspace) REAL array, dimension (LWORK)
LWORK >= NB_A + NQ, with

NQ = NUMROC( N+MOD( IA-1, NB_Y ), NB_Y, MYCOL, IACOL, NPCOL )
IACOL = INDXG2P( JA, NB_A, MYCOL, CSRC_A, NPCOL )

INDXG2P and NUMROC are ScaLAPACK tool functions;
MYROW, MYCOL, NPROW and NPCOL can be determined by calling
the subroutine BLACS_GRIDINFO.
.SH FURTHER DETAILS
The matrices Q and P are represented as products of elementary
reflectors:
.br

   Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)

Each H(i) and G(i) has the form:
.br

   H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'

where tauq and taup are real scalars, and v and u are real vectors.

If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
A(ia+i-1:ia+m-1,ja+i-1); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is
stored on exit in A(ia+i-1,ja+i:ja+n-1); tauq is stored in
TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
.br

If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
A(ia+i+1:ia+m-1,ja+i-1); u(1:i-1) = 0, u(i) = 1, and u(i:n) is
stored on exit in A(ia+i-1,ja+i:ja+n-1); tauq is stored in
TAUQ(ja+i-1) and taup in TAUP(ia+i-1).
.br

The elements of the vectors v and u together form the m-by-nb matrix
V and the nb-by-n matrix U' which are needed, with X and Y, to apply
the transformation to the unreduced part of the matrix, using a block
update of the form:  sub( A ) := sub( A ) - V*Y' - X*U'.
.br

The contents of sub( A ) on exit are illustrated by the following
examples with nb = 2:
.br

m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):

  (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
  (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
  (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
  (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
  (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
  (  v1  v2  a   a   a  )
.br

where a denotes an element of the original matrix which is unchanged,
vi denotes an element of the vector defining H(i), and ui an element
of the vector defining G(i).
.br

