.TH PSLAHRD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAHRD - reduce the first NB columns of a real general N-by-(N-K+1) distributed matrix A(IA:IA+N-1,JA:JA+N-K) so that elements below the k-th subdiagonal are zero
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLAHRD(
N, K, NB, A, IA, JA, DESCA, TAU, T, Y, IY, JY,
DESCY, WORK )
.TP 20
.ti +4
INTEGER
IA, IY, JA, JY, K, N, NB
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCY( * )
.TP 20
.ti +4
REAL
A( * ), T( * ), TAU( * ), WORK( * ), Y( * )
.SH PURPOSE
PSLAHRD reduces the first NB columns of a real general N-by-(N-K+1) distributed matrix A(IA:IA+N-1,JA:JA+N-K) so that elements below the k-th subdiagonal are zero. The reduction is performed by an orthogo- nal similarity transformation Q' * A * Q. The routine returns the
matrices V and T which determine Q as a block reflector I - V*T*V',
and also the matrix Y = A * V * T.
.br

This is an auxiliary routine called by PSGEHRD. In the following
comments sub( A ) denotes A(IA:IA+N-1,JA:JA+N-1).
.br

.SH ARGUMENTS
.TP 8
N       (global input) INTEGER
The number of rows and columns to be operated on, i.e. the
order of the distributed submatrix sub( A ).
N >= 0.
.TP 8
K       (global input) INTEGER
The offset for the reduction. Elements below the k-th
subdiagonal in the first NB columns are reduced to zero.
.TP 8
NB      (global input) INTEGER
The number of columns to be reduced.
.TP 8
A       (local input/local output) REAL pointer into
the local memory to an array of dimension (LLD_A,
LOCc(JA+N-K)). On entry, this array contains the the local
pieces of the N-by-(N-K+1) general distributed matrix
A(IA:IA+N-1,JA:JA+N-K). On exit, the elements on and above
the k-th subdiagonal in the first NB columns are overwritten
with the corresponding elements of the reduced distributed
matrix; the elements below the k-th subdiagonal, with the
array TAU, represent the matrix Q as a product of elementary
reflectors. The other columns of A(IA:IA+N-1,JA:JA+N-K) are
unchanged. See Further Details.
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
TAU     (local output) REAL array, dimension LOCc(JA+N-2)
The scalar factors of the elementary reflectors (see Further
Details). TAU is tied to the distributed matrix A.
.TP 8
T       (local output) REAL array, dimension (NB_A,NB_A)
The upper triangular matrix T.
.TP 8
Y       (local output) REAL pointer into the local memory
to an array of dimension (LLD_Y,NB_A). On exit, this array
contains the local pieces of the N-by-NB distributed
matrix Y. LLD_Y >= LOCr(IA+N-1).
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
WORK    (local workspace) REAL array, dimension (NB)
.SH FURTHER DETAILS
The matrix Q is represented as a product of nb elementary reflectors

   Q = H(1) H(2) . . . H(nb).
.br

Each H(i) has the form
.br

   H(i) = I - tau * v * v'
.br

where tau is a real scalar, and v is a real vector with
.br
v(1:i+k-1) = 0, v(i+k) = 1; v(i+k+1:n) is stored on exit in
A(ia+i+k:ia+n-1,ja+i-1), and tau in TAU(ja+i-1).
.br

The elements of the vectors v together form the (n-k+1)-by-nb matrix
V which is needed, with T and Y, to apply the transformation to the
unreduced part of the matrix, using an update of the form:
A(ia:ia+n-1,ja:ja+n-k) := (I-V*T*V')*(A(ia:ia+n-1,ja:ja+n-k)-Y*V').

The contents of A(ia:ia+n-1,ja:ja+n-k) on exit are illustrated by the
following example with n = 7, k = 3 and nb = 2:
.br

   ( a   h   a   a   a )
.br
   ( a   h   a   a   a )
.br
   ( a   h   a   a   a )
.br
   ( h   h   a   a   a )
.br
   ( v1  h   a   a   a )
.br
   ( v1  v2  a   a   a )
.br
   ( v1  v2  a   a   a )
.br

where a denotes an element of the original matrix
.br
A(ia:ia+n-1,ja:ja+n-k), h denotes a modified element of the upper
Hessenberg matrix H, and vi denotes an element of the vector
defining H(i).
.br

