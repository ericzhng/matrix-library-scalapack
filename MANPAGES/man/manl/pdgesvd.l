.TH PDGESVD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDGESVD - compute the singular value decomposition (SVD) of an M-by-N matrix A, optionally computing the left and/or right singular vectors
.SH SYNOPSIS
.TP 20
SUBROUTINE PDGESVD(
JOBU, JOBVT, M, N, A, IA, JA, DESCA, S, U, IU,
JU, DESCU, VT, IVT, JVT, DESCVT, WORK, LWORK,
INFO )
.TP 20
.ti +4
CHARACTER
JOBU, JOBVT
.TP 20
.ti +4
INTEGER
IA, INFO, IU, IVT, JA, JU, JVT, LWORK, M, N
.TP 20
.ti +4
INTEGER
DESCA( * ), DESCU( * ), DESCVT( * )
.TP 20
.ti +4
DOUBLE
PRECISION A( * ), S( * ), U( * ), VT( * ), WORK( * )
.SH PURPOSE
PDGESVD computes the singular value decomposition (SVD) of an M-by-N matrix A, optionally computing the left and/or right singular vectors. The SVD is written as 
     A = U * SIGMA * transpose(V)
.br

where SIGMA is an M-by-N matrix which is zero except for its
min(M,N) diagonal elements, U is an M-by-M orthogonal matrix, and
V is an N-by-N orthogonal matrix. The diagonal elements of SIGMA
are the singular values of A and the columns of U and V are the
corresponding right and left singular vectors, respectively. The
singular values are returned in array S in decreasing order and
only the first min(M,N) columns of U and rows of VT = V**T are
computed.
.br

Notes
.br
=====
.br
Each global data object is described by an associated description
vector. This vector stores the information required to establish
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

Let K be the number of rows or columns of a distributed matrix, and
assume that its process grid has dimension r x c. LOCr( K ) denotes
the number of elements of K that a process would receive if K were
distributed over the r processes of its process column. Similarly,
LOCc( K ) denotes the number of elements of K that a process would
receive if K were distributed over the c processes of its process
row. The values of LOCr() and LOCc() may be determined via a call
to the ScaLAPACK tool function, NUMROC:
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

MP = number of local rows in A and U
NQ = number of local columns in A and VT
SIZE = min( M, N )
SIZEQ = number of local columns in U
SIZEP = number of local rows in VT
.TP 8
JOBU    (global input) CHARACTER*1
Specifies options for computing U:
.br
= 'V':  the first SIZE columns of U (the left singular
vectors) are returned in the array U;
= 'N':  no columns of U (no left singular vectors) are
computed.
.TP 8
JOBVT   (global input) CHARACTER*1
Specifies options for computing V**T:
.br
= 'V':  the first SIZE rows of V**T (the right singular
vectors) are returned in the array VT;
= 'N':  no rows of V**T (no right singular vectors) are
computed.
.TP 8
M       (global input) INTEGER
The number of rows of the input matrix A.  M >= 0.
.TP 8
N       (global input) INTEGER
The number of columns of the input matrix A.  N >= 0.
.TP 8
A       (local input/workspace) block cyclic DOUBLE PRECISION array,
global dimension (M, N), local dimension (MP, NQ)
On exit, the contents of A are destroyed.
.TP 8
IA      (global input) INTEGER
The row index in the global array A indicating the first
row of sub( A ).
.TP 8
JA      (global input) INTEGER
The column index in the global array A indicating the
first column of sub( A ).
.TP 8
DESCA   (global input) INTEGER array of dimension DLEN_
The array descriptor for the distributed matrix A.
.TP 8
S       (global output) DOUBLE PRECISION array, dimension SIZE
The singular values of A, sorted so that S(i) >= S(i+1).
.TP 8
U       (local output) DOUBLE PRECISION array, local dimension
(MP, SIZEQ), global dimension (M, SIZE)
if JOBU = 'V', U contains the first min(m,n) columns of U
if JOBU = 'N', U is not referenced.
.TP 8
IU      (global input) INTEGER
The row index in the global array U indicating the first
row of sub( U ).
.TP 8
JU      (global input) INTEGER
The column index in the global array U indicating the
first column of sub( U ).
.TP 8
DESCU   (global input) INTEGER array of dimension DLEN_
The array descriptor for the distributed matrix U.
.TP 8
VT      (local output) DOUBLE PRECISION array, local dimension
(SIZEP, NQ), global dimension (SIZE, N).
If JOBVT = 'V', VT contains the first SIZE rows of
V**T. If JOBVT = 'N', VT is not referenced.
.TP 8
IVT     (global input) INTEGER
The row index in the global array VT indicating the first
row of sub( VT ).
.TP 8
JVT     (global input) INTEGER
The column index in the global array VT indicating the
first column of sub( VT ).
.TP 9
DESCVT   (global input) INTEGER array of dimension DLEN_
The array descriptor for the distributed matrix VT.
.TP 8
WORK    (local workspace/output) DOUBLE PRECISION array, dimension
(LWORK)
On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
.TP 8
LWORK   (local input) INTEGER
The dimension of the array WORK.

LWORK > 2 + 6*SIZEB + MAX(WATOBD, WBDTOSVD),

where SIZEB = MAX(M,N), and WATOBD and WBDTOSVD refer,
respectively, to the workspace required to bidiagonalize
the matrix A and to go from the bidiagonal matrix to the
singular value decomposition U*S*VT.

For WATOBD, the following holds:

WATOBD = MAX(MAX(WPDLANGE,WPDGEBRD),
MAX(WPDLARED2D,WPDLARED1D)),

where WPDLANGE, WPDLARED1D, WPDLARED2D, WPDGEBRD are the
workspaces required respectively for the subprograms
PDLANGE, PDLARED1D, PDLARED2D, PDGEBRD. Using the
standard notation

MP = NUMROC( M, MB, MYROW, DESCA( CTXT_ ), NPROW),
NQ = NUMROC( N, NB, MYCOL, DESCA( LLD_ ), NPCOL),

the workspaces required for the above subprograms are

WPDLANGE = MP,
WPDLARED1D = NQ0,
WPDLARED2D = MP0,
WPDGEBRD = NB*(MP + NQ + 1) + NQ,

where NQ0 and MP0 refer, respectively, to the values obtained
at MYCOL = 0 and MYROW = 0. In general, the upper limit for
the workspace is given by a workspace required on
processor (0,0):

WATOBD <= NB*(MP0 + NQ0 + 1) + NQ0.

In case of a homogeneous process grid this upper limit can
be used as an estimate of the minimum workspace for every
processor.

For WBDTOSVD, the following holds:

WBDTOSVD = SIZE*(WANTU*NRU + WANTVT*NCVT) +
MAX(WDBDSQR,
MAX(WANTU*WPDORMBRQLN, WANTVT*WPDORMBRPRT)),

where

1, if left(right) singular vectors are wanted
WANTU(WANTVT) =
0, otherwise

and WDBDSQR, WPDORMBRQLN and WPDORMBRPRT refer respectively
to the workspace required for the subprograms DBDSQR,
PDORMBR(QLN), and PDORMBR(PRT), where QLN and PRT are the
values of the arguments VECT, SIDE, and TRANS in the call
to PDORMBR. NRU is equal to the local number of rows of
the matrix U when distributed 1-dimensional "column" of
processes. Analogously, NCVT is equal to the local number
of columns of the matrix VT when distributed across
1-dimensional "row" of processes. Calling the LAPACK
procedure DBDSQR requires

WDBDSQR = MAX(1, 2*SIZE + (2*SIZE - 4)*MAX(WANTU, WANTVT))

on every processor. Finally,

WPDORMBRQLN = MAX( (NB*(NB-1))/2, (SIZEQ+MP)*NB)+NB*NB,
WPDORMBRPRT = MAX( (MB*(MB-1))/2, (SIZEP+NQ)*MB )+MB*MB,

If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
size for the work array. The required workspace is returned
as the first element of WORK and no error message is issued
by PXERBLA.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit.
.br
< 0:  if INFO = -i, the i-th argument had an illegal value.
.br
> 0:  if SBDSQR did not converge
If INFO = MIN(M,N) + 1, then PDGESVD has detected
heterogeneity by finding that eigenvalues were not
identical across the process grid. In this case, the
accuracy of the results from PDGESVD cannot be
guaranteed.
