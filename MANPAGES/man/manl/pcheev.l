.TH PCHEEV l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PCHEEV - compute selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A by calling the recommended sequence of ScaLAPACK routines
.SH SYNOPSIS
.TP 19
SUBROUTINE PCHEEV(
JOBZ, UPLO, N, A, IA, JA, DESCA, W, Z, IZ, JZ,
DESCZ, WORK, LWORK, RWORK, LRWORK, INFO )
.TP 19
.ti +4
CHARACTER
JOBZ, UPLO
.TP 19
.ti +4
INTEGER
IA, INFO, IZ, JA, JZ, LRWORK, LWORK, N
.TP 19
.ti +4
INTEGER
DESCA( * ), DESCZ( * )
.TP 19
.ti +4
REAL
RWORK( * ), W( * )
.TP 19
.ti +4
COMPLEX
A( * ), WORK( * ), Z( * )
.SH PURPOSE
PCHEEV computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A by calling the recommended sequence of ScaLAPACK routines. 
In its present form, PCHEEV assumes a homogeneous system and makes
only spot checks of the consistency of the eigenvalues across the
different processes.  Because of this, it is possible that a
heterogeneous system may return incorrect results without any error
messages.
.br

Notes
.br
=====
.br
A description vector is associated with each 2D block-cyclicly dis-
tributed matrix.  This vector stores the information required to
establish the mapping between a matrix entry and its corresponding
process and memory location.
.br

In the following comments, the character _ should be read as
"of the distributed matrix".  Let A be a generic term for any 2D
block cyclicly distributed matrix.  Its description vector is DESCA:

NOTATION        STORED IN      EXPLANATION
.br
--------------- -------------- --------------------------------------
DTYPE_A(global) DESCA( DTYPE_) The descriptor type.
.br
CTXT_A (global) DESCA( CTXT_ ) The BLACS context handle, indicating
                               the BLACS process grid A is distribu-
                               ted over. The context itself is glo-
                               bal, but the handle (the integer
                               value) may vary.
.br
M_A    (global) DESCA( M_ )    The number of rows in the distributed
                               matrix A.
.br
N_A    (global) DESCA( N_ )    The number of columns in the distri-
                               buted matrix A.
.br
MB_A   (global) DESCA( MB_ )   The blocking factor used to distribute
                               the rows of A.
.br
NB_A   (global) DESCA( NB_ )   The blocking factor used to distribute
                               the columns of A.
.br
RSRC_A (global) DESCA( RSRC_ ) The process row over which the first
                               row of the matrix A is distributed.
CSRC_A (global) DESCA( CSRC_ ) The process column over which the
                               first column of A is distributed.
LLD_A  (local)  DESCA( LLD_ )  The leading dimension of the local
                               array storing the local blocks of the
                               distributed matrix A.
.br
                               LLD_A >= MAX(1,LOCr(M_A)).

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

.SH ARGUMENTS

NP = the number of rows local to a given process.
NQ = the number of columns local to a given process.
.TP 8
JOBZ    (global input) CHARACTER*1
Specifies whether or not to compute the eigenvectors:
.br
= 'N':  Compute eigenvalues only.
.br
= 'V':  Compute eigenvalues and eigenvectors.
.TP 8
UPLO    (global input) CHARACTER*1
Specifies whether the upper or lower triangular part of the
symmetric matrix A is stored:
.br
= 'U':  Upper triangular
.br
= 'L':  Lower triangular
.TP 8
N       (global input) INTEGER
The number of rows and columns of the matrix A.  N >= 0.
.TP 8
A       (local input/workspace) block cyclic COMPLEX array,
global dimension (N, N), local dimension ( LLD_A,
LOCc(JA+N-1) )

On entry, the symmetric matrix A.  If UPLO = 'U', only the
upper triangular part of A is used to define the elements of
the symmetric matrix.  If UPLO = 'L', only the lower
triangular part of A is used to define the elements of the
symmetric matrix.

On exit, the lower triangle (if UPLO='L') or the upper
triangle (if UPLO='U') of A, including the diagonal, is
destroyed.
.TP 8
IA      (global input) INTEGER
A's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JA      (global input) INTEGER
A's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCA   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
If DESCA( CTXT_ ) is incorrect, PCHEEV cannot guarantee
correct error reporting.
.TP 8
W       (global output) REAL array, dimension (N)
On normal exit, the first M entries contain the selected
eigenvalues in ascending order.
.TP 8
Z       (local output) COMPLEX array,
global dimension (N, N),
local dimension (LLD_Z, LOCc(JZ+N-1))
If JOBZ = 'V', then on normal exit the first M columns of Z
contain the orthonormal eigenvectors of the matrix
corresponding to the selected eigenvalues.
If JOBZ = 'N', then Z is not referenced.
.TP 8
IZ      (global input) INTEGER
Z's global row index, which points to the beginning of the
submatrix which is to be operated on.
.TP 8
JZ      (global input) INTEGER
Z's global column index, which points to the beginning of
the submatrix which is to be operated on.
.TP 8
DESCZ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix Z.
DESCZ( CTXT_ ) must equal DESCA( CTXT_ )
.TP 8
WORK    (local workspace/output) COMPLEX array,
dimension (LWORK)
On output, WORK(1) returns the workspace needed to guarantee
completion.  If the input parameters are incorrect, WORK(1)
may also be incorrect.

If JOBZ='N' WORK(1) = minimal workspace for eigenvalues only.
If JOBZ='V' WORK(1) = minimal workspace required to
generate all the eigenvectors.
.TP 8
LWORK   (local input) INTEGER
See below for definitions of variables used to define LWORK.
If no eigenvectors are requested (JOBZ = 'N') then
LWORK >= MAX( NB*( NP0+1 ), 3 ) +3*N
If eigenvectors are requested (JOBZ = 'V' ) then
the amount of workspace required:
LWORK >= (NP0 + NQ0 + NB)*NB + 3*N + N^2

Variable definitions:
NB = DESCA( MB_ ) = DESCA( NB_ ) =
DESCZ( MB_ ) = DESCZ( NB_ )
NP0 = NUMROC( NN, NB, 0, 0, NPROW )
NQ0 = NUMROC( MAX( N, NB, 2 ), NB, 0, 0, NPCOL )

If LWORK = -1, the LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
size for the WORK array.  The required workspace is returned
as the first element of WORK and no error message is issued
by PXERBLA.
.TP 8
RWORK   (local workspace/output) COMPLEX array,
dimension (LRWORK)
On output RWORK(1) returns the
REAL workspace needed to
guarantee completion.  If the input parameters are incorrect,
RWORK(1) may also be incorrect.
.TP 8
LRWORK  (local input) INTEGER
Size of RWORK array.
If eigenvectors are desired (JOBZ = 'V') then
LRWORK >= 2*N + 2*N-2
If eigenvectors are not desired (JOBZ = 'N') then
LRWORK >= 2*N

If LRWORK = -1, the LRWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
size for the RWORK array.  The required workspace is returned
as the first element of RWORK and no error message is issued
by PXERBLA.
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
> 0:  If INFO = 1 through N, the i(th) eigenvalue did not
converge in CSTEQR2 after a total of 30*N iterations.
If INFO = N+1, then PCHEEV has detected heterogeneity
by finding that eigenvalues were not identical across
the process grid.  In this case, the accuracy of
the results from PCHEEV cannot be guaranteed.

Alignment requirements
======================

The distributed submatrices A(IA:*, JA:*) and C(IC:IC+M-1,JC:JC+N-1)
must verify some alignment properties, namely the following
expressions should be true:

( MB_A.EQ.NB_A.EQ.MB_Z .AND. IROFFA.EQ.IROFFZ .AND. IROFFA.EQ.0 .AND.
IAROW.EQ.IZROW )
where
IROFFA = MOD( IA-1, MB_A ) and ICOFFA = MOD( JA-1, NB_A ).
