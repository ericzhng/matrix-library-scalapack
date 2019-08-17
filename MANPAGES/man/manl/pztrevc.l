.TH PZTREVC l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PZTREVC - compute some or all of the right and/or left eigenvectors of a complex upper triangular matrix T in parallel
.SH SYNOPSIS
.TP 20
SUBROUTINE PZTREVC(
SIDE, HOWMNY, SELECT, N, T, DESCT, VL, DESCVL,
VR, DESCVR, MM, M, WORK, RWORK, INFO )
.TP 20
.ti +4
CHARACTER
HOWMNY, SIDE
.TP 20
.ti +4
INTEGER
INFO, M, MM, N
.TP 20
.ti +4
LOGICAL
SELECT( * )
.TP 20
.ti +4
INTEGER
DESCT( * ), DESCVL( * ), DESCVR( * )
.TP 20
.ti +4
DOUBLE
PRECISION RWORK( * )
.TP 20
.ti +4
COMPLEX*16
T( * ), VL( * ), VR( * ), WORK( * )
.SH PURPOSE
PZTREVC computes some or all of the right and/or left eigenvectors of a complex upper triangular matrix T in parallel. 
The right eigenvector x and the left eigenvector y of T corresponding
to an eigenvalue w are defined by:
.br

             T*x = w*x,     y'*T = w*y'
.br

where y' denotes the conjugate transpose of the vector y.

If all eigenvectors are requested, the routine may either return the
matrices X and/or Y of right or left eigenvectors of T, or the
products Q*X and/or Q*Y, where Q is an input unitary
.br
matrix. If T was obtained from the Schur factorization of an
original matrix A = Q*T*Q', then Q*X and Q*Y are the matrices of
right or left eigenvectors of A.
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
and assume that its process grid has dimension r x c.
.br
LOCr( K ) denotes the number of elements of K that a process
would receive if K were distributed over the r processes of its
process column.
.br
Similarly, LOCc( K ) denotes the number of elements of K that a
process would receive if K were distributed over the c processes of
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
SIDE    (global input) CHARACTER*1
= 'R':  compute right eigenvectors only;
.br
= 'L':  compute left eigenvectors only;
.br
= 'B':  compute both right and left eigenvectors.
.TP 8
HOWMNY  (global input) CHARACTER*1
.br
= 'A':  compute all right and/or left eigenvectors;
.br
= 'B':  compute all right and/or left eigenvectors,
and backtransform them using the input matrices
supplied in VR and/or VL;
= 'S':  compute selected right and/or left eigenvectors,
specified by the logical array SELECT.
.TP 8
SELECT  (global input) LOGICAL array, dimension (N)
If HOWMNY = 'S', SELECT specifies the eigenvectors to be
computed.
If HOWMNY = 'A' or 'B', SELECT is not referenced.
To select the eigenvector corresponding to the j-th
eigenvalue, SELECT(j) must be set to .TRUE..
.TP 8
N       (global input) INTEGER
The order of the matrix T. N >= 0.
.TP 8
T       (global input/output) COMPLEX*16 array, dimension
(DESCT(LLD_),*)
The upper triangular matrix T.  T is modified, but restored
on exit.
.TP 8
DESCT   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix T.
.TP 8
VL      (global input/output) COMPLEX*16 array, dimension
(DESCVL(LLD_),MM)
On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
contain an N-by-N matrix Q (usually the unitary matrix Q of
Schur vectors returned by ZHSEQR).
On exit, if SIDE = 'L' or 'B', VL contains:
if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
if HOWMNY = 'B', the matrix Q*Y;
if HOWMNY = 'S', the left eigenvectors of T specified by
SELECT, stored consecutively in the columns
of VL, in the same order as their
eigenvalues.
If SIDE = 'R', VL is not referenced.
.TP 8
DESCVL  (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix VL.
.TP 8
VR      (global input/output) COMPLEX*16 array, dimension
(DESCVR(LLD_),MM)
On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
contain an N-by-N matrix Q (usually the unitary matrix Q of
Schur vectors returned by ZHSEQR).
On exit, if SIDE = 'R' or 'B', VR contains:
if HOWMNY = 'A', the matrix X of right eigenvectors of T;
if HOWMNY = 'B', the matrix Q*X;
if HOWMNY = 'S', the right eigenvectors of T specified by
SELECT, stored consecutively in the columns
of VR, in the same order as their
eigenvalues.
If SIDE = 'L', VR is not referenced.
.TP 8
DESCVR  (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix VR.
.TP 8
MM      (global input) INTEGER
The number of columns in the arrays VL and/or VR. MM >= M.
.TP 8
M       (global output) INTEGER
The number of columns in the arrays VL and/or VR actually
used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
is set to N.  Each selected eigenvector occupies one
column.
.TP 8
WORK    (local workspace) COMPLEX*16 array,
dimension ( 2*DESCT(LLD_) )
Additional workspace may be required if PZLATTRS is updated
to use WORK.
.TP 8
RWORK   (local workspace) DOUBLE PRECISION array,
dimension ( DESCT(LLD_) )
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
.SH FURTHER DETAILS
The algorithm used in this program is basically backward (forward)
substitution.  It is the hope that scaling would be used to make the
the code robust against possible overflow.  But scaling has not yet
been implemented in PZLATTRS which is called by this routine to solve
the triangular systems.  PZLATTRS just calls PZTRSV.
.br

Each eigenvector is normalized so that the element of largest
magnitude has magnitude 1; here the magnitude of a complex number
(x,y) is taken to be |x| + |y|.
.br

.SH FURTHER DETAILS
Implemented by Mark R. Fahey, June, 2000
.br

