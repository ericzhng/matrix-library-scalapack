.TH PDSTEBZ l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDSTEBZ - compute the eigenvalues of a symmetric tridiagonal matrix in parallel
.SH SYNOPSIS
.TP 20
SUBROUTINE PDSTEBZ(
ICTXT, RANGE, ORDER, N, VL, VU, IL, IU,
ABSTOL, D, E, M, NSPLIT, W, IBLOCK, ISPLIT,
WORK, LWORK, IWORK, LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
ORDER, RANGE
.TP 20
.ti +4
INTEGER
ICTXT, IL, INFO, IU, LIWORK, LWORK, M, N,
NSPLIT
.TP 20
.ti +4
DOUBLE
PRECISION ABSTOL, VL, VU
.TP 20
.ti +4
INTEGER
IBLOCK( * ), ISPLIT( * ), IWORK( * )
.TP 20
.ti +4
DOUBLE
PRECISION D( * ), E( * ), W( * ), WORK( * )
.SH PURPOSE
PDSTEBZ computes the eigenvalues of a symmetric tridiagonal matrix in parallel. The user may ask for all eigenvalues, all eigenvalues in the interval [VL, VU], or the eigenvalues indexed IL through IU. A
static partitioning of work is done at the beginning of PDSTEBZ which
results in all processes finding an (almost) equal number of
eigenvalues.
.br

NOTE : It is assumed that the user is on an IEEE machine. If the user
       is not on an IEEE mchine, set the compile time flag NO_IEEE
       to 1 (in SLmake.inc). The features of IEEE arithmetic that
       are needed for the "fast" Sturm Count are : (a) infinity
       arithmetic (b) the sign bit of a single precision floating
       point number is assumed be in the 32nd bit position
       (c) the sign of negative zero.
.br

See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
Matrix", Report CS41, Computer Science Dept., Stanford
.br
University, July 21, 1966.
.br

.SH ARGUMENTS
.TP 8
ICTXT   (global input) INTEGER
The BLACS context handle.
.TP 8
RANGE   (global input) CHARACTER
Specifies which eigenvalues are to be found.
= 'A': ("All")   all eigenvalues will be found.
.br
= 'V': ("Value") all eigenvalues in the interval
[VL, VU] will be found.
= 'I': ("Index") the IL-th through IU-th eigenvalues (of the
entire matrix) will be found.
.TP 8
ORDER   (global input) CHARACTER
Specifies the order in which the eigenvalues and their block
numbers are stored in W and IBLOCK.
= 'B': ("By Block") the eigenvalues will be grouped by
split-off block (see IBLOCK, ISPLIT) and
ordered from smallest to largest within
the block.
= 'E': ("Entire matrix")
the eigenvalues for the entire matrix
will be ordered from smallest to largest.
.TP 8
N       (global input) INTEGER
The order of the tridiagonal matrix T.  N >= 0.
.TP 8
VL      (global input) DOUBLE PRECISION
If RANGE='V', the lower bound of the interval to be searched
for eigenvalues.  Eigenvalues less than VL will not be
returned.  Not referenced if RANGE='A' or 'I'.
.TP 8
VU      (global input) DOUBLE PRECISION
If RANGE='V', the upper bound of the interval to be searched
for eigenvalues.  Eigenvalues greater than VU will not be
returned.  VU must be greater than VL.  Not referenced if
RANGE='A' or 'I'.
.TP 8
IL      (global input) INTEGER
If RANGE='I', the index (from smallest to largest) of the
smallest eigenvalue to be returned.  IL must be at least 1.
Not referenced if RANGE='A' or 'V'.
.TP 8
IU      (global input) INTEGER
If RANGE='I', the index (from smallest to largest) of the
largest eigenvalue to be returned.  IU must be at least IL
and no greater than N.  Not referenced if RANGE='A' or 'V'.
.TP 8
ABSTOL  (global input) DOUBLE PRECISION
The absolute tolerance for the eigenvalues.  An eigenvalue
(or cluster) is considered to be located if it has been
determined to lie in an interval whose width is ABSTOL or
less.  If ABSTOL is less than or equal to zero, then ULP*|T|
will be used, where |T| means the 1-norm of T.
Eigenvalues will be computed most accurately when ABSTOL is
set to the underflow threshold DLAMCH('U'), not zero.
Note : If eigenvectors are desired later by inverse iteration
( PDSTEIN ), ABSTOL should be set to 2*PDLAMCH('S').
.TP 8
D       (global input) DOUBLE PRECISION array, dimension (N)
The n diagonal elements of the tridiagonal matrix T.  To
avoid overflow, the matrix must be scaled so that its largest
entry is no greater than overflow**(1/2) * underflow**(1/4)
in absolute value, and for greatest accuracy, it should not
be much smaller than that.
.TP 8
E       (global input) DOUBLE PRECISION array, dimension (N-1)
The (n-1) off-diagonal elements of the tridiagonal matrix T.
To avoid overflow, the matrix must be scaled so that its
largest entry is no greater than overflow**(1/2) *
underflow**(1/4) in absolute value, and for greatest
accuracy, it should not be much smaller than that.
.TP 8
M       (global output) INTEGER
The actual number of eigenvalues found. 0 <= M <= N.
(See also the description of INFO=2)
.TP 8
NSPLIT  (global output) INTEGER
The number of diagonal blocks in the matrix T.
1 <= NSPLIT <= N.
.TP 8
W       (global output) DOUBLE PRECISION array, dimension (N)
On exit, the first M elements of W contain the eigenvalues
on all processes.
.TP 8
IBLOCK  (global output) INTEGER array, dimension (N)
At each row/column j where E(j) is zero or small, the
matrix T is considered to split into a block diagonal
matrix.  On exit IBLOCK(i) specifies which block (from 1
to the number of blocks) the eigenvalue W(i) belongs to.
NOTE:  in the (theoretically impossible) event that bisection
does not converge for some or all eigenvalues, INFO is set
to 1 and the ones for which it did not are identified by a
negative block number.
.TP 8
ISPLIT  (global output) INTEGER array, dimension (N)
The splitting points, at which T breaks up into submatrices.
The first submatrix consists of rows/columns 1 to ISPLIT(1),
the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
etc., and the NSPLIT-th consists of rows/columns
ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
(Only the first NSPLIT elements will actually be used, but
since the user cannot know a priori what value NSPLIT will
have, N words must be reserved for ISPLIT.)
.TP 8
WORK    (local workspace) DOUBLE PRECISION array,
dimension ( MAX( 5*N, 7 ) )
.TP 8
LWORK   (local input) INTEGER
size of array WORK must be >= MAX( 5*N, 7 )
If LWORK = -1, then LWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
IWORK   (local workspace) INTEGER array, dimension ( MAX( 4*N, 14 ) )
.TP 8
LIWORK  (local input) INTEGER
size of array IWORK must be >= MAX( 4*N, 14, NPROCS )
If LIWORK = -1, then LIWORK is global input and a workspace
query is assumed; the routine only calculates the minimum
and optimal size for all work arrays. Each of these
values is returned in the first entry of the corresponding
work array, and no error message is issued by PXERBLA.
.TP 8
INFO    (global output) INTEGER
= 0 :  successful exit
.br
< 0 :  if INFO = -i, the i-th argument had an illegal value
.br
> 0 :  some or all of the eigenvalues failed to converge or
.br
were not computed:
.br
= 1 : Bisection failed to converge for some eigenvalues;
these eigenvalues are flagged by a negative block
number.  The effect is that the eigenvalues may not
be as accurate as the absolute and relative
tolerances. This is generally caused by arithmetic
which is less accurate than PDLAMCH says.
= 2 : There is a mismatch between the number of
eigenvalues output and the number desired.
= 3 : RANGE='i', and the Gershgorin interval initially
used was incorrect. No eigenvalues were computed.
Probable cause: your machine has sloppy floating
point arithmetic.
Cure: Increase the PARAMETER "FUDGE", recompile,
and try again.
.SH PARAMETERS
.TP 8
RELFAC  DOUBLE PRECISION, default = 2.0
The relative tolerance.  An interval [a,b] lies within
"relative tolerance" if  b-a < RELFAC*ulp*max(|a|,|b|),
where "ulp" is the machine precision (distance from 1 to
the next larger floating point number.)
.TP 8
FUDGE   DOUBLE PRECISION, default = 2.0
A "fudge factor" to widen the Gershgorin intervals.  Ideally,
a value of 1 should work, but on machines with sloppy
arithmetic, this needs to be larger.  The default for
publicly released versions should be large enough to handle
the worst machine around.  Note that this has no effect
on the accuracy of the solution.
