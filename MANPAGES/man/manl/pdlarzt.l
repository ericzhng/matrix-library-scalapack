.TH PDLARZT l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLARZT - form the triangular factor T of a real block reflector H of order > n, which is defined as a product of k elementary reflectors as returned by PDTZRZF
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLARZT(
DIRECT, STOREV, N, K, V, IV, JV, DESCV, TAU,
T, WORK )
.TP 20
.ti +4
CHARACTER
DIRECT, STOREV
.TP 20
.ti +4
INTEGER
IV, JV, K, N
.TP 20
.ti +4
INTEGER
DESCV( * )
.TP 20
.ti +4
DOUBLE
PRECISION TAU( * ), T( * ), V( * ), WORK( * )
.SH PURPOSE
PDLARZT forms the triangular factor T of a real block reflector H of order > n, which is defined as a product of k elementary reflectors as returned by PDTZRZF. 
If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;

If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.

If STOREV = 'C', the vector which defines the elementary reflector
H(i) is stored in the i-th column of the array V, and
.br

   H  =  I - V * T * V'
.br

If STOREV = 'R', the vector which defines the elementary reflector
H(i) is stored in the i-th row of the array V, and
.br

   H  =  I - V' * T * V
.br

Currently, only STOREV = 'R' and DIRECT = 'B' are supported.

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
DIRECT  (global input) CHARACTER
Specifies the order in which the elementary reflectors are
multiplied to form the block reflector:
.br
= 'F': H = H(1) H(2) . . . H(k) (Forward, not supported yet)
.br
= 'B': H = H(k) . . . H(2) H(1) (Backward)
.TP 8
STOREV  (global input) CHARACTER
Specifies how the vectors which define the elementary
reflectors are stored (see also Further Details):
.br
= 'R': rowwise
.TP 8
N       (global input) INTEGER
The number of meaningful entries of the block reflector H.
N >= 0.
.TP 8
K       (global input) INTEGER
The order of the triangular factor T (= the number of
elementary reflectors). 1 <= K <= MB_V (= NB_V).
.TP 8
V       (input/output) DOUBLE PRECISION pointer into the local memory
to an array of local dimension (LOCr(IV+K-1),LOCc(JV+N-1)).
The distributed matrix V contains the Householder vectors.
See further details.
.TP 8
IV      (global input) INTEGER
The row index in the global array V indicating the first
row of sub( V ).
.TP 8
JV      (global input) INTEGER
The column index in the global array V indicating the
first column of sub( V ).
.TP 8
DESCV   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix V.
.TP 8
TAU     (local input) DOUBLE PRECISION array, dimension LOCr(IV+K-1)
if INCV = M_V, and LOCc(JV+K-1) otherwise. This array
contains the Householder scalars related to the Householder
vectors.  TAU is tied to the distributed matrix V.
.TP 8
T       (local output) DOUBLE PRECISION array, dimension (MB_V,MB_V)
It contains the k-by-k triangular factor of the block
reflector associated with V. T is lower triangular.
.TP 8
WORK    (local workspace) DOUBLE PRECISION array,
dimension (K*(K-1)/2)
.SH FURTHER DETAILS
The shape of the matrix V and the storage of the vectors which define
the H(i) is best illustrated by the following example with n = 5 and
k = 3. The elements equal to 1 are not stored; the corresponding
array elements are modified but restored on exit. The rest of the
array is not used.
.br

DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':

                                            ______V_____
.br
       ( v1 v2 v3 )                        /            \
       ( v1 v2 v3 )                      ( v1 v1 v1 v1 v1 . . . . 1 )
   V = ( v1 v2 v3 )                      ( v2 v2 v2 v2 v2 . . . 1   )
       ( v1 v2 v3 )                      ( v3 v3 v3 v3 v3 . . 1     )
       ( v1 v2 v3 )
.br
          .  .  .
.br
          .  .  .
.br
          1  .  .
.br
             1  .
.br
                1
.br

DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':

                                                      ______V_____
          1                                          /            \
          .  1                           ( 1 . . . . v1 v1 v1 v1 v1 )
          .  .  1                        ( . 1 . . . v2 v2 v2 v2 v2 )
          .  .  .                        ( . . 1 . . v3 v3 v3 v3 v3 )
          .  .  .
.br
       ( v1 v2 v3 )
.br
       ( v1 v2 v3 )
.br
   V = ( v1 v2 v3 )
.br
       ( v1 v2 v3 )
.br
       ( v1 v2 v3 )
.br

