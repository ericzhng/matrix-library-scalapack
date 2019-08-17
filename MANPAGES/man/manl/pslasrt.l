.TH PSLASRT l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLASRT - Sort the numbers in D in increasing order and the corresponding vectors in Q
.SH SYNOPSIS
.TP 20
SUBROUTINE PSLASRT(
ID, N, D, Q, IQ, JQ, DESCQ, WORK, LWORK,
IWORK, LIWORK, INFO )
.TP 20
.ti +4
CHARACTER
ID
.TP 20
.ti +4
INTEGER
INFO, IQ, JQ, LIWORK, LWORK, N
.TP 20
.ti +4
INTEGER
DESCQ( * ), IWORK( * )
.TP 20
.ti +4
REAL
D( * ), Q( * ), WORK( * )
.SH PURPOSE
PSLASRT Sort the numbers in D in increasing order and the corresponding vectors in Q. 
.SH ARGUMENTS
.TP 8
ID      (global input) CHARACTER*1
= 'I': sort D in increasing order;
.br
= 'D': sort D in decreasing order. (NOT IMPLEMENTED YET)
.TP 8
N       (global input) INTEGER
The number of columns to be operated on i.e the number of
columns of the distributed submatrix sub( Q ). N >= 0.
.TP 8
D       (global input/output) REAL array, dimmension (N)
On exit, the number in D are sorted in increasing order.
.TP 8
Q       (local input) REAL pointer into the local memory
to an array of dimension (LLD_Q, LOCc(JQ+N-1) ). This array
contains the local pieces of the distributed matrix sub( A )
to be copied from.
.TP 8
IQ      (global input) INTEGER
The row index in the global array A indicating the first
row of sub( Q ).
.TP 8
JQ      (global input) INTEGER
The column index in the global array A indicating the
first column of sub( Q ).
.TP 8
DESCQ   (global and local input) INTEGER array of dimension DLEN_.
The array descriptor for the distributed matrix A.
.TP 8
WORK    (local workspace/local output) REAL array,
dimension (LWORK)
LWORK   (local or global input) INTEGER
The dimension of the array WORK.
LWORK = MAX( N, NP * ( NB + NQ ))
where
NP = NUMROC( N, NB, MYROW, IAROW, NPROW ),
NQ = NUMROC( N, NB, MYCOL, DESCQ( CSRC_ ), NPCOL )
.TP 8
IWORK   (local workspace/local output) INTEGER array,
dimension (LIWORK)

LIWORK (local or global input) INTEGER
The dimension of the array IWORK.
LIWORK = N + 2*NB + 2*NPCOL
.TP 8
INFO    (global output) INTEGER
= 0:  successful exit
.br
< 0:  If the i-th argument is an array and the j-entry had
an illegal value, then INFO = -(i*100+j), if the i-th
argument is a scalar and had an illegal value, then
INFO = -i.
