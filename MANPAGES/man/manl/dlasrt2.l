.TH DLASRT2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
DLASRT2 - the numbers in D in increasing order (if ID = 'I') or in decreasing order (if ID = 'D' )
.SH SYNOPSIS
.TP 20
SUBROUTINE DLASRT2(
ID, N, D, KEY, INFO )
.TP 20
.ti +4
CHARACTER
ID
.TP 20
.ti +4
INTEGER
INFO, N
.TP 20
.ti +4
INTEGER
KEY( * )
.TP 20
.ti +4
DOUBLE
PRECISION D( * )
.SH PURPOSE
Sort the numbers in D in increasing order (if ID = 'I') or in decreasing order (if ID = 'D' ). 
Use Quick Sort, reverting to Insertion sort on arrays of
.br
size <= 20. Dimension of STACK limits N to about 2**32.
.br

.SH ARGUMENTS
.TP 8
ID      (input) CHARACTER*1
= 'I': sort D in increasing order;
.br
= 'D': sort D in decreasing order.
.TP 8
N       (input) INTEGER
The length of the array D.
.TP 8
D       (input/output) DOUBLE PRECISION array, dimension (N)
On entry, the array to be sorted.
On exit, D has been sorted into increasing order
(D(1) <= ... <= D(N) ) or into decreasing order
(D(1) >= ... >= D(N) ), depending on ID.
.TP 8
KEY     (input/output) INTEGER array, dimension (N)
On entry, KEY contains a key to each of the entries in D()
Typically, KEY(I) = I for all I
On exit, KEY is permuted in exactly the same manner as
D() was permuted from input to output
Therefore, if KEY(I) = I for all I upon input, then
D_out(I) = D_in(KEY(I))
.TP 8
INFO    (output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
