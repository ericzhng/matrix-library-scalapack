.TH DLAPST l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
DLAPST -  Define a permutation INDX that sorts the numbers in D in increasing order (if ID = 'I') or in decreasing order (if ID = 'D' )
.SH SYNOPSIS
.TP 19
SUBROUTINE DLAPST(
ID, N, D, INDX, INFO )
.TP 19
.ti +4
CHARACTER
ID
.TP 19
.ti +4
INTEGER
INFO, N
.TP 19
.ti +4
INTEGER
INDX( * )
.TP 19
.ti +4
DOUBLE
PRECISION D( * )
.SH PURPOSE
Define a permutation INDX that sorts the numbers in D in increasing order (if ID = 'I') or in decreasing order (if ID = 'D' ). 
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
D       (input)  DOUBLE PRECISION array, dimension (N)
The array to be sorted.
.TP 8
INDX    (ouput) INTEGER array, dimension (N).
The permutation which sorts the array D.
.TP 8
INFO    (output) INTEGER
= 0:  successful exit
.br
< 0:  if INFO = -i, the i-th argument had an illegal value
