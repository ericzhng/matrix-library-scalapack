.TH SLAMSH l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
SLAMSH - send multiple shifts through a small (single node) matrix to see how consecutive small subdiagonal elements are modified by subsequent shifts in an effort to maximize the number of bulges that can be sent through
.SH SYNOPSIS
.TP 18
SUBROUTINE SLAMSH
( S, LDS, NBULGE, JBLK, H, LDH, N, ULP )
.TP 18
.ti +4
INTEGER
LDS, NBULGE, JBLK, LDH, N
.TP 18
.ti +4
REAL
ULP
.TP 18
.ti +4
REAL
S(LDS,*), H(LDH,*)
.SH PURPOSE
SLAMSH sends multiple shifts through a small (single node) matrix to see how consecutive small subdiagonal elements are modified by subsequent shifts in an effort to maximize the number of bulges that can be sent through. SLAMSH should only be called when there are multiple shifts/bulges 
   (NBULGE > 1) and the first shift is starting in the middle of an
   unreduced Hessenberg matrix because of two or more consecutive small
   subdiagonal elements.
.br

.SH ARGUMENTS
.TP 8
S       (local input/output) REAL array, (LDS,*)
On entry, the matrix of shifts.  Only the 2x2 diagonal of S is
referenced.  It is assumed that S has JBLK double shifts
(size 2).
On exit, the data is rearranged in the best order for 
applying.
.TP 8
LDS     (local input) INTEGER
On entry, the leading dimension of S.  Unchanged on exit.
1 < NBULGE <= JBLK <= LDS/2
.TP 8
NBULGE  (local input/output) INTEGER
On entry, the number of bulges to send through H ( >1 ).
NBULGE should be less than the maximum determined (JBLK).
1 < NBULGE <= JBLK <= LDS/2
On exit, the maximum number of bulges that can be sent 
through.
.TP 8
JBLK    (local input) INTEGER
On entry, the number of shifts determined for S.
Unchanged on exit.
.TP 8
H       (local input/output) REAL array (LDH,N)
On entry, the local matrix to apply the shifts on.
H should be aligned so that the starting row is 2.
On exit, the data is destroyed.
.TP 8
LDS     (local input) INTEGER
On entry, the leading dimension of S.  Unchanged on exit.
.TP 8
N       (local input) INTEGER
On entry, the size of H.  If all the bulges are expected to
go through, N should be at least 4*NBULGE+2.
Otherwise, NBULGE may be reduced by this routine.
.TP 8
ULP     (local input) REAL            
On entry, machine precision
Unchanged on exit.

Implemented by:  G. Henry, May 1, 1997
