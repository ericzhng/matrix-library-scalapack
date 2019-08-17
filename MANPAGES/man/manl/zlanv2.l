.TH ZLANV2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
ZLANV2 - compute the Schur factorization of a complex 2-by-2 nonhermitian matrix in standard form
.SH SYNOPSIS
.TP 19
SUBROUTINE ZLANV2(
A, B, C, D, RT1, RT2, CS, SN )
.TP 19
.ti +4
DOUBLE
PRECISION CS
.TP 19
.ti +4
COMPLEX*16
A, B, C, D, RT1, RT2, SN
.SH PURPOSE
ZLANV2 computes the Schur factorization of a complex 2-by-2 nonhermitian matrix in standard form: 
     [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
.br
     [ C  D ]   [ SN  CS ] [  0  DD ] [-SN  CS ]
.br

.SH ARGUMENTS
.TP 8
A       (input/output) COMPLEX*16
B       (input/output) COMPLEX*16
C       (input/output) COMPLEX*16
D       (input/output) COMPLEX*16
On entry, the elements of the input matrix.
On exit, they are overwritten by the elements of the
standardised Schur form.
.TP 8
RT1     (output) COMPLEX*16
RT2     (output) COMPLEX*16
The two eigenvalues.
.TP 8
CS      (output) DOUBLE PRECISION
SN      (output) COMPLEX*16
Parameters of the rotation matrix.
.SH FURTHER DETAILS
Implemented by Mark R. Fahey, May 28, 1999
.br

