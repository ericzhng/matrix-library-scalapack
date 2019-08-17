.TH CLANV2 l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
CLANV2 - compute the Schur factorization of a complex 2-by-2 nonhermitian matrix in standard form
.SH SYNOPSIS
.TP 19
SUBROUTINE CLANV2(
A, B, C, D, RT1, RT2, CS, SN )
.TP 19
.ti +4
REAL
CS
.TP 19
.ti +4
COMPLEX
A, B, C, D, RT1, RT2, SN
.SH PURPOSE
CLANV2 computes the Schur factorization of a complex 2-by-2 nonhermitian matrix in standard form: 
     [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
.br
     [ C  D ]   [ SN  CS ] [  0  DD ] [-SN  CS ]
.br

.SH ARGUMENTS
.TP 8
A       (input/output) COMPLEX
B       (input/output) COMPLEX
C       (input/output) COMPLEX
D       (input/output) COMPLEX
On entry, the elements of the input matrix.
On exit, they are overwritten by the elements of the
standardised Schur form.
.TP 8
RT1     (output) COMPLEX
RT2     (output) COMPLEX
The two eigenvalues.
.TP 8
CS      (output) REAL
SN      (output) COMPLEX
Parameters of the rotation matrix.
.SH FURTHER DETAILS
Implemented by Mark R. Fahey, May 28, 1999
.br

