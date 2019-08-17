.TH PSLAMCH l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PSLAMCH - determine single precision machine parameters
.SH SYNOPSIS
.TP 14
REAL FUNCTION
PSLAMCH( ICTXT, CMACH )
.TP 14
.ti +4
CHARACTER
CMACH
.TP 14
.ti +4
INTEGER
ICTXT
.SH PURPOSE
PSLAMCH determines single precision machine parameters. 
.SH ARGUMENTS
.TP 8
ICTXT   (global input) INTEGER
The BLACS context handle in which the computation takes
place.
.TP 8
CMACH   (global input) CHARACTER*1
Specifies the value to be returned by PSLAMCH:
.br
= 'E' or 'e',   PSLAMCH := eps
.br
= 'S' or 's ,   PSLAMCH := sfmin
.br
= 'B' or 'b',   PSLAMCH := base
.br
= 'P' or 'p',   PSLAMCH := eps*base
.br
= 'N' or 'n',   PSLAMCH := t
.br
= 'R' or 'r',   PSLAMCH := rnd
.br
= 'M' or 'm',   PSLAMCH := emin
.br
= 'U' or 'u',   PSLAMCH := rmin
.br
= 'L' or 'l',   PSLAMCH := emax
.br
= 'O' or 'o',   PSLAMCH := rmax

where
.TP 6
eps   = relative machine precision
sfmin = safe minimum, such that 1/sfmin does not overflow
base  = base of the machine
prec  = eps*base
t     = number of (base) digits in the mantissa
rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
emin  = minimum exponent before (gradual) underflow
rmin  = underflow threshold - base**(emin-1)
emax  = largest exponent before overflow
rmax  = overflow threshold  - (base**emax)*(1-eps)
