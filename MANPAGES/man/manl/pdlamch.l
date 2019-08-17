.TH PDLAMCH l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLAMCH - determine double precision machine parameters
.SH SYNOPSIS
.TP 17
DOUBLE PRECISION
FUNCTION PDLAMCH( ICTXT, CMACH )
.TP 17
.ti +4
CHARACTER
CMACH
.TP 17
.ti +4
INTEGER
ICTXT
.SH PURPOSE
PDLAMCH determines double precision machine parameters. 
.SH ARGUMENTS
.TP 8
ICTXT   (global input) INTEGER
The BLACS context handle in which the computation takes
place.
.TP 8
CMACH   (global input) CHARACTER*1
Specifies the value to be returned by PDLAMCH:
.br
= 'E' or 'e',   PDLAMCH := eps
.br
= 'S' or 's ,   PDLAMCH := sfmin
.br
= 'B' or 'b',   PDLAMCH := base
.br
= 'P' or 'p',   PDLAMCH := eps*base
.br
= 'N' or 'n',   PDLAMCH := t
.br
= 'R' or 'r',   PDLAMCH := rnd
.br
= 'M' or 'm',   PDLAMCH := emin
.br
= 'U' or 'u',   PDLAMCH := rmin
.br
= 'L' or 'l',   PDLAMCH := emax
.br
= 'O' or 'o',   PDLAMCH := rmax

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
