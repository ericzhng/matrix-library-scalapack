.TH PDLABAD l "13 August 2001" "ScaLAPACK version 1.7" ")"
.SH NAME
PDLABAD - take as input the values computed by PDLAMCH for underflow and overflow, and returns the square root of each of these values if the log of LARGE is sufficiently large
.SH SYNOPSIS
.TP 20
SUBROUTINE PDLABAD(
ICTXT, SMALL, LARGE )
.TP 20
.ti +4
INTEGER
ICTXT
.TP 20
.ti +4
DOUBLE
PRECISION LARGE, SMALL
.SH PURPOSE
PDLABAD takes as input the values computed by PDLAMCH for underflow and overflow, and returns the square root of each of these values if the log of LARGE is sufficiently large. This subroutine is intended to identify machines with a large exponent range, such as the Crays,
and redefine the underflow and overflow limits to be the square roots
of the values computed by PDLAMCH.  This subroutine is needed because
PDLAMCH does not compensate for poor arithmetic in the upper half of
the exponent range, as is found on a Cray.
.br

In addition, this routine performs a global minimization and maximi-
zation on these values, to support heterogeneous computing networks.

.SH ARGUMENTS
.TP 8
ICTXT   (global input) INTEGER
The BLACS context handle in which the computation takes
place.
.TP 8
SMALL   (local input/local output) DOUBLE PRECISION
On entry, the underflow threshold as computed by PDLAMCH.
On exit, if LOG10(LARGE) is sufficiently large, the square
root of SMALL, otherwise unchanged.
.TP 8
LARGE   (local input/local output) DOUBLE PRECISION
On entry, the overflow threshold as computed by PDLAMCH.
On exit, if LOG10(LARGE) is sufficiently large, the square
root of LARGE, otherwise unchanged.
