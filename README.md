# matrix-library-scalapack

ScaLAPACK, or Scalable LAPACK, is a library of high performance linear
algebra routines for distributed memory computers supporting MPI.

The complete ScaLAPACK package is freely available on netlib and
can be obtained via the World Wide Web or anonymous ftp.

    http://www.netlib.org/scalapack/

ScaLAPACK, version 2.0, includes routines for the solution of dense,
band, and tridiagonal linear systems of equations, condition estimation and
iterative refinement, for LU and Cholesky factorization, matrix inversion,
full-rank linear least squares problems, orthogonal and generalized orthogonal
factorizations, orthogonal transformation routines, reductions to upper
Hessenberg, bidiagonal and tridiagonal form, reduction of a symmetric-definite/
Hermitian-definite generalized eigenproblem to standard form, the
symmetric/Hermitian, divide-and-conquer symmetric/Hermitian, generalized
symmetric/Hermitian and the nonsymmetric eigenproblem, and the singular value
decomposition.  With the exception of the singular value decomposition,
most routines are available in four types: single precision real, double
precision real, single precision complex, and double precision complex.

New in version 2.0:

- ScaLAPACK now only supports MPI.

- The BLACS is now part of ScaLAPACK, and is compiled into the ScaLAPACK
  library.  It is no longer necessary to link against BLACS libraries.

- Building ScaLAPACK using cmake is now supported.

- New MRRR Symmetric Eigenvalue Problem routines are included:
  pssyevr, pdsyevr, pcheevr and pzheevr.

- New Nonsymmetric Eigenvalue Problem QR routines for computing
  eigenvalues of a Hessenberg matrix are included for real matrices:  
  pshseqr and pdhseqr.

Unless otherwise noted, the current scalapack.tgz on netlib contains
all available updates.

Errata for ScaLAPACK (source code and documentation) can be found at:

  http://www.netlib.org/scalapack/errata.html

ScaLAPACK example programs can be found at:

  http://www.netlib.org/scalapack/examples/

A basic example is included in the EXAMPLE directory.

The ScaLAPACK User's Guide for ScaLAPACK version 1.5 is available
from SIAM at:

 http://www.ec-securehost.com/SIAM/SE04.html
 
To view an HTML version of the Users' Guide for version 1.5, see:
 
 http://www.netlib.org/scalapack/slug/
 
A number of technical reports were written during the development of
ScaLAPACK and published as LAPACK Working Notes by the University
of Tennessee.  These working notes are available at:

 http://www.netlib.org/lapack/lawns/

All questions/comments should be directed to scalapack@cs.utk.edu.

----------------------------------------------------------------------

The Makefiles in ScaLAPACK and its support libraries assume the basic
directory structure below:

                               SCALAPACK/
 _____________________________/     |   \______________________________________
/   |        |        | <library.a> | SLmake.inc  |           |           |    \
    |        |        |             |             |           |           |
    |        |        |             |             |           |           |
  BLACS/   EXAMPLE/ PBLAS/       REDIST/         SRC/      TESTING/     TOOLS/
                                                           / Input \                
                                                          / Files & \
                                                         /Executables\
                                   

NOTE:  It is assumed that the BLAS and LAPACK libraries (and MPI)
       are available on your machine.  These libraries are NOT included 
       with this distribution, and may be obtained at the sites below.

       http://www.netlib.org/blas/
       http://www.netlib.org/lapack/
       http://www.mcs.anl.gov/mpi/mpich/
       http://www.lam-mpi.org/
       http://www.open-mpi.org/

All ScaLAPACK routines -- driver, computational, and auxiliary -- can be
found in the SRC/ directory.  Testing routines and input files can be found
in the TESTING/ directory.  All machine-specific parameters and the locations
of BLAS and LAPACK libraries are specified in the SCALAPACK/SLmake.inc file.
This include file is then referenced in all subdirectory Makefiles.  Once 
the include file has been modified, the entire installation process (including
the building of testing executables) can be performed by typing ``make''
in the top-level ScaLAPACK directory.
