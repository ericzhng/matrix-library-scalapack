# Matrix Library (ScaLAPACK)

You can find more details in [this website here](https://ericzhng.github.io/docs-matrix-library/). This repository is only for ScaLAPACK.

## Introduction

ScaLAPACK, or Scalable LAPACK, is a library of high performance linear
algebra routines for distributed memory computers supporting MPI. It is originally hosted here: (http://www.netlib.org/scalapack/)

The distribution contains:

                               SCALAPACK/
 _____________________________/     |   \______________________________________
/   |        |        | <library.a> | SLmake.inc  |           |           |    \
    |        |        |             |             |           |           |
    |        |        |             |             |           |           |
  BLACS/   EXAMPLE/ PBLAS/       REDIST/         SRC/      TESTING/     TOOLS/
                                                           / Input \                
                                                          / Files & \
                                                         /Executables\
                                   

NOTE:  It is assumed that the BLAS and LAPACK libraries
       are available on your machine.  These libraries are NOT included 
       with this distribution, and may be obtained at the sites below.
       http://www.netlib.org/blas/
       http://www.netlib.org/lapack/

	   If BLAS and LAPACK libraries are not provided, ScaLAPACK examples would not work correctly.

## Installation

First you need to compile the code to be installed. Type `make all` to compile.

Make sure gcc/gfortran is installed through `sudo apt install gfortran build-essential`. Also mpif90 is needed, it can be installed using openmpi or mpich packages. MPI can be obtained through following websites.
       http://www.mcs.anl.gov/mpi/mpich/
       http://www.lam-mpi.org/
       http://www.open-mpi.org/

To install ScaLAPACK, simply type `make install PREFIX=/opt/apps/scalapack/2.0.0` to copy all necessary files to destination. The specified directory can be changed to anywhere you want to install ScaLAPACK to. If it is system directories, use `sudo` ahead.

Both the static and dynamic libraries are built and can be linked. Some external tests will be provided to show how to use ScaLAPACK library directly.
