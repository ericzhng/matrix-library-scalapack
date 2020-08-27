############################################################################
#
#  Program:         ScaLAPACK
#
#  Module:          Makefile
#
#  Purpose:         Top-level Makefile
#
#  Creation date:   March 20, 1995
#
#  Modified:        February 15, 2000
#
#  Send bug reports, comments or suggestions to scalapack@cs.utk.edu
#
############################################################################

PREFIX ?= $(HOME)/matrix-library/scalapack/2.0.0

include SLmake.inc

PRECISIONS = single double complex complex16

############################################################################
#
#  The library can be set up to include routines for any combination of the
#  four PRECISIONS.  First, modify the ARCH, ARCHFLAGS, RANLIB, F77, CC,
#  F77FLAGS, CCFLAGS, F77LOADER, CCLOADER, F77LOADFLAGS, CCLOADFLAGS and
#  CDEFS definitions in SLmake.inc to match your library archiver, compiler
#  and the options to be used.
#
#  The command
#       make
#  without any arguments creates the library of precisions defined by the
#  environment variable PRECISIONS as well as the corresponding testing
#  executables,
#       make lib
#  creates only the library,
#       make exe
#  creates only the testing executables.
#       make example
#  creates only the example
#
#  The name of the library is defined in the file called SLmake.inc and
#  is created at this directory level.
#
#  To remove the object files after the library and testing executables
#  are created, enter
#       make clean
#
############################################################################

#all: lib exe example
all: lib shareso

shareso: $(SCALAPACKLIB_SHARED)

$(SCALAPACKLIB_SHARED): lib
	$(CCLOADER) -shared -o $@ -Wl,--whole-archive libscalapack.a -Wl,--no-whole-archive

lib: blacslib toolslib pblaslib redistlib scalapacklib

exe: blacsexe pblasexe redistexe scalapackexe

clean: cleanlib cleanexe cleanexample

blacslib:
	( cd BLACS; $(MAKE) lib )

pblaslib:
	( cd PBLAS/SRC; $(MAKE) $(PRECISIONS) )

redistlib:
	( cd REDIST/SRC; $(MAKE) integer $(PRECISIONS) )

scalapacklib:
	( cd SRC; $(MAKE) $(PRECISIONS) )

toolslib:
	( cd TOOLS; $(MAKE) $(PRECISIONS) )

blacsexe:
	( cd BLACS; $(MAKE) tester )

pblasexe:
	( cd PBLAS/TESTING; $(MAKE) $(PRECISIONS) )
	( cd PBLAS/TIMING; $(MAKE) $(PRECISIONS) )

scalapackexe:
	( cd TESTING/LIN; $(MAKE) $(PRECISIONS) )
	( cd TESTING/EIG; $(MAKE) $(PRECISIONS) )

redistexe:
	( cd REDIST/TESTING; $(MAKE) integer $(PRECISIONS) )

example:
	( cd EXAMPLE; $(MAKE) $(PRECISIONS) )

cleanexe:
	( cd PBLAS/TESTING; $(MAKE) clean )
	( cd PBLAS/TIMING; $(MAKE) clean )
	( cd TESTING/LIN; $(MAKE) clean )
	( cd TESTING/EIG; $(MAKE) clean )
	( cd REDIST/TESTING; $(MAKE) clean )
	( cd BLACS/TESTING; $(MAKE) clean )
	( cd TESTING; rm -f x* )

cleanlib:
	( cd BLACS; $(MAKE) clean )
	( cd PBLAS/SRC; $(MAKE) clean )
	( cd SRC; $(MAKE) clean )
	( cd TOOLS; $(MAKE) clean )
	( cd REDIST/SRC; $(MAKE) clean )
	( rm -f $(SCALAPACKLIB) )

cleanexample:
	( cd EXAMPLE; $(MAKE) clean )

install:
	mkdir -p $(PREFIX)/include
	cp -r BLACS/SRC/*.h $(PREFIX)/include/
	cp -r PBLAS/SRC/*.h $(PREFIX)/include/
	cp -r REDIST/SRC/*.h $(PREFIX)/include/
	
	mkdir -p $(PREFIX)/lib64
	cp libscalapack.so $(PREFIX)/lib64/libscalapack.so
	cp libscalapack.a $(PREFIX)/lib64/libscalapack.a
	
	mkdir -p $(PREFIX)/share/man
	cp -r MANPAGES/man/* $(PREFIX)/share/man/
