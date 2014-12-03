# $Id: GNUmakefile 22732 2012-08-16 16:44:03Z d3y133 $


OBJ = src/c2f.o src/cart2sph.o src/cint1e.o src/cint2e.o \
      src/cint_bas.o src/fblas.o src/g1e.o \
      src/g2e.o src/misc.o src/optimizer.o src/rys_roots.o \
      src/autocode/auto_intor1.o src/autocode/auto_intor2.o \
      src/autocode/auto_intor3.o src/autocode/gaunt1.o \
      src/autocode/grad1.o src/autocode/grad2.o src/autocode/grad3.o

HEADERS = config.h src/c2f.h src/cart2sph.h src/cint1e.h src/cint2e.h \
          src/cint_bas.h src/cint_const.h src/fblas.h src/g1e.h src/g2e.h \
          src/misc.h src/optimizer.h src/rys_roots.h

LIBRARY = libnwints.a
DEFINES += -DI8

include ../../config/makefile.h
include ../../config/makelib.h

justo:	$(OBJ)
	@echo "just objects made"
