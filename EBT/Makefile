# Makefile for the compilation of the EBT programs PSPM-Disease and PSPM-Disease_bif
#==============================================================================

# Default configuration flags

EBTDIR      = $(EBTPATH)/integrator
INCLUDEDIRS = . $(EBTDIR) $(EBTDIR)/Odesolvers $(EXTRAINCLUDEDIRS)
TMPDIR      = /tmp

# Generic compilation and linking flags
CC	 = gcc
CFLAGS	 = -O3 -march=native -Wall -Wpointer-arith -Wcast-qual -Wcast-align $(SPECDEFS) -DPROBLEMFILE="<PSPM-Disease.h>" $(patsubst %, -I%, $(INCLUDEDIRS))
# CFLAGS	 = -g -DDEBUG=1 -march=native -Wall -Wpointer-arith -Wcast-qual -Wcast-align $(SPECDEFS) -DPROBLEMFILE="<PSPM-Disease.h>" $(patsubst %, -I%, $(INCLUDEDIRS))
LD.c     = ld -r
LIBS     = -lm

# Names of the EBT modules linked into the program-independent library
EBTMODS  =  ebtinit ebtmain ebtcohrt ebttint ebtutils ebtstop
EBTOBJS  =  $(TMPDIR)/ebtinit.o $(TMPDIR)/ebtmain.o  $(TMPDIR)/ebtcohrt.o \
	          $(TMPDIR)/ebttint.o $(TMPDIR)/ebtutils.o $(TMPDIR)/ebtstop.o

#==============================================================================

# The dependencies of the executables

PSPM-Disease:     SPECDEFS = -DBIFURCATION=0 -DPROGRAMNAME=RESERVEHABSWITCH-QC
PSPM-Disease_bif: SPECDEFS = -DBIFURCATION=1 -DPROGRAMNAME=RESERVEHABSWITCH-QC

PSPM-Disease: PSPM-Disease.o PSPM-Disease.lib.o
	$(LINK.c) $(LDFLAGS) -o $@ PSPM-Disease.o     PSPM-Disease.lib.o     $(LIBS)

PSPM-Disease_bif: PSPM-Disease_bif.o PSPM-Disease_bif.lib.o
	$(LINK.c) $(LDFLAGS) -o $@ PSPM-Disease_bif.o PSPM-Disease_bif.lib.o $(LIBS)

#==============================================================================

# The dependencies of the problem-specific object files

PSPM-Disease.o: PSPM-Disease.c PSPM-Disease.h $(EXTRADEPENDS)
	$(COMPILE.c) -o $@ PSPM-Disease.c

PSPM-Disease_bif.o: PSPM-Disease.c PSPM-Disease.h $(EXTRADEPENDS)
	$(COMPILE.c) -o $@ PSPM-Disease.c

#==============================================================================

# The dependencies of the problem object library file

PSPM-Disease.lib.o: PSPM-Disease.h
	for I in $(EBTMODS) ; do $(COMPILE.c) -o $(TMPDIR)/$${I}.o $(EBTDIR)/$${I}.c ; done
	$(LD.c) -o $@ $(EBTOBJS)
	@${CLEANCMD}

PSPM-Disease_bif.lib.o: PSPM-Disease.h
	for I in $(EBTMODS) ; do $(COMPILE.c) -o $(TMPDIR)/$${I}.o $(EBTDIR)/$${I}.c ; done
	$(LD.c) -o $@ $(EBTOBJS)
	@${CLEANCMD}

#==============================================================================

# The dependencies of some additional targets

clean:
	@echo "Cleaning up PSPM-Disease and PSPM-Disease_bif"
	@rm -f  PSPM-Disease      PSPM-Disease.o        PSPM-Disease.lib.o     PSPM-Diseasemodule.*
	@rm -f  PSPM-Disease_bif  PSPM-Disease_bif.o    PSPM-Disease_bif.lib.o PSPM-Disease_bifmodule.*
	@for I in $(EBTMODS) ; do rm -f ./$${I}.o $(TMPDIR)/$${I}.o ; done

#==============================================================================
