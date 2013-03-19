# makefile for sparse supernodal LU, implemented in ANSI C, targeted for 
# shared memory machines.
#
include ../make.inc

#######################################################################
#  This makefile creates the example programs for the linear equation
#  routines in SuperLU.
#
#  The command
#       make
#  without any arguments creates all the example programs.
#  The executable files are called
#       pslinsol	
#	pslinsolx
#	psrepeat
#	psspmd
#       pdlinsol
#       pdlinsolx
#       pdrepeat
#       pdspmd
#       pclinsol
#       pclinsolx
#       pcrepeat
#       pcspmd
#       pzlinsol
#       pzlinsolx
#       pzrepeat
#       pzspmd
#
#  To remove the object files after the executable files have been
#  created, enter
#       make clean
#######################################################################

HEADER = ../SRC

all:    single double complex complex16

single:    pslinsol pslinsolx psrepeat psspmd
double:    pdlinsol pdlinsolx pdrepeat pdspmd
complex:   pclinsol pclinsolx pcrepeat pcspmd
complex16: pzlinsol pzlinsolx pzrepeat pzspmd

SLINOBJS 	= pslinsol.o
SLINXOBJS 	= pslinsolx.o
SREPOBJS	= psrepeat.o
SSPMDOBJS	= psspmd.o 

DLINOBJS	= pdlinsol.o sp_ienv.o
DLINXOBJS       = pdlinsolx.o
DREPOBJS	= pdrepeat.o
DSPMDOBJS       = pdspmd.o

CLINOBJS	= pclinsol.o
CLINXOBJS       = pclinsolx.o
CREPOBJS	= pcrepeat.o
CSPMDOBJS       = pcspmd.o

ZLINOBJS	= pzlinsol.o
ZLINXOBJS       = pzlinsolx.o
ZREPOBJS	= pzrepeat.o
ZSPMDOBJS       = pzspmd.o

pslinsol: $(SLINOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pslinsolx: $(SLINXOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SLINXOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

psrepeat: $(SREPOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SREPOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

psspmd: $(SSPMDOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(SSPMDOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pdlinsol: $(DLINOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pdlinsolx: $(DLINXOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DLINXOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pdrepeat: $(DREPOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DREPOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pdspmd: $(DSPMDOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(DSPMDOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pclinsol: $(CLINOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pclinsolx: $(CLINXOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CLINXOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pcrepeat: $(CREPOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CREPOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pcspmd: $(CSPMDOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(CSPMDOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pzlinsol: $(ZLINOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pzlinsolx: $(ZLINXOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZLINXOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pzrepeat: $(ZREPOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZREPOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

pzspmd: $(ZSPMDOBJS) ../lib/$(SUPERLULIB)
	$(LOADER) $(LOADOPTS) $(ZSPMDOBJS) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@ $(VERBOSE)

F77EXM	= f77_main.o hbcode1.o c_bridge_pdgssv.o

f77exm: $(F77EXM) ../lib/$(SUPERLULIB)
	$(FORTRAN) $(LOADOPTS) $(F77EXM) \
	../lib/$(SUPERLULIB) $(BLASLIB) $(MATHLIB) $(MPLIB) -o $@

.c.o:
	$(CC) $(CFLAGS) $(CDEFS) -I$(HEADER) -c $< $(VERBOSE)

.f.o:
	$(FORTRAN) $(FFLAGS) -c $< $(VERBOSE)

clean:	
	rm -f *.o core p*linsol p*linsolx p*repeat p*spmd f77exm

