# ######################################################################
#
# Research group of H. Dal: Please use this structure of the makefile
# in order to make the transfer to C. Miehe's notebook easy.
#
# Makefile for feap_83 with gfortran compiler (Version 31/10/2011)
#
# local changes by group members concern
# 1) Set TOOLPATH, i.e. path to group's fortran library
# 2) Link individual fortran objects FOBJ
# 3) Look for interger4 (32 bit) or integer8 (64 bit standard)
#
# ######################################################################
# global tools folder on Group member's maschines
  TOOLDIR = /home/oem/a_feap_bin/serial/tools_83s

# user defined fortran objects
  FOBJ = \
         $(TOOLDIR)/slcnxd.o \
	 $(TOOLDIR)/feap83.o \
	 $(TOOLDIR)/feap_dif.o \
	 ./prtrea.o \
	 ./umacr0.o \
	 ./mat.o \
	 ./2f-triquad-PF3d.o
# Feap archive directory
  ARFEAP = $(TOOLDIR)/Feap8_3g.a
# include header files (32/64 bit)
  FINCLUDE = $(TOOLDIR)/include -I$(TOOLDIR)/include/integer8
  FF = gfortran
  FFOPTFLAG = -O2
  LDOPTIONS = -L/usr/X11R6/lib -lX11 -lm

  feap: $(FOBJ) $(ARFEAP)
	  ranlib $(ARFEAP)
	$(FF) -o feap $(FOBJ) $(ARFEAP) \
	$(LDOPTIONS)

  .f.o:
	$(FF) -c $(FFOPTFLAG) -I$(FINCLUDE) $*.f -o $*.o

  clean:
	rcsclean -q
	rm -f *.o IX* O* P*  L* M* plots/* nohup.out ./feap $(TOOLDIR)/*.o


