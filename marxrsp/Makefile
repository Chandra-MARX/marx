#---*- sh -*----------------------------------------------------------------
# ANSI C compiler.  The default is to use your CC environment variable
#---------------------------------------------------------------------------
#CC = gcc         # Uncomment this line if you want gcc!!!
#CFLAGS = -g

#---------------------------------------------------------------------------
# Where is your CFITSIO library and include file?  CFITSIO is part of IRAF.
#---------------------------------------------------------------------------
# At MIT, we use:
CFITSIO_LIB = $(ASC)/$(ARCH)/lib
CFITSIO_INC = $(ASC)/$(ARCH)/include
# Perhaps one of these combinations will work for you:
#CFITSIO_LIB = /usr/local/include
#CFITSIO_INC = /usr/local/lib
#
#CFITSIO_INC = /usr/local/src/ftools/SunOS_5.5_sparc/include
#CFITSIO_LIB = /usr/local/src/ftools/SunOS_5.5_sparc/lib

#---------------------------------------------------------------------------
# New versions of CFITSIO require additional networking libraries on some
# systems.  The most well known system with this requirement is solaris.
#---------------------------------------------------------------------------
#NETLIBS = -lsocket -lnsl  # <<<---- Solaris
NETLIBS =                # <<<---- Linux
RPATH = -Wl,-R$(CFITSIO_LIB) #<<--- Linux

#---------------------------------------------------------------------------
# There should be no need to modify anything else below this point
#---------------------------------------------------------------------------

JDSRC_ROOT = ..
OBJDIR = $(ARCH)objs
INSTALL_DIR = ../bin/$(ARCH)
CFITSIO = cfitsio

JDMATHDIR = $(JDSRC_ROOT)/jdmath/src
JDMATHLIB = -L$(JDMATHDIR)/$(ARCH)objs -ljdmath

MARXDIR = $(JDSRC_ROOT)/marx/libsrc
MARXLIB = -L$(MARXDIR)/$(ARCH)objs -lmarx

PFILEINC = $(JDSRC_ROOT)/pfile/src
PFILELIB = -L$(JDSRC_ROOT)/pfile/src/$(ARCH)objs -lpfile

JDFITSINC = -I$(JDSRC_ROOT)/jdfits/src
JDFITSLIB = -L$(JDSRC_ROOT)/jdfits/src/$(ARCH)objs -ljdfits

OBJS = $(OBJDIR)/cfits.o $(OBJDIR)/rsp.o

COMPILE_CMD = $(CC) -c $(CFLAGS) \
 -I$(JDMATHDIR) -I$(PFILEINC) -I$(MARXDIR) -I$(CFITSIO_INC)

EXECINC = -I.
MISCEXEC_OBJS = $(OBJDIR)/argcargv.o

LIBS = -L$(OBJDIR) -lcfits $(MARXLIB) $(PFILELIB) $(JDFITSLIB) $(JDMATHLIB) \
    -L$(CFITSIO_LIB) $(RPATH) -l$(CFITSIO) $(NETLIBS) -lm

all: $(OBJDIR) $(OBJDIR)/marxrsp $(OBJDIR)/rspdump $(OBJDIR)/libcfits.a

$(OBJDIR):
	-/bin/mkdir $(OBJDIR)

$(OBJDIR)/libcfits.a: $(OBJS)
	-rm -f $(OBJDIR)/libcfits.a
	ar cr $(OBJDIR)/libcfits.a $(OBJS)
	ranlib $(OBJDIR)/libcfits.a

$(OBJDIR)/marxrsp: $(OBJDIR)/libcfits.a $(OBJDIR)/marxrsp.o $(MISCEXEC_OBJS)
	$(CC) -o $(OBJDIR)/marxrsp $(OBJDIR)/marxrsp.o $(MISCEXEC_OBJS) $(LIBS)
$(OBJDIR)/marxrsp.o: marxrsp.c
	$(COMPILE_CMD) -I. marxrsp.c -o $(OBJDIR)/marxrsp.o
$(OBJDIR)/marx2fits: $(OBJDIR)/libcfits.a $(OBJDIR)/marx2fits.o $(MISCEXEC_OBJS)
	$(CC) -o $(OBJDIR)/marx2fits $(OBJDIR)/marx2fits.o $(MISCEXEC_OBJS) $(LIBS)
$(OBJDIR)/marx2fits.o: marx2fits.c
	$(COMPILE_CMD) -I. marx2fits.c -o $(OBJDIR)/marx2fits.o
$(OBJDIR)/rspdump: $(OBJDIR)/libcfits.a $(OBJDIR)/rspdump.o $(MISCEXEC_OBJS)
	$(CC) -o $(OBJDIR)/rspdump $(OBJDIR)/rspdump.o $(MISCEXEC_OBJS) $(LIBS)
$(OBJDIR)/rspdump.o: rspdump.c
	$(COMPILE_CMD) -I. rspdump.c -o $(OBJDIR)/rspdump.o
#
$(OBJDIR)/cfits.o: cfits.c
	$(COMPILE_CMD) cfits.c -o $(OBJDIR)/cfits.o
$(OBJDIR)/rsp.o: rsp.c
	$(COMPILE_CMD) rsp.c -o $(OBJDIR)/rsp.o

$(OBJDIR)/argcargv.o: argcargv.c
	$(COMPILE_CMD) argcargv.c -o $(OBJDIR)/argcargv.o

clean:
	-/bin/rm $(OBJDIR)/* *~

install: all
	cp $(OBJDIR)/marxrsp $(OBJDIR)/rspdump $(INSTALL_DIR)
