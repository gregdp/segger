srcdir	= .
TOP	= $(srcdir)/../..
include $(TOP)/mk/common.make

MODULE	= Segger
MODDIR	= $(datadir)/$(MODULE)
SRCS	= *.py
DATA	= volseg.png fitseg.png mapq.png marker.gif
OBJS	= $(SRCS:.py=.pyc)
CLEAN	= $(OBJS)

all: $(OBJS)

install: all
	-mkdir -p $(MODDIR)
	cd "$(srcdir)"; $(RSYNC) $(SRCS) $(DATA) $(MODDIR)
	$(RSYNC) $(OBJS) $(MODDIR)

tarball:
	-mkdir -p temp_tarball/Segger
	$(RSYNC) GNUmakefile $(DATA) $(SRCS) temp_tarball/Segger
	cd temp_tarball ; $(GNUTAR) cf segger.tar Segger ; gzip -f segger.tar
	mv temp_tarball/segger.tar.gz .
	\rm -rf temp_tarball
