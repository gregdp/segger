srcdir	= .
TOP	= $(srcdir)/../../../..
include $(TOP)/mk/common.make

MODULE	= Segger
MODDIR	= $(datadir)/$(MODULE)
SRCS	= __init__.py ChimeraExtension.py attributes.py axes.py \
	  extract_region_dialog.py fit_devel.py fit_dialog.py \
	  graph.py imageviewer.py mousemode.py orthoview.py promod_dialog.py \
	  quaternion.py regions.py regionsfile.py rseg_dialog.py \
	  segcmd.py segfile.py segloop_dialog.py segment_dialog.py
DATA	= volseg.png fitseg.png
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
