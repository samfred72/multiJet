AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
	-L$(OFFLINE_MAIN)/lib64 \
	-L$(ROOTSYS)/lib

pkginclude_HEADERS = \
  DijetTreeMaker.h

ROOTDICTS = 

pcmdir = $(libdir)
nobase_dist_pcm_DATA = 

lib_LTLIBRARIES = \
  libdijettreemaker_io.la \
  libdijettreemaker.la

libdijettreemaker_io_la_SOURCES = \
  $(ROOTDICTS) 

libdijettreemaker_la_SOURCES = \	
  DijetTreeMaker.cc

libdijettreemaker_io_la_LIBADD = \
  -lphool

libdijettreemaker_la_LIBADD = \
	-lcalo_io \
	-lCLHEP \
	-lfun4all \
	-lglobalvertex_io \
	-lmbd_io \
	-lphool \
	-lSubsysReco \
	-lqautils \
	-ljetbase_io \
	-lg4dst \
	-lg4eval \
	-lg4jets \
	-lmbd_io \
	-lcentrality \
	-lcalotrigger \
  -lSubsysReco \
	-lphg4hit \
	-lffaobjects \
	-lRooFit \
	-lRooFitCore


BUILT_SOURCES = testexternals.cc	

noinst_PROGRAMS = \
  testexternals_io \
  testexternals

testexternals_io_SOURCES = testexternals.cc
testexternals_io_LDADD   = libdijettreemaker_io.la

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = libdijettreemaker.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
