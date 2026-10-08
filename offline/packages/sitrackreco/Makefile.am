AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

# Headers installed under $(prefix)/include/sitrackreco
pkgincludedir = $(includedir)/sitrackreco
pkginclude_HEADERS = \
  Si_Chain.h \
  Si_Chainv1.h \
  Si_ChainContainer.h \
  Si_ChainContainerv1.h \
  SiHitSeedData.h \
  SiHitSeedReco.h \
  SiHitSeedQA.h

ROOTDICTS = \
  Si_Chain_Dict.cc \
  Si_Chainv1_Dict.cc \
  Si_ChainContainer_Dict.cc \
  Si_ChainContainerv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

lib_LTLIBRARIES = libsitrackreco.la

libsitrackreco_la_SOURCES = \
  Si_Chainv1.cc \
  Si_ChainContainerv1.cc \
  SiHitSeedReco.cc \
  SiHitSeedQA.cc \
  $(ROOTDICTS)

libsitrackreco_la_LIBADD = \
  -lffaobjects \
  -lphool \
  -lfun4all \
  -lSubsysReco \
  -ltrack_io \
  -ltrackbase_historic_io \
  -lmvtx \
  -lmvtx_io \
  -lintt \
  -lintt_io \
  -lg4detectors_io

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
