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
  Si_Trajectory.h \
  Si_Trajectoryv1.h \
  Si_TrajectoryContainer.h \
  Si_TrajectoryContainerv1.h \
  SiDetectorFrame.h \
  SiTpcBeamAlignment.h \
  SiTpcHelixFit.h \
  SiTpc_Track.h \
  SiTpc_Trackv1.h \
  SiTpc_TrackContainer.h \
  SiTpc_TrackContainerv1.h \
  SiTpcTrackMatcher.h \
  SiTrajectoryFitter.h \
  SiHitSeedData.h \
  SiHitSeedReco.h \
  SiHitSeedQA.h

ROOTDICTS = \
  Si_Chain_Dict.cc \
  Si_Chainv1_Dict.cc \
  Si_ChainContainer_Dict.cc \
  Si_ChainContainerv1_Dict.cc \
  Si_Trajectory_Dict.cc \
  Si_Trajectoryv1_Dict.cc \
  Si_TrajectoryContainer_Dict.cc \
  Si_TrajectoryContainerv1_Dict.cc \
  SiTpc_Track_Dict.cc \
  SiTpc_Trackv1_Dict.cc \
  SiTpc_TrackContainer_Dict.cc \
  SiTpc_TrackContainerv1_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = $(ROOTDICTS:.cc=_rdict.pcm)

lib_LTLIBRARIES = libsitrackreco.la

libsitrackreco_la_SOURCES = \
  Si_Chainv1.cc \
  Si_ChainContainerv1.cc \
  Si_Trajectoryv1.cc \
  Si_TrajectoryContainerv1.cc \
  SiTrajectoryFitter.cc \
  SiTpcHelixFit.cc \
  SiTpc_Trackv1.cc \
  SiTpc_TrackContainerv1.cc \
  SiTpcTrackMatcher.cc \
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
  -lg4detectors_io \
  -ltpctrackreco_io

%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

%_Dict_rdict.pcm: %_Dict.cc ;

clean-local:
	rm -f *Dict* $(BUILT_SOURCES) *.pcm
