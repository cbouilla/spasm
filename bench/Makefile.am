AUTOMAKE_OPTIONS = gnu
AM_CPPFLAGS = -I$(top_builddir)/src
LDADD = $(top_builddir)/src/libspasm.la

bin_PROGRAMS = rank dm pbm pgm pluq kernel transpose permute blocks stack blocks_rank dm_rank schur narrow

if SPASM_HAVE_METIS
bin_PROGRAMS += metis
metis_SOURCES = metis.c
metis_LDFLAGS = $(METIS_LIBS)
metis_CFLAGS = $(METIS_CFLAGS)
endif

if SPASM_HAVE_LINBOX
bin_PROGRAMS += linbox_rank_gauss linbox_rank_wm
linbox_rank_gauss_SOURCES = linbox_rank_gauss.cpp
linbox_rank_gauss_LDFLAGS = $(LINBOX_LIBS)
linbox_rank_gauss_CXXFLAGS = $(LINBOX_CXXFLAGS)

linbox_rank_wm_SOURCES = linbox_rank_wm.cpp
linbox_rank_wm_LDFLAGS = $(LINBOX_LIBS)
linbox_rank_wm_CXXFLAGS = $(LINBOX_CXXFLAGS)
endif

if SPASM_HAVE_FFLAS_FFPACK
bin_PROGRAMS += ffpack_rank ffpack_narrow
ffpack_rank_SOURCES = ffpack_rank.cpp
ffpack_rank_LDFLAGS = $(FFLAS_FFPACK_LIBS) $(GIVARO_LIBS)
ffpack_rank_CXXFLAGS = $(FFLAS_FFPACK_CXXFLAGS) $(GIVARO_CXXFLAGS)

ffpack_narrow_SOURCES = ffpack_narrow.cpp
ffpack_narrow_LDFLAGS = $(FFLAS_FFPACK_LIBS) $(GIVARO_LIBS)
ffpack_narrow_CXXFLAGS = $(FFLAS_FFPACK_CXXFLAGS) $(GIVARO_CXXFLAGS)
endif

%.pbm : %.sms pbm
	./pbm < $< > $@

%.png : %.pbm
	convert $< $@
