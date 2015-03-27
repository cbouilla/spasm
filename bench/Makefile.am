AUTOMAKE_OPTIONS = gnu
AM_CPPFLAGS = -I$(top_builddir)/src
LDADD = $(top_builddir)/src/libspasm.la

bin_PROGRAMS = rank dm pbm pgm pluq kernel transpose permute blocks

if SPASM_HAVE_METIS
bin_PROGRAMS += metis
metis_SOURCES = metis.c
metis_LDFLAGS = $(METIS_LIBS)
metis_CFLAGS = $(METIS_CFLAGS)
endif

%.pbm : %.sms pbm
	./pbm < $< > $@

%.png : %.pbm
	convert $< $@
