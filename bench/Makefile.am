AUTOMAKE_OPTIONS = gnu
AM_CPPFLAGS = -I$(top_builddir)/src
LDADD = $(top_builddir)/src/libspasm.la

bin_PROGRAMS = rank dm pbm pgm pluq kernel transpose permute metis blocks
metis_LDFLAGS = -L/usr/local/lib -lmetis $(top_builddir)/src/libspasm.la
metis_CFLAGS = -I/usr/local/include

%.pbm : %.sms pbm
	./pbm < $< > $@

%.png : %.pbm
	convert $< $@
