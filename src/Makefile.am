AUTOMAKE_OPTIONS = gnu

bin_PROGRAMS = LUstats
LUstats_LDADD = libspasm.la

libspasm_LTLIBRARIES = libspasm.la
libspasmdir = $(libdir)
libspasm_la_SOURCES = spasm_util.c spasm_io.c spasm_triplet.c spasm_GFp.c \
        spasm_scatter.c spasm_gaxpy.c spasm_permutation.c  \
	spasm_reach.c spasm_triangular.c spasm_lu.c spasm_solutions.c

include_HEADERS = spasm.h

LUstats_SOURCES = LUstats.c

