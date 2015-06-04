AM_CPPFLAGS = -I$(top_srcdir)/include

lib_LTLIBRARIES = libsnap.la
libsnap_la_LDFLAGS = 
libsnap_la_SOURCES = 

libsnap_la_LIBADD = \
	$(top_builddir)/src/graph_generation/libgraph_gen.la \
	$(top_builddir)/src/graph_kernels/libgraph_kernels.la \
	$(top_builddir)/src/graph_metrics/libgraph_metrics.la \
	$(top_builddir)/src/graph_partitioning/libgraph_partitioning.la \
	$(top_builddir)/src/misc/libsnap_misc.la \
	$(top_builddir)/src/sprng2.0/libsprng.la
