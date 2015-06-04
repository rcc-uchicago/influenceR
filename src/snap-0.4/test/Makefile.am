AM_CPPFLAGS = -I$(top_srcdir)/include
AM_LDFLAGS = -L$(top_builddir)/lib -lm

bin_PROGRAMS = eval_BFS eval_biconnected_components eval_modularity_spectral \
eval_vertex_cover eval_vertex_betweenness eval_modularity_betweenness \
eval_modularity_greedy_agglomerative save_graph eval_seed_community_detection \
gen_seedfile_random_walk

EXTRA_DIST = test_graph_dimacs_5_8.dim test_graph_karate.gml \
test_graph_rand_100K_1M_undir_unwt.rnd test_graph_rand_1K_10K_undir_unwt.rnd \
test_graph_rmat_100K_1M_undir_unwt.rmat test_graph_snap_5_8_undir_unwt.gr

eval_vertex_betweenness_SOURCES = drive_vertex_betweenness.c
eval_vertex_betweenness_LDADD = $(top_builddir)/lib/libsnap.la

eval_BFS_SOURCES = drive_BFS.c
eval_BFS_LDADD = $(top_builddir)/lib/libsnap.la

eval_vertex_cover_SOURCES = drive_vertex_cover.c
eval_vertex_cover_LDADD = $(top_builddir)/lib/libsnap.la

eval_biconnected_components_SOURCES = drive_biconnected_components.c
eval_biconnected_components_LDADD   = $(top_builddir)/lib/libsnap.la

eval_modularity_spectral_SOURCES = drive_modularity_spectral.c
eval_modularity_spectral_LDADD = $(top_builddir)/lib/libsnap.la

eval_modularity_betweenness_SOURCES = drive_modularity_edge_betweenness.c
eval_modularity_betweenness_LDADD = $(top_builddir)/lib/libsnap.la

eval_modularity_greedy_agglomerative_SOURCES = drive_modularity_greedy_agglomerative.c
eval_modularity_greedy_agglomerative_LDADD = $(top_builddir)/lib/libsnap.la

save_graph_SOURCES = drive_save_graph.c
save_graph_LDADD = $(top_builddir)/lib/libsnap.la

gen_seedfile_random_walk_SOURCES = drive_random_walk.c
gen_seedfile_random_walk_LDADD = $(top_builddir)/lib/libsnap.la

eval_seed_community_detection_SOURCES = drive_seed_community_detection.c
eval_seed_community_detection_LDADD = $(top_builddir)/lib/libsnap.la

# Not building this code for now.
#eval_dyn_ds_SOURCES = test_dyn_ds.c
#eval_dyn_ds_LDADD = $(top_builddir)/lib/libsnap.la

