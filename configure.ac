AC_INIT([influenceR], 0.9)

AC_PROG_CC

AC_LANG(C)

AC_ARG_WITH([mpi],
    AC_HELP_STRING([--with-mpi], [Enable MPI (must have the libraries installed and on the LD_LIBRARY_PATH variable)]),
    AC_CHECK_LIB(mpi, MPI_Init, [AC_SUBST(PKG_LIBS_MPI,-lmpi\ -lutil\ -lpthread) AC_SUBST(USE_MPI,-DUSE_MPI)], [AC_MSG_ERROR("Unable to find MPI library")]))

m4_include([ax_openmp.m4])
AX_OPENMP([AC_SUBST(SNAP_OPENMP,\--enable-openmp)])

AM_CONDITIONAL([ICC], [test x$CC = xicc])
AM_COND_IF([ICC], AC_SUBST(LRT_LIBS,"-lrt"))

AC_CONFIG_FILES([src/Makevars R/graph_metrics.R])

AC_OUTPUT
