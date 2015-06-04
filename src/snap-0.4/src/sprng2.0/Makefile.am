AM_CPPFLAGS = -I$(top_srcdir)/include
noinst_LTLIBRARIES = libsprng.la

libsprng_la_SOURCES = primes_32.c primes_64.c cputime.c makeseed.c \
memory.c checkid.c store.c lcg.c lfg.c lcg64.c \
mlfg.c cmrg.c sprng.c checkid.h cputime.h lcg64.h lfg.h mlfg.h \
primelist_32.h primes_32.h sprng.h store.h cmrg.h int64.h \
lcg.h makeseed.h memory.h multiply.h primelist_64.h primes_64.h sprng_interface.h

SPRNG_CFLAGS = 
libsprng_la_CFLAGS = $(SPRNG_CFLAGS) 
libsprng_la_LDFLAGS = $(SPRNG_CFLAGS)
