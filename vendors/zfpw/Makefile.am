#LDFLAGS=-fPIC -shared

AUTOMAKE_OPTIONS=foreign
include_HEADERS=include/zfp_toolkit.h include/zfpw.h include/zfp_rw.h

lib_LTLIBRARIES=libzfpw.la
libzfpw_la_CFLAGS=-I./include -I./zfp/include -I./zfp/include/zfp
libzfpw_la_LDFLAGS = -version-info  0:1:0
libzfpw_la_LIDADD=./zfp/lib/libzfp.a
libzfpw_la_SOURCES=src/zfp_toolkit.c src/zfpw.c src/zfp_rw.c

libzfpw_la_LINK= $(AM_V_CC)$(LIBTOOL) --tag=CC --mode=link $(CCLD) $(libzfpw_la_CFLAGS) -O3 $(libzfpw_la_LDFLAGS) -o $(lib_LTLIBRARIES)
