AUTOMAKE_OPTIONS=foreign
AM_CPPFLAGS = -I../zfpw/include -I../zfp/include
AM_CFLAGS = -I../zfpw/include -I../zfp/include

bin_PROGRAMS=test_zfp_fixtol test_zfp_fixratio

test_zfp_fixratio_SOURCES=test_zfp_fixratio.c
test_zfp_fixratio_LDADD=../zfpw/.libs/libzfpw.a ../zfp/lib/libzfp.a -lm
test_zfp_fixtol_SOURCES=test_zfp_fixtol.c
test_zfp_fixtol_LDADD=../zfpw/.libs/libzfpw.a ../zfp/lib/libzfp.a -lm
