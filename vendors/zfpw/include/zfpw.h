/**
 *  @file ZFP_Wrapper.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole Wrapper.c.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ZFP_Wrapper_H
#define _ZFP_Wrapper_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zfp_toolkit.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define ZFP_SCES 0
#define ZFP_FERR 1
#define ZFP_TERR 2

size_t computeDataSize(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
unsigned char* zfp_compression_float(float* data, int mode, double tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* outSize);
float* zfp_decompression_float(unsigned char* bytes, int mode, size_t outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

unsigned char* zfp_compression_double(double* data, int mode, double tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* outSize);
double* zfp_decompression_double(unsigned char* bytes, int mode, size_t outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

unsigned char* zfp_compression(void* data, int dataType, int mode, double tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* outSize);
void* zfp_decompression(int dataType, int mode, unsigned char* bytes, size_t outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ZFP_Wrapper_H  ----- */
