/**
 *  @file zfp_rw.h
 *  @author Sheng Di
 *  @date April, 2022
 *  @brief Header file for the whole io interface.
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ZFP_RW_H
#define _ZFP_RW_H

#include <stdio.h>
#include <stdint.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif

#if defined(_MSC_VER) /* MSVC Compiler Case */
#define F_OK    0       /* Test for existence.  */
#define access _access
#endif

float *readFloatData_systemEndian(char *srcFilePath, size_t *nbEle, int *status);
unsigned char *ZFP_readByteData(char *srcFilePath, size_t *byteLength, int *status);
void ZFP_writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status);
float *ZFP_readFloatData(char *srcFilePath, size_t *nbEle, int *status);
void ZFP_writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status);

#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ZFP_RW_H  ----- */
