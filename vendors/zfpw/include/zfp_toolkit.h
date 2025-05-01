/**
 *  @file ZFP_ByteToolkit.h
 *  @author Sheng Di
 *  @date April, 2015
 *  @brief Header file for the whole ByteToolkit.c.
 *  (C) 2015 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#ifndef _ZFP_ByteToolkit_H
#define _ZFP_ByteToolkit_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#define PATH_SEPARATOR ';'
#else
#define PATH_SEPARATOR ':'
#endif

#ifdef __cplusplus
extern "C" {
#endif
#define LITTLE_ENDIAN_DATA 0
#define BIG_ENDIAN_DATA 1 /*big_endian (ppc, max, etc.) ; little_endian (x86, x64, etc.)*/

#define LITTLE_ENDIAN_SYSTEM 0
#define BIG_ENDIAN_SYSTEM 1

#define ZFP_FLOAT 0
#define ZFP_DOUBLE 1

#define ZFP_ABS 0
#define ZFP_REL 1
#define ZFP_FXR 2

extern int sysEndianType;
extern int dataEndianType;

typedef union eclint16
{
	unsigned short usvalue;
	short svalue;
	unsigned char byte[2];
} eclint16;

typedef union eclint32
{
	int ivalue;
	unsigned int uivalue;
	unsigned char byte[4];
} eclint32;

typedef union eclint64
{
	long lvalue;
	unsigned long ulvalue;
	unsigned char byte[8];
} eclint64;


typedef union eclshort
{
        unsigned short svalue;
        unsigned char byte[2];
} eclshort;


typedef union ecldouble
{
    double value;
    unsigned long lvalue;
    unsigned char byte[8];
} ecldouble;

typedef union eclfloat
{
    float value;
    unsigned int ivalue;
    unsigned char byte[4];
} eclfloat;

void ZFP_symTransform_4bytes(unsigned char data[4]);
void ZFP_symTransform_2bytes(unsigned char data[2]);
void ZFP_symTransform_8bytes(unsigned char data[8]);

int ZFP_bytesToInt_bigEndian(unsigned char* bytes);
void ZFP_intToBytes_bigEndian(unsigned char *b, unsigned int num);
long ZFP_bytesToLong_bigEndian(unsigned char* b);
void ZFP_longToBytes_bigEndian(unsigned char *b, unsigned long num);
long ZFP_doubleToOSEndianLong(double value);
int ZFP_floatToOSEndianInt(float value);

short ZFP_bytesToShort(unsigned char* bytes);
int ZFP_bytesToInt(unsigned char* bytes);
long ZFP_bytesToLong(unsigned char* bytes);
float ZFP_bytesToFloat(unsigned char* bytes);
void ZFP_floatToBytes(unsigned char *b, float num);
double ZFP_bytesToDouble(unsigned char* bytes);
void ZFP_doubleToBytes(unsigned char *b, double num);

unsigned short* ZFP_convertByteDataToUShortArray(unsigned char* bytes, size_t byteLength);
void ZFP_convertShortArrayToBytes(short* states, size_t stateLength, unsigned char* bytes);
void ZFP_convertUShortArrayToBytes(unsigned short* states, size_t stateLength, unsigned char* bytes);
void ZFP_convertIntArrayToBytes(int* states, size_t stateLength, unsigned char* bytes);
void ZFP_convertUIntArrayToBytes(unsigned int* states, size_t stateLength, unsigned char* bytes);
void ZFP_convertLongArrayToBytes(long* states, size_t stateLength, unsigned char* bytes);
void ZFP_convertULongArrayToBytes(unsigned long* states, size_t stateLength, unsigned char* bytes);


#ifdef __cplusplus
}
#endif

#endif /* ----- #ifndef _ZFP_ByteToolkit_H  ----- */
