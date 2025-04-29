/**
 *  @file zfp_rw.c
 *  @author Sheng Di
 *  @date April, 2022
 *  @brief io interface for fortrance
 *  (C) 2022 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
# include <unistd.h>

#include "zfp_rw.h"
#include "zfpw.h"
#include "zfp_toolkit.h"


//////////////////////
float *ZFP_readFloatData(char *srcFilePath, size_t *nbEle, int *status)
{
	int state = ZFP_SCES;
	if(dataEndianType==sysEndianType)
	{
		float *daBuf = readFloatData_systemEndian(srcFilePath, nbEle, &state);
		*status = state;
		return daBuf;
	}
	else
	{
		size_t i,j;

		size_t byteLength;
		unsigned char* bytes = ZFP_readByteData(srcFilePath, &byteLength, &state);
		if(state == ZFP_FERR)
		{
			*status = ZFP_FERR;
			return NULL;
		}
		float *daBuf = (float *)malloc(byteLength);
		*nbEle = byteLength/4;

		eclfloat buf;
		for(i = 0;i<*nbEle;i++)
		{
			j = i*4;
			memcpy(buf.byte, bytes+j, 4);
			ZFP_symTransform_4bytes(buf.byte);
			daBuf[i] = buf.value;
		}
		free(bytes);
		return daBuf;
	}
}

void ZFP_writeFloatData_inBytes(float *data, size_t nbEle, char* tgtFilePath, int *status)
{
	size_t i = 0;
	int state = ZFP_SCES;
	eclfloat buf;
	unsigned char* bytes = (unsigned char*)malloc(nbEle*sizeof(float));
	for(i=0;i<nbEle;i++)
	{
		buf.value = data[i];
		bytes[i*4+0] = buf.byte[0];
		bytes[i*4+1] = buf.byte[1];
		bytes[i*4+2] = buf.byte[2];
		bytes[i*4+3] = buf.byte[3];
	}

	size_t byteLength = nbEle*sizeof(float);
	ZFP_writeByteData(bytes, byteLength, tgtFilePath, &state);
	free(bytes);
	*status = state;
}

void ZFP_writeByteData(unsigned char *bytes, size_t byteLength, char *tgtFilePath, int *status)
{
	FILE *pFile = fopen(tgtFilePath, "wb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 3\n");
        *status = ZFP_FERR;
        return;
    }

    fwrite(bytes, 1, byteLength, pFile); //write outSize bytes
    fclose(pFile);
    *status = ZFP_SCES;
}

unsigned char *ZFP_readByteData(char *srcFilePath, size_t *byteLength, int *status)
{
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = ZFP_FERR;
        return 0;
    }
	fseek(pFile, 0, SEEK_END);
    *byteLength = ftell(pFile);
    fclose(pFile);

    unsigned char *byteBuf = ( unsigned char *)malloc((*byteLength)*sizeof(unsigned char)); //sizeof(char)==1

    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = ZFP_FERR;
        return 0;
    }
    fread(byteBuf, 1, *byteLength, pFile);
    fclose(pFile);
    *status = ZFP_SCES;
    return byteBuf;
}

float *readFloatData_systemEndian(char *srcFilePath, size_t *nbEle, int *status)
{
	size_t inSize;
	FILE *pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 1\n");
        *status = ZFP_FERR;
        return NULL;
    }
	fseek(pFile, 0, SEEK_END);
    inSize = ftell(pFile);
    *nbEle = inSize/4;
    fclose(pFile);

    if(inSize<=0)
    {
		printf("Error: input file is wrong!\n");
		*status = ZFP_FERR;
	}

    float *daBuf = (float *)malloc(inSize);

    pFile = fopen(srcFilePath, "rb");
    if (pFile == NULL)
    {
        printf("Failed to open input file. 2\n");
        *status = ZFP_FERR;
        return NULL;
    }
    fread(daBuf, 4, *nbEle, pFile);
    fclose(pFile);
    *status = ZFP_SCES;
    return daBuf;
}
