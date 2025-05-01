/*
 *  Copyright (C) 2021, Argonne National Laboratory
 *  See COPYRIGHT notice in top-level directory.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <zfp.h>
#include <zfp_toolkit.h>

size_t computeDataSize(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
        size_t dataLength;
        if(r1==0)
        {
                dataLength = 0;
        }
        else if(r2==0)
        {
                dataLength = r1;
        }
        else if(r3==0)
        {
                dataLength = r1*r2;
        }
        else if(r4==0)
        {
                dataLength = r1*r2*r3;
        }
        else if(r5==0)
        {
                dataLength = r1*r2*r3*r4;
        }
        else
        {
                dataLength = r1*r2*r3*r4*r5;
        }
        return dataLength;
}


/*ZFP compression/decompression*/
unsigned char* zfp_compression_float(float* data, int mode, double tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* outSize)
{
        int x = 1;
        char *y = (char*)&x;

        if(*y==1)
                sysEndianType = LITTLE_ENDIAN_SYSTEM;
        else //=0
                sysEndianType = BIG_ENDIAN_SYSTEM;

	double realTolerance = tolerance;
	if(mode==ZFP_REL)
	{
		size_t i = 0;
		size_t nbEle = computeDataSize(r5, r4, r3, r2, r1);
		double min = data[0], max = min;
		for(i=1;i<nbEle;i++)
		{
			if(min>data[i]) min = data[i];
			if(max<data[i]) max = data[i];
		}
		realTolerance = tolerance*(max-min);
	}
	zfp_exec_policy exec = zfp_exec_serial;
	unsigned char* buffer = NULL; //to store compressed data
	zfp_field* field = NULL;
	zfp_stream* zfp = NULL;
	bitstream* stream = NULL;
	zfp_type type = zfp_type_float;
//	size_t typesize = zfp_type_size(type);
	
	uint dims = 5;
	if(r5==0) dims--;
	if(r4==0) dims--;
	if(r3==0) dims--;
	if(r2==0) dims--;
	zfp = zfp_stream_open(NULL);
	field = zfp_field_alloc();
	
	zfp_field_set_pointer(field, data);
	
	zfp_field_set_type(field, type);	
	switch (dims) {
		case 1:
			zfp_field_set_size_1d(field, r1);
			break;
		case 2:
			zfp_field_set_size_2d(field, r1, r2);
			break;
		case 3:
			zfp_field_set_size_3d(field, r1, r2, r3);
			break;
	}	
	
	if(mode==ZFP_FXR)
		zfp_stream_set_rate(zfp, 32.0/realTolerance, zfp_type_float, dims, 0);
	else
		zfp_stream_set_accuracy(zfp, realTolerance);		
	
    /* allocate buffer for compressed data */
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    if (!bufsize) {
      fprintf(stderr, "invalid compression parameters\n");
      exit(0);
    }
    buffer = (unsigned char*)malloc(bufsize+sizeof(double));
    ZFP_doubleToBytes(buffer, realTolerance);
    
    if (!buffer) {
      fprintf(stderr, "cannot allocate memory\n");
      exit(0);
    }

	zfp_stream_set_execution(zfp, exec);

    /* associate compressed bit stream with memory buffer */
    stream = stream_open(buffer+sizeof(double), bufsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      exit(0);
    }
    zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);
	
	size_t zfpsize = zfp_compress(zfp, field);	
    if (zfpsize == 0) {
      fprintf(stderr, "compression failed\n");
      exit(0);
    }

	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
    
    *outSize = zfpsize;
    return buffer;
}

float* zfp_decompression_float(unsigned char* bytes, int mode, size_t outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
        int x = 1;
        char *y = (char*)&x;

        if(*y==1)
                sysEndianType = LITTLE_ENDIAN_SYSTEM;
        else //=0
                sysEndianType = BIG_ENDIAN_SYSTEM;

	zfp_exec_policy exec = zfp_exec_serial;
	//unsigned char* buffer = NULL; //to store compressed data
	zfp_field* field = NULL;
	zfp_stream* zfp = NULL;
	bitstream* stream = NULL;
	zfp_type type = zfp_type_float;
	size_t typesize = zfp_type_size(type);;
	uint dims = 5;
	if(r5==0) dims--;
	if(r4==0) dims--;
	if(r3==0) dims--;
	if(r2==0) dims--;
	size_t nbEle = 0;
	switch(dims)
	{
	case 1: nbEle = r1;break;
	case 2: nbEle = r1*r2;break;
	case 3: nbEle = r1*r2*r3;break;
	case 4: nbEle = r1*r2*r3*r4;break;	
	}

	field = zfp_field_alloc();
	zfp_field_set_type(field, type);	
	switch (dims) {
		case 1:
			zfp_field_set_size_1d(field, r1);
		break;
		case 2:
			zfp_field_set_size_2d(field, r1, r2);
		break;
		case 3:
			zfp_field_set_size_3d(field, r1, r2, r3);
		break;
	}	
	zfp = zfp_stream_open(NULL);
    
    double tolerance = ZFP_bytesToDouble(bytes);
    //printf("tolerance = %f\n", tolerance);
    
	if(mode==ZFP_FXR)
		zfp_stream_set_rate(zfp, 32.0/tolerance, zfp_type_float, dims, 0);
	else
		zfp_stream_set_accuracy(zfp, tolerance);	    
    
    /* allocate buffer for compressed data */
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    if (!bufsize) {
      fprintf(stderr, "invalid compression parameters\n");
      exit(0);
    }
    /* associate compressed bit stream with memory buffer */
    stream = stream_open(bytes+sizeof(double), bufsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      exit(0);
    }
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);
    
	//zfp_read_header(zfp, field, ZFP_HEADER_FULL);
    
	zfp_stream_set_execution(zfp, exec);		

    size_t rawsize = typesize * nbEle;
    float* fo = malloc(rawsize);
    if (!fo) {
      fprintf(stderr, "cannot allocate memory\n");
      exit(0);
    }
    zfp_field_set_pointer(field, fo);
    
    //compression
	zfp_decompress(zfp, field);
	
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);	
	
    return fo;
}

unsigned char* zfp_compression_double(double* data, int mode, double tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* outSize)
{
        int x = 1;
        char *y = (char*)&x;

        if(*y==1)
                sysEndianType = LITTLE_ENDIAN_SYSTEM;
        else //=0
                sysEndianType = BIG_ENDIAN_SYSTEM;

	double realTolerance = tolerance;
	if(mode==1)
	{
		size_t i = 0;
		size_t nbEle = computeDataSize(r5, r4, r3, r2, r1);
		double min = data[0], max = min;
		for(i=0;i<nbEle;i++)
		{
			if(min>data[i]) min = data[i];
			if(max<data[i]) max = data[i];
		}
		realTolerance = tolerance*(max-min);
	
	}
	zfp_exec_policy exec = zfp_exec_serial;
	unsigned char* buffer = NULL; //to store compressed data
	zfp_field* field = NULL;
	zfp_stream* zfp = NULL;
	bitstream* stream = NULL;
	zfp_type type = zfp_type_double;
	//size_t typesize = zfp_type_size(type);
	
	uint dims = 5;
	if(r5==0) dims--;
	if(r4==0) dims--;
	if(r3==0) dims--;
	if(r2==0) dims--;
	zfp = zfp_stream_open(NULL);
	field = zfp_field_alloc();
	
	zfp_field_set_pointer(field, data);
	
	zfp_field_set_type(field, type);	
	switch (dims) {
		case 1:
			zfp_field_set_size_1d(field, r1);
			break;
		case 2:
			zfp_field_set_size_2d(field, r1, r2);
			break;
		case 3:
			zfp_field_set_size_3d(field, r1, r2, r3);
			break;
	}	

        if(mode==ZFP_FXR)
                zfp_stream_set_rate(zfp, 64.0/realTolerance, zfp_type_double, dims, 0);
        else
                zfp_stream_set_accuracy(zfp, realTolerance);
	
    /* allocate buffer for compressed data */
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    if (!bufsize) {
      fprintf(stderr, "invalid compression parameters\n");
      exit(0);
    }
    buffer = (unsigned char*)malloc(bufsize+sizeof(double));
    ZFP_doubleToBytes(buffer, realTolerance);
    
    if (!buffer) {
      fprintf(stderr, "cannot allocate memory\n");
      exit(0);
    }

	zfp_stream_set_execution(zfp, exec);

    /* associate compressed bit stream with memory buffer */
    stream = stream_open(buffer+sizeof(double), bufsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      exit(0);
    }
    zfp_stream_set_bit_stream(zfp, stream);
	zfp_stream_rewind(zfp);
	
	size_t zfpsize = zfp_compress(zfp, field);	
    if (zfpsize == 0) {
      fprintf(stderr, "compression failed\n");
      exit(0);
    }

	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);
    
    *outSize = zfpsize;
    return buffer;
}

double* zfp_decompression_double(unsigned char* bytes, int mode, size_t outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
        int x = 1;
        char *y = (char*)&x;

        if(*y==1)
                sysEndianType = LITTLE_ENDIAN_SYSTEM;
        else //=0
                sysEndianType = BIG_ENDIAN_SYSTEM;

	zfp_exec_policy exec = zfp_exec_serial;
//	unsigned char* buffer = NULL; //to store compressed data
	zfp_field* field = NULL;
	zfp_stream* zfp = NULL;
	bitstream* stream = NULL;
	zfp_type type = zfp_type_double;
	size_t typesize = zfp_type_size(type);;
	uint dims = 5;
	if(r5==0) dims--;
	if(r4==0) dims--;
	if(r3==0) dims--;
	if(r2==0) dims--;
	size_t nbEle = 0;
	switch(dims)
	{
	case 1: nbEle = r1;break;
	case 2: nbEle = r1*r2;break;
	case 3: nbEle = r1*r2*r3;break;
	case 4: nbEle = r1*r2*r3*r4;break;	
	}

	field = zfp_field_alloc();
	zfp_field_set_type(field, type);	
	switch (dims) {
		case 1:
			zfp_field_set_size_1d(field, r1);
		break;
		case 2:
			zfp_field_set_size_2d(field, r1, r2);
		break;
		case 3:
			zfp_field_set_size_3d(field, r1, r2, r3);
		break;
	}	
	zfp = zfp_stream_open(NULL);
    
    double tolerance = ZFP_bytesToDouble(bytes);
    //printf("tolerance = %f\n", tolerance);
    
	if(mode==ZFP_FXR)
		zfp_stream_set_rate(zfp, 32.0/tolerance, zfp_type_double, dims, 0);
	else
		zfp_stream_set_accuracy(zfp, tolerance);	   
    
    /* allocate buffer for compressed data */
    size_t bufsize = zfp_stream_maximum_size(zfp, field);
    if (!bufsize) {
      fprintf(stderr, "invalid compression parameters\n");
      exit(0);
    }
    /* associate compressed bit stream with memory buffer */
    stream = stream_open(bytes+sizeof(double), bufsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      exit(0);
    }
    zfp_stream_set_bit_stream(zfp, stream);
    zfp_stream_rewind(zfp);

	//zfp_read_header(zfp, field, ZFP_HEADER_FULL);
    
	zfp_stream_set_execution(zfp, exec);		

    size_t rawsize = typesize * nbEle;
    double* fo = malloc(rawsize);
    if (!fo) {
      fprintf(stderr, "cannot allocate memory\n");
      exit(0);
    }
    zfp_field_set_pointer(field, fo);
    
    //compression
	zfp_decompress(zfp, field);
	
	zfp_field_free(field);
	zfp_stream_close(zfp);
	stream_close(stream);	
	
    return fo;
}

unsigned char* zfp_compression(void* data, int dataType, int mode, double tolerance, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1, size_t* outSize)
{
	if(dataType ==0)
	{
		return zfp_compression_float((float*)data, mode, tolerance, r5, r4, r3, r2, r1, outSize);
	}
	else if(dataType == 1)
	{
		return zfp_compression_double((double*)data, mode, tolerance, r5, r4, r3, r2, r1, outSize);
	}
	return NULL;
}

void* zfp_decompression(int dataType, int mode, unsigned char* bytes, size_t outSize, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
	if(dataType==0)
	{
		return zfp_decompression_float(bytes, mode, outSize, r5, r4, r3, r2, r1);	
	}
	else if(dataType==1)
	{
		return zfp_decompression_double(bytes, mode, outSize, r5, r4, r3, r2, r1);
	}
	return NULL;
}

