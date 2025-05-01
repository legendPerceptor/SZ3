/**
 *  @file ByteToolkit.c
 *  @author Sheng Di
 *  @date April, 2016
 *  @brief Byte Toolkit
 *  (C) 2016 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */
 
#include <stdlib.h>
#include "zfp_toolkit.h"
	
int sysEndianType = 0;	
int dataEndianType = 0;

void ZFP_symTransform_8bytes(unsigned char data[8])
{
        unsigned char tmp = data[0];
        data[0] = data[7];
        data[7] = tmp;

        tmp = data[1];
        data[1] = data[6];
        data[6] = tmp;

        tmp = data[2];
        data[2] = data[5];
        data[5] = tmp;

        tmp = data[3];
        data[3] = data[4];
        data[4] = tmp;
}

void ZFP_symTransform_2bytes(unsigned char data[2])
{
        unsigned char tmp = data[0];
        data[0] = data[1];
        data[1] = tmp;
}

void ZFP_symTransform_4bytes(unsigned char data[4])
{
        unsigned char tmp = data[0];
        data[0] = data[3];
        data[3] = tmp;

        tmp = data[1];
        data[1] = data[2];
        data[2] = tmp;
}
	
	
	
int ZFP_bytesToInt_bigEndian(unsigned char* bytes)
{
	int temp = 0;
	int res = 0;
	
	res <<= 8;
	temp = bytes[0] & 0xff;
	res |= temp;	

	res <<= 8;
	temp = bytes[1] & 0xff;
	res |= temp;

	res <<= 8;
	temp = bytes[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = bytes[3] & 0xff;
	res |= temp;		
	
	return res;
}

/**
 * @unsigned char *b the variable to store the converted bytes (length=4)
 * @unsigned int num
 * */
void ZFP_intToBytes_bigEndian(unsigned char *b, unsigned int num)
{
	b[0] = (unsigned char)(num >> 24);	
	b[1] = (unsigned char)(num >> 16);	
	b[2] = (unsigned char)(num >> 8);	
	b[3] = (unsigned char)(num);	
	
	//note: num >> xxx already considered endian_type...
//if(dataEndianType==LITTLE_ENDIAN_DATA)
//		ZFP_symTransform_4bytes(*b); //change to BIG_ENDIAN_DATA
}

/**
 * @endianType: refers to the endian_type of unsigned char* b.
 * */
long ZFP_bytesToLong_bigEndian(unsigned char* b) {
	long temp = 0;
	long res = 0;
	//int i;

	res <<= 8;
	temp = b[0] & 0xff;
	res |= temp;

	res <<= 8;
	temp = b[1] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[2] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[3] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[4] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[5] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[6] & 0xff;
	res |= temp;
	
	res <<= 8;
	temp = b[7] & 0xff;
	res |= temp;						
	
	return res;
}

void ZFP_longToBytes_bigEndian(unsigned char *b, unsigned long num) 
{
	b[0] = (unsigned char)(num>>56);
	b[1] = (unsigned char)(num>>48);
	b[2] = (unsigned char)(num>>40);
	b[3] = (unsigned char)(num>>32);
	b[4] = (unsigned char)(num>>24);
	b[5] = (unsigned char)(num>>16);
	b[6] = (unsigned char)(num>>8);
	b[7] = (unsigned char)(num);
//	if(dataEndianType==LITTLE_ENDIAN_DATA)
//		ZFP_symTransform_8bytes(*b);
}


long ZFP_doubleToOSEndianLong(double value)
{
	ecldouble buf;
	buf.value = value;
	return buf.lvalue;
}

int ZFP_floatToOSEndianInt(float value)
{
	eclfloat buf;
	buf.value = value;
	return buf.ivalue;
}

/**
 * By default, the endian type is OS endian type.
 * */
short ZFP_bytesToShort(unsigned char* bytes)
{
	eclshort buf;
	memcpy(buf.byte, bytes, 2);
	
	return buf.svalue;
}

int ZFP_bytesToInt(unsigned char* bytes)
{
	eclfloat buf;
	memcpy(buf.byte, bytes, 4);
	return buf.ivalue;
}

long ZFP_bytesToLong(unsigned char* bytes)
{
	ecldouble buf;
	memcpy(buf.byte, bytes, 8);
	return buf.lvalue;
}

//the byte to input is in the big-endian format
float ZFP_bytesToFloat(unsigned char* bytes)
{
	eclfloat buf;
	memcpy(buf.byte, bytes, 4);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		ZFP_symTransform_4bytes(buf.byte);	
	return buf.value;
}

void ZFP_floatToBytes(unsigned char *b, float num)
{
	eclfloat buf;
	buf.value = num;
	memcpy(b, buf.byte, 4);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		ZFP_symTransform_4bytes(b);		
}

//the byte to input is in the big-endian format
double ZFP_bytesToDouble(unsigned char* bytes)
{
	ecldouble buf;
	memcpy(buf.byte, bytes, 8);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		ZFP_symTransform_8bytes(buf.byte);
	return buf.value;
}

void ZFP_doubleToBytes(unsigned char *b, double num)
{
	ecldouble buf;
	buf.value = num;
	memcpy(b, buf.byte, 8);
	if(sysEndianType==LITTLE_ENDIAN_SYSTEM)
		ZFP_symTransform_8bytes(b);
}

void ZFP_convertShortArrayToBytes(short* states, size_t stateLength, unsigned char* bytes)
{
	eclint16 ls;
	size_t i;
	if(sysEndianType==dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			ls.svalue = states[i];
			bytes[i*2] = ls.byte[0];
			bytes[i*2+1] = ls.byte[1];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			ls.svalue = states[i];
			bytes[i*2] = ls.byte[1];
			bytes[i*2+1] = ls.byte[0];
		}			
	}
}

void ZFP_convertUShortArrayToBytes(unsigned short* states, size_t stateLength, unsigned char* bytes)
{
	eclint16 ls;
	size_t i;
	if(sysEndianType==dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			ls.usvalue = states[i];
			bytes[i*2] = ls.byte[0];
			bytes[i*2+1] = ls.byte[1];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			ls.usvalue = states[i];
			bytes[i*2] = ls.byte[1];
			bytes[i*2+1] = ls.byte[0];
		}			
	}
}

void ZFP_convertIntArrayToBytes(int* states, size_t stateLength, unsigned char* bytes)
{
	eclint32 ls;
	size_t index = 0;
	size_t i;
	if(sysEndianType==dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.ivalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.ivalue = states[i];
			bytes[index] = ls.byte[3];
			bytes[index+1] = ls.byte[2];
			bytes[index+2] = ls.byte[1];
			bytes[index+3] = ls.byte[0];
		}			
	}
}

void ZFP_convertUIntArrayToBytes(unsigned int* states, size_t stateLength, unsigned char* bytes)
{
	eclint32 ls;
	size_t index = 0;
	size_t i;
	if(sysEndianType==dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.uivalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 2; //==i*4
			ls.uivalue = states[i];
			bytes[index] = ls.byte[3];
			bytes[index+1] = ls.byte[2];
			bytes[index+2] = ls.byte[1];
			bytes[index+3] = ls.byte[0];
		}			
	}
}

void ZFP_convertLongArrayToBytes(long* states, size_t stateLength, unsigned char* bytes)
{
	eclint64 ls;
	size_t index = 0;
	size_t i;
	if(sysEndianType==dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.lvalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
			bytes[index+4] = ls.byte[4];
			bytes[index+5] = ls.byte[5];
			bytes[index+6] = ls.byte[6];
			bytes[index+7] = ls.byte[7];	
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.lvalue = states[i];
			bytes[index] = ls.byte[7];
			bytes[index+1] = ls.byte[6];
			bytes[index+2] = ls.byte[5];
			bytes[index+3] = ls.byte[4];
			bytes[index+4] = ls.byte[3];
			bytes[index+5] = ls.byte[2];
			bytes[index+6] = ls.byte[1];
			bytes[index+7] = ls.byte[0];	
		}			
	}
}

void ZFP_convertULongArrayToBytes(unsigned long* states, size_t stateLength, unsigned char* bytes)
{
	eclint64 ls;
	size_t index = 0;
	size_t i;
	if(sysEndianType==dataEndianType)
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.ulvalue = states[i];
			bytes[index] = ls.byte[0];
			bytes[index+1] = ls.byte[1];
			bytes[index+2] = ls.byte[2];
			bytes[index+3] = ls.byte[3];
			bytes[index+4] = ls.byte[4];
			bytes[index+5] = ls.byte[5];
			bytes[index+6] = ls.byte[6];
			bytes[index+7] = ls.byte[7];			
		}		
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			index = i << 3; //==i*8
			ls.ulvalue = states[i];
			bytes[index] = ls.byte[7];
			bytes[index+1] = ls.byte[6];
			bytes[index+2] = ls.byte[5];
			bytes[index+3] = ls.byte[4];
			bytes[index+4] = ls.byte[3];
			bytes[index+5] = ls.byte[2];
			bytes[index+6] = ls.byte[1];
			bytes[index+7] = ls.byte[0];	
		}			
	}
}

unsigned short* ZFP_convertByteDataToUShortArray(unsigned char* bytes, size_t byteLength)
{
	eclint16 ls;
	size_t i, stateLength = byteLength/2;
	unsigned short* states = (unsigned short*)malloc(stateLength*sizeof(unsigned short));
	if(sysEndianType==dataEndianType)
	{	
		for(i=0;i<stateLength;i++)
		{
			ls.byte[0] = bytes[i*2];
			ls.byte[1] = bytes[i*2+1];
			states[i] = ls.usvalue;
		}
	}
	else
	{
		for(i=0;i<stateLength;i++)
		{
			ls.byte[0] = bytes[i*2+1];
			ls.byte[1] = bytes[i*2];
			states[i] = ls.usvalue;
		}		
	}
	return states;
} 
