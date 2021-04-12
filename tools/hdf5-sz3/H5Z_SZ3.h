//
// Created by apple on 2021/4/10.
//

#ifndef SZ_H5Z_SZ3_H
#define SZ_H5Z_SZ3_H

#include <hdf5.h>
//#include "sz.h"
#include <H5PLextern.h>
#include "compressor/SZ_General_Compressor.hpp"

#define H5Z_FILTER_SZ 32018
#define MAX_CHUNK_SIZE 4294967295 //2^32-1
static hid_t H5Z_SZ_ERRCLASS = -1;


extern int load_conffile_flag;
extern int init_sz_flag;

extern char cfgFile[256];

/* convenience macro to handle errors */
#define ERROR(FNAME)                                              \
do {                                                              \
    int _errno = errno;                                           \
    fprintf(stderr, #FNAME " failed at line %d, errno=%d (%s)\n", \
        __LINE__, _errno, _errno?strerror(_errno):"ok");          \
    return 1;                                                     \
} while(0)

#define H5Z_SZ_PUSH_AND_GOTO(MAJ, MIN, RET, MSG)     \
do                                                    \
{                                                     \
	H5Epush(H5E_DEFAULT,__FILE__,_funcname_,__LINE__, \
		H5Z_SZ_ERRCLASS,MAJ,MIN,MSG);                \
	retval = RET;                                     \
	goto done;                                        \
} while(0)


enum SZ_TYPE {
    SZ_DOUBLE,
    SZ_FLOAT,
    SZ_UINT8,
    SZ_UINT16,
    SZ_UINT32,
    SZ_UINT64,
    SZ_INT8,
    SZ_INT16,
    SZ_INT32,
    SZ_INT64
};

struct SZ3_config_params {
    bool fallback;
    bool usebitmap;
    bool preserve_sign;
    bool has_bg;
    float bg_value;
    float low, high;
};

void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1);
//void SZ_copymetaDataToCdArray(size_t* cd_nelmts, unsigned int *cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void SZ_metaDataToCdArray(size_t* cd_nelmts, unsigned int** cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

//void SZ_cdArrayToMetaDataErr(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1,
//                             int* error_bound_mode, float* abs_error, float* rel_error, float* pw_rel_error, float* psnr);
//
//void SZ_metaDataErrToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1,
//                             int error_bound_mode, float abs_error, float rel_error, float pw_rel_error, float psnr);
//
//int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[]);
static size_t H5Z_filter_sz(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf);
static herr_t H5Z_sz_set_local(hid_t dcpl_id, hid_t type_id, hid_t space_id);


void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

// SZ Old Functions
int computeDimension(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);
void H5Z_SZ3_Init(SZ::Compressor<float> *& sz, SZ::Compressor<float>*& sz_old, SZ3_config_params & sz3conf, int r3, int r2, int r1);
size_t computeDataLength(size_t r5, size_t r4, size_t r3, size_t r2, size_t r1);

void longToBytes_bigEndian(unsigned char *b, unsigned long num);
int bytesToInt_bigEndian(unsigned char* bytes);

void intToBytes_bigEndian(unsigned char *b, unsigned int num);
long bytesToLong_bigEndian(unsigned char* b);

#endif //SZ_H5Z_SZ3_H
