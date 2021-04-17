/**
 *  @file H5Z_SZ.c
 *  @author Sheng Di
 *  @date July, 2017
 *  @brief SZ filter for HDF5
 *  (C) 2017 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <cstdio>
#include <cstdlib>
#include <cassert>
#include "H5Z_SZ3.h"

#include <compressor/SZ_General_Compressor.hpp>
#include <random>
#include <sstream>

static SZ::Compressor<float> *sz3, *sz_old;
static SZ3_config_params sz3conf;

const H5Z_class2_t H5Z_SZ[1] = {{
                                        H5Z_CLASS_T_VERS,              /* H5Z_class_t version */
                                        (H5Z_filter_t)H5Z_FILTER_SZ, /* Filter id number */
                                        1,              /* encoder_present flag (set to true) */
                                        1,              /* decoder_present flag (set to true) */
                                        "SZ3 compressor/decompressor for floating-point data.", /* Filter name for debugging */
                                        NULL,                          /* The "can apply" callback */
                                        H5Z_sz_set_local,                          /* The "set local" callback */
                                        (H5Z_func_t)H5Z_filter_sz,   /* The actual filter function */
                                }};

H5PL_type_t H5PLget_plugin_type(void) {return H5PL_TYPE_FILTER;}
const void *H5PLget_plugin_info(void) {return H5Z_SZ;}


static herr_t H5Z_sz_set_local(hid_t dcpl_id, hid_t type_id, hid_t chunk_space_id)
{
    printf("start in H5Z_sz_set_local, dcpl_id = %lld\n", dcpl_id);
    size_t r5=0,r4=0,r3=0,r2=0,r1=0, dsize;
    static char const *_funcname_ = "H5Z_sz_set_local";
    int i, ndims, ndims_used = 0;
    hsize_t dims[H5S_MAX_RANK], dims_used[5] = {0,0,0,0,0};
    herr_t retval = 0;
    H5T_class_t dclass;
    H5T_sign_t dsign;
    unsigned int flags = 0;
    //conf_params = H5Z_SZ_Init_Default();
//    if(load_conffile_flag)
//        H5Z_SZ_Init(cfgFile);
//    else
//        H5Z_SZ_Init(NULL);


    SZ_TYPE dataType = SZ_FLOAT;

    size_t mem_cd_nelmts = 0, cd_nelmts = 0;
    unsigned int mem_cd_values[12];
    int freshCdValues = 0;

    if (0 > (dclass = H5Tget_class(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a datatype");

    if (0 == (dsize = H5Tget_size(type_id)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "size is smaller than 0!");

    if (0 > (ndims = H5Sget_simple_extent_dims(chunk_space_id, dims, 0)))
        H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "not a data space");

    for (i = 0; i < ndims; i++)
    {
        if (dims[i] <= 1) continue;
        dims_used[ndims_used] = dims[i];
        ndims_used++;
    }

    //printf("dclass=%d, H5T_FLOAT=%d, H5T_INTEGER=%d\n", dclass, H5T_FLOAT, H5T_INTEGER);
    if (dclass == H5T_FLOAT)
        dataType = dsize==4? SZ_FLOAT: SZ_DOUBLE;
    else if(dclass == H5T_INTEGER)
    {
        if (0 > (dsign = H5Tget_sign(type_id)))
            H5Z_SZ_PUSH_AND_GOTO(H5E_ARGS, H5E_BADTYPE, -1, "Error in calling H5Tget_sign(type_id)....");
        if(dsign == H5T_SGN_NONE) //unsigned
        {
            switch(dsize)
            {
                case 1:
                    dataType = SZ_UINT8;
                    break;
                case 2:
                    dataType = SZ_UINT16;
                    break;
                case 4:
                    dataType = SZ_UINT32;
                    break;
                case 8:
                    dataType = SZ_UINT64;
                    break;
            }
        }
        else
        {
            switch(dsize)
            {
                case 1:
                    dataType = SZ_INT8;
                    break;
                case 2:
                    dataType = SZ_INT16;
                    break;
                case 4:
                    dataType = SZ_INT32;
                    break;
                case 8:
                    dataType = SZ_INT64;
                    break;
            }
        }
    }
    else
    {
        //printf("Error: dclass...\n");
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADTYPE, 0, "datatype class must be H5T_FLOAT or H5T_INTEGER");
    }

    if (0 > H5Pget_filter_by_id(dcpl_id, H5Z_FILTER_SZ, &flags, &mem_cd_nelmts, mem_cd_values, 0, NULL, NULL))
        H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_CANTGET, 0, "unable to get current SZ cd_values");


    switch(ndims_used)
    {
        case 1:
            r1 = dims_used[0];
            if(mem_cd_nelmts<=4)
                freshCdValues = 1;
            break;
        case 2:
            r1 = dims_used[0];
            r2 = dims_used[1];
            if(mem_cd_nelmts<=4)
                freshCdValues = 1;
            break;
        case 3:
            r1 = dims_used[0];
            r2 = dims_used[1];
            r3 = dims_used[2];
            if(mem_cd_nelmts<=5)
                freshCdValues = 1;
            break;
        case 4:
            r1 = dims_used[0];
            r2 = dims_used[1];
            r3 = dims_used[2];
            r4 = dims_used[3];
            if(mem_cd_nelmts<=6)
                freshCdValues = 1;
            break;
        default:
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "requires chunks w/1,2,3 or 4 non-unity dims");
    }

    if(freshCdValues)
    {
        unsigned int* cd_values = NULL;
        SZ_metaDataToCdArray(&cd_nelmts, &cd_values, dataType, r5, r4, r3, r2, r1);
        /* Now, update cd_values for the filter */
        if (0 > H5Pmodify_filter(dcpl_id, H5Z_FILTER_SZ, flags, cd_nelmts, cd_values))
            H5Z_SZ_PUSH_AND_GOTO(H5E_PLINE, H5E_BADVALUE, 0, "failed to modify cd_values");
        free(cd_values);
    }

    retval = 1;
    done:
    return retval;
}

/**
 * to be used in compression, and to be called outside H5Z_filter_sz().
 * */
void SZ_metaDataToCdArray(size_t* cd_nelmts, unsigned int **cd_values, int dataType, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    unsigned char bytes[8] = {0};
    unsigned long size;
    *cd_values = (unsigned int*)malloc(sizeof(unsigned int)*7);
    int dim = computeDimension(r5, r4, r3, r2, r1);
    (*cd_values)[0] = dim;
    (*cd_values)[1] = dataType;	//0: FLOAT ; 1: DOUBLE ; 2,3,4,....: INTEGER....
    switch(dim)
    {
        case 1:
            size = (unsigned long)r1;
            longToBytes_bigEndian(bytes, size);
            (*cd_values)[2] = bytesToInt_bigEndian(bytes);
            (*cd_values)[3] = bytesToInt_bigEndian(&bytes[4]);
            *cd_nelmts = 4;
            break;
        case 2:
            (*cd_values)[2] = (unsigned int) r2;
            (*cd_values)[3] = (unsigned int) r1;
            *cd_nelmts = 4;
            break;
        case 3:
            (*cd_values)[2] = (unsigned int) r3;
            (*cd_values)[3] = (unsigned int) r2;
            (*cd_values)[4] = (unsigned int) r1;
            *cd_nelmts = 5;
            break;
        case 4:
            (*cd_values)[2] = (unsigned int) r4;
            (*cd_values)[3] = (unsigned int) r3;
            (*cd_values)[4] = (unsigned int) r2;
            (*cd_values)[5] = (unsigned int) r1;
            *cd_nelmts = 6;
            break;
        default:
            (*cd_values)[2] = (unsigned int) r5;
            (*cd_values)[3] = (unsigned int) r4;
            (*cd_values)[4] = (unsigned int) r3;
            (*cd_values)[5] = (unsigned int) r2;
            (*cd_values)[6] = (unsigned int) r1;
            *cd_nelmts = 7;
    }
}


int checkCDValuesWithErrors(size_t cd_nelmts, const unsigned int cd_values[])
{
    int result = 0; //0 means no-error-information-in-cd_values; 1 means cd_values contains error information
    int dimSize = cd_values[0];
    //printf("nc_nelmts = %d\n", cd_nelmts);
    switch(dimSize)
    {
        case 1:
            if(cd_nelmts>4)
                result = 1;
            break;
        case 2:
            if(cd_nelmts>4)
                result = 1;
            break;
        case 3:
            if(cd_nelmts>5)
                result = 1;
            break;
        case 4:
            if(cd_nelmts>6)
                result = 1;
            break;
        case 5:
            if(cd_nelmts>7)
                result = 1;
            break;
    }
    return result;
}


void SZ_cdArrayToMetaData(size_t cd_nelmts, const unsigned int cd_values[], int* dimSize, int* dataType, size_t* r5, size_t* r4, size_t* r3, size_t* r2, size_t* r1)
{
    //printf("cd_nelmts=%zu\n", cd_nelmts);
    assert(cd_nelmts >= 4);
    unsigned char bytes[8];
    *dimSize = cd_values[0];
    *dataType = cd_values[1];

    switch(*dimSize)
    {
        case 1:
            intToBytes_bigEndian(bytes, cd_values[2]);
            intToBytes_bigEndian(&bytes[4], cd_values[3]);
            if(sizeof(size_t)==4)
                *r1 = (unsigned int)bytesToLong_bigEndian(bytes);
            else
                *r1 = (unsigned long)bytesToLong_bigEndian(bytes);
            *r2 = *r3 = *r4 = *r5 = 0;
            break;
        case 2:
            *r3 = *r4 = *r5 = 0;
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        case 3:
            *r4 = *r5 = 0;
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        case 4:
            *r5 = 0;
            *r4 = cd_values[5];
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
            break;
        default:
            *r5 = cd_values[6];
            *r4 = cd_values[5];
            *r3 = cd_values[4];
            *r2 = cd_values[3];
            *r1 = cd_values[2];
    }
    //printf("SZ_cdArrayToMetaData: r5=%zu, r4=%zu, r3=%zu, r2=%zu, r1=%zu\n", *r5, *r4, *r3, *r2, *r1);
}

static size_t H5Z_filter_sz(unsigned int flags, size_t cd_nelmts, const unsigned int cd_values[], size_t nbytes, size_t* buf_size, void** buf)
{
    printf("start in H5Z_filter_sz, cd_nelmts=%zu\n", cd_nelmts);
    //H5Z_SZ_Init_Default();


    size_t r1 = 0, r2 = 0, r3 = 0, r4 = 0, r5 = 0;
    int dimSize = 0, dataType = 0;

    if(cd_nelmts==0) //this is special data such as string, which should not be treated as values.
        return nbytes;

    SZ_cdArrayToMetaData(cd_nelmts, cd_values, &dimSize, &dataType, &r5, &r4, &r3, &r2, &r1);

    size_t nbEle = computeDataLength(r5, r4, r3, r2, r1);

    if(nbEle < 200)
        return nbytes;
    printf("In H5Z_SZ3.cpp, dimensions: r5=%zu, r4=%zu, r3=%zu, r2=%zu, r1=%zu\n", r5,r4,r3,r2,r1);
    H5Z_SZ3_Init(sz3, sz_old, sz3conf, r3, r2, r1);

    if (flags & H5Z_FLAG_REVERSE)
    {
        //cost_start();
        /* decompress data */
        if(dataType == SZ_FLOAT)//==0
        {
            float* data;
            if(sz3conf.fallback) {
                data= sz_old->decompress((unsigned char*)(*buf), nbytes);
            }else if(!sz3conf.has_bg){
                data= sz3->decompress((unsigned char*)(*buf), nbytes);
            }else {
                data = sz3->decompress_withBG((unsigned char*)(*buf), nbytes, sz3conf.bg_value,sz3conf.low, sz3conf.high, sz3conf.usebitmap, sz3conf.preserve_sign, sz3conf.has_bg);
            }
            free(*buf);
            *buf = data;
            *buf_size = nbEle*sizeof(float);
        }
        else
        {
            printf("Decompression error: unknown data type: %d\n", dataType);
            exit(0);
        }
        //cost_end();
        //printf("decompression time = %lf, decompression rate = %lf\n", totalCost, 1.0*nbEle*sizeof(float)/totalCost);
    }
    else //compression
    {
        size_t outSize = 0;

        //printf("r5=%d, r4=%d, r3=%d, r2=%d, r1=%d, dataType=%d\n", r5, r4, r3, r2, r1, dataType);
        //cost_start();
        if(dataType == SZ_FLOAT)//==0
        {
            float* data = (float*)(*buf);
            unsigned char *bytes = nullptr;
            if(sz3conf.fallback) {
//                std::random_device rd;  //Will be used to obtain a seed for the random number engine
//                std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
//                std::uniform_real_distribution<> dis(0, 100000000);
//                std::stringstream ss;
//                ss << "QW-"<<dis(gen) << ".szt3";
//                SZ::writefile(ss.str().data(), data, r1);
                bytes = sz_old->compress(data, outSize);
            } else if(!sz3conf.has_bg) {
                bytes = sz3->compress(data, outSize);
            } else {
                bytes = sz3->compress_withBG(data, outSize, sz3conf.bg_value, sz3conf.low, sz3conf.high, sz3conf.usebitmap, sz3conf.preserve_sign, sz3conf.has_bg);
            }
            free(*buf);
            *buf = bytes;
            *buf_size = outSize;
        } else
        {
            printf("Compression error: unknown data type: %d\n", dataType);
            exit(0);
        }
        //cost_end();
        //printf("compression time = %lf, compression rate = %lf\n", totalCost, 1.0*nbEle*sizeof(float)/totalCost);
    }

    //H5Z_SZ_Finalize();
    return *buf_size;
}

void init_dims_chunk(int dim, hsize_t dims[5], hsize_t chunk[5], size_t nbEle, size_t r5, size_t r4, size_t r3, size_t r2, size_t r1)
{
    switch(dim)
    {
        case 1:
            dims[0] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
                chunk[0] = r1;
            else
                chunk[0] = 2147483648;//2^31
            break;
        case 2:
            dims[0] = r2;
            dims[1] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r2;
                chunk[1] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        case 3:
            dims[0] = r3;
            dims[1] = r2;
            dims[2] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r3;
                chunk[1] = r2;
                chunk[2] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        case 4:
            dims[0] = r4;
            dims[1] = r3;
            dims[2] = r2;
            dims[3] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r4;
                chunk[1] = r3;
                chunk[2] = r2;
                chunk[3] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
            break;
        default:
            dims[0] = r5;
            dims[1] = r4;
            dims[2] = r3;
            dims[3] = r2;
            dims[4] = r1;
            if(nbEle <= MAX_CHUNK_SIZE) //2^32-1
            {
                chunk[0] = r5;
                chunk[1] = r4;
                chunk[2] = r3;
                chunk[3] = r2;
                chunk[4] = r1;
            }
            else
            {
                printf("Error: size is too big!\n");
                exit(0);
            }
    }
}


