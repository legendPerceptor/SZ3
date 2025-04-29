/**
 *  @file zfp_compress.c
 *  @author Sheng Di
 *  @date April, 2022
 *  @brief This is an example of using compression interface
 *  (C) 202 by Mathematics and Computer Science (MCS), Argonne National Laboratory.
 *      See COPYRIGHT in top-level directory.
 */

#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include "zfp.h"
#include "zfpw.h"
#include "zfp_rw.h"

struct timeval startTime;
struct timeval endTime;  /* Start and end times */
struct timeval costStart; /*only used for recording the cost*/
double totalCost = 0;


void cost_start()
{
	totalCost = 0;
        gettimeofday(&costStart, NULL);
}

void cost_end()
{
        double elapsed;
        struct timeval costEnd;
        gettimeofday(&costEnd, NULL);
        elapsed = ((costEnd.tv_sec*1000000+costEnd.tv_usec)-(costStart.tv_sec*1000000+costStart.tv_usec))/1000000.0;
        totalCost += elapsed;
}


int main(int argc, char * argv[])
{
    size_t r5=0,r4=0,r3=0,r2=0,r1=0;
    char oriFilePath[640], compressedFilePath[645], outputFilePath[650];
    
    if(argc < 3)
    {
		printf("Test case: test_zfp [srcFilePath] [tolerance] [dimension sizes...]\n");
		printf("Example: test_zfp test.f32 0.01 8 8 128\n");
		exit(0);
    }
   
    sprintf(oriFilePath, "%s", argv[1]);
    
    double tolerance = atof(argv[2]);
    
    if(argc>=4)
		r1 = atoi(argv[3]); //8
    if(argc>=5)
		r2 = atoi(argv[4]); //8
    if(argc>=6)
		r3 = atoi(argv[5]); //128
    if(argc>=7)
        r4 = atoi(argv[6]);
    if(argc>=8)
        r5 = atoi(argv[7]);
   
    sprintf(compressedFilePath, "%s.zfp", oriFilePath);
    sprintf(outputFilePath, "%s.out", compressedFilePath);
  
    int status = 0; 
    size_t nbEle;
    float *data = ZFP_readFloatData(oriFilePath, &nbEle, &status);
    if(status != ZFP_SCES)
    {
		printf("Error: data file %s cannot be read!\n", oriFilePath);
		exit(0);
    }
   
    //compression

    size_t outSize = 0; 
    cost_start();
    unsigned char* bytes = zfp_compression(data, ZFP_FLOAT, ZFP_ABS, tolerance, r5, r4, r3, r2, r1, &outSize);
    cost_end();
    printf("compression time = %f\n",totalCost); 
    ZFP_writeByteData(bytes, outSize, compressedFilePath, &status);

    //decompression
    size_t compressedDataSize = outSize;
    cost_start();
    float* decData = zfp_decompression(ZFP_FLOAT, ZFP_ABS, bytes, compressedDataSize, r5, r4, r3, r2, r1);
    cost_end();
    printf("decompression time = %f\n",totalCost);
    ZFP_writeFloatData_inBytes(decData, nbEle, outputFilePath, &status);   

    printf("done\n");
    free(bytes); 
    free(data);
    free(decData);

    return 0;
}
