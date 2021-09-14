//
// Created by apple on 2021/9/14.
//

#ifndef SZ_QCAT_DATAANALYSIS_H
#define SZ_QCAT_DATAANALYSIS_H


#include <stdio.h>
#include <stdint.h>

#define QCAT_FLOAT 0
#define QCAT_DOUBLE 1
#define QCAT_INT32 2
#define QCAT_INT16 3
#define QCAT_UINT32 4
#define QCAT_UINT16 5


struct QCAT_DataProperty
{
    int dataType; /*DA_DOUBLE or DA_FLOAT*/
    size_t r5;
    size_t r4;
    size_t r3;
    size_t r2;
    size_t r1;

    long numOfElem;
    double minValue;
    double maxValue;
    double valueRange;
    double avgValue;
    double entropy;
    double zeromean_variance;
    size_t totalByteSize;
};

double computeEntropy(int dataType, void* data, size_t nbEle);
QCAT_DataProperty* computeProperty(int dataType, void* data, size_t nbEle);
void printProperty(QCAT_DataProperty* property);

#endif //SZ_QCAT_DATAANALYSIS_H
