#ifndef SZ3_QCAT_DATAANALYSIS_H
#define SZ3_QCAT_DATAANALYSIS_H

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#define QCAT_FLOAT 0
#define QCAT_DOUBLE 1
#define QCAT_INT32 2
#define QCAT_INT16 3
#define QCAT_UINT32 4
#define QCAT_UINT16 5

struct QCAT_DataProperty {
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

double computeEntropy(int dataType, void* data, size_t nbEle) {
    size_t i = 0;
    unsigned char* bytes = (unsigned char*)data;
    size_t totalLen = dataType == QCAT_FLOAT ? nbEle * sizeof(float) : nbEle * sizeof(double);

    double entVal = 0.0;
    unsigned char index = 0;
    size_t table_size = 256;
    long* table = (long*)malloc(table_size * sizeof(long));
    memset(table, 0, table_size * sizeof(long));

    for (i = 0; i < totalLen; i++) {
        index = bytes[i];
        table[index]++;
    }

    size_t sum = nbEle * sizeof(float);
    for (i = 0; i < table_size; i++)
        if (table[i] != 0) {
            double prob = (double)table[i] / sum;
            entVal -= prob * log(prob) / log(2);
        }

    free(table);
    return entVal;
}

QCAT_DataProperty* computeProperty(int dataType, void* data, size_t nbEle) {
    QCAT_DataProperty* property = (QCAT_DataProperty*)malloc(sizeof(QCAT_DataProperty));
    memset(property, 0, sizeof(QCAT_DataProperty));

    property->dataType = dataType;
    property->numOfElem = nbEle;
    size_t i = 0;
    if (dataType == QCAT_FLOAT) {
        float* data_ = (float*)data;
        double min = data_[0], max = data_[0], sum = 0;
        for (i = 0; i < nbEle; i++) {
            if (min > data_[i])
                min = data_[i];
            if (max < data_[i])
                max = data_[i];
            sum += data_[i];
        }

        double med = min + (max - min) / 2;
        double sum_of_square = 0;
        for (i = 0; i < nbEle; i++)
            sum_of_square += (data_[i] - med) * (data_[i] - med);
        property->zeromean_variance = sum_of_square / nbEle;

        property->minValue = min;
        property->maxValue = max;

        property->avgValue = sum / nbEle;
        property->valueRange = max - min;
        property->totalByteSize = nbEle * sizeof(float);

    } else // QCAT_DOUBLE
    {
        double* data_ = (double*)data;
        double min = data_[0], max = data_[0], sum = 0;
        for (i = 0; i < nbEle; i++) {
            if (min > data_[i])
                min = data_[i];
            if (max < data_[i])
                max = data_[i];
            sum += data_[i];
        }

        double med = min + (max - min) / 2;
        double sum_of_square = 0;
        for (i = 0; i < nbEle; i++)
            sum_of_square += (data_[i] - med) * (data_[i] - med);
        property->zeromean_variance = sum_of_square / nbEle;

        property->minValue = min;
        property->maxValue = max;

        property->avgValue = sum / nbEle;
        property->valueRange = max - min;
        property->totalByteSize = nbEle * sizeof(double);
    }

    property->entropy = computeEntropy(dataType, data, nbEle);

    return property;
}

void printProperty(QCAT_DataProperty* property) {
    printf("numOfElem = %zu\n", property->numOfElem);
    printf("totalDataSize = %zu bytes (%f MB)\n", property->totalByteSize,
           property->totalByteSize / (1024.0f * 1024.0f));
    printf("min = %f\n", property->minValue);
    printf("max = %f\n", property->maxValue);
    printf("valueRange = %f\n", property->valueRange);
    printf("avgValue = %f\n", property->avgValue);
    printf("entropy = %f\n", property->entropy);
    printf("zeromean_variance = %f\n", property->zeromean_variance);
}

#endif // SZ3_QCAT_DATAANALYSIS_H
