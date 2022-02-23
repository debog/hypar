/*! @file arrayfunctions_gpu.h
    @brief Contains function definitions for common array operations on GPU.
    @author Youngdae Kim
 */

#ifndef _ARRAYFUNCTIONS_GPU_H_
#define _ARRAYFUNCTIONS_GPU_H_

#include <arrayfunctions.h>

#ifdef __cplusplus
extern "C" {
#endif

enum gpuMemcpyKind {
    gpuMemcpyHostToDevice = 0,
    gpuMemcpyDeviceToHost,
    gpuMemcpyDeviceToDevice
};

void gpuSetDevice(int device);
void gpuMemcpy(void*,const void*,size_t,enum gpuMemcpyKind);
void gpuMalloc(void**,size_t);
void gpuMemset(void*,int,size_t);
void gpuFree(void*);

void gpuArrayCopy1D(const double*,double*,int);
void gpuArraySetValue(double*,int,double);
void gpuArrayAXPY(const double*,double,double*,int);
void gpuArrayBlockMultiply(double*,const double*,int,int);
double gpuArraySumSquarenD(int,int,int*,int,int*,double*);

void gpuArrayCopy1DNewScheme(const double*,double*,int,int);
void gpuArrayCheckEqual(const char*,const double*,const double*,int,int);
void gpuArrayCheckEqual2(const char*,const double*,const double*,int);

#ifdef __cplusplus
}
#endif

#endif
