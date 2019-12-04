/*! @file io_cpp.h
    @author Debojyoti Ghosh
    @brief C++ Function declarations for file I/O functions.
*/

#ifndef _IO_CPP_H_
#define _IO_CPP_H_

extern "C" int WriteBinary   (int,int,int*,double*,double*,char*,int*);
extern "C" int WriteText     (int,int,int*,double*,double*,char*,int*);
extern "C" int WriteTecplot2D(int,int,int*,double*,double*,char*,int*);
extern "C" int WriteTecplot3D(int,int,int*,double*,double*,char*,int*);

#endif
