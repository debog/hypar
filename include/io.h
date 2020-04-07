/*! @file io.h
    @author Debojyoti Ghosh
    @brief Function declarations for file I/O functions.
*/

#ifndef _IO_H_
#define _IO_H_

int ReadArray(int,int,int*,int*,int,void*,void*,double*,double*,char*,int*);

int ReadArraywInterp( int,
                      int,
                      int*,
                      int*,
                      int*,
                      int,
                      void*,
                      void*,
                      double*,
                      double*,
                      char*,
                      int* );

int WriteArray    (int,int,int*,int*,int,double*,double*,void*,void*,char*);

int WriteBinary   (int,int,int*,double*,double*,char*,int*);
int WriteText     (int,int,int*,double*,double*,char*,int*);
int WriteTecplot2D(int,int,int*,double*,double*,char*,int*);
int WriteTecplot3D(int,int,int*,double*,double*,char*,int*);

#endif
