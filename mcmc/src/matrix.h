#ifndef MATRIX_H
#define MATRIX_H 

#include <math.h>

void jacobi(double **a,int n,double d[],double **v,int *nrot);

void eigsrt(double d[],double **v,int n);

void matrixproduct(double **A,double **B,double **AB,int m,int n);

#endif
