#ifndef MATHUTIL_H
#define MATHUTIL_H

#include <math.h>

double gammln(double xx);

double factln(int n);

double bico(int n,int k);

double lnbico(int n,int k);

double beta(double z,double w);

void spline(double x[],double y[],int n,double yp1,double ypn,double y2[]);
void splint(double xa[],double ya[],double y2a[],int n,double x,double *y);

void part_spline(double x[],double y[],int n,int n1,int n2,double ypn1,double ypn2,double y2[]);
void part_splint(double xa[],double ya[],double y2a[],int n,int n1,int n2,double x,double *y,int iexp);

void gaussj(double **a,int n,double **b,int m);
void gaussjinv(double **a,int n);

double zbrent(double (*func)(),double x1,double x2,double tol);

void bilinint(double *x,double *y,double **z,int NDIMx,int NDIMy,double xt,double yt,double *zt);

void bilinint3(double *x,double *y,double **zb,double **zc,double **zd,int NDIMx,int NDIMy,double xt,double yt,double *zbt,double *zct,double *zdt);

#endif

