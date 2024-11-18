#ifndef GET_OBS_H
#define GET_OBS_H

#include <math.h>
#include "typedef.h"

void Get_Covarinv(double **B_covar,int Nobs,double *F_obs);

void Get_Covarinv_E(double **B_covar,int Nobs,double *Err_obs);
void Get_Covarinv_C(double **B_covar,int Nobs);

double Get_chi2(double *F,double **B_covar,double *F_obs,int Nobs);

void atoHODpar(double *a,HOD *HODpar);

void HODpartoa(HOD HODpar,double *a);

void Get_PredObs(double *F, int Nobs, double *Q, int Ndeq, HOD HODpar, table_struct table,int *BinFlag);

#endif

