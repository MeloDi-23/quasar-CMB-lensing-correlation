#ifndef MCMC_UTIL_H
#define MCMC_UTIL_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "mathutil.h"
#include "random.h"
#include "matrix.h"

/*** the \chi^2 ***/
double chi_square(double *F,double **B_covar,double *F_obs,int Nobs);

/*** Inverse Covariance Matrix ***/
void InvCovar(double **B_covar,int Nobs);

/*** random walk step (draw from a gaussian distribution) ***/
void RandomWalkStep(double *dap,double *sig,int Npar,long *idum);

/*** da'=Q*da => da=Q^{-1}da', da' has diagonal covariance matrix ***/
void ParTransfer(double *da,double **Qinv,double *dap,int Npar);

void UpdatePar(double *a,double *da,int Npar);
void VecterCopy(double *x_source,double *x_dest,int n);
void MatrixCopy(double **x_source,double **x_dest,int m,int n);
int  AcceptReject(double chi2_new,double chi2_old,long *idum);
void UpdateParCovar(double **ParCovar,double *a,int Npar);
void ParCovCorrection(double **ParCovar,int *iFreePar,int Npar);
void PrintChain(int N_accept,double AcceptanceRatio,double AcceptanceRatioCov,double chi2,double *a,int Npar,double *F,int Nobs,double *Q,int Ndeq,FILE *stream);

#endif
