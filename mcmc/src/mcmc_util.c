#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#include "mathutil.h"
#include "random.h"
#include "matrix.h"

/*** the \chi^2 ***/
double chi_square(double *F,double **B_covar,double *F_obs,int Nobs)
{
 double chi2;
 int    i,j;

 chi2=0.0;
 for(i=1;i<=Nobs;i++)
    for(j=1;j<=Nobs;j++)
       chi2=chi2+(F[i]-F_obs[i])*B_covar[i][j]*(F[j]-F_obs[j]);
 return chi2;
}

/*** Inverse Covariance Matrix ***/
void InvCovar(double **B_covar,int Nobs)
{
 gaussjinv(B_covar,Nobs);
}

/*** random walk step (draw from a gaussian distribution) ***/
void RandomWalkStep(double *dap,double *sig,int Npar,long *idum)
{
 int i;

 for(i=1;i<=Npar;i++) dap[i]=sig[i]*gasdev(idum);
}

/*** da'=Q*da => da=Q^{-1}da', da' has diagonal covariance matrix ***/
void ParTransfer(double *da,double **Qinv,double *dap,int Npar)
{
 int i,j;

 for(i=1;i<=Npar;i++) {
    da[i]=0.0;
    for(j=1;j<=Npar;j++)
       da[i]=da[i]+Qinv[i][j]*dap[j];
 }  
}

void UpdatePar(double *a,double *da,int Npar)
{
 int i;
 for(i=1;i<=Npar;i++) a[i]=a[i]+da[i];
}

void VecterCopy(double *x_source,double *x_dest,int n)
{
 int i;
 for(i=1;i<=n;i++) x_dest[i]=x_source[i];
}

void MatrixCopy(double **x_source,double **x_dest,int m,int n)
{
 int i,j;
 for(i=1;i<=m;i++) 
    for(j=1;j<=n;j++)
       x_dest[i][j]=x_source[i][j];
}

int AcceptReject(double chi2_new,double chi2_old,long *idum) 
{
 int    ACCEPT=0;
 double ratio,R;

 if(chi2_new<=chi2_old) ACCEPT=1;
 else { 
   ratio=exp(-0.5*(chi2_new-chi2_old));
   R=ran2(idum);
   if(R<ratio) ACCEPT=1;
 }

 return ACCEPT;
}

void UpdateParCovar(double **ParCovar,double *a,int Npar) 
{
 static int n=0;
 static double SS[32][32],S[32];
 int i,j;

 if(n==0) {
   if(Npar>31) {
     fprintf(stderr,"Npar=%d exceeds the dimesion %d set for covariance evaluation!\n",Npar,31);
     exit(-1);
   }
   for(i=1;i<=Npar;i++) {
      S[i]=0.0;
      for(j=1;j<=Npar;j++) 
         SS[i][j]=0.0;
   }
 }
 for(i=1;i<=Npar;i++) {
    S[i]=S[i]+a[i];
    for(j=i;j<=Npar;j++) {
       SS[i][j]=SS[i][j]+a[i]*a[j];
       SS[j][i]=SS[i][j];
    }
 }
 n++;

 for(i=1;i<=Npar;i++) {
    //for(j=1;j<i;j++) printf("%e ",0.);
    for(j=i;j<=Npar;j++) {
       ParCovar[i][j]=SS[i][j]/n-(S[i]/n)*(S[j]/n);
       ParCovar[j][i]=ParCovar[i][j];
       //printf("%e ",ParCovar[i][j]);
    }
    //printf("\n");
 }
}

/** 
 If the i-th parameter is fixed, the ij element of the parameter Covariance
 matrix should be zero, there would be a problem in obtaining eigenvalues 
 of the matrix (signular). The solution here is to use a small number (1e-20) 
 to replace zero for ii element and 0 for ij (j!=i).
 **/
void ParCovCorrection(double **ParCovar,int *iFreePar,int Npar)
{
 int i,j;

 for(i=1;i<=Npar;i++) {
    for(j=1;j<=Npar;j++) {
       if(iFreePar[i]==0 || iFreePar[j]==0) {
         if(i==j) ParCovar[i][j]=1e-20;
         else     ParCovar[i][j]=0.0;
       }
    }
 }
}

void PrintChain(int N_accept,double AcceptanceRatio,double AcceptanceRatioCov,double chi2,double *a,int Npar,double *F,int Nobs,double *Q,int Ndeq,FILE *stream)
{
 int i;
//    fprintf(stream,"%f %f %f ",AcceptanceRatio,AcceptanceRatioCov,chi2);
    fprintf(stream,"%6d %f %f ",N_accept,AcceptanceRatio,chi2);
    for(i=1;i<=Npar;i++) fprintf(stream,"%f ",a[i]);
    for(i=1;i< Nobs;i++) fprintf(stream,"%f ",F[i]);
                         fprintf(stream,"%e ",F[i]);
    for(i=1;i<=Ndeq;i++) fprintf(stream,"%f ",Q[i]);
    fprintf(stream,"\n");
}

