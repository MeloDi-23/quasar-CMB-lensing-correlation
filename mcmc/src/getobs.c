#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "typedef.h"
#include "mathutil.h"
#include "hod.h"
//#include "MedianMass.h"

void Get_Covarinv(double **B_covar,int Nobs,double *F_obs)
{ int i,j;

  /** assume uncorrelated, 10% error **/
  for(i=1;i<=Nobs;i++)
    for(j=1;j<=Nobs;j++)
       { if (i==j) B_covar[i][i]=(0.10*F_obs[i])*(0.10*F_obs[i]);
         else B_covar[i][j]=0.0;
       }
  gaussjinv(B_covar,Nobs);

}

void Get_Covarinv_E(double **B_covar,int Nobs,double *Err_obs)
{ int i,j;

  /** assume uncorrelated, using observational error-bars **/
  for(i=1;i<=Nobs;i++)
    for(j=1;j<=Nobs;j++)
       { if (i==j) B_covar[i][i]=Err_obs[i]*Err_obs[i];
         else B_covar[i][j]=0.0;
       }
  gaussjinv(B_covar,Nobs);

}

/*** if covariance matrix is given, obtain the inverse ***/
void Get_Covarinv_C(double **B_covar,int Nobs)
{
 gaussjinv(B_covar,Nobs);
}

/*** the \chi^2 ***/
double Get_chi2(double *F,double **B_covar,double *F_obs,int Nobs)
{ double chi2;
  int i,j;

  chi2=0.0;
  for(i=1;i<=Nobs;i++)
     for(j=1;j<=Nobs;j++)
        chi2=chi2+(F[i]-F_obs[i])*B_covar[i][j]*(F[j]-F_obs[j]);
//  printf("chi2=%f\n",chi2);
  return chi2;
}

/*** Compute the predicted observables and derived quantities ***/
void Get_PredObs(double *F, int Nobs, double *Q, int Ndeq, HOD HODpar, table_struct table, int *BinFlag)
{
 double ng,ng_sat,b_g,beta_g;
 double Mmed_cen,Mmed_sat; 
 double ni,nj,ninj,ngsq,nisq;
 double Ns,NcNc,NcNs,NsNc,NsNs,lgM,lgM1;
 double wf_cs,wf_ss,wp_cc,wp_cs,wp_sc,wp_ss;
 int    ir,i,j,index,indexM,N_M2,iobs;
 double npair_cs,npair_ss;

 ng    =0.0;
 ng_sat=0.0;
 npair_cs=0.0;
 npair_ss=0.0;
 for(i=0;i<table.N_M;i++) {
    ni          = table.nh[i];
    lgM         = table.lgM[i];
    HODpar.Nc[i]= Nc_avg(lgM,HODpar);
    HODpar.Ns[i]= Ns_avg(lgM,HODpar); 
    HODpar.Ns[i]= HODpar.Ns[i]*HODpar.Nc[i];  /** x Ncen profile  **/
    ng          = ng     + ni*(HODpar.Nc[i]+HODpar.Ns[i]);
    ng_sat      = ng_sat + ni*HODpar.Ns[i];
    npair_cs    = npair_cs + ni*(HODpar.Nc[i]*HODpar.Ns[i]);
    npair_ss    = npair_ss + ni*(0.5*HODpar.Ns[i]*HODpar.Ns[i]);
 }
 ngsq=ng*ng;
 //printf("%e %e npair_cs npair_ss\n",npair_cs,npair_ss);

 iobs=1;
 N_M2=table.N_M*(table.N_M-1)/2;
 for(ir=0;ir<table.N_rp;ir++) {
    if(BinFlag[ir]==0) continue;
    wf_cs=0.0;
    wf_ss=0.0;
    wp_cc=0.0;
    wp_cs=0.0;
    wp_sc=0.0;
    wp_ss=0.0;

    for(i=0;i<table.N_M;i++) {
       index=ir*table.N_M+i;
       ni=table.nh[i];
       nisq=ni*ni;
       NcNc=HODpar.Nc[i]*HODpar.Nc[i];
       NcNs=HODpar.Nc[i]*HODpar.Ns[i];
       //NcNs=             HODpar.Ns[i];
       NsNs=HODpar.Ns[i]*HODpar.Ns[i];
       //Ns  =HODpar.Ns[i];
       wf_cs=wf_cs+2.0*ni  /ngsq*NcNs*table.wf_cs[index];
       wf_ss=wf_ss+    ni  /ngsq*NsNs*table.wf_ss[index];
       wp_cc=wp_cc+    nisq/ngsq*NcNc*table.wp_cc[index];
       wp_cs=wp_cs+2.0*nisq/ngsq*NcNs*table.wp_cs[index];
       wp_ss=wp_ss+    nisq/ngsq*NsNs*table.wp_ss[index];
       //printf("ii %3d %f %f %f %f %f\n",i,wf_cs,wf_ss,wp_cc,wp_cs,wp_ss);
    } 

    indexM=0;
    for(i=0;i<table.N_M;i++) {
       ni=table.nh[i];
       for(j=i+1;j<table.N_M;j++) {
          index=ir*N_M2+indexM;
          //printf("%03d %03d %6d\n",i,j,index);
          nj=table.nh[j];
          ninj=ni*nj;
          NcNc=HODpar.Nc[i]*HODpar.Nc[j];
          NcNs=HODpar.Nc[i]*HODpar.Ns[j];
          NsNc=HODpar.Ns[i]*HODpar.Nc[j];
//          NcNs=0.5*(HODpar.Nc[i]*HODpar.Ns[j]+HODpar.Nc[j]*HODpar.Ns[i]);
          NsNs=HODpar.Ns[i]*HODpar.Ns[j];
          wp_cc=wp_cc+    ninj/ngsq*NcNc*table.wpij_cc[index] *2.0; 
          wp_cs=wp_cs+    ninj/ngsq*NcNs*table.wpij_cs[index] *2.0;
          wp_sc=wp_sc+    ninj/ngsq*NsNc*table.wpij_sc[index] *2.0;
          wp_ss=wp_ss+    ninj/ngsq*NsNs*table.wpij_ss[index] *2.0; 
          /** the 2.0 factor at the end is for i-j symmetry (for i!=j) **/
          indexM++;
       }
       //printf("ij %3d %3d %f %f %f\n",i,j,wp_cc,wp_cs,wp_ss);
    }
    
//    F[ir+1]=wf_cs+wf_ss+wp_cc+wp_cs+wp_ss;
    F[iobs]=wf_cs+wf_ss+wp_cc+wp_cs+wp_sc+wp_ss;
    //printf("%7.4f %7.4e %7.4e %7.4e %7.4e %7.4e %7.4e\n",table.rp[ir],wf_cs,wf_ss,wp_cc,(wp_cs+wp_sc),wp_ss,F[iobs]);
    iobs++;
 }
 F[Nobs]=ng;

// Mmed_cen = MedianMass(cosmology,HODpar,1);
// Mmed_sat = MedianMass(cosmology,HODpar,0);

 Q[1]=b_g;
 Q[2]=ng_sat/ng;
 Q[3]=Mmed_cen;
 Q[4]=Mmed_sat;
 Q[5]=lgM1;

 return;
}

/*** copy the value of a to HODpar **/
/***
     only copy a[], other values (e.g. ng_avg) should be 
     copied from somewhere else. 
     At present, 
     a[1]=lgMmin
     a[2]=sig_lgM
     a[3]=lgM0
     a[4]=lgM1p
     a[5]=alpha
 
     will change when more parameters are involved
 ***/
void atoHODpar(double *a,HOD *HODpar)
{ 
 int i;

 (*HODpar).lgMmin  =a[1];
 (*HODpar).sig_lgM =a[2];
 (*HODpar).lgM0    =a[3];
 (*HODpar).lgM1p   =a[4];
 (*HODpar).alpha   =a[5];
}

/*** copy the value of HODpar to a **/
void HODpartoa(HOD HODpar,double *a)
{ 
 int i;

 a[1]=HODpar.lgMmin;
 a[2]=HODpar.sig_lgM;
 a[3]=HODpar.lgM0;
 a[4]=HODpar.lgM1p;
 a[5]=HODpar.alpha;
}

