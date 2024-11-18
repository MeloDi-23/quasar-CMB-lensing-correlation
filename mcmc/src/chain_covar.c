#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mathutil.h"
#include "random.h"
#include "file_util.h"
#include "hod.h"
#include "getobs.h"
#include "readpar.h"
#include "matrix.h"
#include "mcmc_util.h"
#include "table_util.h"

void helpmessage(char *com)
{
  printf("MCMC code for 2-pt correlation function (by Zheng Zheng, gNFW Ver2.3, 20130424; RRSD+growth+7par_option)\n\n"
         "Usage: %s -seed chain_length model_in.par observ.dat covar.dat [outputfile]\n"
         "       seed          - random seed (positive integer)\n"
         "       chain_length  - length of the chain (positive integer)\n"
         "       model_in.par  - parameter file of the input cosmology + HOD\n"
         "                       (initial guess of the HOD)\n"
         "       observ.dat    - observational data points (rp, wp, err_wp)\n"
         "       outputfile    - outputfile recording the chain\n"
         "                       default name: mcmc.dat\n\n",
         com);
  exit(0);
}

int PriorSatisfied(double *a,int Npar)
{
 /* Priors         *
  * sigma_lgM>0    *
  * 0.5< alpha < 4 *
  * -1 < lgc_0 < 2 *
  *  0 < gamma < 3 */

 int OK;

 if(    a[2]<=0.0  
    || (a[5]<= 0.5 || a[5]>=4.0) 
   ) OK=0;
 else       OK=1;

 //printf("a[1]=%f a[2]=%f a[3]=%f a[4]=%f a[5]=%f OK=%d\n",a[1],a[2],a[3],a[4],a[5],OK);
 return OK;
}

int main(int argc,char **argv)
{
  HOD    HODpar;
  int    i,j;
  double *F_obs,*F,*F0,*a,*da,*a0,**b,**d,**dFda,**B_covar,*rp_data,*Err_obs;
  int    Nobs,Npar;
  double Q[6],Q0[6]; /** array of derived quantities **/
  int    Ndeq=5;

  double ftmp,chi2,chi2_0;
  double incfac0,incfac;
  int    cnt,EOM; /** EOM = End Of Minimization **/
  double kcost;  /** for cost function **/
  int    MONO;      /** force to be monotonicity ? **/
  double incstep,Ran;

  FILE   *fp;
  char   file_par[128],file_dat[128],file_cov[128],file_out[128],file_tab[128];

  double xdum,ydum,errydum;
  double ng_avg;

  double R,ratio;
  int    JUMP;
  long   idum;
  double lgMmin0,lgMmin,lgMc,lgMc0;
  double Mmin_0,M1_0;
  int    node;
  int    Nrandom,Naccept;
  int    ACCEPT;

  double *dap,**ParCovar,**ParCovarTmp,*sig;
  double **Eigenvector,*Eigenvalue;
  double chi2_old,chi2_new,AcceptanceRatio,AcceptanceRatioCov;
  int    burn_in,covar_sampling_length,covar_sampling_start,covar_sampling_end,total_samples;
  int    k,nrot,num_accept,num_accept_cov;
  int    chain_length;
  int    *iFreePar;

  table_struct table;
  int    *BinFlag;

  if(argc<6) helpmessage(argv[0]);
  idum=atoi(argv[1]);
  chain_length=atoi(argv[2]);
  strcpy(file_par,argv[3]);
  strcpy(file_dat,argv[4]);
  strcpy(file_cov,argv[5]);
  if(argc==7) strcpy(file_out,argv[6]);
  else        strcpy(file_out,"mcmc.dat");

  strcpy(file_tab,"wptable.bin");

  /*** Do Files exist? ***/
  CheckFile(file_par);
  CheckFile(file_dat);
  CheckFile(file_cov);
  CheckFile(file_tab);

  /** Load the component CF table **/
  load_table(file_tab,&table);

  /*** read HOD ***/
  readpar(file_par,table.N_M,&HODpar);

  Nobs=DataDimen(file_dat)+1;  /* 1 for number density */
  Npar=5;

  rp_data = dvector(1,Nobs);
  F_obs   = dvector(1,Nobs);
  Err_obs = dvector(1,Nobs);
  F       = dvector(1,Nobs);
  F0      = dvector(1,Nobs);
  B_covar = dmatrix(1,Nobs,1,Nobs);
  a       = dvector(1,Npar);
  da      = dvector(1,Npar);
  dap     = dvector(1,Npar);
  a0      = dvector(1,Npar);
  iFreePar    = ivector(1,Npar);
  sig         = dvector(1,Npar);
  ParCovar    = dmatrix(1,Npar,1,Npar);
  ParCovarTmp = dmatrix(1,Npar,1,Npar);
  Eigenvector = dmatrix(1,Npar,1,Npar);
  Eigenvalue  = dvector(1,Npar);

  /** observations & error-bars **/
  ReadData(rp_data,F_obs,Nobs-1,file_dat);
  F_obs[Nobs]  =HODpar.ng_avg; 
  Err_obs[Nobs]=HODpar.err_ng;
  ng_avg=HODpar.ng_avg;

  ReadCovar(B_covar,Nobs-1,file_cov);

  for(i=1;i<=Nobs-1;i++) { 
     for(j=1;j<=Nobs-1;j++) { 
//        B_covar[i][j]=Err_obs[i]*Err_obs[j]*B_covar[i][j];  
//        /** change to un-normalized one **/
        printf("%f ",B_covar[i][j]);
     }
     printf(" 0 \n");
  }
  for(i=1;i<=Nobs-1;i++) {
     B_covar[i][Nobs]=0.0;
     B_covar[Nobs][i]=0.0;
     printf("%f ",B_covar[Nobs][i]);
  }
  B_covar[Nobs][Nobs]=Err_obs[Nobs]*Err_obs[Nobs];
  printf("%e Err=%f\n",B_covar[Nobs][Nobs],Err_obs[Nobs]);
  B_covar[Nobs][Nobs]=1.0;
 
  Get_Covarinv_C(B_covar,Nobs);
//  Get_Covarinv_E(B_covar,Nobs,Err_obs);
  
  B_covar[Nobs][Nobs]=1./(Err_obs[Nobs]*Err_obs[Nobs]);
  printf("B done\n"); fflush(stdin);

  BinFlag = ivector(0,table.N_rp-1);
  CheckScale(rp_data,Nobs,table,BinFlag); /** data rp bin in the table? **/

  for(i=1;i<=Npar;i++) iFreePar[i]=HODpar.iFreePar[i];

  /** predictions **/
  Get_PredObs(F,Nobs,Q,Ndeq,HODpar,table,BinFlag);
  chi2_old=Get_chi2(F,B_covar,F_obs,Nobs);

  printf("chi^2=%f\n",chi2_old);
  for(i=1;i<=Nobs;i++)
     printf("F_obs[%2d]=%e F[%2d]=%e\n",i,F_obs[i],i,F[i]);
  printf("\nStart the MCMC ...\n");
  
  if(!(fp=fopen(file_out,"wt"))) FO_error(file_out);
  MONO=0; incfac0=0;
  fprintf(fp,"# %s %s %s %s %s %s %s\n",argv[0],argv[1],argv[2],argv[3],argv[4],argv[5],file_out);
  fprintf(fp,"# Npar = %2d ; Nobs = %2d ; ng_avg = %e +/- %e\n",Npar,Nobs,ng_avg,HODpar.err_ng);

  //printf("initial steps in lgM0, lgM1p, alpha_h, lgc_0, gamma");
  //printf("\n");
//  da[1]=0.05; da[2]=0.05; da[3]=0.10; da[4]=0.05; da[5]=0.05;                 //switching step sizes (this line is original)
//  da[1]=.5; da[2]=.5; da[3]=1; da[4]=.5; da[5]=.5;                 		//large step size to find chi^2 minimum
  da[1]=0.05; da[2]=0.05; da[3]=0.05; da[4]=0.05; da[5]=0.05;
  for(i=1;i<=Npar;i++) da[i]=da[i]*iFreePar[i];
//small step size
  fprintf(fp,"# ");
  for(i=1;i<=Npar;i++) {
//   scanf("%lf",&(da[i]));
     fprintf(fp,"%f ",da[i]);
     sig[i]=da[i];
  }
  fprintf(fp," [parameter stepsize; 0->fixed]\n");
//  fprintf(fp,"# %f %f %f %f %e %e %d [lgMmin_a sig_M_a lgMmin_b sig_M_b ng_col ng icen]\n",log10(HODpar.Mmin_a),HODpar.sig_M_a,log10(HODpar.Mmin_b),HODpar.sig_M_b,HODpar.ng_avg_col,ng_avg,HODpar.icen);

  HODpartoa(HODpar,a);
  for(i=1;i<=Npar;i++) a0[i]=a[i];
  for(i=1;i<=Nobs;i++) F0[i]=F[i];
  for(i=1;i<=Ndeq;i++) Q0[i]=Q[i];
  
  AcceptanceRatio=1.;  AcceptanceRatioCov=0.;
  PrintChain(1,AcceptanceRatio,AcceptanceRatioCov,chi2_old,a,Npar,F,Nobs,Q,Ndeq,fp);

  num_accept    =1;
  num_accept_cov=0;
  burn_in       =100;     /** minimum number of accepted points for burn-in **/
  covar_sampling_length=1000; /** length for covar evaluation after burn-in **/
  covar_sampling_start =burn_in+1; /** initial, update below **/
  covar_sampling_end   =burn_in+covar_sampling_length;
  total_samples        =chain_length; /** total length of the chain **/

/** start the MCMC process **/
  for(k=1;k<=total_samples;k++) {
    /** set the range to sample covar after burn-in **/
    if(num_accept==burn_in) {
      covar_sampling_start=k+1;
      covar_sampling_end  =k+covar_sampling_length;
    }

    /** save par and obs from the last step **/
    VecterCopy(a,a0,Npar);
    VecterCopy(F,F0,Nobs);
    VecterCopy(Q,Q0,Ndeq);

    /** during burn-in and covar sampling phases, use only simple covar **/
    /** after the covar is fully sampled, use the full covar            **/
    do  {
      VecterCopy(a0,a,Npar);
      if(k<=covar_sampling_end)
        RandomWalkStep(da,sig,Npar,&idum);
      else {
        RandomWalkStep(dap,sig,Npar,&idum);
        ParTransfer(da,Eigenvector,dap,Npar);
      }
      UpdatePar(a,da,Npar);
    } while (!PriorSatisfied(a,Npar));

    atoHODpar(a,&HODpar);
    Get_PredObs(F,Nobs,Q,Ndeq,HODpar,table,BinFlag);
    chi2_new=Get_chi2(F,B_covar,F_obs,Nobs);
    printf("%5d %8.3f ",k,chi2_new);
    for(i=1;i<=Npar;i++) printf("%f ",a[i]); //printf("\n"); fflush(stdin);

    /** back to the previous point if the trial fails **/
    if(AcceptReject(chi2_new,chi2_old,&idum)==0) {
      VecterCopy(a0,a,Npar);
      VecterCopy(F0,F,Nobs);
      VecterCopy(Q0,Q,Ndeq);
      chi2_new=chi2_old;
      printf(" Reject\n");
    }
    else {
      num_accept++;
      if(k>covar_sampling_end) num_accept_cov++;
      printf(" Accept\n");
    }
    fflush(stdin);
    AcceptanceRatio=num_accept/(k+1.);
    /** also acceptance ratio for the part of chains using the full covar **/
    if(k>covar_sampling_end)
      AcceptanceRatioCov=num_accept_cov/(k-covar_sampling_end+0.0);

    PrintChain(num_accept,AcceptanceRatio,AcceptanceRatioCov,chi2_new,a,Npar,F,Nobs,Q,Ndeq,fp);
    chi2_old=chi2_new;

    /** sample covar and use the covar after 'fully' sampled **/
    if(k>=covar_sampling_start && k<=covar_sampling_end)
      UpdateParCovar(ParCovar,a,Npar);
    if(k==covar_sampling_end) {
    /*** find Eigenvalues and Eigenvectors for the covar
         (and the rotation matrix for diagonalization) ***/
      ParCovCorrection(ParCovar,iFreePar,Npar);
      /** above numerical correction for covar related to fixed parameters **/
      MatrixCopy(ParCovar,ParCovarTmp,Npar,Npar);
      jacobi(ParCovarTmp,Npar,Eigenvalue,Eigenvector,&nrot);
      for(i=1;i<=Npar;i++) sig[i]=0.25*sqrt(Eigenvalue[i]);
    }
  }

  fclose(fp);

  free_dvector(rp_data,1,Nobs);
  free_dvector(F_obs,  1,Nobs);
  free_dvector(Err_obs,1,Nobs);
  free_dvector(F,      1,Nobs);
  free_dvector(F0,     1,Nobs);
  free_dvector(a,      1,Npar);
  free_dvector(da,     1,Npar);
  free_dvector(a0,     1,Npar);
  free_ivector(iFreePar,1,Npar);
  free_dvector(dap,    1,Npar);
  free_dvector(sig,    1,Npar);
  free_dmatrix(ParCovar,   1,Npar,1,Npar);
  free_dmatrix(ParCovarTmp,1,Npar,1,Npar);
  free_dmatrix(Eigenvector,1,Npar,1,Npar);
  free_dvector(Eigenvalue, 1,Npar);
  free_dmatrix(B_covar,    1,Nobs,1,Nobs);
  free_ivector(BinFlag,0,table.N_rp-1);
  free_table(&table);
  free_HOD(&HODpar);
}

