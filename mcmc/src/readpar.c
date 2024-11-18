#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "typedef.h"
#include "nrutil.h"
#include "mathutil.h"
#include "file_util.h"

#define LINE_LEN 512
#define NAME_LEN 32
#define VALUE_LEN 64
#define MAX_PARNUM 64 

#define N_def 9
typedef struct
   {
    char Name[NAME_LEN];
    char Value1[VALUE_LEN];
    char Value2[VALUE_LEN];
    char Value3[VALUE_LEN];
    char Value4[VALUE_LEN];
    char Value5[VALUE_LEN];
   } PARSTRUCT;


static char Valid_Name[N_def][NAME_LEN]={
 "ng_avg",
 "err_ng",
 "lgMmin",
 "sig_lgM",
 "lgM0",
 "lgM1p",
 "alpha",
 "icen",
 "iFreePar"
};

int OK_check_par(PARSTRUCT par[MAX_PARNUM],int N_par)
{ 
  int i,j,count,error,wrong_value,valid;
  char line[LINE_LEN];
  char par0[NAME_LEN];

/* 
 ****** Judge the integerity and completeness of the parameters in .par file
 */

  error=0;

  for(i=0;i<N_par;i++)
     {
      valid=0;
      for(j=0;j<N_def;j++)
         if(strcmp(par[i].Name,Valid_Name[j])==0)
           { valid=1; break;}
      if(!valid) 
        { error=1;
          printf("*Warning: / %s / Invalid parameter.\n",par[i].Name);
        }
      } 

  for(i=0;i<N_def;i++)
     {
      count=0;
      wrong_value=0;
      j=0;
      while(j<N_par)
         {
           if(strcmp(par[j].Name,Valid_Name[i])==0) 
             { count++;
               if((par[j].Value1[0]==' ')||(par[j].Value1[0]=='#')
                ||(par[j].Value2[0]==' ')||(par[j].Value2[0]=='#') 
                ||(par[j].Value3[0]==' ')||(par[j].Value3[0]=='#') 
                ||(par[j].Value4[0]==' ')||(par[j].Value4[0]=='#') 
                ||(par[j].Value5[0]==' ')||(par[j].Value5[0]=='#')
                ) wrong_value=1;
            }
           j++;
         }
      if(wrong_value) 
        { 
          error=1;
          printf("***Error: / %s / - Wrong value given.\n",Valid_Name[i]);
        }
      if(count==0) 
        { 
          error=1;
          printf("***Error: / %s / - NO value given.\n",Valid_Name[i]);
        }
      if(count>1 && (i!=(N_def-1)))
        {
          error=1;
          printf("*Warning: / %s / - Duplicated value given.\n",Valid_Name[i]);
        }
     }

  if(error) return 0;
  else return 1;

}

void readpar(char *filename,int N_M,HOD *HODpar)
{ FILE *fp;
  char line[LINE_LEN];
  PARSTRUCT par0,par[MAX_PARNUM];
  double lgMmin;
 
  int i,j,k,N_par,Nspln,jmin;

  if(!(fp=fopen(filename,"r")))
    { printf("cannot open file %s.\n",filename);
      exit(0);
    }

  N_par=0;
  Nspln=0;

  fscanf(fp,"%s",par0.Name);
  if((par0.Name[0])!='#')
    { fscanf(fp,"%s",par0.Value1);
      if(strcmp(par0.Name,"lgM_lgNs")==0)
        fscanf(fp,"%s",par0.Value2);
      if(strcmp(par0.Name,"iFreePar")==0) {
        fscanf(fp,"%s",par0.Value2);
        fscanf(fp,"%s",par0.Value3);
        fscanf(fp,"%s",par0.Value4);
        fscanf(fp,"%s",par0.Value5);
      }
      fgets(line,LINE_LEN,fp);
    }
  else fgets(line,LINE_LEN,fp);

  while(!feof(fp))
    { 
      if((par0.Name[0])!='#')
        { 
          sprintf(par[N_par].Name,"%s",par0.Name);
          sprintf(par[N_par].Value1,"%s",par0.Value1);
          if(strcmp(par0.Name,"lgM_lgNs")==0)
            { sprintf(par[N_par].Value2,"%s",par0.Value2);
              Nspln=Nspln+1;
            }
          if(strcmp(par0.Name,"iFreePar")==0) {
            sprintf(par[N_par].Value2,"%s",par0.Value2);
            sprintf(par[N_par].Value3,"%s",par0.Value3);
            sprintf(par[N_par].Value4,"%s",par0.Value4);
            sprintf(par[N_par].Value5,"%s",par0.Value5);
          }
          N_par++;
         }
      fscanf(fp,"%s",par0.Name);
      if((par0.Name[0])!='#')
        { fscanf(fp,"%s",par0.Value1);
          if(strcmp(par0.Name,"lgM_lgNs")==0)
            fscanf(fp,"%s",par0.Value2);
          if(strcmp(par0.Name,"iFreePar")==0) {
            fscanf(fp,"%s",par0.Value2);
            fscanf(fp,"%s",par0.Value3);
            fscanf(fp,"%s",par0.Value4);
            fscanf(fp,"%s",par0.Value5);
          }
          fgets(line,LINE_LEN,fp);
        }
      else fgets(line,LINE_LEN,fp);
    }

  if(!OK_check_par(par,N_par)) exit(0);

/*** set values to COSMOS cosmology and HOD HODpar***/
  for(i=0;i<N_par;i++) {
       for(j=0;j<N_def;j++)
         if(strcmp(par[i].Name,Valid_Name[j])==0) break;
       switch (j) { 
           case  0 : (*HODpar).ng_avg     = atof(par[i].Value1); break;
           case  1 : (*HODpar).err_ng     = atof(par[i].Value1); break;
           case  2 : (*HODpar).lgMmin     = atof(par[i].Value1); break;
           case  3 : (*HODpar).sig_lgM    = atof(par[i].Value1); break;
           case  4 : (*HODpar).lgM0       = atof(par[i].Value1); break;
           case  5 : (*HODpar).lgM1p      = atof(par[i].Value1); break;
           case  6 : (*HODpar).alpha      = atof(par[i].Value1); break;
           case  7 : (*HODpar).icen       = atoi(par[i].Value1); break;
           case  8 : (*HODpar).iFreePar[1]= atoi(par[i].Value1);
                     (*HODpar).iFreePar[2]= atoi(par[i].Value2);
                     (*HODpar).iFreePar[3]= atoi(par[i].Value3);
                     (*HODpar).iFreePar[4]= atoi(par[i].Value4);
                     (*HODpar).iFreePar[5]= atoi(par[i].Value5);
                     break;
           default : printf("default");
       }
  }

  (*HODpar).N_M = N_M;
  (*HODpar).Nc  = (double *) calloc(N_M,sizeof(double));
  (*HODpar).Ns  = (double *) calloc(N_M,sizeof(double));
}

void free_HOD(HOD *HODpar)
{
 free((*HODpar).Nc);
 free((*HODpar).Ns);
}

/*
int main()
{
 COSMOS cosmology;
 HOD HODpar;
 int i;

 readpar("cosmos_hod.par",&cosmology,&HODpar);

 printf("Omega_m %4.2f \nOmega_L %4.2f\nOmega_b %4.2f\nh %4.2f\nitrans %d\nn_s %4.2f\nGamma %4.2f\nr_th %4.2f\ndelnorm %4.2f\n",cosmology.Omega_m,cosmology.Omega_L,cosmology.Omega_b,cosmology.h,cosmology.itrans,cosmology.n_s,cosmology.Gamma,cosmology.r_th,cosmology.delnorm);
 printf("ng_avg %4.2f\npomega %4.2f\nDeltaGamma %4.2f\nalpha_v %4.2f\nMmin %5.2e\nnu %4.2f\n",HODpar.ng_avg,HODpar.pomega,HODpar.DeltaGamma,HODpar.alpha_v,HODpar.Mmin,HODpar.nu);
 for(i=1;i<=HODpar.S_lgM.n;i++)
    printf("lgM_lgN %4.2f %f %e\n",HODpar.S_lgM.x[i],HODpar.S_lgM.y[i],HODpar.S_lgM.y2[i]);
 
 free_HOD(&HODpar);

}
*/
