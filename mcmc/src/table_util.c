#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "typedef.h"
#include "file_util.h"

void load_table(char *filename,table_struct *table) 
{
 FILE   *fp;
 int    i,index,N_M,N_M2,N_rp,N1,N2;
 char   string[128];
 double h,Omega_m,Omega_L,Omega_b,n_s,sigma_8,redshift,boxsize,b;

 if(!(fp=fopen(filename,"rb"))) FO_error(filename);

 fread(string,         sizeof(char),128,fp);
 fread(&h,             sizeof(double),1,fp);
 fread(&Omega_m,       sizeof(double),1,fp);
 fread(&Omega_L,       sizeof(double),1,fp);
 fread(&Omega_b,       sizeof(double),1,fp);
 fread(&n_s,           sizeof(double),1,fp);
 fread(&sigma_8,       sizeof(double),1,fp);
 fread(&redshift,      sizeof(double),1,fp);
 fread(&boxsize,       sizeof(double),1,fp);
 fread(&b,             sizeof(double),1,fp);

 strcpy((*table).string,string);
 (*table).h       = h;
 (*table).Omega_m = Omega_m;
 (*table).Omega_L = Omega_L;
 (*table).Omega_b = Omega_b;
 (*table).n_s     = n_s;
 (*table).sigma_8 = sigma_8;
 (*table).redshift= redshift;
 (*table).boxsize = boxsize;
 (*table).b       = b;
 
 fread(&N_rp,sizeof(int),1,fp);
 fread(&N_M, sizeof(int),1,fp);
 fread(&N_M2,sizeof(int),1,fp);
 N1=N_rp*N_M;
 N2=N_rp*N_M2;

 (*table).N_rp=N_rp;
 (*table).N_M =N_M;
 (*table).N1  =N1;
 (*table).N2  =N2;

 (*table).rp      = (double *) calloc(N_rp,sizeof(double));
 (*table).rpmin   = (double *) calloc(N_rp,sizeof(double));
 (*table).rpmax   = (double *) calloc(N_rp,sizeof(double));
 (*table).lgM     = (double *) calloc(N_M, sizeof(double));
 (*table).lgMmin  = (double *) calloc(N_M, sizeof(double));
 (*table).lgMmax  = (double *) calloc(N_M, sizeof(double));
 (*table).nh      = (double *) calloc(N_M, sizeof(double));

 (*table).wf_cs   = (double *) calloc(N1,  sizeof(double));
 (*table).wf_ss   = (double *) calloc(N1,  sizeof(double));
 (*table).wp_cc   = (double *) calloc(N1,  sizeof(double));
 (*table).wp_cs   = (double *) calloc(N1,  sizeof(double));
 (*table).wp_ss   = (double *) calloc(N1,  sizeof(double));
 (*table).wpij_cc = (double *) calloc(N2,  sizeof(double));
 (*table).wpij_cs = (double *) calloc(N2,  sizeof(double));
 (*table).wpij_sc = (double *) calloc(N2,  sizeof(double));
 (*table).wpij_ss = (double *) calloc(N2,  sizeof(double));

 fread((*table).rp,     sizeof(double),N_rp,fp);
 fread((*table).rpmin,  sizeof(double),N_rp,fp);
 fread((*table).rpmax,  sizeof(double),N_rp,fp);

 fread((*table).lgM,    sizeof(double),N_M, fp);
 fread((*table).lgMmin, sizeof(double),N_M, fp);
 fread((*table).lgMmax, sizeof(double),N_M, fp);

 fread((*table).nh,     sizeof(double),N_M, fp);

 fread((*table).wf_cs,  sizeof(double),N1,  fp);
 fread((*table).wf_ss,  sizeof(double),N1,  fp);
 fread((*table).wp_cc,  sizeof(double),N1,  fp);
 fread((*table).wp_cs,  sizeof(double),N1,  fp);
 fread((*table).wp_ss,  sizeof(double),N1,  fp);

 fread((*table).wpij_cc,sizeof(double),N2,  fp);
 fread((*table).wpij_cs,sizeof(double),N2,  fp);
 fread((*table).wpij_sc,sizeof(double),N2,  fp);
 fread((*table).wpij_ss,sizeof(double),N2,  fp);
 
 fclose(fp);

}

void free_table(table_struct *table)
{
 free((*table).rp     );
 free((*table).rpmin  );
 free((*table).rpmax  );
 free((*table).lgM    );
 free((*table).lgMmin );
 free((*table).lgMmax );
 free((*table).nh     );

 free((*table).wf_cs  );
 free((*table).wf_ss  );
 free((*table).wp_cc  );
 free((*table).wp_cs  );
 free((*table).wp_ss  );
 free((*table).wpij_cc);
 free((*table).wpij_cs);
 free((*table).wpij_sc);
 free((*table).wpij_ss);
}

void CheckScale(double *rp_data,int Nobs,table_struct table,int *BinFlag)
{
 int i,j,found,stop;

 for(j=0;j<table.N_rp;j++) BinFlag[j]=0;
 
 stop=0;
 for(i=1;i<Nobs;i++) {
    found=0;
    for(j=0;j<table.N_rp;j++) {
       if(table.rpmin[j]<rp_data[i] && rp_data[i]<=table.rpmax[j]) {
         BinFlag[j]=1;
         found=1;
       }
    }
    if(!found) {
      printf("Sorry, rp scale %f is not in the table range!\n",rp_data[i]);
      stop=1;
    } 
 }
 if(stop) exit(-1); 

}
