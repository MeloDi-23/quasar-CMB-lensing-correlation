#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ftio.h"

#define N_rp 19
#define N_pi 50

typedef struct {
double n_halo;
double f_cs[N_rp][N_pi];
double f_ss[N_rp][N_pi];
double xi2h_cc[N_rp][N_pi];
double xi2h_cs[N_rp][N_pi];
double xi2h_ss[N_rp][N_pi];
} table_struct;

void P_error(char *err_msg)
{
  printf("%s\n",err_msg);
  exit(0);
}

void FO_error(char *filename)
{
  printf("Cannot open file %s\n",filename);
  exit(0);
}

int main(int *argc,char **argv)
{
 FILE   *fp;
 char   filename[128];
 float  n_halo;
 int    nbyte;
 int    N,Nx,Ny;
 double *tmp_array,*wp;
 table_struct table;
 int    i,j,index;
 float  ftmp;
 double binsize_pi;
 //table_struct table1[1765];
 //exit(0);
 
 N = N_rp * N_pi;
 tmp_array = (double *) calloc(N,sizeof(double));
 wp        = (double *) calloc(N_rp,sizeof(double));

 //strcpy(filename,"table2_multidark_167167");
 strcpy(filename,argv[1]);
 if(!(fp=fopen(filename,"rb"))) {
   printf("cannot open file %s\n",filename);
   exit(0);
 }

 ftread(tmp_array,sizeof(double),N,fp);
 for(j=0;j<N_pi;j++)
    for(i=0;i<N_rp;i++) {
       index=i+j*N_rp;
       table.xi2h_cc[i][j]=tmp_array[index];
    }
 printf("array=%e\n",tmp_array[3500]);
 ftread(tmp_array,sizeof(double),N,fp);
 for(j=0;j<N_pi;j++)
    for(i=0;i<N_rp;i++) {
       index=i+j*N_rp;
       table.xi2h_cs[i][j]=tmp_array[index];
    }
 printf("array=%e\n",tmp_array[3500]);
 ftread(tmp_array,sizeof(double),N,fp);
 for(j=0;j<N_pi;j++)
    for(i=0;i<N_rp;i++) {
       index=i+j*N_rp;
       table.xi2h_ss[i][j]=tmp_array[index];
    }
 printf("array=%e\n",tmp_array[3500]);

 fclose(fp);
 free(tmp_array);


}
