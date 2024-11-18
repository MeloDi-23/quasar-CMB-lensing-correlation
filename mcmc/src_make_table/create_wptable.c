#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "ftio.h"

#define N_rp 19
#define N_pi 50
#define N_s  19
#define N_mu 20
#define N_M  199
#define N_M2 19701

double *tmp_array;
int    N;
int    k;

typedef struct {
double n_halo;
/* xi(rp,pi) */
double f_cs[N_rp][N_pi];
double f_ss[N_rp][N_pi];
double xi2h_cc[N_rp][N_pi];
double xi2h_cs[N_rp][N_pi];
double xi2h_sc[N_rp][N_pi];
double xi2h_ss[N_rp][N_pi];
/* xi(s,mu) */
double Sf_cs[N_s][N_mu];
double Sf_ss[N_s][N_mu];
double Sxi2h_cc[N_s][N_mu];
double Sxi2h_cs[N_s][N_mu];
double Sxi2h_sc[N_s][N_mu];
double Sxi2h_ss[N_s][N_mu];
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

double *vector(nl,nh)
int nl,nh;
{
        double *v;

        v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
        if (!v) P_error("allocation failure in vector()");
        return v-nl;
}

double **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
{
        int i;
        double **m;

        m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
        if (!m) P_error("allocation failure 1 in matrix()");
        m -= nrl;

        for(i=nrl;i<=nrh;i++) {
                m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
                if (!m[i]) P_error("allocation failure 2 in matrix()");
                m[i] -= ncl;
        }
        return m;
}

void free_vector(v,nl,nh)
double *v;
int nl,nh;
{
        free((char*) (v+nl));
}

void free_matrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
{
        int i;

        for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
        free((char*) (m+nrl));
}

void read2proj(FILE *fp,double **Arrayw,int iM)
{
 extern double *tmp_array;
 extern int    N,k;
 int           i,j,index;
 double        *w,binsize_pi=2.0;

 w=vector(0,N_rp-1);

 ftread(tmp_array,sizeof(double),N,fp);
 for(i=0;i<N_rp;i++) w[i]=0.0;
 for(j=0;j<N_pi;j++) 
    for(i=0;i<N_rp;i++) {
       index=i+j*N_rp;
       w[i]=w[i]+tmp_array[index];
    }
 for(i=0;i<N_rp;i++) { 
    w[i]=w[i]*binsize_pi*2.0;
    Arrayw[i][iM]=w[i];
 }

 free_vector(w,0,N_rp-1);
}

void read2proj_zero(double **Arrayw,int iM)
{
 extern double *tmp_array;
 extern int    N,k;
 int           i,j,index;
 double        *w,binsize_pi=2.0;

 for(i=0;i<N_rp;i++) Arrayw[i][iM]=0.0;

}



void read2multipole(FILE *fp,double **Array_Sf0,double **Array_Sf2,double **Array_Sf4,int iM)
{
 double xi0,xi2,xi4;
 double P0,P2,P4;
 double mu,dmu=0.05;
 int    i,j,index;
 extern double *tmp_array;
 extern int    N;
 double        **xi;

 xi=matrix(0,N_s-1,0,N_mu-1);

 //printf("N=%d\n",N);
 ftread(tmp_array,sizeof(double),N,fp);
 for(j=0;j<N_mu;j++)
    for(i=0;i<N_s;i++) {
       index=i+j*N_s;
       xi[i][j]=tmp_array[index];
    }

 for(i=0;i<N_s;i++) {
    xi0=0.0;
    xi2=0.0;
    xi4=0.0;
    for(j=0;j<N_mu;j++) {
       mu=(j+0.5)*dmu;
       P0=1.0;
       P2=(3.0*mu*mu-1.0)/2.0;
       P4=((35.0*mu*mu-30.0)*mu*mu+3.0)/8.0;
       xi0=xi0+(2.0*0+1.0)/2.0*P0*xi[i][j];
       xi2=xi2+(2.0*2+1.0)/2.0*P2*xi[i][j];
       xi4=xi4+(2.0*4+1.0)/2.0*P4*xi[i][j];
    }
    xi0=xi0*dmu;
    xi2=xi2*dmu;
    xi4=xi4*dmu;
    Array_Sf0[i][iM]=xi0;
    Array_Sf2[i][iM]=xi2;
    Array_Sf4[i][iM]=xi4;
 }

 free_matrix(xi,0,N_s-1,0,N_mu-1);

}



int main(int *argc,char **argv)
{
 FILE   *fp,*fpout;
 char   filename[128];
 float  n_halo;
 extern double *tmp_array;
 extern int    N;
 double *wf_cs,*wf_ss,*wp_cc,*wp_cs,*wp_sc,*wp_ss;
 double *numden_halo,*lgM,*lgMmin,*lgMmax,*lgMi,*lgMj;
 double *rp,*rpmin,*rpmax;

 double **Array_wf_cs,**Array_wf_ss;
 double **Array_wp_cc,**Array_wp_cs,**Array_wp_ss;
 double **Array_wpij_cc,**Array_wpij_cs,**Array_wpij_sc,**Array_wpij_ss;

 double **Array_f0_cs,**Array_f0_ss;
 double **Array_xi0_cc,**Array_xi0_cs,**Array_xi0_ss;
 double **Array_xiij0_cc,**Array_xiij0_cs,**Array_xiij0_sc,**Array_xiij0_ss;

 double **Array_f2_cs,**Array_f2_ss;
 double **Array_xi2_cc,**Array_xi2_cs,**Array_xi2_ss;
 double **Array_xiij2_cc,**Array_xiij2_cs,**Array_xiij2_sc,**Array_xiij2_ss;

 double **Array_f4_cs,**Array_f4_ss;
 double **Array_xi4_cc,**Array_xi4_cs,**Array_xi4_ss;
 double **Array_xiij4_cc,**Array_xiij4_cs,**Array_xiij4_sc,**Array_xiij4_ss;

 double *skip_array;

 int    i,j,index,iM,jM,itmp,Nskip;
 float  ftmp;
 double binsize_pi;
 char   str[10],strtmp[256];
 extern int k;

 int    Include_wp,Include_zspace;
 char   string[128];
 double h,Omega_m,Omega_L,Omega_b,n_s,sigma_8,redshift,boxsize,b;

 rp        = (double *) calloc(N_rp,sizeof(double));
 rpmin     = (double *) calloc(N_rp,sizeof(double));
 rpmax     = (double *) calloc(N_rp,sizeof(double));
 numden_halo=(double *) calloc(N_M, sizeof(double));
 lgM       = (double *) calloc(N_M, sizeof(double));
 lgMmin    = (double *) calloc(N_M, sizeof(double));
 lgMmax    = (double *) calloc(N_M, sizeof(double));
 lgMi      = (double *) calloc(N_M2,sizeof(double));
 lgMj      = (double *) calloc(N_M2,sizeof(double));

 Array_wf_cs=matrix(0,N_rp-1,0,N_M-1);
 Array_wf_ss=matrix(0,N_rp-1,0,N_M-1);
 Array_wp_cc=matrix(0,N_rp-1,0,N_M-1);
 Array_wp_cs=matrix(0,N_rp-1,0,N_M-1);
 Array_wp_ss=matrix(0,N_rp-1,0,N_M-1);
 Array_wpij_cc=matrix(0,N_rp-1,0,N_M2-1);
 Array_wpij_cs=matrix(0,N_rp-1,0,N_M2-1);
 Array_wpij_sc=matrix(0,N_rp-1,0,N_M2-1);
 Array_wpij_ss=matrix(0,N_rp-1,0,N_M2-1);

 strcpy(filename,"hist_halo_6500_0788.dat");
 if(!(fp=fopen(filename,"rt"))) FO_error(filename);
 fgets(strtmp,256,fp);
 for(iM=0;iM<N_M;iM++) {
    fscanf(fp,"%lf %lf %lf %d",&(lgM[iM]),&(lgMmin[iM]),&(lgMmax[iM]),&itmp);
//    printf("%f %f %f\n",lgM[iM],lgMmin[iM],lgMmax[iM]);
 }
 fclose(fp); 

 strcpy(filename,"rp.dat");
 if(!(fp=fopen(filename,"rt"))) FO_error(filename);
 for(i=0;i<N_rp;i++) {
    fscanf(fp,"%lf %lf %lf",&(rp[i]),&(rpmin[i]),&(rpmax[i]));
    printf("%f %f %f\n",rp[i],rpmin[i],rpmax[i]);
 }
 fclose(fp);

 N = N_rp * N_pi;
 tmp_array = (double *) calloc(N,   sizeof(double));

 for(iM=0;iM<N_M;iM++) {
    sprintf(str,"%03d%03d",iM+1,iM+1);
    strcpy(filename,"table/table_6500_0788_");
    strcat(filename,str);
    printf("%s\n",filename);

    if(!(fp=fopen(filename,"rb"))) FO_error(filename);

    ftread(&n_halo,sizeof(float),1,fp);
    //printf("n_halo=%e\n",n_halo);
    numden_halo[iM]=n_halo;

    read2proj(fp,Array_wf_cs,iM);
    read2proj(fp,Array_wf_ss,iM);
    read2proj(fp,Array_wp_cc,iM);
    read2proj(fp,Array_wp_cs,iM);
    ftread(tmp_array,sizeof(double),N,fp); /** array same as cs */
    read2proj(fp,Array_wp_ss,iM);

    fclose(fp);
 }
 /** Array[rp][M] done **/

 index=0;
 for(iM=0;iM<N_M;iM++) {
    for(jM=iM+1;jM<N_M;jM++) {
       // TODO: This is the part of determining how the table is organized.
       lgMi[index]=lgM[iM];
       lgMj[index]=lgM[jM];
       sprintf(str,"%03d%03d",iM+1,jM+1);
       strcpy(filename,"table/table_6500_0788_");
       strcat(filename,str);
       printf("%s\n",filename);

/*
       if(!(fp=fopen(filename,"rb"))) FO_error(filename);
       read2proj(fp,Array_wpij_cc,index);
       read2proj(fp,Array_wpij_cs,index);
       read2proj(fp,Array_wpij_sc,index);
       read2proj(fp,Array_wpij_ss,index);

       fclose(fp);
*/
       if(!(fp=fopen(filename,"rb"))) {
         read2proj_zero(Array_wpij_cc,index);
         read2proj_zero(Array_wpij_cs,index);
         read2proj_zero(Array_wpij_sc,index);
         read2proj_zero(Array_wpij_ss,index);
         printf("zeros for %s\n",filename);
       } else {
         read2proj(fp,Array_wpij_cc,index);
         read2proj(fp,Array_wpij_cs,index);
         read2proj(fp,Array_wpij_sc,index);
         read2proj(fp,Array_wpij_ss,index);
    
         fclose(fp);
       }
       index++;
    }
 }

 free(tmp_array);

 /** write to table **/
 strcpy(filename,"wptable_Jing16.bin");
 if(!(fpout=fopen(filename,"wb"))) FO_error(filename);

 Include_wp    =1;
 Include_zspace=1;
 strcpy(string,"Average of 16 Jing's Simulations, FOF b=0.2, z=0.512"); 
 h=0.71; Omega_m=0.268; Omega_L=0.732; Omega_b=0.045; n_s=1.00; sigma_8=0.85;
 redshift = 0.512; 
 boxsize=600.0; /* Mpc/h */ 
 b=0.20;
 //fwrite(&Include_wp,    sizeof(int),   1,fpout);
 //fwrite(&Include_zspace,sizeof(int),   1,fpout);
 fwrite(string,         sizeof(char),128,fpout);
 fwrite(&h,             sizeof(double),1,fpout);
 fwrite(&Omega_m,       sizeof(double),1,fpout);
 fwrite(&Omega_L,       sizeof(double),1,fpout);
 fwrite(&Omega_b,       sizeof(double),1,fpout);
 fwrite(&n_s,           sizeof(double),1,fpout);
 fwrite(&sigma_8,       sizeof(double),1,fpout);
 fwrite(&redshift,      sizeof(double),1,fpout);
 fwrite(&boxsize,       sizeof(double),1,fpout);
 fwrite(&b,             sizeof(double),1,fpout);

 i=N_rp; j=N_M; k=N_M*(N_M-1)/2;
 fwrite(&i,   sizeof(int),   1,fpout);
 fwrite(&j,   sizeof(int),   1,fpout);
 fwrite(&k,   sizeof(int),   1,fpout);

 fwrite(rp,   sizeof(double),N_rp,fpout);
 fwrite(rpmin,sizeof(double),N_rp,fpout);
 fwrite(rpmax,sizeof(double),N_rp,fpout);

 fwrite(lgM,   sizeof(double),N_M,fpout);
 fwrite(lgMmin,sizeof(double),N_M,fpout);
 fwrite(lgMmax,sizeof(double),N_M,fpout);

 fwrite(numden_halo,sizeof(double),N_M,fpout);

 for(i=0;i<N_rp;i++) fwrite(Array_wf_cs[i],  sizeof(double),N_M, fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wf_ss[i],  sizeof(double),N_M, fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wp_cc[i],  sizeof(double),N_M, fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wp_cs[i],  sizeof(double),N_M, fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wp_ss[i],  sizeof(double),N_M, fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wpij_cc[i],sizeof(double),N_M2,fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wpij_cs[i],sizeof(double),N_M2,fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wpij_sc[i],sizeof(double),N_M2,fpout);
 for(i=0;i<N_rp;i++) fwrite(Array_wpij_ss[i],sizeof(double),N_M2,fpout);

 fclose(fpout);

 for(i=0;i<N_rp;i++) {
    index=50;
    printf("%7.4f %13.4e %13.4e %13.4e %13.4e %13.4e\n",
      rp[i],Array_wf_cs[i][index],Array_wf_ss[i][index],Array_wp_cc[i][index],Array_wp_cs[i][index],Array_wp_ss[i][index]);
 }

 free(rp         );
 free(rpmin      );
 free(rpmax      );
 free(numden_halo);
 free(lgM        );
 free(lgMmin     );
 free(lgMmax     );
 free(lgMi       );
 free(lgMj       );

 free_matrix(Array_wf_cs,0,N_rp-1,0,N_M-1);
 free_matrix(Array_wf_ss,0,N_rp-1,0,N_M-1);
 free_matrix(Array_wp_cc,0,N_rp-1,0,N_M-1);
 free_matrix(Array_wp_cs,0,N_rp-1,0,N_M-1);
 free_matrix(Array_wp_ss,0,N_rp-1,0,N_M-1);

 free_matrix(Array_wpij_cc,0,N_rp-1,0,N_M2-1);
 free_matrix(Array_wpij_cs,0,N_rp-1,0,N_M2-1);
 free_matrix(Array_wpij_sc,0,N_rp-1,0,N_M2-1);
 free_matrix(Array_wpij_ss,0,N_rp-1,0,N_M2-1);

}
