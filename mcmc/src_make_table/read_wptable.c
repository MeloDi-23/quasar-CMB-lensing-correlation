#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef struct {
  int    Include_wp,Include_zspace;
  char   string[128];
  double h,Omega_m,Omega_L,Omega_b,n_s,sigma_8,redshift,boxsize,b;
  int    N_rp,N_M,N_s,N1,N2;
  double *rp,*rpmin,*rpmax;
  double *lgM,*lgMmin,*lgMmax;
  double *nh;

  double *wf_cs,*wf_ss;
  double *wp_cc,*wp_cs,*wp_ss;
  double *wpij_cc,*wpij_cs,*wpij_sc,*wpij_ss;

  double *f0_cs,*f0_ss;
  double *xi0_cc,*xi0_cs,*xi0_ss;
  double *xiij0_cc,*xiij0_cs,*xiij0_sc,*xiij0_ss;

  double *f2_cs,*f2_ss;
  double *xi2_cc,*xi2_cs,*xi2_ss;
  double *xiij2_cc,*xiij2_cs,*xiij2_sc,*xiij2_ss;

  double *f4_cs,*f4_ss;
  double *xi4_cc,*xi4_cs,*xi4_ss;
  double *xiij4_cc,*xiij4_cs,*xiij4_sc,*xiij4_ss;

} table_struct;

typedef struct {
 double Mmin;
 double sig_lgM;
 double M0;
 double M1p;
 double alpha;
} HOD;

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

void load_table(table_struct *table,char *filename) 
{
 FILE *fp;
 int    i,index,N_M,N_M2,N_rp,N1,N2;
 int    Include_wp,Include_zspace;
 char   string[128];
 double h,Omega_m,Omega_L,Omega_b,n_s,sigma_8,redshift,boxsize,b;

 if(!(fp=fopen(filename,"rb"))) FO_error(filename);

 //fread(&Include_wp,    sizeof(int),   1,fp);
 //fread(&Include_zspace,sizeof(int),   1,fp);
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

 //printf("%d %d\n",Include_wp,Include_zspace);
 printf("%s\n",string);
 printf("%f %f %f %f %f %f %f %f %f\n",h,Omega_m,Omega_L,Omega_b,n_s,sigma_8,redshift,boxsize,b);

 //(*table).Include_wp    =Include_wp;
 //(*table).Include_zspace=Include_zspace;
 strcpy((*table).string,string);
 (*table).h             =h;
 (*table).Omega_m       =Omega_m;
 (*table).Omega_L       =Omega_L;
 (*table).Omega_b       =Omega_b;
 (*table).n_s           =n_s;
 (*table).sigma_8       =sigma_8;
 (*table).redshift      =redshift;
 (*table).boxsize       =boxsize;
 (*table).b             =b;
 fread(&N_rp,sizeof(int),1,fp);
 fread(&N_M, sizeof(int),1,fp);
 fread(&N_M2,sizeof(int),1,fp);
 N1=N_rp*N_M;
 N2=N_rp*N_M2;

 (*table).N_rp=N_rp;
 (*table).N_s =N_rp;
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

//void get_2pcf(double *F,double *Q,int N_obs,int N_deq,table_struct table,HOD HODpar)
//{
// double *Nc,*Ns,*NcNs,*NsNsm1;
 /** the definition of Nc, ... should be out side of this function **/ 
//}
int main()
{
 char   filename[128];
 FILE   *fp;
 int    i,index,N_M,N_M2,N_rp,N1,N2;
 int    ir,j;
 table_struct table;
 double cumu_n;

 strcpy(filename,"wptable.bin");

 load_table(&table,filename);
 N_rp=table.N_rp;
 N_M =table.N_M;

 printf("lgM=%f\n",table.lgM[166]);
 for(i=0;i<N_rp;i++) {
    index=i*N_M+166;
    printf("%7.4f %13.4e %13.4e %13.4e %13.4e %13.4e\n",
      table.rp[i],table.wf_cs[index],table.wf_ss[index],table.wp_cc[index],table.wp_cs[index],table.wp_ss[index]);
 }

 ir=N_rp-1;
 index=0;
 for(i=0;i<N_M;i++)
    for(j=i+1;j<N_M;j++) {
       printf("%3d %3d %5d %e %e %e %e %e\n",i+1,j+1,index,table.nh[i],table.nh[j],table.wpij_cc[index],table.wpij_cs[index],table.wpij_sc[index],table.wpij_ss[index]);
       index++;
    }

 fp=fopen("cumu_mf_table.dat","wt");
 cumu_n=0.0;
 for(i=N_M-1;i>=0;i--) {
    cumu_n=cumu_n+table.nh[i];
    fprintf(fp,"%f %e\n",table.lgMmin[i],cumu_n);
 }
 fclose(fp);
 
 free_table(&table);
}


