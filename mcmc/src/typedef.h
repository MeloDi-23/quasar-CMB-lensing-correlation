#ifndef TYPEDEF_H
#define TYPEDEF_H

/** cosmological model **/
typedef struct {
  double Omega_m;   /** Matter **/
  double Omega_L;   /** Lambda **/
  double Omega_b;   /** baryon **/
  double n_s;       /** index of inflationary power spectrum **/
  double h;         /** Hubble constant = 100h km/s/Mpc **/
  double sig8;      /** rms mass fluctuation at scale 8Mpc/h at z=0 **/
  double redshift;  /** redshift **/
  double M_starz;   /** Nonlinear mass under the given cosmology at z   **/
  double M_star0;   /** Nonlinear mass under the given cosmology at z=0 **/
} COSMOS;

/** HOD parameters **/
typedef struct {
  double lgMmin;
  double sig_lgM;
  double lgM0;
  double lgM1p;
  double alpha;
  int    icen;           /** include central galaxy? icen=0 No; icen!=0 Yes **/
  int    iFreePar[6];    /** is the i-th parameter free? 1=yes, 0=no   **/

  double ng_avg;
  double err_ng;

  int    N_M;
  double *Nc;
  double *Ns;

} HOD;

/** cf table **/
typedef struct {
  char   string[128];
  double h,Omega_m,Omega_L,Omega_b,n_s,sigma_8,redshift,boxsize,b;
  int    N_rp,N_M,N1,N2;
  double *rp,*rpmin,*rpmax;
  double *lgM,*lgMmin,*lgMmax;
  double *nh;
  double *wf_cs,*wf_ss;
  double *wp_cc,*wp_cs,*wp_ss;
  double *wpij_cc,*wpij_cs,*wpij_sc,*wpij_ss;
} table_struct;

#endif
