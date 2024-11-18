#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "typedef.h"
#include "mathutil.h"

double Nc_avg(double lgM,HOD HODpar)
{ double lgMmin,sig_lgM,Nc,x;

  lgMmin  = HODpar.lgMmin;
  sig_lgM = HODpar.sig_lgM;

  x=(lgM-lgMmin)/sig_lgM;
  Nc=0.5*(1.0+erf(x));
//  if(x<=-6.0) Nc=0.0;
//  if(x>=4.0) Nc=1.0;
//  else Nc=0.5*(1.0+erf(x));
 
//  printf("lgM=%f x=%f Nc=%f\n",lgM,x,Nc); 

  return Nc;
}

/** Ns before multiplying the Ncen profile **/
double Ns_avg(double lgM,HOD HODpar)
{ double M,lgMmin,sig_M,lgNum,M0,M1p,alpha,Ns,f,x;

  M     = pow(10.0,lgM);
  M0    = pow(10.0,HODpar.lgM0);
  M1p   = pow(10.0,HODpar.lgM1p);
  alpha = HODpar.alpha;

  x=(M-M0)/M1p;
  if(x>0.0) Ns=pow(x,alpha);
  else      Ns=0.0;

  return Ns;
}

/** Ns after multiplying the Ncen profile **/
double Ns_avg_f(double lgM,HOD HODpar)
{ double M,lgMmin,sig_M,lgNum,M0,M1p,alpha,Ns,f,x;

  return Nc_avg(lgM,HODpar)*Ns_avg(lgM,HODpar);
}

double N_avg(double lgM,HOD HODpar)
{ double Nc,Ns;

  Nc=Nc_avg(lgM,HODpar);
  Ns=Ns_avg(lgM,HODpar);
 
  return Nc+Ns;
}


