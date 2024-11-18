#include <stdio.h>
#include <math.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum)
{
        int j;
        long k;
        static long iy=0;
        static long iv[NTAB];
        double temp;

        if (*idum <= 0 || !iy) {
                if (-(*idum) < 1) *idum=1;
                else *idum = -(*idum);
                for (j=NTAB+7;j>=0;j--) {
                        k=(*idum)/IQ;
                        *idum=IA*(*idum-k*IQ)-IR*k;
                        if (*idum < 0) *idum += IM;
                        if (j < NTAB) iv[j] = *idum;
                }
                iy=iv[0];
        }
        k=(*idum)/IQ;
        *idum=IA*(*idum-k*IQ)-IR*k;
        if (*idum < 0) *idum += IM;
        j=iy/NDIV;
        iy=iv[j];
        iv[j] = *idum;
        if ((temp=AM*iy) > RNMX) return RNMX;
        else return temp;
}
#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)


double ran2(long *idum)
{
  int j;
  long k;
  static long idum2=123456789;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    idum2=(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ1;
      *idum=IA1*(*idum-k*IQ1)-k*IR1;
      if (*idum < 0) *idum += IM1;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ1;
  *idum=IA1*(*idum-k*IQ1)-k*IR1;
  if (*idum < 0) *idum += IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2 += IM2;
  j=iy/NDIV;
  iy=iv[j]-idum2;
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX


double gasdev(long *idum)
{
	static int iset=0;
	static double gset;
	double fac,r,v1,v2;
//	float ran1();

	if  (iset == 0) {
		do {
			v1=2.0*ran2(idum)-1.0;
			v2=2.0*ran2(idum)-1.0;
			r=v1*v1+v2*v2;
		} while (r >= 1.0);
		fac=sqrt(-2.0*log(r)/r);
		gset=v1*fac;
		iset=1;
		return v2*fac;
	} else {
		iset=0;
		return gset;
	}
}


double expdev(long *idum)
{
  /*  double ran2(long *idum);*/
  double dum;
  
  do
    dum=ran2(idum);
  while (dum == 0.0);
  return -log(dum);
}

/*
  random number generator for angular distribution
  p(theta, phi) ~ (1+cos^2 theta) dcos_theta dphi, 
  -1<= cos_theta<=1
  0<=phi<2pi 
  getting p(theta) using rejection method: 
  comparation function is taken as 
  f(x)=2
*/

void get_angle(double *cos_thetap, double *phip, long *idum)
{
  double x, y;
  do{
  x=2*ran2(idum)-1.;
  y=2*ran2(idum);
  }
  while(y>1+x*x);
  *cos_thetap=x;
  *phip=2*M_PI*ran2(idum);
} 

void get_cos_theta(double *cos_theta, long *idum)
{ 
  double x, y;
  do{
     x=2*ran2(idum)-1.;
     y=2*ran2(idum);
    }
  while(y>1+x*x);
  *cos_theta=x;
}

void get_phi(double *phi, long *idum)
{
    *phi=2*M_PI*ran2(idum);
}

/**
 ** Use projection method to produce v_z according to its distribution.
 ** u_z=v_z/v_pT
 ** The comparison function is 1/(1+(u_z-...)**2/a0**2).
 ** The original function to comparison function ratio is exp(-u*u).
 ** When |xiLab|>~3.5, this method becomes somewhat slow.
 ** Test has been made to compare the cumulative distribution with the
 ** one obtained by direct integral.
 **/

void get_uz00(double *uz, double a, double x, long *idum)
{
  double u,v1,v2,r,y;
 
  do {
      do {
          v1=2.0*ran2(idum)-1.0;
          v2=2.0*ran2(idum)-1.0;
         }
      while (v1*v1+v2*v2 > 1.0);
      y=v2/v1;
      u=a*y+x;

      r=ran2(idum);
     }
  while(r > exp(-u*u));
 
  *uz=u;
}

void get_uz0(double *uz, double a, double x, long *idum)
{
  double u,r,d,t,t1,emd2,P,ratio,sign;

  sign=1.0;
  if(x<0) {x=-x; sign=-1.0;}
  d=x/(1.01+x/210.0+x*x/105.0); /*expirence*/
  emd2=exp(-d*d);
  t1=atan((d-x)/a);
  P=(t1+0.5*M_PI)/(0.5*M_PI*(1+emd2)+t1*(1-emd2));

  do {
      if(ran2(idum)<P)
        {
         t=-0.5*M_PI+ran2(idum)*(t1+0.5*M_PI);
         u=a*tan(t)+x;
        }
       else
        {
         t=t1+ran2(idum)*(0.5*M_PI-t1);
         u=a*tan(t)+x;
        }

       r=ran2(idum);
       ratio=exp(-u*u);
       if(u>=d) ratio=ratio/emd2;
     }
  while(r > ratio);

  *uz=sign*u;
}

void get_uz01(double *uz, double a, double x, long *idum)
{
  double u,v1,v2,r,y,b,t,t1,t2;

  b=5.7;
  t1=atan(-(b+x)/a);
  t2=atan((b-x)/a);

  if((x>-b)&&(x<b))
    { do {
        do {
            v1=2.0*ran2(idum)-1.0;
            v2=2.0*ran2(idum)-1.0;
           }
        while (v1*v1+v2*v2 > 1.0);
        y=v2/v1;
        u=a*y+x;

        r=ran2(idum);
         }
       while(r > exp(-u*u));
    }
  else 

    { do {
          t=t1+ran2(idum)*(t2-t1);
          u=a*tan(t)+x;
      
          r=ran2(idum); 
         }
      while(r > exp(-u*u));
    }

  *uz=u;
}

void get_up(double *up, long *idum)
{
  double u2;
  u2=expdev(idum);
  *up=sqrt(u2);
}

void get_up0(double *u_x, double *u_y, long *idum)
{
  double u2,up,phi;

  u2=expdev(idum);
  up=sqrt(u2);
  phi=2*M_PI*ran2(idum);
  
  *u_x=up*cos(phi);
  *u_y=up*sin(phi);

}

/* random number generator for Maxwell-Boltzman distribution, using rejection method. 
   The comparation function is simply a constant, the efficiency is 0.4254 (i.e. for each
   random number generated, 0.4254 of is not rejected */

/* The distribution generated here is 
   f(v)= sqrt(2/pi)/sigma^3 * v^2 exp(-0.5*v^2/sigma^2)
   whose integral from 0 to infinity is 1.
   maxwell(sigma, vmax) returns the desired random number, where sigma is the dispersion defined
   above, and vmax is the maxium velocity one want to generate. 
*/

static double xrandom(double xh, long *idum) /*return a uniform random number between 0 and xh. */
{
  return xh*ran2(idum);
}

static double g_x(double x)
{
  return x*x * exp(-0.5*x*x);
}

double maxwell(double sig, double vmax, long *idum)
{
  static double fmax;
  double f, fx, x;

  f = 1.0;
  fx = 0.0;
  while (f > fx) {
    x = xrandom(vmax/sig, idum);
    f = xrandom(fmax, idum);
    fx = g_x(x);
  }
  return sig*x;
}
  /* The number of points falling in a bin is given by
     10000N f(v) dv=M N f(v) vm/N=M f(v) vm, so 
     f(v)=Number of points in bin/(M vm), so what printed in
     tmp.dat is a list with v vs. f(v), where f(v) is normarlized, as defined
     in the begining part.
  */

/* maxwell_2d produce random numbers with distribution
   f(v) dv = N exp(-v^2/(sigma^2)) v dv
*/


double maxwell_2d(double sigma, long *idum)
{
  double u;
  u=expdev(idum);
  return sigma*sqrt(u);
}

#define PI 3.14159265358979323846264338327950288

double poidev(double xm, long *idum)
{
        double gammln(double xx);
        double ran1(long *idum);
        static double sq,alxm,g,oldm=(-1.0);
        double em,t,y;

        if (xm < 12.0) {
                if (xm != oldm) {
                        oldm=xm;
                        g=exp(-xm);
                }
                em = -1;
                t=1.0;
                do {
                        ++em;
                        t *= ran1(idum);
                } while (t > g);
        } else {
                if (xm != oldm) {
                        oldm=xm;
                        sq=sqrt(2.0*xm);
                        alxm=log(xm);
                        g=xm*alxm-gammln(xm+1.0);
                }
                do {
                        do {
                                y=tan(PI*ran1(idum));
                                em=sq*y+xm;
                        } while (em < 0.0);
                        em=floor(em);
                        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
                } while (ran1(idum) > t);
        }
        return em;
}
#undef PI

