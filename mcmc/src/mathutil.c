#include <math.h>
#include "complex.h"

double gammln(xx)
double xx;
{
        double x,tmp,ser;
        static double cof[6]={76.18009173,-86.50532033,24.01409822,
                -1.231739516,0.120858003e-2,-0.536382e-5};
        int j;

        x=xx-1.0;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.0;
        for (j=0;j<=5;j++) {
                x += 1.0;
                ser += cof[j]/x;
        }
        return -tmp+log(2.50662827465*ser);
}

double factln(n)
int n;
{
        static double a[101];
        double gammln();
        void nrerror();

//        printf("factln> n=%d factln=%e\n",n,gammln(n+1.0));
        if (n < 0) nrerror("Negative factorial in routine FACTLN");
        if (n <= 1) return 0.0;
        if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
        else return gammln(n+1.0);
}

double bico(n,k)
int n,k;
{
        double factln();

        return floor(0.5+exp(factln(n)-factln(k)-factln(n-k)));
}

double lnbico(n,k)
int n,k;
{
        double factln();

        return (factln(n)-factln(k)-factln(n-k));
}

double beta(z,w)
double z,w;
{
        double gammln();

        return exp(gammln(z)+gammln(w)-gammln(z+w));
}

/*** modified -- Beyond the ends of x, linear interpolation is used ***/
void spline(x,y,n,yp1,ypn,y2)
double x[],y[],yp1,ypn,y2[];
int n;
{
	int i,k;
	double p,qn,sig,un,*u,*vector();
	void free_vector();

	u=vector(1,n-1);
	if (yp1 > 0.99e30)
		y2[1]=u[1]=0.0;
	else {
		y2[1] = -0.5;
		u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
	}
	for (i=2;i<=n-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
//                printf("spline loop> i=%d x=%e y=%e u=%e y2=%e\n",i,x[i],y[i],u[i],y2[i]);  
//                printf("spline loop+1> i+1=%d x=%e y=%e u=%e y2=%e\n",i+1,x[i+1],y[i+1],u[i+1],y2[i+1]);  
 
	}
	if (ypn > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
	}
	y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
//        printf("spline loop1> un qn u[n-1]=%e %e %e i=%d y2=%e\n",un,qn,u[n-1],n,y2[n]);
	for (k=n-1;k>=1;k--)
		{ y2[k]=y2[k]*y2[k+1]+u[k];
//                   printf("spline loop2> k=%d y2=%e\n",k,y2[k]);
                }
	free_vector(u,1,n-1);

//        for(k=1;k<=n;k++)
//            printf("spline k x y y2>%d %e %e %e\n",k,x[k],y[k],y2[k]);
}


void splint(xa,ya,y2a,n,x,y)
double xa[],ya[],y2a[],x,*y;
int n;
{
	int klo,khi,k;
	double h,b,a;
	void nrerror();

	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad XA input to routine SPLINT");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;
        b=(x-xa[klo])/h;
        /** linear interpolation outside boundary **/
        if(klo==1 && b<0.0)
             *y=a*ya[klo]+b*ya[khi]-b*h*h*(2.0*y2a[klo]+y2a[khi])/6.0;
        else if(khi==n && a<0.0)
             *y=a*ya[klo]+b*ya[khi]-a*h*h*(y2a[klo]+2.0*y2a[khi])/6.0;
        else
             *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
//        printf("splint> x=%e y=%e\n",x,*y);
//        for(k=1;k<=n;k++)
//            printf("splint k x y y2>%d %e %e %e\n",k,xa[k],ya[k],y2a[k]); 
}

/*** do spline fitting / interpolation only in the range 
     from x[n1] to x[n2] 
***/
void part_spline(x,y,n,n1,n2,ypn1,ypn2,y2)
double x[],y[],ypn1,ypn2,y2[];
int n,n1,n2;
{
	int i,k;
	double p,qn,sig,un,*u,*vector();
	void free_vector();
 

        if(n1<1 || n1>n) n1=1;
        if(n2>n || n2<1) n2=n;
        if(n1==n2) {
          n1=n1-1; 
          n2=n2+1;
          if(n1<1) n1=1;
          if(n2>n) n2=n;
        } 
        for(i=1;i<=n;i++)
           y2[i]=0.0;

	u=vector(1,n-1);
	if (ypn1 > 0.99e30)
		y2[n1]=u[n1]=0.0;
	else {
		y2[n1] = -0.5;
		u[n1]=(3.0/(x[n1+1]-x[n1]))*((y[n1+1]-y[n1])/(x[n2]-x[n1])-ypn1);
	}
	for (i=n1+1;i<=n2-1;i++) {
		sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p=sig*y2[i-1]+2.0;
		y2[i]=(sig-1.0)/p;
		u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
//                printf("spline loop> i=%d x=%e y=%e u=%e y2=%e\n",i,x[i],y[i],u[i],y2[i]);  
//                printf("spline loop+1> i+1=%d x=%e y=%e u=%e y2=%e\n",i+1,x[i+1],y[i+1],u[i+1],y2[i+1]);  
 
	}
	if (ypn2 > 0.99e30)
		qn=un=0.0;
	else {
		qn=0.5;
		un=(3.0/(x[n2]-x[n2-1]))*(ypn2-(y[n2]-y[n2-1])/(x[n2]-x[n2-1]));
	}
	y2[n2]=(un-qn*u[n2-1])/(qn*y2[n2-1]+1.0);
//        printf("spline loop1> un qn u[n-1]=%e %e %e i=%d y2=%e\n",un,qn,u[n-1],n,y2[n]);
	for (k=n2-1;k>=n1;k--)
		{ y2[k]=y2[k]*y2[k+1]+u[k];
//                   printf("spline loop2> k=%d y2=%e\n",k,y2[k]);
                }
	free_vector(u,1,n-1);

//        for(k=1;k<=n;k++)
//            printf("spline k x y y2>%d %e %e %e\n",k,x[k],y[k],y2[k]);
}


/** iexp=1 use exp(y) for interpolation out of the range **/
void part_splint(xa,ya,y2a,n,n1,n2,x,y,iexp)
double xa[],ya[],y2a[],x,*y;
int n,n1,n2,iexp;
{
	int klo,khi,k;
	double h,b,a;
	void nrerror();
       
        if(n1<1 || n1>n) n1=1;
        if(n2>n || n2<1) n2=n;
 
	klo=1;
	khi=n;
	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;
	}
	h=xa[khi]-xa[klo];
	if (h == 0.0) nrerror("Bad XA input to routine PART_SPLINT");
        if((n1>1 && khi<=n1)  || (n2<n && klo>=n2))
          {
            a=(xa[khi]-x)/h;
            b=(x-xa[klo])/h;
            if(iexp==0) *y=a*ya[klo]+b*ya[khi];
            else { *y=a*exp(ya[klo])+b*exp(ya[khi]);
                   if(*y<0.0) *y=-99.0;
                   else *y=log(*y);
                 }
//            printf("x=%f %f %f %f\n",x,ya[klo],*y,ya[khi]);
            /** out of spline range [n1,n2] use linear interpolation 
                if iexp!=0 use exp(y) for interpolation
                but in the range [1,n] still spline3 
             **/
          }
        else
          { a=(xa[khi]-x)/h;
            b=(x-xa[klo])/h;
            /** linear interpolation outside boundary **/
            if(klo==1 && b<0.0)
                 *y=a*ya[klo]+b*ya[khi]-b*h*h*(2.0*y2a[klo]+y2a[khi])/6.0;
            else if(khi==n && a<0.0)
                 *y=a*ya[klo]+b*ya[khi]-a*h*h*(y2a[klo]+2.0*y2a[khi])/6.0;
            else
                 *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
          } 
//            if((klo==1 && b<0.0) || (khi==n && a<0.0)) h=0.0;
//            *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;     }
//        printf("splint> x=%e y=%e\n",x,*y);
//        for(k=1;k<=n;k++)
//            printf("splint k x y y2>%d %e %e %e\n",k,xa[k],ya[k],y2a[k]); 
}

#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}

void gaussj(a,n,b,m)
double **a,**b;
int n,m;
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll,*ivector();
	double big,dum,pivinv;
	void nrerror(),free_ivector();

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

void gaussjinv(a,n)
double **a;
int n;
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll,*ivector();
	double big,dum,pivinv;
	void nrerror(),free_ivector();

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1) nrerror("GAUSSJ: Singular Matrix-1");
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) nrerror("GAUSSJ: Singular Matrix-2");
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}

#undef SWAP

#include <math.h>

#define ITMAX 100
#define EPS 3.0e-8

double zbrent(func,x1,x2,tol)
double x1,x2,tol;
double (*func)();	/* ANSI: float (*func)(float); */
{
	int iter;
	double a=x1,b=x2,c,d,e,min1,min2;
	double fa=(*func)(a),fb=(*func)(b),fc,p,q,r,s,tol1,xm;
	void nrerror();

	if (fb*fa > 0.0) nrerror("Root must be bracketed in ZBRENT");
	fc=fb;
	for (iter=1;iter<=ITMAX;iter++) {
		if (fb*fc > 0.0) {
			c=a;
			fc=fa;
			e=d=b-a;
		}
		if (fabs(fc) < fabs(fb)) {
			a=b;
			b=c;
			c=a;
			fa=fb;
			fb=fc;
			fc=fa;
		}
		tol1=2.0*EPS*fabs(b)+0.5*tol;
		xm=0.5*(c-b);
		if (fabs(xm) <= tol1 || fb == 0.0) return b;
		if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
			s=fb/fa;
			if (a == c) {
				p=2.0*xm*s;
				q=1.0-s;
			} else {
				q=fa/fc;
				r=fb/fc;
				p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
				q=(q-1.0)*(r-1.0)*(s-1.0);
			}
			if (p > 0.0)  q = -q;
			p=fabs(p);
			min1=3.0*xm*q-fabs(tol1*q);
			min2=fabs(e*q);
			if (2.0*p < (min1 < min2 ? min1 : min2)) {
				e=d;
				d=p/q;
			} else {
				d=xm;
				e=d;
			}
		} else {
			d=xm;
			e=d;
		}
		a=b;
		fa=fb;
		if (fabs(d) > tol1)
			b += d;
		else
			b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
		fb=(*func)(b);
	}
	nrerror("Maximum number of iterations exceeded in ZBRENT");
}

/**** Bilinear interpolation ****/
/**** x[] and y[] are required to be monototically increasing ****/
/***
      * z(x[i],y[j+1])         * z(x[i+1],y[j+1])

                * z(xt,yt)

      * z(x[i],y[j])           * z(x[i+1],y[j])
***/

void bilinint(double *x,double *y,double **z,int NDIMx,int NDIMy,double xt,double yt,double *zt)
{
  int indx,indy;
  int i,j;
  double t,u;

/*** find the box region that (xt,yt) locates ***/
  if(xt<x[1]) i=1;
  else if (xt>x[NDIMx]) i=NDIMx-1;
  else { i=2;
         while(xt>x[i]) i=i+1;
         i=i-1;
       }

  if(yt<y[1]) j=1;
  else if (yt>y[NDIMy]) j=NDIMy-1;
  else { j=2;
         while(yt>y[j]) j=j+1;
         j=j-1;
       }

/*** interpolation ***/
  t=(xt-x[i])/(x[i+1]-x[i]);
  u=(yt-y[j])/(y[j+1]-y[j]);

  *zt=(1.0-t)*(1.0-u)*z[i][j]+t*(1.0-u)*z[i+1][j]+(1.0-t)*u*z[i][j+1]+t*u*z[i+1][j+1];

}

void bilinint3(double *x,double *y,double **zb,double **zc,double **zd,int NDIMx,int NDIMy,double xt,double yt,double *zbt,double *zct,double *zdt)
{
  int indx,indy;
  int i,j;
  double t,u;

/*** find the box region that (xt,yt) locates ***/
  if(xt<x[1]) i=1;
  else if (xt>x[NDIMx]) i=NDIMx-1;
  else { i=2;
         while(xt>x[i]) i=i+1;
         i=i-1;
       }

  if(yt<y[1]) j=1;
  else if (yt>y[NDIMy]) j=NDIMy-1;
  else { j=2;
         while(yt>y[j]) j=j+1;
         j=j-1;
       }

/*** interpolation ***/
  t=(xt-x[i])/(x[i+1]-x[i]);
  u=(yt-y[j])/(y[j+1]-y[j]);

  *zbt=(1.0-t)*(1.0-u)*zb[i][j]+t*(1.0-u)*zb[i+1][j]+(1.0-t)*u*zb[i][j+1]+t*u*zb[i+1][j+1];
  *zct=(1.0-t)*(1.0-u)*zc[i][j]+t*(1.0-u)*zc[i+1][j]+(1.0-t)*u*zc[i][j+1]+t*u*zc[i+1][j+1];
  *zdt=(1.0-t)*(1.0-u)*zd[i][j]+t*(1.0-u)*zd[i+1][j]+(1.0-t)*u*zd[i][j+1]+t*u*zd[i+1][j+1];

}


