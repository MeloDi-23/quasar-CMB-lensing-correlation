#ifndef RANDOM_H
#define RANDOM_H
double ran1(long *idum);
double ran2(long *idum);
double expdev(long *idum);
double gasdev(long *idum);
void get_angle(double *cos_thetap, double *phip, long *idum);
void get_cos_theta(double *cos_theta, long *idum);
void get_phi(double *phi, long *idum);
void get_uz0(double *uz, double a, double x, long *idum);
void get_uz(double *uz, double aT, double betapT, double xiLab, long *idum);
void get_uz2(double *uz, double aT, double betapT, double xiLab, long *idum);
void get_uz3(double *uz, double aT, double betapT, double xiLab, long *idum);
void get_up(double *up, long *idum);
void get_up0(double *u_x, double *u_y, long *idum);
double maxwell_2d(double sig, long *idum);
double poidev(double xm, long *idum);

#endif
