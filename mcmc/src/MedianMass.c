#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "constants.h"
#include "typedef.h"
#include "hod.h"
#include "halo.h"
#include "qromo.h"
#include "mathutil.h"
#include "pspectrum.h"

COSMOS cosmology;
HOD HODpar;

double MedianMass(COSMOS cosmology, HOD HODpar, int censat)
{
 int    Nbin=2000,i;
 double lnMmin,lnMmax,dlnM,lnMass,dndlnMass,Navg,Mass;
 double S[2000],halfS,ln10=2.30258509299405;
 
 lnMmin=log(MhaloMin);
 lnMmax=log(MhaloMax);
 dlnM=(lnMmax-lnMmin)/(Nbin-1); 
 
 for(i=0;i<Nbin;i++) {
    lnMass=lnMmin+i*dlnM;
    splint(cosmology.lndndlnM_lnM.x,cosmology.lndndlnM_lnM.y,cosmology.lndndlnM_lnM.y2,cosmology.lndndlnM_lnM.n,lnMass,&dndlnMass);
    dndlnMass=exp(dndlnMass);
    Mass=exp(lnMass);

    if(censat==1) Navg=Nc_avg(Mass,HODpar);
    if(censat==0) Navg=Ns_avg(Mass,HODpar);

    if(i==0) S[i]=dndlnMass*Navg; 
    else     S[i]=S[i-1]+dndlnMass*Navg;
 }

 halfS=0.5*S[Nbin-1]; 
 for(i=0;i<Nbin;i++) {
    if(S[i]>halfS) break;
 }
 lnMass=lnMmin+(i-1+(halfS-S[i-1])/(S[i]-S[i-1]))*dlnM; 
 /** linear interpolation to find the mass **/

 return lnMass/ln10;
}

