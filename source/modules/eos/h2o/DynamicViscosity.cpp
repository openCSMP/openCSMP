#include "DynamicViscosity.h"

using namespace std;

namespace csmp {

DynamicViscosity::DynamicViscosity()
 : kelvin(273.15),
   tstar(647.27), 
   rhostar(317.763)
 {       
    ak[0] =  0.0181583;
    ak[1] =  0.0177624;
    ak[2] =  0.0105287;
    ak[3] = -0.0036744;

    bij[0][0] =  0.501938;
    bij[1][0] =  0.162888;
    bij[2][0] = -0.130356;
    bij[3][0] =  0.907919;
    bij[4][0] = -0.551119;
    bij[5][0] =  0.146543;
    bij[0][1] =  0.235622;
    bij[1][1] =  0.789393;
    bij[2][1] =  0.673665;
    bij[3][1] =  1.207552;
    bij[4][1] =  0.0670665;
    bij[5][1] = -0.084337;
    bij[0][2] = -0.274637;
    bij[1][2] = -0.743539;
    bij[2][2] = -0.959456;
    bij[3][2] = -0.687343;
    bij[4][2] = -0.497089;
    bij[5][2] =  0.195286;
    bij[0][3] =  0.145831;
    bij[1][3] =  0.263129;
    bij[2][3] =  0.347247;
    bij[3][3] =  0.213486;
    bij[4][3] =  0.100754;
    bij[5][3] = -0.032932;
    bij[0][4] = -0.0270448;
    bij[1][4] = -0.0253093;
    bij[2][4] = -0.0267758;
    bij[3][4] = -0.0822904;
    bij[4][4] =  0.0602253;
    bij[5][4] = -0.0202595;
    
 }




DynamicViscosity::~DynamicViscosity()
 { 
 }



double64 DynamicViscosity::ViscosityFromTemperatureAndDensity( double64 t, double64 rho )
 {
   double64 trat,trat1,rhorat,rhorat1,n0,n;
   int i,j,k;

   trat    = (t+kelvin)/tstar;
   trat1   = 1.0/trat-1.0;
   rhorat  = rho/rhostar;
   rhorat1 = rhorat-1.0;

   n0=0;
   for(k=0;k<4;++k){
     n0 += (ak[k]*pow((1.0/trat),k));
   }  
   n0 = 1.0/n0;
   n0 *= (sqrt(trat));

   n=0;
   for(i=0;i<6;++i){
     for(j=0;j<5;++j){
       n += (bij[i][j]*pow(trat1,i)*pow(rhorat1,j));
     }
   }
   n = n0*exp(rhorat*n);

   return n*1.0e-06; // converts to kg/s/m
}
 
 
} 
 

































