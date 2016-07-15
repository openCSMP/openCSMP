#include "ErrorFunction.h"

#include <cmath>

using namespace std;

namespace csmp {

ErrorFunction::ErrorFunction()
 : ITMAX(100), 
   EPS(3.0e-7),
   FPMIN(1.0e-30)
 {
    cof[0] =  76.18009172947146;
    cof[1] = -86.50532032941677;
    cof[2] =  24.01409824083091;
    cof[3] =  -1.231739572450155;
    cof[4] =   0.1208650973866179e-2;
    cof[5] =  -0.5395239384953e-5;
 }


ErrorFunction::~ErrorFunction()
 {
 }
 

/// logarithm of gamma function [6.1]
double ErrorFunction::gammln( double xx )
{
	double x(xx), y(xx);
	double tmp(x+5.5);
	tmp -= (x+0.5)*log(tmp);
	double ser(1.000000000190015);
	for ( int j=0;j<=5;j++) ser += cof[j]/++y;
	return -tmp + log(2.5066282746310005*ser/x);
 }


/// continued fraction used by gammp and gammq [6.2]
void ErrorFunction::gcf( double& gammcf, double a, double x, double& gln )
{
	gln=gammln(a);
	
	double b=x+1.0-a;
	double c=1.0/FPMIN;
	double d=1.0/b;
	double h=d;
	double an, del;
	int    i;
	
	for ( i=1; i<=ITMAX; i++ ) {
		an = -i*(i-a);
		b += 2.0;
		d=an*d+b;
		if ( fabs(d) < FPMIN ) d=FPMIN;
		c=b+an/c;
		if ( fabs(c) < FPMIN ) c=FPMIN;
		d=1.0/d;
		del=d*c;
		h *= del;
		if ( fabs(del-1.0) < EPS ) break;
	}
	if (i > ITMAX) 
	  cout <<"\nErrorFunction::gcf a too large, ITMAX too small in gcf"<< endl;
	  
    gammcf = exp(-x+a*log(x)-gln) * h;
    
} // end gcf




/// series used by gammp and gammq [6.2]
void ErrorFunction::gseries( double& gamser, double a, double x, double& gln )
{
	gln=gammln(a);
	if (x <= 0.0) {
		 if (x < 0.0) 
		   cout <<"\ngseries: x is less than 0."<< endl;
		 gamser=0.0;
		 return;
	  } 
	else 
	  {
		double ap(a);
		double sum(1.0/a);
		double del(sum);
		
		for ( int n=1; n<=ITMAX; n++ ) 
		  {
			 ++ap;
			 del *= x/ap;
			 sum += del;
			 if (fabs(del) < fabs(sum)*EPS) {
			 	gamser=sum*exp(-x+a*log(x)-gln);
			 	return;
			 }
		  }
	  cout<<"\ngseries: a too large, ITMAX too small in routine gseries."<< endl;
      return;
	}
	
} // end gseries


#ifndef erf
double erf( double x ) 
 { 
    static ErrorFunction f; 
    return f.Erf(x); 
 } 
#endif

#ifndef erfc
double erfc( double x ) 
 { 
    static ErrorFunction f; 
    return f.Erfc(x); 
 } 
#endif

} // csp
