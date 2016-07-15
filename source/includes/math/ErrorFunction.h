#ifndef CSMP_ERROR_FUNCTION_H
#define CSMP_ERROR_FUNCTION_H

#include "CSMP_definitions.h"

namespace csmp {


class ErrorFunction {
  public:
    ErrorFunction();
    ~ErrorFunction();
    
    double Erf( double x );
    double Erfc( double x );
    
  private:
    double gammln( double xx ); 
    void   gcf( double& gammcf, double a, double x, double& gln );
    void   gseries( double& gamser, double a, double x, double& gln ); 
    double gammp( double a, double x );
    double gammq( double a, double x );
    
  private:
    const int     ITMAX;
    const double  EPS,
                  FPMIN;
	double  cof[6];
};

#ifndef ERF
double erf( double x ); 
#endif

#ifndef ERFC
double erfc( double x ); 
#endif

// incomplete gamma function [6.2] 
inline double ErrorFunction::gammp( double a, double x )
{
	if ( x < 0. || a <= 0. ) {
	     std::cout <<"\ngammp: Invalid arguments in routine gammp."<< std::endl;
	     throw std::domain_error("ErrorFunction::gammp");
	  }
	  
	double gamser, gammcf, gln;
   
	if (x < (a+1.0)) {
		gseries( gamser, a, x, gln );
		return gamser;
	} else {
		gcf( gammcf, a, x, gln );
		return 1.0-gammcf;
	}
}

// complement of incomplete gamma function [6.2] 
inline double ErrorFunction::gammq( double a, double x )
{
	double gamser, gammcf, gln;
   
	if ( x < 0. || a <= 0. ) {
	     std::cout <<"\ngammq: Invalid arguments in routine gammq."<< std::endl;
	     throw std::domain_error("ErrorFunction::gammq");
	  }
	  
	if (x < (a+1.0)) {
		gseries( gamser, a, x, gln );
		return 1.0-gamser;
	} else {
		gcf( gammcf, a, x, gln );
		return gammcf;
	}
}


inline double ErrorFunction::Erf( double x )
{
   return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}



inline double ErrorFunction::Erfc( double x )
{
   return x < 0.0 ? 1.0+gammp(0.5,x*x) : gammq(0.5,x*x);
}


} // csmp

#endif
