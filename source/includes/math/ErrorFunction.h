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

} // csmp

#endif
