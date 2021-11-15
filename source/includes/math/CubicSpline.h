#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "CSMP_definitions.h"

namespace csmp {

/// implementation of a cubic spline using 3 constraint points and 2 derivatives at the curve ends
class CubicSpline {
  public:
    CubicSpline();
    CubicSpline( const CubicSpline& );
    explicit CubicSpline( const char* datafile );
    ~CubicSpline();
    CubicSpline& operator=( const CubicSpline& );
    void Initialize( const char* datafile );
    void Initialize( const std::vector<double>&,
                     const std::vector<double>&,
                     const double, const double );
    
    double Value( double x ) const;
    double Derivative( double x ) const;
    
    // of the input values
    double MaxDerivative() const;
    double Range_x() const;
    double Range_fx() const;
    
    void Out() const;
    
  private:
    std::vector<double> xa_, ya_, y2a_;
    double  x_range_, y_range_, xa_min, xa_max, ya_min, ya_max;
};

double splint( const std::vector<double>& xa, 
                 const std::vector<double>& ya,
                 const std::vector<double>& y2a,
                 double x );
                      
 } // csmp

#endif
