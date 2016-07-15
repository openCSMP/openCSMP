#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include "CSMP_definitions.h"

namespace csmp {

/// implementation of a cubic splite using 3 constraint points and 2 derivatives at the curve ends
class CubicSpline {
  public:
    CubicSpline();
    CubicSpline( const CubicSpline& );
    explicit CubicSpline( const char* datafile );
    ~CubicSpline();
    CubicSpline& operator=( const CubicSpline& );
    void Initialize( const char* datafile );
    void Initialize( const std::vector<double64>&, const std::vector<double64>&,
                     const double64, const double64 );
    
    double64 Value( double64 x ) const;
    double64 Derivative( double64 x ) const;
    
    // of the input values
    double64 MaxDerivative() const;
    double64 Range_x() const;
    double64 Range_fx() const;
    
    void Out() const;
    
  private:
    std::vector<double64> xa_, ya_, y2a_;
    double64  x_range_, y_range_, xa_min, xa_max, ya_min, ya_max;
};

double64 splint( const std::vector<double64>& xa, 
                 const std::vector<double64>& ya,
                 const std::vector<double64>& y2a,
                 double64 x );
                      


inline double64 CubicSpline::Value( double64 x ) const
 {
    return splint( xa_, ya_, y2a_, x );
 } 

/// two point derivative about x(x-1% of range,x+1% of range
inline double64 CubicSpline::Derivative( double64 x ) const
 {
    return (Value(x+x_range_/100.)-Value(x-x_range_/100.)) / (x_range_/50.);
 } 
 

inline double64 CubicSpline::MaxDerivative() const
 {
    return (*std::max_element( y2a_.begin(), y2a_.end() ));
 } 


inline double64 CubicSpline::Range_x() const
 {
    return x_range_;
 } 


inline double64 CubicSpline::Range_fx() const
 {
    return y_range_;
 } 
 
 } // csmp

#endif
