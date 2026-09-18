// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CUBIC_SPLINE_H
#define CUBIC_SPLINE_H

#include <vector>

namespace csmp {

/** Implementation of a cubic spline using at least 3 constraint points / data values and 2 derivatives at the curve ends.
    
    @note A sufficiently large number of the datapoints need to be supplied to get a satisfactory spline, dependent on the curve represented by the data.
    Graph the spline curve first before using the results in other computations. The data for graphing can be obtained using the Out() function.
    
    CubicSpine is used inside the experimental data-based saturation functions in the CSMP library found in the thmc directory.
*/
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
    
    double Value( const double x ) const;
    
    /// computes derivative using a forward finite difference approach @todo improve by using analytic derivative
    double Derivative( const double x ) const;
    
    /// computes the maximum value of the second derivative of the spline function; use for checking whether there are enough data points to get satisfactory curve
    double MaxDerivative() const;
    double Range_x() const;
    double Range_fx() const;
    
    void Out() const;
    
  private:
    std::vector<double> xa_, ya_, y2a_;
    double  x_range_, y_range_, xa_min, xa_max, ya_min, ya_max;
};

/// function computes second derivative that is needed to define the spline; is called only once during spline construction
double splint( const std::vector<double>& xa, 
               const std::vector<double>& ya,
               const std::vector<double>& y2a,
               double x );
                      
 } // csmp

#endif
