#ifndef BROOKS_COREY_FRONT_VELOCITY_H
#define BROOKS_COREY_FRONT_VELOCITY_H

#include "CSMP_definitions.h"

namespace csmp {

/// @author Insa Neuweiler @date 2002
class BrooksCoreyFrontVelocity {
  public:
    BrooksCoreyFrontVelocity();
    ~BrooksCoreyFrontVelocity();
    
    double ShockVelocityMultiplier( double lambda, double visc_ratio ); 
  
  private:
    double dffbc(double x, double lambda, double mu);
    double fxbc(double x, double lambda, double mu);
    double dfxbc(double x, double lambda, double mu);
    void funcdbc(double x, double *fl, double *df, double lambda, double mu);
    void zbrakbc(double x1,double x2,int n,double *xb1,double *xb2, double lambda, double mu);
    double findrootbc(double x1,double x2,double xacc,double lambda,double mu);
    const int MAXIT;
};

} // end namespace csmp

#endif
