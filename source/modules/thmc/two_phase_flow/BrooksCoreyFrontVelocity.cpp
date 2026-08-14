#include "BrooksCoreyFrontVelocity.h"

using namespace std;

namespace csmp {

BrooksCoreyFrontVelocity::BrooksCoreyFrontVelocity()
 : MAXIT(100)
 {
 }


// destructor
BrooksCoreyFrontVelocity::~BrooksCoreyFrontVelocity()
 {
 }


double BrooksCoreyFrontVelocity::dffbc(double x, double lambda, double mu)
{
  double dfrac,a1,a2,a3;


  a1= pow( 1-x, 2/lambda );
  a2= pow( 1-x, 2+2/lambda);
  a3= pow( 1-x, 3+2/lambda);

	dfrac=(mu*a2*x*(2*x+lambda*(2-2*a1+x+2*a1*x)))/(lambda*(mu*a3+x*x*(1-a1+a1*x))*(mu*a3+x*x*(1-a1+a1*x)));


  return dfrac;
}


double BrooksCoreyFrontVelocity::fxbc(double x, double lambda, double mu)
{
  double nullst,a1,a2,a3;


  a1= pow( 1-x, 2/lambda );
  a2= pow( 1-x, 3+2/lambda);
  a3= pow( 1-x, 1+2/lambda);

	nullst=-x*(1-a3)/(mu*a2+x*x*(1-a3))+
	(mu*(1-x)*a3*x*(2*x+lambda*(2-2*a1+x+2*x*a1)))/
  	(lambda*(mu*a2+x*x*(1-a1+a1*x))*(mu*a2+x*x*(1-a1+a1*x)));


  return nullst;
}



double BrooksCoreyFrontVelocity::dfxbc(double x, double lambda, double mu)
{

  double abl,a1,a2,a3;

  a1= pow( 1-x, 2/lambda );
  a2= pow( 1-x, 3+2/lambda);
  a3= pow( 1-x, 1+2/lambda);
  abl=-(4*mu*a3*x*x*(-mu*a2+x*x*(1+a1-a1*x))+2*lambda*mu*x*a3*(-mu*a2*(3+4*x)-x*x*(-5+5*a1-9*a1*x+4*a1*x*x))-
lambda*lambda*(
  x*x*x*x*pow(1-a1+a1*x,3)+mu*mu*(1-x)*(1-x)*a3*a3*(1-a1+(7-3*a1)*x+(4+3*a1)*x*x+a1*x*x*x)+
  mu*a3*x*x*(
    -6*(-1+a1)*(-1+a1)
    +x*(-3-11*a1+14*a1*a1)
    +x*x*(3-8*a1-8*a1*a1)
    +x*x*x*(7*a1-2*a1*a1)
    +2*x*x*x*x*a1*a1)))/
(lambda*lambda*pow(mu*a2+x*x*(1-a1+a1*x),3));

    return abl;

}



void BrooksCoreyFrontVelocity::funcdbc(double x, double *fl, double *df, double lambda, double mu)

{
(*fl) = fxbc(x,lambda,mu);
(*df) = dfxbc(x,lambda,mu);
}


void BrooksCoreyFrontVelocity::zbrakbc(double x1,double x2,int n,double *xb1,double *xb2,double lambda, double mu)

{
  int i;
  double x,fp,fc,dx;

  dx=(x2-x1)/n;
  fp=fxbc(x=x1,lambda,mu);
  for (i=1;i<=n;i++) {
    fc=fxbc(x += dx,lambda,mu);
    if(fc*fp < 0.0){
      *xb1=x-dx;
      *xb2=x;
    }
    fp=fc;

  }
}



double BrooksCoreyFrontVelocity::findrootbc(double x1,double x2,double xacc,double lambda,double mu)
{
  int j;
  double df,dx,dxold,f,fh,fl;
  double temp,xh,xl,rts;

  funcdbc(x1,&fl,&df,lambda,mu);
  funcdbc(x2,&fh,&df,lambda,mu);
  if (fl*fh >= 0){
    cout <<"\nRoot must be bracketed.\n";
  }
  if (fl < 0.0) {
    xl=x1;
    xh=x2;
  } else {
    xh=x1;
    xl=x2;
  }

  rts=0.5*(x1+x2);
  dxold=fabs(x2-x1);
  dx=dxold;
  funcdbc(rts,&f,&df,lambda,mu);
  for (j=1;j<=MAXIT;j++) {
    if((((rts-xh)*df-f)*((rts-xl)*df-f)>=0.0)
       ||(fabs(2.0*f) > fabs(dxold*df))) {
      dxold=dx;
      dx=0.5*(xh-xl);
      rts=xl+dx;
      if (xl==rts) return rts;
    } else {
      dxold=dx;
      dx=f/df;
      temp=rts;
      rts-=dx;
      if (temp == rts) return rts; 
        }
    if(fabs(dx) < xacc) return rts;
    funcdbc(rts,&f,&df,lambda,mu);
    if(f<0.0)
      xl=rts;
    else
      xh=rts;
  }
  return rts;
}




double BrooksCoreyFrontVelocity::ShockVelocityMultiplier( double lambda, double visc_ratio )
{
     double x1(0.),x2(0.);
     double xacc(1.0e-6);

     zbrakbc(1.0E-6,1.0-1.0E-6,100,&x1,&x2,lambda,visc_ratio);

	 // cout << "Velocity multiplier = " << findrootbc(x1,x2,xacc,lambda,visc_ratio) << endl; 


     return dffbc( findrootbc(x1,x2,xacc,lambda,visc_ratio),lambda, visc_ratio );
}

} // end namespace csmp
