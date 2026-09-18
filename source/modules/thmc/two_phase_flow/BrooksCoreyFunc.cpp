// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "BrooksCoreyFunc.h"

using namespace std;


namespace csmp {

BrooksCoreyFunc::BrooksCoreyFunc()
 {
 }


// destructor
BrooksCoreyFunc::~BrooksCoreyFunc()
 {
 }

/// fractional flow function in terms of nonwetting phase
double BrooksCoreyFunc::ffbc(double x, double lambda, double visc_ratio, double reswet, double resnonwet)
{
  double seff, frac, a1, a2;


  seff = (1-x-reswet)/(1-reswet-resnonwet);
  a1 = (1-seff)*(1-seff)*(1-pow(seff,(2+lambda)/lambda));
  a2 = pow(seff,(2+3*lambda)/lambda);

	frac=a1/(a1+visc_ratio*a2);


  return frac;
}

/// derivative of fractional flow function in terms of nonwetting phase
double BrooksCoreyFunc::ffbcd(double x, double lambda, double visc_ratio, double reswet, double resnonwet)
{
	double seff, a1, a2, a3, a4;
	double dfrac;
	
	seff = (1-x-reswet)/(1-reswet-resnonwet);
	a1 = pow(seff,(2+lambda)/lambda);
	a2 = pow(seff,(2+2*lambda)/lambda);
	a3 = pow(seff,(2+3*lambda)/lambda);
	a4 = ((1-seff)*(1-seff)*(1-a1)+visc_ratio*a3);
	
	dfrac = -visc_ratio/lambda*1/(1-reswet-resnonwet)*
	((1-seff)*a2*(2*(seff-1)+lambda*(seff-3+2*a1)))/(a4*a4);
	
	return dfrac;

}

/// gravity function in terms of nonwetting phase
double BrooksCoreyFunc::gfbc(double x, double lambda, double visc_ratio, double reswet, double resnonwet)
{
	double seff, a1, a2;
	double grav;
	
	seff = (1-x-reswet)/(1-reswet-resnonwet);
  	a1 = (1-seff)*(1-seff)*(1-pow(seff,(2+lambda)/lambda));
  	a2 = pow(seff,(2+3*lambda)/lambda);
  	
  	grav = visc_ratio*a1*a2/(a1+visc_ratio*a2);
  	
  	return grav;

}


/// derivative of gravity function in terms of nonwetting phase
double BrooksCoreyFunc::gfbcd(double x, double lambda, double visc_ratio, double reswet, double resnonwet)
{
	double seff, a1, a2, a3, a4, a5, a6, b1, b2;
	double dgrav;
	
	seff = (1-x-reswet)/(1-reswet-resnonwet);
	
	a1 = pow(seff,(2+lambda)/lambda);
	a2 = pow(seff,(2+2*lambda)/lambda);
	a3 = pow(seff,(2+3*lambda)/lambda);
	a4 = pow(seff,(2+4*lambda)/lambda);
	a5 = pow(seff,(4+3*lambda)/lambda);
	a6 = pow(seff,(4+5*lambda)/lambda);
	
	b1 = (1-seff)*a2*(2*(seff-1)*(-1+2*seff-seff*seff-4*a2+2*a3-a1*a1+2*a5+(visc_ratio-1)*a2*a2+2*a1)
	     +lambda*(3-9*seff+9*seff*seff-3*seff*seff*seff+18*a2-18*a3-2*(visc_ratio-3)*a4+3*a1*a1-9*a5-(visc_ratio-9)*a2*a2+3*(visc_ratio-1)*a6-6*a1));
	b2 = ((1-seff)*(1-seff)*(1-a1)+visc_ratio*a3);
	
	dgrav = -visc_ratio/lambda*1/(1-reswet-resnonwet)*b1/(b2*b2);
	
	return dgrav;
	
}

/// capillary pressure function, x is here x nonwetting
double BrooksCoreyFunc::pcbc(double x, double pcentry, double lambda, double reswet, double resnonwet)
{
	double seff;
	double pc;
	
	seff = (1-x-reswet)/(1-reswet-resnonwet);
  	pc = pcentry*pow(seff,-1/lambda);
  	
  	return pc;

}

/// capillary pressure function, x is here x wetting
double BrooksCoreyFunc::pcwbc(double x, double pcentry, double lambda, double reswet, double resnonwet)
{
	double seff;
	double pc;
	
	seff = (x-resnonwet)/(1-reswet-resnonwet);
  	pc = pcentry*pow(seff,-1/lambda);
  	
  	return pc;

}

/// nonlinear diffusion term for the two-phase flow problem, careful: term is negative, assuming the diffusive flux is counted positive
double BrooksCoreyFunc::diffbc(double x, double pcentry, double lambda,  double visc_ratio, double reswet, double resnonwet)
{
	double seff,a1,a2,a3;
	double diff;
	
	seff = (1-x-reswet)/(1-reswet-resnonwet);
	a1 = (1-seff)*(1-seff)*(1-pow(seff,(2+lambda)/lambda));
  	a2 = pow(seff,(2+3*lambda)/lambda);
  	
  	a3 = visc_ratio*a1*a2/(a1+visc_ratio*a2);
  	
  	diff = (-1)/(1-reswet-resnonwet)*a3*pcentry/lambda*pow(seff,-(1+lambda)/lambda);
  	
  	return diff;

}

/// relative permeability of nonwetting phase (for Richard's equation), x is here x nonwetting
double BrooksCoreyFunc::relp1bc(double x, double lambda, double reswet, double resnonwet)
{
  double seff, rpnwbc;


  seff = (1-x-reswet)/(1-reswet-resnonwet);
  rpnwbc = (1-seff)*(1-seff)*(1-pow(seff,(2+lambda)/lambda));

  return rpnwbc;
}

/// relative permeability of wetting phase (for Richard's equation), careful: x is here x wetting!!!
double BrooksCoreyFunc::relp2bc(double x, double lambda, double reswet, double resnonwet)
{
  double seff, rpwbc;


  seff = (x-resnonwet)/(1-reswet-resnonwet);
  rpwbc = pow(seff,(2+3*lambda)/lambda);

  return rpwbc;
}

/// relative permeability of nonwetting phase (for Richard's equation), careful: x is here x nonwetting
double BrooksCoreyFunc::relp1bcd(double x, double lambda, double reswet, double resnonwet)
{
  double seff, rpnwbcd;


  seff = (1-x-reswet)/(1-reswet-resnonwet);
  rpnwbcd = (1)/(1-reswet-resnonwet)*(2*(1-seff)*(1-pow(seff,(2+lambda)/lambda))+(1-seff)*(1-seff)*(2+lambda)/lambda*(pow(seff,(2)/lambda)));

  return rpnwbcd;
}

/// derivative of relative permeability of wetting phase (for Richard's equation), careful: x is here x wetting!!!
double BrooksCoreyFunc::relp2bcd(double x, double lambda, double reswet, double resnonwet)
{
  double seff, rpwbcd;


  seff = (x-resnonwet)/(1-reswet-resnonwet);
  rpwbcd =(1)/(1-reswet-resnonwet)*(2+3*lambda)/lambda*pow(seff,(2+2*lambda)/lambda);

  return rpwbcd;
}

/// nonlinear diffusion term for Richard's equation if flowing fluid is nonwetting,
/// @warning term is negative, assuming the diffusive flux is counted positive
double BrooksCoreyFunc::diffrichnwbc(double x, double pcentry, double lambda,  double reswet, double resnonwet)
{
	double seff,a1,a2;
	double diffnwr;
	
	seff = (1-x-reswet)/(1-reswet-resnonwet);
	a1 = (1-seff)*(1-seff)*(1-pow(seff,(2+lambda)/lambda));
  	a2 = 1/lambda*pcentry*pow(seff,-(1+lambda)/lambda);
 
  	diffnwr =(-1)/(1-reswet-resnonwet)*a1*a2;
  	
  	return diffnwr;

}

/// nonlinear diffusion term for Richard's equation if flowing fluid is wetting,
/// @warning term is negative, assuming the diffusive flux is counted positive
double BrooksCoreyFunc::diffrichwbc(double x, double pcentry, double lambda,  double reswet, double resnonwet)
{
	double seff,a1,a2;
	double diffwr;
	
	seff = (x-resnonwet)/(1-reswet-resnonwet);
	a1 =  pow(seff,(2+3*lambda)/lambda);
  	a2 = 1/lambda*pcentry*pow(seff,-(1+lambda)/lambda);
 
  	diffwr = (-1)/(1-reswet-resnonwet)*a1*a2;
  	
  	return diffwr;

}



double BrooksCoreyFunc::VelocityMultiplierF( double x, double pcentry, double lambda, double visc_ratio, double reswet, double resnonwet, int cho )

{
	if(cho==1){
     	return ffbc( x, lambda, visc_ratio, reswet, resnonwet );
    }else if(cho==2){
        return ffbcd( x, lambda, visc_ratio, reswet, resnonwet );
	}else if(cho==3){
	    return gfbc( x, lambda, visc_ratio, reswet, resnonwet );
	}else if(cho==4){
	    return gfbcd( x, lambda, visc_ratio, reswet, resnonwet );
	}else if(cho==5){
	    return pcbc( x, pcentry, lambda, reswet, resnonwet );
	}else if(cho==6){
	    return pcwbc(x, pcentry, lambda, reswet, resnonwet);
	}else if(cho==7){
	    return diffbc(x, pcentry, lambda,visc_ratio, reswet, resnonwet);
	}else if(cho==8){
	    return relp1bc(x, lambda, reswet, resnonwet);
	}else if(cho==9){
	    return relp2bc(x, lambda, reswet, resnonwet);
	}else if(cho==10){
	    return relp1bcd(x, lambda, reswet, resnonwet);
    }else if(cho==11){
        return relp2bcd(x, lambda, reswet, resnonwet);
	}else if(cho==12){
	    return diffrichnwbc(x, pcentry, lambda, reswet, resnonwet);
	}else if(cho==13){
	    return diffrichwbc(x, pcentry, lambda, reswet, resnonwet);
	}else if(cho==14){
	    return VelocityMultiplierF(x, pcentry, lambda, visc_ratio, reswet, resnonwet, cho );
	}else{
	    return 0.0;
	}

 
}

} // end namespace csmp
