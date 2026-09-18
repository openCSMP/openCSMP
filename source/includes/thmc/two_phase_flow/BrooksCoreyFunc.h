// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef BROOKS_COREY_FUNC_H
#define BROOKS_COREY_FUNC_H 

#include "CSMP_definitions.h"

namespace csmp {

/// @author Insa Neuweiler @date 2002
class BrooksCoreyFunc {
  public:
    BrooksCoreyFunc();
    ~BrooksCoreyFunc();
 	double VelocityMultiplierF( double x, double pcentry, double lambda, 
 	                            double visc_ratio, double reswet, double resnonwet, int cho);   

  private:
	double ffbc(double x, double lambda, double visc_ratio, double reswet, double resnonwet);
	double ffbcd(double x, double lambda, double visc_ratio, double reswet, double resnonwet);
	double gfbc(double x, double lambda, double visc_ratio, double reswet, double resnonwet);
	double gfbcd(double x, double lambda, double visc_ratio, double reswet, double resnonwet);
	double pcbc(double x, double pcentry, double lambda, double reswet, double resnonwet);
	double pcwbc(double x, double pcentry, double lambda, double reswet, double resnonwet);
	double diffbc(double x, double pcentry, double lambda,  double mu, double reswet, double resnonwet);
	double relp1bc(double x, double lambda, double reswet, double resnonwet);
	double relp2bc(double x, double lambda, double reswet, double resnonwet);
	double relp1bcd(double x, double lambda, double reswet, double resnonwet);
	double relp2bcd(double x, double lambda, double reswet, double resnonwet);
	double diffrichnwbc(double x, double pcentry, double lambda, double reswet, double resnonwet);
	double diffrichwbc(double x, double pcentry, double lambda, double reswet, double resnonwet);

};

} // end namespace csmp

#endif
