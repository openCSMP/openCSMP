#ifndef CVFEM_PRESSURE_GRADIENT_VISITOR_H
#define CVFEM_PRESSURE_GRADIENT_VISITOR_H

/*   Changelog
     February 2015, Philipp Weis:
	 - initial port to CSMP++ and strong simplification as compared to the csp5-version.
*/

#include "Visitor.h"
#include "Model.h"
#include "compareFloats.h"

namespace csmp {

/// Calculating the pressure gradient within the CVFEM scheme (Weis et al., Geofluids, 2014).

template<uint32_t dim>
class CVFEM_PressureGradientVisitor : public Visitor<dim> {
  public:
    CVFEM_PressureGradientVisitor( Model<dim>& model, 
                                   const char* fluid_pressure,        // nodal scalar
                                   const char* lithostatic_pressure,  // nodal scalar
                                   const char* permeability,          // element variable
                                   const char* KgradP,                // element vector variable for
								                                      // fluid pressure gradient x permeability
																	  // for further use in finite volume calculations
																	  // of CVFEM scheme
                                   const char* gradP_scaling);        // this one is only needed for permeability visitor - edit!

    ~CVFEM_PressureGradientVisitor();

    virtual void Visit(Element<dim>* n);   
	virtual void Visit(Region<dim>* n);

  private:
    void ComputeGradient( Element<dim>& e );
    void ComputeGradient2( Element<dim>& e );

    bool gradient_scaling;

    csmp::Index  k_key_, pres_key_, KgradP_key_,
                 lith_pres_key_, scale_key_;

    ScalarVariable  k_, gradP_factor_;
    VectorVariable<dim>  KgradP_;
    std::vector<ScalarVariable> p_, lp_;

    DenseMatrix<DM_MIN> DERIV_;

  };
  /**
     @class CVFEM_PressureGradientVisitor

     @author Philipp Weis, ETH Zuerich
     @section contact Contact
     philipp.weis@erdw.ethz.ch

     @changes changes Latest Changes                                                                                  
  
     @section motivation Motivation
      For the fully-upwinded scheme including fluid viscosity and relative permeability,
	  the CVFEM scheme needs the varibale KgradP instead of velcotities to perform the finite volume calculations.

     @section usage Usage
      Used within the CVFEM scheme (Weis et al., Geofluids, 2014).

     @code
	 The visitor calculates the pressure gradient at the elements barycenter from the nodal variable "fluid pressure".
	 This vector variable is multiplied with the element permeability k.
	 k * grad( P )
          
     @endcode
     
     @section dependencies Dependencies
     
     @section issues Known issues
	 The variable gradp_scaling is not well documented and only needed in combination with the permeability vistor.
	 The bool gradient_scaling can't be changed from outside.
     
     @section testing Testing
     testing was done in the period before publication in 2014.

  */

 } // csmp

#endif
