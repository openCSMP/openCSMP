#ifndef CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H
#define CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H

#include "CSMP_definitions.h"
#include "ArrayVariable.h"

namespace csmp {

template<uint32_t> class Element;
  

 
template<uint32_t dim, template<uint32_t> class USER>
class BrooksCoreySaturationFunctions {
  public:
        BrooksCoreySaturationFunctions();
  
        // pc, kri, and derivative methods that use a user supplied saturation value
        
        double EffectiveSaturation_at( Element<dim>* const, double sw ) const;
  
        double pc_at( Element<dim>* const, double sw ) const;
  
        double dpcds_at( Element<dim>* const, double sw ) const;
  
        double krw_at( Element<dim>* const, double sw ) const;
  
        double krn_at( Element<dim>* const, double sw ) const;
  
        double dkrwds_at( Element<dim>* const, double sw ) const;
  
        double dkrnds_at( Element<dim>* const, double sw ) const;
  
        double dpcds_at_Numerical(  Element<dim>* const, double sw, double h = 0.00001 ) const;
  
        double dkrwds_at_Numerical( Element<dim>* const, double sw, double delta_s=0.001 ) const;

        double dkrnds_at_Numerical( Element<dim>* const, double sw, double delta_s=0.001 ) const;


       // pc, kri, and derivative methods that use saturation values at the element barycentre (for FE mobility calculations etc.)
  
        double EffectiveSaturation( Element<dim>* const ) const;
  
        double pc( Element<dim>* const ) const;
  
        double dpcds( Element<dim>* const ) const;
  
        double krw( Element<dim>* const ) const;
  
        double krn( Element<dim>* const ) const;
  
        double dkrwds( Element<dim>* const ) const;
  
        double dkrnds( Element<dim>* const ) const;
  
        /// maximum value of pc (Pa)
        double MaxCapillaryPressure() const { return max_pc_; }
  
        /// Numerical derivatives of first derivatives of relative permeability of water and CO2
        double dkrwds_Numerical( Element<dim>* const, double delta_s=0.001 ) const ;

        double dkrnds_Numerical( Element<dim>* const, double delta_s=0.001 ) const ;
  
        double dpcds_Numerical(  Element<dim>* const, double h = 0.00001 ) const;
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

    const double max_derivative_ = 1.0e+8; ///< the absolute value of any derivative calculated herein must be less than this value
    const double max_pc_         = 5.0e+7; ///< the absolute value of any capillary pressure must not exceed the tensile strength of the rock
};

  
}  // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H */
