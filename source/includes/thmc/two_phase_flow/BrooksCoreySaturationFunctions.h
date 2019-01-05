#ifndef CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H
#define CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H

#include "CSMP_definitions.h"
#include "ArrayVariable.h"

namespace csmp {

template<size_t> class Element;
  

 
template<size_t dim, template<size_t> class USER>
class BrooksCoreySaturationFunctions {
  public:
        BrooksCoreySaturationFunctions();
  
        // pc, kri, and derivative methods that use a user supplied saturation value
        
        double64 EffectiveSaturation_at( const Element<dim>* const, double64 sw ) const;
  
        double64 pc_at( const Element<dim>* const, double64 sw ) const;
  
        double64 dpcds_at( const Element<dim>* const, double64 sw ) const;
  
        double64 krw_at( const Element<dim>* const, double64 sw ) const;
  
        double64 krn_at( const Element<dim>* const, double64 sw ) const;
  
        double64 dkrwds_at( const Element<dim>* const, double64 sw ) const;
  
        double64 dkrnds_at( const Element<dim>* const, double64 sw ) const;
  
        double64 dpcds_at_Numerical(  const Element<dim>* const, double64 sw, double64 h = 0.00001 ) const;
  
        double64 dkrwds_at_Numerical( const Element<dim>* const, double64 sw, double64 delta_s=0.001 ) const;

        double64 dkrnds_at_Numerical( const Element<dim>* const, double64 sw, double64 delta_s=0.001 ) const;


       // pc, kri, and derivative methods that use saturation values at the element barycentre (for FE mobility calculations etc.)
  
        double64 EffectiveSaturation( const Element<dim>* const ) const;
  
        double64 pc( const Element<dim>* const ) const;
  
        double64 dpcds( const Element<dim>* const ) const;
  
        double64 krw( const Element<dim>* const ) const;
  
        double64 krn( const Element<dim>* const ) const;
  
        double64 dkrwds( const Element<dim>* const ) const;
  
        double64 dkrnds( const Element<dim>* const ) const;
  
        /// maximum value of pc (Pa)
        double64 MaxCapillaryPressure() const { return max_pc_; }
  
        /// Numerical derivatives of first derivatives of relative permeability of water and CO2
        double64 dkrwds_Numerical( const Element<dim>* const, double64 delta_s=0.001 ) const ;

        double64 dkrnds_Numerical( const Element<dim>* const, double64 delta_s=0.001 ) const ;
  
        double64 dpcds_Numerical(  const Element<dim>* const, double64 h = 0.00001 ) const;
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

    const double64 max_derivative_ = 1.0e+8; ///< the absolute value of any derivative calculated herein must be less than this value
    const double64 max_pc_         = 5.0e+7; ///< the absolute value of any capillary pressure must not exceed the tensile strength of the rock
};

  
}  // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H */
