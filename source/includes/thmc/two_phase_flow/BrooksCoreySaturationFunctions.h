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
        
        /// Pseudo functions that do nothing
        void InitialiseBrooksCoreyParameters( const Element<dim>* const );  
        void UpdateBrooksCoreyParameters( Element<dim>* e );   
    
        double64 EffectiveSaturation( const Element<dim>* const ) const;
  
        double64 EffectiveSaturation_at( const Element<dim>* const, double64 s1 ) const;
  
        double64 pc( const Element<dim>* const ) const;
  
        double64 dpcds( const Element<dim>* const ) const;
  
  
        double64 krw( const Element<dim>* const ) const ;
    
        double64 krw_at( const Element<dim>* const , double64 S) const ;
  
        double64 krn( const Element<dim>* const ) const ;
    
        double64 krn_at( const Element<dim>* const, double64 S) const ;
  
        double64 dkrwds( const Element<dim>* const ) const ;
    
        double64 dkrwds_at( const Element<dim>* const, double64 S) const ;
  
        double64 dkrnds( const Element<dim>* const ) const ;
    
        double64 dkrnds_at( const Element<dim>* const, double64 S) const ;
  
        /// maximum value of pc
        double64 MaxCapillaryPressure() const { return 1e7; /* Pa */ }
        /// maximum value of dpcdS
        double64 MaxCapillaryPressureDerivative() const { return 1e6; /* Pa m-1 */ }
  
  
        /// Numerical derivatives of first derivatives of relative permeability of water and CO2
        double64 dkrwds_Numerical( const Element<dim>* const, double64 delta_s=0.001 ) const ;

        double64 dkrwds_at_Numerical( const Element<dim>* const, double64 sw, double64 delta_s=0.001 ) const ;

        double64 dkrnds_Numerical( const Element<dim>* const, double64 delta_s=0.001 ) const ;
  
        double64 dkrnds_at_Numerical( const Element<dim>* const, double64 sw, double64 delta_s=0.001 ) const ;
        
        double64 dpcds_Numerical(  const Element<dim>* const, double64 h = 0.00001 ) const;

  
    
  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }
};

  
}  // end namespace csmp

#endif /* CSMP_BROOKS_COREY_SATURATION_FUNCTIONS_H */
