#ifndef CSMP_HETEROGENEITY_AND_RATE_AWARE_SATURATION_FUNCTIONS_H
#define CSMP_HETEROGENEITY_AND_RATE_AWARE_SATURATION_FUNCTIONS_H

#include "CSMP_definitions.h"
#include "ArrayVariable.h"
//#include "OtwayCRC3_RockTypes.h"
#include "OtwayCRC3_RockTypes_Version_2.h"

namespace csmp {

template<size_t> class Element;
  

 
template<size_t dim, template<size_t> class USER>
class HeterogeneityAndRateAwareSaturationFunctions {
  public:
        HeterogeneityAndRateAwareSaturationFunctions();
        
        virtual ~HeterogeneityAndRateAwareSaturationFunctions();
  
        // pc, kri, and derivative methods that use a user supplied saturation value
        
        double64 EffectiveSaturation_at( Element<dim>* const, double64 sw ) const;
  
        double64 pc_at( Element<dim>* const, double64 sw ) const;
        double64 pc_at_BrookCorey( Element<dim>* const, double64 sw ) const;
        double64 pc_at_VanGenuchten( Element<dim>* const, double64 sw ) const; 
  
        double64 dpcds_at( Element<dim>* const, double64 sw ) const;
        double64 dpcds_at_BC( Element<dim>* const, double64 sw ) const; 
        double64 dpcds_at_VG( Element<dim>* const, double64 sw ) const;
  
        double64 krw_at( Element<dim>* const, double64 sw ) const;
  
        double64 krn_at( Element<dim>* const, double64 sw ) const;
  
        double64 dkrwds_at( Element<dim>* const, double64 sw ) const;
  
        double64 dkrnds_at( Element<dim>* const, double64 sw ) const;
  
        double64 dpcds_at_Numerical(  Element<dim>* const, double64 sw, double64 h = 0.00001 ) const;
        double64 dpcds_at_Numerical_BC(  Element<dim>* const, double64 sw, double64 h = 0.00001 ) const; 
        double64 dpcds_at_Numerical_VG(  Element<dim>* const, double64 sw, double64 h = 0.00001 ) const; 
  
        double64 dkrwds_at_Numerical( Element<dim>* const, double64 sw, double64 delta_s=0.001 ) const;

        double64 dkrnds_at_Numerical( Element<dim>* const, double64 sw, double64 delta_s=0.001 ) const;


       // pc, kri, and derivative methods that use saturation values at the element barycentre (for FE mobility calculations etc.)
  
        double64 EffectiveSaturation( Element<dim>* const ) const;
  
        double64 pc( Element<dim>* const ) const;
        double64 pc_BrookCorey( Element<dim>* const ) const; 
        double64 pc_VanGenuchten( Element<dim>* const ) const; 
  
        double64 dpcds( Element<dim>* const ) const;
        double64 dpcds_BC( Element<dim>* const ) const; 
        double64 dpcds_VG( Element<dim>* const ) const; 
  
        double64 krw( Element<dim>* const ) const;
  
        double64 krn( Element<dim>* const ) const;
  
        double64 dkrwds( Element<dim>* const ) const;
  
        double64 dkrnds( Element<dim>* const ) const;
  
        /// maximum value of pc (Pa)
        double64 MaxCapillaryPressure() const { return max_pc_; }
  
        /// Numerical derivatives of first derivatives of relative permeability of water and CO2
        double64 dkrwds_Numerical( Element<dim>* const, double64 delta_s=0.001 ) const ;

        double64 dkrnds_Numerical( Element<dim>* const, double64 delta_s=0.001 ) const ;
  
        double64 dpcds_Numerical(  Element<dim>* const, double64 h = 0.00001 ) const;
        double64 dpcds_Numerical_BC(  Element<dim>* const, double64 h = 0.00001 ) const; 
        double64 dpcds_Numerical_VG(  Element<dim>* const, double64 h = 0.00001 ) const;

        /// New functions
        void InitilizeRockProperties( Element<dim>* const e ) const;
        double64 GetSwr( Element<dim>* const e ) const;
        double64 GetSnr( Element<dim>* const e ) const;
        double64 GetK( Element<dim>* const e ) const;     
        double64 GetKV( Element<dim>* const e ) const; 
        double64 GetPhi( Element<dim>* const e ) const;     

        double64 GetSwiPc( Element<dim>* const e ) const;
        double64 GetmVG( Element<dim>* const e ) const;
        double64 GetmLow( Element<dim>* const e ) const;
        double64 GetPd( Element<dim>* const e ) const;
        double64 GetPd_VG( Element<dim>* const e ) const;
        double64 GetPdLow( Element<dim>* const e ) const;
        double64 GetBcp( Element<dim>* const e ) const;
        double64 GetKrw( Element<dim>* const e, double64 Sw ) const;
        double64 GetKrn( Element<dim>* const e, double64 Sw ) const;
        double64 GetKrwParallelDrainage( Element<dim>* const e, double64 Sw ) const;
        double64 GetKrwParallelDrainage( Element<dim>* const e, double64 Sw, double64 vt_magnitude) const;
        double64 GetKrwCrossDrainage( Element<dim>* const e, double64 Sw ) const;
        double64 GetKrwCrossDrainage( Element<dim>* const e, double64 Sw, double64 vt_magnitude) const;
        double64 GetKrnParallelDrainage( Element<dim>* const e, double64 Sw ) const;
        double64 GetKrnParallelDrainage( Element<dim>* const e, double64 Sw, double64 vt_magnitude) const;
        double64 GetKrnCrossDrainage( Element<dim>* const e, double64 Sw ) const;
        double64 GetKrnCrossDrainage( Element<dim>* const e, double64 Sw, double64 vt_magnitude) const;
        double64 GetLYLow( Element<dim>* const e ) const;
        double64 GetLYHigh( Element<dim>* const e ) const;
        double64 GetKLow( Element<dim>* const e ) const;
        double64 GetKHigh( Element<dim>* const e ) const;

        int32 RockType( Element<dim>* const e ) const;
        bool IsComposite( Element<dim>* const e ) const;
        double64 Bcp( Element<dim>* const e ) const;
        double64 Bcp_low( Element<dim>* const e ) const;
        VectorVariable<dim> InitializeVelocity( Element<dim>* const e ) const;
        double64 krw_parallel( Element<dim>* const e, size_t direction ) const;
        double64 krw_crossflow( Element<dim>* const e ) const;
        double64 krw_parallel_at( Element<dim>* const e, double64 Sw, size_t direction ) const;
        double64 krw_crossflow_at( Element<dim>* const e, double64 Sw  ) const;
        double64 PermeabilityPerpendicularToLaminations( Element<dim>* const e ) const;
        double64 PermeabilityParallelToLaminations( Element<dim>* const e ) const;
        double64 PermeabilityInFlowDirection( Element<dim>* const e, const TensorVariable<dim>& KK) const;
        double64 K_reduction_in_flow_direction( Element<dim>* const e ) const;
        double64 krn_parallel( Element<dim>* const e, size_t direction ) const;
        double64 krn_crossflow( Element<dim>* const e ) const;
        double64 krn_parallel_at( Element<dim>* const e, double64 Sw, size_t direction ) const;
        double64 krn_crossflow_at( Element<dim>* const e, double64 Sw ) const;
        
        void WriteRelativePermeabilityTable (const char* filename, long RT, Element<dim>* const e);
        void OutputTestingResults(Element<dim>* const e);


  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

    const double64 max_derivative_ = 1.0e+8; ///< the absolute value of any derivative calculated herein must be less than this value
    const double64 max_pc_         = 5.0e+7; ///< the absolute value of any capillary pressure must not exceed the tensile strength of the rock
    
    OtwayRockTypes Otway_;
};

  
}  // end namespace csmp

#endif /* CSMP_HETEROGENEITY_AND_RATE_AWARE_SATURATION_FUNCTIONS_H */
