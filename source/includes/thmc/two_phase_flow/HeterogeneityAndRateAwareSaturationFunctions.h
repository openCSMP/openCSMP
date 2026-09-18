// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_HETEROGENEITY_AND_RATE_AWARE_SATURATION_FUNCTIONS_H
#define CSMP_HETEROGENEITY_AND_RATE_AWARE_SATURATION_FUNCTIONS_H

#include "CSMP_definitions.h"
#include "ArrayVariable.h"
//#include "OtwayCRC3_RockTypes.h"
#include "OtwayCRC3_RockTypes_Version_2.h"

namespace csmp {

template<uint32_t> class Element;
  

 
template<uint32_t dim, template<uint32_t> class USER>
class HeterogeneityAndRateAwareSaturationFunctions {
  public:
        HeterogeneityAndRateAwareSaturationFunctions();
        
        virtual ~HeterogeneityAndRateAwareSaturationFunctions();
  
        // pc, kri, and derivative methods that use a user supplied saturation value
        
        double EffectiveSaturation_at( Element<dim>* const, double sw ) const;
  
        double pc_at( Element<dim>* const, double sw ) const;
        double pc_at_BrookCorey( Element<dim>* const, double sw ) const;
        double pc_at_VanGenuchten( Element<dim>* const, double sw ) const; 
  
        double dpcds_at( Element<dim>* const, double sw ) const;
        double dpcds_at_BC( Element<dim>* const, double sw ) const; 
        double dpcds_at_VG( Element<dim>* const, double sw ) const;
  
        double krw_at( Element<dim>* const, double sw ) const;
  
        double krn_at( Element<dim>* const, double sw ) const;
  
        double dkrwds_at( Element<dim>* const, double sw ) const;
  
        double dkrnds_at( Element<dim>* const, double sw ) const;
  
        double dpcds_at_Numerical(  Element<dim>* const, double sw, double h = 0.00001 ) const;
        double dpcds_at_Numerical_BC(  Element<dim>* const, double sw, double h = 0.00001 ) const; 
        double dpcds_at_Numerical_VG(  Element<dim>* const, double sw, double h = 0.00001 ) const; 
  
        double dkrwds_at_Numerical( Element<dim>* const, double sw, double delta_s=0.001 ) const;

        double dkrnds_at_Numerical( Element<dim>* const, double sw, double delta_s=0.001 ) const;


       // pc, kri, and derivative methods that use saturation values at the element barycentre (for FE mobility calculations etc.)
  
        double EffectiveSaturation( Element<dim>* const ) const;
  
        double pc( Element<dim>* const ) const;
        double pc_BrookCorey( Element<dim>* const ) const; 
        double pc_VanGenuchten( Element<dim>* const ) const; 
  
        double dpcds( Element<dim>* const ) const;
        double dpcds_BC( Element<dim>* const ) const; 
        double dpcds_VG( Element<dim>* const ) const; 
  
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
        double dpcds_Numerical_BC(  Element<dim>* const, double h = 0.00001 ) const; 
        double dpcds_Numerical_VG(  Element<dim>* const, double h = 0.00001 ) const;

        /// New functions
        void InitilizeRockProperties( Element<dim>* const e ) const;
        double GetSwr( Element<dim>* const e ) const;
        double GetSnr( Element<dim>* const e ) const;
        double GetK( Element<dim>* const e ) const;     
        double GetKV( Element<dim>* const e ) const; 
        double GetPhi( Element<dim>* const e ) const;     

        double GetSwiPc( Element<dim>* const e ) const;
        double GetmVG( Element<dim>* const e ) const;
        double GetmLow( Element<dim>* const e ) const;
        double GetPd( Element<dim>* const e ) const;
        double GetPd_VG( Element<dim>* const e ) const;
        double GetPdLow( Element<dim>* const e ) const;
        double GetBcp( Element<dim>* const e ) const;
        double GetKrw( Element<dim>* const e, double Sw ) const;
        double GetKrn( Element<dim>* const e, double Sw ) const;
        double GetKrwParallelDrainage( Element<dim>* const e, double Sw ) const;
        double GetKrwParallelDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const;
        double GetKrwCrossDrainage( Element<dim>* const e, double Sw ) const;
        double GetKrwCrossDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const;
        double GetKrnParallelDrainage( Element<dim>* const e, double Sw ) const;
        double GetKrnParallelDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const;
        double GetKrnCrossDrainage( Element<dim>* const e, double Sw ) const;
        double GetKrnCrossDrainage( Element<dim>* const e, double Sw, double vt_magnitude) const;
        double GetLYLow( Element<dim>* const e ) const;
        double GetLYHigh( Element<dim>* const e ) const;
        double GetKLow( Element<dim>* const e ) const;
        double GetKHigh( Element<dim>* const e ) const;

        int32_t RockType( Element<dim>* const e ) const;
        bool IsComposite( Element<dim>* const e ) const;
        double Bcp( Element<dim>* const e ) const;
        double Bcp_low( Element<dim>* const e ) const;
        VectorVariable<dim> InitializeVelocity( Element<dim>* const e ) const;
        double krw_parallel( Element<dim>* const e, int direction ) const;
        double krw_crossflow( Element<dim>* const e ) const;
        double krw_parallel_at( Element<dim>* const e, double Sw, int direction ) const;
        double krw_crossflow_at( Element<dim>* const e, double Sw  ) const;
        double PermeabilityPerpendicularToLaminations( Element<dim>* const e ) const;
        double PermeabilityParallelToLaminations( Element<dim>* const e ) const;
        double PermeabilityInFlowDirection( Element<dim>* const e, const TensorVariable<dim>& KK) const;
        double K_reduction_in_flow_direction( Element<dim>* const e ) const;
        double krn_parallel( Element<dim>* const e, int direction ) const;
        double krn_crossflow( Element<dim>* const e ) const;
        double krn_parallel_at( Element<dim>* const e, double Sw, int direction ) const;
        double krn_crossflow_at( Element<dim>* const e, double Sw ) const;
        
        void WriteRelativePermeabilityTable (const char* filename, long RT, Element<dim>* const e);
        void OutputTestingResults(Element<dim>* const e);


  private:
    USER<dim>* User() { return static_cast<USER<dim>*>(this); }
    USER<dim> const* User() const { return static_cast<const USER<dim>*>(this); }

    const double max_derivative_ = 1.0e+8; ///< the absolute value of any derivative calculated herein must be less than this value
    const double max_pc_         = 5.0e+7; ///< the absolute value of any capillary pressure must not exceed the tensile strength of the rock
    
    OtwayRockTypes Otway_;
};

  
}  // end namespace csmp

#endif /* CSMP_HETEROGENEITY_AND_RATE_AWARE_SATURATION_FUNCTIONS_H */
