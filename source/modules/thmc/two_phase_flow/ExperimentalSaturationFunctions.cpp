
#include "ExperimentalSaturationFunctions.h"
#include "CO2H2O_FunctionsModule1.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"

using namespace std;

namespace csmp {
  
/**
   The default constructor of Exoperimental Saturation Functions With Hysteresis class
*/
template<size_t dim, template<size_t> class USER>
ExperimentalSaturationFunctions<dim,USER>::ExperimentalSaturationFunctions( const char* filename )
  {
    ConstructRTs(filename);
  }
  



/**
 
 read from Rock Type file
 
 */
 template<size_t dim, template<size_t> class USER>
 void ExperimentalSaturationFunctions<dim,USER>::ConstructRTs(const char* rt_file_name)
  {
    std::cout << "\nExperimentalRT::ConstructRTs()\n";
    
    std::ifstream rt_file;
    rt_file.open(rt_file_name);
    
    if (!rt_file)
    {
      std::cout << "\nExperimentalRT::ConstructRTs(): " << rt_file_name << " file doesn't exist.\n";
      std::exit(1);
    }
    else
    {
      unsigned int number_of_tables;
      rt_file >> number_of_tables;
      
      std::cout << "\nExperimentalRT::ConstructRTs(): Number of tables: " << number_of_tables << "\n";
      
      kr1_.resize(number_of_tables);
      kr2_.resize(number_of_tables);
      pc_.resize(number_of_tables);
      
      for (unsigned int i = 0; i < number_of_tables; i++)
        {
          std::cout << "\nTable " << i << "\n";
          
          double64 kro_start_derivative, kro_end_derivative, krw_start_derivative, krw_end_derivative, pc_start_derivative, pc_end_derivative;
          std::vector<double64> sw, kro, krw, pc;
          
          double64 sw_value, kro_value, krw_value, pc_value;
          
          unsigned int number_of_entries;
          
          rt_file >> number_of_entries;
          
          rt_file >> kro_start_derivative >> kro_end_derivative >> krw_start_derivative >> krw_end_derivative >> pc_start_derivative >> pc_end_derivative;
          
          std::cout << "Derivaties:\t" << kro_start_derivative << "\t" << kro_end_derivative << "\t" << krw_start_derivative << "\t" << krw_end_derivative
          << "\t" << pc_start_derivative << "\t" << pc_end_derivative << "\n";
          
          std::cout << "sw\tkro\tkrw\tpc\n";
          
          for ( size_t n = 0; n < number_of_entries; n++ )
          {
            rt_file >> sw_value >> kro_value >> krw_value >> pc_value;
            
            std::cout << sw_value << "\t" << kro_value << "\t" << krw_value << "\t" << pc_value <<"\n";
            
            sw.push_back(sw_value);
            kro.push_back(kro_value);
            krw.push_back(krw_value);
            pc.push_back(pc_value);
          }
          
          std::cout << "\n"; 
          
          kr1_[i].Initialize(sw, krw, krw_start_derivative, krw_end_derivative);
          kr2_[i].Initialize(sw, kro, kro_start_derivative, kro_end_derivative);
          pc_[i].Initialize(sw, pc, pc_start_derivative, pc_end_derivative);
       }
    }
    
  } // end ConstructRTs

 
  

  
  
  
  
  // Numerical derivative added
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds_Numerical( const TARGET_PLACEMENT& p, double64 h ) const
  {
    const double64 dSedSw(1./(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(p));
    
    if ( seff < 0.+h )
      return ( this->krw_at(p,seff + h) - this->krw_at(p,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(p,seff) - this->krw_at(p,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(p,seff + h) - this->krw_at(p,seff - h) ) / (2. * h) * dSedSw;
  }
  
  template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrwds_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrwds_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrwds_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
  template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrwds_Numerical( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrwds_Numerical( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrwds_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds_at_Numerical( const TARGET_PLACEMENT& p, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(p,sw));
    
    if ( seff < 0.+h )
      return ( this->krw_at(p,seff + h) - this->krw_at(p,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at(p,seff) - this->krw_at(p,seff - h)) / h * dSedSw;
    
    return ( this->krw_at(p,seff + h) - this->krw_at(p,seff - h) ) / (2. * h) * dSedSw;
  }
  
  template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrwds_at_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64 , double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrwds_at_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64 , double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrwds_at_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64, double64  ) const;
  
  template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrwds_at_Numerical( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, double64 , double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrwds_at_Numerical( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, double64 , double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrwds_at_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64, double64  ) const;
  
  
  
  
  
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds_Numerical( const TARGET_PLACEMENT& p, double64 h ) const
  {
    const double64 dSedSw(1./(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(p));
    
    if ( seff < 0.+h )
      return ( this->krn_at(p,seff + h) - this->krn_at(p,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(p,seff) - this->krn_at(p,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(p,seff + h) - this->krn_at(p,seff - h) ) / (2. * h) * dSedSw;
  }
  
  template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrnds_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrnds_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrnds_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
   template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrnds_Numerical( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrnds_Numerical( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrnds_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds_at_Numerical( const TARGET_PLACEMENT& p, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at(p,sw));
    
    if ( seff < 0.+h )
      return ( this->krn_at(p,seff + h) - this->krn_at(p,seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at(p,seff) - this->krn_at(p,seff - h)) / h * dSedSw;
    
    return ( this->krn_at(p,seff + h) - this->krn_at(p,seff - h) ) / (2. * h) * dSedSw;
  }
  
  
  template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrnds_at_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64, double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrnds_at_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64, double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrnds_at_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64, double64  ) const;
  
   template double64 ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>::dkrnds_at_Numerical( const FiniteElementPlacement<1U,FACET_INTEGRATION_POINT>&, double64, double64 ) const;
  template double64 ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>::dkrnds_at_Numerical( const FiniteElementPlacement<2U,FACET_INTEGRATION_POINT>&, double64, double64 ) const;
  template double64 ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>::dkrnds_at_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64, double64  ) const;
 
  
  
  
  
  
  
  
  

  
template class ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>;
template class ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>;
template class ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>;
  
  
} // csmp
