
#include "ExperimentalSaturationFunctions.h"
#include "FlowFunctions1.h"
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
  


template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
size_t ExperimentalSaturationFunctions<dim,USER>::RockType( const TARGET_PLACEMENT& p ) const
 {
    const size_t rocktype = static_cast<size_t>(p.Obtain( User()->key_RRT ));
    assert( rocktype < 254 );
    return rocktype;
 }

template size_t ExperimentalSaturationFunctions<1U,FlowFunctions2>::RockType( const FiniteElementPlacement<1U,ELEMENT>&  ) const;
template size_t ExperimentalSaturationFunctions<2U,FlowFunctions2>::RockType( const FiniteElementPlacement<2U,ELEMENT>&  ) const;
template size_t ExperimentalSaturationFunctions<3U,FlowFunctions2>::RockType( const FiniteElementPlacement<3U,ELEMENT>&  ) const;


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

 
  

  
  
/**
   
  This calculate the effective saturation..
   
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT>
double64 ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation( const TARGET_PLACEMENT& p ) const
{
  double64 seff = (p.Obtain(User()->key_sH2O) - p.Obtain(User()->key_srH2O)) /
                  (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
  
  return std::min( std::max( seff, 0. ), 1. );
}

  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::EffectiveSaturation( const FiniteElementPlacement<1U,ELEMENT>&  ) const ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::EffectiveSaturation( const FiniteElementPlacement<2U,ELEMENT>&  ) const ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteElementPlacement<3U,ELEMENT>&  ) const ;

template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<1U,NODE>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<2U,NODE>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteVolumePlacement<3U,NODE>& ) const;

template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteElementPlacement<3U,NODE>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;


  
  
  
  
  
  
  
  
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT>
  double64 ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation_at( const TARGET_PLACEMENT& p, double64 sw ) const
  {
    double64 seff =  (sw - p.Obtain(User()->key_srH2O)) /
    (1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
  }
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::EffectiveSaturation_at( const FiniteElementPlacement<1U,ELEMENT>&  , double64 sw) const ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::EffectiveSaturation_at( const FiniteElementPlacement<2U,ELEMENT>&  , double64 sw) const ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation_at( const FiniteElementPlacement<3U,ELEMENT>&  , double64 sw) const ;

template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::EffectiveSaturation_at( const FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::EffectiveSaturation_at( const FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation_at( const FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::EffectiveSaturation_at( const FiniteVolumePlacement<1U,NODE>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::EffectiveSaturation_at( const FiniteVolumePlacement<2U,NODE>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::EffectiveSaturation_at( const FiniteVolumePlacement<3U,NODE>&, double64 sw ) const;



  


  
  
  
  
  
  
  
/**
    TODO: capillary pressure also exists outside of the effective saturation range
*/
template<size_t dim, template<size_t> class USER>
template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::pc( const TARGET_PLACEMENT& p )
{
    return pc_[ RockType(p) ].Value( EffectiveSaturation(p) );
}

template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::pc( const FiniteElementPlacement<1U,ELEMENT>&  )  ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::pc( const FiniteElementPlacement<2U,ELEMENT>&  )  ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::pc( const FiniteElementPlacement<3U,ELEMENT>&  )  ;

template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::pc( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& );
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::pc( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& );
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::pc( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& );  

  
  
  
  
  
 
  
  
  
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::dpcds( const TARGET_PLACEMENT& p) const
  {
    
    double64 srH2O  = p.Obtain(User()->key_srH2O);
    double64 srCO2  = p.Obtain(User()->key_srCO2);
    
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    return pc_[ RockType(p) ].Derivative( EffectiveSaturation(p) ) * seff_mult;
}
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dpcds( const FiniteElementPlacement<1U,ELEMENT>&  ) const ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dpcds( const FiniteElementPlacement<2U,ELEMENT>&  ) const ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dpcds( const FiniteElementPlacement<3U,ELEMENT>&  ) const ;
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dpcds( const FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dpcds( const FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dpcds( const FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;

template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dpcds( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dpcds( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dpcds( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  


  
  
  
  
 
  
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::krw( const TARGET_PLACEMENT& p ) const
  {
     return kr1_[ RockType(p) ].Value( EffectiveSaturation(p) );
}
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krw( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krw( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw( const FiniteElementPlacement<3U,ELEMENT>& ) const ;
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krw( const FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krw( const FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw( const FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krw( const FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krw( const FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw( const FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krw( const FiniteVolumePlacement<1U,NODE>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krw( const FiniteVolumePlacement<2U,NODE>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw( const FiniteVolumePlacement<3U,NODE>& ) const;

template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  

  
  
  
  
  
  
  
  
  
  
  
/**
   
   Second form of the water relative permeability.. the same as previous function for any saturation..
   
*/
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::krw_at( const TARGET_PLACEMENT& p , double64 S) const
  {
    
      double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
      return kr1_[ RockType(p) ].Value( seff );
}
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krw_at( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krw_at( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw_at( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krw_at( const FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krw_at( const FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw_at( const FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krw_at( const FiniteVolumePlacement<1U,NODE>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krw_at( const FiniteVolumePlacement<2U,NODE>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krw_at( const FiniteVolumePlacement<3U,NODE>&, double64 ) const;
  

  
  
  
  
  
  
  
/**
   
   For calculating the oil relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
   The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
   the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
   
*/
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::krn( const TARGET_PLACEMENT& p) const
  {
    
    return kr2_[ RockType(p) ].Value( EffectiveSaturation(p) );
}

template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krn( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krn( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn( const FiniteElementPlacement<3U,ELEMENT>& ) const ;
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krn( const FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krn( const FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn( const FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krn( const FiniteVolumePlacement<1U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krn( const FiniteVolumePlacement<2U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn( const FiniteVolumePlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krn( const FiniteVolumePlacement<1U,NODE>& ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krn( const FiniteVolumePlacement<2U,NODE>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn( const FiniteVolumePlacement<3U,NODE>& ) const;

template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  

  
/**
   
   Second form of the oil relative permeability.. the same as previous function for any saturation..
   
*/
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::krn_at( const TARGET_PLACEMENT& p, double64 S) const
  {
    
    double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    return kr2_[ RockType(p) ].Value( seff );
}
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krn_at( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krn_at( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn_at( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krn_at( const FiniteVolumePlacement<1U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krn_at( const FiniteVolumePlacement<2U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn_at( const FiniteVolumePlacement<3U,FACET_INTEGRATION_POINT>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::krn_at( const FiniteVolumePlacement<1U,NODE>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::krn_at( const FiniteVolumePlacement<2U,NODE>&, double64 sw ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::krn_at( const FiniteVolumePlacement<3U,NODE>&, double64 sw ) const;
  

  
  
/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds( const TARGET_PLACEMENT& p ) const
  {
    
    double64 srH2O  = p.Obtain(User()->key_srH2O);
    double64 srCO2  = p.Obtain(User()->key_srCO2);
    
    double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    return kr1_[ RockType(p) ].Derivative( EffectiveSaturation(p) )*seff_mult;
}
  
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrwds( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrwds( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds( const FiniteElementPlacement<3U,ELEMENT>& ) const ;
  
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  

  
  
  
  
  
  
  
  
/**
   
   calculating the 1st derivative of water relative permeability for any water saturation
   
*/
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds_at( const TARGET_PLACEMENT& p , double64 S) const
  {
    
    double64 srH2O  = p.Obtain(User()->key_srH2O);
    double64 srCO2  = p.Obtain(User()->key_srCO2);
    
    double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));

    return kr1_[ RockType(p) ].Derivative( seff )*seff_mult;
}
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrwds_at( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrwds_at( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrwds_at( const FiniteVolumePlacement<1U,NODE>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrwds_at( const FiniteVolumePlacement<2U,NODE>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at( const FiniteVolumePlacement<3U,NODE>&, double64 ) const;


  
  
  
  
  
  
  
  
  
  
  
  
/**
   
  calculating the 1st derivative of oil relative permeability
   
*/
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds( const TARGET_PLACEMENT& p) const
  {
    
    double64 srH2O  = p.Obtain(User()->key_srH2O);
    double64 srCO2  = p.Obtain(User()->key_srCO2);
    
    double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    return kr2_[ RockType(p) ].Derivative( EffectiveSaturation(p) )*seff_mult;
}
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrnds( const FiniteElementPlacement<1U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrnds( const FiniteElementPlacement<2U,ELEMENT>& ) const ;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds( const FiniteElementPlacement<3U,ELEMENT>& ) const ;
  
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>& ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>& ) const;
  

  
  
  
  
  
  
  
/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<size_t dim, template<size_t> class USER>
  template<class TARGET_PLACEMENT> double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds_at( const TARGET_PLACEMENT& p, double64 S ) const
  {
    
    double64 srH2O  = p.Obtain(User()->key_srH2O);
    double64 srCO2  = p.Obtain(User()->key_srCO2);
    
    double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    double64 seff =  (S - p.Obtain(User()->key_srH2O)) /(1. - p.Obtain(User()->key_srH2O) - p.Obtain(User()->key_srCO2));
    
    return kr2_[ RockType(p) ].Derivative( seff )*seff_mult;
    
}
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrnds_at( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrnds_at( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;

template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrnds_at( const FiniteVolumePlacement<1U,NODE>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrnds_at( const FiniteVolumePlacement<2U,NODE>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at( const FiniteVolumePlacement<3U,NODE>&, double64 ) const;


  
  
  
  
  
  
  
  
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
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrwds_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrwds_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_Numerical( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_Numerical( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;  
  
  
  
  
  
  
  
  
  
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
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64 , double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64 , double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64, double64  ) const;
  
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteVolumePlacement<1U,NODE>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteVolumePlacement<2U,NODE>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrwds_at_Numerical( const FiniteVolumePlacement<3U,NODE>&, double64, double64 ) const;    
  
  
  
  
  
  
  
  
  
  
  
  
  
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
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrnds_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrnds_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64 ) const;
  
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_Numerical( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_Numerical( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64 ) const;  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
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
  
  
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteElementPlacement<1U,ELEMENT>&, double64 , double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteElementPlacement<2U,ELEMENT>&, double64 , double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteElementPlacement<3U,ELEMENT>&, double64, double64  ) const;
  
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteElementPlacement<3U,ELEMENT_INTEGRATION_POINT>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteElementPlacement<3U,FACET_INTEGRATION_POINT>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteElementPlacement<3U,SECTOR_INTEGRATION_POINT>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<1U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteVolumePlacement<1U,NODE>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<2U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteVolumePlacement<2U,NODE>&, double64, double64 ) const;
template double64 ExperimentalSaturationFunctions<3U,FlowFunctions2>::dkrnds_at_Numerical( const FiniteVolumePlacement<3U,NODE>&, double64, double64 ) const;    
  
  
  
  
  
  
  
  
  

  
template class ExperimentalSaturationFunctions<1U,FlowFunctions2>;
template class ExperimentalSaturationFunctions<2U,FlowFunctions2>;
template class ExperimentalSaturationFunctions<3U,FlowFunctions2>;
  
  
} // csmp
