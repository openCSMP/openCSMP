
#include "ExperimentalSaturationFunctions.h"
#include "CO2H2O_FunctionsModule1.h"
#include "Element.h"

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
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds_Numerical( const Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krw_at( e, seff + h) - this->krw_at( e, seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at( e, seff) - this->krw_at( e, seff - h)) / h * dSedSw;
    
    return ( this->krw_at( e, seff + h) - this->krw_at( e, seff - h) ) / (2. * h) * dSedSw;
  }
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at( e, sw));
    
    if ( seff < 0.+h )
      return ( this->krw_at( e, seff + h) - this->krw_at( e, seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krw_at( e, seff) - this->krw_at( e, seff - h)) / h * dSedSw;
    
    return ( this->krw_at( e, seff + h) - this->krw_at( e, seff - h) ) / (2. * h) * dSedSw;
  }
  
  
  
  
  
  
  
  
  
  
  template<size_t dim, template<size_t> class USER>
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds_Numerical( const Element<dim>* const e, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation(e));
    
    if ( seff < 0.+h )
      return ( this->krn_at( e, seff + h) - this->krn_at( e, seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at( e, seff) - this->krn_at( e, seff - h)) / h * dSedSw;
    
    return ( this->krn_at( e, seff + h) - this->krn_at( e, seff - h) ) / (2. * h) * dSedSw;
  }
  
  
  
  
  

  
  
  
  template<size_t dim, template<size_t> class USER>
  double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds_at_Numerical( const Element<dim>* const e, double64 sw, double64 h ) const
  {
    const double64 dSedSw(1./(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2)));
    const double64 seff(EffectiveSaturation_at( e, sw));
    
    if ( seff < 0.+h )
      return ( this->krn_at( e, seff + h) - this->krn_at( e, seff) ) / h * dSedSw;
    if ( seff > 1.-h )
      return ( this->krn_at( e, seff) - this->krn_at( e, seff - h)) / h * dSedSw;
    
    return ( this->krn_at( e, seff + h) - this->krn_at( e, seff - h) ) / (2. * h) * dSedSw;
  }
  
  
  
  
  
// INLINE FUNCTIONS
  
template<size_t dim, template<size_t> class USER>
size_t ExperimentalSaturationFunctions<dim,USER>::RockType( const Element<dim>* const e ) const
 {
    const size_t rocktype = static_cast<size_t>(e->Read( User()->key_RRT ));
    assert( rocktype < 254 );
    return rocktype;
 }








template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation( const Element<dim>* const e ) const
{
  double64 sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );

  double64 seff = (sw - e->Read(User()->key_srH2O)) /
                  (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
  
  return std::min( std::max( seff, 0. ), 1. );
}
  
  
  



template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation_at( const Element<dim>* const e, double64 sw ) const
  {
    double64 seff = (sw - e->Read(User()->key_srH2O)) / (1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    
    return std::min( std::max( seff, 0. ), 1. );
  }
  


/**
    TODO: capillary pressure also exists outside of the effective saturation range
*/
template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::pc( const Element<dim>* const e )
{
    return pc_[ RockType(e) ].Value( EffectiveSaturation(e) );
}




template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::dpcds( const Element<dim>* const e) const
  {
    double64 srH2O  = e->Read(User()->key_srH2O);
    double64 srCO2  = e->Read(User()->key_srCO2);
    
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    return pc_[ RockType(e) ].Derivative( EffectiveSaturation(e) ) * seff_mult;
}
  


template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::krw( const Element<dim>* const e ) const
{
   return kr1_[ RockType(e) ].Value( EffectiveSaturation(e) );
}
  


template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::krw_at( const Element<dim>* const e , double64 S) const
{
   double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
   return kr1_[ RockType(e) ].Value( seff );
}


/**
   
   For calculating the oil relative permeability, we use notation which it has been inherated from the Skjaeveland et al. 2000
   The relative permeability has been calculated from the Brooks Corey Capillary pressure. These formula has been called as
   the Corey-Burdine realative permeablity equations... for further information see the page 65 in Skjaeveland et al. 2000.
   
*/
template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::krn( const Element<dim>* const e) const
  {
    return kr2_[ RockType(e) ].Value( EffectiveSaturation(e) );
  }

  
template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::krn_at( const Element<dim>* const e, double64 S) const
{
    double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    return kr2_[ RockType(e) ].Value( seff );
}
  


/**
 
 calculating the 1st derivative of water relative permeability
 
*/
template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds( const Element<dim>* const e ) const
  {
    const double64 srH2O  = e->Read(User()->key_srH2O);
    const double64 srCO2  = e->Read(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    return kr1_[ RockType(e) ].Derivative( EffectiveSaturation(e) ) * seff_mult;
}


template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::dkrwds_at( const Element<dim>* const e , double64 S) const
  {
    const double64 srH2O  = e->Read(User()->key_srH2O);
    const double64 srCO2  = e->Read(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));

    return kr1_[ RockType(e) ].Derivative( seff ) * seff_mult;
}




template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds( const Element<dim>* const e) const
  {
    const double64 srH2O  = e->Read(User()->key_srH2O);
    const double64 srCO2  = e->Read(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    return kr2_[ RockType(e) ].Derivative( EffectiveSaturation(e) ) * seff_mult;
}



/**
   
   calculating the 1st derivative of oil relative permeability for any water saturations
   
*/
template<size_t dim, template<size_t> class USER>
double64 ExperimentalSaturationFunctions<dim,USER>::dkrnds_at( const Element<dim>* const e, double64 S ) const
  {
    const double64 srH2O  = e->Read(User()->key_srH2O);
    const double64 srCO2  = e->Read(User()->key_srCO2);
    const double64 seff_mult( 1.0/ (1.0 - srH2O - srCO2 ) );
    
    const double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    
    return kr2_[ RockType(e) ].Derivative( seff ) * seff_mult;
    
}


  
  

  
template class ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>;
template class ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>;
template class ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>;
  
  
} // csmp
