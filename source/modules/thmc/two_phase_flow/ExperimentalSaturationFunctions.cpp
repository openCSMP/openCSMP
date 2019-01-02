
#include "ExperimentalSaturationFunctions.h"
#include "CO2H2O_FunctionsModule1.h"
#include "Element.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {
  
/**
   The default constructor of Exoperimental Saturation Functions With Hysteresis class
*/
template<size_t dim, template<size_t> class USER>
ExperimentalSaturationFunctions<dim,USER>::ExperimentalSaturationFunctions( const char* filename )
  {
     InitialiseReservoirRockTypes(filename);
  }
  


/**
 
 reading the rocktype file
 
 */
template<size_t dim, template<size_t> class USER>
size_t ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes( const char* rt_file_name )
  {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
    std::cout << "\nExperimentalSaturationFunctions::InitialiseReservoirRockTypes: reading data from file: '"<< rt_file_name <<"' ...\n";
    
    std::ifstream rt_file;
    rt_file.open(rt_file_name);
    
    if ( !rt_file.is_open() )
      throw csmp::Exception( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes", "input file not found." );

    int number_of_tables(0);
    rt_file >> number_of_tables;
  
    std::cout << "\nExperimentalSaturationFunctions::InitialiseReservoirRockTypes: Number of tables: " << number_of_tables << "\n";
    if ( number_of_tables <= 0 )
      throw csmp::Exception( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes", "no data tables contained in input file." );

    kr1_.resize(number_of_tables);
    kr2_.resize(number_of_tables);
    pc_.resize(number_of_tables);
  
    for ( int i=0; i < number_of_tables; i++ )
      {
        std::cout << "\nReading data table for RRT: " << i << "\n";
      
        unsigned int number_of_entries;
        rt_file >> number_of_entries;
      
        double64 kro_start_derivative, kro_end_derivative, krw_start_derivative, krw_end_derivative, pc_start_derivative, pc_end_derivative;
        rt_file >> kro_start_derivative >> kro_end_derivative >> krw_start_derivative >> krw_end_derivative >> pc_start_derivative >> pc_end_derivative;
        cout << "Derivatives (krnw0,krnw1,krw0,krw1,dpcds0,dpcds1):\t" << kro_start_derivative << "\t" << kro_end_derivative;
        cout << "\t" << krw_start_derivative << "\t" << krw_end_derivative << "\t" << pc_start_derivative << "\t" << pc_end_derivative << "\n";
 
        if ( fabs(pc_start_derivative) > max_derivative_ ) {
             cerr <<"\n\tsaturation derivative of pc at sw_min: "<< pc_start_derivative;
             csmp_error.notice( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
          }
        if ( fabs(pc_end_derivative) > max_derivative_ ) {
             cerr <<"\n\tsaturation derivative of pc at sw_max: "<< pc_end_derivative;
             csmp_error.notice( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
          }
        if ( fabs(krw_start_derivative) > max_derivative_ ) {
             cerr <<"\n\tsaturation derivative of krw at sw_min: "<< krw_start_derivative;
             csmp_error.notice( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
          }
        if ( fabs(krw_end_derivative) > max_derivative_ ) {
             cerr <<"\n\tsaturation derivative of krw at sw_max: "<< krw_end_derivative;
             csmp_error.notice( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
          }
        if ( fabs(kro_start_derivative) > max_derivative_ ) {
             cerr <<"\n\tsaturation derivative of krn at sw_min: "<< kro_start_derivative;
             csmp_error.notice( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
          }
        if ( fabs(kro_end_derivative) > max_derivative_ ) {
             cerr <<"\n\tsaturation derivative of krn at sw_max: "<< kro_end_derivative;
             csmp_error.notice( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
          }
      
        std::cout << "sw\tkro\tkrw\tpc\n";
      
        std::vector<double64> sw, kro, krw, pc;
        for ( size_t n = 0; n < number_of_entries; n++ )
          {
            double64 sw_value, kro_value, krw_value, pc_value;
            rt_file >> sw_value >> kro_value >> krw_value >> pc_value;
          
            std::cout << sw_value << "\t" << kro_value << "\t" << krw_value << "\t" << pc_value <<"\n";
          
            sw.push_back(sw_value);
            kro.push_back(kro_value);
            krw.push_back(krw_value);
            pc.push_back(pc_value);
          }
      
        std::cout << "\n";
      
        kr1_[i].Initialize( sw, krw, krw_start_derivative, krw_end_derivative );
        kr2_[i].Initialize( sw, kro, kro_start_derivative, kro_end_derivative );
        pc_[i].Initialize( sw, pc, pc_start_derivative, pc_end_derivative );
     }

  //Out();
  // limiting the rocktype number range in the property database to the actual maximum value

  cout << "\nExperimentalSaturationFunctions::InitialiseReservoirRockTypes: rock types initialised successfully from file '"<< rt_file_name <<".'"<< endl;
  return pc_.size();
  
} // end InitialiseReservoirRockTypes

 



template<size_t dim, template<size_t> class USER>
size_t ExperimentalSaturationFunctions<dim,USER>::RockTypes() const
 { return pc_.size(); }


  
  
  
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
   //double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
   double64 seff =  EffectiveSaturation_at(e, S);
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
    //double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    double64 seff =  EffectiveSaturation_at(e, S);
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
    
    //double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    double64 seff =  EffectiveSaturation_at(e, S);

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
    
    //const double64 seff =  (S - e->Read(User()->key_srH2O)) /(1. - e->Read(User()->key_srH2O) - e->Read(User()->key_srCO2));
    double64 seff =  EffectiveSaturation_at(e, S);
    
    return kr2_[ RockType(e) ].Derivative( seff ) * seff_mult;
    
}


  
template<size_t dim, template<size_t> class USER>
void ExperimentalSaturationFunctions<dim,USER>::Out() const
 {
    cout <<"\nExperimentalSaturationFunctions<dim,USER>::Out: current rock type data:\n";
    cout <<"\nrelative permeability of phase 1 (krw(sw)), for all rock types:\n";
    int number(0);
    for ( auto it=kr1_.begin(); it!=kr1_.end(); it++ ) {
        cout <<"\nrock type "<< number++;
        (*it).Out();
      }
    cout <<"\nrelative permeability of phase 2 (krnw(sw)), for all rock types:\n";
    number = 0;
    for ( auto it=kr2_.begin(); it!=kr2_.end(); it++ ) {
        cout <<"\nrock type "<< number++;
        (*it).Out();
      }
    cout <<"\ncapillary pressure curve, pc(sw), for all rock types:\n";
    number = 0;
    for ( auto it=kr2_.begin(); it!=kr2_.end(); it++ ) {
        cout <<"\nrock type "<< number++;
        (*it).Out();
      }

    cout <<"\nmaximum slope of derivative curves: "<< max_derivative_ << endl;
 }

  
template class ExperimentalSaturationFunctions<1U,CO2H2O_FunctionsModule2>;
template class ExperimentalSaturationFunctions<2U,CO2H2O_FunctionsModule2>;
template class ExperimentalSaturationFunctions<3U,CO2H2O_FunctionsModule2>;
  
  
} // csmp
