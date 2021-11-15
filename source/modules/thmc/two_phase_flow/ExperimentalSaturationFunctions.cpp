
#include "ExperimentalSaturationFunctions.h"
#include "FlowFunctionsModule.h"
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
        
        double kro_start_derivative, kro_end_derivative, krw_start_derivative, krw_end_derivative, pc_start_derivative, pc_end_derivative;
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
      
        std::vector<double> sw, kro, krw, pc;
        for ( size_t n = 0; n < number_of_entries; n++ )
        {
            double sw_value, kro_value, krw_value, pc_value;
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

 


/**
    returns for how many rocktypes information is stored in the underlying data arrays.
*/
template<size_t dim, template<size_t> class USER>
size_t ExperimentalSaturationFunctions<dim,USER>::RockTypes() const
 { return pc_.size(); }


  
  
/**
    returns the rocktype associated with the current Element
*/
template<size_t dim, template<size_t> class USER>
size_t ExperimentalSaturationFunctions<dim,USER>::RockType( Element<dim>* const e ) const
 {
    const size_t rocktype = static_cast<size_t>(e->Read( User()->key_RRT ));
    assert( rocktype < pc_.size() );
    return rocktype;
 }







/**
    (sw-swr) / (1-swr-snr); use only in 2-phase flow simulations
 
    TODO: this needs to be computed from the curves in the input file
*/
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation( Element<dim>* const e ) const
 {
    if ( e->Read(User()->key_srH2O) > 0. || e->Read(User()->key_srCO2) > 0. )
      throw csmp::Exception( ERROR, "ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation",
                             "Do not set end-point saturation values in conjunction with experimental flow functions as they are defined implicitly");
   
    return e->PropertyValueAtBaryCenter( User()->key_sH2O );
 }
  
  
  
 
    /// (sw-swr) / (1-swr-snr); use only in 2-phase flow simulations
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation_at( Element<dim>* const e, double sw ) const
 {
    if ( e->Read(User()->key_srH2O) > 0. || e->Read(User()->key_srCO2) > 0. )
      throw csmp::Exception( ERROR, "ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation_at",
                             "Do not set end-point saturation values in conjunction with experimental flow functions as they are defined implicitly");
    return sw;
 }







/**
    @note Note that capillary pressure also exists outside of the effective saturation range
*/
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::pc( Element<dim>* const e ) const
 {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return min( pc_[ RockType(e) ].Value(sw), max_capillary_pressure_ );
 }





template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::pc_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return min( pc_[ RockType(e) ].Value(sw), max_capillary_pressure_ );
 }



/**
    @attention this treatment assumes that the capillary pressure derivative with regard to the
    wetting phase saturation is always negative.
*/
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dpcds( Element<dim>* const e ) const
  {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( pc_[ RockType(e) ].Derivative(sw), -max_derivative_ );
}
  



template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dpcds_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( pc_[ RockType(e) ].Derivative(sw), -max_derivative_ );
}





template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krw( Element<dim>* const e ) const
{
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( kr1_[ RockType(e) ].Value( sw ), 0. );
}
  


template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krw_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );
   
    // as defined by spline curves over full saturation range
    return max( kr1_[ RockType(e) ].Value( sw ), 0. );
 }



  
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krn( Element<dim>* const e ) const
  {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( kr2_[ RockType(e) ].Value( sw ), 0. );
  }

 
 
 
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krn_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( kr2_[ RockType(e) ].Value( sw ), 0. );
 }
  


/**
 
    calculating the 1st derivative of water relative permeability at the element barycentre.
 
*/
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrwds( Element<dim>* const e ) const
  {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );
  
    return min( kr1_[ RockType(e) ].Derivative( sw ), max_derivative_ );
}





template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrwds_at( Element<dim>* const e, double sw ) const
{
    assert( sw >= 0. );
    assert( sw <= 1. );

    return min( kr1_[ RockType(e) ].Derivative( sw ), max_derivative_ );
 }




/**
    @note by contrast with the wetting phase relative permeability, the non-wetting phase
    kri is expected to have a negative slope w.r.t. water saturation.
*/
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrnds( Element<dim>* const e ) const
 {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
 
    return max( kr2_[ RockType(e) ].Derivative( sw ), -max_derivative_ );
 }




/**
    1st derivative of non-wetting phase relative permeability as a function of saturation.
*/
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrnds_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );
  
    return max( kr2_[ RockType(e) ].Derivative( sw ), -max_derivative_ );
 }



// NUMERICAL DERIVATIVE CALCULATIONS


template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dpcdsw_Numerical( Element<dim>* const e, double h ) const
  {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );

    // deal with the more common case of a high water saturation first
    if ( sw > (1. - h) )
      return (pc_at( e, 1. ) - pc_at( e, 1. - h)) / h;
  
    // low water saturation
    if ( sw < h )
      return ( pc_at( e, h ) - pc_at( e, 0. ) ) / h;

    // water saturation between the endpoints
    return ( pc_at( e, sw + h ) - pc_at( e, sw - h ) ) / (2. * h);
}





template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dpcdsw_at_Numerical( Element<dim>* const e, double sw, double h ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    // deal with the more common case of a high water saturation first
    if ( sw > (1. - h) )
      return ( pc_at( e, 1. ) - pc_at( e, 1. - h)) / h;
  
    // low water saturation
    if ( sw < h )
      return ( pc_at( e, h ) - pc_at( e, 0. ) ) / h;

    // water saturation between the endpoints
    return ( pc_at( e, sw + h ) - pc_at( e, sw - h ) ) / (2. * h);
}





/**
    Computes derivative of the wetting phase relative permeability using central finite difference method.
*/
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrwds_Numerical( Element<dim>* const e, double h ) const
  {
    assert( h > 0. );
    assert( h <= 0.01 );
  
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
  
    // deal with the more common case of a high water saturation first
    if ( sw > (1. - h) )
      return ( krw_at( e, 1. ) - krw_at( e, 1. - h)) / h;
  
    // low water saturation
    if ( sw < h )
      return ( krw_at( e, h ) - krw_at( e, 0. ) ) / h;

    // water saturation between the endpoints
    return ( krw_at( e, sw + h ) - krw_at( e, sw - h ) ) / (2. * h);
  }
  
  
  
  
  
  
  
  
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrwds_at_Numerical( Element<dim>* const e, double sw, double h ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );
    assert( h > 0. );
    assert( h <= 0.01 );

    // deal with the more common case of a high water saturation first
    if ( sw > (1. - h))
      return ( krw_at( e, 1. ) - krw_at( e, 1. - h)) / h;
  
    // low water saturation
    if ( sw < h )
      return ( krw_at( e, h ) - krw_at( e, 0. ) ) / h;

    // water saturation between the endpoints
    return ( krw_at( e, sw + h ) - krw_at( e, sw - h ) ) / (2. * h);
 }
  
  
  
  
  
  
  
  
  
  
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrnds_Numerical( Element<dim>* const e, double h ) const
 {
    assert( h > 0. );
    assert( h <= 0.01 );
  
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );
  
    // deal with the more common case of a high water saturation first
    if ( sw > (1. - h) )
      return ( krn_at( e, 1. ) - krn_at( e, 1. - h)) / h;
  
    // low water saturation

    if ( sw < h )
      return ( krn_at( e, h ) - krn_at( e, 0. ) ) / h;

    // water saturation between the endpoints
    return ( krn_at( e, sw + h ) - krn_at( e, sw - h ) ) / (2. * h);
 }
  
  
  
  
  

  
  
  
template<size_t dim, template<size_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrnds_at_Numerical( Element<dim>* const e, double sw, double h ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );

    assert( h > 0. );
    assert( h <= 0.01 );
  
    // deal with the more common case of a high water saturation first
    if ( sw > (1. - h) )
      return ( krn_at( e, 1. ) - krn_at( e, 1. - h)) / h;
  
    // low water saturation
    if ( sw < h )
      return ( krn_at( e, h ) - krn_at( e, 0. ) ) / h;

    // water saturation between the endpoints
    return ( krn_at( e, sw + h ) - krn_at( e, sw - h ) ) / (2. * h);
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

  
template class ExperimentalSaturationFunctions<1U,FlowFunctionsModule3>;
template class ExperimentalSaturationFunctions<2U,FlowFunctionsModule3>;
template class ExperimentalSaturationFunctions<3U,FlowFunctionsModule3>;

template class ExperimentalSaturationFunctions<1U,FlowFunctionsModule6>;
template class ExperimentalSaturationFunctions<2U,FlowFunctionsModule6>;
template class ExperimentalSaturationFunctions<3U,FlowFunctionsModule6>;
  
  
} // csmp
