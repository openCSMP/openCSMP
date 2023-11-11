
#include "ExperimentalSaturationFunctions.h"
#include "FlowFunctionsModule.h"
#include "Element.h"
#include "ErrorHandler.h"
#include "CSMP_highLevelUtilities.h"
#include "Model.h"
#include "CSMP_mathUtilities.h"
#include "Region.h"

using namespace std;

namespace csmp {
  
/**
   The default constructor of Exoperimental Saturation Functions With Hysteresis class
*/
template<uint32_t dim, template<uint32_t> class USER>
ExperimentalSaturationFunctions<dim,USER>::ExperimentalSaturationFunctions( const char* filename )
  {
     InitialiseReservoirRockTypes(filename);
  }
  


/**
 
 reading the rocktype file
 
 */
template<uint32_t dim, template<uint32_t> class USER>
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

    for ( int rocktype=0; rocktype < number_of_tables; rocktype++ )
      {
   /* TODO: rocktype number is not read but inferred from the number of tables read!
          size_t rocktype;
          rt_file >> rocktype; // check that this is indeed a rocktype number
          assert( rocktype < number_of_tables );
          std::cout << "\nReading data table for RRT: " << rocktype << "\n";
  */
          uint32_t number_of_entries{0U};
          rt_file >> number_of_entries;
          
          double kro_start_derivative, kro_end_derivative, krw_start_derivative, krw_end_derivative, pc_start_derivative, pc_end_derivative;
          rt_file >> kro_start_derivative >> kro_end_derivative >> krw_start_derivative >> krw_end_derivative >> pc_start_derivative >> pc_end_derivative;
          cout << "Derivatives (krnw0,krnw1,krw0,krw1,dpcds0,dpcds1):\t" << kro_start_derivative << "\t" << kro_end_derivative;
          cout << "\t" << krw_start_derivative << "\t" << krw_end_derivative << "\t" << pc_start_derivative << "\t" << pc_end_derivative << "\n";
   
          if ( fabs(pc_start_derivative) > max_derivative_ ) {
               cerr <<"\n\tsaturation derivative of pc at sw_min: "<< pc_start_derivative;
               csmp_error.Note( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
            }
          if ( fabs(pc_end_derivative) > max_derivative_ ) {
               cerr <<"\n\tsaturation derivative of pc at sw_max: "<< pc_end_derivative;
               csmp_error.Note( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
            }
          if ( fabs(krw_start_derivative) > max_derivative_ ) {
               cerr <<"\n\tsaturation derivative of krw at sw_min: "<< krw_start_derivative;
               csmp_error.Note( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
            }
          if ( fabs(krw_end_derivative) > max_derivative_ ) {
               cerr <<"\n\tsaturation derivative of krw at sw_max: "<< krw_end_derivative;
               csmp_error.Note( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
            }
          if ( fabs(kro_start_derivative) > max_derivative_ ) {
               cerr <<"\n\tsaturation derivative of krn at sw_min: "<< kro_start_derivative;
               csmp_error.Note( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
            }
          if ( fabs(kro_end_derivative) > max_derivative_ ) {
               cerr <<"\n\tsaturation derivative of krn at sw_max: "<< kro_end_derivative;
               csmp_error.Note( ERROR, "ExperimentalSaturationFunctions<dim,USER>::InitialiseReservoirRockTypes:", "derivative out of range; check input table.");
            }
        
          std::cout << "sw\tkro\tkrw\tpc\n";
        
          std::vector<double> sw, kro, krw, pc;
          double swr(0.), snr(0.);
          double pre_sw_value(0.);
          double sw_value(0.), kro_value, krw_value, pc_value;
          bool swr_found(false), snr_found(false);
          for ( uint32_t n{0U}; n < number_of_entries; n++ )
            {
                pre_sw_value = sw_value;
                rt_file >> sw_value >> kro_value >> krw_value >> pc_value;

                std::cout << sw_value << "\t" << kro_value << "\t" << krw_value << "\t" << pc_value <<"\n";
              
                sw.push_back(sw_value);
                kro.push_back(kro_value);
                krw.push_back(krw_value);
                pc.push_back(pc_value);

                if(!swr_found && krw_value > 0.0) {swr = pre_sw_value; swr_found = true;}
                if(!snr_found && kro_value == 0.0) {snr = 1.0 - sw_value; snr_found = true;}
            }

          std::cout << "swr = " <<swr<<", snr = "<<snr<<endl;
          std::cout << "\n";

          input_sw_[rocktype] = sw;
          input_pc_[rocktype] = pc;

          csmp::CubicSpline kr1_spline, kr2_spline, pc_spline;
          kr1_spline.Initialize( sw, krw, krw_start_derivative, krw_end_derivative );
          kr2_spline.Initialize( sw, kro, kro_start_derivative, kro_end_derivative );
          pc_spline.Initialize( sw, pc, pc_start_derivative, pc_end_derivative );

          kr1_[rocktype] = kr1_spline;
          kr2_[rocktype] = kr2_spline;
          pc_[rocktype] = pc_spline;
          swr_[rocktype] = swr;
          snr_[rocktype] = snr;
       }

  //Out();
  // limiting the rocktype number range in the property database to the actual maximum value

  cout << "\nExperimentalSaturationFunctions::InitialiseReservoirRockTypes: rock types initialised successfully from file '"<< rt_file_name <<".'"<< endl;
  return pc_.size();
  
} // end InitialiseReservoirRockTypes

 


/**
    returns for how many rocktypes information is stored in the underlying data arrays.
*/
template<uint32_t dim, template<uint32_t> class USER>
size_t ExperimentalSaturationFunctions<dim,USER>::RockTypes() const
 { return pc_.size(); }


  
  
/**
    returns the rocktype associated with the current Element
*/
template<uint32_t dim, template<uint32_t> class USER>
size_t ExperimentalSaturationFunctions<dim,USER>::RockType( Element<dim>* const e ) const
 {
    const size_t rocktype = static_cast<uint32_t>(e->Read( User()->key_RRT ));
    return rocktype;
 }







/**
    (sw-swr) / (1-swr-snr); use only in 2-phase flow simulations
 
    TODO: this needs to be computed from the curves in the input file
*/
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation( Element<dim>* const e ) const
 {
   /*
    if ( e->Read(User()->key_srH2O) > 0. || e->Read(User()->key_srCO2) > 0. )
      throw csmp::Exception( ERROR, "ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation",
                             "Do not set end-point saturation values in conjunction with experimental flow functions as they are defined implicitly");
    */
    return e->PropertyValueAtBaryCenter( User()->key_sH2O );
 }
  
  
  
 
    /// (sw-swr) / (1-swr-snr); use only in 2-phase flow simulations
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation_at( Element<dim>* const e, double sw ) const
 {
    /*
    if ( e->Read(User()->key_srH2O) > 0. || e->Read(User()->key_srCO2) > 0. )
      throw csmp::Exception( ERROR, "ExperimentalSaturationFunctions<dim,USER>::EffectiveSaturation_at",
                             "Do not set end-point saturation values in conjunction with experimental flow functions as they are defined implicitly");
    */
    return sw;
 }







/**
    @note Note that capillary pressure also exists outside of the effective saturation range
*/
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::pc( Element<dim>* const e ) const
 {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return min( pc_.find( RockType(e) )->second.Value(sw), max_capillary_pressure_ );
 }





template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::pc_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return min( pc_.find( RockType(e) )->second.Value(sw), max_capillary_pressure_ );
 }



/**
    @attention this treatment assumes that the capillary pressure derivative with regard to the
    wetting phase saturation is always negative.
*/
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dpcds( Element<dim>* const e ) const
  {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( pc_.find( RockType(e) )->second.Derivative(sw), -max_derivative_ );
}
  



template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dpcds_at( Element<dim>* const e, double sw ) const
  {
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( pc_.find( RockType(e) )->second.Derivative(sw), -max_derivative_ );
}





template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krw( Element<dim>* const e ) const
{
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( kr1_.find( RockType(e) )->second.Value( sw ), 0. );
}
  


template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krw_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );
   
    // as defined by spline curves over full saturation range
    return max( kr1_.find( RockType(e) )->second.Value( sw ), 0. );
 }



  
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krn( Element<dim>* const e ) const
  {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( kr2_.find( RockType(e) )->second.Value( sw ), 0. );
  }

 
 
 
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::krn_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );

    // as defined by spline curves over full saturation range
    return max( kr2_.find( RockType(e) )->second.Value( sw ), 0. );
 }
  


/**
 
    calculating the 1st derivative of water relative permeability at the element barycentre.
 
*/
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrwds( Element<dim>* const e ) const
  {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O  );
    assert( sw >= 0. );
    assert( sw <= 1. );

    return min( kr1_.find( RockType(e) )->second.Derivative( sw ), max_derivative_ );
  }





template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrwds_at( Element<dim>* const e, double sw ) const
{
    assert( sw >= 0. );
    assert( sw <= 1. );

    return min( kr1_.find( RockType(e) )->second.Derivative( sw ), max_derivative_ );
 }




/**
    @note by contrast with the wetting phase relative permeability, the non-wetting phase
    kri is expected to have a negative slope w.r.t. water saturation.
*/
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrnds( Element<dim>* const e ) const
 {
    const double sw = e->PropertyValueAtBaryCenter( User()->key_sH2O );

    return max( kr2_.find( RockType(e) )->second.Derivative( sw ), -max_derivative_ );
 }




/**
    1st derivative of non-wetting phase relative permeability as a function of saturation.
*/
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::dkrnds_at( Element<dim>* const e, double sw ) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );

    return max( kr2_.find( RockType(e) )->second.Derivative( sw ), -max_derivative_ );
 }



// NUMERICAL DERIVATIVE CALCULATIONS


template<uint32_t dim, template<uint32_t> class USER>
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





template<uint32_t dim, template<uint32_t> class USER>
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
template<uint32_t dim, template<uint32_t> class USER>
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
  
  
  
  
  
  
  
  
template<uint32_t dim, template<uint32_t> class USER>
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
  
  
  
  
  
  
  
  
  
  
template<uint32_t dim, template<uint32_t> class USER>
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
  
  
  
  
  

  
  
  
template<uint32_t dim, template<uint32_t> class USER>
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
  
  
  
  
  



  
template<uint32_t dim, template<uint32_t> class USER>
void ExperimentalSaturationFunctions<dim,USER>::Out() const
 {
    cout <<"\nExperimentalSaturationFunctions<dim,USER>::Out: current rock type data:\n";
    cout <<"\nrelative permeability of phase 1 (krw(sw)), for all rock types:\n";
    int number(0);
    for ( auto it=kr1_.begin(); it!=kr1_.end(); it++ ) {
        cout <<"\nrock type "<< it->first;
        it->second.Out();
      }
    cout <<"\nrelative permeability of phase 2 (krnw(sw)), for all rock types:\n";
    number = 0;
    for ( auto it=kr2_.begin(); it!=kr2_.end(); it++ ) {
        cout <<"\nrock type "<< it->first;
        it->second.Out();
      }
    cout <<"\ncapillary pressure curve, pc(sw), for all rock types:\n";
    number = 0;
    for ( auto it=pc_.begin(); it!=pc_.end(); it++ ) {
        cout <<"\nrock type "<< it->first;
        it->second.Out();
      }

    cout <<"\nmaximum slope of derivative curves: "<< max_derivative_ << endl;
 }




/// return wetting-phase saturation based on effective saturation
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::seff_to_sw( Element<dim>* const e, double seff ) const
{
    double swr = e->Read(User()->key_srH2O);
    double snr = e->Read(User()->key_srCO2);
    return seff * (1. - swr - snr) + swr;

}



/// inverse capillary pressure function
template<uint32_t dim, template<uint32_t> class USER>
double ExperimentalSaturationFunctions<dim,USER>::sw_from_pc_at( Element<dim>* const e, double pc, double sw ) const
{
    // only for the min saturation of water precautions are needed
    // the actual saturation is used instead of the effective saturation
    if ( pc >= MaxCapillaryPressure() )
      return seff_to_sw(e, 0.);

    if ( pc == 0.0 )
      return seff_to_sw(e, 1.);

    std::vector<double> input_pc = input_pc_.find( RockType(e) )->second;
    std::vector<double> input_sw = input_sw_.find( RockType(e) )->second;
    size_t i(0);
    for( i = 1; i < input_pc.size() ; ++i )
      if( ( pc <=input_pc[i-1] ) && (pc > input_pc[i]))
        break;

    double x = pc, x1 = input_pc[i-1], x2 = input_pc[i], y1 = input_sw[i-1], y2 = input_sw[i];
    double k1 = 1.0 / dpcds_at(e,y1);
    double k2 = 1.0 / dpcds_at(e,y2);
    return splineValue( x, x1, x2, y1, y2, k1, k2);
  }



///initialise swr and snr from input table data
template<uint32_t dim, template<uint32_t> class USER>
void ExperimentalSaturationFunctions<dim,USER>::initialiseResidualSaturations( Model<dim>& model )
{
    Region<dim>& mref = model.Region("Model");
    for ( auto eit = mref.CellsBegin(); eit!= mref.CellsEnd(); eit++ ) {
      double swr = swr_.find( RockType(*eit) )->second;
      double snr = snr_.find( RockType(*eit) )->second;
      (*eit)->Store( User()->key_srH2O, makeScalar(PLAIN, swr) );
      (*eit)->Store( User()->key_srCO2, makeScalar(PLAIN, snr) );
    }
    cout<<"ExperimentalSaturationFunctions<dim,USER>:finished assigning residual saturations to all elements"<<endl;
    printRangeOfVariable( model, "residual saturation carbonic phase" );
    printRangeOfVariable( model, "residual saturation aqueous phase" );
  }

  
template class ExperimentalSaturationFunctions<1U,FlowFunctionsModule3>;
template class ExperimentalSaturationFunctions<2U,FlowFunctionsModule3>;
template class ExperimentalSaturationFunctions<3U,FlowFunctionsModule3>;

template class ExperimentalSaturationFunctions<1U,FlowFunctionsModule6>;
template class ExperimentalSaturationFunctions<2U,FlowFunctionsModule6>;
template class ExperimentalSaturationFunctions<3U,FlowFunctionsModule6>;
  
  
} // csmp
