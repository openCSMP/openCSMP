#include "BrooksCoreyWithHysteresis.h"
#include "BrooksCoreyFrontVelocity.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {


template<size_t dim>
BrooksCoreyWithHysteresis<dim>::BrooksCoreyWithHysteresis()
 : acc_gravity_(9.8066)
 {
 }
 
 
/**
   we also need to initialize the material parameter key in the base class
   in this case the material parameter is used as the lambda parameter
   of the Brooks-Corey model
*/
template<size_t dim>
BrooksCoreyWithHysteresis<dim>::BrooksCoreyWithHysteresis( const PropertyDatabase<dim>& database,
                                                           const char* permeability,
                                                           double64 viscosity_nw, double64 viscosity_w,
                                                           double64 density_nw, double64 density_w, 
                                                           const char* lamda, const char* pc_entry )
 : pd_key(database.StorageKey(pc_entry)),
   lamda_key(database.StorageKey(lamda)),
   sormax_key(database.StorageKey("maximum residual oil saturation")),
   sat_previous_key(database.StorageKey("previous water saturation barycenter")),
   sat_inflection_key(database.StorageKey("inflection water saturation barycenter")),
   acc_gravity_(9.8066),
   on_scaning_curve(false),
   imbibing(false),
   sat_previous(std::numeric_limits<double64>::quiet_NaN()),
   sat_inflection(std::numeric_limits<double64>::quiet_NaN()),
   sormax(std::numeric_limits<double64>::quiet_NaN()),
   TwoPhaseModel<dim>(database, viscosity_nw, viscosity_w, density_nw, density_w,
                         permeability, "saturation water", 
                        "residual saturation non-wetting phase",
                        "residual saturation wetting phase")
 {
  if ( sat_previous_key.place != ELEMENT || sat_previous_key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "BrooksCoreyWithHysteresis(constructor)", 
                           "previous water saturation barycenter", 
                           " variable must be scalar element property" );
                           
  if ( sormax_key.place != NODE || sormax_key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "BrooksCoreyWithHysteresis(constructor)", 
                           "maximum residual oil saturation", 
                           " variable must be scalar node property" ); 

  if ( sat_inflection_key.place != ELEMENT || sat_inflection_key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "BrooksCoreyWithHysteresis(constructor)", 
                           "inflection water saturation barycenter", 
                           " variable must be scalar element property" ); 
 }



template<size_t dim>
BrooksCoreyWithHysteresis<dim>::BrooksCoreyWithHysteresis( const PropertyDatabase<dim>& database,
                                                           const char* lamda, const char* pc_entry )
 : pd_key(database.StorageKey(pc_entry)),
   lamda_key(database.StorageKey(lamda)),
   sormax_key(database.StorageKey("maximum residual oil saturation")),
   sat_previous_key(database.StorageKey("previous water saturation barycenter")),
   sat_inflection_key(database.StorageKey("inflection water saturation barycenter")),
   acc_gravity_(9.8066),
   on_scaning_curve(false),
   imbibing(false),
   sat_previous(std::numeric_limits<double64>::quiet_NaN()),
   sat_inflection(std::numeric_limits<double64>::quiet_NaN()),
   sormax(std::numeric_limits<double64>::quiet_NaN()),
   TwoPhaseModel<dim>(database, "permeability",  
                        "viscosity oil", "viscosity water",
                        "density oil", "density water", "saturation water",
                        "residual saturation non-wetting phase",
                        "residual saturation wetting phase" )
 {
 
  if ( sat_previous_key.place != ELEMENT || sat_previous_key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "BrooksCoreyWithHysteresis(constructor)", 
                           "previous water saturation barycenter", 
                           " variable must be scalar element property" );
                           
  if ( sormax_key.place != NODE || sormax_key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "BrooksCoreyWithHysteresis(constructor)", 
                           "maximum residual oil saturation", 
                           " variable must be scalar node property" ); 

  if ( sat_inflection_key.place != ELEMENT || sat_inflection_key.type != SCALAR )
     throw csmp::Exception( FATAL_ERROR, "BrooksCoreyWithHysteresis(constructor)", 
                           "inflection water saturation barycenter", 
                           " variable must be scalar element property" ); 
 }



template<size_t dim>
BrooksCoreyWithHysteresis<dim>::~BrooksCoreyWithHysteresis()
 {
 }
 
 
 

// Testing
// linear & BC relperms, f's dfds's are O.K.


template<size_t dim>
void BrooksCoreyWithHysteresis<dim>::Initialize( Element<dim>& e )
 {
    TwoPhaseModel<dim>::k_   = e.Read( TwoPhaseModel<dim>::perm_key_ );
    TwoPhaseModel<dim>::swr_ = e.Read( TwoPhaseModel<dim>::swr_key_ );
    TwoPhaseModel<dim>::snr_ = e.Read( TwoPhaseModel<dim>::snr_key_ );
    pm1 = e.Read( pd_key );
    pm2 = e.Read( lamda_key );
    sormax = e.Read( sormax_key );
    sat_previous = e.Read( sat_previous_key );
    sat_inflection = e.Read( sat_inflection_key );
    
    imbibing = ( this->sat_ > sat_previous );
    if ( this->sat_ < this->swr_ ) sat_inflection = this->swr_;
    on_scaning_curve = !( sat_inflection == -1. );
           
    if ( imbibing & !on_scaning_curve )
    {
       // set sat_inflection = sw_previous
       e.Store(sat_inflection_key, makeScalar(e.Status(sat_previous_key),e.Read(sat_previous_key)) );
    }

    if ( !(!imbibing & !on_scaning_curve) )
    {
      CC = 1. / ( sormax / ( 1. - this->swr_ ) ) - 1.;
      CK = 1. / ( sormax - this->snr_ ) - 1. / ( 1. - this->swr_ - this->snr_ );
      Sot = this->snr_ + ( 1. - sat_inflection - this->snr_ ) / ( 1. + CK * ( 1. - sat_inflection - this->snr_ ) );
      Snorm = sormax + ( 1. - this->sat_ - Sot ) * ( 1. - this->swr_ - sormax ) / ( 1. - sat_inflection - Sot );
      kroDswinflection = std::pow( 1. - ( sat_inflection - this->swr_ ) / ( 1. - this->swr_ - this->snr_ ), 2.) *
                         ( 1. - std::pow( ( sat_inflection - this->swr_ ) / ( 1. - this->swr_ - this->snr_ ),
                         ( ( 2. + this->pm2 ) / this->pm2 ) ) );
      Sof_1_Snorm = 0.5 * ( Snorm / ( 1. - this->swr_ ) - sormax / ( 1. - this->swr_ ) + std::pow(
                    std::pow( Snorm / ( 1. - this->swr_ ) - sormax / ( 1. - this->swr_ ), 2. ) + 4. / CC *
                    ( Snorm / ( 1. - this->swr_ ) - sormax / ( 1. - this->swr_ ) ), 0.5 ) );
      kroI_1_Snorm = std::pow( Sof_1_Snorm, 2.) * ( 1. - std::pow( 1. - Sof_1_Snorm,
                     ( 1. + 2. / this->pm2 ) ) );
    }
    
    if ( !imbibing & on_scaning_curve )
    {
       if ( this->sat_ < sat_inflection )
       {
       // set sat_inflection = -1
          ScalarVariable  sc( PLAIN, -1. );
          e.Store(sat_inflection_key, sc);
       }
    }
 } // end Initialize




template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::krw_Phase() const
 { 
    if ( TwoPhaseModel<dim>::seff_ <= 0. ) return static_cast<double64>(0.);
    if ( TwoPhaseModel<dim>::seff_ >= 1. ) return static_cast<double64>(1.);

    // pm2 = lambda, the Brooks-Corey parameter
    return std::pow( TwoPhaseModel<dim>::seff_,
                    (2. + 3. * pm2) / pm2 ); 
 }


template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::krn_Phase() const
 { 
    if ( imbibing )
    {
       if ( !on_scaning_curve )
       {
          if ( this->sat_ >= ( 1. - Sot ) ) return static_cast<double64>(0.);
          if ( this->sat_ < this->swr_ ) return static_cast<double64>(1.);
          return kroI_1_Snorm * kroDswinflection;
       }
       else
       { 
          if ( this->sat_ >= ( 1. - Sot ) ) return static_cast<double64>(0.);
          if ( this->sat_ < this->swr_ ) return static_cast<double64>(1.);
          return kroI_1_Snorm * kroDswinflection;
       }
    }
    else
    {
       if ( !on_scaning_curve )
       {
          if ( this->sat_ <= this->swr_ ) return static_cast<double64>(1.);
          if ( this->sat_ >= ( 1 - this->snr_ ) ) return static_cast<double64>(0.);
          return std::pow( 1. - ( this->sat_ - this->swr_ ) / ( 1. - this->swr_ - this->snr_ ), 2.) *
                 ( 1. - std::pow( ( this->sat_ - this->swr_ ) / ( 1. - this->swr_ - this->snr_ ),
                 ( ( 2. + this->pm2 ) / this->pm2 ) ) );
       }
       else
       {
          if ( this->sat_ <= this->swr_ ) return static_cast<double64>(1.);
          if ( this->sat_ >= ( 1. - Sot ) ) return static_cast<double64>(0.);
          if ( this->sat_ < sat_inflection )
          {
             if ( this->sat_ >= ( 1. - this->snr_ ) ) return static_cast<double64>(0.);
             return std::pow( 1. - ( this->sat_ - this->swr_ ) / ( 1. - this->swr_ - this->snr_ ), 2.) *
                    ( 1. - std::pow( ( this->sat_ - this->swr_ ) / ( 1. - this->swr_ - this->snr_ ),
                    ( ( 2. + this->pm2 ) / this->pm2 ) ) );
          }
          return kroI_1_Snorm * kroDswinflection;
       }
    } 
 }




// for the wetting phase
template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::MaxFractionalFlowDerivative() const
 {
    return static_cast<double64>(5.6); // as computed with dfds method
 } 
  

// pc covers the full saturation range, pc is capped if sw<swr     
// tested: O.K.     
template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::pc_Phase( size_t ) const
{
   // compute sw_eff for which pc = 40MPa, seff_min = (pc/pd)^-lamda
   // applying the limit on capillary pressure 
   if ( TwoPhaseModel<dim>::seff_ <= std::pow( TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ / pm1, -pm2 ) )
     return TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
     
   return pm1 * std::pow( TwoPhaseModel<dim>::seff_, -1. / pm2 );
}



// dpcdS covers the full saturation range, dpcdS is capped if sw<swr
// tested: O.K.    
template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::dpcds_Phase( size_t phase ) const
{
   // the wetting phase has no capilllary pressure
   if ( phase == 1U ) return static_cast<double64>(0.);

   // compute Se for which pc = 40MPa, seff_min = (pc/pd)^-lamda
   const double64 Se_min = std::pow( TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ / pm1, -pm2 );
   // applying the limit on capillary pressure and its derivative
   if ( TwoPhaseModel<dim>::seff_ <= Se_min )
     return pm1 * std::pow( Se_min, -1. / pm2 ) / (Se_min * pm2);

   // computing the unlimited capillary pressure derivative
   return pm1 * std::pow( TwoPhaseModel<dim>::seff_, -1. / pm2 ) /
          (TwoPhaseModel<dim>::seff_ * pm2);
}
    


 

template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::dfds() const 
{
   const double64 seff = TwoPhaseModel<dim>::seff_;

   if ( seff <= TwoPhaseModel<dim>::swr_ or seff >= (1. - TwoPhaseModel<dim>::snr_) ) return static_cast<double64>(0.);

   // the linear model case
   if ( pm2 == static_cast<double64>(0.) )
     return 1. / (1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_);
   
	const double64 a1 = std::pow ( seff , ( 2. + pm2 ) / pm2 );
	const double64 a2 = std::pow ( seff , ( 2. + 2. * pm2 ) / pm2 );
	const double64 a3 = std::pow ( seff , ( 2. + 3. * pm2 ) / pm2 );
    const double64 a4 = ( ((1.-seff) * (1.-seff)) * (1. - a1 ) + TwoPhaseModel<dim>::mun_ / TwoPhaseModel<dim>::muw_ * a3);
	
    return -TwoPhaseModel<dim>::mun_ / TwoPhaseModel<dim>::muw_ * 1. / pm2 * 1. /
           (1. -  TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_) *
	       ((1.-seff) * a2 * (2. * -(1.-seff) + pm2 * (-(1.-seff) - 2. + 2. * a1))) / (a4 * a4);
	
}  // end dfdS_Phase




// tested: O.K.
template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::dGds() const 
{
   // 2. compute dGdS for the Brooks-Corey model
   double64 seffw(TwoPhaseModel<dim>::sat_ / (1.-TwoPhaseModel<dim>::swr_-TwoPhaseModel<dim>::snr_));
   if ( TwoPhaseModel<dim>::sat_ > 1.-TwoPhaseModel<dim>::snr_ ) seffw = 1.;
   if ( TwoPhaseModel<dim>::sat_ < TwoPhaseModel<dim>::swr_ )    return 0.; // phase is immobile, no flow

   // see Maple worksheet 'BrooksCorey_G&dGdS.mw'
   const double64 t1 = 1. - seffw;
   const double64 t3 = 1. /  pm2; // pm2 = Brooks-Corey lambda
   const double64 t5 = std::pow( TwoPhaseModel<dim>::sat_,  ((2. + pm2) * t3) ); // sat = sw
   const double64 t6 = 0.1e1 - t5;
   const double64 t8 = 0.1e1 / TwoPhaseModel<dim>::mun_;
   const double64 t9 =  t1 * t6 * t8;
   const double64 t11 = 2. + 3. * pm2;
   const double64 t12 = t11 * t3;
   const double64 t13 = std::pow( seffw,  t12);
   const double64 t14 = 1. /  TwoPhaseModel<dim>::muw_;
   const double64 t15 = t13 * t14;
   const double64 t16 = t1 * t1;
   const double64 t17 = t16 * t6;
   const double64 t18 = t17 * t8;
   const double64 t19 = t18 +  t15;
   const double64 t20 = 0.1e1 / t19;
   const double64 t26 = 1. / seffw;
   const double64 t31 = t19 * t19;

   return -0.2e1 * t9 * t15 * t20 + t17 * t8 * t13 * t12 * t26 * t14 * t20 - t18 *  t15 / 
           t31 * (-0.2e1 * t9 + (t13 * t11 * t3 * t26 * t14));

}  // end dGdS_Phase

//cout <<"\nsw "<< TwoPhaseModel<dim>::sat <<", bcp "<< pm2 <<", dGdS "<< dGdS <<" "; cout.flush();




template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::ShockSpeed() const
 {
    return BrooksCoreyFrontVelocity().ShockVelocityMultiplier( pm2, 
                                                 TwoPhaseModel<dim>::mun_/TwoPhaseModel<dim>::muw_ );
 }


// look this one up in the book by Randy LeVeque
template<size_t dim>
double64 BrooksCoreyWithHysteresis<dim>::ShockHeight() const
 {
    cerr <<"\nBrooksCorey<dim>::ShockHeight: not implemented yet."<< endl;
    return -1.; 
 }




template<size_t dim>
void BrooksCoreyWithHysteresis<dim>::Out( std::ostream& os, size_t phase ) const
 {
    TwoPhaseModel<dim>::Out(os, phase);
    os <<"\nBrooksCorey<"<< dim << ">::Out: Additional properties: "<< endl;
    os <<"\nelement properties:";
    os <<"\n       capillary entry pressure, pd: "<< pm1;
    os <<"\n      Brooks-Corey lambda parameter: "<< pm2;
    os <<"\n             MAX_CAPILLARY_PRESSURE: "<< TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ << endl << endl;

 } // end Out
 
 

template class BrooksCoreyWithHysteresis<1U>;
template class BrooksCoreyWithHysteresis<2U>;
template class BrooksCoreyWithHysteresis<3U>;

} // end namespace csmp


























