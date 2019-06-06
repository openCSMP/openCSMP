#include "HeterogeneityAndRateAwareModel.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {

template<size_t dim>
HeterogeneityAndRateAwareModel<dim>::HeterogeneityAndRateAwareModel( const PropertyDatabase<dim>& database,
                                                                     const char* rocktype, const char* total_velocity,
                                                                     const char* pc_entry,
                                                                     const bool sw_ro_mu_placement )
  
 : TwoPhaseModel<dim>(database, "permeability",  
                     "viscosity oil", "viscosity water",
                     "density oil", "density water", "saturation water",
                     "residual saturation non-wetting phase",
                     "residual saturation wetting phase",
                      sw_ro_mu_placement ),

   RRT_key_(database.StorageKey(rocktype)),
   vt_key_(database.StorageKey(total_velocity)),
   bcp_key_(database.StorageKey("brooks corey parameter")),
   pd_key_(database.StorageKey(pc_entry))
{
}



template<size_t dim>
HeterogeneityAndRateAwareModel<dim>::~HeterogeneityAndRateAwareModel()
 {
 }
 
 
 

/// Reading parameters for the linear & BC relperm models
template<size_t dim>
void HeterogeneityAndRateAwareModel<dim>::Initialize( const Element<dim>& e )
 {
    TwoPhaseModel<dim>::swr_ = e.Read( TwoPhaseModel<dim>::swr_key_ );
    TwoPhaseModel<dim>::snr_ = e.Read( TwoPhaseModel<dim>::snr_key_ );

    if( TwoPhaseModel<dim>::tensor_permeability_){
        e.Read( TwoPhaseModel<dim>::perm_key_, TwoPhaseModel<dim>::K_);
        TwoPhaseModel<dim>::k_ = TwoPhaseModel<dim>::K_.Trace()/static_cast<double64>(dim);
    }else{
        TwoPhaseModel<dim>::k_ = e.Read( TwoPhaseModel<dim>::perm_key_ );
        TwoPhaseModel<dim>::K_.operator=( VectorVariable<dim>(PLAIN, TwoPhaseModel<dim>::k_ ) );
    }

    entry_pressure_ = e.Read( pd_key_ );
    lambda_ = e.Read( bcp_key_ );

    // if properties are discretized on element
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
        TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
        TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
        TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
        TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
        TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );     
    }

 } // end Initialize





/*
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::dkrwds_Phase() const
 {
    const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    //  switch to linear relperm model if lambda = 0
    if ( lambda_ == 0.0 )
        return seff_mult;

    return ( 2./lambda_ + 3.0 )*std::pow( TwoPhaseModel<dim>::seff_, 2./lambda_ + 2.0 )*seff_mult;
 }


template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::dkrnds_Phase() const
 {
    const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    //  switch to linear relperm model if lambda = 0
    if ( lambda_ == 0. )
        return -seff_mult;

    return 2.0*(TwoPhaseModel<dim>::seff_ - 1.0)*
            ( 1.0 + pow(TwoPhaseModel<dim>::seff_, 2.0/lambda_)*( 1.0/lambda_ + 1./2. - TwoPhaseModel<dim>::seff_*(1.0/lambda_ + 3./2.)))*
            seff_mult;

 }
*/





/// pc covers the full saturation range, pc is capped based on maximum dpcds 
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::pc_Phase() const
{
    // linear relperm model
    if ( lambda_ == 0. ) {

        if ( TwoPhaseModel<dim>::seff_ <= 0. )
            return TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;

        if ( TwoPhaseModel<dim>::seff_ >= 1. )
            return entry_pressure_;

        return entry_pressure_ + ( 1. - TwoPhaseModel<dim>::seff_ ) * ( TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ - entry_pressure_ );
   }

   // for zero entry pressure capillary pressure always is zero
   if ( entry_pressure_ == 0. ) return 0.;

} // end pc_Phase





/// dpcdS covers the full saturation range, dpcds is capped based on maximum dpcds
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::dpcds_Phase() const
{

    // linear relperm model
     if ( lambda_ == 0. ){

         if ( ( TwoPhaseModel<dim>::seff_ < 0. ) || ( TwoPhaseModel<dim>::seff_ > 1. ) )
             return 0.;

         const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

         return -(TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ - entry_pressure_ )*seff_mult;

    }

    // for zero entry pressure dpcds is zero
    if ( entry_pressure_ == 0. )
       return 0.;

 } // end
 
 
 
 

/*
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::Sw_Phase( double64 pc ) const
{

    // for linear relperm model
    if ( lambda_ == 0. )
    {
        // not unique solution
        if ( pc_max_ == entry_pressure_ )
            return std::min( std::max( TwoPhaseModel<dim>::sat_, TwoPhaseModel<dim>::swr_), 1.0 - TwoPhaseModel<dim>::snr_);

        if ( pc >= pc_max_ ){
            TwoPhaseModel<dim>::seff_ = 0.0;
            return TwoPhaseModel<dim>::SeffToSw();
        }

        if ( pc <= entry_pressure_ ){
            TwoPhaseModel<dim>::seff_ = 1.0;
            return TwoPhaseModel<dim>::SeffToSw();
        }

        TwoPhaseModel<dim>::seff_ =  std::min( std::max( 1.0 - ( pc - entry_pressure_ )/( pc_max_ - entry_pressure_ ), 0.0), 1.0 );
        return TwoPhaseModel<dim>::SeffToSw();
   }

   // for zero entry pressure capillary pressure always is zero, so the saturation is not unique
   if ( entry_pressure_ == 0. )
       return TwoPhaseModel<dim>::sat_;

   const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double64 pcmax = entry_pressure_ * pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                           ( -1. / ( 1. + lambda_ ) ) );

   if ( pc >= pcmax )
   {
      // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
       const double64 Se_min =  pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                ( lambda_ / ( 1. + lambda_ ) ) );
       // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
       TwoPhaseModel<dim>::seff_ = Se_min - ( pc - pcmax )/TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_*seff_mult ;

       return TwoPhaseModel<dim>::SeffToSw();
   }

   TwoPhaseModel<dim>::seff_ =  std::pow(pc/entry_pressure_, -lambda_);

   return TwoPhaseModel<dim>::SeffToSw();

}
*/



/*
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::dsdpc_Phase( double64 pc ) const
{
    // linear relperm model
    if( lambda_ == 0.0 )
    {
       const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
        return -1.0/( TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE - entry_pressure_ )/seff_mult;

    }

    if ( entry_pressure_ == 0. )
        return 0.0;

    const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double64 Se_min =  pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                ( lambda_ / ( 1. + lambda_ ) ) );

    const double64 Seff = (TwoPhaseModel<dim>::Sw_at( pc )- TwoPhaseModel<dim>::swr_ ) *seff_mult;

    if( Seff<= Se_min )
        return -1.0/TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_;

    return -lambda_/entry_pressure_ * pow(pc/entry_pressure_, -lambda_ - 1.0 )/seff_mult;

}
*/




/**

   Calculates the derivative of water fractional flow respect to water saturation for 
   Brooks-Corey and Linear relative permeability models. The derivative is conducted by
   analytical differentiation.

*/

template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::dfds() const
{
	
}  // end dfdS





template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::dGds() const
{

   const double64 seff = TwoPhaseModel<dim>::seff_;
   
   //if ( ( seff <= static_cast<double64>(0.) ) || ( seff >= static_cast<double64>(1.) ) )
   //    return static_cast<double64>(0.);

   const double64 lambda = lambda_;

   const double64 mob_w = krw_Phase() / TwoPhaseModel<dim>::muw_;
   const double64 mob_n = krn_Phase() / TwoPhaseModel<dim>::mun_;
   
   // linear model
   if ( lambda == static_cast<double64>(0.) )
   {
       return ( ( mob_n * mob_n / TwoPhaseModel<dim>::muw_ ) - ( mob_w * mob_w / TwoPhaseModel<dim>::mun_ ) ) / 
              ( ( mob_w + mob_n ) * ( mob_w + mob_n ) ) / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ );
   }
   
   // Brooks-Corey 
   const double64 dmob_wdsw = ( 3. + 2. / lambda ) * std::pow ( seff , 2. + 2. / lambda ) / 
                             ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) / TwoPhaseModel<dim>::muw_;
   const double64 dmob_ndsw = - ( 2 * ( 1 - seff ) * ( 1 - std::pow ( seff , 1. + 2. / lambda ) ) + ( 1 - seff ) * ( 1 - seff ) 
                             * ( 1. + 2. / lambda ) * std::pow ( seff ,2. / lambda ) ) / 
                             ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) / TwoPhaseModel<dim>::mun_;
 
   return ( dmob_wdsw * ( mob_n * mob_n ) + dmob_ndsw * ( mob_w * mob_w ) ) / ( ( mob_w + mob_n ) * ( mob_w + mob_n ) );

}  // end dGdS_Phase


/// for the wetting phase
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::MaxFractionalFlowDerivative() const
 {
    // linear relperm model
    // in the linear relperm model maximum dfds belongs to seff=0 or seff=1 depending on the viscosity of water and oil
    // if mu_oil>mu_water maximum dfds is dfds@seff=0. and if mu_oil<mu_water maximum dfds is dfds@seff=1.
    if ( lambda_ == static_cast<double64>(0.) )
        return ( TwoPhaseModel<dim>::mun_ >= TwoPhaseModel<dim>::muw_ ) ?
               ( TwoPhaseModel<dim>::mun_ / TwoPhaseModel<dim>::muw_ / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) ) :
               ( TwoPhaseModel<dim>::muw_ / TwoPhaseModel<dim>::mun_ / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    return static_cast<double64>(5.6); // as computed with dfds method
 }



template<size_t dim>
void HeterogeneityAndRateAwareModel<dim>::Out( size_t phase ) const
 {
    TwoPhaseModel<dim>::Out(phase);
    cout <<"\nHeterogeneityAndRateAwareModel<"<< dim << ">::Out: Additional properties: "<< endl;
    cout <<"\nelement properties:";
    cout <<"\n       capillary entry pressure, pd: "<< entry_pressure_;
    cout <<"\n      Brooks-Corey lambda parameter: "<< lambda_;
    cout <<"\n  MAX_CAPILLARY_PRESSURE_DERIVATIVE: "<< TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_ << endl << endl;

 } // end Out
 
 



// ===============================================================================
//
//      NEW FUNCTIONS - STANFORD VISIT
//
// ===============================================================================

// CO2 RATE-DEPENDENT RELPERM



/**
    returns CL water saturation value below which water will enter the fine-grained layers.
*/
template<size_t dim>
double64 csmp::HeterogeneityAndRateAwareModel<dim>::SwStarKinkCO2( double64 sw, double64 Nc ) const
 {
    // TODO: read from model
    const double64  sw_star_kink_CL_(0.7891), sw_star_kink_VL_(0.6805), b_kink_(8000.);
   
    return (sw_star_kink_CL_ - sw_star_kink_VL_) * exp(- b_kink_ * Nc) + SwStar(sw) * sw_star_kink_VL_;
 }



/**
*/
template<size_t dim>
double64 csmp::HeterogeneityAndRateAwareModel<dim>::SwStarKinkH2O( double64 sw, double64 Nc ) const
 {
    // TODO: read from model
    const double64  sw_star_kink_CL_(UNSPECIFIED), sw_star_kink_VL_(UNSPECIFIED), b_kink_(UNSPECIFIED);

    return (sw_star_kink_CL_ - sw_star_kink_VL_) * exp(- b_kink_ * Nc) + SwStar(sw) * sw_star_kink_VL_;
 }



/// standard form: Nc = vt mu_CO2 / sigma
template<size_t dim>
double64 csmp::HeterogeneityAndRateAwareModel<dim>::Nc_vtmuCO2_Version() const
 {
     return (vt_magnitude_ * TwoPhaseModel<dim>::mun_) / IFT_;
 }



/// pressure gradient form: Nc = k ||grad p|| / sigma
template<size_t dim>
double64 csmp::HeterogeneityAndRateAwareModel<dim>::Nc_kgradP_Version() const
 {
     double64 grad_p_magnitude(UNSPECIFIED);
   
     return k_ * grad_p_magnitude / IFT_;
 }


/**
    Maartje, slide 6
*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krw_Phase() const
 {
    //  return linear relperm model if lambda = 0
    if ( lambda_ == static_cast<double64>(0.) ) return TwoPhaseModel<dim>::seff_;
   
    double64 sw_star = SwStar(sw_);
    double64 Nc      = Nc_vtmuCO2_Version();
   
    // if we are at a pc above the entry pressure of the low permeability layer
    if ( sw_star < SwStarKinkH2O( sw_, Nc ) )
      {
         double64 krw_end_CL(0.7); // at sw* = 1.
         double64 a = -(1. - krw_end_CL) * exp(-1500. * Nc) + 1.;
         double64 b(6.), c(0.);
         return c + a * pow( sw_star, b );
      }
   
    // for a lowr capillary pressure
    else {
         double64 krw_begin_CL(0.1222); // at sw* = 0.
         double64 c = krw_begin_CL * exp(-37670. * Nc);
         double64 a(1. - c);
         double64 b(6. + 6. * exp(-20450. * Nc));
         // TODO: is this the right form?
         return c + a * pow( sw_star, b );
      }
   
 } // end krw_Phase





/**
*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krn_Phase() const
 {
    //  switch to linear relperm model if lambda = 0
    if ( lambda_ == static_cast<double64>(0.) ) return static_cast<double64>( 1. - TwoPhaseModel<dim>::seff_ );

    double64 sw_star = SwStar(sw_);
    double64 Nc      = Nc_vtmuCO2_Version();
   
    // if we are at a pc above the entry pressure of the low permeability layer
    if ( sw_star < SwStarKinkCO2( sw_, Nc ) )
      {
         double64 a = -0.34 * exp(-8000. * Nc) + 1.25;
         double64 b = -0.79 * exp(-8000. * Nc) + 2.4;
         double64 c(0.);
         return c + a * pow( (1. - sw_star), b );
      }
   
    // for a lowr capillary pressure
    else {
         double64 a = 6.12 * exp(-8000. * Nc) + 1.63;
         double64 b(2.8);
         return a * pow( (1. - sw_star), b );
      }
 }



template class HeterogeneityAndRateAwareModel<1U>;
template class HeterogeneityAndRateAwareModel<2U>;
template class HeterogeneityAndRateAwareModel<3U>;

} // end namespace csmp


























