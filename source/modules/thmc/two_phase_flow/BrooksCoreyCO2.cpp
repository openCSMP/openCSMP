#include "BrooksCoreyCO2.h"
#include "BrooksCoreyFrontVelocity.h"
#include "PropertyDatabase.h"

using namespace std;

namespace csmp {


template<uint32_t dim>
BrooksCoreyCO2<dim>::BrooksCoreyCO2()
:default_capillary_pressure_max_(false)
{
 pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}

/** 
    One also has to initialize the material parameter key in the base class
    in this case the material parameter is used as the lambda parameter
    of the Brooks-Corey model
*/
template<uint32_t dim>
BrooksCoreyCO2<dim>::BrooksCoreyCO2( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               double viscosity_nw, double viscosity_w,
                               double density_nw, double density_w, 
                               const char* lamda, const char* pc_entry,
                               const bool sw_ro_mu_placement )

 : TwoPhaseModel<dim>(database, viscosity_nw, viscosity_w, density_nw, density_w,
                      permeability, "saturation aqueous phase", 
                     "residual saturation carbonic phase",
                     "residual saturation aqueous phase",
                      sw_ro_mu_placement),

   pd_key_(database.StorageKey(pc_entry)),
   lamda_key_(database.StorageKey(lamda)),
   pc_max_key_(pd_key_),
   default_capillary_pressure_max_(false)
{
  pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}


template<uint32_t dim>
BrooksCoreyCO2<dim>::BrooksCoreyCO2( const PropertyDatabase<dim>& database,
                               const char* lamda, const char* pc_entry,
                               const bool sw_ro_mu_placement )
                               
 : TwoPhaseModel<dim>(database, "permeability",  
                     "viscosity carbonic phase", "viscosity aqueous phase",
                     "density carbonic phase", "density aqueous phase", "saturation aqueous phase",
                     "residual saturation carbonic phase",
                     "residual saturation aqueous phase",
                      sw_ro_mu_placement ),

   pd_key_(database.StorageKey(pc_entry)),
   lamda_key_(database.StorageKey(lamda)),
   pc_max_key_(pd_key_),
   default_capillary_pressure_max_(false)
{
  pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}

template<uint32_t dim>
BrooksCoreyCO2<dim>::BrooksCoreyCO2( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               const char* viscosity_nw, const char* viscosity_w,
                               const char* density_nw, const char* density_w,
                               const char* lamda,
                               const char* pc_entry,
                               const char* sat_w, const char* res_sat_nw,
                               const char* res_sat_w,
                               const bool sw_ro_mu_placement )

 : TwoPhaseModel<dim>(database, permeability,
                      viscosity_nw, viscosity_w,
                      density_nw, density_w, sat_w,
                      res_sat_nw,
                      res_sat_w,
                      sw_ro_mu_placement ),

   pd_key_(database.StorageKey(pc_entry)),
   lamda_key_(database.StorageKey(lamda)),
   pc_max_key_(pd_key_),
   default_capillary_pressure_max_(false)
{
  pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}

template<uint32_t dim>
BrooksCoreyCO2<dim>::BrooksCoreyCO2( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               const char* viscosity_nw, const char* viscosity_w,
                               const char* density_nw, const char* density_w,
                               const char* lamda,
                               const char* pc_entry,
                               const char* maximum_pc,
                               const char* sat_w, const char* res_sat_nw,
                               const char* res_sat_w,
                               const bool sw_ro_mu_placement )

 : TwoPhaseModel<dim>(database, permeability,
                      viscosity_nw, viscosity_w,
                      density_nw, density_w, sat_w,
                      res_sat_nw,
                      res_sat_w,
                      sw_ro_mu_placement ),

   pd_key_(database.StorageKey(pc_entry)),
   lamda_key_(database.StorageKey(lamda)),
   pc_max_key_((maximum_pc==NULL)? Index(): database.StorageKey(maximum_pc)),
   default_capillary_pressure_max_( (maximum_pc==NULL)? true : false)
{
  pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}

template<uint32_t dim>
BrooksCoreyCO2<dim>::~BrooksCoreyCO2()
 {
 }
 
 
 

/// Reading parameters for the linear & BC relperm models
template<uint32_t dim>
void BrooksCoreyCO2<dim>::Initialize( const Element<dim>& e )
 {
    TwoPhaseModel<dim>::swr_ = e.Read( TwoPhaseModel<dim>::swr_key_ );
    TwoPhaseModel<dim>::snr_ = e.Read( TwoPhaseModel<dim>::snr_key_ );

    if( TwoPhaseModel<dim>::tensor_permeability_){
        e.Read( TwoPhaseModel<dim>::perm_key_, TwoPhaseModel<dim>::K_);
        TwoPhaseModel<dim>::k_ = TwoPhaseModel<dim>::K_.Trace()/static_cast<double>(dim);
    }else{
        TwoPhaseModel<dim>::k_ = e.Read( TwoPhaseModel<dim>::perm_key_ );
        TwoPhaseModel<dim>::K_.operator=( VectorVariable<dim>(PLAIN, TwoPhaseModel<dim>::k_ ) );
    }

    entry_pressure_ = e.Read( pd_key_ );
    lambda_ = e.Read( lamda_key_ );

    if( lambda_ == 0.0 )
    {
        if(!default_capillary_pressure_max_)
            pc_max_ = e.Read( pc_max_key_ );
        else
            pc_max_ = TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
    }

    // if properties are discretized on element
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
        TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
        TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
        TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
        TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
        TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );     
    }

 } // end Initialize




template<uint32_t dim>
double BrooksCoreyCO2<dim>::krw_Phase() const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( lambda_ == static_cast<double>(0.) )
        return TwoPhaseModel<dim>::seff_;

    // pm2 = lambda, the Brooks-Corey parameter
    return std::pow( TwoPhaseModel<dim>::seff_, 2./lambda_ + 3.0 );
 }


/// as in Helmig 97, eqn.2.57, p. 75
template<uint32_t dim>
double BrooksCoreyCO2<dim>::krn_Phase() const
 { 
    //  switch to linear relperm model if lambda = 0
    if ( lambda_ == static_cast<double>(0.) )
        return static_cast<double>( 1. - TwoPhaseModel<dim>::seff_ );

    const double  seffn(1. - TwoPhaseModel<dim>::seff_);
    return (seffn * seffn) * (1. - pow( this->seff_, 2./lambda_ + 1.0) );
 }

template<uint32_t dim>
double BrooksCoreyCO2<dim>::dkrwds_Phase() const
 {
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    //  switch to linear relperm model if lambda = 0
    if ( lambda_ == 0.0 )
        return seff_mult;

    return ( 2./lambda_ + 3.0 )*std::pow( TwoPhaseModel<dim>::seff_, 2./lambda_ + 2.0 )*seff_mult;
 }


template<uint32_t dim>
double BrooksCoreyCO2<dim>::dkrnds_Phase() const
 {
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    //  switch to linear relperm model if lambda = 0
    if ( lambda_ == 0. )
        return -seff_mult;

    return 2.0*(TwoPhaseModel<dim>::seff_ - 1.0)*
            ( 1.0 + pow(TwoPhaseModel<dim>::seff_, 2.0/lambda_)*( 1.0/lambda_ + 1./2. - TwoPhaseModel<dim>::seff_*(1.0/lambda_ + 3./2.)))*
            seff_mult;

 }


/// pc covers the full saturation range, pc is capped based on maximum dpcds 
template<uint32_t dim>
double BrooksCoreyCO2<dim>::pc_Phase( ) const
{
    // linear relperm model
    if ( lambda_ == 0. ){

        if ( TwoPhaseModel<dim>::seff_ <= 0. )
            return pc_max_;

        if ( TwoPhaseModel<dim>::seff_ >= 1. )
            return entry_pressure_;

		if (entry_pressure_ == pc_max_)
			return entry_pressure_;

        return entry_pressure_ + ( 1. - TwoPhaseModel<dim>::seff_ ) * ( pc_max_ - entry_pressure_ );
   }

   // for zero entry pressure capillary pressure always is zero
   if ( entry_pressure_ == 0. )
       return 0.;

   const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double pcmax = entry_pressure_ * pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                           ( -1. / ( 1. + lambda_ ) ) );

   if ( TwoPhaseModel<dim>::seff_ <= std::pow( pcmax / entry_pressure_, -lambda_ ) )
   {
      // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
       const double Se_min =  pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                ( lambda_ / ( 1. + lambda_ ) ) );
       // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
       return pcmax + ( Se_min - TwoPhaseModel<dim>::seff_ ) * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult;
   }

   return entry_pressure_ * std::pow( TwoPhaseModel<dim>::seff_, -1. / lambda_ );

}





/// dpcdS covers the full saturation range, dpcds is capped based on maximum dpcds
template<uint32_t dim>
double BrooksCoreyCO2<dim>::dpcds_Phase( ) const
{

    // linear relperm model
     if ( lambda_ == 0. ){

         if ( ( TwoPhaseModel<dim>::seff_ < 0. ) || ( TwoPhaseModel<dim>::seff_ > 1. ) )
             return 0.;

		 if (entry_pressure_ == pc_max_)
			 return 0.0;

         const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

         return -( pc_max_ - entry_pressure_ )*seff_mult;

    }

    // for zero entry pressure dpcds is zero
    if ( entry_pressure_ == 0. )
       return 0.;

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double Se_min = pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                            ( lambda_ / ( 1. + lambda_ ) ) );

    // below Se_min, capillary pressure's derivative is constant and equal to MAXIMUM_DPCDS
    if ( TwoPhaseModel<dim>::seff_ <= Se_min )
        return -TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_;

    // computing the capillary pressure derivative for seff > Se_min
    return -entry_pressure_ * std::pow( TwoPhaseModel<dim>::seff_, -1. - 1. / lambda_ ) /
           lambda_* seff_mult;

}


template<uint32_t dim>
double BrooksCoreyCO2<dim>::Sw_Phase( double pc ) const
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

   const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

   // compute maximum capillary pressure based on maximum dpcds of MAXIMUM_DPCDS
   // applying the limit on capillary pressure
   double pcmax = entry_pressure_ * pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                           ( -1. / ( 1. + lambda_ ) ) );

   if ( pc >= pcmax )
   {
      // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
       const double Se_min =  pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                ( lambda_ / ( 1. + lambda_ ) ) );
       // assuming linear changes in capillary pressure below Se_min with slope of MAXIMUM_DPCDS
       TwoPhaseModel<dim>::seff_ = Se_min - ( pc - pcmax )/TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_*seff_mult ;

       return TwoPhaseModel<dim>::SeffToSw();
   }

   TwoPhaseModel<dim>::seff_ =  std::pow(pc/entry_pressure_, -lambda_);

   return TwoPhaseModel<dim>::SeffToSw();

}

template<uint32_t dim>
double BrooksCoreyCO2<dim>::dsdpc_Phase( double pc ) const
{

    // linear relperm model
    if( lambda_ == 0.0 )
    {
        // not unique solution
        if ( pc_max_ == entry_pressure_ )
            return 0.0;

        const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
        return -1.0/( pc_max_ - entry_pressure_ )/seff_mult;

    }

    if ( entry_pressure_ == 0. )
        return 0.0;

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    // compute minimun effective saturation for which dpcds = MAXIMUM_DPCDS
    const double Se_min =  pow( ( entry_pressure_ / ( lambda_ * TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_/seff_mult ) ),
                                ( lambda_ / ( 1. + lambda_ ) ) );

    const double Seff = (TwoPhaseModel<dim>::Sw_at( pc )- TwoPhaseModel<dim>::swr_ ) *seff_mult;

    if( Seff<= Se_min )
        return -1.0/TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_;

    return -lambda_/entry_pressure_ * pow(pc/entry_pressure_, -lambda_ - 1.0 )/seff_mult;

}


/**

   Calculates the derivative of water fractional flow respect to water saturation for 
   Brooks-Corey and Linear relative permeability models. The derivative is conducted by
   analytical differentiation.

*/

template<uint32_t dim>
double BrooksCoreyCO2<dim>::dfds() const
{

   const double seff = TwoPhaseModel<dim>::seff_;
   
   //if ( ( seff <= static_cast<double>(0.) ) || ( seff >= static_cast<double>(1.) ) )
   //    return static_cast<double>(0.);

   const double lambda = lambda_;
   const double mob_w = krw_Phase() / TwoPhaseModel<dim>::muw_;
   const double mob_n = krn_Phase() / TwoPhaseModel<dim>::mun_;

   // linear model
   if ( lambda == static_cast<double>(0.) )
   {
       return 1. / ( ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) * TwoPhaseModel<dim>::muw_ * TwoPhaseModel<dim>::mun_
              * ( mob_w + mob_n ) * ( mob_w + mob_n ) );
   }
   
   // Brooks-Corey
   const double dmob_wdsw = ( 3. + 2. / lambda ) * std::pow ( seff , 2. + 2. / lambda ) / 
                             ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) / TwoPhaseModel<dim>::muw_;
   const double dmob_ndsw = - ( 2 * ( 1 - seff ) * ( 1 - std::pow ( seff , 1. + 2. / lambda ) ) + ( 1 - seff ) * ( 1 - seff ) 
                             * ( 1. + 2. / lambda ) * std::pow ( seff ,2. / lambda ) ) / 
                             ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) / TwoPhaseModel<dim>::mun_;
 
   return ( dmob_wdsw *  mob_n - dmob_ndsw * mob_w ) / ( ( mob_w + mob_n ) * ( mob_w + mob_n ) );
	
}  // end dfdS


template<uint32_t dim>
double BrooksCoreyCO2<dim>::dGds() const
{

   const double seff = TwoPhaseModel<dim>::seff_;
   
   //if ( ( seff <= static_cast<double>(0.) ) || ( seff >= static_cast<double>(1.) ) )
   //    return static_cast<double>(0.);

   const double lambda = lambda_;

   const double mob_w = krw_Phase() / TwoPhaseModel<dim>::muw_;
   const double mob_n = krn_Phase() / TwoPhaseModel<dim>::mun_;
   
   // linear model
   if ( lambda == static_cast<double>(0.) )
   {
       return ( ( mob_n * mob_n / TwoPhaseModel<dim>::muw_ ) - ( mob_w * mob_w / TwoPhaseModel<dim>::mun_ ) ) / 
              ( ( mob_w + mob_n ) * ( mob_w + mob_n ) ) / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ );
   }
   
   // Brooks-Corey 
   const double dmob_wdsw = ( 3. + 2. / lambda ) * std::pow ( seff , 2. + 2. / lambda ) / 
                             ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) / TwoPhaseModel<dim>::muw_;
   const double dmob_ndsw = - ( 2 * ( 1 - seff ) * ( 1 - std::pow ( seff , 1. + 2. / lambda ) ) + ( 1 - seff ) * ( 1 - seff ) 
                             * ( 1. + 2. / lambda ) * std::pow ( seff ,2. / lambda ) ) / 
                             ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) / TwoPhaseModel<dim>::mun_;
 
   return ( dmob_wdsw * ( mob_n * mob_n ) + dmob_ndsw * ( mob_w * mob_w ) ) / ( ( mob_w + mob_n ) * ( mob_w + mob_n ) );

}  // end dGdS_Phase


/// for the wetting phase
template<uint32_t dim>
double BrooksCoreyCO2<dim>::MaxFractionalFlowDerivative() const
 {
    // linear relperm model
    // in the linear relperm model maximum dfds belongs to seff=0 or seff=1 depending on the viscosity of water and oil
    // if mu_oil>mu_water maximum dfds is dfds@seff=0. and if mu_oil<mu_water maximum dfds is dfds@seff=1.
    if ( lambda_ == static_cast<double>(0.) )
        return ( TwoPhaseModel<dim>::mun_ >= TwoPhaseModel<dim>::muw_ ) ?
               ( TwoPhaseModel<dim>::mun_ / TwoPhaseModel<dim>::muw_ / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) ) :
               ( TwoPhaseModel<dim>::muw_ / TwoPhaseModel<dim>::mun_ / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    return static_cast<double>(5.6); // as computed with dfds method
 }

template<uint32_t dim>
double BrooksCoreyCO2<dim>::ShockSpeed() const
 {
    if ( lambda_ == static_cast<double>(0.) )
    {
        cout <<"\nBrooksCoreyCO2<dim>::ShockHeight: not implemented for linear relative permeability model yet."<< endl;
        return -1.;
    }
   
    return BrooksCoreyFrontVelocity().ShockVelocityMultiplier( lambda_,
                                                               TwoPhaseModel<dim>::mun_/TwoPhaseModel<dim>::muw_ ); 
 }


/// look this one up in the book by Randy LeVeque
template<uint32_t dim>
double BrooksCoreyCO2<dim>::ShockHeight() const
 {

    cout <<"\nBrooksCoreyCO2<dim>::ShockHeight: not implemented yet."<< endl;
    return -1.; 

 }



template<uint32_t dim>
void BrooksCoreyCO2<dim>::Out( uint32_t phase ) const
 {
    TwoPhaseModel<dim>::Out(phase);
    cout <<"\nBrooksCoreyCO2<"<< dim << ">::Out: Additional properties: "<< endl;
    cout <<"\nelement properties:";
    cout <<"\n       capillary entry pressure, pd: "<< entry_pressure_;
    cout <<"\n      Brooks-Corey lambda parameter: "<< lambda_;
    cout <<"\n  MAX_CAPILLARY_PRESSURE_DERIVATIVE: "<< TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_ << endl << endl;

 } // end Out
 
 

template<uint32_t dim>
double csmp::BrooksCoreyCO2<dim>::Pd() const
  {
    return entry_pressure_;
  }


template<uint32_t dim>
void csmp::BrooksCoreyCO2<dim>::Pd( double pd )
  {
    entry_pressure_ = pd;
  }


template<uint32_t dim>
void csmp::BrooksCoreyCO2<dim>::Lambda( double lambda )
  {
    lambda_ = lambda;
  }

template<uint32_t dim>
double csmp::BrooksCoreyCO2<dim>::Lambda() const
  {
    return lambda_;
  }


template class BrooksCoreyCO2<1U>;
template class BrooksCoreyCO2<2U>;
template class BrooksCoreyCO2<3U>;

} // end namespace csp


























