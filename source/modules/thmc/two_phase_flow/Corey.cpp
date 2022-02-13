#include "Corey.h"

using namespace std;

namespace csmp {


template<uint32_t dim>
Corey<dim>::Corey()
 {
 }
 


template<uint32_t dim>
Corey<dim>::Corey( const PropertyDatabase<dim>& database )
 : expn_key_(database.StorageKey("corey exponent oil")),
   expw_key_(database.StorageKey("corey exponent water")),
   krn_key_(database.StorageKey("relperm endpoint oil")),
   krw_key_(database.StorageKey("relperm endpoint water")),
   lambda_key_(database.StorageKey("corey exponent capillary pressure")),
   pd_key_(database.StorageKey("capillary entry pressure")),
   TwoPhaseModel<dim>(database, "permeability",
                        "viscosity oil", "viscosity water",
                        "density oil", "density water", "saturation water",
                        "residual saturation oil",
                        "residual saturation water" ),
   pc_max_key_( pd_key_),
   default_capillary_pressure_max_(false)
{
   pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}


template<uint32_t dim>
Corey<dim>::Corey( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               const char* viscosity_nw, const char* viscosity_w,
                               const char* density_nw, const char* density_w,
                               const char* lambda,
                               const char* pc_entry,
                               const char* exp_nw, const char* exp_w,
                               const char* krnw, const char* krw,
                               const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                               const bool sw_ro_mu_placement)
: TwoPhaseModel<dim>(database, permeability,
                     viscosity_nw, viscosity_w,
                     density_nw, density_w, sat_w,
                     res_sat_nw,
                     res_sat_w,
                     sw_ro_mu_placement ),
   expn_key_(database.StorageKey(exp_nw)),
   expw_key_(database.StorageKey(exp_w)),
   krn_key_(database.StorageKey(krnw)),
   krw_key_(database.StorageKey(krw)),
   lambda_key_(database.StorageKey(lambda)),
   pd_key_(database.StorageKey(pc_entry)),
   pc_max_key_( pd_key_),
   default_capillary_pressure_max_(false)
{
   pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}

template<uint32_t dim>
Corey<dim>::Corey( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               const char* viscosity_nw, const char* viscosity_w,
                               const char* density_nw, const char* density_w,
                               const char* lambda,
                               const char* pc_entry,
                               const char* maximum_pc,
                               const char* exp_nw, const char* exp_w,
                               const char* krnw, const char* krw,
                               const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                               const bool sw_ro_mu_placement)
: TwoPhaseModel<dim>(database, permeability,
                     viscosity_nw, viscosity_w,
                     density_nw, density_w, sat_w,
                     res_sat_nw,
                     res_sat_w,
                     sw_ro_mu_placement ),
   expn_key_(database.StorageKey(exp_nw)),
   expw_key_(database.StorageKey(exp_w)),
   krn_key_(database.StorageKey(krnw)),
   krw_key_(database.StorageKey(krw)),
   lambda_key_(database.StorageKey(lambda)),
   pd_key_(database.StorageKey(pc_entry)),
   pc_max_key_((maximum_pc==NULL)? Index(): database.StorageKey(maximum_pc)),
   default_capillary_pressure_max_( (maximum_pc==NULL)? true : false)
 {
    pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
 }


template<uint32_t dim>
Corey<dim>::~Corey()
 {
 }
 
 
 



template<uint32_t dim>
void Corey<dim>::Initialize( const Element<dim>& e )
 {
    TwoPhaseModel<dim>::swr_   = e.Read( TwoPhaseModel<dim>::swr_key_ ); // irreducible water
    TwoPhaseModel<dim>::snr_   = e.Read( TwoPhaseModel<dim>::snr_key_ ); // residual oil

    if( TwoPhaseModel<dim>::tensor_permeability_){
        e.Read( TwoPhaseModel<dim>::perm_key_, TwoPhaseModel<dim>::K_);
        TwoPhaseModel<dim>::k_ = TwoPhaseModel<dim>::K_.Trace()/static_cast<double>(dim);
    }else{
        TwoPhaseModel<dim>::k_ = e.Read( TwoPhaseModel<dim>::perm_key_ );
        TwoPhaseModel<dim>::K_.operator=( VectorVariable<dim>(PLAIN, TwoPhaseModel<dim>::k_ ) );
    }

    entry_pressure_    = e.Read( pd_key_ );
    lambda_ = e.Read( lambda_key_ );

    if( lambda_ == 0.0 )
    {
        if(!default_capillary_pressure_max_)
            pc_max_ = e.Read( pc_max_key_ );
        else
            pc_max_ = TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
    }

    expw_  = e.Read( expw_key_ );
    expn_  = e.Read( expn_key_ );
    krw_   = e.Read( krw_key_ );
    krn_   = e.Read( krn_key_ );


 } // end Initialize




template<uint32_t dim>
double Corey<dim>::krw_Phase() const
 { 
 
    if ( TwoPhaseModel<dim>::seff_ <= 0. ) return static_cast<double>(0.);
    if ( TwoPhaseModel<dim>::seff_ >= 1. ) return krw_;

    return krw_ * std::pow( TwoPhaseModel<dim>::seff_, expw_ );
 }



template<uint32_t dim>
double Corey<dim>::krn_Phase() const
 { 
    if ( TwoPhaseModel<dim>::seff_ <= 0. ) return krn_;
    if ( TwoPhaseModel<dim>::seff_ >= 1. ) return static_cast<double>(0.);
   	   	
    return krn_ * std::pow( 1.0 - TwoPhaseModel<dim>::seff_, expn_ ) ;
 }

template<uint32_t dim>
double Corey<dim>::dkrwds_Phase() const
 {
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    return krw_ * expw_* std::pow( TwoPhaseModel<dim>::seff_, expw_ - 1.0)*seff_mult;
 }



template<uint32_t dim>
double Corey<dim>::dkrnds_Phase() const
 {
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    return -krn_ * expn_ * std::pow( 1.0 - TwoPhaseModel<dim>::seff_ , expn_ -1.0 )*seff_mult ;
 }

/// pc covers the full saturation range, pc is capped based on maximum dpcds
template<uint32_t dim>
double Corey<dim>::pc_Phase( ) const
{
    // linear relperm model
    if ( lambda_ == 0. ){

        if ( TwoPhaseModel<dim>::seff_ <= 0. )
            return pc_max_;

        if ( TwoPhaseModel<dim>::seff_ >= 1. )
            return entry_pressure_;

        return entry_pressure_ + ( 1. - TwoPhaseModel<dim>::seff_ ) * ( TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ - entry_pressure_ );
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
double Corey<dim>::dpcds_Phase( ) const
{

    // linear relperm model
     if ( lambda_ == 0. ){

         if ( ( TwoPhaseModel<dim>::seff_ < 0. ) || ( TwoPhaseModel<dim>::seff_ > 1. ) )
             return 0.;

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
double Corey<dim>::Sw_Phase( double pc ) const
{

     // linear relperm model
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
double Corey<dim>::dsdpc_Phase( double pc ) const
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

    return -(lambda_/entry_pressure_) * pow( pc/entry_pressure_, -lambda_ - 1.0 )/seff_mult;

}


template<uint32_t dim>
double Corey<dim>::MaxFractionalFlowDerivative() const
{
    // this is NOT the maximum of dfdS !!
    // since this function is used to compute the CFL criterion, which
    // needs dfdS for the Welge tangent (at the front saturation Swf) and
    // NOT the maximum of dfdS to compute the shock speed, 1 is returned

    return 1.0;
}


template<uint32_t dim>
void Corey<dim>::Out( size_t phase ) const
 {
    TwoPhaseModel<dim>::Out(phase);
    cout << "\nCorey<csp_float," << dim << ">::Out: Additional properties: " << endl;
    cout << "\nElement properties:";
    cout << "\nCapillary entry pressure (Pa): " << entry_pressure_;
    cout << "\nCorey parameter brine: " << expw_;
    cout << "\nCorey parameter oil: " << expn_;
    cout << "\nRelative permeability endpoint brine: " << krw_;
    cout << "\nRelative permeability endpoint oil: " << krn_;
    cout << "\nCorey exponent capillary pressure: " << lambda_;
	cout << endl << endl;

 } // end Out
 
 

template class Corey<1U>;
template class Corey<2U>;
template class Corey<3U>;

} // end namespace csp


























