#include "LinearTwoPhaseModel.h"

using namespace std;

namespace csmp {


// LINEAR RELATIVE PERMEABILITY MODEL
template<size_t dim>
LinearTwoPhaseModel<dim>::LinearTwoPhaseModel(const PropertyDatabase<dim>& database,
                                              const char* permeability, 
                                              double viscosity_nw, double viscosity_w,
                                              double density_nw, double density_w,
                                              const char* pc_entry, const char* sw, const char* rsnw, const char* rsw,
                                              const bool sw_ro_mu_placement )
  :TwoPhaseModel<dim>(database, viscosity_nw, viscosity_w, density_nw, density_w,
   permeability, sw,   rsnw,  rsw,  sw_ro_mu_placement),
   pd_key_(database.StorageKey(pc_entry)),
   pc_max_key_(database.StorageKey(pc_entry)),
   default_capillary_pressure_max_(false),
   pc_max_( TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_)
  {
  }

template<size_t dim>
LinearTwoPhaseModel<dim>::LinearTwoPhaseModel(const PropertyDatabase<dim>& database,
                                         const char* permeability,
                                         double viscosity_nw, double viscosity_w,
                                         double density_nw, double density_w,
                                         const char* pc_entry,
                                         const bool sw_ro_mu_placement)
 :TwoPhaseModel<dim>(database, viscosity_nw, viscosity_w, density_nw, density_w,
                      permeability, "saturation water",
                     "residual saturation non-wetting phase",
                     "residual saturation wetting phase",
                      sw_ro_mu_placement),
  pd_key_(database.StorageKey(pc_entry)),
  pc_max_key_(database.StorageKey(pc_entry)),
  default_capillary_pressure_max_( false )
 {
    pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
 }


template<size_t dim>
LinearTwoPhaseModel<dim>::LinearTwoPhaseModel( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               const char* viscosity_nw, const char* viscosity_w,
                               const char* density_nw, const char* density_w,
                               const char* pc_entry,
                               const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                               const bool sw_ro_mu_placement )

 : TwoPhaseModel<dim>(database, permeability,
                      viscosity_nw, viscosity_w,
                      density_nw, density_w,
                      sat_w, res_sat_nw, res_sat_w,
                      sw_ro_mu_placement ),
   pd_key_(database.StorageKey(pc_entry)),
   pc_max_key_(database.StorageKey(pc_entry)),
   default_capillary_pressure_max_(false)
 {
    pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
 }

// LINEAR RELATIVE PERMEABILITY MODEL
template<size_t dim>
LinearTwoPhaseModel<dim>::LinearTwoPhaseModel(const PropertyDatabase<dim>& database,
                                         const char* permeability,
                                         double viscosity_nw, double viscosity_w,
                                         double density_nw, double density_w,
                                         const char* pc_entry,
                                         const char* maximum_pc,
                                         const bool sw_ro_mu_placement)
 :TwoPhaseModel<dim>(database, viscosity_nw, viscosity_w, density_nw, density_w,
                      permeability, "saturation water",
                     "residual saturation non-wetting phase",
                     "residual saturation wetting phase",
                      sw_ro_mu_placement),
  pd_key_(database.StorageKey(pc_entry)),
  pc_max_key_((maximum_pc==NULL)? Index(): database.StorageKey(maximum_pc)),
  default_capillary_pressure_max_( (maximum_pc==NULL)? true : false)
{
    pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
}

template<size_t dim>
LinearTwoPhaseModel<dim>::LinearTwoPhaseModel( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               const char* viscosity_nw, const char* viscosity_w,
                               const char* density_nw, const char* density_w,
                               const char* pc_entry,
                               const char* maximum_pc,
                               const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                               const bool sw_ro_mu_placement )

 : TwoPhaseModel<dim>(database, permeability,
                      viscosity_nw, viscosity_w,
                      density_nw, density_w,
                      sat_w, res_sat_nw, res_sat_w,
                      sw_ro_mu_placement ),
   pd_key_(database.StorageKey(pc_entry)),
   pc_max_key_((maximum_pc==NULL)? Index(): database.StorageKey(maximum_pc)),
   default_capillary_pressure_max_( (maximum_pc==NULL)? true : false)
 {
    pc_max_ =  TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
 }



template<size_t dim> 
LinearTwoPhaseModel<dim>::~LinearTwoPhaseModel()
 {
 }


/// Reading parameters for the linear relperm model
template<size_t dim>
void LinearTwoPhaseModel<dim>::Initialize( const Element<dim>& e )
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

    if(!default_capillary_pressure_max_)
        pc_max_ = e.Read( pc_max_key_ );
    else
        pc_max_ = TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;

    // if properties are discretized on element
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
        TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
        TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
        TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
        TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
        TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );     
    }

 } // end Initialize
 
 
 
template<size_t dim>
double LinearTwoPhaseModel<dim>::krw_Phase() const
 {
    return TwoPhaseModel<dim>::seff_;
 }



template<size_t dim>
double LinearTwoPhaseModel<dim>::krn_Phase() const
 {
    return 1. - TwoPhaseModel<dim>::seff_;
 } 
 
template<size_t dim>
double LinearTwoPhaseModel<dim>::dkrwds_Phase() const
 {
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    return seff_mult;
 }



template<size_t dim>
double LinearTwoPhaseModel<dim>::dkrnds_Phase() const
 {
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    return -seff_mult;
 }


template<size_t dim>
double LinearTwoPhaseModel<dim>::pc_Phase( ) const
{
    if ( TwoPhaseModel<dim>::seff_ <= 0. )
        return pc_max_;

    if ( TwoPhaseModel<dim>::seff_ >= 1. )
        return entry_pressure_;
        
	if (entry_pressure_ == pc_max_)
		return entry_pressure_;

    return entry_pressure_ + ( 1. - TwoPhaseModel<dim>::seff_ ) * ( pc_max_ - entry_pressure_ );

}

template<size_t dim>
double LinearTwoPhaseModel<dim>::dpcds_Phase( ) const
{
    if ( ( TwoPhaseModel<dim>::seff_ < 0. ) || ( TwoPhaseModel<dim>::seff_ > 1. ) )
        return 0.;

	if (entry_pressure_ == pc_max_)
		return 0.0;

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    return -( pc_max_ - entry_pressure_ )*seff_mult;

}

template<size_t dim>
double LinearTwoPhaseModel<dim>::Sw_Phase( double pc ) const
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

template<size_t dim>
double LinearTwoPhaseModel<dim>::dsdpc_Phase( double ) const
{
    // not unique solution
    if ( pc_max_ == entry_pressure_ )
        return 0.0;

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    return -1.0/( pc_max_ - entry_pressure_ )/seff_mult;

}


/// for the wetting phase
template<size_t dim>
double LinearTwoPhaseModel<dim>::MaxFractionalFlowDerivative() const
{
    // linear relperm model
    // in the linear relperm model maximum dfds belongs to seff=0 or seff=1 depending on the viscosity of water and oil
    // if mu_oil>mu_water maximum dfds is dfds@seff=0. and if mu_oil<mu_water maximum dfds is dfds@seff=1.

    return ( TwoPhaseModel<dim>::mun_ >= TwoPhaseModel<dim>::muw_ ) ?
           ( TwoPhaseModel<dim>::mun_ / TwoPhaseModel<dim>::muw_ / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) ) :
           ( TwoPhaseModel<dim>::muw_ / TwoPhaseModel<dim>::mun_ / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

}


template<size_t dim>
double LinearTwoPhaseModel<dim>::dfds() const
 {   
   //if ( ( TwoPhaseModel<dim>::seff_ <= static_cast<double>(0.) ) || ( TwoPhaseModel<dim>::seff_ >= static_cast<double>(1.) ) )
   //    return static_cast<double>(0.);
       
   const double mob_w = krw_Phase() / TwoPhaseModel<dim>::muw_;
   const double mob_n = krn_Phase() / TwoPhaseModel<dim>::mun_;

   return 1. / ( ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) * TwoPhaseModel<dim>::muw_ * TwoPhaseModel<dim>::mun_
              * ( mob_w + mob_n ) * ( mob_w + mob_n ) );
             
 } // end DerivativeOfFractionalFlowOfWettingPhase
 
 
 

template<size_t dim>
double LinearTwoPhaseModel<dim>::dGds( ) const
 {
   //if ( ( TwoPhaseModel<dim>::seff_ <= static_cast<double>(0.) ) || ( TwoPhaseModel<dim>::seff_ >= static_cast<double>(1.) ) )
   //    return static_cast<double>(0.);

   const double mob_w = krw_Phase() / TwoPhaseModel<dim>::muw_;
   const double mob_n = krn_Phase() / TwoPhaseModel<dim>::mun_;
   
   return ( ( mob_n * mob_n / TwoPhaseModel<dim>::muw_ ) - ( mob_w * mob_w / TwoPhaseModel<dim>::mun_ ) ) /
            ( ( mob_w + mob_n ) * ( mob_w + mob_n ) ) / ( 1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ );

 } // end dGds


template class LinearTwoPhaseModel<1U>;
template class LinearTwoPhaseModel<2U>;
template class LinearTwoPhaseModel<3U>;


} // end namespace csmp
