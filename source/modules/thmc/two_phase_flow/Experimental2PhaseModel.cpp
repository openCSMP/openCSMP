#include "Experimental2PhaseModel.h"
#include "PropertyDatabase.h"
#include "FiniteVolumeStencil.h"

using namespace std;

namespace csmp {


template<size_t dim>
Experimental2PhaseModel<dim>::Experimental2PhaseModel()
 : acc_gravity_(9.8066),
   dfds_max_(std::numeric_limits<double64>::quiet_NaN()),
   dfds_shock_(std::numeric_limits<double64>::quiet_NaN()),
   s1_shock_(std::numeric_limits<double64>::quiet_NaN())
 {
    cout <<"\nExperimental2PhaseModel: Error: this default constructor does not suffice to get this model to work.\n";
 }
 
 
 
 
/**

@warning material parameter key in the base class must be initialised
as it is used to distinguish betweeen linear and experimental
relative permeability model (linear=0., experimental=any other)*/
template<size_t dim>
Experimental2PhaseModel<dim>::Experimental2PhaseModel( 
                                  const PropertyDatabase<dim>& database,
                                  const char* permeability,
                                  double64 viscosity_nw, double64 viscosity_w,
                                  double64 density_nw, double64 density_w,
                                  const char* kr1_data_file,
                                  const char* kr2_data_file,
                                  const char* pc_data_file,
                                  const bool sw_ro_mu_placement )
 : acc_gravity_(9.8066),
   dfds_max_(std::numeric_limits<double64>::quiet_NaN()),
   dfds_shock_(std::numeric_limits<double64>::quiet_NaN()),
   s1_shock_(std::numeric_limits<double64>::quiet_NaN()),
   // NB: in this base class the names of the saturation variables must be changed if so desired
   TwoPhaseModel<dim>(database, viscosity_nw, viscosity_w, density_nw, density_w,
                         permeability, "saturation water", 
                        "residual saturation non-wetting phase",
                        "residual saturation wetting phase",
                         sw_ro_mu_placement),
    kr1_(kr1_data_file),
    kr2_(kr2_data_file),
    pc_(pc_data_file),
    model_key_(database.StorageKey("brooks corey parameter")), 
    pd_key_(database.StorageKey("entry pressure")) 
 {
 }




template<size_t dim>
Experimental2PhaseModel<dim>::Experimental2PhaseModel( 
                                  const PropertyDatabase<dim>& database,
                                  const char* kr1_data_file,
                                  const char* kr2_data_file,
                                  const char* pc_data_file,
                                  const bool sw_ro_mu_placement )
 : acc_gravity_(9.8066),
   dfds_max_(std::numeric_limits<double64>::quiet_NaN()),
   dfds_shock_(std::numeric_limits<double64>::quiet_NaN()),
   s1_shock_(std::numeric_limits<double64>::quiet_NaN()),
   // NB: in this base class the names of the saturation variables must be changed if so desired
   TwoPhaseModel<dim>(database, "permeability",  
                        "viscosity oil", "viscosity water",
                        "density oil", "density water", "saturation water",
                        "residual saturation non-wetting phase",
                        "residual saturation wetting phase",
                         sw_ro_mu_placement ),
    kr1_(kr1_data_file),
    kr2_(kr2_data_file),
    pc_(pc_data_file),
    model_key_(database.StorageKey("brooks corey parameter")), 
    pd_key_(database.StorageKey("entry pressure")) 
 {
 }




template<size_t dim>
Experimental2PhaseModel<dim>::~Experimental2PhaseModel()
 {
 }
 
 
 

template<size_t dim>
void Experimental2PhaseModel<dim>::Initialize( const Element<dim>& e )
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
    lambda_ = e.Read( model_key_ );

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
void Experimental2PhaseModel<dim>::ComputeShockSpeedAndHeight() const
 {
    const double64  seff_original(TwoPhaseModel<dim>::seff_);
    double64        dfds_s_max(0.); 
    dfds_shock_ = dfds_max_ = 0.;    
    // loop over the saturation interval finding the maximum value of the fractional flow derivative
    // note the bounds! - only within these dfds is actually defined
    for ( TwoPhaseModel<dim>::seff_=0.005;
          TwoPhaseModel<dim>::seff_<=0.995;
          TwoPhaseModel<dim>::seff_+=0.005 ) {
         double64 dfds1 = Experimental2PhaseModel<dim>::dfds();
         // max fractional flow derivative
         dfds_max_ = std::max( dfds_max_, dfds1 );
         // dfds at shock front and shock height
         double64  dfds_s = dfds1 * TwoPhaseModel<dim>::seff_;
         if ( dfds_s > dfds_s_max ) {
              dfds_s_max = dfds_s;
              dfds_shock_= dfds1;
              s1_shock_ = TwoPhaseModel<dim>::seff_;
         }
      }
   TwoPhaseModel<dim>::seff_ = seff_original;
     
 } // end ComputeShockSpeedAndHeight


/// for the wetting phase
template<size_t dim>
double64 Experimental2PhaseModel<dim>::MaxFractionalFlowDerivative() const
 {

    ComputeShockSpeedAndHeight();
    return dfds_max_;
 }

/// dfds(s1) which will maximise the transport, i.e., the product: vt * dfds * s1
template<size_t dim>
double64 Experimental2PhaseModel<dim>::ShockSpeed() const
 {

    ComputeShockSpeedAndHeight();
    return dfds_shock_;
 }


/// saturation of phase 1 at the shock front
template<size_t dim>
double64 Experimental2PhaseModel<dim>::ShockHeight() const
 {
    ComputeShockSpeedAndHeight();
    return s1_shock_;
 }


template<size_t dim>
double64 Experimental2PhaseModel<dim>::krw_Phase() const
 { 

    if ( TwoPhaseModel<dim>::seff_ <= 0. ) return static_cast<double64>(0.);
    if ( TwoPhaseModel<dim>::seff_ >= 1. ) return static_cast<double64>(1.);

    return kr1_.Value( TwoPhaseModel<dim>::seff_ );
 }



template<size_t dim>
double64 Experimental2PhaseModel<dim>::krn_Phase() const
 { 

    if ( TwoPhaseModel<dim>::seff_ <= 0. ) return static_cast<double64>(1.);
    if ( TwoPhaseModel<dim>::seff_ >= 1. ) return static_cast<double64>(0.);
    
    return kr2_.Value( TwoPhaseModel<dim>::seff_ );
 }

template<size_t dim>
double64 Experimental2PhaseModel<dim>::dkrwds_Phase() const
 {
    const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    return kr1_.Derivative( TwoPhaseModel<dim>::seff_ )*seff_mult;
 }



template<size_t dim>
double64 Experimental2PhaseModel<dim>::dkrnds_Phase() const
 {
    const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
    return kr2_.Derivative( TwoPhaseModel<dim>::seff_ )*seff_mult;
 }

/// pc covers the full saturation range, pc is capped if sw<swr
template<size_t dim>
double64 Experimental2PhaseModel<dim>::pc_Phase( ) const
{
   return pc_.Value( TwoPhaseModel<dim>::seff_ );
}



/// dpcdS covers the full saturation range, dpcdS is capped if sw<swr
template<size_t dim>
double64 Experimental2PhaseModel<dim>::dpcds_Phase( ) const
{
   const double64 seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
   return pc_.Derivative( TwoPhaseModel<dim>::seff_ )*seff_mult;
}
    


 

template<size_t dim>
double64 Experimental2PhaseModel<dim>::dfds() const
{

    return TwoPhaseModel<dim>::dfds_numerical();

}  // end dfdS_Phase




template<size_t dim>
double64 Experimental2PhaseModel<dim>::dGds( ) const
{
    return TwoPhaseModel<dim>::dGds_numerical();

}  // end dGdS_Phase



template<size_t dim>
void Experimental2PhaseModel<dim>::Out( size_t phase ) const
 {
    TwoPhaseModel<dim>::Out(phase);
    cout <<"\nExperimental2PhaseModel<"<< dim << ">::Out: Input data for relperm and pc calculations:"<< endl;
    cout <<"\nrelative permeability curve kr1:";
    kr1_.Out();
    cout <<"\nrelative permeability curve kr2:";
    kr2_.Out();
    cout <<"\ncapillary pressure saturation curve:";
    pc_.Out();

 } // end Out
 
 

template class Experimental2PhaseModel<1U>;
template class Experimental2PhaseModel<2U>;
template class Experimental2PhaseModel<3U>;

} // end namespace csp


























