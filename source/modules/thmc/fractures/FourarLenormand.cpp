#include "FourarLenormand.h"
#include "PropertyDatabase.h"
#include "CSMP_mathUtilities.h"
#include "FiniteVolumeStencil.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
FourarLenormand<dim>::FourarLenormand()
 {
 }


template<uint32_t dim>
FourarLenormand<dim>::FourarLenormand( const PropertyDatabase<dim>& database,
                                       const char* fractureAperture )
  :  TwoPhaseModel<dim>(database, "permeability",
                        "viscosity oil", "viscosity water",
                        "density oil", "density water", "saturation water",
                        "residual saturation non-wetting phase",
                        "residual saturation wetting phase" ),
     fractureApertureKey_( database.StorageKey( fractureAperture ) )
 {
     assert( fractureApertureKey_.place == ELEMENT );
 }

template<uint32_t dim>
FourarLenormand<dim>::FourarLenormand( const PropertyDatabase<dim>& database,
                                       const char* fractureAperture,
                                       const char* permeability,
                                       const char* my_non_wetting,
                                       const char* my_wetting,
                                       const char* rho_non_wetting,
                                       const char* rho_wetting,
                                       const char* sat_wetting,
                                       const char* res_sat_non_wetting,
                                       const char* res_sat_wetting )
  :  TwoPhaseModel<dim>(database, permeability,
                        my_non_wetting, my_wetting,
                        rho_non_wetting, rho_wetting, sat_wetting,
                        res_sat_non_wetting,
                        res_sat_wetting ),
     fractureApertureKey_( database.StorageKey( fractureAperture ) )
 {
     assert( fractureApertureKey_.place == ELEMENT );
 }
 
 
template<uint32_t dim>
FourarLenormand<dim>::~FourarLenormand()
 {
 }



template<uint32_t dim>
void FourarLenormand<dim>::Initialize( const Element<dim>& e )
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

    fractureAperture_       = e.Read( fractureApertureKey_ );

 } // end Initialize



template<uint32_t dim>
double FourarLenormand<dim>::krw_Phase() const
{
  double se = TwoPhaseModel<dim>::seff_;
  if ( se <= 0.) return static_cast<double>(0.);
  if ( se >= 1.) return static_cast<double>(1.);

  return pow( se, 2 )/2 * ( 3-se );
}

template<uint32_t dim>
double FourarLenormand<dim>::krn_Phase() const
{
  double se = TwoPhaseModel<dim>::seff_;
  if ( se <= 0.) return static_cast<double>(1.);
  if ( se >= 1.) return static_cast<double>(0.);

  return pow( ( 1.-se ), 3. ) + 3./2. * TwoPhaseModel<dim>::ViscosityRatio() * se * ( 1.-se ) * ( 2.-se );

}

template<uint32_t dim>
double FourarLenormand<dim>::dkrwds_Phase() const
{
  double se = TwoPhaseModel<dim>::seff_;

  const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
  return 3.0*se*(1.0 - 0.5*se)*seff_mult;

}

template<uint32_t dim>
double FourarLenormand<dim>::dkrnds_Phase() const
{
  double se = TwoPhaseModel<dim>::seff_;
  const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );
  return (-3.0*pow( ( 1.-se ), 2. ) + 3./2. * TwoPhaseModel<dim>::ViscosityRatio() * ( 2.0 - 6.0*se + 3.0*se*se ))*seff_mult;

}

/// uses dummy aperture if key is not specified in according constructor
template<uint32_t dim>
double FourarLenormand<dim>::pc_Phase(  ) const
{
  double se = TwoPhaseModel<dim>::seff_;

  // only for the min saturation of water precautions are needed
  // the actual saturation is used instead of the effective saturation
  if ( se == 0. ) return TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_;
  if ( se == 1. ) return static_cast<double>(0.);

  return 0.001 * ( 4 * 75 * cos( PI * 45. / 180 )  /  fractureAperture_ ); //for single fracture, water-air in [Pa]

  // formulation for DFN REV
  /*
    double se = TwoPhaseModel<dim>::seff;
    double pcEntry = 0.001 * ( 4 * 75 * cos(45)  / fractureAperture_  ); //for single fracture, water-air in [Pa]
    double pcMax = pcEntry * 10;
    double parameter = 2.;

    if( se <= pow( pcMax / pcEntry, -parameter))
      return pcMax;

    return std::min( pcMax, pcEntry * pow( se, -1/parameter));
    */

}

//dummy as above
template<uint32_t dim>
double FourarLenormand<dim>::dpcds_Phase( ) const
{
  double se = TwoPhaseModel<dim>::seff_;

  if ( se == 0. ) return -TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_;
  if ( se == 1. ) return -TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_;

  return 0.0;

  //numerical difference
  //return TwoPhaseModel<dim>::dpcds_numeric();

}

template<uint32_t dim>
double FourarLenormand<dim>::Sw_Phase( double pc ) const
{
    // not unique solution
    return std::min( std::max( TwoPhaseModel<dim>::sat_, TwoPhaseModel<dim>::swr_), 1.0 - TwoPhaseModel<dim>::snr_);
}


template<uint32_t dim>
double FourarLenormand<dim>::dsdpc_Phase( double pc ) const
{
    return 0.0;
}


template class FourarLenormand<1U>;
template class FourarLenormand<2U>;
template class FourarLenormand<3U>;

} // csmp namespace
