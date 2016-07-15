#include "FractureMatrixUpscaled.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "ErrorFunction.h"

using namespace std;

namespace csmp {

template<size_t dim>
FractureMatrixUpscaled<dim>::FractureMatrixUpscaled()
 {
 }


template<size_t dim>
FractureMatrixUpscaled<dim>::FractureMatrixUpscaled( const PropertyDatabase<dim>& database,
                                                     const char* permeability,
                                                     const char* visc_nw, const char* visc_w,
                                                     const char* rho_nw, const char* rho_w,
                                                     const char* lamda, const char* pc_entry,
                                                     const char* specific_fracture_matrix_interface_area,
                                                     const char* phi_f, // fracture porosity
                                                     const char* phi_m, // matrix~total porosity
                                                     const char* qfqm_ratio,
                                                     const char* volume_flux,
                                                      // for the GTF
                                                     const char* sw_initial,  // new
                                                     const char* block_radius,
                                                     const bool sw_ro_mu_placement )
 : pd_key_(database.StorageKey(pc_entry)),
   lambda_key_(database.StorageKey(lamda)),
   phim_key_(database.StorageKey(phi_m)),
   phif_key_(database.StorageKey(phi_f)),
   Af_key_(database.StorageKey(specific_fracture_matrix_interface_area)),
   qfqm_key_(database.StorageKey(qfqm_ratio)),
   flux_key_(database.StorageKey(volume_flux)),
   swi_key_(database.StorageKey(sw_initial)),
   r_key_(database.StorageKey(block_radius)),
   TwoPhaseModel<dim>(database, permeability,
                         visc_nw, visc_w,
                         rho_nw, rho_w,
                        "saturation water",
                        "residual saturation non-wetting phase",
                        "residual saturation wetting phase",
                         sw_ro_mu_placement)
 {
 }

template<size_t dim>
FractureMatrixUpscaled<dim>::FractureMatrixUpscaled( const PropertyDatabase<dim>& database,
                                                     const char* permeability,
                                                     const char* visc_nw, const char* visc_w,
                                                     const char* rho_nw, const char* rho_w,
                                                     const char* lamda, const char* pc_entry,
                                                     const char* specific_fracture_matrix_interface_area,
                                                     const char* phi_f, // fracture porosity
                                                     const char* phi_m, // matrix~total porosity
                                                     const char* qfqm_ratio,
                                                     const char* volume_flux,
                                                      // for the GTF
                                                     const char* sw_initial,  // new
                                                     const char* block_radius,
                                                     const char* rsnw, const char* rsw,
                                                     const bool sw_ro_mu_placement )
 : pd_key_(database.StorageKey(pc_entry)),
   lambda_key_(database.StorageKey(lamda)),
   phim_key_(database.StorageKey(phi_m)),
   phif_key_(database.StorageKey(phi_f)),
   Af_key_(database.StorageKey(specific_fracture_matrix_interface_area)),
   qfqm_key_(database.StorageKey(qfqm_ratio)),
   flux_key_(database.StorageKey(volume_flux)),
   swi_key_(database.StorageKey(sw_initial)),
   r_key_(database.StorageKey(block_radius)),
   TwoPhaseModel<dim>(database, permeability,
                      visc_nw, visc_w,
                      rho_nw, rho_w,
                      "saturation water",
                      rsnw,
                      rsw,
                      sw_ro_mu_placement)
 {
 }

template<size_t dim>
FractureMatrixUpscaled<dim>::FractureMatrixUpscaled( double permeability,
                                                     double visc_nw, double visc_w,
                                                     double rho_nw, double rho_w,
                                                     double lambda, double pd,
                                                     double Af,
                                                     double phif, 
                                                     double phim, 
                                                     double qfqm,
                                                     double qv,
                                                     double swi,  
                                                     double radius,
                                                     double satwr,
                                                     double satnr )
 : TwoPhaseModel<dim>(),
    Af_( Af ), pd_( pd ), lambda_( lambda ), phim_ ( phim ), phif_( phif ), qfqm_( qfqm ), qv_( qv ), radius_( radius ),
    swi_( PLAIN, swi )
 {
    TwoPhaseModel<dim>::k_     = permeability;
    TwoPhaseModel<dim>::mun_   = visc_nw;
    TwoPhaseModel<dim>::muw_   = visc_w;
    TwoPhaseModel<dim>::rhn_   = rho_nw;
    TwoPhaseModel<dim>::rhw_   = rho_w;
    TwoPhaseModel<dim>::swr_   = satwr;
    TwoPhaseModel<dim>::snr_   = satnr;
 }



template<size_t dim>
FractureMatrixUpscaled<dim>::~FractureMatrixUpscaled()
 {
 }

template<size_t dim>
void FractureMatrixUpscaled<dim>::InitializeForSw( double sw )
 {
   TwoPhaseModel<dim>::sat_ = sw;
 }

/// @note only the element properties
template<size_t dim>
void FractureMatrixUpscaled<dim>::Initialize( const Element<dim>& e )
 {
    if( TwoPhaseModel<dim>::tensor_permeability_){
        e.Read( TwoPhaseModel<dim>::perm_key_, TwoPhaseModel<dim>::K_);
        TwoPhaseModel<dim>::k_ = TwoPhaseModel<dim>::K_.Trace()/static_cast<double64>(dim);
    }else{
        TwoPhaseModel<dim>::k_ = e.Read( TwoPhaseModel<dim>::perm_key_ );
        TwoPhaseModel<dim>::K_.operator=( VectorVariable<dim>(PLAIN, TwoPhaseModel<dim>::k_ ) );
    }

    TwoPhaseModel<dim>::swr_ = e.Read( TwoPhaseModel<dim>::swr_key_ );
    TwoPhaseModel<dim>::snr_ = e.Read( TwoPhaseModel<dim>::snr_key_ );
    // the scalar local flow velocity
    qv_                        = e.Read( flux_key_ );
    // relative permeability parameters (Brooks-Corey)
    pd_                        = e.Read( pd_key_ );
    lambda_                    = e.Read( lambda_key_ );
    assert( lambda_ > 0. );
    phim_                      = e.Read( phim_key_ );
    // fracture matrix parameters
    Af_                        = e.Read( Af_key_ );
    phif_                      = e.Read( phif_key_ );
    qfqm_                      = e.Read( qfqm_key_ );
    radius_                    = e.Read( r_key_ );
    // as interpolated from nodes
    e.PropertyValueAtBaryCenter( swi_key_, swi_ );

    // if properties are discretized on element
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
        TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
        TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
        TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
        TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
        TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );     
    }

 } // end Initialize




/// @note node properties
template<size_t dim>
void FractureMatrixUpscaled<dim>::InitializeForNode( const Element<dim>& e,
                                                        size_t fem_node )
 {
    this->sat_  = e. N(fem_node)->Read( this->sat_key_ );
    swi_       = e. N(fem_node)->Read( swi_key_ );
    this->mun_  = e. N(fem_node)->Read( this->mun_key_ );
    this->muw_  = e. N(fem_node)->Read( this->muw_key_ );
    this->rhn_  = e. N(fem_node)->Read( this->rhn_key_ );
    this->rhw_  = e. N(fem_node)->Read( this->rhw_key_ );
 }




/// wetting phase relative permeability
/// @note these functions are all based on the absolute saturation
template<size_t dim>
double64 FractureMatrixUpscaled<dim>::krw_Phase() const
 {
    // 1. Brooks-Corey based upscaled model using the GTF without boost factor
    GenericTransferFunction  gtf( phim_, pd_, lambda_, this->muw_, this->mun_ );
    // matrix relperm
    const double64 krw_m = gtf.krw_BC( this->seff_ );
    // fracture relperm from SKM's new formulation
    const double64 krw_f = (erf(this->sat_/phif_) * this->sat_ * qfqm_) / ((1.-this->sat_) / qfqm_ + this->sat_ * qfqm_);

    // krw weighted by qfqm ratio
    return ( qfqm_ > 1. ) ? krw_f : krw_m;
  }





/// nonwetting phase relative permeability with transfer term
template<size_t dim>
double64 FractureMatrixUpscaled<dim>::krn_Phase() const
 {

    // 1. upscaled relative permeability
    GenericTransferFunction  gtf( phim_, pd_, lambda_, this->muw_, this->mun_ );
    // cout <<"\ngtf "<< gtf.Transfer( this->sat, swi_.Value(), swi_.Value(), this->k, radius_ );

    // matrix relperm
    const double64 krn_m = gtf.krn_BC( this->seff_ );

    // fracture relperm
    double64 krn_f = (erfc(this->sat_/phif_) * (1.-this->sat_)) / (qfqm_ * ((1.-this->sat_)/qfqm_ + this->sat_ * qfqm_));

    // contribution to kro due to capillary transfer
    krn_f += Af_sw() * gtf.Transfer( this->sat_, swi_.Value(), swi_.Value(), this->k_, radius_ ) / (1. + qv_ * qfqm_);

    // 2. weighting krn by qfqm ratio and limiting it so that sum of relperms is <= 1.
    return ( qfqm_ > 1. ) ? std::max( 1. - krw_Phase(), krn_f ) : krn_m;
  }

template<size_t dim>
double64 FractureMatrixUpscaled<dim>::pc_Phase( ) const
{
   GenericTransferFunction  gtf( phim_, pd_, lambda_, this->muw_, this->mun_ );
   return gtf.Pc( this->sat_ );
}


/// for the wetting phase
template<size_t dim>
double64 FractureMatrixUpscaled<dim>::MaxFractionalFlowDerivative() const
 {
    return static_cast<double64>(5.6); // as computed with dfds method
 }

template<size_t dim>
double64 FractureMatrixUpscaled<dim>::dfds() const
{

    return TwoPhaseModel<dim>::dfds_numerical();

}  // end dfdS_Phase



template<size_t dim>
double64 FractureMatrixUpscaled<dim>::dGds( ) const
{

    return TwoPhaseModel<dim>::dGds_numerical();

}  // end dGdS_Phase

//cout <<"\nsw "<< TwoPhaseModel<dim>::sat <<", bcp "<< TwoPhaseModel<dim>::pm2 <<", dGdS "<< dGdS <<" "; cout.flush();


template<size_t dim>
void FractureMatrixUpscaled<dim>::Out( size_t phase ) const
 {
    TwoPhaseModel<dim>::Out(phase);
    cout <<"\nFractureMatrixUpscaled<"<< dim << ">::Out: Additional properties: "<< endl;
    cout <<"\nelement properties:";
    cout <<"\n       capillary entry pressure, pd: "<< pd_;
    cout <<"\n      Brooks-Corey lambda parameter: "<< lambda_;
    cout <<"\n             MAX_CAPILLARY_PRESSURE: "<< TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ << endl << endl;

 } // end Out




template class FractureMatrixUpscaled<1U>;
template class FractureMatrixUpscaled<2U>;
template class FractureMatrixUpscaled<3U>;

} // end namespace csmp

