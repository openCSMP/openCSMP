#include "CSMP_mathUtilities.h"
#include "VanGenuchten.h"
#include "PropertyDatabase.h"
#include "CSMP_mathUtilities.h"
#include "FiniteVolumeStencil.h"

using namespace std;

namespace csmp {


template<uint32_t dim>
VanGenuchten<dim>::VanGenuchten()
 : acc_gravity_(9.80665), 
   PC_LOW_SW_LIMIT_ (0.01),
   PC_HIGH_SW_LIMIT_(0.99),
   KRW_HIGH_SW_LIMIT_(0.9),
   KRN_LOW_SW_LIMIT_(0.1),
   matrix_tensor_permeability_(false)
 {
 }
 
 
template<uint32_t dim>
VanGenuchten<dim>::VanGenuchten( const PropertyDatabase<dim>& database,
                                 const char* model_parameter_n,
                                 const char* model_parameter_alpha,
                                 bool withFractureMatrixTransfer,
                                 const bool sw_ro_mu_placement )
 : n_key_(database.StorageKey(model_parameter_n)),
   alpha_key_(database.StorageKey(model_parameter_alpha)),
   residualMatrixSaturationWKey_((withFractureMatrixTransfer)? database.StorageKey( "residual saturation water matrix" ) : csmp::Index()),
   residualMatrixSaturationOKey_((withFractureMatrixTransfer)? database.StorageKey( "residual saturation oil matrix" ) : csmp::Index()),
   matrixPermeabilityKey_((withFractureMatrixTransfer)? database.StorageKey( "matrix permeability" ) : csmp::Index()),
   matrix_tensor_permeability_(false),
   acc_gravity_(9.80665),
   PC_LOW_SW_LIMIT_ (0.01),
   PC_HIGH_SW_LIMIT_(0.99),
   KRW_HIGH_SW_LIMIT_(0.9),
   KRN_LOW_SW_LIMIT_(0.1),
   TwoPhaseModel<dim>(database, "permeability",
                        "viscosity oil", "viscosity water",
                        "density oil", "density water", "saturation water",
                        "residual saturation non-wetting phase",
                        "residual saturation wetting phase",
                         sw_ro_mu_placement )
 {

    if(withFractureMatrixTransfer)
        if( matrixPermeabilityKey_.type == TENSOR )
            matrix_tensor_permeability_ = true;

 }

template<uint32_t dim>
VanGenuchten<dim>::VanGenuchten( const PropertyDatabase<dim>& database,
                               const char* permeability,
                               const char* viscosity_nw, const char* viscosity_w,
                               const char* density_nw, const char* density_w,
                               const char* model_parameter_n,
                               const char* model_parameter_alpha,
                               const char* sat_w, const char* res_sat_nw, const char* res_sat_w,
                               const char* permeability_matrix, const char* res_sat_nw_matrix, const char* res_sat_w_matrix,
                               bool withFractureMatrixTransfer,
                               const bool sw_ro_mu_placement )

 : n_key_(database.StorageKey(model_parameter_n)),
   alpha_key_(database.StorageKey(model_parameter_alpha)),
   residualMatrixSaturationWKey_((withFractureMatrixTransfer)? database.StorageKey( res_sat_w_matrix ) : csmp::Index()),
   residualMatrixSaturationOKey_((withFractureMatrixTransfer)? database.StorageKey( res_sat_nw_matrix ) : csmp::Index()),
   matrixPermeabilityKey_((withFractureMatrixTransfer)? database.StorageKey( permeability_matrix ) : csmp::Index()),
   matrix_tensor_permeability_(false),
   acc_gravity_(9.80665),
   PC_LOW_SW_LIMIT_ (0.01),
   PC_HIGH_SW_LIMIT_(0.99),
   KRW_HIGH_SW_LIMIT_(0.9),
   KRN_LOW_SW_LIMIT_(0.1),
   TwoPhaseModel<dim>(database, permeability,
                         viscosity_nw, viscosity_w,
                         density_nw, density_w,
                         sat_w, res_sat_nw, res_sat_w,
                         sw_ro_mu_placement )
{

    if( withFractureMatrixTransfer )
        if( matrixPermeabilityKey_.type == TENSOR )
            matrix_tensor_permeability_ = true;

}
template<uint32_t dim>
VanGenuchten<dim>::~VanGenuchten()
 {
 }



template<uint32_t dim>
void VanGenuchten<dim>::Initialize( const Element<dim>& e )
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

    n_     = e.Read( n_key_ );
    alpha_ = e.Read( alpha_key_ );

    // if properties are discretized on element
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
        TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
        TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
        TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
        TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
        TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );     
    }

 } // end Initialize


/// initializer for DFN model with virtual matrix transfer function
template<uint32_t dim>
void VanGenuchten<dim>::InitializeVirtualMatrix( const Element<dim>& e )
 {

    TwoPhaseModel<dim>::swr_ = e.Read( residualMatrixSaturationWKey_ );
    TwoPhaseModel<dim>::snr_ = e.Read( residualMatrixSaturationOKey_ );

    if( matrix_tensor_permeability_){
        e.Read( matrixPermeabilityKey_, TwoPhaseModel<dim>::K_);
        TwoPhaseModel<dim>::k_ = TwoPhaseModel<dim>::K_.Trace()/static_cast<double>(dim);
    }else{
        TwoPhaseModel<dim>::k_ = e.Read( matrixPermeabilityKey_ );
        TwoPhaseModel<dim>::K_.operator=( VectorVariable<dim>(PLAIN, TwoPhaseModel<dim>::k_ ) );
    }

    n_     = e.Read( n_key_ );
    alpha_ = e.Read( alpha_key_ );

 } // end Initialize


// NB: All these functions assume that the saturation of the wetting phase is used in the computations


/// calculates van Genuchten parameter m = 1 - 1/n
template<uint32_t dim>
double VanGenuchten<dim>::m_from_n( double n ) const
 {
    return 1. - 1. / n;
 }



/// Helmig 97', p. 75 eqn. 2.58
template<uint32_t dim>
double VanGenuchten<dim>::krw_Phase() const
 {
    if ( TwoPhaseModel<dim>::seff_ <= 0.) return static_cast<double>(0.);
    if ( TwoPhaseModel<dim>::seff_ >= 1.) return static_cast<double>(1.);

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    // regularization
    if ( TwoPhaseModel<dim>::seff_ > KRW_HIGH_SW_LIMIT_ ){

        const double y = TwoPhaseModel<dim>::krw_at( KRW_HIGH_SW_LIMIT_ );
        const double k = TwoPhaseModel<dim>::dkrwds_at( KRW_HIGH_SW_LIMIT_)/seff_mult;

        return splineValue( TwoPhaseModel<dim>::seff_, KRW_HIGH_SW_LIMIT_, 1.0 , y, 1.0, k, 0.0);
    }

    //        krw = [ 1 - (1 - Se^1/m)^m ]^2 * Se^eps
    double m = m_from_n(n_);
    double term = 1. - std::pow( 1. - std::pow( TwoPhaseModel<dim>::seff_, 1./m ), m );

    return std::pow( TwoPhaseModel<dim>::seff_, 0.5 ) * (term * term);
 }




/// Helmig 97', p. 75 eqn. 2.59
template<uint32_t dim>
double VanGenuchten<dim>::krn_Phase() const
 {
    if ( TwoPhaseModel<dim>::seff_ <= 0.) return static_cast<double>(1.);
    if ( TwoPhaseModel<dim>::seff_ >= 1.) return static_cast<double>(0.);

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    // regularization
    if ( TwoPhaseModel<dim>::seff_ < KRN_LOW_SW_LIMIT_ ){

        const double y = TwoPhaseModel<dim>::krn_at( KRN_LOW_SW_LIMIT_ );
        const double k = TwoPhaseModel<dim>::dkrnds_at( KRN_LOW_SW_LIMIT_)/seff_mult;

        return splineValue( TwoPhaseModel<dim>::seff_, 0.0, KRN_LOW_SW_LIMIT_, 1.0, y, 0.0, k);
    }

    //        krn = [ 1 - Se^1/m]^2m * (1 - Se)^gamma
    double m = m_from_n(n_);
    double term = std::pow( 1. - std::pow( TwoPhaseModel<dim>::seff_, 1./m ), 2. * m );

    return std::pow( 1. - TwoPhaseModel<dim>::seff_, 1./3. ) * term;
 }

template<uint32_t dim>
double VanGenuchten<dim>::dkrwds_Phase() const
 {

    if ( TwoPhaseModel<dim>::seff_ <= 0.) return static_cast<double>(0.);
    if ( TwoPhaseModel<dim>::seff_ >= 1.) return static_cast<double>(0.);

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    // regularization
    if ( TwoPhaseModel<dim>::seff_ > KRW_HIGH_SW_LIMIT_ ){

        const double y = TwoPhaseModel<dim>::krw_at( KRW_HIGH_SW_LIMIT_ );
        const double k = TwoPhaseModel<dim>::dkrwds_at( KRW_HIGH_SW_LIMIT_)/seff_mult;

        return splineDerivative( TwoPhaseModel<dim>::seff_, KRW_HIGH_SW_LIMIT_, 1.0 , y, 1.0, k, 0.0)*seff_mult;
    }

    const double m = m_from_n(n_);
    const double term  = 1. - std::pow( TwoPhaseModel<dim>::seff_, 1./m );
    const double termM = pow( term , m);

    return (1.0-termM)*( (1.0 - termM)/2. + 2.0*termM*(1.0-term)/term )/sqrt(TwoPhaseModel<dim>::seff_)*seff_mult;
 }


template<uint32_t dim>
double VanGenuchten<dim>::dkrnds_Phase() const
 {

    if ( TwoPhaseModel<dim>::seff_ <= 0.) return static_cast<double>(0.);
    if ( TwoPhaseModel<dim>::seff_ >= 1.) return static_cast<double>(0.);

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    // regularization
    if ( TwoPhaseModel<dim>::seff_ < KRN_LOW_SW_LIMIT_ ){

        const double y = TwoPhaseModel<dim>::krn_at( KRN_LOW_SW_LIMIT_ );
        const double k = TwoPhaseModel<dim>::dkrnds_at( KRN_LOW_SW_LIMIT_)/seff_mult;

        return splineDerivative( TwoPhaseModel<dim>::seff_, 0.0, KRN_LOW_SW_LIMIT_, 1.0, y, 0.0, k)*seff_mult;
    }

    const double m = m_from_n(n_);
    const double term  = std::pow( TwoPhaseModel<dim>::seff_, 1./m );
    const double ratio1  = (1.0 - TwoPhaseModel<dim>::seff_)/(1.0 - term);
    const double ratio2  = term/TwoPhaseModel<dim>::seff_;

    return -pow(1.0-term, 2.0*m )*pow(1.0 - TwoPhaseModel<dim>::seff_, -2.0/3.0)*( 1./3. + 2.0*ratio1*ratio2)*seff_mult;

 }

template<uint32_t dim>
double VanGenuchten<dim>::pc_Phase( ) const
{
    /// Not consistent way of regularization
    /// Roman, 2013
    // only for the min saturation of water precautions are needed
    // the actual saturation is used instead of the effective saturation
    // if ( TwoPhaseModel<dim>::seff_ == 0. ) return TwoPhaseModel<dim>::pc_max_;
    // if ( TwoPhaseModel<dim>::seff_ == 1. ) return static_cast<double>(0.);

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    /// New way of regularization
    // linear regularization for lower part of sw range
    if ( TwoPhaseModel<dim>::seff_ < PC_LOW_SW_LIMIT_ ){

        const double pc_lim = TwoPhaseModel<dim>::pc_at( PC_LOW_SW_LIMIT_ );
        const double dpcds_lim = TwoPhaseModel<dim>::dpcds_at( PC_LOW_SW_LIMIT_ )/seff_mult;
        return pc_lim + dpcds_lim*( TwoPhaseModel<dim>::seff_ - PC_LOW_SW_LIMIT_ );

    // linear regularization for higher part of sw range
    }else if ( TwoPhaseModel<dim>::seff_ > PC_HIGH_SW_LIMIT_ ){

        const double pc_lim = TwoPhaseModel<dim>::pc_at( PC_HIGH_SW_LIMIT_ );
        const double dpcds_lim ( (0.0 - pc_lim)/( 1.0 - PC_HIGH_SW_LIMIT_ ) );
        return dpcds_lim*( TwoPhaseModel<dim>::seff_ - 1.0 );

    }


    double m = m_from_n( n_ );
    double term = std::pow( TwoPhaseModel<dim>::seff_, -1. / m ) - 1.;
    return std::pow( term, 1. / n_ ) / alpha_;

}


template<uint32_t dim>
double VanGenuchten<dim>::dpcds_Phase( ) const
{
    /// Not consistent way of regularization
    /// Roman, 2013
    //if ( TwoPhaseModel<dim>::seff_ < PC_LOW_SW_LIMIT_   ) return -MAX_CAPILLARY_PRESSURE_SLOPE_;
    //if ( TwoPhaseModel<dim>::seff_ > PC_HIIGH_SW_LIMIT_ ) return -MAX_CAPILLARY_PRESSURE_SLOPE_;

    /// New way of regularization
    // linear regularization for lower part of sw range
    if ( TwoPhaseModel<dim>::seff_ < PC_LOW_SW_LIMIT_   )
        return TwoPhaseModel<dim>::dpcds_at( PC_LOW_SW_LIMIT_ );

    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    // linear regularization for higher part of sw range
    if ( TwoPhaseModel<dim>::seff_ > PC_HIGH_SW_LIMIT_ )
        return (0.0 - TwoPhaseModel<dim>::pc_at(PC_HIGH_SW_LIMIT_))/( 1.0 - PC_HIGH_SW_LIMIT_ ) * seff_mult;

    const double m  = m_from_n(  n_ );
    const double term0 = std::pow( TwoPhaseModel<dim>::seff_, -1. / m );
    const double term1 = term0 - 1.;
    const double term2 = std::pow( term1, 1. / n_ );

    return -(term2 * term0) / (alpha_ * n_ * m * TwoPhaseModel<dim>::seff_ * term1) * seff_mult;
}




template<uint32_t dim>
double VanGenuchten<dim>::Sw_Phase( double pc ) const
{

    double m = m_from_n( n_ );
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    // useful for numerical calculations
    if( pc<0.0){

        const double pc_lim = TwoPhaseModel<dim>::pc_at( PC_HIGH_SW_LIMIT_ );
        const double dpcds_lim ( (0.0 - pc_lim)/( 1.0 - PC_HIGH_SW_LIMIT_ ) );
        TwoPhaseModel<dim>::seff_ =  pc/dpcds_lim + 1.0;
        return TwoPhaseModel<dim>::SeffToSw();

    }else{

        double term = std::pow( pc*alpha_, n_ ) + 1.0;
        TwoPhaseModel<dim>::seff_ =  pow( term, -m );
    }

    // linear regularization for low part of sw range
    if ( TwoPhaseModel<dim>::seff_ <= PC_LOW_SW_LIMIT_ ){

        const double pc_lim = TwoPhaseModel<dim>::pc_at( PC_LOW_SW_LIMIT_ );
        const double dpcds_lim = TwoPhaseModel<dim>::dpcds_at( PC_LOW_SW_LIMIT_ )/seff_mult;
        TwoPhaseModel<dim>::seff_ =  PC_LOW_SW_LIMIT_ + (pc - pc_lim)/dpcds_lim;

    // linear regularization for higher part of sw range
    }else if ( TwoPhaseModel<dim>::seff_ > PC_HIGH_SW_LIMIT_ ){

        const double pc_lim = TwoPhaseModel<dim>::pc_at( PC_HIGH_SW_LIMIT_ );
        const double dpcds_lim ( (0.0 - pc_lim)/( 1.0 - PC_HIGH_SW_LIMIT_ ));

        TwoPhaseModel<dim>::seff_ = pc/dpcds_lim + 1.0;

    }

    return TwoPhaseModel<dim>::SeffToSw();

}


template<uint32_t dim>
double VanGenuchten<dim>::dsdpc_Phase( double pc ) const
{

    double m = m_from_n( n_ );
    double Seff;
    const double seff_mult( 1.0/ (1.0 - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_ ) );

    if( pc<0.0)
        Seff = 1.5; // just means that it is more than 1
    else
        Seff = (TwoPhaseModel<dim>::Sw_at( pc )- TwoPhaseModel<dim>::swr_ ) *seff_mult;

    // linear regularization for lower part of sw range
    if ( Seff < PC_LOW_SW_LIMIT_   )
        return 1.0/TwoPhaseModel<dim>::dpcds_at( PC_LOW_SW_LIMIT_ );

    // linear regularization for higher part of sw range
    if ( Seff > PC_HIGH_SW_LIMIT_ )
        return ( 1.0 - PC_HIGH_SW_LIMIT_ )/(0.0 - TwoPhaseModel<dim>::pc_at(PC_HIGH_SW_LIMIT_))/seff_mult;


    double term = std::pow( pc*alpha_, n_ );
    return -m*n_*term*pow( term + 1.0, -m -1.0)/pc/seff_mult;

}


template<uint32_t dim>
double VanGenuchten<dim>::MaxFractionalFlowDerivative() const
 {
    return static_cast<double>(7.); // as computed with dfds_Phase method
 }


// from Sebastian's capillary pressure model

/**

Member function that computes the capillary pressure from the
pore parameters (pore and pore throat radii) and interfacial tension during the
visitation. Reference saturation is the non-wetting phase saturation.

@section arguments Input Arguments

Reference saturation, pore throat radius, pore radius.

@return The capillary pressure
*/
template<uint32_t dim>
double VanGenuchten<dim>::PoreParametersCapillaryPressure( double pt, double pr, double tension ) const
 {
    // saturation function f(satw) = arccos(2satw - 1) / pi
    double pc = std::acos( 2. * TwoPhaseModel<dim>::seff_ - 1. ) / PI;

    // capillary pressure according to Berg 75
    pc *= tension * ( 1. / pr + 1. / pt );

    return pc;
 }


template class VanGenuchten<1U>;
template class VanGenuchten<2U>;
template class VanGenuchten<3U>;

} // end namespace csp








