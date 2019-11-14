#include "HeterogeneityAndRateAwareModel.h"
#include "PropertyDatabase.h"
#include "ErrorHandler.h"
#include "OtwayCRC3_RockTypes.h"


// #define DEBUG_HETEROGENEITY_AWARE_MODEL


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
//   RRT_key_(database.StorageKey(rocktype)),
   pf_key_(database.StorageKey("fluid pressure")), // to calculate fluid pressure gradient
   vt_key_(database.StorageKey("velocity"))  // to calculate flow direction
{
}



template<size_t dim>
HeterogeneityAndRateAwareModel<dim>::~HeterogeneityAndRateAwareModel()
 {
 }
 
 



/**
    Computing the magnitude of the capillary pressure gradient.
    and reading rocktype and parameters for the van Genuchten relperm model.
*/
template<size_t dim>
void HeterogeneityAndRateAwareModel<dim>::Initialize( const Element<dim>& e )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // inferring input parameters from the rocktype
    // rocktype_ = static_cast<int>(e.Read( RRT_key_ ));


    // 1. base class parameter initialisation
    // --------------------------------------
    // if saturation is an element variable
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
         sw_ = TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
         TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
         TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
         TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
         TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );
      }
    // else warning is issued and values are interpolated to barycentre
    else {
         csmp_error.notice( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize", "variable placement set to nodal" );
         sw_ = TwoPhaseModel<dim>::sat_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::sat_key_ );
         TwoPhaseModel<dim>::mun_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::mun_key_ );
         TwoPhaseModel<dim>::muw_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::muw_key_ );
         TwoPhaseModel<dim>::rhn_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::rhn_key_ );
         TwoPhaseModel<dim>::rhw_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::rhw_key_ );
      }
   
    // permeability anisotropy
    const double64 kv = PermeabilityPerpendicularToLaminations();
    const double64 kh = PermeabilityParallelToLaminations();
    KK_ = 0.;
    KK_(0,0)          = kh;
    KK_(1,1)          = kv;

    // communicating assigned values to base class to get information for testing
    TwoPhaseModel<dim>::k_   = (flow_direction_ == HORIZONTAL) ? kh : kv;
    TwoPhaseModel<dim>::tensor_permeability_ = false; // TODO: Needs to change later
    TwoPhaseModel<dim>::swr_  = Swr_Composite();
    TwoPhaseModel<dim>::snr_  = 0.;
    TwoPhaseModel<dim>::ift_  = 0.035; // 35 mN/m water - CO2


    // 2. Flow-rate and direction dependent parameters
    // -----------------------------------------------
    // computing the magnitude of the pressure gradient
    e.Read( vt_key_, vt_ );
    flow_direction_   = ProminentFlowDirection();
    K_flow_direction_ = PermeabilityInFlowDirection( vt_ );
    grad_p_magnitude_ = PressureGradientMagnitude(e);
    Nc_               = Nc_kgradP_Version( grad_p_magnitude_ ); // capillary number
    RVC_              = RVC( Nc_ );
    Sw_VL_            = Sw_VL( RVC_ );
    // default values (VL case, where the sample saturation is uniform)
    sw_low_k_star_ = sw_high_k_star_ = sw_;
    // VERY IMPORTANT as the whole model is based on this
    FlowRateDependentLayerSaturations( Sw_VL_, RVC_, sw_low_k_star_, sw_high_k_star_ );
    // ATTENTION: do not move: seff depends on Nc, RVC and flow direction
    TwoPhaseModel<dim>::seff_ = EffectiveSaturation();

// DEBUGGING
/*
if ( isnan(TwoPhaseModel<dim>::TotalMobility()) )
  throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                        "mobt = NaN.");

if ( TwoPhaseModel<dim>::TotalMobility() < 0. )
  throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                        "mobt is negative.");

if ( TwoPhaseModel<dim>::TotalMobility() <= numeric_limits<double64>::epsilon() )
  throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                        "mobt is zero.");
*/
#ifdef DEBUG_HETEROGENEITY_AWARE_MODEL
Out(1);
// no flow
vt_(0) = 0.;
vt_(1) = 0.;
flow_direction_   = ProminentFlowDirection();
K_flow_direction_ = PermeabilityInFlowDirection( vt_ );
Nc_               = 0.;
RVC_              = RVC( Nc_ );
Sw_VL_            = Sw_VL( RVC_ );
grad_p_magnitude_ = 0.;
WriteRelativePermeabilityTable( "Maartje2_layer_parallel_Nc0", 0. );
// horizontal flow
vt_(0) = 1.0e-6;
vt_(1) = 0.;
flow_direction_   = ProminentFlowDirection();
K_flow_direction_ = PermeabilityInFlowDirection( vt_ );
Nc_               = 1.0e-6;
RVC_              = RVC( Nc_ );
Sw_VL_            = Sw_VL( RVC_ );
grad_p_magnitude_ = 1.0e-9;
WriteRelativePermeabilityTable( "Maartje2_layer_parallel_Nc-6", 1.0e-6 );
// vertical flow
vt_(0) = 0.;
vt_(1) = 1.0e-6;
flow_direction_   = ProminentFlowDirection();
K_flow_direction_ = PermeabilityInFlowDirection( vt_ );
WriteRelativePermeabilityTable( "Maartje2_layer_perpendicular_Nc-6", 1.0e-6 );
throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                       "written relative permeability curves to file for testing.");
#endif

 } // end Initialize








/// weighted permeability average for flow along layers
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::PermeabilityParallelToLaminations() const
 {
    // layer 1,3,5 are high_k layers, layer 2, 4 are low_k layers
    double64 LY_1(0.015),  LY_2 (0.01),  LY_3 (0.005), LY_4 (0.015), LY_5 (0.005);
    return (LY_1*k_high_ + LY_2*k_low_ + LY_3*k_high_ + LY_4*k_low_ + LY_5*k_high_)/(LY_1+LY_2+LY_3+LY_4+LY_5);
 }
 
 
 
/// weigthed harmonic mean of permeability for flow across layers
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::PermeabilityPerpendicularToLaminations() const
 {
    // layer 1,3,5 are high_k layers, layer 2, 4 are low_k layers
    double64 LY_1(0.015),  LY_2 (0.01),  LY_3 (0.005), LY_4 (0.015), LY_5 (0.005);
    const double64 sum_of_weights = LY_1+LY_2+LY_3+LY_4+LY_5;
    return sum_of_weights / (LY_1/k_high_ + LY_2/k_low_ + LY_3/k_high_ + LY_4/k_low_ + LY_5/k_high_);
 }



/**
    Assuming that the layers are horizontal, this method returns the permeability in the flow direction
    which by interpolation between the layer parallel and layer-perpendicular permeability.
*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::PermeabilityInFlowDirection( const VectorVariable<dim>& mv ) const
 {
    // copying the velocity and normalising it to one
    vc_ = vt_;
    vc_.EuclideanNormalize();

    // finding the permeability in the direction of the velocity vector
    vc_ = KK_ * vc_;

    // the length of vc_ is equal to the magnitude of permeability in the tarfet direction
    return vc_.Length();
   
 } // end PermeabilityInFlowDirection



/**
    Computes prominent direction of the flow of the fluid mixture.
 
    Uses total velocity, vt, to establish whether vertical component of flow is greater than horizontal one.
 
    @note if the flow velocity magnitude is below threshold value, HORIZONAL is returned.
*/
template<size_t dim>
typename HeterogeneityAndRateAwareModel<dim>::FLOW_DIRECTION HeterogeneityAndRateAwareModel<dim>::ProminentFlowDirection() const
 {
     if ( fabs(vt_[0]+vt_[1]) <= numeric_limits<double64>::epsilon() ) return HORIZONTAL;
   
     // comparing the vertical component with the horizontal magnitude of the flow
     const double64 horizontal_magnitude = ( dim == 2U )  ? vt_[0] : sqrt( vt_[0]*vt_[0] + vt_[2]*vt_[2] );
     FLOW_DIRECTION direction = ( fabs(vt_[1]) > horizontal_magnitude ) ? VERTICAL : HORIZONTAL;
     return direction;
 }



/**
    Returns the magnitude of the pressure gradient
*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::PressureGradientMagnitude( const Element<dim>& e ) const
 {
    array<double64, 3> gradP{ 0., 0., 0. };
    e.dN_AtBaryCenter( DN_ );
    for ( size_t i=0U; i<e.Nodes(); ++i ) {
         const double64 pf = e.N(i)->Read( pf_key_ );
         for ( size_t j=0U; j<dim; ++j )
         gradP[j] += DN_(j,i) * pf;
      }
    double64 grad_p_magnitude = 0.;
    for ( size_t j=0U; j<dim; ++j )
      grad_p_magnitude += gradP[j] * gradP[j];
    grad_p_magnitude = sqrt(grad_p_magnitude);
   
    return grad_p_magnitude;

} // end PressureGradientMagnitude



/// pressure gradient form: Nc = k ||grad p|| / sigma
template<size_t dim>
double64 csmp::HeterogeneityAndRateAwareModel<dim>::Nc_kgradP_Version( double64 pf_gradient_magnitude ) const
 {
     assert( pf_gradient_magnitude != UNSPECIFIED );
     return PermeabilityParallelToLaminations() * pf_gradient_magnitude / TwoPhaseModel<dim>::ift_;
 }



/**
    Viscous - capillary force balance (RVC)
 
    Uses the capillary number and the interfacial tension as input parameters.

    RVC = ((RVC_high_k * (Ly1+Ly3+Ly5)) + (RVC_low_k * (Ly2+Ly4)))/(Ly1+Ly2+Ly3+Ly4+Ly5);

*/
template<size_t dim>
double64 csmp::HeterogeneityAndRateAwareModel<dim>::RVC( double64 Ncap ) const
 {
    const double64 rvc_low_k  = ((Ncap * TwoPhaseModel<dim>::ift_) * L_low_/2.) / (dPc_ * k_low_);
    const double64 rvc_high_k = ((Ncap * TwoPhaseModel<dim>::ift_) * L_high_/3.) / (dPc_ * k_high_);
   
    return (rvc_high_k * L_high_ + rvc_low_k * L_low_) / (L_low_ + L_high_);
 }




 /**
    @attention This is where the actual water saturation in the cell enters the computation.
 
    Since Sw_CL can be written as a function of the Sw_VL (curve fit) and cell average water saturation
    (Sw(Nc)) is known, Sw_VL can be calculated. Once you know Sw_VL you can calculate Sw_CL.
    So, after substitution of

    𝑆𝑤_𝐶L(𝑐𝑒𝑙𝑙)=0.3365+0.6711∗𝑆𝑤_𝑉𝐿(𝑐𝑒𝑙𝑙)

    into

    SwNC := RVC*Sw_CL/(RVC+1)

    we get

    Sw_VL = 1.490090896*SwNC*RVC-.5014155864*RVC+1.490090896*SwNC

    where SwNc   Sw_cell(Nc)
 
    @attention Sw_VL must be limited to [0,1], else kri will have NaN values.

*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::Sw_VL( double64 RVC ) const
 {
    //return ((sw_ * (RVC + 1.)) - 0.3365) / (RVC + 0.6711);
    return min( max( ((sw_ * (RVC + 1.)) - 0.3365) / (RVC + 0.6711), 0. ), 1. );
 }


/// by correlation
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::Sw_CL( double64 sw_VL ) const
 {
    // limited version: return min( 0.3365 + 0.6711 * Sw_VL_, 1. );
    return 0.3365 + 0.6711 * Sw_VL_;
 }




/**
    From the average cell saturation at the viscous limit, and the given capillary number and
    viscous-to-capillary force ratio, this method computes the effective saturations in the individual laminations.
 
    @note the effective saturations are bracketed to 0..1 to avoid error propagation when the model is applied outside of its
    established saturation range.
 
    TODO: extend this method to contain both, the horizontal and vertical flow cases
    TODO: include the VL and CL cases: in vertical VL case sw=s_ave and is the same in all layers
*/
template<size_t dim>
void HeterogeneityAndRateAwareModel<dim>::FlowRateDependentLayerSaturations( double64 sw_VL, double64 RVC, double64& sw_low_k_star, double64& sw_high_k_star ) const
 {
    // ESSENTIAL: using max to eliminate negative values
    const double64 Sw_CL_high_k = max( -1. + 2. * sqrt( sw_VL ), 0. );
    const double64 Sw_CL_low_k  = 0.55 + 0.45 * pow( Sw_CL_high_k, 0.29 );

    // limiting these water saturations to the residual saturation of water
    // TODO: how about capillary desaturation in highly permeable layers; put at least a warning here for sufficienctly high Nc
    const double64 sw_high_k = max( (RVC * sw_VL + Sw_CL_high_k) / (RVC + 1.), Swi_high_ );
    const double64 sw_low_k  = max( (RVC * sw_VL + Sw_CL_low_k) / (RVC + 1.), Swi_low_ );

    // Sw_star_FSst=(Sw-Swi_FSst)/(1-Swi_FSst)
    sw_high_k_star = (sw_high_k - Swi_high_) / (1. - Swi_high_);
    sw_low_k_star  = (sw_low_k  - Swi_low_)  / (1. - Swi_low_);

 } // end FlowRateDependentLayerSaturations






/**
    Maartje, write-up page 3.
 
    6.  Use saturations found in step 5 to calculate relative permeability using the VanGenuchten relationship
 
    Sw_star_FSst=(Sw-Swi_FSst)/(1-Swi_FSst)
    Krw_ave=((Krw_high_k.*k_high.*(Ly1+Ly3+Ly5))+(Krw_low_k.*k_low.*(Ly2+Ly4)))./(k_ave.*(Ly1+Ly2+Ly3+Ly4+Ly5));

*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krw_Phase() const
 {
    if ( sw_ <= Swr_Composite() ) return 0.;
   
    // layer-parallel case
    // -------------------
    if ( flow_direction_ == HORIZONTAL ) {
         // Krw_high_k=(Sw_high_k_star^0.5)*(1-(1-Sw_high_k_star^(1/m_high))^m_high)^2
         double64 term = 1. - pow( 1. - pow( sw_high_k_star_, (1./m_high_)), m_high_ );
         const double64 Krw_high_k = sqrt( sw_high_k_star_ ) * (term * term);
 
         // Krw_low_k=(Sw_low_k^0.5)*(1-(1-Sw_low_k^(1/m_low))^m_low)^2
         term = 1. - pow( 1. - pow( sw_low_k_star_, (1./m_low_)), m_low_ );
         const double64 Krw_low_k = sqrt( sw_low_k_star_ ) * (term * term);
   
         // Krw_ave=((Krw_high_k.*k_high.*(Ly1+Ly3+Ly5))+(Krw_low_k.*k_low.*(Ly2+Ly4)))./(k_ave.*(Ly1+Ly2+Ly3+Ly4+Ly5));
         return min( max( (L_high_ * Krw_high_k * k_high_ + L_low_ * Krw_low_k * k_low_) / (k_high_ * L_high_ + k_low_ * L_low_), 0. ), 1. );
      }
   
    // layer-perpendicular case
    // ------------------------
    const double64 KvKh_ratio = K_flow_direction_ / KK_(0,0);
    // Maartje 2/9/2019
    double64 u = vt_.Length() * KvKh_ratio;
    //If ux=>5e-5 m/s use viscous limit relative permeability curves.
    if ( u >= 5.5e-5 ) return krw_VL_LayerPerpendicular();
    //If ux=<5e-7 m/s use capillary limit relative permeablity curves.
    if ( u <= 5.0e-7 ) return krw_CL_LayerPerpendicular();
   
    const double64 krw = (0.6554e37 * pow(u,0.8e1) - 0.1456e34 * pow(u,0.7e1) + 0.1352e30 * pow(u,0.6e1) - 0.6797e25 * pow(u,0.5e1) + 0.2007e21 * pow(u,0.4e1) - 0.3526e16 * pow(u,0.3e1) + 0.3547e11 * u * u - 0.1815e6 * u + 0.9447e0) * pow(sw_,-0.5214e19 * pow(u,0.4e1) + 0.6725e15 * pow(u,0.3e1) - 0.3122e11 * u * u + 0.6112e6 * u + 0.4933e1) - 0.1177e37 * pow(u,0.8e1) + 0.2601e33 * pow(u,0.7e1) - 0.2404e29 * pow(u, 0.6e1) + 0.1205e25 * pow(u,0.5e1) - 0.3558e20 * pow(u,0.4e1) + 0.6307e15 * pow(u,0.3e1) - 0.6549e10 * u * u + 0.3644e5 * u - 0.5009e-1;
   
    return max( krw * KvKh_ratio, 0. );
   
 } // end krw_Phase




/**
    CO2 relative permeability
 
    Krnw_ave=((Krnw_high_k.*k_high.*(Ly1+Ly3+Ly5))+(Krnw_low_k.*k_low.*(Ly2+Ly4)))./(k_ave.*(Ly1+Ly2+Ly3+Ly4+Ly5));

*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krn_Phase() const
 {
    // introducing a percolation threshold for the non-wetting phase
    const double64 percolation_threshold_nw(0.01);
    if ( (1. - sw_) <= percolation_threshold_nw ) return 0.;

    // layer-parallel case
    // ------------------------
    if ( flow_direction_ == HORIZONTAL ) {
         // Krnw_high_k=((1-Sw_high_k)^2).*((1-Sw_high_k^2))
         const double64 Krnw_high_k = ((1. - sw_high_k_star_)*(1. - sw_high_k_star_)) * (1. - sw_high_k_star_ * sw_high_k_star_);

         // Krnw_low_k=((1-Sw_low_k)^2)*((1-Sw_low_k^2))
         const double64 Krnw_low_k = ((1. - sw_low_k_star_)*(1. - sw_low_k_star_)) * (1. - sw_low_k_star_ * sw_low_k_star_);

         // Krnw_ave=((Krnw_high_k.*k_high.*(Ly1+Ly3+Ly5))+(Krnw_low_k.*k_low.*(Ly2+Ly4)))./(k_ave.*(Ly1+Ly2+Ly3+Ly4+Ly5));
         return min( max( (L_high_ * Krnw_high_k * k_high_ + L_low_ * Krnw_low_k * k_low_) / (k_high_ * L_high_ + k_low_ * L_low_), 0. ), 1. );
      }
   
    // layer-perpendicular case
    // ------------------------
    const double64 KvKh_ratio = K_flow_direction_ / KK_(0,0);
    // Maartje 2/9/2019
    double64 u = vt_.Length() * KvKh_ratio;
    //If ux=>5e-5 m/s use viscous limit relative permeability curves.
    if ( u >= 5.5e-5 ) return krw_VL_LayerPerpendicular();
    //If ux=<5e-7 m/s use capillary limit relative permeablity curves.
    if ( u <= 5.0e-7 ) return krw_CL_LayerPerpendicular();

    const double64 krn = (-0.3228e9 * u * u + 0.3561e5 * u + 0.2139e0) * pow(-0.2838e37 * pow(u,0.8e1) + 0.6254e33 * pow(u,0.7e1) - 0.5748e29 * pow(u,0.6e1) + 0.2856e25 * pow(u,0.5e1) - 0.8316e20 * pow(u,0.4e1) + 0.1438e16 * pow(u,0.3e1) - 0.1422e11 * u * u + 0.707e5 * u + 0.8596e0 - sw_, 0.2e1);
   
    return min( krn * KvKh_ratio, 1. );

 } // end krn_Phase









/// for the wetting phase
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::MaxFractionalFlowDerivative() const
 {
    return static_cast<double64>(5.6); // as computed with dfds method
 }




/**

Pc_drain_high = Pd_high*((Sw_star_high)^(-1/m_high)-1)^(1-m_high); % capillary pressure curve high permeable layer

Pc_drain_low = Pd_low*((Sw_star_low)^(-1/m_low)-1)^(1-m_low);      % capillary pressure curve low permeable layer

Apply averaging to get vertical flow.

@attention pc covers the full saturation range, pc(sw) is capped based on maximum dpcds.

TODO: extend capillary pressure curve beyond saturation endpoints.

*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::pc_Phase() const
 {
    double64 sw_low_k_star, sw_high_k_star;
    FlowRateDependentLayerSaturations( Sw_VL_, RVC_, sw_low_k_star, sw_high_k_star );

    if ( flow_direction_ == HORIZONTAL ) {
         // Pc_drain_high = Pd_high * ((Sw_star_high)^(-1/m_high) - 1)^(1-m_high)
         const double64 Pc_drain_high = pd_high_ * pow( (pow( max(Swi_high_,sw_high_k_star), -1./m_high_ ) - 1.), 1. - m_high_ );
         // Pc_drain_low=Pd_low.*((Sw_star_low).^(-1/m_low)-1).^(1-m_low)
         const double64 Pc_drain_low  = pd_low_ * pow( (pow( max(Swi_low_,sw_low_k_star), -1./m_low_ ) - 1.), 1. - m_low_ );
         // thickness-weighted averaging, value must be greater than entry pressure of the more permeable layer
         return max( min( (L_high_ * Pc_drain_high + L_low_ * Pc_drain_low) / (L_high_ + L_low_), TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ ), pd_high_ );
      }

    // layer-perpendicular case (just using the properties of the low_k layer
    // ----------------------------------------------------------------------
// HAS NO ENTRY PRESSURE EFFECT:    return pd_low_ * pow( (pow( max(Swi_low_,sw_low_k_star), -1./m_low_ ) - 1.), 1. - m_low_ );
    // capping pc_max by dissallowing Sw values below the irreducible saturation of composite
    const double64 sw = max( sw_, Swr_Composite() );
    const double64 lambda = (m_low_ / (1. - m_low_)) * (1. - pow( sw, 1./m_low_ ));
    assert( lambda <= 10. );
    // Brooks-Corey for low-k layer instead (checked against Helmig, 97)
    return ( sw > 1. - numeric_limits<double64>::epsilon() ) ? pd_low_ : pd_low_ * pow( sw, -1. / lambda );

 } // end pc_Phase




/**
    SKM: thickness-weighted average of the effective saturations in the high and the low k laminations.
 
    @attention must always be called first because Seff must be communicated with the base class.
*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::EffectiveSaturation() const
 {
    if ( flow_direction_ == HORIZONTAL )
      return (L_high_ * sw_high_k_star_ + L_low_ * sw_low_k_star_) / (L_high_ + L_low_);
   
    return (sw_ - TwoPhaseModel<dim>::swr_) / (1. - TwoPhaseModel<dim>::swr_ - TwoPhaseModel<dim>::snr_);
   
 } // end EffectiveSaturation



/**
    Capillary pressure - saturation derivative estimated numerically.
 
    @attention dpcdsw is zero as soon as one of the phases is immobile
    to avoid that CFL is influenced while capillary spreading is not possible.
*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::dpcds_Phase() const
 {
    assert( sw_ >= 0. );
    assert( sw_ <= 1. );
    const double64  h(0.01), sw(sw_);
//    assert( h > 0. );
//    assert( h <= 0.01 );

    // deal with the more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         sw_ = 1.;
         double64 pc1 = pc_Phase();
         sw_ = 1. - h;
         double64 pc2 = pc_Phase();
         return (pc1 - pc2) / h;
      }
  
    // low water saturation
    if ( sw <= h ) {
         sw_ = h;
         double64 pc1 = pc_Phase();
         sw_ = 0.;
         double64 pc2 = pc_Phase();
         return (pc1 - pc2) / h;
      }

    // water saturation between the endpoints
    sw_ = sw + h;
    double64 pc1 = pc_Phase();
    sw_ = sw - h;
    double64 pc2 = pc_Phase();
    // resetting sw value
    sw_ = sw;
   
    return (pc1 - pc2) / (2. * h);
}



/**
   Irreducible water saturation of the facies association of rocktypes which form the composite.
*/
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::Swr_Composite() const
 {
     return (Swi_low_ * phi_low_ * L_low_ + Swi_high_ * phi_high_ * L_high_) / (L_low_ * phi_low_ + L_high_ * phi_high_);
 }





    // limit approximations for vertical flow arrived at by curve fitting (6/9/19)
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krw_VL_LayerPerpendicular() const
 {
    const double64 C1(6.342), C2(0.4549);
    // Krw_VL = Sw(ux)^0.5 * (1 - (1 - Sw(ux)^C1)^C2
    return max( sqrt(sw_) * pow( 1. - (1. - pow( sw_, C1 ) ), C2 ), 0.);
 }



// piecewise defined functions 
template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krn_VL_LayerPerpendicular() const
 {
    const double64 C1(1.491 ), C2(1.945), C3(2.472);
    // Krnw_VL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
    return max( C1 * (1. - pow( sw_, C2 )) * (1. - pow( sw_, C3 )), 0.);
 }



template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krw_CL_LayerPerpendicular() const
 {
    const double64 C1 = (sw_ > 0.83) ? 2.717  : 7.39;
    const double64 C2 = (sw_ > 0.83) ? 0.5099 : 1.557;
   
    // Krw_CL = (Sw(ux)^0.5)*(1-(1-Sw(ux)^C1)^C2)
    return max( sqrt(sw_) * pow( 1. - (1. - pow( sw_, C1 ) ), C2 ), 0.);
 }



template<size_t dim>
double64 HeterogeneityAndRateAwareModel<dim>::krn_CsL_LayerPerpendicular() const
 {
    const double64 C1 = (sw_ > 0.3) ? 3.1   : 2.199;
    const double64 C2 = (sw_ > 0.3) ? 5.033 : 4.148;
    const double64 C3 = (sw_ > 0.3) ? 2.758 : 30.1;
   
    // Krnw_CL = C1 * ((1-Sw(ux))^C2) * ((1-Sw(ux)^C3))
    return max( C1 * (1. - pow( sw_, C2 )) * (1. - pow( sw_, C3 )), 0.);
}











/**
    writes textfile with sw, krw(sw,Nc), krn(sw,Nc), and pc(sw) values computed for (composite) rocktype in 0.05 saturation increments
*/
template<size_t dim>
void HeterogeneityAndRateAwareModel<dim>::WriteRelativePermeabilityTable( const char* filename, double64 Ncap )
 {
    ofstream  ofs( string(filename) + ".txt" );
   
    const string rocktype("laminated-sand-silt");
    const double64 original_sw(sw_);
    const double64 original_Nc(Nc_);
    Nc_ = Ncap;

    ofs <<"sw(Nc="<< Nc_ <<")\t krw(sw,Nc="<< Nc_ <<",rocktype="<< rocktype <<")\t krnw(sw,Nc) \t pc(sw,Nc)\n";
    for ( double64 sw(0.); sw<1.05; sw+=0.05 )
      {
         // dynamic parameters
         sw_    = sw;
         RVC_   = RVC( Nc_ );
         Sw_VL_ = Sw_VL( RVC_ );
         FlowRateDependentLayerSaturations( Sw_VL_, RVC_, sw_low_k_star_, sw_high_k_star_ );
         ofs << sw << "\t"<< krw_Phase();
         ofs <<"\t"<< krn_Phase();
         ofs <<"\t"<< pc_Phase();
         ofs << endl;
      }
   
    // restoring current saturation
    sw_    = original_sw;
    Nc_    = original_Nc;
    RVC_   = RVC( Nc_ );
    Sw_VL_ = Sw_VL( RVC_ );
    FlowRateDependentLayerSaturations( Sw_VL_, RVC_, sw_low_k_star_, sw_high_k_star_ );

 } // end WriteRelativePermeabilityTable
  





template<size_t dim>
void HeterogeneityAndRateAwareModel<dim>::Out( size_t phase ) const
 {
    TwoPhaseModel<dim>::Out(phase);
    cout <<"\nHeterogeneityAndRateAwareModel<"<< dim << ">::Out: Additional properties and return values of functions: "<< endl;
    cout <<"\nelement properties:";
    cout <<"\n                                 velocity (m/s): "<< vt_;
    if ( ProminentFlowDirection() == HORIZONTAL ) cout <<", dominantly horizontal flow.";
    else cout <<", dominantly vertical flow.";
    cout <<"\n                           capillary number, Nc: "<< Nc_;
    cout <<"\n               layer-parellel permeability (m2): "<< PermeabilityParallelToLaminations();
    cout <<"\n          layer-perpendicular permeability (m2): "<< PermeabilityPerpendicularToLaminations();
    cout <<"\nratio between viscous and capillary forces, RVC: "<< RVC( Nc_ );
    cout <<"\n      average water saturation in composite, sw: "<< sw_;
    cout <<"\n                                        krw(sw): "<< krw_Phase();
    cout <<"\n                                        krn(sw): "<< krn_Phase();
    cout <<"\n                                         pc(sw): "<< pc_Phase();

    cout <<"\n  MAX_CAPILLARY_PRESSURE_DERIVATIVE: "<< TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_ << endl << endl;

 } // end Out
 
 


template class HeterogeneityAndRateAwareModel<1U>;
template class HeterogeneityAndRateAwareModel<2U>;
template class HeterogeneityAndRateAwareModel<3U>;

} // end namespace csmp


























