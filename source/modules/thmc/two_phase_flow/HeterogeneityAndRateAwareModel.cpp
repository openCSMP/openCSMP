#include "HeterogeneityAndRateAwareModel.h"
#include "PropertyDatabase.h"
#include "ErrorHandler.h"

//#define DEBUG_HETEROGENEITY_AWARE_MODEL
//#define TURN_OFF_RATE_AWARE


using namespace std;

namespace csmp {

template<uint32_t dim>
HeterogeneityAndRateAwareModel<dim>::HeterogeneityAndRateAwareModel( const PropertyDatabase<dim>& database,
                                                                     const bool sw_ro_mu_placement )
  
 : TwoPhaseModel<dim>(database, "permeability",  
                     "viscosity oil", "viscosity water",
                     "density oil", "density water", "saturation water",
                     "residual saturation non-wetting phase",
                     "residual saturation wetting phase",
                      sw_ro_mu_placement ),
   RRT_key_(database.StorageKey("rock type")),
   pf_key_(database.StorageKey("fluid pressure")), // to calculate fluid pressure gradient
   kh_key_(database.StorageKey("permeability")), // scalar representing horizontal direction
   kv_key_(database.StorageKey("vertical permeability")), // scalar representing vertical direction
   phi_key_(database.StorageKey("porosity")),
   vt_key_(database.StorageKey("velocity")),  // to calculate flow direction
   rocktype_(0), // WELL
   is_composite_(false),
   K_flow_direction_(HORIZONTAL), 
   k_low_(numeric_limits<double>::quiet_NaN()), 
   k_high_(numeric_limits<double>::quiet_NaN()), 
   K_reduction_in_flow_direction_(numeric_limits<double>::quiet_NaN()),
   L_low_(numeric_limits<double>::quiet_NaN()), 
   L_high_(numeric_limits<double>::quiet_NaN()),       ///< permeability in flow direction; smallest over highest permeability
   vt_magnitude_(numeric_limits<double>::quiet_NaN()),
   vt_magnitude_x_(numeric_limits<double>::quiet_NaN()),
   vt_magnitude_y_(numeric_limits<double>::quiet_NaN()), 
   krw_(numeric_limits<double>::quiet_NaN()), 
   krn_(numeric_limits<double>::quiet_NaN()),
   krw_parallel_(numeric_limits<double>::quiet_NaN()), 
   krw_crossflow_(numeric_limits<double>::quiet_NaN()),
   krn_parallel_(numeric_limits<double>::quiet_NaN()), 
   krn_crossflow_(numeric_limits<double>::quiet_NaN()),
   phi_(numeric_limits<double>::quiet_NaN()), 
   pd_(numeric_limits<double>::quiet_NaN()), 
   pd_high_(numeric_limits<double>::quiet_NaN()), 
   pd_low_(numeric_limits<double>::quiet_NaN()), 
   m_VG_(numeric_limits<double>::quiet_NaN()), 
   bcp_(numeric_limits<double>::quiet_NaN()), 
   bcp_high_(numeric_limits<double>::quiet_NaN()), 
   bcp_low_(numeric_limits<double>::quiet_NaN()), 
   Swi_pc_(numeric_limits<double>::quiet_NaN()), 
   dPc_(numeric_limits<double>::quiet_NaN()),
   pd_flow_direction_(numeric_limits<double>::quiet_NaN()), 
   pc_flow_direction_(numeric_limits<double>::quiet_NaN()),
   grad_p_magnitude_(numeric_limits<double>::quiet_NaN()) 
{
}



template<uint32_t dim>
HeterogeneityAndRateAwareModel<dim>::~HeterogeneityAndRateAwareModel()
 {
 }
 
 

template<uint32_t dim>
int32_t HeterogeneityAndRateAwareModel<dim>::RockType( const Element<dim>& e ) const {
     const double rock_type = e.Read(RRT_key_);
     assert( !isnan(rock_type) );
     return static_cast<int32_t>( rock_type );
  }



/** 
    Reads velocity variable, vt_, normalises it (vt_normalised_) and computes magnitude (vt_magnitude_).
    If the velocity is zero, horizontal flow will be indicated.
*/
template<uint32_t dim>
void HeterogeneityAndRateAwareModel<dim>::InitializeVelocity( const Element<dim>& e )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    e.Read( vt_key_, vt_ );
    if ( isnan(vt_[0]) || isnan(vt_[1]) )
      csmp_error.notice( ERROR, "HeterogeneityAndRateAwareModel<dim>::Initialize", "the velocity variable has not been initialised.");
    
    vt_magnitude_ = vt_.Length();
    vt_magnitude_x_ = fabs(vt_[0]);
    vt_magnitude_y_ = fabs(vt_[1]);
    // in zero velocity case, the horizontal relative permeability is set to dominate
    if ( vt_magnitude_ <= numeric_limits<double>::epsilon() * 100. ) {
         vt_normalised_(0) = 1.; 
         vt_normalised_(1) = 0.; 
      }
    else vt_normalised_ = vt_ / vt_magnitude_;
    
#ifdef TURN_OFF_RATE_AWARE
    vt_magnitude_ = vt_magnitude_x_ = vt_magnitude_y_ = 1.0e-12; //set to below capillary limit
#endif      
 }



/**
    Computing the magnitude of the capillary pressure gradient.
    and reading rocktype and parameters for the van Genuchten relperm model.
    
    @attention this method ignores any end-point saturations provided in the input dataset (discretised on the model)
    
    @attention this method only deals with horizontal or vertical laminations; oblique directions are not considered.
*/
template<uint32_t dim>
void HeterogeneityAndRateAwareModel<dim>::Initialize( const Element<dim>& e )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 1. base class parameter initialisation
    // --------------------------------------
    // NB: residual saturations are ignored because they are redefined here on the basis of the rocktype
    // if saturation is an element variable
    if (!TwoPhaseModel<dim>::sw_ro_mu_placement_) {
         Sw_ = TwoPhaseModel<dim>::sat_ = e.Read( TwoPhaseModel<dim>::sat_key_ );
         TwoPhaseModel<dim>::mun_ = e.Read( TwoPhaseModel<dim>::mun_key_ );
         TwoPhaseModel<dim>::muw_ = e.Read( TwoPhaseModel<dim>::muw_key_ );
         TwoPhaseModel<dim>::rhn_ = e.Read( TwoPhaseModel<dim>::rhn_key_ );
         TwoPhaseModel<dim>::rhw_ = e.Read( TwoPhaseModel<dim>::rhw_key_ );
      }
    // else warning is issued and values are interpolated to barycentre
    else {
         csmp_error.notice( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize", "variable placement set to nodal" );
         Sw_ = TwoPhaseModel<dim>::sat_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::sat_key_ );
         TwoPhaseModel<dim>::mun_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::mun_key_ );
         TwoPhaseModel<dim>::muw_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::muw_key_ );
         TwoPhaseModel<dim>::rhn_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::rhn_key_ );
         TwoPhaseModel<dim>::rhw_ = e.PropertyValueAtBaryCenter( TwoPhaseModel<dim>::rhw_key_ );
      }
          

    // 2. rocktype and flow velocity
    // ----------------------------------------------------------------
    rocktype_ = RockType(e);
    assert( rocktype_ >= 0 && rocktype_ <= 16 );
    InitializeVelocity(e);
    // default: standard rocktype
    is_composite_ = false;
    TwoPhaseModel<dim>::tensor_permeability_ = false;
    TwoPhaseModel<dim>::ift_  = 0.035; // 35 mN/m water - CO2
    dPc_ = 0.;

    switch( rocktype_ ) {
         case 0: // WELL
             k_low_ = k_high_         = get<0>(Otway_.rocktype_).k_;
             // reading single-petro-type rocktypes; scalar permeability is read from the elements
             phi_                     = get<0>(Otway_.rocktype_).phi_;
             pd_                      = get<0>(Otway_.rocktype_).pd_;
             m_VG_                    = get<0>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<0>(Otway_.rocktype_).Swi_pc_; // swr used in fitting the pc curve
             bcp_                     = 0.; // linear model
             // these are not used in any of the underlying functions of the 2-phase model
             TwoPhaseModel<dim>::swr_ = get<0>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<0>(Otway_.rocktype_).Sgr_;
             krw_                     = get<0>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<0>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 1: // H_Mst
             k_low_ = k_high_ = get<1>(Otway_.rocktype_).k_;
             phi_                     = get<1>(Otway_.rocktype_).phi_;
             pd_                      = get<1>(Otway_.rocktype_).pd_;
             m_VG_                    = get<1>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<1>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<1>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<1>(Otway_.rocktype_).Sgr_;
             krw_                     = get<1>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<1>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 2: // M_CbSst_Mst (composite but not rate-dependent)
             is_composite_ = true;
             L_low_         = get<2>(Otway_.rocktype_).LY_low_; L_high_ = get<2>(Otway_.rocktype_).LY_high_;
             k_low_         = get<2>(Otway_.rocktype_).k_low_;  k_high_ = get<2>(Otway_.rocktype_).k_high_;
             pd_high_       = get<2>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<2>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<2>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<2>(Otway_.rocktype_).m_low_ );
             phi_           = get<2>(Otway_.rocktype_).phi_; 
             pd_            = get<2>(Otway_.rocktype_).pd_;
             m_VG_          = get<2>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<2>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<2>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<2>(Otway_.rocktype_).Swi_pc_; // for this type Swi=Swi_pc
             TwoPhaseModel<dim>::snr_ = get<2>(Otway_.rocktype_).Sgr_;
             // no rate dependence for this carbonate cemented rocktype
             krw_parallel_  = get<2>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<2>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<2>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<2>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 3: // H_CbSst
             k_low_ = k_high_         = get<3>(Otway_.rocktype_).k_;
             phi_                     = get<3>(Otway_.rocktype_).phi_;
             pd_                      = get<3>(Otway_.rocktype_).pd_;
             m_VG_                    = get<3>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<3>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<3>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<3>(Otway_.rocktype_).Sgr_;
             krw_                     = get<3>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<3>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 4: // P_Mst_Slt: COMPOSITE
             is_composite_ = true;
             // layer-specific properties
             L_low_         = get<4>(Otway_.rocktype_).LY_low_; L_high_ = get<4>(Otway_.rocktype_).LY_high_;
             k_low_         = get<4>(Otway_.rocktype_).k_low_;  k_high_ = get<4>(Otway_.rocktype_).k_high_;
             pd_high_       = get<4>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<4>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<4>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<4>(Otway_.rocktype_).m_low_ );
             // composite properties
             phi_           = get<4>(Otway_.rocktype_).phi_; // average porosity
             pd_            = get<4>(Otway_.rocktype_).Pd( vt_magnitude_x_ ); 
             m_VG_          = get<4>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<4>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<4>(Otway_.rocktype_).dPc_;
             // saturation endpoints        
             TwoPhaseModel<dim>::swr_ = get<4>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<4>(Otway_.rocktype_).Sgr_; // make rate dependent?
             // relative permeability parallel and perpendicular to layers
             krw_parallel_  = get<4>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<4>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<4>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<4>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );   
           break;
         case 5: // M_CbSst_Slt: COMPOSITE - not rate dependent
             is_composite_ = true;
             L_low_         = get<5>(Otway_.rocktype_).LY_low_; L_high_ = get<5>(Otway_.rocktype_).LY_high_;
             k_low_         = get<5>(Otway_.rocktype_).k_low_;  k_high_ = get<5>(Otway_.rocktype_).k_high_;
             pd_high_       = get<5>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<5>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<5>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<5>(Otway_.rocktype_).m_low_ );
             phi_           = get<5>(Otway_.rocktype_).phi_; 
             pd_            = get<5>(Otway_.rocktype_).pd_;
             m_VG_          = get<5>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<5>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<5>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<5>(Otway_.rocktype_).Swi_pc_; // for this type Swi=Swi_pc
             TwoPhaseModel<dim>::snr_ = get<5>(Otway_.rocktype_).Sgr_;
             // no rate dependence for this carbonate cemented rocktype
             krw_parallel_  = get<5>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<5>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<5>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<5>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 6: // P_Mst_FSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<6>(Otway_.rocktype_).LY_low_; L_high_ = get<6>(Otway_.rocktype_).LY_high_;
             k_low_ = get<6>(Otway_.rocktype_).k_low_;  k_high_ = get<6>(Otway_.rocktype_).k_high_;
             pd_high_       = get<6>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<6>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<6>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<6>(Otway_.rocktype_).m_low_ );
             phi_           = get<6>(Otway_.rocktype_).phi_; 
             pd_            = get<6>(Otway_.rocktype_).Pd( vt_magnitude_x_ );
             m_VG_          = get<6>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<6>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<6>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<6>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<6>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<6>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<6>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<6>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<6>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );
           break;
         case 7: // M_CbSst_FSst: COMPOSITE - not rate dependent
             is_composite_ = true;
             L_low_ = get<7>(Otway_.rocktype_).LY_low_; L_high_ = get<7>(Otway_.rocktype_).LY_high_;
             k_low_ = get<7>(Otway_.rocktype_).k_low_;  k_high_ = get<7>(Otway_.rocktype_).k_high_;
             pd_high_       = get<7>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<7>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<7>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<7>(Otway_.rocktype_).m_low_ );
             phi_           = get<7>(Otway_.rocktype_).phi_; 
             pd_            = get<7>(Otway_.rocktype_).pd_;
             m_VG_          = get<7>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<7>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<7>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<7>(Otway_.rocktype_).Swi_pc_; // Swi = Sw_pc
             TwoPhaseModel<dim>::snr_ = get<7>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<7>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<7>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<7>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<7>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 8: // P_Mst_CSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<8>(Otway_.rocktype_).LY_low_; L_high_ = get<8>(Otway_.rocktype_).LY_high_;
             k_low_ = get<8>(Otway_.rocktype_).k_low_;  k_high_ = get<8>(Otway_.rocktype_).k_high_;
             pd_high_       = get<8>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<8>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<8>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<8>(Otway_.rocktype_).m_low_ );
             phi_           = get<8>(Otway_.rocktype_).phi_; 
              pd_            = get<8>(Otway_.rocktype_).Pd( vt_magnitude_x_ );
             m_VG_          = get<8>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<8>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<8>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<8>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<8>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<8>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<8>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<8>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<8>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );
           break;
         case 9: // M_CbSst_CSst: COMPOSITE - not rate dependent 
             is_composite_ = true;
             L_low_ = get<9>(Otway_.rocktype_).LY_low_; L_high_ = get<9>(Otway_.rocktype_).LY_high_;
             k_low_ = get<9>(Otway_.rocktype_).k_low_;  k_high_ = get<9>(Otway_.rocktype_).k_high_;
             pd_high_       = get<9>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<9>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<9>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<9>(Otway_.rocktype_).m_low_ );
             phi_           = get<9>(Otway_.rocktype_).phi_; 
             pd_            = get<9>(Otway_.rocktype_).pd_;
             m_VG_          = get<9>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<9>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<9>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<9>(Otway_.rocktype_).Swi_pc_; // Swi = Swi_pc for this type
             TwoPhaseModel<dim>::snr_ = get<9>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<9>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<9>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<9>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<9>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 10: // H_Slt
             k_low_ = k_high_         = get<10>(Otway_.rocktype_).k_;
             phi_                     = get<10>(Otway_.rocktype_).phi_;
             pd_                      = get<10>(Otway_.rocktype_).pd_;
             m_VG_                    = get<10>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<10>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<10>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<10>(Otway_.rocktype_).Sgr_;
             krw_                     = get<10>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<10>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 11: // P_Slt_FSst: COMPOSITE 
             is_composite_ = true;
             L_low_ = get<11>(Otway_.rocktype_).LY_low_; L_high_ = get<11>(Otway_.rocktype_).LY_high_;
             k_low_ = get<11>(Otway_.rocktype_).k_low_;  k_high_ = get<11>(Otway_.rocktype_).k_high_;
             pd_high_       = get<11>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<11>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<11>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<11>(Otway_.rocktype_).m_low_ );
             phi_           = get<11>(Otway_.rocktype_).phi_; 
             pd_            = get<11>(Otway_.rocktype_).Pd( vt_magnitude_x_ ); 
             m_VG_          = get<11>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<11>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<11>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<11>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<11>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<11>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<11>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<11>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<11>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ ); 
           break;
         case 12: // P_Slt_CSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<12>(Otway_.rocktype_).LY_low_; L_high_ = get<12>(Otway_.rocktype_).LY_high_;
             k_low_ = get<12>(Otway_.rocktype_).k_low_;  k_high_ = get<12>(Otway_.rocktype_).k_high_;
             pd_high_       = get<12>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<12>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<12>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<12>(Otway_.rocktype_).m_low_ );
             phi_           = get<12>(Otway_.rocktype_).phi_; 
             pd_            = get<12>(Otway_.rocktype_).Pd( vt_magnitude_x_ );
             m_VG_          = get<12>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<12>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<12>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<12>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = 0.; // get<12>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<12>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<12>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<12>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<12>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ ); 
           break;
         case 13: // H_FSst
             k_low_ = k_high_ = get<13>(Otway_.rocktype_).k_;
             phi_                     = get<13>(Otway_.rocktype_).phi_;
             pd_                      = get<13>(Otway_.rocktype_).pd_;
             m_VG_                    = get<13>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<13>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<13>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<13>(Otway_.rocktype_).Sgr_;
             krw_                     = get<13>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<13>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 14: // X_CSst_FSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<14>(Otway_.rocktype_).LY_low_; L_high_ = get<14>(Otway_.rocktype_).LY_high_;
             k_low_ = get<14>(Otway_.rocktype_).k_low_;  k_high_ = get<14>(Otway_.rocktype_).k_high_;
             pd_high_       = get<14>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<14>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<14>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<14>(Otway_.rocktype_).m_low_ );
             phi_           = get<14>(Otway_.rocktype_).phi_; 
             pd_            = get<14>(Otway_.rocktype_).pd_;
             m_VG_          = get<14>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<14>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<14>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<14>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<14>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<14>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<14>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<14>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<14>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );  
           break;
         case 15: // H_CSst  
             k_low_ = k_high_ = get<15>(Otway_.rocktype_).k_;
             phi_                     = get<15>(Otway_.rocktype_).phi_;
             pd_                      = get<15>(Otway_.rocktype_).pd_;
             m_VG_                    = get<15>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<15>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<15>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<15>(Otway_.rocktype_).Sgr_;
             krw_                     = get<15>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<15>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 16: // BAFFLE
             k_low_ = k_high_ = get<16>(Otway_.rocktype_).k_;
             phi_                     = get<16>(Otway_.rocktype_).phi_;
             pd_                      = get<16>(Otway_.rocktype_).pd_;
             m_VG_                    = get<16>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<16>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<16>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<16>(Otway_.rocktype_).Sgr_;
             krw_                     = get<16>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<16>(Otway_.rocktype_).Krn(Sw_);
           break;
         default:
           cerr <<"\n\t rocktype: "<< rocktype_;
           csmp_error.notice( ERROR, "HeterogeneityAndRateAwareModel<dim>::Initialize", "rocktype not recognized");
      }
     
    // uses the poroperms discretised on the simulation model 
    const bool k_from_rocktypes(false);  
    InitialisePermeability( e, k_from_rocktypes );

    // determine permeability anisotropy: assumption: all laminations are horizontal
    if ( k_low_ != k_high_ ) {
         // permeability
         K_flow_direction_              = PermeabilityInFlowDirection( vt_normalised_ );
         // TwoPhaseModel<dim>::k_ stores the magnitude of the permeability in the x,z plane
         K_reduction_in_flow_direction_ = K_flow_direction_ / TwoPhaseModel<dim>::k_; 
      }

    // 2. Computing the force balances
    // -----------------------------------------------
    grad_p_magnitude_ = PressureGradientMagnitude(e); 
    Nc_               = Nc_kgradP_Version( grad_p_magnitude_ ); // capillary number

    const double lambda_t = TwoPhaseModel<dim>::TotalMobility();
    if ( isnan(lambda_t) ) {
         Out(1); // 1=wetting phase
         throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                                      "mobt = NaN.");
      }

    if ( lambda_t < 0. ) {
         Out(1);
         throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                                      "mobt is negative.");
      }
    if ( lambda_t <= numeric_limits<double>::epsilon() ) {
         Out(1);
         throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                                      "mobt is zero.");
      }

/*
 #ifdef DEBUG_HETEROGENEITY_AWARE_MODEL





// no flow = Ncap=0
VectorVariable<dim> vt; vt=0.;
// for all rocktypes
for ( long i=0; i<=15; i++ )
  WriteRelativePermeabilityTable( "Maartje2_layer_parallel_Nc0", i, vt );

// horizontal flow, Ncap=1.0e-6
vt(0) = 1.0e-4;
vt(1) = 0.;
for ( long i=0; i<=15; i++ )
  WriteRelativePermeabilityTable( "Maartje2_layer_parallel_Nc6", i, vt  );

// vertical flow, Ncap=1.0e-6
vt(0) = 0.;
vt(1) = 1.0e-4;
for ( long i=0; i<=15; i++ )
  WriteRelativePermeabilityTable( "Maartje2_layer_perpendicular_Nc6", i, vt );
  
// MISC tests
cerr <<"\npoly2C: "<< polyC2( 1.05, 11.73, 0.3316 );
cerr <<"\npow() with negative value: "<< pow( -0.05, 2 );    
  
throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                       "written relative permeability curves to file for testing.");
#endif
*/

 } // end Initialize



/**
      Intialises model for testing and plotting
*/
template<uint32_t dim>
void HeterogeneityAndRateAwareModel<dim>::Initialize( long rocktype, double Sw, const VectorVariable<dim>& vt )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // 1. base class parameter initialisation
    // --------------------------------------
    Sw_ = TwoPhaseModel<dim>::sat_ = Sw;
    TwoPhaseModel<dim>::mun_ = 0.000022959;
    TwoPhaseModel<dim>::muw_ = 0.00054843338222415;
    TwoPhaseModel<dim>::rhn_ = 283.;
    TwoPhaseModel<dim>::rhw_ = 992.;

    // 2. rocktype and flow velocity
    // ----------------------------------------------------------------
    rocktype_ = rocktype;
    assert( rocktype_ >= 0 && rocktype_ <= 16 );
    vt_ = vt_normalised_ = vt;
    vt_magnitude_        = vt.Length();
    vt_magnitude_x_        = fabs(vt[0]); 
    vt_magnitude_y_        = fabs(vt[1]);
    if ( vt_magnitude_ > numeric_limits<double>::epsilon() )
      vt_normalised_    /= vt_magnitude_;
    else {
         vt_normalised_(0) = 1.;
         vt_normalised_(1) = 0.;
      }
      
#ifdef TURN_OFF_RATE_AWARE
    vt_magnitude_ = vt_magnitude_x_ = vt_magnitude_y_ = 1.0e-12; //set to below capillary limit
#endif        
 
#ifdef DEBUG_HETEROGENEITY_AWARE_MODEL

  // no flow = Ncap=0
  VectorVariable<dim> vt; 
  vt=0.;

  //1. horizontal flows
  //1.1 vt(0) < CL
  vt(0) = 4.0e-8; vt(1) = 0.;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Horizontal_CL_4e-8", i, vt );
    
  //1.2 CL < vt(0) < VL: 3 velocities
  vt(0) = 5.0e-5; vt(1) = 0.;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Horizontal_CL_VL_5e-5", i, vt );
    
  vt(0) = 5.0e-4; vt(1) = 0.;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Horizontal_CL_VL_5e-4", i, vt );
    
  vt(0) = 5.0e-3; vt(1) = 0.;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Horizontal_CL_VL_5e-3", i, vt );
    
  //1.3 vt(0) > VL
  vt(0) = 6.0e-1; vt(1) = 0.;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Horizontal_VL_6e-1", i, vt );    

  //2. vertical flows
  //2.1 vt(1) < CL
  vt(0) = 0; vt(1) = 4.0e-8;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Vertical_CL_4e-8", i, vt );
    
  //2.2 CL < vt(1) < VL: 3 velocities
  vt(0) = 0.; vt(1) = 5.0e-6;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Vertical_CL_VL_5e-6", i, vt );
    
  vt(0) = 0.; vt(1) = 5.0e-5;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Vertical_CL_VL_5e-5", i, vt );
    
  vt(0) = 0.; vt(1) = 5.0e-4;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Vertical_CL_VL_5e-4", i, vt );
    
  //2.3 vt(1) > VL
  vt(0) = 0.; vt(1) = 6.0e-2;
  for ( long i=0; i<=16; i++ ) // for all rocktypes
      WriteRelativePermeabilityTable( "Vertical_VL_6e-2", i, vt );    
    
  throw csmp::Exception( INFO, "HeterogeneityAndRateAwareModel<dim>::Initialize:",
                         "written relative permeability curves to file for testing.");
#endif  
 
    // default: standard rocktype
    is_composite_ = false;
    TwoPhaseModel<dim>::tensor_permeability_ = false;
    TwoPhaseModel<dim>::ift_  = 0.035; // 35 mN/m water - CO2

    switch( rocktype_ ) {
         case 0: // WELL
             // reading single-petro-type rocktypes; scalar permeability is read from the elements
             k_low_ = k_high_ = get<0>(Otway_.rocktype_).k_;
             TwoPhaseModel<dim>::k_   = get<0>(Otway_.rocktype_).k_;
             phi_                     = get<0>(Otway_.rocktype_).phi_;
             pd_                      = get<0>(Otway_.rocktype_).pd_;
             m_VG_                    = get<0>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<0>(Otway_.rocktype_).Swi_pc_; // swr used in fitting the pc curve
             bcp_                     = 0.; // linear model
             // these are not used in any of the underlying functions of the 2-phase model
             TwoPhaseModel<dim>::swr_ = get<0>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<0>(Otway_.rocktype_).Sgr_;
             krw_                     = get<0>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<0>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 1: // H_Mst
             k_low_ = k_high_ = get<1>(Otway_.rocktype_).k_;
             TwoPhaseModel<dim>::k_   = get<1>(Otway_.rocktype_).k_;
             phi_                     = get<1>(Otway_.rocktype_).phi_;
             pd_                      = get<1>(Otway_.rocktype_).pd_;
             m_VG_                    = get<1>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<1>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<1>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<1>(Otway_.rocktype_).Sgr_;
             krw_                     = get<1>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<1>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 2: // M_CbSst_Mst (composite but not rate-dependent)
             is_composite_ = true;
             L_low_ = get<2>(Otway_.rocktype_).LY_low_; L_high_ = get<2>(Otway_.rocktype_).LY_high_;
             k_low_ = get<2>(Otway_.rocktype_).k_low_;  k_high_ = get<2>(Otway_.rocktype_).k_high_;
             pd_high_       = get<2>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<2>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<2>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<2>(Otway_.rocktype_).m_low_ );
             phi_           = get<2>(Otway_.rocktype_).phi_; 
             pd_            = get<2>(Otway_.rocktype_).pd_;
             m_VG_          = get<2>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<2>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<2>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<2>(Otway_.rocktype_).Swi_pc_; // for this type Swi=Swi_pc
             TwoPhaseModel<dim>::snr_ = get<2>(Otway_.rocktype_).Sgr_;
             // no rate dependence for this carbonate cemented rocktype
             krw_parallel_  = get<2>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<2>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<2>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<2>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 3: // H_CbSst
             k_low_ = k_high_ = get<3>(Otway_.rocktype_).k_;
             TwoPhaseModel<dim>::k_   = get<3>(Otway_.rocktype_).k_;
             phi_                     = get<3>(Otway_.rocktype_).phi_;
             pd_                      = get<3>(Otway_.rocktype_).pd_;
             m_VG_                    = get<3>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<3>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<3>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<3>(Otway_.rocktype_).Sgr_;
             krw_                     = get<3>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<3>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 4: // P_Mst_Slt: COMPOSITE
             is_composite_ = true;
             // layer-specific properties
             L_low_ = get<4>(Otway_.rocktype_).LY_low_; L_high_ = get<4>(Otway_.rocktype_).LY_high_;
             k_low_ = get<4>(Otway_.rocktype_).k_low_;  k_high_ = get<4>(Otway_.rocktype_).k_high_;
             pd_high_       = get<4>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<4>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<4>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<4>(Otway_.rocktype_).m_low_ );
             // composite properties
             phi_           = get<4>(Otway_.rocktype_).phi_; // average porosity
             pd_            = get<4>(Otway_.rocktype_).Pd( vt_magnitude_x_ );
             m_VG_          = get<4>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<4>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<4>(Otway_.rocktype_).dPc_;
             // saturation endpoints        
             TwoPhaseModel<dim>::swr_ = get<4>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<4>(Otway_.rocktype_).Sgr_; 
             // relative permeability parallel and perpendicular to layers
             krw_parallel_  = get<4>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<4>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<4>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<4>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );
           break;
         case 5: // M_CbSst_Slt: COMPOSITE - not rate dependent
             is_composite_ = true;
             L_low_ = get<5>(Otway_.rocktype_).LY_low_; L_high_ = get<5>(Otway_.rocktype_).LY_high_;
             k_low_ = get<5>(Otway_.rocktype_).k_low_;  k_high_ = get<5>(Otway_.rocktype_).k_high_;
             pd_high_       = get<5>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<5>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<5>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<5>(Otway_.rocktype_).m_low_ );
             phi_           = get<5>(Otway_.rocktype_).phi_; 
             pd_            = get<5>(Otway_.rocktype_).pd_;
             m_VG_          = get<5>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<5>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<5>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<5>(Otway_.rocktype_).Swi_pc_; // for this type Swi=Swi_pc
             TwoPhaseModel<dim>::snr_ = get<5>(Otway_.rocktype_).Sgr_;
             // no rate dependence for this carbonate cemented rocktype
             krw_parallel_  = get<5>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<5>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<5>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<5>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 6: // P_Mst_FSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<6>(Otway_.rocktype_).LY_low_; L_high_ = get<6>(Otway_.rocktype_).LY_high_;
             k_low_ = get<6>(Otway_.rocktype_).k_low_;  k_high_ = get<6>(Otway_.rocktype_).k_high_;
             pd_high_       = get<6>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<6>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<6>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<6>(Otway_.rocktype_).m_low_ );
             phi_           = get<6>(Otway_.rocktype_).phi_; 
             pd_            = get<6>(Otway_.rocktype_).Pd( vt_magnitude_x_ );
             m_VG_          = get<6>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<6>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<6>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<6>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<6>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<6>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<6>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<6>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<6>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ ); 
           break;
         case 7: // M_CbSst_FSst: COMPOSITE - not rate dependent
             is_composite_ = true;
             L_low_ = get<7>(Otway_.rocktype_).LY_low_; L_high_ = get<7>(Otway_.rocktype_).LY_high_;
             k_low_ = get<7>(Otway_.rocktype_).k_low_;  k_high_ = get<7>(Otway_.rocktype_).k_high_;
             pd_high_       = get<7>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<7>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<7>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<7>(Otway_.rocktype_).m_low_ );
             phi_           = get<7>(Otway_.rocktype_).phi_; 
             pd_            = get<7>(Otway_.rocktype_).pd_;
             m_VG_          = get<7>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<7>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<7>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<7>(Otway_.rocktype_).Swi_pc_; // Swi = Sw_pc
             TwoPhaseModel<dim>::snr_ = get<7>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<7>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<7>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<7>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<7>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 8: // P_Mst_CSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<8>(Otway_.rocktype_).LY_low_; L_high_ = get<8>(Otway_.rocktype_).LY_high_;
             k_low_ = get<8>(Otway_.rocktype_).k_low_;  k_high_ = get<8>(Otway_.rocktype_).k_high_;
             pd_high_       = get<8>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<8>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<8>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<8>(Otway_.rocktype_).m_low_ );
             phi_           = get<8>(Otway_.rocktype_).phi_; 
             pd_            = get<8>(Otway_.rocktype_).Pd( vt_magnitude_x_ );
             m_VG_          = get<8>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<8>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<8>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<8>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<8>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<8>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<8>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<8>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<8>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );  
           break;
         case 9: // M_CbSst_CSst: COMPOSITE - not rate dependent 
             is_composite_ = true;
             L_low_ = get<9>(Otway_.rocktype_).LY_low_; L_high_ = get<9>(Otway_.rocktype_).LY_high_;
             k_low_ = get<9>(Otway_.rocktype_).k_low_;  k_high_ = get<9>(Otway_.rocktype_).k_high_;
             pd_high_       = get<9>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<9>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<9>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<9>(Otway_.rocktype_).m_low_ );
             phi_           = get<9>(Otway_.rocktype_).phi_; 
             pd_            = get<9>(Otway_.rocktype_).pd_;
             m_VG_          = get<9>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<9>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<9>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<9>(Otway_.rocktype_).Swi_pc_; // Swi = Swi_pc for this type
             TwoPhaseModel<dim>::snr_ = get<9>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<9>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_ );
             krw_crossflow_ = get<9>(Otway_.rocktype_).Krw_CrossDrainage( Sw_ );
             krn_parallel_  = get<9>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_ );
             krn_crossflow_ = get<9>(Otway_.rocktype_).Krn_CrossDrainage( Sw_ );
           break;
         case 10: // H_Slt
             k_low_ = k_high_ = get<10>(Otway_.rocktype_).k_;
             TwoPhaseModel<dim>::k_   = get<10>(Otway_.rocktype_).k_;
             phi_                     = get<10>(Otway_.rocktype_).phi_;
             pd_                      = get<10>(Otway_.rocktype_).pd_;
             m_VG_                    = get<10>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<10>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<10>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<10>(Otway_.rocktype_).Sgr_;
             krw_                     = get<10>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<10>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 11: // P_Slt_FSst: COMPOSITE 
             is_composite_ = true;
             L_low_ = get<11>(Otway_.rocktype_).LY_low_; L_high_ = get<11>(Otway_.rocktype_).LY_high_;
             k_low_ = get<11>(Otway_.rocktype_).k_low_;  k_high_ = get<11>(Otway_.rocktype_).k_high_;
             pd_high_       = get<11>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<11>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<11>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<11>(Otway_.rocktype_).m_low_ );
             phi_           = get<11>(Otway_.rocktype_).phi_; 
             pd_            = get<11>(Otway_.rocktype_).Pd( vt_magnitude_x_ ); 
             m_VG_          = get<11>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<11>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<11>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<11>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<11>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<11>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<11>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<11>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<11>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ ); 
           break;
         case 12: // P_Slt_CSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<12>(Otway_.rocktype_).LY_low_; L_high_ = get<12>(Otway_.rocktype_).LY_high_;
             k_low_ = get<12>(Otway_.rocktype_).k_low_;  k_high_ = get<12>(Otway_.rocktype_).k_high_;
             pd_high_       = get<12>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<12>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<12>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<12>(Otway_.rocktype_).m_low_ );
             phi_           = get<12>(Otway_.rocktype_).phi_; 
             pd_            = get<12>(Otway_.rocktype_).Pd( vt_magnitude_x_ );
             m_VG_          = get<12>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<12>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<12>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<12>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = 0.; // get<12>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<12>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<12>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<12>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<12>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );  
           break;
         case 13: // H_FSst
             k_low_ = k_high_ = get<13>(Otway_.rocktype_).k_;
             TwoPhaseModel<dim>::k_   = get<13>(Otway_.rocktype_).k_;
             phi_                     = get<13>(Otway_.rocktype_).phi_;
             pd_                      = get<13>(Otway_.rocktype_).pd_;
             m_VG_                    = get<13>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<13>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<13>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<13>(Otway_.rocktype_).Sgr_;
             krw_                     = get<13>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<13>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 14: // X_CSst_FSst: COMPOSITE
             is_composite_ = true;
             L_low_ = get<14>(Otway_.rocktype_).LY_low_; L_high_ = get<14>(Otway_.rocktype_).LY_high_;
             k_low_ = get<14>(Otway_.rocktype_).k_low_;  k_high_ = get<14>(Otway_.rocktype_).k_high_;
             pd_high_       = get<14>(Otway_.rocktype_).pd_high_;
             pd_low_        = get<14>(Otway_.rocktype_).pd_low_;
             bcp_high_      = lambdaFrom_VG( get<14>(Otway_.rocktype_).m_high_ );
             bcp_low_       = lambdaFrom_VG( get<14>(Otway_.rocktype_).m_low_ );
             phi_           = get<14>(Otway_.rocktype_).phi_; 
             pd_            = get<14>(Otway_.rocktype_).pd_;
             m_VG_          = get<14>(Otway_.rocktype_).m_; 
             Swi_pc_        = get<14>(Otway_.rocktype_).Swi_pc_;
             dPc_           = get<14>(Otway_.rocktype_).dPc_;
             TwoPhaseModel<dim>::swr_ = get<14>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<14>(Otway_.rocktype_).Sgr_;
             krw_parallel_  = get<14>(Otway_.rocktype_).Krw_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krw_crossflow_ = get<14>(Otway_.rocktype_).Krw_CrossDrainage( Sw_, vt_magnitude_y_ );
             krn_parallel_  = get<14>(Otway_.rocktype_).Krn_ParallelDrainage( Sw_, vt_magnitude_x_ );
             krn_crossflow_ = get<14>(Otway_.rocktype_).Krn_CrossDrainage( Sw_, vt_magnitude_y_ );  
           break;
         case 15: // H_CSst  
             k_low_ = k_high_ = get<15>(Otway_.rocktype_).k_;
             TwoPhaseModel<dim>::k_   = get<15>(Otway_.rocktype_).k_;
             phi_                     = get<15>(Otway_.rocktype_).phi_;
             pd_                      = get<15>(Otway_.rocktype_).pd_;
             m_VG_                    = get<15>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<15>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<15>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<15>(Otway_.rocktype_).Sgr_;
             krw_                     = get<15>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<15>(Otway_.rocktype_).Krn(Sw_);
           break;
         case 16: // BAFFLE
             k_low_ = k_high_ = get<16>(Otway_.rocktype_).k_;
             phi_                     = get<16>(Otway_.rocktype_).phi_;
             pd_                      = get<16>(Otway_.rocktype_).pd_;
             m_VG_                    = get<16>(Otway_.rocktype_).m_; 
             Swi_pc_                  = get<16>(Otway_.rocktype_).Swi_pc_;
             bcp_                     = lambdaFrom_VG( m_VG_);
             TwoPhaseModel<dim>::swr_ = get<16>(Otway_.rocktype_).Swi_;
             TwoPhaseModel<dim>::snr_ = get<16>(Otway_.rocktype_).Sgr_;
             krw_                     = get<16>(Otway_.rocktype_).Krw(Sw_);
             krn_                     = get<16>(Otway_.rocktype_).Krn(Sw_);
           break;         
         default:
           cerr <<"\n\t rocktype: "<< rocktype_;
           csmp_error.notice( ERROR, "HeterogeneityAndRateAwareModel<dim>::Initialize", "rocktype not recognized");
      }
   
    // determine permeability anisotropy TODO: here we assume that laminations are horizontal
    if ( is_composite_ ) {
         // permeability
         // TwoPhaseModel<dim>::tensor_permeability_ = true; // anisotropy is not communicated to base class
         const double k_crossflow_ = PermeabilityPerpendicularToLaminations();
         const double k_parallel_  = PermeabilityParallelToLaminations();
         // assuming that layers are horizontal and that the stored K is the horizontal one
         TwoPhaseModel<dim>::K_ = 0.;
         TwoPhaseModel<dim>::K_(0,0)    = k_parallel_;
         TwoPhaseModel<dim>::K_(1,1)    = k_crossflow_;
         K_flow_direction_ = PermeabilityInFlowDirection( vt_normalised_ );
         TwoPhaseModel<dim>::k_         = k_parallel_;
         K_reduction_in_flow_direction_ = K_flow_direction_ / k_parallel_; 
      }

} // end Initialise (testing & plotting)




/**
       Initialises K_ in TwoPhaseModel base class (2D version)
       @attention TwoPhaseModel<3U>::k_   is used to store the magnitude of the maximum permeability in the horizontal plane
*/
template<>
void HeterogeneityAndRateAwareModel<2U>::InitialisePermeability( const Element<2U>& e, bool k_from_rocktypes )
 {
     TwoPhaseModel<2U>::K_ = 0.;
     if ( k_from_rocktypes ) {
          assert( k_low_  >= 1.0e-18 ); 
          assert( k_high_ >= 1.0e-16 ); 
          TwoPhaseModel<2U>::K_(0,0) = k_high_;
          TwoPhaseModel<2U>::K_(1,1) = k_low_; 
          TwoPhaseModel<2U>::k_      = k_high_;      
          return;
       }
     TwoPhaseModel<2U>::K_(0,0) = e.Read( kh_key_ );
     TwoPhaseModel<2U>::K_(1,1) = e.Read( kv_key_ );
     TwoPhaseModel<2U>::k_      = TwoPhaseModel<2U>::K_(0,0);      
 } 


/**
       Initialises K_ in TwoPhaseModel base class (3D version)
       @attention TwoPhaseModel<3U>::k_   is used to store the magnitude of the maximum permeability in the horizontal plane
*/
template<>
void HeterogeneityAndRateAwareModel<3U>::InitialisePermeability( const Element<3U>& e, bool k_from_rocktypes )
 {
     TwoPhaseModel<3U>::K_ = 0.;
     if ( k_from_rocktypes ) {
          assert( k_low_  >= 1.0e-18 ); 
          assert( k_high_ >= 1.0e-16 ); 
          TwoPhaseModel<3U>::K_(0,0) = k_high_; // x
          TwoPhaseModel<3U>::K_(1,1) = k_low_;  // y = vertical        
          TwoPhaseModel<3U>::K_(2,2) = k_high_; // z = horizontal
          TwoPhaseModel<3U>::k_      = k_high_;      
          return;
       }
     TwoPhaseModel<3U>::K_(0,0) = e.Read( kh_key_ );
     TwoPhaseModel<3U>::K_(1,1) = e.Read( kv_key_ );
     TwoPhaseModel<3U>::K_(2,2) = TwoPhaseModel<3U>::K_(0,0);
     TwoPhaseModel<3U>::k_      = sqrt( TwoPhaseModel<3U>::K_(0,0)*TwoPhaseModel<3U>::K_(0,0) +
                                        TwoPhaseModel<3U>::K_(2,2)*TwoPhaseModel<3U>::K_(2,2) );      
 } 

// 1D stub
template<>
void HeterogeneityAndRateAwareModel<1U>::InitialisePermeability( const Element<1U>& e, bool k_from_rocktypes )
 {
     TwoPhaseModel<1U>::K_ = 0.;
     if ( k_from_rocktypes) {
          assert( k_low_  >= 1.0e-18 ); 
          assert( k_high_ >= 1.0e-16 ); 
          TwoPhaseModel<1U>::K_(0,0) = k_high_;
          return;
       }
     TwoPhaseModel<1U>::K_(0,0) = TwoPhaseModel<1U>::k_ = e.Read( kh_key_ );
 } 




/// weighted permeability average for flow along layers
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::PermeabilityParallelToLaminations() const
 {
    assert( is_composite_ );
    return (L_high_ * k_high_ + L_low_ * k_low_) / (L_high_ + L_low_);
 }
 
 
 
/// weigthed harmonic mean of permeability for flow across layers
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::PermeabilityPerpendicularToLaminations() const
 {
    assert( is_composite_ );
    const double sum_of_weights = L_low_ + L_high_;
    return sum_of_weights / (L_high_ / k_high_ + L_low_ / k_low_);
 }



/**
    From the permeability tensor (in TwoPhaseModel) and the flow direction, this method calculates and returns the permeability in the flow direction.
    
        @attention the normalized velocity vt_normalised must be initiaised.
        
        @attention if the flow velocity is zero, horizontal flow is assumed.
    
*/
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::PermeabilityInFlowDirection( const VectorVariable<dim>& vt_normalised ) const
 {
    assert( is_composite_ );
    assert( TwoPhaseModel<dim>::K_(0,0) > 0. );
    assert( TwoPhaseModel<dim>::K_(1,1) > 0. );
    
    if ( fabs(vt_normalised.Length() - 1.) <= numeric_limits<double>::epsilon() )
      // horizontal permeability
      return TwoPhaseModel<dim>::K_(0,0);

    // finding the permeability in the direction of the velocity vector
    vc_ = TwoPhaseModel<dim>::K_ * vt_normalised;

    // the length of vc_ is equal to the magnitude of permeability in the target direction
    return vc_.Length();
   
 } // end PermeabilityInFlowDirection



/**
    Computes prominent direction of the flow of the fluid mixture.
 
    Uses total velocity, vt, to establish whether vertical component of flow is greater than horizontal one.
 
    @note if the flow velocity magnitude is below threshold value, HORIZONAL is returned.
*/
template<uint32_t dim>
typename HeterogeneityAndRateAwareModel<dim>::FLOW_DIRECTION 
HeterogeneityAndRateAwareModel<dim>::ProminentFlowDirection( const VectorVariable<dim>& vt ) const
 {
     if ( fabs(vt[0])+fabs(vt[1]) <= numeric_limits<double>::epsilon() ) return HORIZONTAL;
     // comparing the vertical component with the horizontal magnitude of the flow
     const double horizontal_magnitude = ( dim == 2U ) ? vt[0] : sqrt( vt[0]*vt[0] + vt[2]*vt[2] );
     FLOW_DIRECTION direction = ( fabs(vt[1]) > horizontal_magnitude ) ? VERTICAL : HORIZONTAL;
     return direction;
 }




/**
    Returns the magnitude of the pressure gradient corrected for the hydrostatic component which causes no flow
*/
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::PressureGradientMagnitude( const Element<dim>& e ) const
 {
    array<double,dim> gradP = {0.}; // the unspecified elements are initialised to zero
    e.dN_AtBaryCenter( DN_ );
    for ( auto i{0U}; i<e.Nodes(); ++i ) {
         const double pf = e.N(i)->Read( pf_key_ );
         for ( size_t j{0U}; j<dim; ++j )
         gradP[j] += DN_(j,i) * pf;
      }
    // elimination of the hydrostatic pressure gradient
    gradP[dim-1U] -= TwoPhaseModel<dim>::acc_gravity_ * (Sw_*TwoPhaseModel<dim>::rhw_ + (1.-Sw_)*TwoPhaseModel<dim>::rhn_);
    // magnitude of the reduced pressure gradient
    double grad_p_magnitude(0.);
    for ( size_t j{0U}; j<dim; ++j ) grad_p_magnitude += gradP[j] * gradP[j];
    return sqrt(grad_p_magnitude);

} // end PressureGradientMagnitude



/**
     pressure gradient form: Nc = k ||grad p|| / sigma; 
     @attention uses horizontal permeability only
*/
template<uint32_t dim>
double csmp::HeterogeneityAndRateAwareModel<dim>::Nc_kgradP_Version( double pf_gradient_magnitude ) const
 {
     assert( pf_gradient_magnitude != UNSPECIFIED );
     if ( !is_composite_ ) return TwoPhaseModel<dim>::K_(0,0) * pf_gradient_magnitude / TwoPhaseModel<dim>::ift_;
     return PermeabilityParallelToLaminations() * pf_gradient_magnitude / TwoPhaseModel<dim>::ift_;
 }




/**
    SKM: thickness-weighted average of the effective saturations in the high and the low k laminations.
 
    @attention must always be called first because Seff must be communicated with the base class.
    @attention this method just transfers the saturation because the effective saturation is calculated later
*/
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::EffectiveSaturation() const
 {
    TwoPhaseModel<dim>::sat_  = Sw_;
    TwoPhaseModel<dim>::seff_ = Sw_;
    return seffL(Sw_,Swi_pc_);
   
 } // end EffectiveSaturation





/**
    Agreed weighting procedure in which the flow-direction based average of the relative permeability is weighted
    lumping in the Kv/Kh ratio.
    The result is limited to greater than or equal to zero.
*/
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::krw_Phase() const
 {
    if ( !is_composite_ ) return krw_; // krw_VG( Sw_, m_VG_ );
   
    // to get the ensemble krw for the composite, the parallel and perpendicular values are blended
    // taking into account the flow direction 
    // --------------------------------------
    const double krw = krw_parallel_ * (vt_normalised_[0]*vt_normalised_[0]) + krw_crossflow_ * (vt_normalised_[1]*vt_normalised_[1]);

    // scaling the relative permeability by the vertical permeability
    return max( krw * K_reduction_in_flow_direction_, 0. );
   
 } // end krw_Phase




/**
    CO2 relative permeability
*/
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::krn_Phase() const
 {
     if ( !is_composite_ ) return krn_; // krn_BC( Sw_ );

     const double krn = krn_parallel_ * (vt_normalised_[0]*vt_normalised_[0]) + krn_crossflow_ * (vt_normalised_[1]*vt_normalised_[1]);

     // scaling the relative permeability by the vertical permeability
     return max( krn * K_reduction_in_flow_direction_, 0. );

} // end krn_Phase









/// for the wetting phase
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::MaxFractionalFlowDerivative() const
 {
    return static_cast<double>(5.6); // as computed with dfds method
 }




/**

Not directionally dependent.

@attention pc covers the full saturation range, pc(sw) is capped based on maximum dpcds.

@attention this pc function takes into account the effect of flow rate, which also shifts the irreducible saturations.

*/
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::pc_Phase() const
 {
    //rocktype 16 uses BC
    if (rocktype_==16) {
        // if the BC parameter is set to zero, a linear capillary pressure model should be used
        if ( bcp_ <= numeric_limits<double>::epsilon() ) {
             return pd_ + ( 1. - Sw_) * (TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ - pd_); 
          }
        return min( pc_BC( Sw_, Swi_pc_, pd_, bcp_ ), TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_ );    
     }
    
    const double Pc_MAX(1.0e+7);
    const double PC_LOW_SW_LIMIT  = Swi_pc_ + 0.01;
    const double PC_HIGH_SW_LIMIT = 1.0 - 0.01;
    // linear regularization for lower part of sw range
    if ( Sw_ < PC_LOW_SW_LIMIT ) {
        const double pc_lim = pc_Phase_at( PC_LOW_SW_LIMIT );
        const double seff_mult( 1.0/ (1.0 - Swi_pc_ ) );
        const double dpcds_lim = dpcds_Phase_at( PC_LOW_SW_LIMIT )/seff_mult;
        const double pc = pc_lim + dpcds_lim*( Sw_ - PC_LOW_SW_LIMIT );
        return std::min(pc, Pc_MAX);    
    
      }
    // linear regularization for higher part of sw range
    else if ( Sw_ > PC_HIGH_SW_LIMIT ) {
        const double pc_lim = pc_Phase_at( PC_HIGH_SW_LIMIT );
        const double dpcds_lim ( (0.0 - pc_lim)/( 1.0 - PC_HIGH_SW_LIMIT ) );
        const double pc = dpcds_lim*( Sw_ - 1.0 );
        return std::max(pc, 0.);
     } 

    return std::min( pc_VG( Sw_, pd_, m_VG_, Swi_pc_ ), Pc_MAX);       

 } // end pc_Phase



template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::pc_Phase_at(double sw) const
{
    return pc_VG( sw, pd_, m_VG_, Swi_pc_ );   
}



/**
    Capillary pressure - saturation derivative estimated numerically.
    
    @attention the capillary pressure saturation derivative is capped to a minimum value of 10 bars/m to limit slowdown of simulator during pc spreading.
    
    @note this method is never called by Shao's FECFVM simulator as it only uses the absolute values of pc.
 
    @todo dpcdsw is zero as soon as one of the phases is immobile
    to avoid that CFL is influenced while capillary spreading is not possible.
*/
template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::dpcds_Phase() const
 {
    assert( Sw_ >= 0. );
    assert( Sw_ <= 1. );
    
    //use constant derivative at PC_LOW_SW_LIMIT for lower part of sw range
    const double PC_LOW_SW_LIMIT   = Swi_pc_ + 0.01;
    const double MIN_PC_DERIVATIVE = -1.0e6;    
    
    if(Sw_ < PC_LOW_SW_LIMIT) return max( MIN_PC_DERIVATIVE, dpcds_Phase_at(PC_LOW_SW_LIMIT ) );
        
    const double  h(0.01), sw(Sw_); // backing up the water saturation

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         Sw_ = 1.;
         double pc1 = pc_Phase();
         Sw_ = 1. - h;
         double pc2 = pc_Phase();
         return max( MIN_PC_DERIVATIVE, (pc1 - pc2) / h );
      }
  
    // low water saturation
    if ( sw <= h ) {
         Sw_ = h;
         double pc1 = pc_Phase();
         Sw_ = 0.;
         double pc2 = pc_Phase();
         return max( MIN_PC_DERIVATIVE, (pc1 - pc2) / h );
      }

    // water saturation between the endpoints
    Sw_ = sw + h;
    double pc1 = pc_Phase();
    Sw_ = sw - h;
    double pc2 = pc_Phase();
    // resetting sw value
    Sw_ = sw;
   
    return max( MIN_PC_DERIVATIVE, (pc1 - pc2) / (2. * h) );
}



template<uint32_t dim>
double HeterogeneityAndRateAwareModel<dim>::dpcds_Phase_at(double sw) const
 {
    assert( sw >= 0. );
    assert( sw <= 1. );
    const double  h(0.005);

    // more common case of a high water saturation first
    if ( sw >= (1. - h) ) {
         double pc1 = pc_Phase_at(1.);
         double pc2 = pc_Phase_at(1. - h);
         return (pc1 - pc2) / h;
      }
  
    // low water saturation
    if ( sw <= h ) {
         double pc1 = pc_Phase_at(h);
         double pc2 = pc_Phase_at(0.);
         return (pc1 - pc2) / h;
      }

    // water saturation between the endpoints
    double pc1 = pc_Phase_at(sw + h);
    double pc2 = pc_Phase_at(sw - h);
   
    return (pc1 - pc2) / (2. * h);
}





/**
    writes textfile with sw, krw(sw,Nc), krn(sw,Nc), and pc(sw) values computed for (composite) rocktype in 0.05 saturation increments
*/
template<uint32_t dim>
void HeterogeneityAndRateAwareModel<dim>::WriteRelativePermeabilityTable( const char* filename, 
                                                                          long RT, 
                                                                          const VectorVariable<dim>& vt )
 {
    // computing the capillary number from V, mu & k, assuming mu = muw at test conditions
    const double dyn_visc(0.00054843338222415),
                   Nc = (dyn_visc * vt.Length()) / TwoPhaseModel<dim>::ift_;

    std::string text_file(filename);
    text_file += "_RT";
    text_file += to_string(RT);
    text_file += "_Ncap";
    text_file += to_string(Nc);

    ofstream  ofs( string(text_file) + ".txt" );
   
    const string   rocktype(parseRockType(RT));
    const double original_sw(Sw_);
    
    ofs <<"rocktype="<< rocktype <<"="<< RT <<"\n";
    ofs <<"sw\t krw(sw,Nc="<< Nc <<"),\t krnw(sw,Nc) \t pc(sw,Nc)\n";
    //for ( double sw(0.); sw<1.05; sw+=0.05 )
    for ( double sw(0.); sw<1.005; sw+=0.005 )
      {
         // dynamic parameters
         Initialize(RT,sw,vt);
         assert( !isnan(Sw_) );
         ofs << sw << "\t"<< krw_Phase();
         ofs <<"\t"<< krn_Phase();
         ofs <<"\t"<< pc_Phase();
         ofs << endl;
      }
   
    // restoring current saturation
    Sw_ = original_sw;

 } // end WriteRelativePermeabilityTable
  


/**
    Viscous - capillary force balance (RVC)
 
    Uses the capillary number and the interfacial tension as input parameters.

    RVC = ((RVC_high_k * (L_high)) + (RVC_low_k * (L_low)))/(total thicknes=L);

*/
template<uint32_t dim>
double csmp::HeterogeneityAndRateAwareModel<dim>::RVC( double Ncap ) const
 {
    const double rvc_low_k  = ((Ncap * TwoPhaseModel<dim>::ift_) * L_low_) / (dPc_ * k_low_);
    const double rvc_high_k = ((Ncap * TwoPhaseModel<dim>::ift_) * L_high_) / (dPc_ * k_high_);
    // 
    return (rvc_high_k * L_high_ + rvc_low_k * L_low_) / (L_low_ + L_high_);
 }




template<uint32_t dim>
void HeterogeneityAndRateAwareModel<dim>::Out( uint32_t phase ) const
 {
    TwoPhaseModel<dim>::Out(phase);
    cout <<"\nHeterogeneityAndRateAwareModel<"<< dim << ">::Out(rocktype="<< rocktype_ <<"): return values of functions: "<< endl;
    cout <<"\nelement properties:";
    cout <<"\n                                 velocity (m/s): "<< scientific << vt_;
    if ( ProminentFlowDirection( vt_ ) == HORIZONTAL ) 
      cout <<"- dominantly horizontal flow.";
    else 
      cout <<"- dominantly vertical flow.";
    cout <<"\n                           capillary number, Nc: "<< scientific << Nc_;
    if ( is_composite_ ) {
         cout <<"\n      average water saturation in composite, sw: "<< defaultfloat << Sw_;
         cout <<"\n               layer-parallel permeability (m2): "<< scientific << PermeabilityParallelToLaminations();
         cout <<"\n          layer-perpendicular permeability (m2): "<< scientific << PermeabilityPerpendicularToLaminations();
         cout <<"\nratio between viscous and capillary forces, RVC: "<< defaultfloat << RVC( Nc_ );
      }
    else {
        cout <<"\n       average water saturation in rocktype, sw: "<< defaultfloat << Sw_;
        cout <<"\n                              permeability (m2): "<< scientific << TwoPhaseModel<dim>::K_;
      }
    cout <<"\n                                        krw(sw): "<< scientific << krw_Phase();
    cout <<"\n                                        krn(sw): "<< scientific << krn_Phase();
    cout <<"\n                                         pc(sw): "<< scientific << pc_Phase();
    cout <<"\n                                     dpc/ds_max: "<< scientific << TwoPhaseModel<dim>::MAX_CAPILLARY_PRESSURE_SLOPE_ << endl << endl;

 } // end Out
 
 


template class HeterogeneityAndRateAwareModel<1U>;
template class HeterogeneityAndRateAwareModel<2U>;
template class HeterogeneityAndRateAwareModel<3U>;

} // end namespace csmp


























