#ifndef HETEROGENEITY_AND_RATE_AWARE_MODEL_H
#define HETEROGENEITY_AND_RATE_AWARE_MODEL_H

#include "TwoPhaseModel.h"
//#include "OtwayCRC3_RockTypes.h"
#include "OtwayCRC3_RockTypes_Version_2.h" 

namespace csmp {

/**
    Maartje Boons & Sally Benson (Stanford visit May 16-23, 2019),
    model for a layered composite. Material parameters are inferred from rocktype.
    The composite is a dual of averages for the high-k laminations and the w-k ones.
 
    Model Concept
 
    - compute the flow rate dependent saturations in the sublayers of the composite
      - end-member cases:
        - viscous limit (VL) where all saturations are the same, i.e. the saturation is uniform and equal to that of the cell
        - CL where pc(sw_cell_average) is the same everywhere, but the saturations in the layers reflect their pc curves
      these limiting cases are not directionaly dependent
    - compute the relperms, layer by layer
    - apply averaging to find the horizontal (weighted mean) and vertical (harmonic mean) relative permeability pairs
    - since this model does drainage only, relperm effects of spontaneous imbibition can be ignored
 
    Implementation
 
    - extended curve fitting and interpolation between the curves
    - there are only 2 flow directions: horizontal and vertical, switch occurs at 45o
    - the anisotropy of the absolute permeability, is treated by scaling the vertical relative permeabilities with Kv/Kh ratio
    - capillary pressure is also averaged, but for the vertical CL case a Brooks-Corey model is used,
      parameterised with the properties of the low permeability layer.
    - for primary drainage, the capillary entry pressure, pd of the lowest permeability layer is used
 
    see Maple worksheet 'BHP/Stanford-visit/IMPLEMENTATION/HeterogeneityAwareSaturationFunctions'
 
    @note thus far, this is a drainage only model
    @note an alternative capillary number: k ||grad p|| / sigma   is used to get around the viscosity problem
 
    @todo MUST WE BLEND BETWEEN HORIZONTAL AND VERTICAL FLOW?
 
    @todo: perhaps define irreducible water saturation and use it in relperm model
    @todo: check how this can be handled efficiently for a suite of composite rocktypes
    @todo: think about what parameters to keep rather than using functions to compute them on the fly when needed
 
    @todo: change FlowRateDependentLayerSaturations that it deals with horizontal and vertical flow and all funky cases

*/
template<size_t dim>
class HeterogeneityAndRateAwareModel : public TwoPhaseModel<dim> {
  public:
    enum FLOW_DIRECTION { HORIZONTAL, VERTICAL } flow_direction_;
  public:
    HeterogeneityAndRateAwareModel( const PropertyDatabase<dim>& database,
                                    const bool sw_ro_mu_placement = false ); // NODE=true ELEMENT=false
  
    virtual ~HeterogeneityAndRateAwareModel();
    
    /// return integer key of rocktype
    int32 RockType( const Element<dim>& ) const;
  
    /// reads vt and makes some initiaisations
    void InitializeVelocity( const Element<dim>& );
    /**
        Computes all relevant values for the single-valued average saturation inside of the element.
        Dependent on flow direction and kv/kh ratio, averaging and scaling is done later inside of the saturation functions. 
    */
    virtual void Initialize( const Element<dim>& );

    /// updating model parameters for testing and plotting
    void Initialize( long rocktype, double64 Sw, const VectorVariable<dim>& vt );

    /// average effective saturation as required by 2-phase model
    virtual double64 EffectiveSaturation() const;
  
   // relative permeabilities
    /// heterogeneity-aware, rate-dependent versions, Nc is calculated in the background
    virtual double64 krw_Phase() const;
    virtual double64 krn_Phase() const;
  
    /// capillary pressure of the non-wetting phase; cap value is applied
    virtual double64 pc_Phase() const;
    double64 pc_Phase_at(double64 sw) const; 
  
    /// numeric implementation of capillary pressure derivative; cap value is applied
    virtual double64 dpcds_Phase() const;
    double64 dpcds_Phase_at(double64 sw) const;
    
    // maximum absolute value returned by dfdS
    virtual double64 MaxFractionalFlowDerivative() const;

    /// Outputs the relperms for the current model initialisation to a plot for visual examination; filename appends rocktype and Ncap calculated
    void WriteRelativePermeabilityTable( const char* filename, long rocktype, const VectorVariable<dim>& vt );
  
    virtual void Out( size_t phase ) const;
  
  private:
    /// takes permeability values (from Model or Rocktypes) and initialises scalar permeability k and tensor K in TwoPhaseModel base class
    void InitialisePermeability( const Element<dim>& e, bool k_from_rocktypes );
    /// weighted average
    double64  PermeabilityParallelToLaminations() const;
    /// harmonic mean
    double64  PermeabilityPerpendicularToLaminations() const;
    /// kv, kh tensor decomposition
    double64  PermeabilityInFlowDirection( const VectorVariable<dim>& normalised_mixture_velocity ) const;
    /// Prominent direction of flow
    FLOW_DIRECTION ProminentFlowDirection( const VectorVariable<dim>& vt ) const;
    /// volume averaged irreducible water saturation
    double64  Swr_Composite() const;
    /// Magnitude of the pressure gradient
    double64  PressureGradientMagnitude( const Element<dim>& ) const;
    /// Capillary number, pressure gradient form: Nc = k ||grad p|| / sigma
    double64  Nc_kgradP_Version( double64 pf_gradient_magnitude ) const;
    /// Ratio between viscous and capillary forces
    double64  RVC( double64 Ncap ) const;
     
    // TODO: add function that assesses whether we are dealing with imbibition or drainage
  
  private:
    const csmp::Index  RRT_key_,         ///<  reservoir rock type
                       pf_key_,          ///<  (absolute) fluid pressure
                       kh_key_, kv_key_, ///<  horizontal and vertical permeability values
                       phi_key_,         ///< porosity (only used to write values if the RT values are to be used)
                       vt_key_;          ///<  total velocity of the fluid mixture
  
    mutable DenseMatrix<DM_MIN>  DN_; ///<  for gradient computations
    VectorVariable<dim>          vt_; ///<  velocity
    VectorVariable<dim>          vt_normalised_; ///< to unit length
    mutable VectorVariable<dim>  vc_; ///< for all kinds of purposes
    
    // dynamic variables
    mutable double64 Sw_;

    // composite is modelled as a dual of two rock types
    OtwayRockTypes Otway_;
    long     rocktype_;
    bool     is_composite_;
    double64 K_flow_direction_, k_low_, k_high_, K_reduction_in_flow_direction_,
             L_low_, L_high_,                ///< permeability in flow direction; smallest over highest permeability
             vt_magnitude_,
             vt_magnitude_x_, vt_magnitude_y_, 
             krw_, krn_,                     ///< for standard rocktypes 
             krw_parallel_, krw_crossflow_,  ///< for composites 
             krn_parallel_, krn_crossflow_,
             phi_, pd_, pd_high_, pd_low_,
             m_VG_, bcp_, bcp_high_, bcp_low_, Swi_pc_, dPc_,
             pd_flow_direction_, pc_flow_direction_;         

    // derived quantities
    double64  grad_p_magnitude_; ///< magnitude of the fluid pressure gradient reduced by hydrostatic gradient
    double64  Nc_; // capillary mumber
};

} // end namespace csmp

#endif
