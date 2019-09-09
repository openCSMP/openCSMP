#ifndef HETEROGENEITY_AND_RATE_AWARE_MODEL_H
#define HETEROGENEITY_AND_RATE_AWARE_MODEL_H

#include "TwoPhaseModel.h"

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
                                    const char* rocktype, const char* total_velocity,
                                    const char* pc_entry,
                                    const bool sw_ro_mu_placement = false ); // NODE=true ELEMENT=false
  
    virtual ~HeterogeneityAndRateAwareModel();
  
    virtual void Initialize( const Element<dim>& e );

    /// average effective saturation as required by 2-phase model
    virtual double64  EffectiveSaturation() const;
  
   // relative permeabilities
    /// heterogeneity-aware, rate-dependent version, Nc is calculated in the background
    virtual double64 krw_Phase() const;
    virtual double64 krn_Phase() const;
  
    /// capillary pressure of the non-wetting phase
    virtual double64 pc_Phase( ) const;
  
    /// numeric implementation of capillary pressure derivative
    virtual double64 dpcds_Phase() const;
    
    // maximum absolute value returned by dfdS
    virtual double64 MaxFractionalFlowDerivative() const;

    /// Outputs the relperms for the current model initialisation to a plot for visual examination
    void WriteRelativePermeabilityTable( const char* filename, double64 Ncap );
  
    virtual void Out( size_t phase ) const;
  
  private:
    /// weighted average
    double64 PermeabilityParallelToLaminations() const;
    /// harmonic mean
    double64 PermeabilityPerpendicularToLaminations() const;
    /// prominent direction of flow, determined from vt
    HeterogeneityAndRateAwareModel::FLOW_DIRECTION ProminentFlowDirection() const;
    /// volume averaged irreducible water saturation
    double64 Swr_Composite() const;
    /// Magnitude of the pressure gradient
    double64 PressureGradientMagnitude( const Element<dim>& ) const;
    /// Capillary number, pressure gradient form: Nc = k ||grad p|| / sigma
    double64  Nc_kgradP_Version( double64 pf_gradient_magnitude ) const;
    /// tensor decomposition
    double64 PermeabilityInFlowDirection( const VectorVariable<dim>& mixture_velocity ) const;
  
    // limit approximations for vertical flow arrived at by curve fitting (Maartje 6/9/19)
    double64 krw_VL_LayerPerpendicular() const;
    double64 krn_VL_LayerPerpendicular() const;
    // piecewise defineed functions
    double64 krw_CL_LayerPerpendicular() const;
    double64 krn_CsL_LayerPerpendicular() const;

    // TODO: add function that assesses whether we are dealing with imbibition or drainage
  
    // LAYER-PARALLEL FLOW

    /// Maartje: viscous - capillary force balance (RVC) from capillary number and interfacial tension between wetting and non-wetting phase
    double64 RVC( double64 Ncap ) const;
     /// Maartje: water saturation in the cell at the viscous limit
    double64 Sw_VL( double64 RVC ) const;
    /// Maartje: water saturation in the cell at the capillary limit
    double64 Sw_CL( double64 sw_VL ) const;
  
    /// from the average cell saturation and the viscous-to-capillary force ratio, compute the effective saturations in the laminations at the given capillary number
    void FlowRateDependentLayerSaturations( double64 sw_VL, double64 RVC, double64& sw_low_k_star, double64& sw_high_k_star ) const;
  
    // (BUOYANCY-DRIVEN) FLOW PERPENDICULAR TO THE LAYERS

  private:
    const csmp::Index  RRT_key_, pf_key_, vt_key_;
  
    mutable DenseMatrix<DM_MIN>  DN_; ///< for gradient computations
    TensorVariable<dim>          KK_; ///< tensor permeability
    VectorVariable<dim>          vt_; ///< velocity
    mutable VectorVariable<dim>  vc_; ///< for all kinds of purposes

    // Composite1 - hypothetical sample based on Achyut's rocktypes
    const int      rocktype_ = 1;           ///< FSst-Slt (Fine Sandstone - Silt) Planar Bedding
    // composite is modelled as a dual of two rock types
    const double64 L_low_    = 0.025;       ///< siltstone: cumulative thickness of low-k layers
    const double64 L_high_   = 0.025;       ///< fine sandstone: cumulative thickness of high-k layers
    const double64 k_low_    = 3.4759e-14;  ///< low-permeability lamination (default)
    const double64 k_high_   = 3.6075e-13;  ///< high permeability lamination
    const double64 phi_low_  = 0.19;        ///< porosity of low-perm layer
    const double64 phi_high_ = 0.28;        ///< porosity of high_perm layer
    const double64 m_low_    = 0.5;         ///< van Genuchten exponents for the 2 different layers
    const double64 m_high_   = 0.6;
    const double64 pd_low_   = 3000.;       ///< entry pressures of layers
    const double64 pd_high_  = 1000.;
    const double64 dPc_      = 1.3061e+06;  ///< capillary pressure difference between high_k and low_k layer at connate water saturation
    const double64 Swi_low_  = 0.18;        ///< irreducible saturations of the 2 different layers (CL)
    const double64 Swi_high_ = 0.159;
    const double64 Swc_      = 0.4;         ///< irreducible water saturation of composite measured in lab @todo check

    // water saturation
    mutable double64  sw_;

    // derived quantities
    double64  grad_p_magnitude_ = UNSPECIFIED;
    double64  K_flow_direction_;             ///< permeability in flow direction
    double64  X_PV_low_   =  (phi_low_ * L_low_) / (phi_low_ * L_low_+ phi_high_ * L_high_); ///< pore volume fraction of low perm layer
    double64  Nc_; // capillary mumber

  
    // limits and ratios for composite
    double64 Sw_VL_, RVC_;                     ///< as computed from the correlation between Sw_CL amd Sw_VL (not sure however why that should exist)
    double64 sw_low_k_star_, sw_high_k_star_;  ///< effective saturations in the low and high k layers at the given Nc and flow direction
  
};

} // end namespace csmp

#endif
