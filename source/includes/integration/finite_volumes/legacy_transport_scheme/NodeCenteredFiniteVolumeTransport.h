#ifndef NODE_CENTERED_FINITE_VOLUME_TRANSPORT_H
#define NODE_CENTERED_FINITE_VOLUME_TRANSPORT_H

#include "CSMP_definitions.h"
#include "Box.h"
#include "StencilProcessor.h"
#include "NodeCenteredFiniteVolumeAlgorithm.h"
#include "FiniteVolumeStencilManager.h"
#include "FV_Parameter.h"
#include "GenericNodePropertyGradient.h"
#include "GenericNodePropertyGradientLimiter.h"

// #define DEBUG_NodeCenteredFiniteVolumeTransport

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t>  class Node;
template<size_t>  class Model;
template<size_t>  class TwoPhaseModel;
template<size_t>  class Region;
template<size_t>  class NodeCenteredFiniteVolumeAlgorithm;

/// base class for node-centered generic finite volume calculations on an entire Model object
template<size_t dim>
class NodeCenteredFiniteVolumeTransport {
  public:
    /// single-phase solute advection-ONLY constructor, 1st-order method
    NodeCenteredFiniteVolumeTransport( const char* group_name,
                                       Model<dim>&,
                                       const char* porosity,
                                       const char* advected_prop,
                                       const char* transp_velocity,
                                       const char* nodal_source,
                                       bool second_order_in_space,
                                       bool second_order_in_time,
                                       const char* elmt_thickness_attribute=NULL,
                                       const char* velocity_multiplier=NULL);

    /// single-phase solute advection-diffusion constructor, 1st-order method
    NodeCenteredFiniteVolumeTransport( const char* group_name,
                                       Model<dim>&,
                                       const char* porosity,
                                       const char* diffusivity,
                                       const char* advected_prop,
                                       const char* transp_velocity,
                                       const char* nodal_source,
                                       bool second_order_in_space,
                                       bool second_order_in_time,
                                       const char* elmt_thickness_attribute=NULL,
                                       const char* velocity_multiplier=NULL);

    virtual ~NodeCenteredFiniteVolumeTransport();

    /// single-phase passive advection, returns Courant increment
    virtual double AdvectVariable( double time_interval,
                                     double cfl_multiplication_factor=1.,
                                     bool apply_flux_balance_correction=true, // if there are poro-elastic sources or sinks
                                     bool update_pore_volumes=false);
  
    /// single-phase passive advection, does NOT return courant increment, single timestep calculation
    /// no checks are made for courant condition.  Assumes external checks.
    virtual void AdvectVariableSingleStep( double time_increment,
                                           bool apply_flux_balance_correction=true,
                                           bool update_pore_volumes=false);

    /// advection of one of two phases in two-phase flow handled by subclasses
    virtual double TransportPhase( TwoPhaseModel<dim>&, double time_interval );

    // bijective mapping
    // void MapFiniteVolumeVariableToFiniteElementSpace( const char* variable );

    /// evaluation of Courant-Friedrich-Levy condition, i.e., grid Courant number
    //virtual double CourantIncrement();

    /// single phase, takes into account element shape in 3D
    virtual double AnisotropicCourantIncrement();

    /// advection, capillary diffusion and gravitational flow
    //virtual double CourantIncrement( TwoPhaseModel<dim>& );

    /// two-phase, takes into account element shape in 3D
    virtual double AnisotropicCourantIncrement( TwoPhaseModel<dim>&,
                                                  double max_time_increment=3153600000. ); // 100 years

    /// to set scaling (default=1)
    void     CFL_Multiplier( double desired_value );
    double CFL_Multiplier() const;    

    /// to apply the gradient limiter
    void  WithLsmGradientLimiter(Model<dim>&);
    void  MaxNonlinear2ndOrderIterations(size_t max_nonlinear_iterations);
    void  TargetNonlinear2ndOrderResidual(double target_residual);

    double ModelInflow() const;
    double ModelOutflow() const;

    /// returns inflow, outflow, and mismatch between these; calculation can FVs where variable is flagged with exclude_flag
    double BoundaryFluxes( double& inflow,
                             double& out_flow,
                             bool box_shaped_model                      = true,
                             bool use_advected_variable                 = false,
                             const char* advected_variable              = "no advected variable",
                             const char* variable_to_determine_no_flow  = "fluid pressure",
                             VARIABLE_FLAG exclude_flag                 = ANY
                           ) const;
  
    /// for each finite volume it sums up incoming and outgoing fluxes reporting the extrema
    void FluxBalance( double& fmin, double& fmax ) const;

    /// as previous method, but assigns surface integrals to the result property
    void Divergence( const char* div_property, const char* result_prop ) const;

    /// integrates nodal scalar properties finite volume by finite volume
    void MultiplyScalarNodePropertyByFiniteVolume( const char* property );

    /// integrates scalar node and element variables
    double VolumeIntegrateScalarFiniteVolumeVariable( const char* property, bool consider_porosity ) const;
    double VolumeIntegrateScalarFiniteElementVariable( const char* property, bool take_porosity_into_account = false, const char *region=NULL ) const;

    /// integrates target property over a region of interest, result can be scaled by porosity
    double VolumeIntegrateScalarFiniteVolumeVariable( const char* group, const Model<dim>&,
                                                        const char* property, bool consider_porosity ) const;
    /// integrals for finite-element computations
    void     VolumeIntegrate( const char* integrand_property, const char* result_prop ) const;
    void     VolumeIntegrate( const char* integrand_property, const char* integrand_multiplier, const char* result_prop ) const;

    /// conversion of values already assigned to the nodes (only possible if nodes comprise model boundary)
    void MultiplyScalarBoundaryValuesByFiniteVolumeCrossSectionalArea( Model<dim>&,BOX_BOUNDARY, const char* property );

    /// computation and assignment of Neumann boundary gradients to target boundary
    void TransformScalarBoundaryValuesIntoNeumannConditions( BOX_BOUNDARY boundary,
                                                             const char* node_property,
                                                             const char* propertionality_constant,
                                                             bool distribute_total_amount );

    /// assignment of a single specific value that is projected on outward point normal
    void AssignScalarBoundaryValues( BOX_BOUNDARY, const char* property, VARIABLE_FLAG bcond,
                                     double val, bool distribute_total_amount );

    /// returns the actual volume of the finite volume i
    double FiniteVolume( const char* volume_property ) const;
    double PoreVolume( size_t iFv_cell ) const;

    /// Adjust flexible samg solver settings
    virtual void AdjustSolverSettings();

    /// Access to advector's solver settings
    virtual CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS& GetSolverSettings();
    virtual CSMP_DEFAULT_LINEAR_SOLVER* GetSolver();

    virtual void  Out() const;

    void Verbose(bool verbose){this->verbose_=verbose; if (this->baseAdvector_) this->baseAdvector_->Verbose(verbose);}
    bool Verbose(){return this->verbose_;}

  protected:

    NodeCenteredFiniteVolumeTransport();
    NodeCenteredFiniteVolumeTransport( const NodeCenteredFiniteVolumeTransport<dim>& );
    void   InitializeFiniteVolumeData( bool multiply_pore_volumes_with_thickness_attribute=false );
    void   InitializeFiniteVolumeDataParametricToPhysical( bool store_normals );
    void   InitializeSectorPoreVolumeData( bool multiply_pore_volumes_with_thickness_attribute );
    bool   InitializeArraysForFirstOrderMethod();
    bool   InitializeArraysForSecondOrderMethod( bool second_order_in_space, bool second_order_in_time );
    virtual void   CheckTransportVariables() const;
    virtual size_t MeasureAllocatedMemory() const;

    bool SecondOrderInTime() const;
    bool SecondOrderInSpace() const;

    void InitialAdvectedPropertyValues( const csmp::Index& adv_key );

    void InitialAdvectedPropertyValues( const csmp::Index& adv_key,
                                        double& smin, double& smax );

    void MinMaxAdvectedProperty();

    void BackupFacetFluxes();

    void UpdateProjectedVelocitiesAndFluxBalances();

    bool FluxThroughBoundaryFiniteVolume( const Node<dim>* fv_ptr, double& inflow, double& flux_balance ) const;

    void AssignFluxBoundaryConditions( NodeCenteredFiniteVolumeAlgorithm<dim>& advection_algorithm,
                                       bool use_t0_saturation ) const;

    // outputs the flux balance in each FV cell to the variable 'balance var' and integrates it
    void RecordFluxBalances( const char* balance_var,
                             double& total_surplus, double& total_deficit ) const;

    void AdvectVariable1stOrder( NodeCenteredFiniteVolumeAlgorithm<dim>&,
                                 double time_interval, bool with_flux_balance_correction );
    // 2nd-order in space
    void AdvectVariable2ndOrder( NodeCenteredFiniteVolumeAlgorithm<dim>&,
                                 double time_increment, bool with_flux_balance_correction );

    void AdvectVariable2ndOrderInSpaceAndTime( NodeCenteredFiniteVolumeAlgorithm<dim>&,
                                                 double time_increment, bool with_flux_balance_correction );

protected:
    const PropertyDatabase<dim>& pref_;
    Region<dim>&                 gref_;
    Model<dim>&                  mref_;

    bool                         verbose_;

    std::vector<FV_Parameter>           STENCIL_DATA;  ///< finite-volume data for fast calculations
    std::vector<double>               FVPOREVOL,     ///< pore volumes of finite-volume cells for theta calculation etc.
                                        SAT0,          ///< advected variable at t0
                                        FLUX_BALANCE;  ///< for each FV, but excluding truncated FV at global model boundaries

    // arrays for second-order accurate calculations
    std::vector<std::pair<double,double> >  SMINMAX; ///< min, max values of advected property at neighbor FV's of current FV    
    std::vector<std::vector<double> > FACETFLUXES0,  ///< at time, t (at last time step)
                                        LTDSATS0;      ///< ltd values of advected property

    const csmp::Index                   phi_key_;       ///< porosity
    const csmp::Index                   diff_key_;      ///< diffusivity of advected variable (only in 2nd-order method)
    const csmp::Index                   adv1_key_;      ///< transported variable
    const csmp::Index                   vel_key_;       ///< transport velocity
    const csmp::Index                   src_key_;       ///< nodal fluid volume source
    const csmp::Index                   thi_key_;       ///< thickness multiplier for lower dimensional elements
    const csmp::Index                   velo_mult_key_; ///< velocity multiplier key
    csmp::Index                         mass_center_key_;			///< center of mass for each CV
    csmp::Index                         grad_advprop_key_;			///< gradient of the advacted property
    csmp::Index                         grad_advprop_limiter_key_;	///< limiter for the gradient of the advected property
    std::string                         advected_variable_;
    double                            cfl_multiplier_;
    size_t                              var_ncomponents_; /// default size is 1 if only a scalar is being advected. (Julian)

    StencilProcessor<dim>               stencil_;
    bool                                firstCall_;
    bool                                with_lsmgrad_limiter_;

    double target_nonlinear_limiting_case_residual_; ///< targed residual in order to get correct 2nd order approximations
    size_t max_nonlinear_limiting_case_iterations_;    ///< maximum number of nonlinear iterations in 2nd order scheme

    NodeCenteredFiniteVolumeAlgorithm<dim>*  baseAdvector_ = nullptr;
    GenericNodePropertyGradientLimiter<dim>*  grad_advprop_limiter_ = nullptr;
};



template<size_t dim>
inline double NodeCenteredFiniteVolumeTransport<dim>::PoreVolume( size_t fv_cell ) const
 {
    assert( fv_cell < FVPOREVOL.size() );
    return FVPOREVOL[fv_cell];
 }


// verify underlying generic FV scheme
template<size_t dim>
bool testFiniteVolumeStencil( const PropertyDatabase<dim>&, const Region<dim>&,
                              NodeCenteredFiniteVolumeTransport<dim>&, long data_precision=3 );


/**
@class NodeCenteredFiniteVolumeTransport  NodeCenteredFiniteVolumeTransport "generic_node_centered_finite_volumes\NodeCenteredFiniteVolumeTransport.h"
@author S.K. Matthaei
@date 2008


Generic transport module for CVFE computations on poly-element-type finite
element (-finite-volume) meshes. Computations are done using an implicit
or semi-implicit approach that is only weakly depedent on the
grid-Courant number.

A range of constructors allows to build the module for first-order,
second-order conservative and non-conservative transport simulations
including source and sink terms.

@todo (1) Bring theta-limited space time second-order accurate scheme back to work (A)
*/

} // end namespace csmp



#endif
