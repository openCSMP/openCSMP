#include "TwoPhaseImplicitNodeCenteredFVTransport.h"
#include "NodeCenteredFiniteVolumeAlgorithm.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "PropertyDatabase.h"
#include "TwoPhaseModel.h"
#include "Exception.h"
#include "StencilProcessor.h"

using namespace std;

namespace csmp {

/**

Constructor for the modelling of two-phase slightly compressible flow 
with fluid density and viscosity variations. The properties 
fluid density and viscosity are expected to be node variables and
are interpolated to the FV facet integration points to carry out the
computations. The transport scheme is setup for the non-wetting 
phase. The onus is on the user to update the second phase saturations
after each advection step, i.e., calculate sw=1-sn.  

The TwoPhaseImplicitNodeCenteredFVTransport can be constructed for first- and second
order calculations.  

Diffusion is only taken into account if second-order is chosen, because
the first-order scheme already is so diffusive that adding an extra
diffusion term makes little sense.   

@section arguments Input Arguments 

The TwoPhaseImplicitNodeCenteredFVTransport gains access via reference to the 
model to which it is applied to. The user needs to specify the piecewise
constant element variables porosity, diffusivity and transport velocity.
These must be defined in the CSP_variables.txt file. The variables 
fluid density and viscosity of the wetting and non-wetting phases must
be nodal variables as well. They are interpolated to the facet and sector
integration points using the finite element basis functions.  

The transported variable must be nodal since this explicit transport scheme 
is based on node-centered finite volumes.  

@section implementation Implementation

This constructor creates a minor storage array for repetitively used
variables like the finite (pore) volumes, the flux balance for each
FV cell and a vector<double64> of FV_parameters called FV_stencil_data that holds
the FV sector volumes, facet areas and facet normal fluxes which are 
always updated before the first advection step, when the CFL criterion
is computed. This vector<double64> will become redundant when the generic FV scheme
has been optimized for speed.  

@section application Application

Use TwoPhaseImplicitNodeCenteredFVTransport as constructed with this particular
constructor for two-phase flow in a rapidly evolving
flow field that necessitates frequent updates of the velocity fields.  

For a weak coupling between phase saturations and total mobility, use 
the implicit transport module FiniteVolumeTransport.  

@section messages Messages 

The memory that is consumed during the construction process is reported.

Constructor for first-order scheme.

*/
template<size_t dim, template<size_t> class STP>
TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::TwoPhaseImplicitNodeCenteredFVTransport( // for entire model
                                                                    const char* group_name,
                                                                    Model<dim>& sg,
                                                                    const char* porosity,
                                                                    const char* viscosity_n,
                                                                    const char* viscosity_w,
                                                                    const char* density_n,
                                                                    const char* density_w,
                                                                    const char* phase1_to_update,
                                                                    const char* phase2_to_advect,
                                                                    const char* transp_velocity,
                                                                    const char* nodal_source,
                                                                    bool with_capillary_spreading,
                                                                    bool with_gravitational_forces,
                                                                    const char* reference_variable_to_reference_variable_to_no_flow_bc,
                                                                    bool nonlinear_scheme )
 : NodeCenteredFiniteVolumeTransport<dim>( group_name, sg,
                                           porosity,
                                           phase2_to_advect, transp_velocity,
                                           nodal_source, false, false ),
   // general flow-related variables
   rhn_key(this->pref_.StorageKey(density_n)),         // density oil
   rhw_key(this->pref_.StorageKey(density_w)),         // density water
   mun_key(this->pref_.StorageKey(viscosity_n)),       // viscosity oil
   muw_key(this->pref_.StorageKey(viscosity_w)),       // viscosity water
   ph1_key(this->pref_.StorageKey(phase1_to_update)),  // the first phase that is not advected but updated
   reference_variable_to_no_flow_bc_key_(this->pref_.StorageKey(reference_variable_to_reference_variable_to_no_flow_bc)),
   with_capillary_spreading_(with_capillary_spreading),
   with_gravitational_forces_(with_gravitational_forces),
   advector_(this->gref_),
   nonlinear_scheme_(nonlinear_scheme),
   target_newton_raphson_residual_(1.0e-3),
   max_newton_raphson_iterations_(5U),
   max_line_search_iterations_(5U)
 {
    // testing whether all variables are placed and initialized as expected
    CheckTransportVariables();

    // establishing the halo stencils
    for ( typename vector<Node<dim>*>::const_iterator
          nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
      for ( size_t i=0U; i<(*nit)->Parents(); i++ )
        if ( !IsInteriorStencil( (*nit)->Parent(i) ) )
          halo_stencils_.insert( (*nit)->Parent(i) );

 } // end constructor (multiphase fluid flow, no diffusion)





template<size_t dim, template<size_t> class STP>
TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::TwoPhaseImplicitNodeCenteredFVTransport( 
                                                                    const char* group_name,
                                                                    Model<dim>& sg,
                                                                    bool second_order_in_space,
                                                                    bool second_order_in_time,
                                                                    const char* porosity,
                                                                    const char* viscosity_n,
                                                                    const char* viscosity_w,
                                                                    const char* density_n,
                                                                    const char* density_w,
                                                                    const char* phase1_to_update,
                                                                    const char* phase2_to_advect,
                                                                    const char* transp_velocity,
                                                                    const char* nodal_source,
                                                                    bool with_capillary_spreading,
                                                                    bool with_gravitational_forces,
                                                                    const char* reference_variable_to_no_flow_bc,
                                                                    bool nonlinear_scheme )
 : NodeCenteredFiniteVolumeTransport<dim>( group_name, sg,
                                           porosity,
                                           phase2_to_advect, transp_velocity,
                                           nodal_source, second_order_in_space, second_order_in_time),
   // general flow-related variables
   rhn_key(this->pref_.StorageKey(density_n)),         // density oil
   rhw_key(this->pref_.StorageKey(density_w)),         // density water
   mun_key(this->pref_.StorageKey(viscosity_n)),       // viscosity oil
   muw_key(this->pref_.StorageKey(viscosity_w)),       // viscosity water
   ph1_key(this->pref_.StorageKey(phase1_to_update)),  // the first phase that is not advected but updated
   reference_variable_to_no_flow_bc_key_(this->pref_.StorageKey(reference_variable_to_no_flow_bc)),
   with_capillary_spreading_(with_capillary_spreading),
   with_gravitational_forces_(with_gravitational_forces),
   advector_(this->gref_),
   nonlinear_scheme_(nonlinear_scheme),
   target_newton_raphson_residual_(1.0e-3),
   max_newton_raphson_iterations_(5U),
   max_line_search_iterations_(5U)
 {
    // testing whether all variables are placed and initialized as expected
    CheckTransportVariables();

    // establishing the halo stencils
    for ( typename vector<Node<dim>*>::const_iterator
          nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
      for ( size_t i=0U; i<(*nit)->Parents(); i++ )
        if ( !IsInteriorStencil( (*nit)->Parent(i) ) )
          halo_stencils_.insert( (*nit)->Parent(i) );

 } // end constructor (multiphase fluid flow, no diffusion)



/**

Constructor for the modelling of two-phase slightly compressible flow
with fluid density and viscosity variations. The properties
fluid density and viscosity are expected to be node variables and
are interpolated to the FV facet integration points to carry out the
computations. The transport scheme is setup for the non-wetting
phase. The onus is on the user to update the second phase saturations
after each advection step, i.e., calculate sw=1-sn.

The TwoPhaseImplicitNodeCenteredFVTransport can be constructed for first- and second
order calculations.

Diffusion is only taken into account if second-order is chosen, because
the first-order scheme already is so diffusive that adding an extra
diffusion term makes little sense.

@section arguments Input Arguments

The TwoPhaseImplicitNodeCenteredFVTransport gains access via reference to the
model to which it is applied to. The user needs to specify the piecewise
constant element variables porosity, diffusivity and transport velocity.
These must be defined in the CSP_variables.txt file. The variables
fluid density and viscosity of the wetting and non-wetting phases must
be nodal variables as well. They are interpolated to the facet and sector
integration points using the finite element basis functions.

The transported variable must be nodal since this explicit transport scheme
is based on node-centered finite volumes.

@section implementation Implementation

This constructor creates a minor storage array for repetitively used
variables like the finite (pore) volumes, the flux balance for each
FV cell and a vector<double64> of FV_parameters called FV_stencil_data that holds
the FV sector volumes, facet areas and facet normal fluxes which are
always updated before the first advection step, when the CFL criterion
is computed. This vector<double64> will become redundant when the generic FV scheme
has been optimized for speed.

@section application Application

Use TwoPhaseImplicitNodeCenteredFVTransport as constructed with this particular
constructor for two-phase flow in a rapidly evolving
flow field that necessitates frequent updates of the velocity fields.

For a weak coupling between phase saturations and total mobility, use
the implicit transport module FiniteVolumeTransport.

@section messages Messages

The memory that is consumed during the construction process is reported.
*/
template<size_t dim, template<size_t> class STP>
TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::TwoPhaseImplicitNodeCenteredFVTransport( // for entire model
                                                              const char* group_name,
                                                              Model<dim>& sg,
                                                              const char* porosity,
                                                              const char* diffusivity,
                                                              const char* viscosity_n,
                                                              const char* viscosity_w,
                                                              const char* density_n,
                                                              const char* density_w,
                                                              const char* phase1_to_update,
                                                              const char* phase2_to_advect,
                                                              const char* transp_velocity,
                                                              const char* nodal_source,
                                                              bool with_capillary_spreading,
                                                              bool with_gravitational_forces,
                                                              const char* reference_variable_to_no_flow_bc,
                                                              bool nonlinear_scheme)
 : NodeCenteredFiniteVolumeTransport<dim>( group_name, sg,
                                           porosity, diffusivity,
                                           phase2_to_advect, transp_velocity,
                                           nodal_source, false, false ),
   // general flow-related variables
   rhn_key(this->pref_.StorageKey(density_n)),         // density oil
   rhw_key(this->pref_.StorageKey(density_w)),         // density water
   mun_key(this->pref_.StorageKey(viscosity_n)),       // viscosity oil
   muw_key(this->pref_.StorageKey(viscosity_w)),       // viscosity water
   ph1_key(this->pref_.StorageKey(phase1_to_update)),  // the first phase that is not advected but updated
   reference_variable_to_no_flow_bc_key_(this->pref_.StorageKey(reference_variable_to_no_flow_bc)),
   with_capillary_spreading_(with_capillary_spreading),
   with_gravitational_forces_(with_gravitational_forces),
   advector_(this->gref_),
   nonlinear_scheme_(nonlinear_scheme),
   target_newton_raphson_residual_(1.0e-3),
   max_newton_raphson_iterations_(5U),
   max_line_search_iterations_(5U)
 {
    // testing whether all variables are placed and initialized as expected
    CheckTransportVariables();

    // establishing the halo stencils
    for ( typename vector<Node<dim>*>::const_iterator
          nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
      for ( size_t i=0U; i<(*nit)->Parents(); i++ )
        if ( !IsInteriorStencil( (*nit)->Parent(i) ) )
          halo_stencils_.insert( (*nit)->Parent(i) );

 } // end constructor (multiphase fluid flow)

/** Second-order constructor, with capillary spreading
*/

template<size_t dim, template<size_t> class STP>
TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::TwoPhaseImplicitNodeCenteredFVTransport( // for entire model
                                                              const char* group_name,
                                                              Model<dim>& sg,
                                                              bool second_order_in_space,
                                                              bool second_order_in_time,
                                                              const char* porosity,
                                                              const char* diffusivity,
                                                              const char* viscosity_n,
                                                              const char* viscosity_w,
                                                              const char* density_n,
                                                              const char* density_w,
                                                              const char* phase1_to_update,
                                                              const char* phase2_to_advect,
                                                              const char* transp_velocity,
                                                              const char* nodal_source,
                                                              bool with_capillary_spreading,
                                                              bool with_gravitational_forces,
                                                              const char* reference_variable_to_no_flow_bc,
                                                              bool nonlinear_scheme)

 : NodeCenteredFiniteVolumeTransport<dim>( group_name, sg,
                                           porosity, diffusivity,
                                           phase2_to_advect, transp_velocity,
                                           nodal_source, second_order_in_space, second_order_in_time ),
   // general flow-related variables
   rhn_key(this->pref_.StorageKey(density_n)),         // density oil
   rhw_key(this->pref_.StorageKey(density_w)),         // density water
   mun_key(this->pref_.StorageKey(viscosity_n)),       // viscosity oil
   muw_key(this->pref_.StorageKey(viscosity_w)),       // viscosity water
   ph1_key(this->pref_.StorageKey(phase1_to_update)),  // the first phase that is not advected but updated
   reference_variable_to_no_flow_bc_key_(this->pref_.StorageKey(reference_variable_to_no_flow_bc)),
   with_capillary_spreading_(with_capillary_spreading),
   with_gravitational_forces_(with_gravitational_forces),
   advector_(this->gref_),
   nonlinear_scheme_(nonlinear_scheme),
   target_newton_raphson_residual_(1.0e-3),
   max_newton_raphson_iterations_(5U),
   max_line_search_iterations_(5U)
 {
    // testing whether all variables are placed and initialized as expected
    CheckTransportVariables();

    // establishing the halo stencils
    for ( typename vector<Node<dim>*>::const_iterator
          nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
      for ( size_t i=0U; i<(*nit)->Parents(); i++ )
        if ( !IsInteriorStencil( (*nit)->Parent(i) ) )
          halo_stencils_.insert( (*nit)->Parent(i) );

 } // end constructor (multiphase fluid flow)



/**
     Constructor for field-scale simulation:
    
     with standard variables, gravity and thickness attribute for lower-dim elements, no-pc
*/
template<size_t dim, template<size_t> class STP>
TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::TwoPhaseImplicitNodeCenteredFVTransport(
                                                                     const char* group_name,
                                                                     Model<dim>& sg,
                                                                     bool second_order_in_space,
                                                                     bool second_order_in_time,
                                                                     const char* transp_velocity,
                                                                     bool with_capillary_forces,
                                                                     const char* thickness_low_dim_elmt,
                                                                     bool nonlinear_scheme )
 : NodeCenteredFiniteVolumeTransport<dim>( group_name, sg,
                                          "porosity",
                                          "saturation carbonic phase", transp_velocity,
                                          "nodal fluid volume source",
                                           second_order_in_space, second_order_in_time,
                                           thickness_low_dim_elmt ),
   // general flow-related variables
   rhn_key(this->pref_.StorageKey("density carbonic phase")),         
   rhw_key(this->pref_.StorageKey("density aqueous phase")),
   mun_key(this->pref_.StorageKey("viscosity carbonic phase")),
   muw_key(this->pref_.StorageKey("viscosity aqueous phase")),
   ph1_key(this->pref_.StorageKey("saturation aqueous phase")),  // the first phase that is not advected but updated
   reference_variable_to_no_flow_bc_key_(this->pref_.StorageKey("fluid pressure")),
   with_capillary_spreading_(with_capillary_forces),
   with_gravitational_forces_(true),
   advector_(this->gref_),
   nonlinear_scheme_(nonlinear_scheme),
   target_newton_raphson_residual_(1.0e-3),
   max_newton_raphson_iterations_(5U),
   max_line_search_iterations_(5U)
 {
    // testing whether all variables are placed and initialized as expected
    CheckTransportVariables();

    // establishing the halo stencils
    for ( typename vector<Node<dim>*>::const_iterator
          nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
      for ( size_t i=0U; i<(*nit)->Parents(); i++ )
        if ( !IsInteriorStencil( (*nit)->Parent(i) ) )
          halo_stencils_.insert( (*nit)->Parent(i) );

 } // end constructor (multiphase fluid flow, no diffusion or capillary transfer)






template<size_t dim, template<size_t> class STP>
inline bool TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::IsInteriorStencil(
                                                      const Element<dim>* const eptr ) const
 {
    return this->gref_.IsPerimeterElement( eptr );
 }



/** Frees all the allocated memory.

  */
template<size_t dim, template<size_t> class STP>
TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::~TwoPhaseImplicitNodeCenteredFVTransport()
 {
 } 




template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::DisableCapillarySpreading()
  {
     with_capillary_spreading_ = false;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::EnableCapillarySpreading()
  {
     with_capillary_spreading_ = true;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::DisableGravitationalForces()
  {
     with_gravitational_forces_ = false;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::EnableGravitationalForces()
  {
     with_gravitational_forces_ = true;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SetNoFlowBoundaryConditionKey(Model<dim>& sg,const char* reference_variable_to_no_flow_bc)
  {
     reference_variable_to_no_flow_bc_key_ = sg.Database().StorageKey(reference_variable_to_no_flow_bc);
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::ApplyLinearScheme()
  {
     nonlinear_scheme_ = false;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::ApplyNonlinearScheme()
  {
    nonlinear_scheme_ = true;
  }



template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::MaxNewtonRahsonIterations(size_t max_newton_raphson_iteratons)
{
    max_newton_raphson_iterations_ = max_newton_raphson_iteratons;
}


template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::MaxLineSearchIterations(size_t max_line_search_iterations)
{
    max_line_search_iterations_ = max_line_search_iterations;
}

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::TargetNewtonRaphsonResidual(double64 target_residual)
{
    target_newton_raphson_residual_ = target_residual;
}



/**
 
Computes the Courant time increment (CFL citerion) taking into account
viscous, gravitational and capillary fluid displacements using the 
contraints from the provided relative permeability model. The CFL
criterion is calculated using the element diameter in the direction
of the flow but not account for the deviation from this vector<double64> of 
the flow of the considered phase.  

For the viscous flow the shock speed is used as a multiplier for the 
total velocity. This may lead to a too tight constraint but is safe
for the case where CFL is computed only once at the onset of series
of advection steps with an explicit scheme.   

@section arguments Input Arguments 

First, a const reference to the current model, a reference to the 
TwoPhaseModel used (this will be a subclass), and, third, the maximum
time-increment imposed by external constraints.  

@return  The CFL criterion.

@section implementation Implementation

The criterion is computed sector by sector for all the finite
elements in the entire model.  

@section application Application

Method has been designed primarily to give stability to an explicit
transport scheme but it also serves as a guide for the timestepping
using an implicit approach.

@section messages Messages 

The method alerts the user to the case where CFL is larger than the
input max_time_increment which may be the case if there is no flow at all
in the domain. Equally, the user is informed if the CFL increment is 
less than a millisecond (usually a prohibitively small increment).  
*/
template<size_t dim, template<size_t> class STP>
double64 TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::AnisotropicCourantIncrement( TwoPhaseModel<dim>& relperm,
                                                                                        double64 max_time_increment )
 {
    //assert( dim != 1U );
    assert( max_time_increment > 0. );
    this->UpdateProjectedVelocitiesAndFluxBalances();
    
    static DenseMatrix<DM_MIN>  DN;
    vector<double64>            gradPc(dim);
    VectorVariable<dim>         vc;
    double64                    velocity,
                                courant_increment(max_time_increment);
    const double64              millisecond(1.0e-3);
    const bool                  multiply_with_cell_thickess = (this->thi_key_ == csmp::Index()) ? false : true;
    const bool                  unless_has_equal_dimension(dim!=1U);

    for ( typename vector<Element<dim>*>::const_iterator 
          eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++ )
      if ( !((*eit)->FE()->IsLineElement() && unless_has_equal_dimension) )
        {
           // 0. relative permeability model computed at element barycenter
           // -------------------------------------------------------------
           relperm.Initialize( *(*eit) );
           relperm.InitializeForBaryCenter( *(*eit) );
           relperm.EffectiveSaturation();

           // 1. limit imposed by advection
           // -----------------------------
           (*eit)->Read( this->vel_key_, vc );
           velocity = vc.Length();

           // NB: This may be a too conservative estimate for the implicit scheme but is necessary for the explicit one
           // velocity = std::max( velocity, velocity * relperm.MaxFractionalFlowDerivative() );
           // thus for implicit scheme is better:
           velocity *= relperm.dfds();

           double64 cell_diameter = (*eit)->LengthInDirection(vc) * (*eit)->Read( this->phi_key_ );
           if ( multiply_with_cell_thickess ) cell_diameter *= (*eit)->Read( this->thi_key_ );

           // 2. additional buoyancy-related flow
           // -------------------------------------
           if(with_gravitational_forces_)
               velocity += fabs(relperm.GravityMultiplier_dGds());

           // 3. limit due to potential capillary spreading
           // ---------------------------------------------
           if(with_capillary_spreading_){

               const double64 k_lambda_overbar(relperm.Permeability() * relperm.G());
               // computing the capillary pressure gradient
               if ( k_lambda_overbar > numeric_limits<double64>::epsilon() )
               {
                   fill( gradPc.begin(), gradPc.end(), 0. );
                   (*(*eit)).dN_AtBaryCenter( DN );
                   for ( size_t j=0U; j<(*eit)->Nodes(); j++ ) {
                       double64 sn = (*eit)->N(j)->Read( this->adv1_key_ );
                       relperm.SaturationWettingPhase( 1. - sn );
                       relperm.EffectiveSaturation();
                       double64 pc = relperm.pc_Phase( );
                       for ( size_t k=0U; k<dim; k++ ) gradPc[k] += DN(k,j) * pc;
                   }
                   // getting the maximum capillary flux (G= lambda overbar)
                   double64  magnitude_grad_pc(gradPc[0]); // 1D
                   if      ( dim == 3U ) magnitude_grad_pc = sqrt(gradPc[0]*gradPc[0]+gradPc[1]*gradPc[1]+gradPc[2]*gradPc[2]);
                   else if ( dim == 2U ) magnitude_grad_pc = hypot(gradPc[0],gradPc[1]);

                   // use data from barycenter: O.K. as long as grad_pc does not increase during iterations
                   velocity += fabs( magnitude_grad_pc * k_lambda_overbar );
               }
           }

           // 4. calculating the CFL criterion from the cell diameter
           // -------------------------------------------------------
           // guarding against degenerate cases
           if ( velocity > numeric_limits<double64>::epsilon() and cell_diameter > numeric_limits<double64>::epsilon() ) {
                courant_increment = std::min( courant_increment, cell_diameter / velocity );
             }
        }

    if ( courant_increment >= max_time_increment ) { 
         cout <<"\nINFO, TwoPhaseImplicitNodeCenteredFVTransport::AnisotropicCourantIncrement (2-phase flow):";
         cout <<" calculated courant increment is larger than maximum permitted increment, there may be no flow in the model domain !";
         cout <<"\n\tCFL is set to "<< fixed << setprecision(0) << max_time_increment << endl;
         return max_time_increment;
      }
    if ( courant_increment <= millisecond ) {
         throw csmp::Exception( WARNING, "TwoPhaseImplicitNodeCenteredFVTransport::AnisotropicCourantIncrement (2-phase flow):",
                                         "courant increment is smaller than a millisecond. Check your boundary conditions." );
      }
    else {
         cout <<"\nTwoPhaseImplicitNodeCenteredFVTransport::AnisotropicCourantIncrement (2-phase flow): ";
         cout << fixed << setprecision(0) << courant_increment <<" and CFL multiplier: "<< this->cfl_multiplier_ << endl;
      }
      
    return courant_increment;

 } // end AnisotropicCourantIncrement (2-phase flow - all)







/**
 
Checks whether the key variables used in the transport calculations have 
the correct placement, i.e. node or element. If not, FATAL_ERRORs are 
raised to stop the user from proceeding.  

@section arguments Input Arguments 

The method checks either the basic variables that the user specified 
to denote the element variables porosity, diffusivity, transport velocity
and transported variable, or it also checks the whether fluid densities
and viscosities for the two phases have been specified correctly
(two_phase_flow=true).   

@section messages Messages 

Falsely specified variables are reported and need to be fixed in the 
specific CSP_variables.txt file used for the simulation.  
 */
template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::CheckTransportVariables() const
  {
     NodeCenteredFiniteVolumeTransport<dim>::CheckTransportVariables();
   
     // if two-phase flow is enabled 
     if ( mun_key.index != muw_key.index )
	     // fluid viscosities and densities
	     if ( mun_key.place != NODE   || muw_key.place != NODE ||
	          rhn_key.place != NODE   || rhw_key.place != NODE || 
	          mun_key.type  != SCALAR || muw_key.type  != SCALAR ||
	          rhn_key.type  != SCALAR || rhw_key.type  != SCALAR )
	       throw csmp::Exception( FATAL_ERROR, "TwoPhaseImplicitNodeCenteredFVTransport::CheckTransportVariables", 
	                                           "Fluid viscosities and densities must be scalar node properties" );

 } // end CheckTransportVariables

/**

The method is the equivalent to AdvectVariable() for 2-phase flow. Thus,
it advects the transport variable for the user-specified time interval
and writes the results back to the current Model object.  

@section arguments Input Arguments 

The method takes a reference to the Model object which contains 
the transport variable field that shall be advected. The relative 
permeabilities and other 2-phase properties are obtained via the 
supplied 'AdvectionDiffusionTwoPhaseModel'. The onus is on the user
to make sure that the necessary variables are defined and initialized.
The third method argument is the time interval over which the 
transport variable shall be advected.  

@return The advected variable field is written back into the supplied Model.
The method returns the CFL condition that it determined before the 
onset of the advection.
 
@section implementation Implementation

The CFL_FACTOR is a multiplier that represents the maximum factor by
which the diameter of the idealized spherical FV can be multiplied
while the transport scheme is stable.  
*/

template<size_t dim, template<size_t> class STP>
double64  TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::TransportPhase( TwoPhaseModel<dim>& relperm,
                                                                            double64 time_interval )
 {
    this->gref_.UpdateMemberIndexes();

    // 1. compute the CFL condition to identify value for overstepping
    const double64 courant_increment = AnisotropicCourantIncrement( relperm, time_interval );

    // 2. applying CFL multiplication factor for overstepping
    double64 time(0.), time_increment = courant_increment * this->cfl_multiplier_;

    cout <<"\nTwoPhaseImplicitNodeCenteredFVTransport<"<< fixed << setprecision(0) << dim <<",STP>::TransportPhase:";
    cout<<"\n\tTime interval     = "<< time_interval;
    cout<<"\n\tCourant increment = "<< courant_increment;
    cout<<"\n\tCFL multiplier    = "<< this->cfl_multiplier_;
    cout<<"\n\tTime increment    = "<< time_increment << endl;

    cout <<"\n\n\nTwoPhaseImplicitNodeCenteredFVTransport::TransportPhase: Advection initiated...\n";
    size_t  substep(1);
    while ( time < time_interval ) {
         cout <<"\n\n\tadvection (sub)step: "<< substep;
         if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
         if ( nonlinear_scheme_ )
           {
              if ( this->SecondOrderInSpace() )
                  SolveTransportEquation2ndOrderInSpace_NonlinearNewtonRaphson( this->mref_, relperm, time_increment );
              else
                  SolveTransportEquation1stOrder_NonlinearNewtonRaphson( relperm, time_increment );

           } else SolveTransportEquation1stOrder( relperm, time_increment );

         time += time_increment;
         substep++;
      }
    cout <<"\nTwoPhaseImplicitNodeCenteredFVTransport::TransportPhase: 'Advection completed.\n";
    return courant_increment;

 } // end TransportPhase

#ifdef CSMP_WITH_SAMG_SOLVER
template<size_t dim, template<size_t> class STP>
SAMG_Settings& TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::GetSolverSettings()
 {
   return advector_.GetSolverSettings();
 }
#else
/// add extra functionality for laternative solver
#endif

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::AdjustSolverSettings()
 {
#if CSMP_WITH_SAMG_SOLVER
     GetSolverSettings().Set_iout1( 0 );
     GetSolverSettings().Set_iout2( 0 );

     /// SAMG solution criterion
     GetSolverSettings().Set_eps(0.);
     GetSolverSettings().Set_rel_eps(1.E-10);

     /// SAMG output to file
     //GetSolverSettings().Set_idmp( 8 );        // define SAMG command and file output
     //GetSolverSettings().Set_ioform( "f" );    // define SAMG file output format for reduced file size, idmp > 1 is required
     //GetSolverSettings().Set_filnam_dump( "SAMG_Transport" ); // set filename for SAMG file output other than default "level", idmp > 1 is required
#else
/// add extra functionality for laternative solver
#endif

 } // end AdjustSolverSettings


// Not needed when sector balances are obtained automatically and are compensated for:

// 4. correct RHS for inflow contributions only at model boundaries
//    advector_.SubtractNonConservativeOutFluxesFromRHS( sg, relperm, this->FLUX_BALANCE, this->ad1_key );






template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::MinMaxAdvectedProperty()
{
   typename vector<pair<double64,double64> >::iterator  sit(this->SMINMAX.begin());

   for ( typename vector<Node<dim>*>::const_iterator
         nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++, sit++ ) {
         // 1. the advected property value at the current node is assigned to min-max pair
        (*sit).first = (*sit).second = (*nit)->Read( this->adv1_key_ );
        for ( size_t i=0U; i<(*nit)->Neighbors(); i++ ) {
             const double64 adv_var((*nit)->Neighbor(i)->Read( this->adv1_key_ ));
             // if element value is smaller the current minimum is assigned etc.
             (*sit).first  = std::min( (*sit).first,  adv_var );
             (*sit).second = std::max( (*sit).second, adv_var );
          }
     } // end for all nodes

 } // end MinMaxAdvectedProperty


template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::MinMaxAdvectedPropertyIncludingTheCurrentNode()
{
   typename vector<pair<double64,double64> >::iterator  sit(this->SMINMAX.begin());
   Element<dim>     current_el;
   size_t  current_n_id,global_neighb_el_id;
   for ( typename vector<Node<dim>*>::const_iterator
         nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++, sit++ ) {

         current_n_id= (*(*nit)).Idx();
         // 1. the advected property value at the current node is assigned to min-max pair
         (*sit).first = (*sit).second = (*nit)->Read( this->adv1_key_ );

         for(  size_t p = 0; p< (*nit)->Parents() ; p++ ){
               // get the global parent id:
               global_neighb_el_id = (*nit)->Parent( p)->Idx();
               // get the corresponding element:
               current_el = *this->gref_.E( global_neighb_el_id );
               for(size_t i=0;i<current_el.Nodes();i++){
                   //ids[i]=current_el.N(i)->Idx();
                   if(current_n_id!=current_el.N(i)->Idx()){
                       const double64 adv_var(current_el.N(i)->Read( this->adv1_key_ ));
                       (*sit).first  = std::min( (*sit).first,  adv_var );
                       (*sit).second = std::max( (*sit).second, adv_var );
                   }
               }
          }

    } // end for all nodes

 } // end MinMaxAdvectedPropertyIncludingTheCurrentNode

template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::MinMaxAdvectedPropertyExceptTheCurrentNode()
{
   typename vector<pair<double64,double64> >::iterator  sit(this->SMINMAX.begin());
   Element<dim>     current_el;
   size_t  current_n_id,global_neighb_el_id;
   for ( typename vector<Node<dim>*>::const_iterator
         nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++, sit++ ) {

         current_n_id= (*(*nit)).Idx();
         // 1. the advected property value at the current node is assigned to min-max pair
         (*sit).first =std::numeric_limits<double64>::max();
         (*sit).second = std::numeric_limits<double64>::min();

         for(  size_t p = 0; p< (*nit)->Parents() ; p++ ){
               // get the global parent id:
               global_neighb_el_id = (*nit)->Parent( p)->Idx();
               // get the corresponding element:
               current_el = *this->gref_.E( global_neighb_el_id );
               for(size_t i=0;i<current_el.Nodes();i++){
                   //ids[i]=current_el.N(i)->Idx();
                   if(current_n_id!=current_el.N(i)->Idx()){
                       const double64 adv_var(current_el.N(i)->Read( this->adv1_key_ ));
                       (*sit).first  = std::min( (*sit).first,  adv_var );
                       (*sit).second = std::max( (*sit).second, adv_var );
                   }
               }
          }

    } // end for all nodes

 } // end MinMaxAdvectedPropertyExceptTheCurrentNode






/**
     Rather than using the sector integrals of the source term interpolated externally, this method integrates
     sources by multiplication with finite volume (the pore volume)
     
     @author who?
*/
template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::ComputePiecewiseConstantNodalSource()
 {
     size_t          n(0U);
     const double64  tolerance(1000. * numeric_limits<double64>::epsilon());

     for ( typename vector<Node<dim>*>::const_iterator
         nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++, n++ ) {
         /// Reading the source term
         double64 nodal_source = (*nit)->Read( this->src_key_ );
         /// Integrating it by multiplication with the finite volume volume
         nodal_source *= this->FVPOREVOL[n];
         /// Assigning it to diagonal of solution matrix if this a sink we just put it into diagonal of lefthandside
         if ( nodal_source < 0. ) advector_.AddToLHS( n, n, -nodal_source );
         else {
             /// Check whether we could fit it into the diagonal of LHS without loosing positive definiteness
             const double64 lhs_diagonal_val = advector_.LHS_Value( n, n );
             if ( nodal_source < (lhs_diagonal_val - tolerance) ) advector_.AddToLHS( n, n, -nodal_source );
             else {
             /// Put as much of the source as possible into the LHS and the rest into the RHS
             const double64 lhs_contribution = nodal_source - (lhs_diagonal_val - tolerance);
             advector_.AddToLHS( n, n, -lhs_contribution );
             /// What is left goes into RHS
             const double64 rhs_contribution = (nodal_source - lhs_contribution) * (*nit)->Read( this->adv1_key_ );
             advector_.AddToRHS( n, rhs_contribution );
             }
         }
     }
 } // ComputePiecewiseConstantNodalSource





template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::AddSourceTerm()
 {
     size_t  n(0U);
     for ( typename vector<Node<dim>*>::const_iterator
         nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++, n++ ) {
         /// Reading the source term
         double64 nodal_source = (*nit)->Read( this->src_key_ );
         /// Integrating it by multiplication with the finite volume volume
         nodal_source *= this->FVPOREVOL[n];
         if(nodal_source!=0.){
             advector_.AddToRHS( n,-nodal_source );
         }
     }
 } // AddSourceTerm









/**

SolveTransportEquation1stOrderInNonConservativeForm() advects the non-wetting
phase in a two-phase flow using the user-supplied relative permeability
model. The upstream values of the advected variable are used so that
the scheme is unconditionally stable as long as the CFL condition is
obeyed.

The derivatives of the fractional flow function are computed at each
facet integration point. Thus, there is a dependency of mesh refinement
for the position of the saturation front, i.e. the scheme is more
diffusive where the mesh is coarser.

@section arguments Input Arguments

The method needs access via reference to the current simulation model,
the ExplicitFiniteVolumeAlgorithm that manages the result vector, and the
specific relative permeability model which will be a subclass of
AdvectionDiffusionTwoPhaseModel.

The input time-increment must satisfy the CFL condition.

@section implementation Implementation

The divergence of the flow field as a consequence of advecting the
non-wetting phase only, is treated as a source - sink term that is
compensated for in each FV. Transient flow -related sources and sinks
are also considered.

@section application Application

Use the second-order accurate scheme wherever you can since it takes
almost the same effor to compute in contrast to the conservative case,
where no interpolations of the transported variable have to be
performed by the first-order accurate scheme.

Takes into account the "nodal fluid volume source" as a source of the transported phase
distributed over the node-centered finite volume.
 */
template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SolveTransportEquation1stOrder( TwoPhaseModel<dim>& relperm,
                                                                                       double64 time_increment,
                                                                                       bool account_for_nodal_sources )
 {
    vector<FV_Parameter>::const_iterator  fvt=this->STENCIL_DATA.begin();

    for ( typename vector<Element<dim>*>::const_iterator
          eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
      {
         // 2.1 get necessary data from each element
         this->stencil_.AccumulateImplicitTwoPhaseSolution1( (*fvt), *(*eit), relperm );

         // 2.2 first-order saturations are accumulated into lefthand side
         advector_.AccumulateLHS( this->stencil_, time_increment );

         // 2.3 source terms due a divergence of vt
         advector_.AccumulateSectorSourceTermsInLHS( this->stencil_ );

         // 2.4 conductance matrix for diffusion
         if ( with_capillary_spreading_ ) advector_.AccumulateIntegral_DNT_op_DN_dV_LHS( this->stencil_ );

         // 2.5 previous solution multiplied by storage and divided by time increment is
         //     accumulated into righthandside for Backward-Euler time stepping
         advector_.AccumulateRHS( this->stencil_, time_increment );

      } // end of accumulation

    // 3. compute nodal source terms adding them to the diagonal of the LHS of the solution matrix
    if ( account_for_nodal_sources ) ComputePiecewiseConstantNodalSource();

    // 4. algebraic multigrid solver is applied to compute FV saturations
    advector_.SolveMatrixEquation();

    // 5. Saving the computed new saturations at the FV centers
    const size_t  n_result_phase(2U); // always phase is the one that gets advected
    advector_.OutputResults( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true );

    // 6. empty the sparse matrix and righthand vector
    advector_.ResetLHS( this->gref_.Nodes() );
    advector_.ResetRHS( this->gref_.Nodes() );

 } // end SolveTransportEquation1stOrder


template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SolveTransportEquation1stOrder_NonlinearNewtonRaphson( TwoPhaseModel<dim>& relperm,
                                                                                                              double64 time_increment )
 {
    const size_t  n_result_phase(2U); // always phase is the one that gets advected

    size_t iter(1U);
    double64 LineSearchIter(1.0),ls_res_prev(1.0),ls_res(1.0);
    double64 res(numeric_limits<double64>::max());
    bool LineSearchContinue(false);
    vector<double64> scale_vec;

    // -------------------------------------------------------
    // Setting for Nonlinear loop

    double64 ls_res_target(0.99*this->target_newton_raphson_residual_);
    // -------------------------------------------------------

    //Initialization section
    double64 smin(0.0), smax(0.0);
    this->SAT0.resize( this->gref_.Nodes() );
    vector<double64>(this->SAT0).swap(this->SAT0);
    this->SN_.resize( this->gref_.Nodes() );
    vector<double64>(this->SN_).swap(this->SN_);
    this->DS_.resize( this->gref_.Nodes() );
    vector<double64>(this->DS_).swap(this->DS_);

    this->InitialAdvectedPropertyValues(this->adv1_key_,smin,smax);

    while( (res > this->target_newton_raphson_residual_) && (iter < this->max_newton_raphson_iterations_)){

        // 1. empty the sparse matrix and righthand vector
        advector_.ResetLHS( this->gref_.Nodes() );
        advector_.ResetRHS( this->gref_.Nodes() );

        //#######################################################################################################################################

        vector<FV_Parameter>::const_iterator  fvt=this->STENCIL_DATA.begin();

        for ( typename vector<Element<dim>*>::const_iterator
              eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
        {

            // 2.1 get necessary data from each element
            this->stencil_.AccumulateImplicitTwoPhaseSolution1_NonlinearNewtonRaphson((*fvt), *(*eit), relperm, with_gravitational_forces_, with_capillary_spreading_);

            // 2.2 first-order saturations are accumulated into lefthand side
            advector_.AccumulateMatrix_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment);

            // 2.3 first-order saturations are accumulated into lefthand side
            advector_.AccumulateResidual_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment);

        } // end of accumulation


        // Calculations at Boundary
        // loop over the boundary cells and adjust fluxes

        for ( size_t i=this->gref_.InteriorNodes(); i<this->gref_.Nodes(); i++ )
        {
             const Node<dim>* const nd_ptr=this->gref_.N(i);
             size_t nid=this->gref_.N(i)->Idx();
             double64        flux_balance(0.0),inflow(0.0);

             if((nd_ptr->Status(reference_variable_to_no_flow_bc_key_) == DIRICH)&&(nd_ptr->Status(this->adv1_key_)!=DIRICH)){

                 // for all SECTORS of the FE_FV-stencils which contribute to boundary finite volume (surrounding the node)
                 for ( size_t t=0U; t<nd_ptr->Parents(); t++ )
                 {
                     Element<dim>* const eptr(nd_ptr->Parent(t));
                     size_t pnid(this->gref_.N(i)->ParentNodeNumber(t));
                     //if(((dim==3) && (!eptr->FE()->IsSurfaceElement()) && (!eptr->FE()->IsLineElement())) ||((dim==2) && (!eptr->FE()->IsLineElement()))){
                     //if((dim>1) && (!eptr->FE()->IsLineElement())){

                          double64 flux(0.0);
                          this->stencil_.CorrectImplicitTwoPhaseSolutionAtBoundary_NonlinearNewtonRaphson(this->STENCIL_DATA[eptr->Idx()], *eptr, relperm, pnid, flux, with_gravitational_forces_, with_capillary_spreading_);
                          // inflow and outflow are measured using the stencils in the interior of the computational region
                          // thus inflows to model originate as positive and outflows as negative
                          inflow += flux; // +to satisfy convention above (that inflow is positive)
                          // the balance can only be evaluated if there is a halo stencil
                          if(!IsInteriorStencil(eptr)){
                                if (halo_stencils_.find(eptr)!= halo_stencils_.end()) flux_balance -= flux;
                                else throw csmp::Exception( ERROR, "TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SolveTransportEquation1stOrder_NonlinearNewtonRaphson",
                                                         "Attempt to access a finite volume stencil that was not initialized" );
                          }

                  //} //end if statement for low dimensional elements

                 } //end for loop for parent elements

             } // end if statement for DIRICHLET BC

             advector_.CompensateInflowOutFlowBoundaries(nid,inflow,flux_balance);

        }//end for loop for interior nodes

        // 3. compute nodal source term and add it to the RHS
        AddSourceTerm();

        //#######################################################################################################################################

        // 4. algebraic multigrid solver is applied to compute FV saturations
        advector_.SolveMatrixEquation();

        // 5. Saving the computed new saturations at the FV centers

        advector_.OutputResults_NonlinearNewtonRaphson( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true , this->DS_, this->SN_);

        res=vector_norm2( advector_.RHSVector(), scale_vec );

        // ------------------------------------------------
        // Line Search Step
        // ------------------------------------------------

        // 6. Do line search step

        LineSearchContinue=true;
        LineSearchIter=0.0;

        while(LineSearchContinue && (LineSearchIter < this->max_line_search_iterations_))
        {
            for (  typename vector<Node<dim>*>::const_iterator nit  = this->gref_.NodesBegin();
                nit != this->gref_.NodesEnd(); nit++ )
            {
                advector_.AssignRESULT((*nit)->Idx(), std::min(std::max(0.0,this->SN_[(*nit)->Idx()]+DS_[(*nit)->Idx()]*1/pow(2.0,LineSearchIter)),1.0));
            }
            LineSearchIter ++;

            // 1. Saving the computed new saturations at the FV centers=nodes
            advector_.OutputResults_NonlinearNewtonRaphson( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true );

            // 2. Recalculate residual in the new point

            //#######################################################################################################################################

            // Reset previos data in RHS
            advector_.ResetRHS( this->gref_.Nodes() );

            // Calculate fluxes and residual inside the region

            fvt=this->STENCIL_DATA.begin();
            for ( typename vector<Element<dim>*>::const_iterator
                eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
            {
                // 2.1 get necessary data from each element
                this->stencil_.AccumulateImplicitTwoPhaseSolution1_NonlinearNewtonRaphson((*fvt), *(*eit), relperm, with_gravitational_forces_, with_capillary_spreading_);

                // 2.2 first-order saturations are accumulated into lefthand side
                advector_.AccumulateResidual_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment );

            }

            // Calculate fluxes and residual at the boundary of the region
            // loop over the boundary cells and adjust fluxes

            for ( size_t i=this->gref_.InteriorNodes(); i<this->gref_.Nodes(); i++ )
            {
                 const Node<dim>* const nd_ptr=this->gref_.N(i);
                 size_t nid=this->gref_.N(i)->Idx();
                 double64        flux_balance(0.0),inflow(0.0);

                 if((nd_ptr->Status(reference_variable_to_no_flow_bc_key_) == DIRICH)&&(nd_ptr->Status(this->adv1_key_)!=DIRICH)){

                     // for all SECTORS of the FE_FV-stencils which contribute to boundary finite volume (surrounding the node)
                     for ( size_t t=0U; t<nd_ptr->Parents(); t++ )
                     {
                         Element<dim>* const eptr(nd_ptr->Parent(t));
                         size_t pnid(this->gref_.N(i)->ParentNodeNumber(t));
                         //if(((dim==3) && (!eptr->FE()->IsSurfaceElement()) && (!eptr->FE()->IsLineElement())) ||((dim==2) && (!eptr->FE()->IsLineElement()))){
                         //if((dim>1) && (!eptr->FE()->IsLineElement())){

                              double64 flux(0.0);
                              this->stencil_.CorrectImplicitTwoPhaseSolutionAtBoundary_NonlinearNewtonRaphson(this->STENCIL_DATA[eptr->Idx()], *eptr, relperm, pnid, flux, with_gravitational_forces_, with_capillary_spreading_);
                              // inflow and outflow are measured using the stencils in the interior of the computational region
                              // thus inflows to model originate as positive and outflows as negative
                              inflow += flux; // +to satisfy convention above (that inflow is positive)
                              // the balance can only be evaluated if there is a halo stencil
                              if(!IsInteriorStencil(eptr)){
                                    if (halo_stencils_.find(eptr)!= halo_stencils_.end()) flux_balance -= flux;
                                    else throw csmp::Exception( ERROR, "TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SolveTransportEquation1stOrder_NonlinearNewtonRaphson",
                                                             "Attempt to access a finite volume stencil that was not initialized" );
                              }

                      //} //end if statement for low dimensional elements

                     } //end for loop for parent elements

                 } // end if statement for DIRICHLET BC

                 advector_.CompensateInflowOutFlowBoundaries(nid,inflow,flux_balance);

            }//end for loop for interior nodes

            AddSourceTerm();

            //#######################################################################################################################################

            // 3. Do extra line search step if needed

            ls_res = vector_norm2( advector_.RHSVector(), scale_vec );
            LineSearchContinue = (ls_res > ls_res_target);

            if ((ls_res > ls_res_prev) && (LineSearchIter > 1))
            {
                for (  typename vector<Node<dim>*>::const_iterator nit  = this->gref_.NodesBegin();
                    nit != this->gref_.NodesEnd(); nit++ )
                {
                    advector_.AssignRESULT((*nit)->Idx(), std::min(std::max(0.0,SN_[(*nit)->Idx()]+DS_[(*nit)->Idx()]*1/pow(2.0,LineSearchIter-1)),1.0));
                }
                LineSearchContinue = false;
            }
            ls_res_prev = ls_res;
        }

        advector_.OutputResults_NonlinearNewtonRaphson( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true );
        res = vector_norm2( advector_.RHSVector(), scale_vec );

        res=ls_res;
        iter++;
    }

    // 7. Output Results
    advector_.OutputResults( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true );

    // empty the sparse matrix and righthand vector
    advector_.ResetLHS( this->gref_.Nodes() );
    advector_.ResetRHS( this->gref_.Nodes() );

 } // end SolveTransportEquation1stOrder_NonlinearSolver



















template<size_t dim, template<size_t> class STP>
void TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SolveTransportEquation2ndOrderInSpace_NonlinearNewtonRaphson( Model<dim>& sg,
                                                                                                                     TwoPhaseModel<dim>& relperm,
                                                                                                                     double64 time_increment )
 {
    const size_t  n_result_phase(2U); // always phase is the one that gets advected
    size_t iter(1U);
    double64 LineSearchIter(1.0),ls_res_prev(1.0),ls_res(1.0);
    double64 res(numeric_limits<double64>::max());
    bool LineSearchContinue(false);
    vector<double64> scale_vec;

    // -------------------------------------------------------
    // Setting for Nonlinear loop
    double64 ls_res_target(0.99*this->target_newton_raphson_residual_);
    // -------------------------------------------------------

    //Initialization section
    double64 smin(0.0), smax(0.0);
    this->SAT0.resize( this->gref_.Nodes() );
    vector<double64>(this->SAT0).swap(this->SAT0);
    this->SN_.resize( this->gref_.Nodes() );
    vector<double64>(this->SN_).swap(this->SN_);
    this->DS_.resize( this->gref_.Nodes() );
    vector<double64>(this->DS_).swap(this->DS_);
    this->InitialAdvectedPropertyValues(this->adv1_key_,smin,smax);
    //this->InitialAdvectedPropertyValues(this->ad1_key_);


    while( (res > this->target_newton_raphson_residual_) && (iter < this->max_newton_raphson_iterations_)){

        // 1. empty the sparse matrix and righthand vector
        advector_.ResetLHS( this->gref_.Nodes() );
        advector_.ResetRHS( this->gref_.Nodes() );

        //#######################################################################################################################################

        // 1. finding smin/smax of previous solution in the neighborhood of each FV
        this->MinMaxAdvectedPropertyIncludingTheCurrentNode();

        // 2. fully-upstream weighted version of finite-volume scheme
        // looping over all FE-FV stencils

        vector<FV_Parameter>::const_iterator  fvt=this->STENCIL_DATA.begin();

        if(this->with_lsmgrad_limiter_){
            /*
            //Element gradient of saturation
            ElementToNodePropertyVisitor<VectorVariable<dim>,dim> ElmToNode(this->pref_,gradsn_elmt_name_,gradsn_node_name_,this->gref_.Nodes());
            this->gref_.CopyGradientOfProperty_A_To_B(sn_name_,gradsn_elmt_name_);
            this->gref_.Accept(ElmToNode);
            */
            this->grad_advprop_limiter_->CalculateGenericNodalGradient();
            this->grad_advprop_limiter_->CalculateSlopeLimiter(this->SMINMAX);
            for ( typename vector<Element<dim>*>::const_iterator
                  eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
            {

                // 2.1 get necessary data from each element
                this->stencil_.AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson(this->SMINMAX,(*fvt), *(*eit), relperm, this->with_gravitational_forces_, this->with_capillary_spreading_,this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);

                // 2.2 first-order saturations are accumulated into lefthand side
                advector_.AccumulateMatrix_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment  );

                // 2.3 first-order saturations are accumulated into lefthand side
                advector_.AccumulateResidual_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment );

            } // end of accumulation
        }else{
            for ( typename vector<Element<dim>*>::const_iterator
                  eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
            {

                // 2.1 get necessary data from each element
                this->stencil_.AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson(this->SMINMAX,(*fvt), *(*eit), relperm, this->with_gravitational_forces_, this->with_capillary_spreading_);

                // 2.2 first-order saturations are accumulated into lefthand side
                advector_.AccumulateMatrix_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment );

                // 2.3 first-order saturations are accumulated into lefthand side
                advector_.AccumulateResidual_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment );

            } // end of accumulation

        }

        // Calculations at Boundary
        // loop over the boundary cells and adjust fluxes

        for ( size_t i=this->gref_.InteriorNodes(); i<this->gref_.Nodes(); i++ )
        {
             const Node<dim>* const nd_ptr=this->gref_.N(i);
             size_t nid=this->gref_.N(i)->Idx();
             double64        flux_balance(0.0),inflow(0.0);

             if((nd_ptr->Status(reference_variable_to_no_flow_bc_key_) == DIRICH)&&(nd_ptr->Status(this->adv1_key_)!=DIRICH)){

                 // for all SECTORS of the FE_FV-stencils which contribute to boundary finite volume (surrounding the node)
                 for ( size_t t=0U; t<nd_ptr->Parents(); t++ )
                 {
                     Element<dim>* const eptr(nd_ptr->Parent(t));
                     size_t pnid(this->gref_.N(i)->ParentNodeNumber(t));
                     //if(((dim==3) && (!eptr->FE()->IsSurfaceElement()) && (!eptr->FE()->IsLineElement())) ||((dim==2) && (!eptr->FE()->IsLineElement()))){
                     //if((dim>1) && (!eptr->FE()->IsLineElement())){

                          double64 flux(0.0);
                          this->stencil_.CorrectImplicitTwoPhaseSolutionAtBoundary_NonlinearNewtonRaphson(this->STENCIL_DATA[eptr->Idx()], *eptr, relperm, pnid, flux, this->with_gravitational_forces_, this->with_capillary_spreading_);
                          // inflow and outflow are measured using the stencils in the interior of the computational region
                          // thus inflows to model originate as positive and outflows as negative
                          inflow += flux; // +to satisfy convention above (that inflow is positive)
                          // the balance can only be evaluated if there is a halo stencil
                          if(!IsInteriorStencil(eptr)){
                                if (halo_stencils_.find(eptr)!= halo_stencils_.end()) flux_balance -= flux;
                                else throw csmp::Exception( ERROR, "TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SolveTransportEquation2ndOrder_NonlinearNewtonRaphson",
                                                         "Attempt to access a finite volume stencil that was not initialized" );
                          }

                  //} //end if statement for low dimensional elements

                 } //end for loop for parent elements

             } // end if statement for DIRICHLET BC

             advector_.CompensateInflowOutFlowBoundaries(nid,inflow,flux_balance);

        }//end for loop for interior nodes


        // 3. compute nodal source term and add it to the RHS
        AddSourceTerm();

        //#######################################################################################################################################

        // 4. algebraic multigrid solver is applied to compute FV saturations
        advector_.SolveMatrixEquation();

        // 5. Saving the computed new saturations at the FV centers

        advector_.OutputResults_NonlinearNewtonRaphson( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true , this->DS_, this->SN_);

        res=vector_norm2( advector_.RHSVector(), scale_vec );

        // ------------------------------------------------
        // Line Search Step
        // ------------------------------------------------

        // 6. Do line earch step

        LineSearchContinue=true;
        LineSearchIter=0.0;

        while(LineSearchContinue && (LineSearchIter < this->max_line_search_iterations_))
        {
            for (  typename vector<Node<dim>*>::const_iterator nit  = this->gref_.NodesBegin();
                nit != this->gref_.NodesEnd(); nit++ )
            {
                advector_.AssignRESULT((*nit)->Idx(), std::min(std::max(0.0,this->SN_[(*nit)->Idx()]+DS_[(*nit)->Idx()]*1/pow(2.0,LineSearchIter)),1.0));
            }
            LineSearchIter ++;

            // 1. Saving the computed new saturations at the FV centers=nodes
            advector_.OutputResults_NonlinearNewtonRaphson( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true );

            // 2. Recalculate residual in the new point

            //#######################################################################################################################################

            // Reset previos data in RHS
            advector_.ResetRHS( this->gref_.Nodes() );

            // Calculate fluxes and residual inside the region

            // 1. finding smin/smax of previous solution in the neighborhood of each FV
            this->MinMaxAdvectedPropertyIncludingTheCurrentNode();
            fvt=this->STENCIL_DATA.begin();

            // 2. fully-upstream weighted version of finite-volume scheme
            // looping over all FE-FV stencils
            if(this->with_lsmgrad_limiter_){
                /*
                //Element gradient of saturation
                ElementToNodePropertyVisitor<VectorVariable<dim>,dim> ElmToNode(this->pref_,gradsn_elmt_name_,gradsn_node_name_,this->gref_.Nodes());
                this->gref_.CopyGradientOfProperty_A_To_B(sn_name_,gradsn_elmt_name_);
                this->gref_.Accept(ElmToNode);
                */
                this->grad_advprop_limiter_->CalculateGenericNodalGradient();
                this->grad_advprop_limiter_->CalculateSlopeLimiter(this->SMINMAX);

                for ( typename vector<Element<dim>*>::const_iterator
                    eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
                {
                    // 2.1 get necessary data from each element
                    this->stencil_.AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson(this->SMINMAX,(*fvt), *(*eit), relperm, this->with_gravitational_forces_, this->with_capillary_spreading_,this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);

                    // 2.2 first-order saturations are accumulated into lefthand side
                    advector_.AccumulateResidual_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment );

                }
            }else{
                for ( typename vector<Element<dim>*>::const_iterator
                    eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
                {
                    // 2.1 get necessary data from each element
                    this->stencil_.AccumulateImplicitTwoPhaseSolution2_NonlinearNewtonRaphson(this->SMINMAX,(*fvt), *(*eit), relperm, this->with_gravitational_forces_, this->with_capillary_spreading_);

                    // 2.2 first-order saturations are accumulated into lefthand side
                    advector_.AccumulateResidual_NonlinearNewtonRaphson( this->stencil_, this->SAT0, time_increment );

                }
            }


            // Calculate fluxes and residual at the boundary of the region
            // loop over the boundary cells and adjust fluxes

            for ( size_t i=this->gref_.InteriorNodes(); i<this->gref_.Nodes(); i++ )
            {
                 const Node<dim>* const nd_ptr=this->gref_.N(i);
                 size_t nid=this->gref_.N(i)->Idx();
                 double64        flux_balance(0.0),inflow(0.0);

                 if((nd_ptr->Status(reference_variable_to_no_flow_bc_key_) == DIRICH)&&(nd_ptr->Status(this->adv1_key_)!=DIRICH)){

                     // for all SECTORS of the FE_FV-stencils which contribute to boundary finite volume (surrounding the node)
                     for ( size_t t=0U; t<nd_ptr->Parents(); t++ )
                     {
                         Element<dim>* const eptr(nd_ptr->Parent(t));
                         size_t pnid(this->gref_.N(i)->ParentNodeNumber(t));
                         //if(((dim==3) && (!eptr->FE()->IsSurfaceElement()) && (!eptr->FE()->IsLineElement())) ||((dim==2) && (!eptr->FE()->IsLineElement()))){
                         //if((dim>1) && (!eptr->FE()->IsLineElement())){

                              double64 flux(0.0);
                              this->stencil_.CorrectImplicitTwoPhaseSolutionAtBoundary_NonlinearNewtonRaphson(this->STENCIL_DATA[eptr->Idx()], *eptr, relperm, pnid, flux, this->with_gravitational_forces_, this->with_capillary_spreading_);
                              // inflow and outflow are measured using the stencils in the interior of the computational region
                              // thus inflows to model originate as positive and outflows as negative
                              inflow += flux; // +to satisfy convention above (that inflow is positive)
                              // the balance can only be evaluated if there is a halo stencil
                              if(!IsInteriorStencil(eptr)){
                                    if (halo_stencils_.find(eptr)!= halo_stencils_.end()) flux_balance -= flux;
                                    else throw csmp::Exception( ERROR, "TwoPhaseImplicitNodeCenteredFVTransport<dim,STP>::SolveTransportEquation2ndOrder_NonlinearNewtonRaphson",
                                                             "Attempt to access a finite volume stencil that was not initialized" );
                              }

                      //} //end if statement for low dimensional elements

                     } //end for loop for parent elements

                 } // end if statement for DIRICHLET BC

                 advector_.CompensateInflowOutFlowBoundaries(nid,inflow,flux_balance);

            }//end for loop for interior nodes


            AddSourceTerm();

            //#######################################################################################################################################

            // 3. Do extra line search step if needed

            ls_res = vector_norm2( advector_.RHSVector(), scale_vec );
            LineSearchContinue = (ls_res > ls_res_target);

            if ((ls_res > ls_res_prev) && (LineSearchIter > 1))
            {
                for (  typename vector<Node<dim>*>::const_iterator nit  = this->gref_.NodesBegin();
                    nit != this->gref_.NodesEnd(); nit++ )
                {
                    advector_.AssignRESULT((*nit)->Idx(), std::min(std::max(0.,SN_[(*nit)->Idx()]+DS_[(*nit)->Idx()]*1/pow(2.,LineSearchIter-1)),1.));
                }
                LineSearchContinue = false;
            }
            ls_res_prev = ls_res;
        }

        advector_.OutputResults_NonlinearNewtonRaphson( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true );
        res = vector_norm2( advector_.RHSVector(), scale_vec );

        res=ls_res;
        iter++;
    }

    // 7. Output Results
    advector_.OutputResults( this->pref_, n_result_phase, ph1_key, this->adv1_key_, true );

    // empty the sparse matrix and righthand vector
    advector_.ResetLHS( this->gref_.Nodes() );
    advector_.ResetRHS( this->gref_.Nodes() );

} // end SolveTransportEquation2ndOrderInSpace_NonlinearSolver



template class TwoPhaseImplicitNodeCenteredFVTransport<1U,StencilProcessor>;
template class TwoPhaseImplicitNodeCenteredFVTransport<2U,StencilProcessor>;
template class TwoPhaseImplicitNodeCenteredFVTransport<3U,StencilProcessor>;
 
} // end namespace csmp
