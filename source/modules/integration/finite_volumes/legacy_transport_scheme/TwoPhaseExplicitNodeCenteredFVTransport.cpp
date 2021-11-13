#include "TwoPhaseExplicitNodeCenteredFVTransport.h"
#include "NodeCenteredFiniteVolumeAlgorithm.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "TwoPhaseModel.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "StencilProcessor.h"
#include "ExplicitStencilProcessor.h"
#include "STL_utilities.h"
#include "compute2PhaseMobilityAtBaryCenter.h"

using namespace std;

namespace csmp {

/**
 
Constructor for modelling transport by two-phase slightly compressible flow 
with fluid density and viscosity variations. The properties 
fluid density and viscosity are expected to be node variables and
are interpolated to the FV facet integration points to carry out the
computations. The transport scheme is setup for the non-wetting 
phase. The onus is on the user to update the second phase saturations
after each advection step, i.e., calculate sw=1-sn.  

The TwoPhaseExplicitNodeCenteredFVTransport can be constructed for first- 
and second order calculations.  

Diffusion is only taken into account if second-order accuracy is chosen because
the first-order scheme already is so diffusive that adding an extra
diffusion term makes little sense.   

The TwoPhaseExplicitNodeCenteredFVTransport gains access via reference to the 
model to which it is applied to. The user needs to specify the piecewise
constant element variables porosity, capillary diffusivity and transport velocity.
These must be defined in the CSP_variables.txt file. The variables 
fluid density and viscosity of the wetting and non-wetting phases must
be nodal variables as well. They are interpolated to the facet and sector
integration points using the finite element basis functions.  

The transported variable must be nodal since this explicit transport scheme 
is based on node-centered finite volumes.  

@section implementation Implementation 

This constructor creates a minor storage array for repetitively used
variables like the finite (pore) volumes, the flux balance for each
FV cell and a vector of FV_parameters called FV_stencil_data that holds
the FV sector volumes, facet areas and facet normal fluxes which are 
always updated before the first advection step, when the CFL criterion
is computed. This vector will become redundant when the generic FV scheme
has been optimized for speed.  

@section application  Application 

Use TwoPhaseExplicitNodeCenteredFVTransport as constructed with this particular
constructor for two-phase flow in a rapidly evolving
flow field that necessitates frequent updates of the velocity fields.  

For a weak coupling between phase saturations and total mobility, use 
the implicit transport module FiniteVolumeTransport.  

Capillary diffusion will not correctly model transport across material
interfaces unless these are split.
 
*/
template<size_t dim, template<size_t> class STP>
TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::TwoPhaseExplicitNodeCenteredFVTransport( // 2-phase
                                                                                const char* group_name,
                                                                                Model<dim>& sg,
                                                                                const char* porosity,
                                                                                const char* diffusivity,
                                                                                const char* advected_prop1,
                                                                                const char* advected_prop2,
                                                                                const char* transp_velocity,
                                                                                const char* nodal_source,
                                                                                bool second_order_accurate,
                                                                                bool with_capillary_spreading,
                                                                                bool with_gravitational_forces,
                                                                                const char* reference_variable_to_no_flow_bc,
                                                                                const char* thickness_attribute )
 : NodeCenteredFiniteVolumeTransport<dim>( group_name, sg,
                                           porosity, diffusivity,
                                           advected_prop2, transp_velocity,
                                           nodal_source, second_order_accurate, false, thickness_attribute ),
   stencil_( this->adv1_key_, this->vel_key_, this->diff_key_ ),
   RESULT(this->gref_.Nodes()),
   sw_key_(sg.Database().StorageKey(advected_prop1)),
   with_capillary_spreading_(with_capillary_spreading),
   with_gravitational_forces_(with_gravitational_forces),
   no_flow_bc_key_(sg.Database().StorageKey(reference_variable_to_no_flow_bc))
 {
    cout <<"\n\nTwoPhaseExplicitNodeCenteredFVTransport<"<< dim;
    cout <<">(constructor - multiphase fluid flow): Constructed successfully."<< endl;

    // establishing the halo stencils
    for ( typename vector<Node<dim>*>::const_iterator
          nit=this->gref_.NodesBegin(); nit!=this->gref_.NodesEnd(); nit++ )
      for ( size_t i=0U; i<(*nit)->Parents(); i++ )
        if ( !IsInteriorStencil( (*nit)->Parent(i) ) )
          halo_stencils_.insert( (*nit)->Parent(i) );
    
 } // end constructor (multiphase fluid flow)


/// reports whether an element stencil has a face at the boundary of the advection region
template<size_t dim, template<size_t> class STP>
inline bool TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::IsInteriorStencil(
                                                      const Element<dim>* const eptr ) const
 {
    return this->gref_.IsPerimeterElement( eptr );
 }


template<size_t dim, template<size_t> class STP>
TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::~TwoPhaseExplicitNodeCenteredFVTransport()
 {
 }




/**
    Capillary spreading is modelled by default.
    This method disables it speeding up computation for large-scale models.
*/
template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::DisableCapillarySpreading()
  {
     with_capillary_spreading_ = false;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::EnableCapillarySpreading()
  {
     with_capillary_spreading_ = true;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::DisableGravitationalForces()
  {
     with_gravitational_forces_ = false;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::EnableGravitationalForces()
  {
     with_gravitational_forces_ = true;
  }

template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::SetNoFlowBoundaryConditionKey(Model<dim>& sg,const char* no_flow_bc)
  {
     no_flow_bc_key_ = sg.Database().StorageKey(no_flow_bc);
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
double64 TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::AnisotropicCourantIncrement( TwoPhaseModel<dim>& relperm,
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
           velocity = std::max( velocity, velocity * relperm.MaxFractionalFlowDerivative() );
           // for implicit scheme is better:
           //velocity *= relperm.dfds();

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
         cout <<"\nINFO, TwoPhaseImplicitNodeCenteredFVTransport::AnisotropicCourantIncrement (2-phase flow)";
         cout <<" calculated courant increment is larger than maximum permitted increment, there may be no flow in the model domain";
         cout <<"\nCFL is set to "<< max_time_increment << endl;
         return max_time_increment;
      }
    if ( courant_increment <= millisecond ) {
         throw csmp::Exception( WARNING, "TwoPhaseImplicitNodeCenteredFVTransport::AnisotropicCourantIncrement (2-phase flow)",
                                         "courant increment is smaller than a millisecond. Check your boundary conditions" );
      }
    else {
         cout <<"\nTwoPhaseImplicitNodeCenteredFVTransport::AnisotropicCourantIncrement (2-phase flow): ";
         cout << courant_increment <<" and CFL multiplier: "<< this->cfl_multiplier_ << endl;
      }

    return courant_increment;

 } // end AnisotropicCourantIncrement (2-phase flow - all)




/**

The method is the equivalent to AdvectVariable() for 2-phase flow. 
It advects the non-wetting phase saturation for the user-specified time interval
and writes the results back to the current Model object.  
The advected variable field is written back into the supplied Model.

@return The method returns the CFL condition that it determined before the 
onset of the advection.

@attention After every sub-CFL transport step, the method prints a dot on the screen. 
  */
template<size_t dim, template<size_t> class STP>
double64  TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::TransportPhase( TwoPhaseModel<dim>& relperm, 
                                                                            double64 time_interval )
 {
    this->gref_.RenumberNodes();

    // true = uses only the result vector
    NodeCenteredFiniteVolumeAlgorithm<dim>  advector( this->gref_ );

    // 1. compute the CFL condition to identify value for overstepping
    // ---------------------------------------------------------------

    const double64 MAX_TIME(86400. * 365.);
    const double64 courant_increment = this->AnisotropicCourantIncrement( relperm, MAX_TIME );
                                 
    this->cfl_multiplier_ = (with_capillary_spreading_==true) ? 0.1 : 0.4;

    double64  time(0.), time_increment = std::min( this->cfl_multiplier_ * courant_increment, time_interval );

    cout<<"\n Time interval="<<time_interval<<endl;
    cout<<" Courant increment="<<courant_increment<<endl;
    cout<<" Time increment="<<time_increment<<endl;

    // 2. compute solution
    // -------------------
    cout <<"\n\n\nTwoPhaseExplicitNodeCenteredFVTransport<"<< dim;
    cout <<">::TransportPhase: Advecting phase";
    while ( time < time_interval ) {

         if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;

         if ( this->SecondOrderInSpace() )
               AdvectVariable2ndOrder( relperm, time_increment );
         else
               AdvectVariable1stOrder( relperm, time_increment );

         time += time_increment;

         //cout <<".";
         cout.flush();

    }

    return courant_increment;

 } // end TransportPhase


/// Copies the result vector to model performing range checks against the PropertyDatabase.
template<size_t dim, template<size_t> class STP>
double64 TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::OutputResults( const csmp::Index& adv_key,
                                                                          bool show_range ) const
 {
    double64         rmin, rmax,
                     amin = RESULT[0],
                     amax = RESULT[0],
                     difference_to_last_output(0.);
    size_t           error_counter(0);
    ScalarVariable   sc;

    this->pref_.RangeOf( this->pref_.Name(adv_key), rmin, rmax );
    for ( size_t i=0U; i<RESULT.size(); i++ ) {
         // recording output range
         amin = std::min( amin, RESULT[i] );
         amax = std::max( amax, RESULT[i] );
           if ( this->gref_.N(i)->Status( adv_key ) != DIRICH ) {
             // reading the pre-existing value and calculating the maximum change per node
             this->gref_.N(i)->Read( adv_key, sc );
             difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc()) );
             // result checking
             if ( RESULT[i] <= rmax && RESULT[i] >= rmin )
               this->gref_.N(i)->Store( adv_key, sc=RESULT[i] );
             else {
                  cout <<"\nvalue: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
                  if ( RESULT[i] > rmax )
                    this->gref_.N(i)->Store( adv_key, sc=rmax );
                  else if ( RESULT[i] < rmin )
                    this->gref_.N(i)->Store( adv_key, sc=rmin );
                  error_counter++;
               }
            }
      }

    // reporting problems
    if ( error_counter > 20U ) {
          throw csmp::Exception( ERROR, "TwoPhaseExplicitNodeCenteredFVTransport::OutputResults",
                                "Output property was out of range, legal (min/max) was stored instead");
      }
    if ( error_counter > (this->gref_.Nodes() / 20U) )
      throw out_of_range("TwoPhaseExplicitNodeCenteredFVTransport::OutputResults: Advected variable out of range");

    if ( show_range ) {
         cout <<"\n\nTwoPhaseExplicitNodeCenteredFVTransport<"<< dim;
         cout <<">::OutputResults: Variable range after advection: ";
         cout << amin <<" to "<< amax << endl << endl;
      }

    return difference_to_last_output / std::max( amax - amin, 1.0e-20 );

 } // end OutputResults



/// Copies the result vector and the saturation of the other phase to model performing range checks against the PropertyDatabase.
template<size_t dim, template<size_t> class STP>
double64 TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::OutputResults( const csmp::Index& adv1_key,
                                                                          const csmp::Index& adv2_key,
                                                                          bool show_range,
                                                                          bool do_range_check) const
 {
    double64        rmin, rmax,
                    amin = RESULT[0],
                    amax = RESULT[0],
                    difference_to_last_output(0.);
    size_t          error_counter(0);
    ScalarVariable  sc;

    this->pref_.RangeOf( this->pref_.Name(adv2_key), rmin, rmax );
      for ( size_t i=0U; i<RESULT.size(); i++ )
         {
            // recording output range
            amin = std::min( amin, RESULT[i] );
            amax = std::max( amax, RESULT[i] );
              if ( this->gref_.N(i)->Status( adv2_key ) != DIRICH ) {
                // reading the pre-existing value and calculating the maximum change per node
                this->gref_.N(i)->Read( adv2_key, sc );
                difference_to_last_output = std::max( difference_to_last_output, fabs(RESULT[i]-sc()) );
                // result checking
                if ( RESULT[i] <= rmax && RESULT[i] >= rmin ) {
                     this->gref_.N(i)->Store( adv2_key, sc=RESULT[i] );
                     this->gref_.N(i)->Store( adv1_key, sc=1.-RESULT[i] );
                    }
                else {
                     cout <<"\nvalue: "<< RESULT[i] <<" versus range from PropertyDatabase: "<< rmin <<"-"<< rmax << endl;
                     if ( RESULT[i] > rmax ) {
                          this->gref_.N(i)->Store( adv2_key, sc=rmax );
                          this->gref_.N(i)->Store( adv1_key, sc=1.-rmax );
                         }
                     else if ( RESULT[i] < rmin ) {
                           this->gref_.N(i)->Store( adv2_key, sc=rmin );
                           this->gref_.N(i)->Store( adv1_key, sc=1.-rmin );
                       }
                     error_counter++;
                  }
              }
         }

    if( do_range_check){

        if ( error_counter > 20U ) {
              throw csmp::Exception( ERROR, "TwoPhaseExplicitNodeCenteredFVTransport::OutputResults",
                                         "Output property was out of range, legal (min/max) was stored instead");
          }
        if ( error_counter > (this->gref_.Nodes() / 20U) )
          throw out_of_range("TwoPhaseExplicitNodeCenteredFVTransport::OutputResults: Advected variable out of range");

    }

    if ( show_range ) {
         cout <<"\n\nTwoPhaseExplicitNodeCenteredFVTransport<"<< dim;
         cout <<">::OutputResults: Variable range after advection: ";
         cout << amin <<" to "<< amax << endl << endl;
      }

    return difference_to_last_output / std::max( amax - amin, 1.0e-20 );

 } // end OutputResults



template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::MinMaxAdvectedProperty()
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
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::MinMaxAdvectedPropertyExceptTheCurrentNode()
{
   typename vector<pair<double64,double64> >::iterator  sit(this->SMINMAX.begin());
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
               Element<dim>* current_el = this->gref_.E( global_neighb_el_id );
               for(size_t i=0;i<current_el->Nodes();i++){
                   //ids[i]=current_el.N(i)->Idx();
                   if(current_n_id!=current_el->N(i)->Idx()){
                       const double64 adv_var(current_el->N(i)->Read( this->adv1_key_ ));
                       (*sit).first  = std::min( (*sit).first,  adv_var );
                       (*sit).second = std::max( (*sit).second, adv_var );
                   }
               }
          }

    } // end for all nodes

 } // end MinMaxAdvectedPropertyExceptTheCurrentNode


template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::MinMaxAdvectedPropertyIncludingTheCurrentNode()
{
   typename vector<pair<double64,double64> >::iterator  sit(this->SMINMAX.begin());
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
               Element<dim>* current_el = this->gref_.E( global_neighb_el_id );
               for(size_t i=0;i<current_el->Nodes();i++){
                   //ids[i]=current_el.N(i)->Idx();
                   if(current_n_id!=current_el->N(i)->Idx()){
                       const double64 adv_var(current_el->N(i)->Read( this->adv1_key_ ));
                       (*sit).first  = std::min( (*sit).first,  adv_var );
                       (*sit).second = std::max( (*sit).second, adv_var );
                   }
               }
          }

    } // end for all nodes

 } // end MinMaxAdvectedPropertyIncludingTheCurrentNode


/**

Computes the non-wetting phase saturation at the new time-level 
for 2-phase flow.  

Solves  S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]

- a single loop over all nodes is required
- this steps also considers in and outflow of finite volume cells
- (distributed) sources and sinks will be considered in the future

Outfluxes at outflow boundaries 
and poroelastic sources and sinks during transient fluid flow are  
compensated for through the consideration of a flux balance vector that 
is supplied as a method argument.  

@param time_increment  the current time-increment for
the advection step needs to be specified in order to compose the 
solution. This time increment must satisfy the CFL condition because 
this is an explicit transport scheme.  

@section application  Application 

This method is used by the higher order-in space accurate transport scheme.  
 */
template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::Compose2PhaseSolution(
                                                               TwoPhaseModel<dim>& relperm,
                                                               double64 time_interval,
                                                               bool with_divergence_free_correction)
 {
     for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ ) 
       {
         if( with_divergence_free_correction){

          // 1. A FLUX_BALANCE computation is performed but only if the node is not at the model boundary
          if ( this->gref_.IsPerimeterNode(nidx) and 
               fabs(this->FLUX_BALANCE[nidx]) > numeric_limits<double64>::epsilon() )
            {
              // computing the average fractional flow for the current finite volume
              double64 fn_avg = static_cast<double64>(0.);
              // for each finite volume, f is evaluated on a sector by sector basis
              for ( size_t t=0U; t<this->gref_.N(nidx)->Parents(); t++ ) {
                   const size_t nid(this->gref_.N(nidx)->ParentNodeNumber(t));
                   // the relperm model is initialised for the saturation of the current finite volume
                   relperm.Initialize( *this->gref_.N(nidx)->Parent(t) );
                   relperm.InitializeForNode( *this->gref_.N(nidx)->Parent(t), nid );
                   relperm.EffectiveSaturation();
                   // weighted for each specific FV sector
                   fn_avg += relperm.f_Phase(2U);
                }
               fn_avg /= static_cast<double64>(this->gref_.N(nidx)->Parents());
               
               // to compensate for any non-physical or physical divergence, the flux balance must be substracted
               // (note that is opposite to implicit scheme where this enters the LHS) 
               RESULT[nidx] -= this->FLUX_BALANCE[nidx] * fn_avg;
           }

         }

         // 2. ACCUMULATION: subtracting the flux time-interval product
         RESULT[nidx] = this->gref_.N(nidx)->Read( this->adv1_key_ ) -
                        (time_interval / this->FVPOREVOL[nidx]) * RESULT[nidx];

         // 3. Nodal fluid SOURCE/SINK terms are considered now
         RESULT[nidx] += time_interval * this->gref_.N(nidx)->Read( this->src_key_ );


       }
       
} // end Compose2PhaseSolution                                                 



/**

For finite volumes corresponding at the model boundary,
the flow across the boundary is calculated also establishing
a flux balance if the finite volume extends beyond the boundary
so that this balance can be used to split of a source term if
necessary.

The flux balance is equal to the divergence of the flux for the current
finite volume (if there is more outflow, the flux balance and divergence
are greater than zero).

The flux balance is calculated from the element stencils that the node
is connected to. These give respective sector contributions.

To get model in and outfluxes the interior stencils are used.

@return The method returns the inflow into the model domain and whether the flux
balance was evaluated as indicated by the state of a bool. The convention adopted
here is that inflow is negative while outflow is positive.

The flux_balance can be evaluated only if the finite volume lies inside
the model domain so that it is not truncated and the exterior part can also be
used. In this case the method returns true, else false.

@attention Use this method only on finite volumes that are located on the boundary
of the computational domain.

*/
template<size_t dim, template<size_t> class STP>
bool TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::FractionalFlowThroughBoundaryFiniteVolume(
                                                           const Node<dim>* nd_ptr,
                                                           TwoPhaseModel<dim>& relperm,
                                                           double64& inflow, double64& flux_balance ) const
 {
    inflow = flux_balance = static_cast<double64>(0.);

    // for all SECTORS of the FE_FV-stencils which contribute to boundary finite volume (surrounding the node)
    for ( size_t t=0U; t<nd_ptr->Parents(); t++ ) {
         Element<dim>* const eptr(nd_ptr->Parent(t));
         const size_t nid(nd_ptr->ParentNodeNumber(t));
         double64  flux(0.);

         // k, swr, snr initialisation precedes all other steps
         relperm.Initialize( *eptr );
         // now the saturation dependent properties are computed
         relperm.InitializeForNode( *eptr, nid );
         relperm.EffectiveSaturation();

         // for all FACETS per SECTOR surrounding the finite volume at the boundary
         for ( size_t i=0U; i<eptr->FV()->FacetsPerSector(nid); i++ )
           {
              size_t iFacet( eptr->FV()->FacetSurroundingSector(nid,i) );
              size_t inside_node = eptr->FV()->InsideNode( iFacet );
              // get velocity across FV facet
              double64 velo = this->STENCIL_DATA[eptr->Idx()].FacetNormalVelocity(iFacet) *
                              this->STENCIL_DATA[eptr->Idx()].FacetArea(iFacet);
              // using upstream weighted first-order fluxes (gravity is not considered)
              if ( nid == inside_node )
                flux += velo * relperm.f_Phase(2U);
              else  // outside node - outgoing flux
                flux -= velo * relperm.f_Phase(2U);
           }
         // inflow and outflow are measured using the stencils in the interior of the computational region
         // thus inflows to model originate as positive and outflows as negative
         // if this is an interior stencil
         if ( IsInteriorStencil(eptr) ) inflow += flux; // +to satisfy convention above (that inflow is positive)
         // the balance can only be evaluated if there is a halo stencil
         else if ( halo_stencils_.find(eptr) != halo_stencils_.end() ) flux_balance -= flux;
         else
         throw csmp::Exception( ERROR, "TwoPhaseExplicitNodeCenteredFVTransport<dim>::FractionalFlowThroughBoundaryFiniteVolume",
                                        "Attempt to access a finite volume stencil that was not initialized" );
      }

    // give an indication whether the flux balance was evaluated
    if ( this->gref_.IsPerimeterNode(nd_ptr->Idx()) ) return false;

    return true;

 } // end FractionalFlowThroughBoundaryFiniteVolume



/// Compensation of in- and outflows for finite volumes that are truncated by region or model boundaries.
template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::AssignFractionalFlowBoundaryConditions( TwoPhaseModel<dim>& relperm )
 {
    double64        inflow, flux_balance;
    const double64  zero(0.);

    // loop over the boundary finite volumes and adjust fluxes
    for ( size_t i=this->gref_.InteriorNodes(); i<this->gref_.Nodes(); i++ )
      {
         const Node<dim>* const nit=this->gref_.N(i);
         // if the flux balance cannot be evaluated because the node sits at a model boundary
         if ( !FractionalFlowThroughBoundaryFiniteVolume( nit, relperm, inflow, flux_balance ) )
           {
              // inflow compensation (both in and outflow compensations are necessary in explicit scheme)
              if ( inflow > zero ) RESULT[ i ] += inflow;
              // else: outflow compensation
              else                 RESULT[ i ] -= inflow;

           }
         // if we have a flux balance because this is an internal boundary
         else {
               if ( flux_balance > zero ) // counter balancing
                  RESULT[ i ] -= flux_balance;
               else if ( flux_balance < zero )
                  RESULT[ i ] += flux_balance;
           }
      }

 } // end AssignFractionalFlowBoundaryConditions




/// Compensation of in- and outflows for finite volumes that are truncated by region or model boundaries
template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::AssignGenericFlowBoundaryConditions(TwoPhaseModel<dim>& relperm )
 {
    double64        flux_balance,inflow;
    Point<dim>      dsdn;
    DenseMatrix<DM_MIN>    DN;

    const size_t v( (dim==1u) ? 0u : 1u );
    VectorVariable<dim> velo;
    double64 facetArea;

    // loop over the boundary finite volumes and adjust fluxes
    for ( size_t i=this->gref_.InteriorNodes(); i<this->gref_.Nodes(); i++ )
      {
         const Node<dim>* const nd_ptr=this->gref_.N(i);

         if((nd_ptr->Status(no_flow_bc_key_) == DIRICH)&&(nd_ptr->Status(this->adv1_key_)!=DIRICH)){
             // if the flux balance cannot be evaluated because the node sits at a model boundary
             flux_balance = static_cast<double64>(0.);
             inflow = static_cast<double64>(0.);

             // for all SECTORS of the FE_FV-stencils which contribute to boundary finite volume (surrounding the node)
             for ( size_t t=0U; t<nd_ptr->Parents(); t++ )
                 {
                      double64  flux(0.);
                      Element<dim>* e(nd_ptr->Parent(t));
                      const size_t pnid(nd_ptr->ParentNodeNumber(t));
                      // k, swr, snr initialisation precedes all other steps
                      relperm.Initialize( *e );
                      relperm.InitializeForNode( *e, pnid );
                      relperm.EffectiveSaturation();

                      if( this->with_capillary_spreading_ ){

                          e->dN_AtBaryCenter( DN );
                          dsdn = 0.;
                          for ( size_t j=0U; j<e->Nodes(); j++ ) {
                               const double64 sn = e->N(j)->Read( this->adv1_key_);
                               for ( size_t k=0U; k<dim; k++ ) dsdn[k] += DN(k,j) * sn;
                          }

                      }

                      // now the saturation dependent properties are computed
                      // for all FACETS per SECTOR surrounding the finite volume at the boundary
                      for ( size_t i=0U; i<e->FV()->FacetsPerSector(pnid); i++ )
                        {
                           size_t iFacet( e->FV()->FacetSurroundingSector(pnid,i) );
                           size_t inside_node,outside_node;
                           e->FV()->FacetEdgeNodes( iFacet, inside_node, outside_node );

                           // get velocity across FV facet
                           //double64 viscous_vel_component = this->STENCIL_DATA[eptr->Idx()].FacetNormalVelocity(iFacet) *this->STENCIL_DATA[eptr->Idx()].FacetArea(iFacet);
                           //double64 grav_vel_component = this->STENCIL_DATA[eptr->Idx()].FacetNormalComponent(iFacet,v) *this->STENCIL_DATA[eptr->Idx()].FacetArea(iFacet);
                           Point<dim>  n =e->FacetNormal(iFacet);
                           e->Read( this->vel_key_, velo );
                           facetArea = e->FacetArea(iFacet);

                           double64 viscous_vel_component = velo.DotProduct(n);

                           if( this->with_gravitational_forces_){

                               double64 gravity_vel_component = n[v];

                               if( this->with_capillary_spreading_ ){

                                   const double64 capillary_vel_component( -dotProduct( dsdn, n) *relperm.dpcds_Phase( ));

                                   if ( pnid == inside_node )
                                       flux += (relperm.f_Phase(2U) * viscous_vel_component  - relperm.GravityMultiplier_G()*gravity_vel_component - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;
                                   else
                                       flux -= (relperm.f_Phase(2U) * viscous_vel_component  - relperm.GravityMultiplier_G()*gravity_vel_component - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;


                               }else{

                                   if ( pnid == inside_node )
                                       flux += (relperm.f_Phase(2U) * viscous_vel_component  - relperm.GravityMultiplier_G()*gravity_vel_component )*facetArea;
                                   else
                                       flux -= (relperm.f_Phase(2U) * viscous_vel_component  - relperm.GravityMultiplier_G()*gravity_vel_component )*facetArea;

                               }

                           }else{

                               if( this->with_capillary_spreading_ ){

                                   const double64 capillary_vel_component( -dotProduct( dsdn, n)*relperm.dpcds_Phase( ));

                                   if ( pnid == inside_node )
                                       flux += (relperm.f_Phase(2U) * viscous_vel_component  - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;
                                   else
                                       flux -= (relperm.f_Phase(2U) * viscous_vel_component  - relperm.Permeability()*relperm.G()*capillary_vel_component)*facetArea;


                               }else{

                                   if ( pnid == inside_node )
                                       flux += (relperm.f_Phase(2U) * viscous_vel_component  )*facetArea;
                                   else
                                       flux -= (relperm.f_Phase(2U) * viscous_vel_component  )*facetArea;

                               }

                           }

                        }

                        // inflow and outflow are measured using the stencils in the interior of the computational region
                        // thus inflows to model originate as positive and outflows as negative
                        inflow += flux; // +to satisfy convention above (that inflow is positive)
                        // the balance can only be evaluated if there is a halo stencil
                        if(!IsInteriorStencil(e)){
                                if (halo_stencils_.find(e)!= halo_stencils_.end()) flux_balance -= flux;
                                else throw csmp::Exception( ERROR, "TwoPhaseExplicitNodeCenteredFVTransport<dim>::AssignGenericFlowBoundaryConditions",
                                                       "Attempt to access a finite volume stencil that was not initialized" );
                        }

                 }//end for loop for parent elements


                 if(this->gref_.IsPerimeterNode(nd_ptr->Idx())){
                     RESULT[ i ] -= inflow;
                     // inflow compensation (both in and outflow compensations are necessary in explicit scheme)
                     //if ( inflow > zero ) {
                     //    RESULT[ i ] += inflow;
                     //}
                     // else: outflow compensation
                     //else{
                     //    RESULT[ i ] -= inflow;
                     //}
                 }
                 // if we have a flux balance because this is an internal boundary
                 else {
                     RESULT[ i ] -= flux_balance;
                     //if ( flux_balance > zero ) // counter balancing
                     //   RESULT[ i ] -= flux_balance;
                     //else if ( flux_balance < zero )
                     //    RESULT[ i ] += flux_balance;
                 }

         }//end if statement for DIRICHLET BC

    }//end for loop for interior nodes

 } // end AssignGenericFlowBoundaryConditions


template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::DivergenceFreeCorrection(
                                                               TwoPhaseModel<dim>& relperm,
                                                               double64 time_interval )
 {
    for ( size_t nidx=0U; nidx<this->gref_.Nodes(); nidx++ ){

          const Node<dim>* const nd_ptr=this->gref_.N(nidx);
          // 1. A FLUX_BALANCE computation is performed but only if the node is not at the model boundary
          if ( !this->gref_.IsPerimeterNode(nidx) )
          {
              double64 div(0.0);
              // for each finite volume, f is evaluated on a sector by sector basis
              for ( size_t t=0U; t<this->gref_.N(nidx)->Parents(); t++ ) {
                   Element<dim>* const eptr(nd_ptr->Parent(t));
                   const size_t nid(nd_ptr->ParentNodeNumber(t));
                   relperm.Initialize( *eptr );
                   // now the saturation dependent properties are computed
                   relperm.InitializeForNode( *eptr, nid );
                   relperm.EffectiveSaturation();

                   // for all FACETS per SECTOR surrounding the finite volume at the boundary
                   for ( size_t i=0U; i<eptr->FV()->FacetsPerSector(nid); i++ ){
                        size_t iFacet( eptr->FV()->FacetSurroundingSector(nid,i) );
                        size_t inside_node,outside_node;
                        eptr->FV()->FacetEdgeNodes( iFacet, inside_node, outside_node );

                        double64 velo = this->STENCIL_DATA[eptr->Idx()].FacetNormalVelocity(iFacet) *this->STENCIL_DATA[eptr->Idx()].FacetArea(iFacet);
                        // using upstream weighted first-order fluxes (gravity is not considered)
                        if ( nid == inside_node )
                          div+=velo;
                        else
                          div-=velo;
                    }
               }
              RESULT[nidx] += div;

         }//end Section for Perimeter Nodes
   }// end for cycle for nodes

} // end DivergenceFreeCorrection


/**

SolveNonConservativeTransportEquation1stOrder() advects the non-wetting
phase in a two-phase flow using the user-supplied relative permeability
model. The upstream values of the advected variable are used so that
the scheme is unconditionally stable as long as the CFL condition is
obeyed.

@param time_increment  The input time-increment must satisfy the CFL condition.

@section implementation  Implementation

The divergence of the flow field as a consequence of advecting the
non-wetting phase only, is treated as a source - sink term that is
compensated for in each FV. Transient flow -related sources and sinks
are also considered.
 */
template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::AdvectVariable1stOrder( TwoPhaseModel<dim>& relperm,
                                                                               double64 time_increment )
 {
    fill( RESULT.begin(), RESULT.end(), static_cast<double64>(0.) );

    vector<FV_Parameter>::const_iterator  fvt = this->STENCIL_DATA.begin();

    // 1. Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]

    // fully-upstream weighted version of finite-volume scheme
    if( !this->with_gravitational_forces_ && !this->with_capillary_spreading_){

        for ( typename vector<Element<dim>*>::const_iterator
              eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
          stencil_.AccumulateExplicitTwoPhaseSolution1_Visc( (*fvt), *(*eit), relperm, RESULT );

    }else{

        for ( typename vector<Element<dim>*>::const_iterator
              eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
            stencil_.AccumulateExplicitTwoPhaseSolution1((*fvt), *(*eit), relperm, RESULT, this->with_gravitational_forces_,this->with_capillary_spreading_);

    }

    // compensate for the inflow and the outflow boundaries
    //AssignFractionalFlowBoundaryConditions( relperm );
    AssignGenericFlowBoundaryConditions( relperm );

    // Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]
    //     - a single loop over all nodes is required
    //     - this steps also considers in and outflow of finite volume cells
    //     - (distributed) sources and sinks will be considered in the future
    Compose2PhaseSolution( relperm, time_increment, false );

    // 2. Saving the computed new saturations at the FV centers
    OutputResults( sw_key_, this->adv1_key_, false, false);


    // indicating that something has been done
    //cout <<".";
    //cout.flush();

 } // end AdvectVariable1stOrder


/**

SolveNonConservativeTransportEquation2ndOrder() advects the non-wetting
phase in a two-phase flow using the user-supplied relative permeability
model. The MINMOD slope-limited values of the advected variable are interpolated
to the FV facet integration points using the FE interpolation functions
are used to obtain second-order accuracy in space of this scheme.
The scheme is TVD and robust as long as the CFL condition is
obeyed.

The derivatives of the fractional flow function are computed at each
facet integration point. Thus, there is a dependency on mesh refinement
for the position of the saturation front, i.e. the scheme is more
diffusive where the mesh is coarser.

The divergence of the flow field is treated as a source - sink term that is
compensated for in each FV. Transient flow -related sources and sinks
are also considered.

@param time_increment  the input time-increment must satisfy the
CFL condition.

@section application Application

Use this scheme in models that are not critically dependent on CFL and
when the velocity field needs to be updated in very short time-intervals
so that little is gained from large transport steps that can be performed
using the implicit scheme.
*/
template<size_t dim, template<size_t> class STP>
void TwoPhaseExplicitNodeCenteredFVTransport<dim,STP>::AdvectVariable2ndOrder( TwoPhaseModel<dim>& relperm,
                                                                               double64 time_increment )
 {
    fill( RESULT.begin(), RESULT.end(), static_cast<double64>(0.) );

    // 1. finding smin/smax of previous solution in the neighborhood of each FV
    this->MinMaxAdvectedProperty();
    //this->MinMaxAdvectedPropertyIncludingTheCurrentNode();

    vector<FV_Parameter>::const_iterator  fvt(this->STENCIL_DATA.begin());

    // 2. fully-upstream weighted version of finite-volume scheme
    if ( !this->with_gravitational_forces_ && !this->with_capillary_spreading_ ){

        if(!this->with_lsmgrad_limiter_)
            for ( typename vector<Element<dim>*>::const_iterator
            eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
                stencil_.AccumulateExplicitTwoPhaseSolution2_Visc( this->SMINMAX, (*fvt), *(*eit), relperm, RESULT);
        else{
            this->grad_advprop_limiter_->CalculateGenericNodalGradient();
            this->grad_advprop_limiter_->CalculateSlopeLimiter(this->SMINMAX);
            for ( typename vector<Element<dim>*>::const_iterator
            eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++ )
                stencil_.AccumulateExplicitTwoPhaseSolution2_Visc( this->SMINMAX, (*fvt), *(*eit), relperm, RESULT, this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);
        }

    }else{

        if(this->with_lsmgrad_limiter_){
            /*
            //Element gradient of saturation
            // Alternative
            ElementToNodePropertyVisitor<VectorVariable<dim>,dim> ElmToNode(this->pref_,gradsn_elmt_name_,gradsn_node_name_,this->gref_.Nodes());
            this->gref_.CopyGradientOfProperty_A_To_B(sn_name_,gradsn_elmt_name_);
            this->gref_.Accept(ElmToNode);
            */
            this->grad_advprop_limiter_->CalculateGenericNodalGradient();
            this->grad_advprop_limiter_->CalculateSlopeLimiter(this->SMINMAX);
            for ( typename vector<Element<dim>*>::const_iterator
            eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++)
                stencil_.AccumulateExplicitTwoPhaseSolution2( this->SMINMAX, (*fvt), *(*eit), relperm, RESULT, this->with_gravitational_forces_,this->with_capillary_spreading_,this->mass_center_key_,this->grad_advprop_key_,this->grad_advprop_limiter_key_);

        }else{

            for ( typename vector<Element<dim>*>::const_iterator
            eit=this->gref_.ElementsBegin(); eit!=this->gref_.ElementsEnd(); eit++, fvt++)
                stencil_.AccumulateExplicitTwoPhaseSolution2( this->SMINMAX, (*fvt), *(*eit), relperm, RESULT, this->with_gravitational_forces_, this->with_capillary_spreading_);

        }

    }

    // 3. compensate for the inflow and the outflow boundaries (O.K.)
    //AssignFractionalFlowBoundaryConditions( relperm );
    AssignGenericFlowBoundaryConditions( relperm );

    // 4. Solve S^t+1 = S^t - dt/(phi Vi) * sum_j^faces Aj n . [f vt]
    //     - a single loop over all nodes is required
    //     - this steps also considers in and outflow of finite volume cells
    //     - (distributed) sources and sinks will be considered in the future
    Compose2PhaseSolution( relperm, time_increment, false );

    // 5. Saving the computed new saturations at the FV centers
    OutputResults( sw_key_, this->adv1_key_, false, false );


    // indicating that something has been done
    //cout <<".";
    //cout.flush();

 } // end AdvectVariable2ndOrder





template class TwoPhaseExplicitNodeCenteredFVTransport<1U,ExplicitStencilProcessor>;
template class TwoPhaseExplicitNodeCenteredFVTransport<2U,ExplicitStencilProcessor>;
template class TwoPhaseExplicitNodeCenteredFVTransport<3U,ExplicitStencilProcessor>;
 
} // end namespace csmp
 







   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
   
