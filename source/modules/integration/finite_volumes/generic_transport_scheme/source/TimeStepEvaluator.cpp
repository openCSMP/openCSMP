//
//  TimeStepEvaluator.cpp
//  CSMP_ReservoirSimulator
//
//  Created by Stephan Matthai on 2/21/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "TimeStepEvaluator.h"
#include "ExplicitTransport.h"
#include "ImplicitTransport.h"
#include "Exception.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "TwoPhaseModel.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class USER>
TimeStepEvaluator<dim,USER>::TimeStepEvaluator( double timestep_reduction_factor, double max_time_increment )
 : max_time_increment_(max_time_increment),
   step_size_reduction_factor_(timestep_reduction_factor)
 {
    assert( timestep_reduction_factor > 0. );
    assert( max_time_increment > 0. );
    const double one_million_years(1e6 * 365. * 86400.);
    assert( max_time_increment < one_million_years );
 }






/// adjust setting for the solve (default=0.5)
template<uint32_t dim, template<uint32_t> class USER>
void TimeStepEvaluator<dim,USER>::StepSizeReductionFactor( double factor )
 {
    step_size_reduction_factor_ = factor;
 }



template<uint32_t dim, template<uint32_t> class USER>
double TimeStepEvaluator<dim,USER>::StepSizeReductionFactor() const
 {
    return step_size_reduction_factor_;
 }



template<uint32_t dim, template<uint32_t> class USER>
double TimeStepEvaluator<dim,USER>::MaxTimeIncrement() const
 {
    return max_time_increment_;
 }





/** 
    Computes 1) volumetric flow-based time-increment that ascertains that the amount of tracer that leaves the FV
    is less than is stored in its pore volume.
 
    This is robust criterion in the presence of fluid sources and sinks.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TimeStepEvaluator<dim,USER>::OutFlowLessThanContentIncrement( Node<dim>* const nptr, double outflow ) const
 {
     // 1. (phi * V) / q_out = dt
     double time_increment = nptr->Read( User()->Notation.key_FVPV ) / outflow;
     assert( time_increment > 0. );

     return std::min( time_increment, max_time_increment_ ) * step_size_reduction_factor_;
 
 } // OutFlowLessThanContentIncrement




/**
    Special case where inflow or outflow needs to be established, differently because the FV is sliced
    by the model boundary.
 
    2 boundary cases are considered:
 
      1. No-flow boundary where flux balance is known because no in- or outflow should occur
 
      2. Dirichlet in- or outflow boundary where balance has to be guessed and compensation
         is necessary, else there is a saturation buildup at out-flow boundaries.
*/
template<uint32_t dim, template<uint32_t> class USER>
double TimeStepEvaluator<dim,USER>::OutFlowLessThanContentIncrementBoundary( const Node<dim>* const nptr, double outflow ) const
  {
     assert( nptr != NULL );
 
     // 0. COMPUTING THE FLUX BALANCE
     const double pore_volume  = nptr->Read( User()->Notation.key_FVPV );
     double       flux_balance = nptr->Read( User()->Notation.key_FB );

    // 1. FINITE VOLUMES TRUNCATED BY MODEL BOUNDARY
    //    if we are at the model boundary we either have in- or outflow; this flow is given by the flux balance.
    //    however, if we are at a no-flow boundary of the model, there should be no inflow or outflow through the missing facets
    //    and any flux balance will therefore be correctly computed.
    if ( nptr->AtBoundary() != NOT ) {
         if ( nptr->Status( User()->Notation.key_PF ) == DIRICH ) {
               // if we are at an outflow boundary, no outflow was recorded as all FV facet normals point into the model domaim
               // in this the flux-balance is equivalent to the outflow
               if ( flux_balance < 0. ) outflow = fabs(flux_balance);
           }
         // for no-flow boundary cases, the outflow and the flux balances are recorded correctly
      }
 
    // 2. INTACT PERIMETER FINITE VOLUMES
    //    if the perimeter of the FV at the boundary is intact because it lies inside of the model,
    //    the volumetric outflow is calculated and used to limit the time increment

     // (phi * V) / q_out = dt
     return std::min( pore_volume / outflow, max_time_increment_ ) * step_size_reduction_factor_;
 
 } // OutFlowLessThanContentIncrementBoundary







/** CFL for the strictly hyperbolic case in multiphase flow

Computes the Courant time increment (CFL citerion) taking into account
viscous, gravitational and capillary fluid displacements using the 
contraints from the provided relative permeability model. The CFL
criterion is calculated using the element length in the direction
of the flow but not account for the deviation from this vector of
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
template<uint32_t dim, template<uint32_t> class USER>
double  TimeStepEvaluator<dim,USER>::StreamlineCFL( const Element<dim>* const eit ) const
 {
    vector<double>            gradPc(dim);
    VectorVariable<dim>         vc;
    double                    velocity, courant_increment(max_time_increment_);
    const double              millisecond(1.0e-3);
    const bool                  multiply_with_cell_thickess = (User()->Notation.key_THI == csmp::Index()) ? false : true;

   // 1. limit imposed by advection
   // -----------------------------
   eit->Read( User()->Notation.key_V, vc );
   velocity = vc.Length();


   double cell_diameter = eit->LengthInDirection(vc) * eit->Read( User()->Notation.key_PHI );
   if ( multiply_with_cell_thickess ) cell_diameter *= eit->Read( User()->Notation.key_THI );


   // 4. calculating the CFL criterion from the cell diameter
   // -------------------------------------------------------
   // guarding against degenerate cases
   if ( velocity > numeric_limits<double>::epsilon() and cell_diameter > numeric_limits<double>::epsilon() ) {
        courant_increment = std::min( courant_increment, cell_diameter / velocity );
     }

    if ( courant_increment >= max_time_increment_ ) {
         cout <<"\nINFO, TimeStepEvaluator<dim>::StreamlineCFL (2-phase flow):";
         cout <<" calculated courant increment is larger than maximum permitted increment, there may be no flow in the model domain !";
         cout <<"\n\tCFL is set to "<< fixed << setprecision(0) << max_time_increment_ << endl;
         return max_time_increment_;
      }
    if ( courant_increment <= millisecond ) {
         throw csmp::Exception( WARNING, "TimeStepEvaluator<dim>::StreamlineCFL (2-phase flow):",
                                         "courant increment is smaller than a millisecond. Check your boundary conditions." );
      }
    else {
         cout <<"\nTimeStepEvaluator<dim>::StreamlineCFL (2-phase flow): ";
         cout << fixed << setprecision(0) << courant_increment <<" secs.\n";
      }

    return courant_increment * step_size_reduction_factor_;
   
 } // end StreamlineCFL (multiphase case)






/** multiphase flow version

*/
/* TODO: port this to new relperm framework
template<uint32_t dim, template<uint32_t> class USER>
double  TimeStepEvaluator<dim,USER>::StreamlineCFL( Node<dim>* const, double max_time_increment ) const
 {
    assert( max_time_increment > 0. );
    static DenseMatrix<DM_MIN>  DN;
    vector<double>            gradPc(dim);
    VectorVariable<dim>         vc;
    double                    velocity,
                                courant_increment(max_time_increment);
    const double              millisecond(1.0e-3);
    const bool                  multiply_with_cell_thickess = (User()->key_THI == csmp::Index()) ? false : true;
    const bool                  unless_has_equal_dimension(dim!=1U);

    for ( typename vector<Element<dim>*>::const_iterator 
          eit=gref_.ElementsBegin(); eit!=gref_.ElementsEnd(); ++eit )
      if ( !((*eit)->FE()->IsLineElement() && unless_has_equal_dimension) )
        {
           // 0. relative permeability model computed at element barycenter
           // -------------------------------------------------------------
           relperm.Initialize( *(*eit) );
           relperm.InitializeForBaryCenter( *(*eit) );
           relperm.EffectiveSaturation();

           // 1. limit imposed by advection
           // -----------------------------
           (*eit)->Read( User()->vD_key, vc );
           velocity = vc.Length();

           // NB: This may be a too conservative estimate for the implicit scheme but is necessary for the explicit one
           // velocity = std::max( velocity, velocity * relperm.MaxFractionalFlowDerivative() );
           // thus for implicit scheme is better:
           velocity *= relperm.dfds();

           double cell_diameter = (*eit)->LengthInDirection(vc) * (*eit)->Read( User()->phi_key );
           if ( multiply_with_cell_thickess ) cell_diameter *= (*eit)->Read( User()->thi_key );

           // 2. additional buoyancy-related flow
           // -------------------------------------
           if ( with_gravitational_forces_ )
             velocity += fabs(relperm.GravityMultiplier_dGds());

           // 3. limit due to potential capillary spreading
           // ---------------------------------------------
           if ( with_capillary_spreading_ ) {
               const double k_lambda_overbar(relperm.Permeability() * relperm.G());
               // computing the capillary pressure gradient
               if ( k_lambda_overbar > numeric_limits<double>::epsilon() )
                 {
                     fill( gradPc.begin(), gradPc.end(), 0. );
                     (*(*eit)).dN_AtBaryCenter( DN );
                     for ( size_t j=0U; j<(*eit)->Nodes(); j++ ) {
                         const double sn = (*eit)->N(j)->Read( User()->C0_key );
                         relperm.SaturationWettingPhase( 1. - sn );
                         relperm.EffectiveSaturation();
                         const double pc = relperm.pc_Phase();
                         for ( size_t k=0U; k<dim; k++ ) gradPc[k] += DN(k,j) * pc;
                     }
                     // getting the maximum capillary flux (G= lambda overbar)
                     double  magnitude_grad_pc(gradPc[0]); // 1D
                     if      ( dim == 3U ) magnitude_grad_pc = sqrt(gradPc[0]*gradPc[0]+gradPc[1]*gradPc[1]+gradPc[2]*gradPc[2]);
                     else if ( dim == 2U ) magnitude_grad_pc = sqrt(gradPc[0]*gradPc[0]+gradPc[1]*gradPc[1]);

                     // use data from barycenter: O.K. as long as grad_pc does not increase during iterations
                     velocity += fabs( magnitude_grad_pc * k_lambda_overbar );
                 }
           }

           // 4. calculating the CFL criterion from the cell diameter
           // -------------------------------------------------------
           // guarding against degenerate cases
           if ( velocity > numeric_limits<double>::epsilon() and cell_diameter > numeric_limits<double>::epsilon() ) {
                courant_increment = std::min( courant_increment, cell_diameter / velocity );
             }
        }

    if ( courant_increment >= max_time_increment ) { 
         cout <<"\nINFO, TimeStepEvaluator<dim>::StreamlineCFL (2-phase flow):";
         cout <<" calculated courant increment is larger than maximum permitted increment, there may be no flow in the model domain !";
         cout <<"\n\tCFL is set to "<< fixed << setprecision(0) << max_time_increment << endl;
         return max_time_increment;
      }
    if ( courant_increment <= millisecond ) {
         throw csmp::Exception( WARNING, "TimeStepEvaluator<dim>::StreamlineCFL (2-phase flow):",
                                         "courant increment is smaller than a millisecond. Check your boundary conditions." );
      }
    else {
         cout <<"\nTimeStepEvaluator<dim>::StreamlineCFL (2-phase flow): ";
         cout << fixed << setprecision(0) << courant_increment <<" secs.\n";
      }
      
    return courant_increment;
   
 } // end StreamlineCFL (multiphase case)
*/

template class TimeStepEvaluator<2U,ExplicitTransport>;
template class TimeStepEvaluator<2U,ImplicitTransport>;

template class TimeStepEvaluator<3U,ExplicitTransport>;
template class TimeStepEvaluator<3U,ImplicitTransport>;

} // end csmp


