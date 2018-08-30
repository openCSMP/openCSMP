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

template<size_t dim, template<size_t> class USER>
TimeStepEvaluator<dim,USER>::TimeStepEvaluator( double64 timestep_reduction_factor, double64 max_time_increment )
 : max_time_increment_(max_time_increment),
   step_size_reduction_factor_(timestep_reduction_factor)
 {
    assert( timestep_reduction_factor > 0. );
    assert( max_time_increment > 0. );
    const double64 one_million_years(1e6 * 365. * 86400.);
    assert( max_time_increment < one_million_years );
 }






/// adjust setting for the solve (default=0.5)
template<size_t dim, template<size_t> class USER>
void TimeStepEvaluator<dim,USER>::StepSizeReductionFactor( double64 factor )
 {
    step_size_reduction_factor_ = factor;
 }



template<size_t dim, template<size_t> class USER>
double64 TimeStepEvaluator<dim,USER>::StepSizeReductionFactor() const
 {
    return step_size_reduction_factor_;
 }



template<size_t dim, template<size_t> class USER>
double64 TimeStepEvaluator<dim,USER>::MaxTimeIncrement() const
 {
    return max_time_increment_;
 }


/** 
    Computes 1) volumetric flow-based time-increment that ascertains that the amount of tracer that leaves the FV
    is less than is stored in its pore volume.
 
    This is robust criterion in the presence of fluid sources and sinks
    
    The method also computes and stores the 2) flux balance (=divergence of flow field) on the current
    FV for later used for correction.
    
    3) the flux concentration products are accumulated into the new concentration variable
*/
template<size_t dim, template<size_t> class USER>
double64 TimeStepEvaluator<dim,USER>::OutFlowLessThanContentIncrement( Node<dim>& n ) const
 {
     double64 flux_balance(0.), outflow(0.), flux_concentration_products(0.);

   if (n.Idx() == 186) {
     std::cerr << "This is the interesting case\n";
   }
   size_t i = 0;
     for (auto fip : n.AllFacetIntegrationPoints()) {
       const double64 ffc = fip.Read(User()->key_ffC);
       const double64 sign = fip.FromInside() ? 1. : -1.;
       // accumulation of volumetric facet flow into flux balance
       const double64 facet_flux = sign * fip.Read( User()->key_ff );
       auto w = fip.IntegrationWeight();
       if ( facet_flux > 0. ) outflow += facet_flux * w;
       flux_balance += facet_flux * w;
       // temporary accumulation of flux-concentration products into the variable 'new concentration'
       flux_concentration_products += sign * ffc * w;
     }

     // 1. (phi * V) / q_out = dt
     double64 time_increment = n.Read( User()->key_FVPV ) / outflow;
     assert( time_increment > 0. );

     // 2. recording the flux balance
     n.Store( User()->key_FB, makeScalar(n.Status(User()->key_FB),flux_balance) );
   
     // 3. recording the flux concentration product balance
     n.Store( User()->key_NC, makeScalar(n.Status(User()->key_NC),flux_concentration_products) );

     return std::min( time_increment, max_time_increment_ ) * step_size_reduction_factor_;
 
 } // OutFlowLessThanContentIncrement




/**
    Special case where the inflow or outflow needs to be established but the FV may be sliced
    by the model boundary.
*/
template<size_t dim, template<size_t> class USER>
double64 TimeStepEvaluator<dim,USER>::OutFlowLessThanContentIncrementBoundary( Node<dim>& n ) const
 {
    const double64 pore_volume = n.Read( User()->key_FVPV );
   
    // 1. if we are at the model boundary we either have in- or outflow; this flow is given by the flux balance
    if ( n.AtBoundary() != NOT ) {
         // 1.0 (Q) checking sources because these will not be picked up in the fluxes correctly if the FV is truncated
         const double64 fluid_source = ( n.Read(User()->key_NQV) > 0. ) ? n.Read(User()->key_NQV) : 0.;
         // 1.1 (phi * V) / q_out = dt
         const double64 flux_balance = std::max( fabs( n.Read( User()->key_FB )), fluid_source );
         if ( fabs(flux_balance) < numeric_limits<double64>::epsilon() ) return max_time_increment_ * step_size_reduction_factor_;
         return std::min( fabs(pore_volume / flux_balance), max_time_increment_ ) * step_size_reduction_factor_;
      }
 
    // 2. for a perimeter FV that is intact, the volumetric outflow needs to be calculated 
     double64 outflow(0.);

     for (auto fip : n.AllFacetIntegrationPoints()) {
       const double64 sign = fip.FromInside() ? 1. : -1.;
       // accumulation of volumetric outflow
       const double64 facet_flux = sign * fip.Read( User()->key_ff );
       if ( facet_flux > 0. ) outflow += facet_flux * fip.IntegrationWeight();
     }

     // (phi * V) / q_out = dt
     if ( outflow < numeric_limits<double64>::epsilon() ) return max_time_increment_ * step_size_reduction_factor_;
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

TODO: implement streamline CFL

*/
template<size_t dim, template<size_t> class USER>
double64  TimeStepEvaluator<dim,USER>::StreamlineCFL( Node<dim>& node ) const
 {
    static DenseMatrix<DM_MIN>  DN;
    vector<double64>            gradPc(dim);
    VectorVariable<dim>         vc;
    double64                    velocity,
                                courant_increment(max_time_increment_);
    const double64              millisecond(1.0e-3);
    const bool                  multiply_with_cell_thickess = (User()->key_THI == csmp::Index()) ? false : true;
    const bool                  unless_has_equal_dimension(dim!=1U);

 /*          // 1. limit imposed by advection
           // -----------------------------
           (*eit)->Read( User()->key_V, vc );
           velocity = vc.Length();


           double64 cell_diameter = (*eit)->LengthInDirection(vc) * (*eit)->Read( User()->key_PHI );
           if ( multiply_with_cell_thickess ) cell_diameter *= (*eit)->Read( User()->key_THI );


           // 4. calculating the CFL criterion from the cell diameter
           // -------------------------------------------------------
           // guarding against degenerate cases
           if ( velocity > numeric_limits<double64>::epsilon() and cell_diameter > numeric_limits<double64>::epsilon() ) {
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
*/
    return courant_increment * step_size_reduction_factor_;

 } // end StreamlineCFL (multiphase case)






/** multiphase flow version

template<size_t dim, template<size_t> class USER>
double64  TimeStepEvaluator<dim,USER>::StreamlineCFL( Node<dim>* const, double64 max_time_increment ) const
 {
    assert( max_time_increment > 0. );
    static DenseMatrix<DM_MIN>  DN;
    vector<double64>            gradPc(dim);
    VectorVariable<dim>         vc;
    double64                    velocity,
                                courant_increment(max_time_increment);
    const double64              millisecond(1.0e-3);
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
           (*eit)->Read( User()->key_V, vc );
           velocity = vc.Length();

           // NB: This may be a too conservative estimate for the implicit scheme but is necessary for the explicit one
           // velocity = std::max( velocity, velocity * relperm.MaxFractionalFlowDerivative() );
           // thus for implicit scheme is better:
           velocity *= relperm.dfds();

           double64 cell_diameter = (*eit)->LengthInDirection(vc) * (*eit)->Read( User()->key_PHI );
           if ( multiply_with_cell_thickess ) cell_diameter *= (*eit)->Read( User()->key_THI );

           // 2. additional buoyancy-related flow
           // -------------------------------------
           if ( with_gravitational_forces_ )
             velocity += fabs(relperm.GravityMultiplier_dGds());

           // 3. limit due to potential capillary spreading
           // ---------------------------------------------
           if ( with_capillary_spreading_ ) {
               const double64 k_lambda_overbar(relperm.Permeability() * relperm.G());
               // computing the capillary pressure gradient
               if ( k_lambda_overbar > numeric_limits<double64>::epsilon() )
                 {
                     fill( gradPc.begin(), gradPc.end(), 0. );
                     (*(*eit)).dN_AtBaryCenter( DN );
                     for ( size_t j=0U; j<(*eit)->Nodes(); j++ ) {
                         const double64 sn = (*eit)->N(j)->Read( User()->key_C );
                         relperm.SaturationWettingPhase( 1. - sn );
                         relperm.EffectiveSaturation();
                         const double64 pc = relperm.pc_Phase();
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

template class TimeStepEvaluator<1U,ExplicitTransport>;
template class TimeStepEvaluator<2U,ExplicitTransport>;
template class TimeStepEvaluator<3U,ExplicitTransport>;

template class TimeStepEvaluator<1U,ImplicitTransport>;
template class TimeStepEvaluator<2U,ImplicitTransport>;
template class TimeStepEvaluator<3U,ImplicitTransport>;

} // end csmp


