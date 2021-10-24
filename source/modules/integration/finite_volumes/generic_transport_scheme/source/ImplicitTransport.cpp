#include "ImplicitTransport.h"
#include "Region.h"
#include "Model.h"
#include "finiteVolumeAuxiliaryFunctions.h"
#include "CSMP_mathUtilities.h"
#include "VTK_Interface.h"

using namespace std;
using namespace csmp::variables;

namespace csmp {

template<size_t dim>
ImplicitTransport<dim>::ImplicitTransport( Model<dim>& model, const std::string& target_region, 
                                           const GoverningEquation<dim>& eqn, 
                                           bool second_order_in_space )
  : Integrator<dim,ImplicitTransport>(model.Mesh().Nodes(),
                                             model.Database().LowerLimitOf("concentration"),
                                             model.Database().UpperLimitOf("concentration")),
    GoverningEquation<dim>(eqn),
    subdomain_(model.Region(target_region)),
    linalg_sys_(subdomain_.Nodes()),
    second_order_(second_order_in_space)
 {
    const bool initialise_flux(false);
    initializeFiniteVolumeProperties( model, subdomain_, initialise_flux );
    
    // 0. model-wide initialisation: results will be accumulated into this variable
 }
  
  
/**
     stub to FV auxiliary function
*/
template<size_t dim>
void ImplicitTransport<dim>::UpdateFluxesAndFluxBalances()
 {
    this->VolumetricFlowAndTransportVariableFluxBalances( ComputationDomain(), halo_cells_ );

 } // end UpdateFluxesAndFluxBalances
 
  
  
/**
    Computation of time increment, flux balance, and temporary new concentration.
*/
template<size_t dim>
double64 ImplicitTransport<dim>::TimeIncrementAndFluxBalance( double64 max_time_increment ) const
 {
     double64 dt_min(max_time_increment);
 
     // 1. processing interior and FVs for which all facet fluxes have been initialised
     const typename vector<Node<dim>*>::const_iterator interior_nodes_end(subdomain_.PerimeterNodesBegin());
     for ( typename vector<Node<dim>*>::const_iterator nit=subdomain_.InteriorNodesBegin(); nit!=interior_nodes_end; ++nit )
       {
          assert( (*nit)->AtBoundary() == NOT );
          const double64 out_flow = this->FluxBalanceAndOutFlow( (*nit) );
          // computes time-increment, flux balance, and flux-concentration product balance
          const double64 time_increment = this->OutFlowLessThanContentIncrement( *nit, out_flow );
          dt_min = std::min( dt_min, time_increment );
       }

     // 2. collecting time-stepping constraints from FVs on region perimeter
     const typename vector<Node<dim>*>::const_iterator nodes_end(ComputationDomain().PerimeterNodesEnd());
     for ( typename vector<Node<dim>*>::const_iterator
           nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
          const double64 out_flow = this->FluxBalanceAndOutFlow( (*nit) );
          // boundary fluxes must be part of the time-increment calculation
          dt_min = std::min( dt_min, this->OutFlowLessThanContentIncrementBoundary( *nit, out_flow ) );
       }
   
    return dt_min;
 }





/**
    Computation of time increment (default limited to 1 year).
*/
template<size_t dim>
double64 ImplicitTransport<dim>::TimeIncrement() const
 {
    return TimeIncrementAndFluxBalance( 356. * 86400. );
 }


  
  
  
/**
    Reports the volumetric flow into the computational region.
*/
template<size_t dim>
double64 ImplicitTransport<dim>::IncomingVolumetricFlow() const
 {
     double64 inflow(0.);
 
     const typename vector<Node<dim>*>::const_iterator nodes_end(subdomain_.NodesEnd());
     for ( typename vector<Node<dim>*>::const_iterator
           nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
           if ( HasHaloStencils() ) {
                // away from model boundaries, the inflow is computed correctly and can be used
                if ( (*nit)->AtBoundary() == NOT ) inflow += this->InFlow( (*nit) );           
             }
           // at model boundaries, FV's are truncated and only the their inside part exist; inflows are outflows (+) for these
           if ( (*nit)->AtBoundary() != NOT ) {
               // to get accurate predictions and without wasting computational cost, no-flow boundaries are excluded
               const VARIABLE_FLAG pf_flag = (*nit)->Status( this->Notation.key_PF );
               if ( pf_flag == DIRICH || pf_flag == NEUMANN ) {
                    const double64 flow = this->VolumetricFlowBalance( (*nit) );
                    if ( flow > 0. )
                      inflow += flow;
                 }
             }
       }
 
    // incoming flux should be positive
    return inflow;

 } // end IncomingVolumetricFlow
  
  
/**
    Reports the volumetric flow outside of the computational region.
 
    @todo SKM check whether the variable 'outflow' can be used for this, saving some computations.
*/
template<size_t dim>
double64 ImplicitTransport<dim>::OutgoingVolumetricFlow() const
 {
     double64 outflow(0.);
 
     const typename vector<Node<dim>*>::const_iterator nodes_end(subdomain_.NodesEnd());
     for ( typename vector<Node<dim>*>::const_iterator
           nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
           if ( HasHaloStencils() ) {
                // away from model boundaries, the inflow is computed correctly and can be used
                if ( (*nit)->AtBoundary() == NOT ) outflow += this->OutFlow( (*nit) );           
             }
           // at model boundaries, FV's are truncated and only the their inside part exist; inflows are outflows (+) for these
           if ( (*nit)->AtBoundary() != NOT ) {
               // to get accurate predictions and without wasting computational cost, no-flow boundaries are excluded
               const VARIABLE_FLAG pf_flag = (*nit)->Status( this->Notation.key_PF );
               if ( pf_flag == DIRICH || pf_flag == NEUMANN ) {
                    const double64 flow = this->VolumetricFlowBalance( (*nit) );
                    if ( flow < 0. )
                      outflow += flow;
                 }
             }
       }
 
    // since the accumulated outflow is negative, it needs to be inverted 
    return -outflow;
 
 } // end OutgoingVolumetricFlow






/**
    1. Computation of (velocity and) facet fluxes
    
    2. Evaluation of time increment, flux balances and flux balance concentration products
 
    3. Assembly of solution

    4. Assignment of boundary conditions (treatment of truncated FVs at the boundary)
 
    5. tranfer of results into 'concentration', vacating 'new concentration' for 
       next assembly.
*/
template<size_t dim>
void ImplicitTransport<dim>::AdvectVariable( double64 time_interval )
 {
   // 1. evaluation of time increment
   double64 time_increment = TimeIncrementAndFluxBalance( this->MaxTimeIncrement() );
   
   const double64 one(1.);
   cout <<"\nImplicitTransport<"<< fixed << setprecision(0) << dim <<">::AdvectVariable:";
   cout <<"\n\tTime interval         = "<< time_interval;
   cout <<"\n\tScaled time increment = "<< time_increment;
   cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);
   
   cout <<"\n\n\nImplicitTransport::AdvectVariable: simulating "<< time_interval/60. <<" minutes of transport...\n";
   size_t   substep(1);
   double64 time(0.);
   
   // 2. setup matrix and vectors; initialise indices in subdomain
   this->SetUp(); 

   // 3. time stepping loop
   while ( time < time_interval ) 
     {
        cout <<"\n\n\tadvection (sub)step: "<< substep;
        if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
       
        this->Accumulate( time_increment );
        
        // assignment of initial and essential conditions, SAMG solution, postprocessing, solution verification and return of results
        this->IntegrateOver( time_increment );
       
        time += time_increment;
       
        time_increment = TimeIncrementAndFluxBalance( this->MaxTimeIncrement() );
       
        substep++;
     }
   
   cout <<"\nImplicitTransport::AdvectVariable: 'transport completed.\n";   

 } // end AdvectVariable



//  template class ImplicitTransport<1U>;
//  template class ImplicitTransport<2U>;
  template class ImplicitTransport<3U>;

} // end csmp

