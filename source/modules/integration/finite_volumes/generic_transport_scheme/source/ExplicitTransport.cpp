//
//  ExplicitTransport.cpp
//
//  Created by Stephan Matthai on 2/21/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "ExplicitTransport.h"
#include "Region.h"
#include "Model.h"
#include "finiteVolumeAuxiliaryFunctions.h"

#include "VTK_Interface.h"

using namespace std;
using namespace csmp::variables;

namespace csmp {

/**
    Custom constructor using an existing model and the default values of the policy classes.
    Halo stencils are created when the computaitonal region is smaller than the model.
    Finite volume and sector volumes are initialised.
    Facet fluxes and flux balances are initialised if the 'velocity' is initialised.
*/
template<size_t dim>
ExplicitTransport<dim>::ExplicitTransport( Model<dim>& m, const char* target_region )
  : VariableSet_TracerTransfer(m.Database()),
    subdomain_(m.Region(target_region)),
    key_acc_(m.CreateProperty("accumulation","node")),
    key_out_(m.CreateProperty("outflow","node"))
  {
    m.InstantiateFiniteVolumes();
  
    // 0. diagnostics: has velocity been initialised, if so, it is used in initialisation.
    double64 vmin, vmax;
    subdomain_.MinMaxOf( "velocity", vmin, vmax );
    const bool initialise_flux = ( !isnan(vmin) && !isnan(vmax) ) ? true : false;
    
    // NOTE: this method does not consider the halo stencils in the flux calculations
    initializeFiniteVolumeProperties( m, m.Region(target_region), initialise_flux );

    // 1. model-wide initialisation: results will be accumulated into this variable
    m.Region("Model").InputPropertyValue( "accumulation", makeScalar(PLAIN,0.), COMPLETE );
   
    // 2. retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->key_C0), lower_limit_, upper_limit_ );
  
    // 3. recording potential halo elements that will need to be initialised
    CollectHaloStencils();
  
    // 4. initialising facet fluxes and flux balances, but only if transport has reasonable values
    if ( initialise_flux ) TransportVariableFluxBalances();
 }
  
  

/**
    Initialises the halo stencils vector and sorts it so that it can be searched.
 
    @attention memory management is left entirely to the vector.
*/
template<size_t dim>
size_t ExplicitTransport<dim>::CollectHaloStencils()
 {
    const typename vector<Node<dim>*>::const_iterator nds_end(subdomain_.NodesEnd());
 
    // checking whether halo elements exist, and estimating how many there may be
    size_t potential_halo_elements(0U);
    for ( typename vector<Node<dim>*>::const_iterator
          nit=subdomain_.NodesBegin(); nit!=nds_end; ++nit )
      // only where perimeter nodes are not at the model boundary, halo elements may exist
      if ( (*nit)->AtBoundary() == NOT )
        potential_halo_elements++;
 
    if ( potential_halo_elements == 0U ) return 0U;
 
    // collecting the halo elements into the vector
    halo_elmts_.reserve( potential_halo_elements * 3 );
    for ( typename vector<Node<dim>*>::const_iterator
          nit=subdomain_.NodesBegin(); nit!=nds_end; ++nit ) {
         // looping over the parent elements of perimeter nodes that are not situated
         // at the model boundary,  adding those parent elements to the halo
         // which are not contained in the subdomain
         if ( (*nit)->AtBoundary() == NOT )
           for ( size_t i=0U; i<(*nit)->Parents(); ++i )
             if ( !subdomain_.Contains( (*nit)->Parent(i) ) )
               halo_elmts_.push_back( (*nit)->Parent(i) );
      }
 
    sort( halo_elmts_.begin(), halo_elmts_.end() );
 
    return halo_elmts_.size();
 
 } // end CollectHaloStencils



/**
   initialises transport class for current pressure/velocity/transport variable field.
*/
template<size_t dim>
void ExplicitTransport<dim>::UpdateFluxesAndFluxBalances()
 {
    VolumetricFlowAndTransportVariableFluxBalances();
 }



/**
    Volumetric facet flow and (chemical) flux calculation (element by element), for all elements in the domain
    
    Facet fluxes in the FV sectors outside the domain are considered.
    by Advective_O1_FluxesBoundary().
    This method also writes the correct volumetric flow balances onto the nodes on the perimeter.
*/
template<size_t dim>
void ExplicitTransport<dim>::VolumetricFlowAndTransportVariableFluxBalances()
 {
     // 1. setting 'flux balance' and 'accumulation' variables to be accumulated on the nodes=FVs to zero
     const typename vector<Node<dim>*>::iterator nodes_end(subdomain_.NodesEnd());
     for ( typename vector<Node<dim>*>::iterator
           nit=subdomain_.NodesBegin(); nit!=nodes_end; ++nit ) {
           (*nit)->Store( key_FB, makeScalar(ANY,0.) );
           (*nit)->Store( key_acc_, makeScalar(ANY,0.) );
           (*nit)->Store( key_out_, makeScalar(ANY,0.) );
       }
     // including the halo elements
     if ( HasHaloElements() )
       for ( typename vector<Element<dim>*>::iterator it=halo_elmts_.begin(); it!=halo_elmts_.end(); ++it ) {
            const size_t nodes((*it)->Nodes());
            for ( size_t i=0U; i<nodes; ++i ) {
                 (*it)->N(i)->Store( key_FB, makeScalar(ANY,0.) );
                 (*it)->N(i)->Store( key_acc_, makeScalar(ANY,0.) );
                 (*it)->N(i)->Store( key_out_, makeScalar(ANY,0.) );
              }
         }

     // 2. element-by-element processing of facet fluxes and flux balances
     const typename vector<Element<dim>*>::iterator elements_end(subdomain_.ElementsEnd());
     for ( typename vector<Element<dim>*>::iterator
           eit=subdomain_.ElementsBegin(); eit!=elements_end; ++eit )
       {
          // 1.1 computation of transport velocity from fluid pressure gradient
       
          // 1.2 computes facet fluxes (including upstream concentration- and flux balances
          //     on all elements of the domain, but not the halo stencils
          this->Advective_O1_FluxesAndBalances( (*eit) );
       }
 
     // 3. processing potential halo stencils
     if ( HasHaloElements() ) {
         const typename vector<Element<dim>*>::iterator elements_end(halo_elmts_.end());
         for ( typename vector<Element<dim>*>::iterator
               eit=halo_elmts_.begin(); eit!=elements_end; ++eit )
           {
              // 1.1 computation of transport velocity from fluid pressure gradient
           
              // 1.2 computes facet fluxes (including upstream concentration- and flux balances
              //     on all elements of the domain, but not the halo stencils
              this->Advective_O1_FluxesAndBalances( (*eit) );
           }
       }

     // 4. processing fluxes through the FVs on the regions perimeter computing outside facet fluxes as necessary
     for ( typename vector<Node<dim>*>::iterator
           nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
          // computing volumetric flux balances and concentration-facet flux product balances.
          // At sliced boundaries 3-typed of conditions are applied:
          //   1) prescribed concentration value at inflow boundaries,
          //   2) prescribed flux (has consequence only where there is inflow),
          //   3) free outflow found from the flux balance and the FV cell's current concentration
          this->FluxBalancesAtBoundary( (*nit) );
       }


 } // end VolumetricFlowAndTransportVariableFluxBalances

 
 

// TESTING: checking the flux balances
//for ( typename vector<Node<dim>*>::const_iterator
//      nit=subdomain_.NodesBegin(); nit!=subdomain_.PerimeterNodesBegin(); ++nit )
//  if ( fabs((*nit)->Read(key_FB)) > numeric_limits<double64>::epsilon()*3. )
//    cerr << ios::scientific << (*nit)->Read(key_FB) <<" ";


 
 
 
 

/**
    Computes FV balances of products of the transported variable with the facet fluxes.
    The volumetric flux balance is not updated, aussuming that the volumetric facet flows
    did not change.
*/
template<size_t dim>
void ExplicitTransport<dim>::TransportVariableFluxBalances()
 {
     // 1. setting 'accumulation' variable to zero
     const typename vector<Node<dim>*>::iterator nodes_end(subdomain_.NodesEnd());
     for ( typename vector<Node<dim>*>::iterator
           nit=subdomain_.NodesBegin(); nit!=nodes_end; ++nit ) {
          (*nit)->Store( key_acc_, makeScalar(ANY,0.) );
          (*nit)->Store( key_out_, makeScalar(ANY,0.) );
       }

     // including the halo elements
     if ( HasHaloElements() )
       for ( typename vector<Element<dim>*>::iterator it=halo_elmts_.begin(); it!=halo_elmts_.end(); ++it ) {
            const size_t nodes((*it)->Nodes());
            for ( size_t i=0U; i<nodes; ++i ) {
                 (*it)->N(i)->Store( key_acc_, makeScalar(ANY,0.) );
                 (*it)->N(i)->Store( key_out_, makeScalar(ANY,0.) );
              }
         }
         
     // 2. computing volumetric flow - transport variable products, storing them in variable 'accumulation'
     const typename vector<Element<dim>*>::iterator elements_end(subdomain_.ElementsEnd());
     for ( typename vector<Element<dim>*>::iterator
           eit=subdomain_.ElementsBegin(); eit!=elements_end; ++eit )
       {
          this->TransportVariableFluxesAndBalances( (*eit) );
       }
 
     // 3. same process applied to the halo stencils
     if ( HasHaloElements() ) {
         const typename vector<Element<dim>*>::iterator elements_end(halo_elmts_.end());
         for ( typename vector<Element<dim>*>::iterator
               eit=halo_elmts_.begin(); eit!=elements_end; ++eit )
           {
              this->TransportVariableFluxesAndBalances( (*eit) );
           }
       }

     // 4. processing fluxes through the FVs on the regions perimeter computing outside facet fluxes as necessary
     for ( typename vector<Node<dim>*>::iterator
           nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
          // computing volumetric flux balances and concentration-facet flux product balances.
          // At sliced boundaries 3-types of conditions are considered:
          //   1) prescribed concentration value at inflow boundaries,
          //   2) TODO: prescribed flux (has consequence only where there is inflow),
          //   3) free outflow found from the flux balance and the FV cell's current concentration
          this->FluxBalancesAtBoundary( (*nit) );
       }

 } // end TransportVariableFluxBalances






/**
    Computation of time increment, flux balance, and temporary new concentration.
 
    @attention default is the most stringent time increment: outflux < PV.
*/
template<size_t dim>
double64 ExplicitTransport<dim>::TimeIncrement_CFL_Outflow( double64 max_time_increment )
 {
     double64 dt_min(max_time_increment);
 
     // processing interior and perimeter FVs for which facet fluxes have been initialised
     const typename vector<Node<dim>*>::const_iterator nodes_end(subdomain_.NodesEnd());
     for ( typename vector<Node<dim>*>::const_iterator
           nit=subdomain_.NodesBegin(); nit!=nodes_end; ++nit )
       {
          // using pre-computed facet fluxes
          // this recomputes facet fluxes, avoid this
//          const double64 out_flow = (*nit)->Read( key_out_ ); // this->OutFlow( (*nit) ); TODO: breaks the code
          const double64 out_flow = this->OutFlow( (*nit) );
          // computes time-increment, flux balance, and flux-concentration product balance
          const double64 time_increment = this->OutFlowLessThanContentIncrement( (*nit), out_flow );
          dt_min = std::min( dt_min, time_increment );
       }
 
    return dt_min;
 }

// tested OK: cerr <<"\noutflow: node "<< (*nit)->Idx() <<": "<< out_flow <<" vs. "<< (*nit)->Read( key_out_ );




/**
    Computation of time increment (default limited to 1 year).
 
    @attention Assumes that the facet flows and tracer fluxes are initialised.
*/
template<size_t dim>
double64 ExplicitTransport<dim>::TimeIncrement()
 {
    const double64 default_max_time_increment( 356. * 86400. ); // 1 year
    return TimeIncrement_CFL_Outflow( default_max_time_increment );
 }








/**  AssembleSolution()

      C1^t+1 = C0^t - dt/(phi Vi) * (sum_j^faces Aj n . [C0 vt] + sources/sinks)
 
    - a single loop over all nodes is required
    - this steps also considers in and outflow of finite volume cells
    - (distributed) sources and sinks will be considered in the future
 
    @attention diffusion is not considered because scheme is highly diffusive anyway
    
    @todo build second-order in space version, considering diffusion

         // TODO: account for finite-element 'fluid volume source' terms by distributing contributions to the nodes (divide by # nodes)
         // TODO: account for FV source terms ('nodal fluid volume source' directly
         // TODO: make this more accurate use fractional step method where source is accounted for at 2 time levels using dt/2 and C0 and C1
*/
template<size_t dim>
void ExplicitTransport<dim>::Assemble1stOrderSolution( double64 delta_t,
                                                       bool enforce_divergence_free_vt_field )
 {
    const typename vector<Node<dim>*>::const_iterator  nodes_end(subdomain_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=subdomain_.NodesBegin(); nit!=nodes_end; ++nit )
      {
         // 1. reading the fluxes accumulated in 'accumulation'
         //    - flux upwind CO products
         //    - boundary conditions applied to perimeter volumes and truncated model boundary volumes
         double64 accumulation((*nit)->Read(this->key_acc_));
      
         // 2. correcting for div vD non-zero 'flux balance' term except at Dirichlet boundaries
         //    or where nodal sources were present
         if ( enforce_divergence_free_vt_field ) {
              // taking care of 'nodal fluid volume source' or sink terms
              const double64 source((*nit)->Read(this->key_NQV));
              if ( fabs(source) < numeric_limits<double64>::epsilon() &&
                   (*nit)->Status( key_PF ) != DIRICH &&
                   (*nit)->AtBoundary() != NOT )
                // multiply because the error relates to the volumetric flow part of  C0 * volumetric flow products
                accumulation += (*nit)->Read(this->key_C0) * -(*nit)->Read(this->key_FB);
           }

         // compute time-increment - pore volume product
         const double64 dt_divided_by_PV  = (delta_t/(*nit)->Read(this->key_FVPV));

         // 3. starting with the sum of facet flux-concentration products temporarily stored in 'accumulation'
          //   solution is assembled by subtracting flux time-interval products from value at previous time level
         double64 solution = (*nit)->Read(this->key_C0) - dt_divided_by_PV * accumulation;
      
         // 4. TODO: element 'fluid volume source' terms have to be distributed over the nodes of the finite element

         // 5. storing the new concentration
         if ( (*nit)->Status( this->key_C0 ) != DIRICH )
           (*nit)->Store( this->key_C0, makeScalar( (*nit)->Status(this->key_C0), solution ) );
     }

} // end AssembleSolution


/*
 
    // flux compensation at outflow boundaries
    for ( typename vector<Node<dim>*>::const_iterator nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
      if ( (*nit)->AtBoundary() != NOT and
           (*nit)->Status(key_PF) == DIRICH and
           (*nit)->Read(key_FB) < 0. )
        {
           // 1. computing the compensation
           double64 corrected_accumulation = (*nit)->Read( key_C0 ) - (*nit)->Read( key_acc_ );
           // 2. storing the new concentration
           (*nit)->Store( this->key_C0, makeScalar( (*nit)->Status(this->key_C0), corrected_accumulation ) );
        }
*/



/**
    For all nodes on the perimeter of the target region, this method
    loops over the finite-element sectors of the associated finite volume
    integrating facet fluxes based on total velocity.
    
    This flux balance is added to the node variable 'new saturation oil'
*/
template<size_t dim>
void ExplicitTransport<dim>::AdjustResultsAssumingDivergenceFreeVelocityField( double64 time_interval )
 {
    const typename vector<Node<dim>*>::iterator  nodes_end(subdomain_.NodesEnd());
    for ( typename vector<Node<dim>*>::iterator nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
      {
          double64 div(0.);
          // for each finite volume, f is evaluated on a sector by sector basis
          const size_t parents((*nit)->Parents());
          for ( size_t t=0U; t<parents; t++ ) {
               Element<dim>* const eptr((*nit)->Parent(t));
               const size_t nid((*nit)->ParentNodeNumber(t));

               // for all FACETS per SECTOR surrounding the finite volume at the boundary
               for ( size_t i=0U; i<eptr->FV()->FacetsPerSector(nid); i++ ) {
                    size_t iFacet( eptr->FV()->FacetSurroundingSector(nid,i) );
                    double64 velo = eptr->ProjectionOnFacetNormal( iFacet, this->key_V );
                    if ( nid == eptr->FV()->InsideNode(iFacet) )div += velo;
                    else div -= velo;
                }
           }
          ScalarVariable result( makeScalar( (*nit)->Status(this->key_C1), (*nit)->Read(this->key_C1) ) );
          result += div;
          (*nit)->Store( this->key_C1, result );
      
      } // end for cycle for nodes

} // end AdjustResultsAssumingDivergenceFreeVelocityField





/** 
    Checks the range of 'new saturation oil' against that specified in the property database.
    Where the new values comply, they are used to replace the previous ones stored as 'saturation oil'.
    If not, the deviations are reported and the nearest maximum or minimum permitted values
    of saturation are stored.
*/
template<size_t dim>
double64 ExplicitTransport<dim>::VerifyAndAssignResults( bool show_range, bool do_range_check ) const
  {
    const typename vector<Node<dim>*>::iterator  nodes_end(subdomain_.NodesEnd());
    typename vector<Node<dim>*>::iterator nit = subdomain_.NodesBegin();

    double64        Cmin(upper_limit_),
                    Cmax(0.),
                    difference_to_last_output(0.);
    size_t          error_counter(0);
    VARIABLE_FLAG  C0_status(ANY);
   
    while ( nit != nodes_end )
       {
           // reading the newly computed saturation values
           double64 result = (*nit)->Read( this->key_C0 );
           Cmin = std::min( Cmin, result );
           Cmax = std::max( Cmax, result );

           // reading the previous values and determining whether nothing needs to be done because this is a Dirichlet node
           if ( (C0_status=(*nit)->Status(this->key_C0)) != DIRICH )
             {
               // calculating the maximum change per node
//               difference_to_last_output = std::max( difference_to_last_output, fabs(result - C0()) );

               // result checking and assignment
               if ( result <= upper_limit_ && result >= lower_limit_ ) (*nit)->Store( this->key_C0, makeScalar( C0_status, result ) );
               else {
                    cerr <<"\nExplicitTransport<dim>::VerifyAndAssignResults: ";
                    cerr <<"value: "<< result <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
                    // bracketing result
                    if ( result > upper_limit_ ) (*nit)->Store( this->key_C0, makeScalar( C0_status, upper_limit_ ) );
                    else if ( result < lower_limit_ ) (*nit)->Store( this->key_C0, makeScalar( C0_status, lower_limit_ ) );
                    error_counter++;
                 }
             }
          nit++;
       }

    if ( do_range_check ) {
          if ( error_counter > (this->subdomain_.Nodes()/10) )
            throw out_of_range("ExplicitTransport<dim>::VerifyAndAssignResults: Advected variable out of range.");
      }

    if ( show_range ) {
         cout <<"\n\nExplicitTransport<"<< dim;
         cout <<">::VerifyAndAssignResults: Variable range after advection: ";
         cout << Cmin <<" to "<< Cmax << endl << endl;
      }

    return difference_to_last_output / std::max( Cmax - Cmin, 1.0e-20 );

 } // end VerifyAndAssignResults





/**
    1. Element by element computation of (velocity), facet fluxes, flux balances and tracer balances
    
    2. Evaluation of time increment
 
    3. Assembly of solution

    4. Assignment of boundary conditions (treatment of truncated FVs at the boundary)
 
    5. tranfer of results into 'concentration', vacating 'new concentration' for 
       next assembly.
*/
template<size_t dim>
void ExplicitTransport<dim>::AdvectVariable( double64 time_interval )
 {
    // 1. element-by-element computation of volumetric and transport variable fluxes on FV facets
    //    within the computational domain and the facets of halo elements
    //
    // 2. accumulation of upstream variable flux products, and
    //    computation of the divergence of the volumetric flow on the finite volumes
    VolumetricFlowAndTransportVariableFluxBalances();

    // 2. node-by-node evaluation of time increment (no new variable values are stored)
    double64 time_increment = TimeIncrement_CFL_Outflow( this->MaxTimeIncrement() );

    cout <<"\nExplicitTransport<"<< fixed << setprecision(0) << dim <<">::AdvectVariable:";
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< std::max( rint(floor(time_interval/time_increment)), 1. );

    cout <<"\n\n\nExplicitTransport::EvolveSolution: Finite volume transport simulation initiated...\n";
    size_t   substep(1);
    double64 time(0.);

    // time incrementation loop
    while ( time < time_interval ) {
          cout <<"\n\n\tadvection (sub)step: "<< substep;
          if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;

          // 3. Solve C^t+1 = C^t - dt/(phi Vi) * sum_j^faces Aj n . [C vD]
          //     - a single loop over all nodes is required
          //     - this steps also considers in and outflow of finite volume cells
          //     - (distributed) sources and sinks will be considered
          const bool with_divergence_correction(true);
          // node loop: FV by FV
          Assemble1stOrderSolution( time_increment, with_divergence_correction );
      
// TESTING
double64 so1_min, so1_max;
subdomain_.MinMaxOf( "concentration", so1_min, so1_max );
cerr <<"\n\ttime-increment: "<< time_increment <<": range of assembled solution: "<< so1_min <<" to "<< so1_max << endl;

          // 4. node loop: 'new concentration' is used to replace 'concentration' performing a range check
          const bool range_check(true);
          const bool show_range(false);
          VerifyAndAssignResults( show_range, range_check );

          // updating facet flow - concentration products only (since velocity did not change)
          TransportVariableFluxBalances();

          time_increment = TimeIncrement_CFL_Outflow( this->MaxTimeIncrement() );

          time += time_increment;
          substep++;
      }

    cout <<"\nExplicitTransport::AdvectVariable: 'transport completed.\n";

 } // end AdvectVariable


// TESTING
/*
double64 gmin, gmax;
subdomain_.MinMaxOf( "facet flux", gmin, gmax );
cout <<"\nExplicitTransport<"<< fixed << setprecision(0) << dim <<">::AdvectVariable: range of 'facet flux': ";
cout << scientific << setprecision(5) << gmin <<" to "<< gmax;

{
  double64 vmin, vmax;
  subdomain_.MinMaxOf( "accumulation", vmin, vmax );
  cout <<"\n\t range of accumulation: "<< vmin <<" to "<< vmax << endl;
  subdomain_.MinMaxOf( "flux balance", vmin, vmax );
  cout <<"\n\t range of flux balance: "<< vmin <<" to "<< vmax << endl;
}

*/

/**
    Reports the volumetric flow into the computational region.
*/
template<size_t dim>
double64 ExplicitTransport<dim>::IncomingVolumetricFlow() const
 {
     double64 influx(0.);
 
     const typename vector<Node<dim>*>::const_iterator nodes_end(subdomain_.NodesEnd());
     for ( typename vector<Node<dim>*>::const_iterator
           nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
           if ( HasHaloElements() ) {
                // away from model boundaries, the inflow is computed correctly and can be used
                if ( (*nit)->AtBoundary() == NOT ) influx += this->InFlow( (*nit) );
                // for a FV truncated by the inflow boundary outflow is recorded
                else influx += this->OutFlow( (*nit) );
             }
           // at the model boundaries, only the inside part of the FV exists, registering inflows as outflows
           else {
                assert( (*nit)->AtBoundary() != NOT );
                influx += this->OutFlow( (*nit) );
             }
       }
 
    // incoming flux should be positive
    return influx;

 } // end IncomingVolumetricFlow
  
  
/**
    Reports the volumetric flow outside of the computational region.
 
    @todo SKM check whether the variable 'outflow' can be used for this, saving some computations.
*/
template<size_t dim>
double64 ExplicitTransport<dim>::OutgoingVolumetricFlow() const
 {
     double64 outflux(0.);
 
     const typename vector<Node<dim>*>::const_iterator nodes_end(subdomain_.NodesEnd());
     for ( typename vector<Node<dim>*>::const_iterator
           nit=subdomain_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
           if ( HasHaloElements() ) {
                // away from model boundaries, the inflow is computed correctly and can be used
                if ( (*nit)->AtBoundary() == NOT ) outflux += this->OutFlow( (*nit) );
                // for a FV truncated by the inflow boundary outflow is recorded
                else outflux += this->InFlow( (*nit) );
             }
           // at the model boundaries, only the inside part of the FV exists, registering inflows as outflows
           else {
                assert( (*nit)->AtBoundary() != NOT );
                outflux += this->InFlow( (*nit) );
             }
       }
 
    return outflux;
 
 } // end OutgoingVolumetricFlow




template class ExplicitTransport<3U>;

} // end csmp

