//
//  ExplicitTransport.cpp
//
//  Created by Stephan Matthai on 2/21/13.
//  Copyright (c) 2013 Stephan Matthai. All rights reserved.
//

#include "ExplicitTransport.h"
#include "Region.h"
#include "Model.h"
#include "finiteVolumeFunctions.h"
#include "CSMP_mathUtilities.h"
#include "VTK_Interface.h"

using namespace std;

namespace csmp {

template<size_t dim>
ExplicitTransport<dim>::ExplicitTransport( Model<dim>& m, const char* target_region )
  : VariableSet_TracerTransferExplicit<dim>(m.Database()),
    gref_(m.Region(target_region)),
    upper_limit_(1.), lower_limit_(0.)
 {
    m.InstantiateFiniteVolumes();
    initializeFiniteVolumeProperties( m, m.Region(target_region) );
    // 0. model-wide initialisation: results will be accumulated into this variable
    m.Region("Model").InputPropertyValue( "new concentration", makeScalar(PLAIN,0.), COMPLETE );
   
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->C0_key), lower_limit_, upper_limit_ );
 }
  



/**
    Velocity and flux calculation (element by element), for all elements in the domain
    
    @attention this means that facet fluxes in FV sectors outside the domain are not considered.
    This is done in TimeIncrementAndFluxBalance() by Advective_O1_FluxesBoundary().
    This method writes the correct balances onto the nodes on the perimeter.
*/
template<size_t dim>
void ExplicitTransport<dim>::UpdateFacetFluxes()
 {
     // 1. element-by-element processing of the facet fluxes
     const typename vector<Element<dim>*>::iterator elements_end(gref_.ElementsEnd());
     for ( typename vector<Element<dim>*>::iterator
           eit=gref_.ElementsBegin(); eit!=elements_end; ++eit )
       {
          // 1.1 computation of transport velocity from fluid pressure gradient
       
          // 1.2 computation of facet fluxes (including upstream concentrations, but no-time increment yet)
          this->Advective_O1_FluxesInterior( (*eit) );
       }

     // 2. processing fluxes through the FVs on regions perimeter computing outside facet fluxes as necessary
     const typename vector<Node<dim>*>::iterator nodes_end(gref_.NodesEnd());
     for ( typename vector<Node<dim>*>::iterator
           nit=gref_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
          // computing flux balances and concentration-facet flux products where possible,
          // at sliced boundaries 3-typed of conditions are applied: 1) prescribed value (only at inflow),
          // 2) prescribed flux (has consequence only where there is inflow), 3) free outflow (outflow)
          // in this case the influx is found from the flux balance and FV cell's concentration
          this->Advective_O1_FluxesAtBoundary( (*nit) );
       }
   
 } // end UpdateFlowTerms

  








/**
    Computation of time increment, flux balance, and temporary new concentration.
*/
template<size_t dim>
double64 ExplicitTransport<dim>::TimeIncrementAndFluxBalance( double64 max_time_increment )
 {
     double64 dt_min(max_time_increment);
 
     // 1. processing interior and FVs for which all facet fluxes have been initialised
     const typename vector<Node<dim>*>::iterator interior_nodes_end(gref_.PerimeterNodesBegin());
     for ( typename vector<Node<dim>*>::iterator
           nit=gref_.NodesBegin(); nit!=interior_nodes_end; ++nit )
       {
          assert( (*nit)->AtBoundary() == NOT );
          // computes time-increment, flux balance, and flux-concentration product balance
          const double64 time_increment = this->OutFlowLessThanContentIncrement( (*nit) );
          dt_min = std::min( dt_min, time_increment );
       }

     // 2. collecting time-stepping constraints from FVs on region perimeter
     const typename vector<Node<dim>*>::const_iterator nodes_end(gref_.NodesEnd());
     for ( typename vector<Node<dim>*>::const_iterator
           nit=gref_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
       {
          // boundary fluxes must be part of the time-increment calculation
          dt_min = std::min( dt_min, this->OutFlowLessThanContentIncrementBoundary( (*nit) ) );
       }
   
    return dt_min;
 }





/**
    Computation of time increment (default limited to 1 year).
*/
template<size_t dim>
double64 ExplicitTransport<dim>::TimeIncrement()
 {
    UpdateFacetFluxes();
    return TimeIncrementAndFluxBalance( 356. * 86400. );
 }








/**  AssembleSolution()

      S^t+1 = S^t - dt/(phi Vi) * (sum_j^faces Aj n . [f vt] + sources/sinks)
 
    - a single loop over all nodes is required
    - this steps also considers in and outflow of finite volume cells
    - (distributed) sources and sinks will be considered in the future
    
    @todo update water saturations after every step
*/
template<size_t dim>
void ExplicitTransport<dim>::AssembleSolution( double64 delta_t,
                                               bool enforce_divergence_free_vt_field )
 {
    const typename vector<Node<dim>*>::const_iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::const_iterator nit=gref_.NodesBegin(); nit!=nodes_end; ++nit )
      {
         // 1. starting with the sum of facet flux-concentration products stored in 'new concentration'
         double64 accumulation = (*nit)->Read(this->C1_key);
        
         // 2. correcting this sum for div vD using 'flux balance' except for at model boundary
         //if ( (*nit)->AtBoundary() == NOT ) accumulation -= accumulation * (*nit)->Read( this->fb_key );
         accumulation -= accumulation * (*nit)->Read( this->fb_key );
        
         // 3. ACCUMULATION: subtracting flux time-interval products from concentration at previous time level
         accumulation = (*nit)->Read(this->C0_key) - (delta_t/(*nit)->Read(this->PV_key)) * accumulation;

         // 4. accounting for absolute 'nodal fluid volume source' terms or sinks after the advection step
         // TODO: make this more accurate using a fractional step method where the source is accounted for at 2 time levels using dt/2 and C0 and C1
         const double64 source((*nit)->Read(this->nsrc_key));
         //                               new concentration
         accumulation += source * (*nit)->Read(this->C1_key) * delta_t;
        
         // 5. storing the new concentration
         (*nit)->Store( this->C1_key, makeScalar( (*nit)->Status(this->C1_key), accumulation ) );
    }
   
} // end AssembleSolution






/**
    For all nodes on the perimeter of the target region, this method
    loops over the finite-element sectors of the associated finite volume
    integrating facet fluxes based on total velocity.
    
    This flux balance is added to the node variable 'new saturation oil'
*/
template<size_t dim>
void ExplicitTransport<dim>::AdjustResultsAssumingDivergenceFreeVelocityField( double64 time_interval )
 {
    const typename vector<Node<dim>*>::iterator  nodes_end(gref_.NodesEnd());
    for ( typename vector<Node<dim>*>::iterator nit=gref_.PerimeterNodesBegin(); nit!=nodes_end; ++nit )
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
                    double64 velo = eptr->ProjectionOnFacetNormal( iFacet, this->vD_key );
                    if ( nid == eptr->FV()->InsideNode(iFacet) )div += velo;
                    else div -= velo;
                }
           }
          ScalarVariable result( makeScalar( (*nit)->Status(this->C1_key), (*nit)->Read(this->C1_key) ) );
          result() += div;
          (*nit)->Store( this->C1_key, result );
      
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
    const typename vector<Node<dim>*>::iterator  nodes_end(gref_.NodesEnd());
    typename vector<Node<dim>*>::iterator nit = gref_.NodesBegin();

    double64        amin(1.),  // = (*nit)->Read( this->so1_key ),
                    amax(0.),  // = (*nit)->Read( this->so1_key ),
                    difference_to_last_output(0.);
    size_t          error_counter(0);
    ScalarVariable  C1;
   
    while ( nit != nodes_end )
       {
          const VARIABLE_FLAG status((*nit)->Status( this->C0_key ));
          if ( status != DIRICH )
            {
               // reading the newly computed saturation values
               (*nit)->Read( this->C1_key, C1 );
               amin = std::min( amin, C1() );
               amax = std::max( amax, C1() );

               // reading the previous values and calculating the maximum change per node
               const double64 C0 = (*nit)->Read( this->C0_key );
               difference_to_last_output = std::max( difference_to_last_output, fabs(C1() - C0) );
  
               // result checking and assignment
               if ( C1() <= upper_limit_ && C1() >= lower_limit_ ) (*nit)->Store( this->C0_key, C1 );
               else {
                    cerr <<"\nExplicitTransport<dim>::VerifyAndAssignResults: ";
                    cerr <<"value: "<< C1() <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
                    if ( C1() > upper_limit_ ) (*nit)->Store( this->C0_key, makeScalar( status, upper_limit_ ) );
                    else if ( C1() < lower_limit_ ) (*nit)->Store( this->C0_key, makeScalar( status, lower_limit_ ) );
                    error_counter++;
                 }
            }
          nit++;
       }

    if ( do_range_check ) {
          if ( error_counter > (this->gref_.Nodes() * 1000000U) )
            throw out_of_range("ExplicitTransport<dim>::VerifyAndAssignResults: Advected variable out of range");
      }

    if ( show_range ) {
         cout <<"\n\nExplicitTransport<"<< dim;
         cout <<">::VerifyAndAssignResults: Variable range after advection: ";
         cout << amin <<" to "<< amax << endl << endl;
      }

    return difference_to_last_output / std::max( amax - amin, 1.0e-20 );

 } // end VerifyAndAssignResults





/**
    1. Computation of (velocity and) facet fluxes
    
    2. Evaluation of time increment, flux balances and flux balance concentration products
 
    3. Assembly of solution

    4. Assignment of boundary conditions (treatment of truncated FVs at the boundary)
 
    5. tranfer of results into 'concentration', vacating 'new concentration' for 
       next assembly.
*/
template<size_t dim>
void ExplicitTransport<dim>::AdvectVariable( double64 time_interval )
 {
    // 1. computing (velocity and) facet fluxes as necessary
    UpdateFacetFluxes();
   
    // 2. evaluation of time increment
    double64 time_increment = TimeIncrementAndFluxBalance( this->MaxTimeIncrement() );
   
    cout <<"\nExplicitTransport<"<< fixed << setprecision(0) << dim <<">::AdvectVariable:";
    cout <<"\n\tTime interval         = "<< time_interval;
    cout <<"\n\tScaled time increment = "<< time_increment;
    cout <<"\n\tSolution steps needed = "<< std::max(rint(floor(time_interval/time_increment)),1);

    cout <<"\n\n\nExplicitTransport::AdvectVariable: FV transport simulation initiated...\n";
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
          const bool with_divergence_correction(false);
          AssembleSolution( time_increment, with_divergence_correction );
      
// testing
double64 so1_min, so1_max;
gref_.MinMaxOf( "new concentration", so1_min, so1_max );
cerr <<"\n\ttime-increment: "<< time_increment <<": range of assembled solution: "<< so1_min <<" to "<< so1_max << endl;

          // 4. 'new saturation oil' is used to replace 'saturation oil' performing a range check
          const bool range_check(true);
          const bool show_range(true);
          VerifyAndAssignResults( show_range, range_check );

          UpdateFacetFluxes();
      
          time_increment = TimeIncrementAndFluxBalance( this->MaxTimeIncrement() );

          time += time_increment;
          substep++;
      }

    cout <<"\nExplicitTransport::AdvectVariable: 'transport completed.\n";

 } // end AdvectVariable









template class ExplicitTransport<3U>;

} // end csmp

