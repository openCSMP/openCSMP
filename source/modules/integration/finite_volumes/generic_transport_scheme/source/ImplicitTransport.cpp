//
//  ImplicitTransport.h
//
//  Created by Andrew J. Bromage on 5/12/2017.
//  Copyright (c) 2017 The University of Melbourne. All rights reserved.
//

#include "ImplicitTransport.h"
#include "Region.h"
#include "Model.h"
#include "finiteVolumeFunctions.h"
#include "CSMP_mathUtilities.h"
#include "VTK_Interface.h"
#include "LinearSystemAccumulator.h"

using namespace std;

namespace csmp {


template<size_t dim>
const Model<dim>&
ImplicitTransport<dim>::GetModel() const
{
    return model_;
}


template<size_t dim>
ImplicitTransport<dim>::ImplicitTransport( Solver& solver, Model<dim>& m, const char* target_region, bool second_order )
  : variables::Variables_TracerTransfer(m.Database()),
    Equation_TracerTransferImplicit<dim>( m, target_region ),
    solver_(solver),
    model_(m),
    gref_(m.Region(target_region)),
    upper_limit_(1.), lower_limit_(0.),
    second_order_(second_order)
 {
    m.InstantiateFiniteVolumes();
     m.Region(target_region).InputPropertyValue( "new concentration", makeScalar(PLAIN,0.), COMPLETE );

    initializeFiniteVolumeProperties( m, m.Region(target_region) );
    // 0. model-wide initialisation: results will be accumulated into this variable
   
    // retrieving the physically meaningful upper and lower solution limit from database
    m.Database().RangeOf( m.Database().Name(this->key_C), lower_limit_, upper_limit_ );
 }
  
/**
    Computation of time increment, flux balance, and temporary new concentration.
*/
template<size_t dim>
double64 ImplicitTransport<dim>::TimeIncrementAndFluxBalance( double64 max_time_increment )
 {
     double64 dt_min(max_time_increment);
 
     // 1. processing interior and FVs for which all facet fluxes have been initialised
     auto interior_nodes_end(gref_.InteriorNodesEnd());
     for ( auto nit=gref_.InteriorNodesBegin(); nit!=interior_nodes_end; ++nit )
       {
          assert( (*nit)->AtBoundary() == NOT );
          // computes time-increment, flux balance, and flux-concentration product balance
          const double64 time_increment = this->OutFlowLessThanContentIncrement( (*nit) );
          dt_min = std::min( dt_min, time_increment );
       }

     // 2. collecting time-stepping constraints from FVs on region perimeter
     const typename vector<Node<dim>*>::const_iterator nodes_end(gref_.PerimeterNodesEnd());
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
double64 ImplicitTransport<dim>::TimeIncrement()
 {
    return TimeIncrementAndFluxBalance( 356. * 86400. );
 }

  
  
  
  template<size_t dim>
  void ImplicitTransport<dim>::AccumulateSystem( double64 dt )
  {
    // Reset the linear system
    const size_t nodes = gref_.Nodes();
    LHS_.Erase();
    LHS_.Resize( nodes );
    
    if ( RHS_.size() != nodes ) {
      RHS_.resize( nodes );
      vector<double64>(RHS_).swap(RHS_);
    }
    fill( RHS_.begin(), RHS_.end(), static_cast<double64>(0.) );
    
    if ( RESULT_.size() != nodes ) {
      RESULT_.resize( nodes );
      vector<double64>(RESULT_).swap(RESULT_);
    }
    fill( RESULT_.begin(), RESULT_.end(), static_cast<double64>(0.) );
    
    this->EquationTimeIncrement(dt);
    this->accumulator.AccumulateByStencil(LHS_, RHS_);
} // end AccumulateSystem

  
  template<size_t dim>
  void ImplicitTransport<dim>::AssignBoundaryConditions()
  {
    const auto nodes_end(gref_.PerimeterNodesEnd());
    for ( auto nit = gref_.PerimeterNodesBegin(); nit != nodes_end; ++nit ) {
      const auto i = (*nit)->Idx();
      
      const auto status = (*nit)->Status( key_C );
      
      // Dirichlet boundary condition: concentration should be unaltered.
      if (status == DIRICH) {
        const auto c0 = (*nit)->Read( key_C );
        LHS_.ZeroRow(i);
        LHS_.Assign(i, i, 1.0);
        RHS_[i] = c0;
        continue;
      }
      
      // Calculate flow through boundary
      double64 inflow = 0.0;
      double64 influx = 0.0;
      const auto c0 = (*nit)->Read( key_C );

      const size_t iNrParents = (*nit)->Parents();
      for ( size_t iParent = 0; iParent < iNrParents; ++iParent ) {
        const size_t pnid = (*nit)->ParentNodeNumber( iParent );
        const auto eptr = (*nit)->Parent( iParent );
        
        const double64 K(eptr->Read(key_k));
        if (isnan(K)) {
          // XXX AJB HACK
          // Deleting boundaries during model creation means you can't set
          // properties during configuration. Retaining the boundaries means
          // that they are still included in the parents of a node.
          //
          // For now, we skip over any element which doesn't have a
          // permeability. They are not the flow domain.
          continue;
        }

        const size_t iNrSectorFacets(eptr->FV()->FacetsPerSector(pnid));
        for ( size_t iSectorFacet=0U; iSectorFacet<iNrSectorFacets; ++iSectorFacet ) {
          const size_t iFacet( eptr->FV()->FacetSurroundingSector(pnid,iSectorFacet) );
          const size_t inside_node(eptr->FV()->InsideNode(iFacet));
          const size_t outside_node(eptr->FV()->OutsideNode(iFacet));
          const double64 ff = eptr->Read( iFacet, 0u, key_ff );
          
          const double64 C_upstream = (ff < 0.) ? eptr->N(outside_node)->Read( key_C ) :
          eptr->N(inside_node)->Read( key_C );
          if ( pnid == inside_node ) {
            inflow += ff;
            influx += ff * C_upstream;
          }
          else {
            inflow -= ff;
            influx -= ff * C_upstream;
          }
        }
      }

      if ( inflow >= 0 ) {
        // inflow compensation
        RHS_[i] += inflow * c0;
      }
      else {
        // outflow compensation
        LHS_.Add( i, i, -inflow );
      }
    }
  }
  
  template<size_t dim>
  void ImplicitTransport<dim>::Solve()
  {
    solver_.Solve(LHS_, RHS_, RESULT_);
  } // end Solve
  


/** 
    Checks the range of 'new saturation oil' against that specified in the property database.
    Where the new values comply, they are used to replace the previous ones stored as 'saturation oil'.
    If not, the deviations are reported and the nearest maximum or minimum permitted values
    of saturation are stored.
*/
template<size_t dim>
double64 ImplicitTransport<dim>::VerifyAndAssignResults( bool show_range, bool do_range_check ) const
  {
    double64 amin(+std::numeric_limits<double64>::max());
    double64 amax(-std::numeric_limits<double64>::max());
    double64 difference_to_last_output(0.);
    size_t          error_counter(0);
    
    const typename vector<Node<dim>*>::iterator  nodes_end(gref_.NodesEnd());
    for ( auto nit = gref_.NodesBegin(); nit != nodes_end; ++nit )
    {
      const VARIABLE_FLAG status((*nit)->Status( this->key_C ));
      if ( status != DIRICH )
      {
        // reading the newly computed saturation values
        const double64 c1 = RESULT_[ (*nit)->Idx() ];
        amin = std::min( amin, c1 );
        amax = std::max( amax, c1 );
        
        // reading the previous values and calculating the maximum change per node
        const double64 C0 = (*nit)->Read( this->key_C );

        difference_to_last_output = std::max( difference_to_last_output, fabs(c1 - C0) );
        
        // result checking and assignment
        if ( c1 <= upper_limit_ && c1 >= lower_limit_ ) (*nit)->Store( this->key_C, makeScalar(status, c1) );
        else {
          cerr <<"\nExplicitTransport<dim>::VerifyAndAssignResults: ";
          cerr <<"value: "<< c1 <<" versus range from PropertyDatabase: "<< lower_limit_ <<"-"<< upper_limit_ << endl;
          if ( c1 > upper_limit_ ) (*nit)->Store( this->key_C, makeScalar( status, upper_limit_ ) );
          else if ( c1 < lower_limit_ ) (*nit)->Store( this->key_C, makeScalar( status, lower_limit_ ) );
          error_counter++;
        }
      }
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
void ImplicitTransport<dim>::AdvectVariable( double64 time_interval )
 {
   // 1. computing (velocity and) facet fluxes as necessary
   const bool reuse_previous_velocity = false;
   this->FacetFluxes(gref_, reuse_previous_velocity, second_order_);
   
   // 2. evaluation of time increment
   double64 time_increment = TimeIncrementAndFluxBalance( this->MaxTimeIncrement() );
   
   const double64 one(1.);
   cout <<"\nImplicitTransport<"<< fixed << setprecision(0) << dim <<">::AdvectVariable:";
   cout <<"\n\tTime interval         = "<< time_interval;
   cout <<"\n\tScaled time increment = "<< time_increment;
   cout <<"\n\tSolution steps needed = "<< max(floor(time_interval/time_increment),one);
   
   cout <<"\n\n\nImplicitTransport::AdvectVariable: FV transport simulation initiated...\n";
   size_t   substep(1);
   double64 time(0.);
   
   // time incrementation loop
   while ( time < time_interval ) {
     cout <<"\n\n\tadvection (sub)step: "<< substep;
     
     if ( (time_interval - time) < time_increment ) time_increment = time_interval - time;
     
     AccumulateSystem(time_increment);
     AssignBoundaryConditions();
     Solve();
     
#if 1
     // testing
     double64 so1_min = +std::numeric_limits<double64>::max();
     double64 so1_max = -std::numeric_limits<double64>::max();
     
     for (auto so : RESULT_) {
       so1_min = std::min(so1_min, so);
       so1_max = std::max(so1_max, so);

     }
     cerr <<"\n\ttime-increment: "<< time_increment <<": range of assembled solution: "<< so1_min <<" to "<< so1_max << endl;
#endif
     
     // 4. 'new saturation oil' is used to replace 'saturation oil' performing a range check
     const bool range_check(true);
     const bool show_range(true);
     VerifyAndAssignResults( show_range, range_check );
     
     time += time_increment;
     
     const bool reuse_previous_velocity = true;
     this->FacetFluxes(gref_, reuse_previous_velocity, second_order_);
     
     time_increment = TimeIncrementAndFluxBalance( this->MaxTimeIncrement() );
     
     substep++;
   }
   
   cout <<"\nExplicitTransport::AdvectVariable: 'transport completed.\n";
   

 } // end AdvectVariable



  template class ImplicitTransport<1U>;
  template class ImplicitTransport<2U>;
  template class ImplicitTransport<3U>;

} // end csmp

