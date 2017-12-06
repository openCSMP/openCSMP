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

using namespace std;

namespace csmp {


template<size_t dim>
const Model<dim>&
ImplicitTransport<dim>::GetModel() const
{
    return model_;
}


template<size_t dim>
ImplicitTransport<dim>::ImplicitTransport( Model<dim>& m, const char* target_region, bool second_order )
  : VariableSet_TracerTransferImplicit(m.Database()),
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
    m.Database().RangeOf( m.Database().Name(this->C0_key), lower_limit_, upper_limit_ );
 }
  



/**
    Velocity and flux calculation (element by element), for all elements in the domain
*/
template<size_t dim>
void ImplicitTransport<dim>::UpdateFacetFluxes(bool reuse_previous_velocity)
 {

 } // end UpdateFacetFluxes



/**
    Velocity and flux calculation (element by element), for all elements in the domain
    1st order accurate in space and time
    
    @attention this means that facet fluxes in FV sectors outside the domain are not considered.
    This is done in TimeIncrementAndFluxBalance() by Advective_O1_FluxesBoundary().
    This method writes the correct balances onto the nodes on the perimeter.
*/
template<size_t dim>
void ImplicitTransport<dim>::UpdateFacetFluxes_O1(bool reuse_previous_velocity)
 {
 } // end UpdateFacetFluxes_O1


/**
    Velocity and flux calculation (element by element), for all elements in the domain
    2nd order accurate in space
    
    @attention this means that facet fluxes in FV sectors outside the domain are not considered.
    This is done in TimeIncrementAndFluxBalance() by Advective_O1_FluxesBoundary().
    This method writes the correct balances onto the nodes on the perimeter.
*/
template<size_t dim>
void ImplicitTransport<dim>::UpdateFacetFluxes_O2(bool reuse_previous_velocity)
 { 
 } // end UpdateFacetFluxes_O2

  








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



/** 
    Checks the range of 'new saturation oil' against that specified in the property database.
    Where the new values comply, they are used to replace the previous ones stored as 'saturation oil'.
    If not, the deviations are reported and the nearest maximum or minimum permitted values
    of saturation are stored.
*/
template<size_t dim>
double64 ImplicitTransport<dim>::VerifyAndAssignResults( bool show_range, bool do_range_check ) const
  {
      return 0.;
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

 } // end AdvectVariable



  template class ImplicitTransport<1U>;
  template class ImplicitTransport<2U>;
  template class ImplicitTransport<3U>;

} // end csmp

