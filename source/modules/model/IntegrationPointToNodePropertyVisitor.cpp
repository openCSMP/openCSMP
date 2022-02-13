#include "IntegrationPointToNodePropertyVisitor.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "TensorVariable.h"
#include "Exception.h"
#include <limits>

using namespace std;

namespace csmp {

/**

Steps:

- Visitor requires all nodes in the mesh to be numbered consecutively

- the target node property has to be set to zero prior to the accumulation

- Values are interpolated to nodes where these are accumulated 

- the nodal contributions are weighted by the distance of the element barycenter to each
  node (all weights are summed so that the absolute value can be rectified once the 
  accumulation is complete).

*/
template<typename Var, uint32_t dim>
IntegrationPointToNodePropertyVisitor<Var,dim>::IntegrationPointToNodePropertyVisitor( const PropertyDatabase<dim>& p, 
                                                                                     const char* cpoint_prop,
                                                                                     const char* node_prop,
                                                                                     size_t nodes )
    : cprop_key_(p.StorageKey(cpoint_prop)),
      nprop_key_(p.StorageKey(node_prop)),
      summed_weights_(nodes,0.),
      weighting_completed_(nodes)
  { 
     this->ApplicationLevel(REGION);
     this->ApplicationTarget(ELEMENT);
     weighting_completed_.SetAll(false);
  }



template<typename Var, uint32_t dim>
void IntegrationPointToNodePropertyVisitor<Var,dim>::ApplyWeightingToExtrapolatedValues()
 {
     
    if ( summed_weights_[0] > numeric_limits<double>::epsilon() )
       throw csmp::Exception( ERROR, "IntegrationPointToNodePropertyVisitor::ApplyWeightingToExtrapolatedValues",
                          "this visitor has to be applied beforehand to collect the weights for the nodes");

    if ( weighting_completed_.GetBit(0U) )
       throw csmp::Exception( WARNING, "IntegrationPointToNodePropertyVisitor::ApplyWeightingToExtrapolatedValues",
                      "weighting appears to have been carried out already");
       
    this->ApplicationTarget(NODE);
 }


template<typename Var, uint32_t dim>
IntegrationPointToNodePropertyVisitor<Var,dim>::~IntegrationPointToNodePropertyVisitor() 
 {  
 }


/// first application cycle (extrapolation)
template<typename Var, uint32_t dim>
void IntegrationPointToNodePropertyVisitor<Var,dim>::Visit( Element<dim>* eptr ) 
 { 
    // storing the constraint point property into a vector
    eptr->IntegrationPointPropertyVector( cprop_key_, vars_vector_ );
    
    const size_t  vcomponents(vars_vector_[0].Size());
    cp_vars_.resize( eptr->IntegrationPoints() * vcomponents );
    nd_vars_.resize( eptr->Nodes() * vcomponents );
    
    for ( auto i{0}; i<vars_vector_.size(); i++ )
      for ( size_t j=0U; j<vcomponents; j++ ) 
        cp_vars_[ i * vcomponents + j ] = vars_vector_[i].Component(j);
    
    // extrapolating constraint point properties to nodes
    eptr->ExtrapolateIntegrationPointVariableToNodes( vcomponents, cp_vars_, nd_vars_ );
    
    // finding the location of the Element barycenter
    csmp::Point<dim>  bc = eptr->BaryCenter();
    
    // accumulating results into nodes vector for later averaging
    for ( auto i{0}; i<eptr->Nodes(); i++ ) {
         eptr->N(i)->Read( nprop_key_, variable_ );
         
         // using 1 / (distance from barycenter to node)  as a weight
         const double weight = 1. / (bc - eptr->N(i)->Coordinate()).Length();
         assert( eptr->N(i)->Idx() < summed_weights_.size() );
         summed_weights_[ eptr->N(i)->Idx() ] += weight;
         
         // adding extrapolated contributions
         for ( size_t j=0U; j<vcomponents; j++ ) 
           variable_.Component( j, variable_.Component(j) + weight * nd_vars_[ i * vcomponents + j ] );
         
         eptr->N(i)->Store( nprop_key_, variable_ );
      }
      
 } // end Visit Element



/// second application cycle (weighting)
template<typename Var, uint32_t dim>
void IntegrationPointToNodePropertyVisitor<Var,dim>::Visit( Node<dim>* nptr ) 
 { 
    if ( !weighting_completed_.GetBit( nptr->Idx() ) ) {
          nptr->Read( nprop_key_, variable_ );
          variable_ /= summed_weights_[ nptr->Idx() ];
          nptr->Store( nprop_key_, variable_ );
          weighting_completed_.SetBit( nptr->Idx(), true );
      }
 }



template class IntegrationPointToNodePropertyVisitor<ScalarVariable,1U>;
template class IntegrationPointToNodePropertyVisitor<VectorVariable<1U>,1U>;
template class IntegrationPointToNodePropertyVisitor<TensorVariable<1U>,1U>;

template class IntegrationPointToNodePropertyVisitor<ScalarVariable,2U>;
template class IntegrationPointToNodePropertyVisitor<VectorVariable<2U>,2U>;
template class IntegrationPointToNodePropertyVisitor<TensorVariable<2U>,2U>;

template class IntegrationPointToNodePropertyVisitor<ScalarVariable,3U>;
template class IntegrationPointToNodePropertyVisitor<VectorVariable<3U>,3U>;
template class IntegrationPointToNodePropertyVisitor<TensorVariable<3U>,3U>;

} // end namespace csmp
