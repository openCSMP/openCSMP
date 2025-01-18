#include "ElementToNodePropertyVisitor.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "TensorVariable.h"
#include "Exception.h"
#include "ErrorHandler.h"

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
ElementToNodePropertyVisitor<Var,dim>::ElementToNodePropertyVisitor( 
                                                             const PropertyDatabase<dim>& p,
                                                             const char* elmt_prop,
                                                             const char* node_prop,
                                                             size_t nodes )
    : eprop_key_(p.StorageKey(elmt_prop)),
      nprop_key_(p.StorageKey(node_prop)),
      summed_weights_(nodes,0.),
      weighting_completed_(nodes)
  { 
     weighting_completed_.SetAll(false);
  }



template<typename Var, uint32_t dim>
void ElementToNodePropertyVisitor<Var,dim>::ApplyWeightingToExtrapolatedValues()
 {
    
    
     if ( summed_weights_[0] > numeric_limits<double>::epsilon() )
       throw csmp::Exception( ERROR, "ElementToNodePropertyVisitor::ApplyWeightingToExtrapolatedValues",
                      "this visitor has to be applied beforehand to collect the weights for the nodes");
       
     this->ApplicationTarget(NODE);
 }


template<typename Var, uint32_t dim>
ElementToNodePropertyVisitor<Var,dim>::~ElementToNodePropertyVisitor() 
 {  
 }


/// first application cycle (extrapolation)
template<typename Var, uint32_t dim>
void ElementToNodePropertyVisitor<Var,dim>::Visit( Element<dim>* eptr ) 
 { 
    eptr->Read( eprop_key_, evariable_ );

    // finding the location of the Element barycenter
    csmp::Point<dim>  bc = eptr->BaryCenter();
    
    // accumulating results into nodes vector for later averaging
    for ( auto i{0}; i<eptr->Nodes(); i++ ) {
         // using 1 / (distance from barycenter to node)  as a weight
         const double weight = 1. / (bc - eptr->N(i)->Coordinate()).Length();
         assert( eptr->N(i)->Idx() < summed_weights_.size() );
         summed_weights_[ eptr->N(i)->Idx() ] += weight;
         
         // adding extrapolated contribution to node variable
         eptr->N(i)->Read( nprop_key_, nvariable_ );
         nvariable_ += evariable_ * weight;
         eptr->N(i)->Store( nprop_key_, nvariable_ );
      }
      
 } // end Visit Element



/// second application cycle (weighting)
template<typename Var, uint32_t dim>
void ElementToNodePropertyVisitor<Var,dim>::Visit( Node<dim>* nptr ) 
 { 
    if ( !weighting_completed_.GetBit( nptr->Idx() ) ) {
          nptr->Read( nprop_key_, nvariable_ );
          nvariable_ /= summed_weights_[ nptr->Idx() ];
          nptr->Store( nprop_key_, nvariable_ );
          weighting_completed_.SetBit( nptr->Idx(), true );
      }
 }



template class ElementToNodePropertyVisitor<ScalarVariable,1U>;
template class ElementToNodePropertyVisitor<VectorVariable<1U>,1U>;
template class ElementToNodePropertyVisitor<TensorVariable<1U>,1U>;

template class ElementToNodePropertyVisitor<ScalarVariable,2U>;
template class ElementToNodePropertyVisitor<VectorVariable<2U>,2U>;
template class ElementToNodePropertyVisitor<TensorVariable<2U>,2U>;

template class ElementToNodePropertyVisitor<ScalarVariable,3U>;
template class ElementToNodePropertyVisitor<VectorVariable<3U>,3U>;
template class ElementToNodePropertyVisitor<TensorVariable<3U>,3U>;

} // end namespace csmp
