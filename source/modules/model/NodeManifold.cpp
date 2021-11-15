#include "NodeManifold.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<size_t dim>
NodeManifold<dim>::~NodeManifold()
{
   // disconnecting the node pointers so that no damage is done
   for ( auto& nit : branches_ )
     nit.first = static_cast<Node<dim>*>(nullptr);
}



template<size_t dim>
NodeManifold<dim>::NodeManifold( const std::vector<Node<dim>*>& nodes,
                                 const std::vector<INTERFACE_SIDE>& sides,
                                 ManifoldType classifier )
 : parent_geometry_(classifier)
{
    const size_t n_branches(nodes.size());
    assert( n_branches >= 2 );
    assert( n_branches == sides.size() );
    branches_.reserve( n_branches );
    for( size_t i=0U; i<n_branches; ++i ) {
         branches_.push_back( make_pair( nodes[i], sides[i] ) );
         nodes[i]->Assign(this);
      }
    // sorting branches using the node pointers as keys (default of sort)
    sort( branches_.begin(), branches_.end() );
}




template<size_t dim>
NodeManifold<dim>::NodeManifold( const manifold& nodes, ManifoldType geometry )
 : branches_(nodes),
   parent_geometry_(geometry)
{
    assert( branches_.size() >= 2 );
    for( auto& nd : branches_ )
      nd.first->Assign(this);
}







/**
   Check whether the current topologic classification of the NodeManifold (still) makes sense
   
      Idea: dependent on topology there are only so many options
      
      @todo use more diagnostics, e.g., element types (line, surface, volume)
      @todo use also the material IDs that come together at that node
      @todo use std::initialiser_list<>  or something to create little comparitor functions that make the comparisons more readable
*/
template<size_t dim>
ManifoldType  consistencyCheck( const NodeManifold<dim>* const nmf  )
 {
    // diagnostics: go over all possible cases excluding the ones that are not possible
    const Node<dim>* const mnd_ptr = nmf->MIDDLE_Node();
    const bool   with_middle_node = (mnd_ptr == nullptr) ? false : true;
    const size_t collocated_nodes = nmf->Branches();
    const ManifoldType te = nmf->GeometricClassifier();

    // extended diagnostics: checking the parent elements for their types
    const size_t connected_elmts(mnd_ptr->Parents());
    size_t       line_counter(0U), surf_counter(0U);
//                 vol_counter(0U); // track inner and outer nodes
    // in the case of a fracture manifold (we can have an endpoint...)
    if ( with_middle_node )
      for ( size_t i=0U; i<connected_elmts; ++i ) {
           if ( mnd_ptr->Parent(i)->IsLineElement() ) line_counter++;
           else surf_counter++;
        }
    const bool line_elmt_manifold = (line_counter == connected_elmts) ? true : false;
    const bool surf_elmt_manifold = (surf_counter == connected_elmts) ? true : false;

    // most common: somewhere on a split boundary
    if ( collocated_nodes == 2 && !with_middle_node ) {
         // classification is plausible
         if ( te == ManifoldType::INTERFACE || te == ManifoldType::BEDGE || te == ManifoldType::BINTERSECTION_POINT ) return te;
         // current classification is probably not correct
         return ManifoldType::NOT_CLASSIFIED;
      }
      
    // if there is a middle node which is only connected to line elements
    if ( collocated_nodes == 2 && with_middle_node ) {
         // classification is plausible
         if ( line_elmt_manifold ) {
              if ( te == ManifoldType::POINT2 && connected_elmts == 2 ) return te;
              else return ManifoldType::NOT_CLASSIFIED;
              if ( te == ManifoldType::POINT3 && connected_elmts == 3 ) return te;
              else return ManifoldType::NOT_CLASSIFIED;
           }
         else if ( te == ManifoldType::INTERFACE || te == ManifoldType::BEDGE || te == ManifoldType::BINTERSECTION_POINT ) return te;
         // current classification is probably not correct
         return ManifoldType::NOT_CLASSIFIED;
      }
      
    // perimeter of SplitBoundary with fracture inside that terminates in a volume
    if ( collocated_nodes == 2 && with_middle_node ) {
         if ( te == ManifoldType::EDGE || te == ManifoldType::END_POINT || te == ManifoldType::INTERSECTION_POINT ) return te;
         else return ManifoldType::NOT_CLASSIFIED;
      }
    // split boundary T junction
    if ( collocated_nodes == 3 && !with_middle_node ) {
         if ( te == ManifoldType::INTERSECTION || te == ManifoldType::BINTERSECTION_POINT || te == ManifoldType::BPOINT2 || te == ManifoldType::BPOINTX ) return te;
         else return ManifoldType::NOT_CLASSIFIED;
      }
    // simple split boundary intersection
    if ( collocated_nodes >= 4 && !with_middle_node ) {
         if ( te == ManifoldType::INTERSECTION || te == ManifoldType::POINTX || te == ManifoldType::BINTERSECTION_POINT ) return te;
         else return ManifoldType::NOT_CLASSIFIED;
      }
    // multiple intersecting split boundaries
    if ( collocated_nodes >= 5 && !with_middle_node ) {
         if ( te == ManifoldType::MULTI_INTERSECTION || te == ManifoldType::POINTX || te == ManifoldType::BPOINTX ) return te;
         else return ManifoldType::NOT_CLASSIFIED;
      }
      
    // all other cases
    return te;
 }



// accessor
template<size_t dim>
ManifoldType NodeManifold<dim>::GeometricClassifier() const
 {return parent_geometry_;}

// mutator
template<size_t dim>
void NodeManifold<dim>::GeometricClassifier( ManifoldType te )
 { parent_geometry_ = te; }



// output of enum
string parse( ManifoldType topology )
 {
   switch( topology ) {
        case ManifoldType::POINT2: return "POINT2";                           ///<  intersection or contact point between 2 line SplitBoundaries
        case ManifoldType::POINT3: return "POINT3";                           ///<  3 line SplitBoundary objects in contact
        case ManifoldType::POINTX: return "POINTX";                           ///<  multiple SB in contact
        case ManifoldType::END_POINT: return "END_POINT";                     ///< where a line SplitBoundary terminates in a volumetric region
        case ManifoldType::EDGE: return "EDGE";                               ///<  outer edge of a SplitBoundary terminating within volume
        case ManifoldType::BPOINT1:  return "BPOINT1";                         ///< intersection between line SB and model boundary
        case ManifoldType::BEDGE:  return "BEDGE";                             ///< edge of SB located on model boundary
        case ManifoldType::BPOINT2:  return "BPOINT2";                         ///< intersection of boundary edge with model boundary
        case ManifoldType::BPOINTX:  return "BPOINTX";
        case ManifoldType::INTERSECTION:  return "INTERSECTION";               ///<  line where 2 SplitBoundary objects meet
        case ManifoldType::MULTI_INTERSECTION:  return "INTERSECTION";         ///<  more than 2 SplitBoundary objects meet
        case ManifoldType::INTERSECTION_POINT:  return "INTERSECTION_POINT";   ///< edge point due the contact between multiple SB objects
        case ManifoldType::BINTERSECTION_POINT: return "BINTERSECTION_POINT"; ///< SB's touch each other on model boundary
        case ManifoldType::INTERFACE:  return "INTERFACE";                 ///<  on SplitBoundary
        default: return "NOT_CLASSIFIED";
      }
    return "NOT_CLASSIFIED";                                           ///<  an isolated Point that has therefore no classification
 }


template<size_t dim>
void NodeManifold<dim>::SortByVariableValue( const Index& index )
{
    assert(index.type==SCALAR);
    assert(index.place==NODE || index.place==ELEMENT);

    if(index.type!=SCALAR)
      throw csmp::Exception( INFO, "NodeManifold<dim>::SortByVariableValue",
                                   "variable type not supported, must be SCALAR");
                                   
    const size_t n_branches(branches_.size());
    
    // if the variable is associated with the node manifold
    // (different values on each node of the manifold)
    if(index.place==NODE){
        for (size_t i = 0; i < n_branches; ++i) {
            double first = branches_[i].first->Read(index);
            assert(!isnan(first));
            for (size_t j = i + 1; j < n_branches; ++j) {
                double second = branches_[j].first->Read(index);
                assert(!isnan(second));
                if (first > second) std::swap(branches_[i], branches_[j]);
            }
        }
    // if the variable is associated with the parent elements
    } else if(index.place==ELEMENT) {
        for (size_t i = 0; i < branches_.size(); ++i) {
            //double first = std::numeric_limits<double>::quiet_NaN();
            double first_value(0.);
            size_t first_count(0);
            for(size_t e = 0; e < branches_[i].first->Parents(); e++) {
                //if( (dim==2 && nodes_[i]->Parent(e)->IsSurfaceElement()) || (dim==3 && nodes_[i]->Parent(e)->IsVolumeElement()) ) {
                if( !isnan(branches_[i].first->Parent(e)->Read(index)) ) {
                    //first = nodes_[i]->Parent(e)->Read(index);
                    //break;
                    first_value += branches_[i].first->Parent(e)->Read(index);
                    first_count ++;
                }
            }
            //assert(!isnan(first));
            assert(first_count != 0);
            first_value /= first_count;
            assert(!isnan(first_value));

            for (size_t j = i + 1; j < n_branches; ++j) {
                //double second = std::numeric_limits<double>::quiet_NaN();
                double second_value(0.);
                size_t second_count(0);                
                for(size_t e = 0; e < branches_[j].first->Parents(); e++) {
                    //if( (dim==2 && nodes_[j]->Parent(e)->IsSurfaceElement()) || (dim==3 && nodes_[j]->Parent(e)->IsVolumeElement()) ) {
                    if( !isnan(branches_[j].first->Parent(e)->Read(index)) ) {
                        //second = nodes_[j]->Parent(e)->Read(index);
                        //break;
                        second_value += branches_[j].first->Parent(e)->Read(index);
                        second_count ++;
                    }
                }
                //assert(!isnan(second));
                assert(second_count != 0);
                second_value /= second_count;
                assert(!isnan(second_value));                
                //if (first > second) std::swap(nodes_[i], nodes_[j]);
                if (first_value > second_value) std::swap(branches_[i], branches_[j]);
            }
        } 
    } else {
        throw csmp::Exception( INFO, "NodeManifold<dim>::SortByVariableValue",
                                     "variable placement not supported");      

    }
    
} // end SortByVariableValue





template<size_t dim>
size_t NodeManifold<dim>::Branches() const
{
  return branches_.size();
}




template<size_t dim>
Node<dim>* const NodeManifold<dim>::N( size_t index )
{
  return branches_[index].first;
}


/// where the node resides
template<size_t dim>
INTERFACE_SIDE NodeManifold<dim>::InterFaceSide( size_t branch ) const
 {
   return branches_[branch].second;
 }


template<size_t dim>
Node<dim>* const  NodeManifold<dim>::MIDDLE_Node() const
 {
    for ( auto& nit : branches_ )
      if ( nit.second == MIDDLE )
        return nit.first;
        
    return static_cast<Node<dim>* const>(nullptr);
 }




/** retrieve node(s) with specific flags (manifold should only contain one intervening Node (MIDDLE)
         @note due to the move() semantics there won't be redundant copying of the vector
 */
template<size_t dim>
vector<Node<dim>*>  NodeManifold<dim>::NodesLocatedAt( INTERFACE_SIDE side ) const
 {
    vector<Node<dim>*> temp;
    temp.reserve( branches_.size() );
    for ( auto& nit : branches_ )
      if ( nit.second == side )
        temp.push_back( nit.first );
    
    return temp;
 }




template<size_t dim>
bool NodeManifold<dim>::Add( Node<dim>* nd, INTERFACE_SIDE side, ManifoldType geometry )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    //a node can only exist in one manifold
    if(nd->Manifold() != nullptr) {
        csmp_error.notice( INFO, "NodeManifold<dim>::Add(Node<dim>*)",
                                 "Node already in another manifold. Nothing was done.");         
        return false;
    }

    //the node is already in current manifold
    for ( auto& nit : branches_ )
      if ( nit.first == nd ) {
          csmp_error.notice( INFO, "NodeManifold<dim>::Add(Node<dim>*)",
                                   "Node already part of current manifold. Nothing was done.");
          return false;
        }

    // making sure that the vector does not grow by multiples of 2
    if ( branches_.capacity() == branches_.size() ) branches_.reserve( branches_.size() + 1 );
    branches_.push_back( make_pair(nd,side) );
    nd->Assign(this);
    
    // reviewing topology of node after insertion to see whether change is necessary
    parent_geometry_ = consistencyCheck( this );

    return true;

}




template<size_t dim>
bool NodeManifold<dim>::Remove( Node<dim>* nd )
{
   for ( auto& nit : branches_ )
     if ( nit.first == nd ) {
         // disconnecting the node from the manifold
         nd->Assign( static_cast<NodeManifold<dim>*>(nullptr) );
         // removing the Node entry from the branch list
         branches_.erase( remove( branches_.begin(), branches_.end(), nit ) );
         return true;
       }

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.notice( WARNING, "NodeManifold<dim>::Remove(Node<dim>*)",
                                "Node does not exist in manifold. Nothing was done.");
    return false;

}




template<size_t dim>
void NodeManifold<dim>::Out() const
{
  std::cout << "NodeManifold: "<< parse(parent_geometry_) <<" with nodes with the IDs:\t";
  for (auto const& node : branches_ )
    {
      std::cout <<"\n\t"<< parseSide(node.second) <<": "<< node.first->Idx() << "\t";
    }
  std::cout << std::endl;
}

template class NodeManifold<1U>;
template class NodeManifold<2U>;
template class NodeManifold<3U>;

}// csmp
