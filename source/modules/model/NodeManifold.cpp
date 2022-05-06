#include "NodeManifold.h"
#include "plf_colony.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
NodeManifold<dim>::~NodeManifold()
{
   // before destruction you should disconnect the nodes from the manifold, like so
   // for ( auto& nit : branches_ )
   //  nit.first = static_cast<Node<dim>*>(nullptr);
}


// tested: OK
template<uint32_t dim>
NodeManifold<dim>::NodeManifold( plf::colony<Node<dim> >& nodes,
                                 const set<pair<size_t,INTERFACE_SIDE> >& manifold_nodes,
                                 ManifoldType classifier )
 : parent_geometry_(classifier)
{
    const size_t n_branches(manifold_nodes.size());
    assert( n_branches >= 2 );
    branches_.reserve( n_branches );

    for ( const auto& nit : manifold_nodes ) {
         // connecting the new manifold to its nodes
         branches_.emplace_back( make_pair( &(*next(nodes.begin(),nit.first)), nit.second ) );
         // connecting the nodes to this new manifold
         // DOES NOT WORK because Manifold not in node manifold manager yet:  (*next(nodes.begin(),nit.first)).Assign( *this );
      }

    // sorting branches using the node pointers as keys (default of sort)
    sort( branches_.begin(), branches_.end() );
}




/**
    Creates manifold assigning pointers to the supplied nodes.
    
    @attention the nodes must be assigned to the new Manifold separately once it has been constructed.
*/
template<uint32_t dim>
NodeManifold<dim>::NodeManifold( Node<dim>& inside_node, Node<dim>& outside_node, ManifoldType classifier )
 : parent_geometry_(classifier)
 {
     branches_.push_back( make_pair( &inside_node, INSIDE ) );
     branches_.push_back( make_pair( &outside_node, OUTSIDE ) );
 }




template<uint32_t dim>
NodeManifold<dim>::NodeManifold( const manifold& nodes, ManifoldType geometry )
 : branches_(nodes),
   parent_geometry_(geometry)
{
    assert( branches_.size() >= 2 );
    for( auto& nd : branches_ )
      nd.first->Assign( *this );
}



template<uint32_t dim>
NodeManifold<dim>::NodeManifold( const NodeManifold& nmf )
 : branches_(nmf.branches_),
   parent_geometry_(nmf.parent_geometry_)
 {
 }


template<uint32_t dim>
NodeManifold<dim>::NodeManifold( NodeManifold&& nmf )
 : branches_( move(nmf.branches_) ),
   parent_geometry_( move(nmf.parent_geometry_) )
 {
 }



template<uint32_t dim>
NodeManifold<dim>& NodeManifold<dim>::operator=( const NodeManifold<dim>& nmf )
 {
    if ( this != &nmf ) {
         branches_        = nmf.branches_;
         parent_geometry_ = nmf.parent_geometry_;
      }
    return *this;
 }
 
 

template<uint32_t dim>
NodeManifold<dim>& NodeManifold<dim>::operator=( NodeManifold<dim>&& nmf )
 {
    if ( this != &nmf ) {
         branches_        = move(nmf.branches_);
         parent_geometry_ = move(nmf.parent_geometry_);
      }
    return *this;
 }



// accessor
template<uint32_t dim>
ManifoldType NodeManifold<dim>::GeometricClassifier() const
 {return parent_geometry_;}

// mutator
template<uint32_t dim>
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


template<uint32_t dim>
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
        for (auto i = 0; i < n_branches; ++i) {
            double first = branches_[i].first->Read(index);
            assert(!isnan(first));
            for (size_t j = i + 1; j < n_branches; ++j) {
                double second = branches_[j].first->Read(index);
                assert(!isnan(second));
                if (first > second) swap(branches_[i], branches_[j]);
            }
        }
    // if the variable is associated with the parent elements
    } else if(index.place==ELEMENT) {
        for (auto i = 0; i < branches_.size(); ++i) {
            //double first = numeric_limits<double>::quiet_NaN();
            double first_value(0.);
            size_t first_count(0);
            for( auto e = 0; e < branches_[i].first->Parents(); e++) {
                //if( (dim==2 && nodes_[i]->Parent(e)->IsSurface()) || (dim==3 && nodes_[i]->Parent(e)->IsVolume()) ) {
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

            for ( auto j = i + 1; j < n_branches; ++j) {
                //double second = numeric_limits<double>::quiet_NaN();
                double second_value(0.);
                uint32_t second_count(0);
                for( auto e = 0; e < branches_[j].first->Parents(); e++) {
                    //if( (dim==2 && nodes_[j]->Parent(e)->IsSurface()) || (dim==3 && nodes_[j]->Parent(e)->IsVolume()) ) {
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
                //if (first > second) swap(nodes_[i], nodes_[j]);
                if (first_value > second_value) swap(branches_[i], branches_[j]);
            }
        } 
    } else {
        throw csmp::Exception( INFO, "NodeManifold<dim>::SortByVariableValue",
                                     "variable placement not supported");      

    }
    
} // end SortByVariableValue





template<uint32_t dim>
size_t NodeManifold<dim>::Branches() const
{
  return branches_.size();
}




template<uint32_t dim>
Node<dim>* const NodeManifold<dim>::N( size_t branch ) const
{
  assert( branch < branches_.size() );
  return branches_[branch].first;
}


/// where the node resides
template<uint32_t dim>
INTERFACE_SIDE NodeManifold<dim>::InterFaceSide( size_t branch ) const
 {
   assert( branch < branches_.size() );
   return branches_[branch].second;
 }


template<uint32_t dim>
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
template<uint32_t dim>
vector<Node<dim>*>  NodeManifold<dim>::NodesLocatedAt( INTERFACE_SIDE side ) const
 {
    vector<Node<dim>*> temp;
    temp.reserve( branches_.size() );
    for ( auto& nit : branches_ )
      if ( nit.second == side )
        temp.push_back( nit.first );
    
    return temp;
 }




template<uint32_t dim>
bool NodeManifold<dim>::Add( Node<dim>* nd, INTERFACE_SIDE side )
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    //a node can only exist in one manifold
    if(nd->Manifold() != nullptr) {
        csmp_error.Note( INFO, "NodeManifold<dim>::Add(Node<dim>*)",
                                 "Node already in another manifold. Nothing was done.");         
        return false;
    }

    //the node is already in current manifold
    for ( auto& nit : branches_ )
      if ( nit.first == nd ) {
          csmp_error.Note( INFO, "NodeManifold<dim>::Add(Node<dim>*)",
                                   "Node already part of current manifold. Nothing was done.");
          return false;
        }

    // making sure that the vector does not grow by multiples of 2
    if ( branches_.capacity() == branches_.size() ) branches_.reserve( branches_.size() + 1 );
    branches_.push_back( make_pair(nd,side) );
    nd->Assign( *this );
    
    // reviewing topology of node after insertion to see whether change is necessary
    // Do this after multiple additions: parent_geometry_ = consistencyCheck( *this );

    return true;

}




template<uint32_t dim>
bool NodeManifold<dim>::Remove( const Node<dim>* const nd )
{
   for ( auto& nit : branches_ )
     if ( nit.first == nd ) {
         // disconnecting the node from the manifold
         nit.first = nullptr;
         // removing the Node entry from the branch list
         branches_.erase( remove( branches_.begin(), branches_.end(), nit ) );
         return true;
       }

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.Note( WARNING, "NodeManifold<dim>::Remove(Node<dim>*)",
                                "Node does not exist in manifold. Nothing was done.");
    return false;

}




template<uint32_t dim>
void NodeManifold<dim>::Out() const
{
  cout << "\nNodeManifold: "<< parse(parent_geometry_) <<" with nodes (indices):\t";
  for ( auto i{0U}; i<Branches(); i++ ) {
      if ( NodeManifold<dim>::N(i) )
        cout <<"\n\t"<< NodeManifold<dim>::N(i)->Idx() <<": "<< parseSide( InterFaceSide(i) ) << "\t";
      else
        cerr <<"\nNodeManifold<dim>::Out: branch node pointer is a null pointer.";
    }
  cout << endl;
}

template class NodeManifold<1U>;
template class NodeManifold<2U>;
template class NodeManifold<3U>;


/**
   Check whether the current topologic classification of the NodeManifold (still) makes sense
   
      Idea: dependent on topology there are only so many options
      
      @todo use more diagnostics, e.g., element types (line, surface, volume)
      @todo use also the material IDs that come together at that node
      @todo use initialiser_list<>  or something to create little comparitor functions that make the comparisons more readable
*/
template<uint32_t dim>
ManifoldType  consistencyCheck( const NodeManifold<dim>& nmf  )
 {
    // diagnostics: go over all possible cases excluding the ones that are not possible
    const ManifoldType te       = nmf.GeometricClassifier();
    const auto collocated_nodes = nmf.Branches();
    
    // Manifolds in 2D models
    // ----------------------
    if constexpr ( dim == 2 )
      {
         // anywhere along a split boundary
         if ( collocated_nodes == 2 ) {
             // if the classification is plausible, nothing is done
             if ( te == ManifoldType::INTERFACE ) {
                  bool isBPOINT1{true};
                  for ( auto i{0U}; i<collocated_nodes; ++i ) {
                       if ( nmf.N(i)->AtBoundary() == NOT ||
                            nmf.N(i)->AtBoundary() == INTERNAL ) isBPOINT1 = false;
                       break;
                    }
                  if ( isBPOINT1 ) return ManifoldType::BPOINT1;
                  return te;
               }
             // current classification is probably not correct
             return ManifoldType::NOT_CLASSIFIED;
          }
        // intersections of split boundaries
        if ( collocated_nodes >= 3 ) {
             // special case: an isolated manifold point (POINTX) like for a well
             if ( nmf.N(0)->Neighbors() == 0 ) return ManifoldType::POINTX;
             return ManifoldType::INTERSECTION_POINT;
          }
         // TODO: deal with single-point manifolds
      } // end 2D models

    // Manifolds in 3D models
    // ----------------------
    // extended diagnostics: checking intervening element to determine manifold type
    if constexpr ( dim == 3 )
      {
        // 1. most common case: manifold is located somewhere on a split boundary
        if ( collocated_nodes == 2 ) {
             // classification is plausible
             if ( te == ManifoldType::INTERFACE || te == ManifoldType::BEDGE || te == ManifoldType::BINTERSECTION_POINT ) return te;
             // current classification is probably not correct
             return ManifoldType::NOT_CLASSIFIED;
          }
          
        // 2. if there is a middle node which is only connected to line elements
        if ( collocated_nodes == 2 && nmf.MIDDLE_Node() )
          {
             if ( te == ManifoldType::INTERFACE || te == ManifoldType::BEDGE || te == ManifoldType::BINTERSECTION_POINT ) return te;
             if ( te == ManifoldType::EDGE || te == ManifoldType::END_POINT || te == ManifoldType::INTERSECTION_POINT ) return te;

             const Node<dim>* const mnd_ptr = nmf.MIDDLE_Node();
             bool line_elmt_manifold{true}, surf_elmt_manifold{true};
             const uint32_t connected_elmts = mnd_ptr->Parents();
             
             // in the case of a fracture manifold (we can have an endpoint...)
             for ( auto i{0U}; i<connected_elmts; ++i ) {
                  assert( mnd_ptr->Parent(i) );
                  if ( !mnd_ptr->Parent(i)->IsLine() ) line_elmt_manifold = false;
                  if ( !mnd_ptr->Parent(i)->IsSurface() ) surf_elmt_manifold = false;
               }
             // classification is plausible
             if ( line_elmt_manifold ) {
                  if ( te == ManifoldType::POINT2 && connected_elmts == 2 ) return te;
                  if ( te == ManifoldType::POINT3 && connected_elmts == 3 ) return te;
               }
             // current classification is probably not correct
             return ManifoldType::NOT_CLASSIFIED;
          }
          
        // 3. split boundary T junction
        if ( collocated_nodes == 3 ) {
             if ( te == ManifoldType::INTERSECTION || te == ManifoldType::BINTERSECTION_POINT ||
                  te == ManifoldType::BPOINT2 || te == ManifoldType::BPOINTX ) return te;
             else return ManifoldType::NOT_CLASSIFIED;
          }
        // 4. simple split boundary intersection
        if ( collocated_nodes >= 4 ) {
             if ( te == ManifoldType::INTERSECTION || te == ManifoldType::POINTX || te == ManifoldType::BINTERSECTION_POINT ) return te;
             else return ManifoldType::NOT_CLASSIFIED;
          }
        // 5. multiple intersecting split boundaries
        if ( collocated_nodes >= 5 ) {
             if ( te == ManifoldType::MULTI_INTERSECTION || te == ManifoldType::POINTX || te == ManifoldType::BPOINTX ) return te;
             else return ManifoldType::NOT_CLASSIFIED;
          }
      
      } // end 3D manifolds
      
    // all other cases
    return te;
 }

template ManifoldType  consistencyCheck( const NodeManifold<3>& );
template ManifoldType  consistencyCheck( const NodeManifold<2>& );
template ManifoldType  consistencyCheck( const NodeManifold<1>& );

}// csmp
