#include "NodeManifold.h"
#include "plf_colony.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/**
    Transfers manifold information directly from supplied data to new manifold,
    connecting the nodes from the container in the MeshManager with the newly created manifolds.
    
    @attention the manifold vector does not get sorted!
 */
template<uint32_t dim>
NodeManifold<dim>::NodeManifold( plf::colony<Node<dim> >& nodes, const vector<size_t>& manifold_nodes )
{
    const size_t n_branches(manifold_nodes.size());
    assert( n_branches >= 2 );
    branches_.reserve( n_branches );

    for ( const auto& nit : manifold_nodes ) {
         // connecting the new manifold to its nodes
         branches_.emplace_back( &(*next(nodes.begin(),static_cast<long>(nit)) ) );
         // connecting the nodes back to this new manifold DOES NOT WORK because
         // Manifold is not in node manifold manager yet.
         // To do so after the manifold is properly constructed:  (*next(nodes.begin(),nit.first)).Assign( *this );
      }

    // sorting branches using the node pointers as keys (default of sort)
    // sort( branches_.begin(), branches_.end() );
}




/**
    Creates manifold assigning pointers to the supplied nodes.
    
    @attention the nodes must be assigned to the new Manifold separately once it has been constructed.
*/
template<uint32_t dim>
NodeManifold<dim>::NodeManifold( Node<dim>& inside_node, Node<dim>& outside_node )
 {
     branches_.reserve(2U);
     branches_.push_back( &inside_node );
     branches_.push_back( &outside_node );
 }




template<uint32_t dim>
NodeManifold<dim>::NodeManifold( const manifold& nodes )
 : branches_(nodes)
{
    assert( branches_.size() >= 2 );
    for( auto& nd : branches_ )
      nd->Assign( *this );
}





struct TopoCounts {
    uint8_t interior_surface = 0;
    uint8_t perimeter_surface = 0;
    uint8_t exterior_surface = 0;
    uint8_t interior_line = 0;
    uint8_t perimeter_line = 0;
    uint8_t exterior_line = 0;
    uint8_t interior_point = 0;
    uint8_t perimeter_point = 0;
    uint8_t exterior_point = 0;

    uint8_t total() const noexcept {
        return interior_surface + perimeter_surface + exterior_surface +
               interior_line + perimeter_line + exterior_line +
               interior_point + perimeter_point + exterior_point;
    }
};


template<uint32_t dim>
ManifoldType NodeManifold<dim>::Classify() const noexcept
{
    if constexpr ( dim == 3 )
      {
        const std::size_t n = branches_.size();

        if (n == 0) return ManifoldType::STAND_ALONE;

        TopoCounts c;

        for (const Node<dim>* node : branches_) {
            if (!node) continue;

            switch (node->Attribute()) {
                case INTERIOR_SURFACE:  ++c.interior_surface;  break;
                case PERIMETER_SURFACE: ++c.perimeter_surface; break;
                case EXTERIOR_SURFACE:  ++c.exterior_surface;  break;

                case INTERIOR_LINE:     ++c.interior_line;     break;
                case PERIMETER_LINE:    ++c.perimeter_line;    break;
                case EXTERIOR_LINE:     ++c.exterior_line;     break;

                case INTERIOR_POINT:    ++c.interior_point;    break;
                case PERIMETER_POINT:   ++c.perimeter_point;   break;
                case EXTERIOR_POINT:    ++c.exterior_point;    break;

                case MESH_VERTEX:
                default:
                    // deliberately ignored — does not constrain topology
                    break;
            }
        }

        /* ------------------------------------------------------------
           1. Pure interior / duplicated nodes
           ------------------------------------------------------------ */

        if (c.total() == 0 || n == 1)
            return ManifoldType::STAND_ALONE;

        /* ------------------------------------------------------------
           2. Split-boundary interior (most common)
           ------------------------------------------------------------ */

        if (c.interior_surface == 2 &&
            c.perimeter_line == 0 &&
            c.perimeter_point == 0)
            return ManifoldType::SPLIT_BOUNDARY;

        /* ------------------------------------------------------------
           3. Split-boundary with internal mesh constraint
           ------------------------------------------------------------ */

        if (c.interior_surface >= 2 &&
            (c.interior_line > 0 || c.interior_point > 0))
            return ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH;

        /* ------------------------------------------------------------
           4. Split-boundary crossings
           ------------------------------------------------------------ */

        if (c.interior_surface == 4)
            return ManifoldType::SPLIT_BOUNDARY_CROSSING;

        if (c.interior_surface >= 6)
            return ManifoldType::MULTI_SB_CROSSING;

        /* ------------------------------------------------------------
           5. Split-boundary terminations
           ------------------------------------------------------------ */

        if (c.interior_surface >= 2 &&
            (c.perimeter_line > 0 || c.perimeter_point > 0))
            return ManifoldType::SPLIT_BOUNDARY_END;
            
        /* ------------------------------------------------------------
           Fallback
           ------------------------------------------------------------ */

        return ManifoldType::STAND_ALONE;
    }

  // for 2 dimensional models
  else if constexpr ( dim == 2 )
      {
        const std::size_t n = branches_.size();

        if (n == 0) return ManifoldType::STAND_ALONE;

        TopoCounts c;

        for (const Node<dim>* node : branches_) {
            if (!node) continue;

            switch (node->Attribute()) {
                case INTERIOR_LINE:     ++c.interior_line;     break;
                case PERIMETER_LINE:    ++c.perimeter_line;    break;
                case EXTERIOR_LINE:     ++c.exterior_line;     break;

                case INTERIOR_POINT:    ++c.interior_point;    break;
                case PERIMETER_POINT:   ++c.perimeter_point;   break;
                case EXTERIOR_POINT:    ++c.exterior_point;    break;

                case MESH_VERTEX:
                default:
                    // deliberately ignored — does not constrain topology
                    break;
            }
        }

        /* ------------------------------------------------------------
           1. Pure interior / duplicated nodes
           ------------------------------------------------------------ */

        if (c.total() == 0 || n == 1) return ManifoldType::STAND_ALONE;

        /* ------------------------------------------------------------
           2. Split-boundary interior (most common)
           ------------------------------------------------------------ */

        if (c.interior_line == 2 &&
            c.perimeter_line == 0 &&
            c.perimeter_point == 0)
            return ManifoldType::SPLIT_BOUNDARY;

        /* ------------------------------------------------------------
           3. Split-boundary with internal mesh constraint
           ------------------------------------------------------------ */

        if (c.interior_line >= 2 &&
            (c.interior_line > 0 || c.interior_point > 0))
            return ManifoldType::SPLIT_BOUNDARY_WITH_INTERNAL_MESH;

        /* ------------------------------------------------------------
           4. Split-boundary crossings
           ------------------------------------------------------------ */

        if (c.interior_line == 4)
            return ManifoldType::SPLIT_BOUNDARY_CROSSING;

        if (c.interior_line >= 6)
            return ManifoldType::MULTI_SB_CROSSING;

        /* ------------------------------------------------------------
           5. Split-boundary terminations
           ------------------------------------------------------------ */

        if (c.interior_line >= 2 &&
            (c.perimeter_line > 0 || c.perimeter_point > 0))
            return ManifoldType::SPLIT_BOUNDARY_END;

        /* ------------------------------------------------------------
           Fallback
           ------------------------------------------------------------ */

        return ManifoldType::STAND_ALONE;
    }

    // 1D models
     return ManifoldType::SPLIT_BOUNDARY_END;
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
        for ( size_t i = 0; i < n_branches; ++i) {
            double first = branches_[i]->Read(index);
            assert(!isnan(first));
            for (size_t j = i + 1; j < n_branches; ++j) {
                double second = branches_[j]->Read(index);
                assert(!isnan(second));
                if (first > second) swap(branches_[i], branches_[j]);
            }
        }
    // if the variable is associated with the parent elements
    } else if(index.place==ELEMENT) {
        for ( size_t i = 0; i < branches_.size(); ++i) {
            //double first = numeric_limits<double>::quiet_NaN();
            double first_value(0.);
            size_t first_count(0);
            for( uint32_t e = 0; e < branches_[i]->Parents(); e++) {
                //if( (dim==2 && nodes_[i]->Parent(e)->IsSurface()) || (dim==3 && nodes_[i]->Parent(e)->IsVolume()) ) {
                if ( !isnan(branches_[i]->Parent(e)->Read(index)) ) {
                    first_value += branches_[i]->Parent(e)->Read(index);
                    first_count ++;
                }
            }
            //assert(!isnan(first));
            assert(first_count != 0);
            first_value /= first_count;
            assert(!isnan(first_value));

            for ( size_t j = i + 1; j < n_branches; ++j) {
                //double second = numeric_limits<double>::quiet_NaN();
                double second_value(0.);
                uint32_t second_count(0);
                for( uint32_t e = 0; e < branches_[j]->Parents(); e++) {
                    //if( (dim==2 && nodes_[j]->Parent(e)->IsSurface()) || (dim==3 && nodes_[j]->Parent(e)->IsVolume()) ) {
                    if( !isnan(branches_[j]->Parent(e)->Read(index)) ) {
                        second_value += branches_[j]->Parent(e)->Read(index);
                        second_count ++;
                    }
                }
                //assert(!isnan(second));
                assert(second_count != 0);
                second_value /= second_count;
                assert(!isnan(second_value));                
                if (first_value > second_value) swap(branches_[i], branches_[j]);
            }
        } 
    } else {
        throw csmp::Exception( INFO, "NodeManifold<dim>::SortByVariableValue",
                                     "variable placement not supported");      

    }
    
} // end SortByVariableValue




///Clears existing InterFaces of assigned to node and assigns a new set to them.
/// @attention This overwrites existing interfaces assigned to node if they are assigned.
template<uint32_t dim>
void NodeManifold<dim>::Assign(Node<dim>* n, std::set<std::pair<InterFace<dim>*,std::pair<uint32_t,INTERFACE_SIDE>>> interface_indexes )
{
  assert(!interface_indexes.empty());

  auto nmap_it = node_parent_interface_map_.find(n);
  //If new node entry
  if ( nmap_it==node_parent_interface_map_.end()){
    //create interface_index_vector from set
    vector<pair<InterFace<dim>*,pair<uint32_t,INTERFACE_SIDE>>> if_idx_vec;
    for (auto& pair:interface_indexes){
      if_idx_vec.push_back(pair);
    }
    //insert into map
    node_parent_interface_map_.insert(make_pair(n,if_idx_vec));
    return;
  } else {
    auto if_indexes = nmap_it->second;
    //if vector is already assigned - we overwrite
    if (!if_indexes.empty())
      if_indexes.clear(); //clear current vector

    for ( auto pair : interface_indexes)
      if_indexes.push_back(pair);           //insert interface index pair into map

    return;
  }
}//end of Assign







template<uint32_t dim>
uint32_t NodeManifold<dim>::Branches() const  noexcept
{
  return static_cast<uint32_t>(branches_.size());
}


template<uint32_t dim>
uint32_t NodeManifold<dim>::NodeMapSize() const
{
  return static_cast<uint32_t>(node_parent_interface_map_.size());
}




template<uint32_t dim>
uint32_t NodeManifold<dim>::InterFaces( Node<dim>* const n ) const noexcept
{
  //check map has been calibrated
  assert( !node_parent_interface_map_.empty() );

  return static_cast<uint32_t>(node_parent_interface_map_.at(n).size());
}


template<uint32_t dim>
Node<dim>* const NodeManifold<dim>::N( size_t branch ) const noexcept
{
  assert( branch < branches_.size() );
  return branches_[branch];
}



template<uint32_t dim>
InterFace<dim>* NodeManifold<dim>::I( Node<dim>* const n, uint32_t i )
{
  assert(!node_parent_interface_map_.empty());

  auto interface_indexes = node_parent_interface_map_.at(n);
  assert(i < interface_indexes.size()); //out of scope
  return interface_indexes[i].first;
}


template<uint32_t dim>
std::pair<InterFace<dim>*, std::pair<uint32_t,INTERFACE_SIDE>> NodeManifold<dim>::InterFaceIndex( Node<dim>* const n, uint32_t i )
{
  assert(!node_parent_interface_map_.empty());
  auto interface_indexes = node_parent_interface_map_.at(n);
  assert(i < interface_indexes.size()); //out of scope

  return interface_indexes[i];
}


template<uint32_t dim>
std::vector<std::pair<InterFace<dim>*, std::pair<uint32_t,INTERFACE_SIDE>>> NodeManifold<dim>::InterFaceIndexVector( Node<dim>* const n )
{
  assert(!node_parent_interface_map_.empty());
  return  node_parent_interface_map_.at(n);
}




template<uint32_t dim>
bool NodeManifold<dim>::Add( Node<dim>* nd )
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
      if ( nit == nd ) {
          csmp_error.Note( INFO, "NodeManifold<dim>::Add(Node<dim>*)",
                                   "Node already part of current manifold. Nothing was done.");
          return false;
        }

    // making sure that the vector does not grow by multiples of 2
    if ( branches_.capacity() == branches_.size() ) branches_.reserve( branches_.size() + 1 );
    branches_.push_back( nd );
    nd->Assign( *this );
    
    // reviewing topology of node after insertion to see whether change is necessary
    // Do this after multiple additions: parent_geometry_ = consistencyCheck( *this );

    return true;

}




template<uint32_t dim>
bool NodeManifold<dim>::Remove( const Node<dim>* const nd )
{
   // in this is true, this should not be a manifold anymore
   if ( branches_.empty() ) return false;
   
   for ( auto& nit : branches_ )
     if ( nit == nd ) {
         // disconnecting the node from the manifold
         nit = nullptr;
         // removing the Node entry from the branch list
         branches_.erase( remove( branches_.begin(), branches_.end(), nit ) );
         return true;
       }

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.Note( WARNING, "NodeManifold<dim>::Remove(Node<dim>*)",
                                "Node does not exist in manifold. Nothing was done.");
    return false;
}



      /// checks whether all nodes in the  manifold have the same location using operator< of point
template<uint32_t dim>
bool NodeManifold<dim>::AreNodesCollocated() const noexcept
 {
    const Point<dim> pt = branches_[0]->Coordinate();
    for ( const auto& nit : branches_ )
      if ( pt != nit->Coordinate() )
        return false;
        
    return true;
 }


/**
    Returns manifold node indices and type to data structure used to initialise VData
*/
template<uint32_t dim>
pair<vector<size_t>,ManifoldType>  NodeManifold<dim>::Data() const noexcept
 {
    // there must be at least 2 nodes
    vector<size_t> node_ids{ branches_[0]->Idx(), branches_[1]->Idx() };
    if ( branches_.size() > 2 ) {
         node_ids.reserve( branches_.size() );
         for ( auto i{2U}; i<branches_.size(); i++ )
           node_ids.push_back( branches_[i]->Idx() );
      }
    return make_pair( node_ids, Classify() );
 }




template<uint32_t dim>
void NodeManifold<dim>::Out() const
{
  cout << "\nNodeManifold: "<< parse(Classify()) <<" with nodes (indices):\t";
  for ( auto i{0U}; i<Branches(); i++ ) {
      if ( NodeManifold<dim>::N(i) )
        cout <<"\n\t"<< NodeManifold<dim>::N(i)->Idx() <<": ";
      else
        cerr <<"\nNodeManifold<dim>::Out: branch node pointer is a null pointer.";
    }
  cout << endl;
}

template class NodeManifold<1U>;
template class NodeManifold<2U>;
template class NodeManifold<3U>;



/**
   Checks whether the current topologic classification of the combination of nodes within one NodeManifold makes sense,
   using the number of collocated nodes and their geometric attributes as evaluation criteria.
   Assigns classification that adds SplitBoundary specific features.

   What can be distinguished here? - To start with the nodes have the attributes, see CSMP_global_enumerations.h:
   
    VERTEX_POINT,      ///< a point within the model volume
    INTERSECTION_POINT,///< a point where lines cross or multiple surfaces intersect
    PERIMETER_POINT,   ///< point at the end of a line inside a 2D model
    EXTERIOR_POINT,    ///< on an outside surface of the model
    INTERIOR_LINE,     ///< a line on the interior of the model
    PERIMETER_LINE,    ///< on a surface edge inside of the model
    EXTERIOR_LINE,     ///< an edge of the model
    INTERSECTION_LINE, ///<  belonging to multiple surfaces
    INTERIOR_SURFACE,  ///< a surface withing the model
    PERIMETER_SURFACE, ///< a surface forming the hull of an object inside of the model
    EXTERIOR_SURFACE  ///< a surface delimiting the model
 
    In which way does a NodeManifoldType go beyond these classifications?
    
         - classification that applies to the group of nodes united in the manifold
         - indication of a SplitBoundary (but this is already done by having a node manifold)
         - specifics lie the distinction of a point on a line ending in 3D space from a perimeter point om a line boundary in 2D
         - intersection between (multiple) split boundaries
         - intersection between boundary and SplitBoundary
         - intersection points between a line and a splitboundary LINE_SPLITBOUNDARY_INTERSECTION
    
    What are the possible combinations of BREP TOPOTYPES distinguishing different types of manifolds?
    
   0. collocated point (2D, not a splitboundary), and line splitboundaries (3D) = as many nodes as parent elements of the line-element node
   1. crossing points -> 4 nodes for 2 lines in 2D & 3D, 5 nodes for 3, 6 nodes for 4...
   2. surface intersections -> only in 3D: 4 nodes  where 2 splitboundaries cross
   3. three nodes -> T intersection in 2D
   
   Consistency Checks performed, using set<TOPOTYPE>:
   ==========================================
   
   1. If all the nodes in the set have the same TOPOTYPE, consistency is achieved and the ManifoldType is determined
   
   2. If the nodes have different types, each node is checked to see whether the TOPOTYPE clashes with the actual geometry that it forms part of. For the potentially corrected classifications a re-evaluation of the ManifoldType is performed.
   
   3. SplitBoundary specific criteria are considered to distinguish basic intersections from boundary intersections from splitboundary intersections
   
   @test refactored 9/07/2022 after removal of the side information
   
*/
template<uint32_t dim>
ManifoldType  consistencyCheck( const NodeManifold<dim>& nmf, bool verbose  )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    // 0. getting the TOPOTYPES of the nodes
    const auto collocated_nodes = nmf.Branches();
    set<TOPOTYPE> node_attributes;
    for ( uint32_t i{0U}; i<collocated_nodes; i++ )
      node_attributes.insert( nmf.N(i)->Attribute() );

    // 1. diagnostics: if there is only one topotype
    // ---------------------------------------------
    if ( node_attributes.size() == 1U ) {
         if ( verbose ) {
           switch( (*node_attributes.begin()) ) {
                case MESH_VERTEX: {
                     if ( nmf.Classify() != ManifoldType::STAND_ALONE ) {
                          cerr <<"\n\t"<< parse(nmf.Classify()) <<" vs. 'STAND_ALONE'";
                          csmp_error.Note( WARNING, "consistencyCheck", "manifold type possible incorrect, resetting." );
                       }
                    }
                  break;
                case INTERIOR_POINT: {
                     if ( nmf.Classify() != ManifoldType::SPLIT_BOUNDARY_CROSSING ) {
                          cerr <<"\n\t"<< parse(nmf.Classify()) <<" vs. 'SPLIT_BOUNDARY_CROSSING'";
                          csmp_error.Note( WARNING, "consistencyCheck", "manifold type possible incorrect, resetting." );
                       }
                    }
                  break;
                case EXTERIOR_LINE:
                case PERIMETER_POINT:
                case EXTERIOR_POINT: {
                     if ( nmf.Classify() != ManifoldType::SPLIT_BOUNDARY_END ) {
                          cerr <<"\n\t"<< parse(nmf.Classify()) <<" vs. 'SPLIT_BOUNDARY_END'";
                          csmp_error.Note( WARNING, "consistencyCheck", "manifold type possible incorrect, resetting." );
                       }
                    }
                  break;
                case INTERIOR_SURFACE:
                case PERIMETER_SURFACE:
                case INTERIOR_LINE: {
                     if ( nmf.Classify() != ManifoldType::SPLIT_BOUNDARY ) {
                          cerr <<"\n\t"<< parse(nmf.Classify()) <<" vs. 'SPLIT_BOUNDARY'";
                          csmp_error.Note( WARNING, "consistencyCheck", "manifold type possible incorrect, resetting." );
                       }
                    }
                  break;
                case EXTERIOR_SURFACE:
                     csmp_error.Note( ERROR, "consistencyCheck", "NodeManifolds cannot contain Nodes classified as EXTERIOR_SURFACE." );
                  break;
                default:
                  cerr <<"\nconsistencyCheck(NodeManifold): Node TOPOTYPE not resolved."<< endl;
                  return nmf.Classify();
             }
           }
      }

    
    // 2. if all diagnostics fails the originally assigned qualifier is returned
    return nmf.Classify();

 } // end consistency check

template ManifoldType  consistencyCheck( const NodeManifold<3>&, bool );
template ManifoldType  consistencyCheck( const NodeManifold<2>&, bool );
template ManifoldType  consistencyCheck( const NodeManifold<1>&, bool );

}// csmp
