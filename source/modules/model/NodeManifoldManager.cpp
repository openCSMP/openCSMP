// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "NodeManifoldManager.h"
#include "NodeManifold.h"
#include "ErrorHandler.h"
#include "Node.h"
#include "binaryReadWrite.h"


using namespace std;

namespace csmp {

/**
    Preferred constructor of  node manifolds, including the case where these are read back from a CSMP native binary fileset
    
    @attention do not use this constructor after nodes were deleted from the MeshManager.
    
    @todo make sure that the map:key nodes are ideed contained in the manifolds.
*/
template<uint32_t dim>
NodeManifoldManager<dim>::NodeManifoldManager( plf::colony<Node<dim>>& mesh_nodes,
                                               VData::manifoldContainer::const_iterator first,
                                               VData::manifoldContainer::const_iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( distance(first,last) == 0 ) {
        csmp_error.Note( WARNING, "NodeManifoldManager::constructor:",
                          "no manifold information contained in vertexManifoldIndices map; no manifolds were constructed" );
        return;
      }
    if ( mesh_nodes.empty() ) {
        csmp_error.Note( ERROR, "NodeManifoldManager::constructor:",
                          "no nodes available to create manifolds from" );
        return;
      }
    // checking that the deque does indeed have ther required node entries
    if ( distance(first,last) > static_cast<long>(mesh_nodes.size()) )
      csmp_error.Note( WARNING, "NodeManifoldManager::constructor:",
                        "it appears that more manifold indices are supplied than nodes" );

    // creating the node manifolds
    while ( first != last ) {
        // checking the manifold construction data
        const size_t  n_manifold_nodes((*first).first.size());
        // a manifold requires at least two nodes
        assert( n_manifold_nodes >= 2 );
        if ( n_manifold_nodes >= 2 ) {
            typename plf::colony<NodeManifold<dim>>::iterator mit =
              node_manifolds_.emplace( NodeManifold<dim>( mesh_nodes, (*first).first ) );
            // applying a consistency check
            if ( (*first).second != (*mit).Classify() ) {
                 cout <<"\n\t"<< parse((*first).second) <<" vs "<< parse((*mit).Classify()) << endl;
                 csmp_error.Note( WARNING, "NodeManifoldManager::constructor:",
                                 "manifold classification stored in input VSet does not match the new manifolds own diagnostic" );
              }
            // assigning the new NodeManifold to the nodes it points to
            for ( uint32_t i{0U}; i<(*mit).Branches(); i++ )
              (*mit).N(i)->Assign( (*mit) );
          }
         first++;
      }
      
    cout <<"\nNodeManifoldManager(custom ctor): constructed "<< node_manifolds_.size();
    cout <<" node manifolds from the input data.\n";
    // Out();
    
 } // end custom constructor



template<uint32_t dim>
typename NodeManifoldManager<dim>::manifoldIterator NodeManifoldManager<dim>::ManifoldsBegin()
  { return node_manifolds_.begin(); }


template<uint32_t dim>
typename NodeManifoldManager<dim>::manifoldIterator NodeManifoldManager<dim>::ManifoldsEnd()
  { return node_manifolds_.end(); }

template<uint32_t dim>
typename NodeManifoldManager<dim>::manifoldConstIterator NodeManifoldManager<dim>::ManifoldsBegin() const
  { return node_manifolds_.begin(); }


template<uint32_t dim>
typename NodeManifoldManager<dim>::manifoldConstIterator NodeManifoldManager<dim>::ManifoldsEnd() const
  { return node_manifolds_.end(); }


template<uint32_t dim>
size_t  NodeManifoldManager<dim>::Manifolds() const
  { return node_manifolds_.size(); }




/**
    Creates a node manifold from the 2 pointers to co-located nodes.
    Additional nodes can be added to it later.
    
    @param nodes reference to the primary storage of the nodes in the MeshManager
    @param inside pointer to the node that will assumed to be on the INSIDE of the interface that that manifold is part of
    @param outside opposing yet collocated node in the manifold
    @return an iterator to the new manifold stored in a colony in the NodeManifoldManager.
    
    @attention this method does not assign the nodes to the new manifold.; this must be done by calling  Node::Assign( manifold ) after its successful creation
    which only happens on the return from this method.
 */
template<uint32_t dim>
typename plf::colony<NodeManifold<dim> >::iterator NodeManifoldManager<dim>::AddManifold( plf::colony<Node<dim> >& nodes,
                                                                                          Node<dim>* const inside,
                                                                                          Node<dim>* const outside )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // both nodes must be manifolds
     if ( inside == nullptr || outside == nullptr ) {
          csmp_error.Note( ERROR, "NodeManifoldManager<dim>::AddManifold",
                            "input parameters contain nullptr nodes; nothing can be done" );
          return node_manifolds_.end();
       }
       
    // constructing the manifold supplying references to the actual Node objects stored in the colony via iterators
    auto inside_it  = nodes.get_iterator( inside );
    auto outside_it = nodes.get_iterator( outside );
    
    // return node_manifolds_.insert( NodeManifold( (*inside_it), (*outside_it), manifold_type ) );
    return node_manifolds_.insert( NodeManifold( (*inside_it), (*outside_it) ) );

 } // end AddManifold




/**
    Tries to replace the two manifolds by a single one that connects all of  their nodes;
    succeeds if the two share nodes. Will return true in this case;
    Fixes all node connections.
    
    @param nmf1  if the operation is successful, the enlarged manifold is returned into the first manifold pointer, else both manifolds remain untouched.
    @param nmf2  if the operation is successful, the enlarged manifold is returned into the first manifold pointer, else both manifolds remain untouched.

   @return true if it was possible to merge the manifolds because they did share nodes.
    
*/
template<uint32_t dim>
bool NodeManifoldManager<dim>::MergeManifolds( NodeManifold<dim>* nmf1, NodeManifold<dim>* nmf2 )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // both nodes must be manifolds
     if ( nmf1 == nullptr || nmf2 == nullptr ) {
          csmp_error.Note( ERROR, "NodeManifoldManager<dim>::MergeManifolds",
                            "input contains nullptr manifolds; nothing can be done" );
          return false;
       }
    
    // 1. do the manifolds share nodes
    // -------------------------------
    // (the set of their node pointers must be smaller than the sum of their branches)
    set<const Node<dim>*>  connected_nodes;
    const auto n_nodes_nmf1{ nmf1->Branches() };
    for ( uint32_t i{0U}; i<n_nodes_nmf1; ++i ) {
         assert( nmf1->N(i) != nullptr );
         connected_nodes.insert( nmf1->N(i) );
      }
    const auto n_nodes_nmf2{ nmf2->Branches() };
    for ( uint32_t i{0U}; i<n_nodes_nmf2; ++i ) {
         assert( nmf2->N(i) != nullptr );
         connected_nodes.insert( nmf2->N(i) );
      }
    // if there are no shared nodes, merging is not possible
    if ( (n_nodes_nmf1 + n_nodes_nmf2) >=  connected_nodes.size() ) return false;
    
    // 2. merging the manifolds into nmf1, deleting nmf2
    // -------------------------------------------------
    for ( auto i{0U}; i<n_nodes_nmf2; ++i )
      nmf1->Add( nmf2->N(i) );
     
    // 3. reclassifying the manifold geometry
    // --------------------------------------
    // (this gets done when that classification is needed)
   
    // 4. deleting the merged manifold 2
    // ---------------------------------
    node_manifolds_.erase( node_manifolds_.get_iterator(nmf2) );

    return true;
     
 } // end MergeManifolds







    /// puts the nodes inside of the manifolds into the ascending order of values of the user specified  variable
template<uint32_t dim>
void NodeManifoldManager<dim>::SortManifoldsByVariableValue( std::string var_name, const csmp::Index& var_index )
 {
    for ( auto& nmf : node_manifolds_ ) {
         nmf.SortByVariableValue( var_index );
      }
    current_sort_variable_ = var_name;
   
 } // end SortManifoldsByVariableValue

 
 
 



template<uint32_t dim>
void NodeManifoldManager<dim>::Delete( NodeManifold<dim>* const md )
{
   node_manifolds_.erase( node_manifolds_.get_iterator(md) );
}




/**
  If there are single-node manifolds, their nodes are deconnected and the manifolds are deleted.
  Later, the manager is searched for the remaining null pointers and these are removed.
*/
template<uint32_t dim>
size_t NodeManifoldManager<dim>::DeleteSingleNodeManifolds()
  {
     size_t n_single_node_manifolds(0U);
     for ( auto nit=node_manifolds_.begin(); nit!=node_manifolds_.end(); ++nit )
       if ( (*nit).Branches() < 2 ) {
            // disconnecting the remaining node
            (*nit).Remove( (*nit).N(0) );
            // deleting the manifold
            nit = node_manifolds_.erase( nit );
            n_single_node_manifolds++;
         }
     
      return n_single_node_manifolds;
  }




/**
    writes manifolds to binary file; relies on unique node indexes
    
    @attention call DeleteSingleNodeManifolds() before to avoid storing any defunct manifolds
    @attention make sure that the nodes have unique numbers that match those in the other binary files
*/
template<uint32_t dim>
void NodeManifoldManager<dim>::OutputNodeManifoldsToBinary( const char* file_name,
                                                            const char* node_sorting_variable ) const
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( node_manifolds_.empty() ) {
         csmp_error.Note( ERROR, "NodeManifoldManager::OutputNodeManifoldsToBinary:",
                           "no node manifolds found; call this method only if there are SplitBoundaries in the model" );
         return;
      }

    std::string bin_file( file_name );
    std::fstream fp( bin_file.c_str(), std::ios::out | std::ios::binary );
    if ( !fp.is_open() ) {
        csmp_error.Note( ERROR, "NodeManifoldManager::OutputNodeManifoldsToBinary:",
                           bin_file, "file could not be opened; nothing was done." );
        return;
      }

    // -----------------------------
    // 1. File header
    // -----------------------------
    {
      BinaryFileSectionWrite hdr( fp, "NDMFHEDR" ); // node manifold header

      std::string heading( "NodeManifoldManager::OutputNodeManifoldsToBinary: " );
      heading += "topologically collocated nodes in Model interconnected by manifolds";
      heading += " written to file: ";
      heading += bin_file;
      heading += "'.";

      binaryFileWrite( fp, heading.c_str() );
    }

    // -----------------------------
    // 2. writing the manifolds
    // -----------------------------
    {
      BinaryFileSectionWrite hdr( fp, "MANIFLDS" );

      // writing their number
      size_t records = this->Manifolds();
      fp.write( reinterpret_cast<const char*>(&records), sizeof( size_t ) );

      // 1. record of topologic entity classifiers
      {
        vector<int8_t>  topo_classifiers;
        topo_classifiers.reserve( records );
        for ( auto& nmf : node_manifolds_ )
          topo_classifiers.push_back( static_cast<int8_t>(nmf.Classify()) );
        binaryFileWrite( fp, topo_classifiers );
      }

      // 2. record of nodes per manifold (size, values)
      size_t n_manifold_node_entries(0U);
      {
        vector<int8_t>  nodes_per_manifold;
        nodes_per_manifold.reserve( records );
        for ( auto& nmf : node_manifolds_ ) {
             nodes_per_manifold.push_back( static_cast<int8_t>(nmf.Branches()) );
             n_manifold_node_entries += nmf.Branches();
          }
        binaryFileWrite( fp, nodes_per_manifold );
      }

      // 3. 'plist'-like record of node indices
      //   (the entries will have been sorted by variable values)
      {
        // name of variable that was used for sorting the nodes
        const string sort_variable{ node_sorting_variable };
        binaryFileWrite( fp, sort_variable );
        // sorted nodes
        vector<size_t>  manifold_node_list;
        manifold_node_list.reserve( n_manifold_node_entries );
        // 'plist' like record of nodes per manifold
        for ( auto& nmf : node_manifolds_ ) {
             const size_t entries(nmf.Branches());
             for ( size_t i{0U}; i<entries; ++i )
               manifold_node_list.push_back( nmf.N(i)->Idx() );
          }
        binaryFileWrite( fp, manifold_node_list );
      }
  } // end writing the manifold records
  
} // end OutputNodeManifoldsToBinary
 
 
 
 
 

/** Reads manifolds from binary file, using indices to create pointer connections
 */
template<uint32_t dim>
string NodeManifoldManager<dim>::InputNodeManifoldsFromBinary( plf::colony<Node<dim>>& mesh_nodes,
                                                               const char* file_name )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    std::string bin_file( file_name );
    std::fstream fp( bin_file.c_str(), std::ios::out | std::ios::binary );
    if ( !fp.is_open() ) {
        csmp_error.Note( ERROR, "NodeManifoldManager::InputNodeManifoldsFromBinary:",
                           bin_file, "file could not be opened; nothing was done." );
        return "unknown";
      }

    // -----------------------------
    // 1. File header
    // -----------------------------
    {
      BinaryFileSectionRead hdr( fp, "NDMFHEDR" ); // node manifold header

      std::string heading;
      binaryFileRead( fp, heading );
      cout <<"\nNodeManifoldManager<dim>::InputNodeManifoldsFromBinary: reading: ";
      cout << heading << endl;
    }

    // -----------------------------
    // 2. read the manifolds
    // -----------------------------
    BinaryFileSectionRead hdr( fp, "MANIFLDS" );

    // reading their number
    size_t manifolds;
    fp.read( reinterpret_cast<char*>(&manifolds), sizeof( size_t ) );
    assert( manifolds > 0 ); // this file should only be written if there are manifolds

    // 1. record of topologic entity classifiers
    vector<ManifoldType>  manifold_topology; // int8_t
    binaryFileRead( fp, manifold_topology );

    // 2. record of nodes per manifold (size, values)
    vector<int8_t>  nodes_per_manifold;
    binaryFileRead( fp, nodes_per_manifold );

    // 3. 'plist'-like record of node indices
    //   (the entries will have been sorted by variable values)
    // name of variable that was used for sorting the nodes
    binaryFileRead( fp, current_sort_variable_ );
    // sorted nodes
    vector<uint32_t>  nodes_of_manifolds;
    binaryFileRead( fp, nodes_of_manifolds );

    // ------------------------------------
    // 3. reconstructing the node manifolds 
    // ------------------------------------
    vector<size_t>  manifold_nodes;
    size_t          counter(0U);
    for ( uint32_t i{0}; i<manifolds; ++i ) {
         const size_t n_branches = nodes_per_manifold[i];
         manifold_nodes.reserve( n_branches );
         for ( uint32_t j{0U}; j<n_branches; ++j ) {
              assert( nodes_of_manifolds[counter] < mesh_nodes.size() );
              manifold_nodes.push_back( nodes_of_manifolds[counter] );
              counter++;
           }
         sort( manifold_nodes.begin(), manifold_nodes.end() );
         manifold_nodes.erase( unique( manifold_nodes.begin(), manifold_nodes.end() ), manifold_nodes.end() );
         assert( manifold_topology[i] < ManifoldType::SPLIT_BOUNDARY_END );
         node_manifolds_.emplace( NodeManifold<dim>( mesh_nodes, manifold_nodes ) );
      }
    
    return current_sort_variable_;
    
 } // end InputNodeManifoldsFromBinary
    
    
    
    

template<uint32_t dim>
void NodeManifoldManager<dim>::Out() const
 {
    if ( node_manifolds_.empty() )
      cout <<"\n\nNodeManifoldManager<"<< dim <<">::Out: no manifolds in store.\n";
    else {
         cout <<"\n\nNodeManifoldManager<"<< dim <<">::Out: "<< node_manifolds_.size() <<" manifolds in store.\n";
         for ( const auto& nit : node_manifolds_ ) {
              nit.Out();
           }
      }
 } // end Out



template class NodeManifoldManager<1U>;
template class NodeManifoldManager<2U>;
template class NodeManifoldManager<3U>;

} // csmp
