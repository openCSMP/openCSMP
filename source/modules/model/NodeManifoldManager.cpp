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
NodeManifoldManager<dim>::NodeManifoldManager( const vertexManifoldIndices& indices,
                                               plf::colony<Node<dim>>& mesh_nodes )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( indices.empty() ) {
        csmp_error.notice( WARNING, "NodeManifoldManager::constructor:",
                          "no manifold information contained in vertexManifoldIndices map; no manifolds were constructed" );
        return;
      }
    if ( mesh_nodes.empty() ) {
        csmp_error.notice( ERROR, "NodeManifoldManager::constructor:",
                          "no nodes available to create manifolds from" );
        return;
      }
    // checking that the deque does indeed have ther required node entries
    if ( indices.size() > mesh_nodes.size() )
      csmp_error.notice( WARNING, "NodeManifoldManager::constructor:",
                        "it appears that more manifold indices are supplied than nodes" );

#ifdef DEBUG
// checking that the key nodes in vertexManifoldIndices map are also contained in the corresponding sets
for ( auto& nit : indices ) {
    set<size_t> mnodes;
    for ( const auto& mf_nodes : nit.second ) mnodes.insert( mf_nodes.first );
    // if the key node is not contained this is reported
    if ( mnodes.find(nit.first) == mnodes.end() )
        csmp_error.notice( WARNING, "NodeManifoldManager::constructor:",
                          "manifold does not contain key node: ", to_string(nit.first) );
  }
#endif

    // creating the node manifolds
    // NB: node_id, connected nodes and their INTERFACE_SIDE identifiers
    //     map<size_t,set<pair<size_t,int8_t> > >
    for ( const auto& nit : indices ) {
        const size_t  n_nodes(nit.second.size());
        // a manifold requires at least two nodes
        if ( n_nodes >= 2U ) {  // indices: map<size_t,set<pair<size_t,int8_t> > >
            typename plf::colony<NodeManifold<dim>>::iterator mit =
              node_manifolds_.emplace( NodeManifold<dim>( mesh_nodes, nit.second, ManifoldType::INTERFACE ) );
            // geometric qualifier is determined through a consistency check once the manifold is in place
            (*mit).GeometricClassifier( consistencyCheck( (*mit) ) );
            // assigning the new NodeManifold to the nodes it points to
            for ( auto i{0}; i<(*mit).Branches(); i++ )
              (*mit).N(i)->Assign( (*mit) );
          }
      }
      
    cout <<"\nNodeManifoldManager(custom ctor): constructed "<< node_manifolds_.size();
    cout <<" node manifolds from the input data.\n";
    // Out();
    
 } // end custom constructor






// disconnecting and deleting the node manifolds
template<uint32_t dim>
NodeManifoldManager<dim>::~NodeManifoldManager<dim>()
 {
 }

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
    Tries to replace the two manifolds by a single one that connects all of  their nodes;
    succeeds if the two share nodes. Will return true in this case;
    Fixes all node connections.
    
    @return true if it was possible to merge the manifolds because they did share nodes.
    
    @param nmf1  if the operation is successful, the enlarged manifold is returned into the first manifold pointer, else both manifolds remain untouched.
*/
template<uint32_t dim>
bool NodeManifoldManager<dim>::MergeManifolds( NodeManifold<dim>* nmf1, NodeManifold<dim>* nmf2 )
 {
    assert( nmf1 != nullptr );
    assert( nmf2 != nullptr );
    
    // 1. do the manifolds share nodes
    // -------------------------------
    // (the set of their node pointers must be smaller than the sum of their branches)
    //set<const Node<dim>* const>  connected_nodes;
    set<const Node<dim>*>  connected_nodes;
    const auto n_nodes_nmf1{ nmf1->Branches() };
    for ( auto i{0}; i<n_nodes_nmf1; ++i ) {
         assert( nmf1->N(i) != nullptr );
         connected_nodes.insert( nmf1->N(i) );
      }
    const auto n_nodes_nmf2{ nmf2->Branches() };
    for ( auto i{0}; i<n_nodes_nmf2; ++i ) {
         assert( nmf2->N(i) != nullptr );
         connected_nodes.insert( nmf2->N(i) );
      }
    // if there are no shared nodes, merging is not possible
    if ( (n_nodes_nmf1 + n_nodes_nmf2) >=  connected_nodes.size() ) return false;
    
    // 2. merging the manifolds into nmf1, deleting nmf2
    // -------------------------------------------------
    for ( auto i{0}; i<n_nodes_nmf2; ++i )
      nmf1->Add( nmf2->N(i), nmf2->InterFaceSide(i) );
     
     // 3. reclassifying the manifold geometry
     // --------------------------------------
     nmf1->GeometricClassifier( consistencyCheck( *nmf1 ) );
     
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

 
 
/**
   replace two separate node manifolds by a single one that contains the union of their nodes
   
   @return bool: if the manifolds share nodes, they can indeed be merged, else they are kept separate.
   
   @attention the new Manifold does not get sorted
*/
/* REFACTOR
template<uint32_t dim>
bool NodeManifoldManager<dim>::MergeManifolds( NodeManifold<dim>* mnf1, NodeManifold<dim>* mnf2 )
  {
     assert( mnf1 != nullptr );
     assert( mnf2 != nullptr );
     
     // 1. making sure that the two manifolds actually share nodes
     size_t sum_nodes = mnf1->Branches() + mnf2->Branches();
     vector<pair<Node<dim>*,INTERFACE_SIDE> > combined_manifolds;
     combined_manifolds.reserve(sum_nodes);
     for ( auto i{0}; i<mnf1->Branches(); ++i )
       combined_manifolds.push_back( make_pair( mnf1->N(i), mnf1->InterFaceSide(i) ) );

     // making the vector unique
     sort( combined_manifolds.begin(), combined_manifolds.end() );
     combined_manifolds.erase( unique(combined_manifolds.begin(), combined_manifolds.end()), combined_manifolds.end() );

     // checking - if there are no shared nodes, the two manifolds cannot be merged
     if ( combined_manifolds.size() == sum_nodes ) return false;
     
     // 2. creating a new temp manifold that gets assigned to mnf1
     vector<Node<dim>*>     nodes; nodes.reserve( combined_manifolds.size() );
     vector<INTERFACE_SIDE> sides; sides.reserve( combined_manifolds.size() );
     for ( auto& it : combined_manifolds ) {
          nodes.push_back( it.first );
          sides.push_back( it.second );
       }
       
     NodeManifold<dim> merged_manifold( nodes, sides, ManifoldType::INTERFACE );
     // making sure it is appropriately classified in terms of the manifold geometry
     ManifoldType mtype = consistencyCheck( merged_manifold );
     merged_manifold.GeometricClassifier( mtype );
     
     // 3. replacing the first manifold with the new one
     *mnf1 = merged_manifold;
     // deleting the second one
     Delete( mnf2 );

     return true;
     
  } // end MergeManifolds
*/
 



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
         csmp_error.notice( ERROR, "NodeManifoldManager::OutputNodeManifoldsToBinary:",
                           "no node manifolds found; call this method only if there are SplitBoundaries in the model" );
         return;
      }

    std::string bin_file( file_name );
    std::fstream fp( bin_file.c_str(), std::ios::out | std::ios::binary );
    if ( !fp.is_open() ) {
        csmp_error.notice( ERROR, "NodeManifoldManager::OutputNodeManifoldsToBinary:",
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
          topo_classifiers.push_back( static_cast<int8_t>(nmf.GeometricClassifier()) );
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
        vector<uint32_t>  manifold_node_list;
        manifold_node_list.reserve( n_manifold_node_entries );
        for ( auto& nmf : node_manifolds_ ) {
             const size_t entries(nmf.Branches());
             for ( auto i{0}; i<entries; ++i )
               manifold_node_list.push_back( nmf.N(i)->Idx() );
          }
        binaryFileWrite( fp, manifold_node_list );
      }
      // 4. like previous record but of INTERFACE_SIDE specifiers
      {
        vector<int8_t>  manifold_node_topo_list;
        manifold_node_topo_list.reserve( n_manifold_node_entries );
        for ( auto& nmf : node_manifolds_ ) {
             const size_t entries(nmf.Branches());
             for ( auto i{0}; i<entries; ++i )
               manifold_node_topo_list.push_back( nmf.InterFaceSide(i) );
          }
        binaryFileWrite( fp, manifold_node_topo_list );
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
        csmp_error.notice( ERROR, "NodeManifoldManager::InputNodeManifoldsFromBinary:",
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
    vector<int8_t>  manifold_topology;
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

    // 4. like previous record but of INTERFACE_SIDE specifiers
    vector<int8_t>  topo_of_nodes;
    binaryFileRead( fp, topo_of_nodes );


    // ------------------------------------
    // 3. reconstructing the node manifolds 
    // ------------------------------------
    set<pair<size_t,INTERFACE_SIDE> > manifold_nodes;
    size_t counter(0U);
    for ( auto i{0}; i<manifolds; ++i ) {
         const auto n_branches( nodes_per_manifold[i] );
         for ( auto j=0U; j<n_branches; ++j ) {
              assert( nodes_of_manifolds[counter] < mesh_nodes.size() );
              assert( topo_of_nodes[counter] <= OUTSIDE );
              manifold_nodes.insert( make_pair( nodes_of_manifolds[counter],
                                                static_cast<INTERFACE_SIDE>(topo_of_nodes[counter]) ) );
              counter++;
           }
         assert( manifold_topology[i] < static_cast<int8_t>(ManifoldType::NOT_CLASSIFIED) );
         const ManifoldType topology = static_cast<ManifoldType>(manifold_topology[i]) ;
         node_manifolds_.emplace( NodeManifold<dim>( mesh_nodes, manifold_nodes, topology ) );
      }
    
    return current_sort_variable_;
 }
    
    

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
