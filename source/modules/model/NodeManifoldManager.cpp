#include "NodeManifoldManager.h"
#include "ErrorHandler.h"
#include "Node.h"
#include "binaryReadWrite.h"


using namespace std;

namespace csmp {

// disconnecting and deleting the node manifolds
template<size_t dim>
NodeManifoldManager<dim>::~NodeManifoldManager<dim>()
 {
    // detaching the nodes from the manifolds
    for ( auto& nit : node_manifolds_ ) {
         delete nit;
         nit = nullptr;
      }
 }

template<size_t dim>
typename NodeManifoldManager<dim>::manifoldIterator NodeManifoldManager<dim>::ManifoldsBegin()
  { return node_manifolds_.begin(); }


template<size_t dim>
typename NodeManifoldManager<dim>::manifoldIterator NodeManifoldManager<dim>::ManifoldsEnd()
  { return node_manifolds_.end(); }

template<size_t dim>
typename NodeManifoldManager<dim>::manifoldConstIterator NodeManifoldManager<dim>::ManifoldsBegin() const
  { return node_manifolds_.begin(); }


template<size_t dim>
typename NodeManifoldManager<dim>::manifoldConstIterator NodeManifoldManager<dim>::ManifoldsEnd() const
  { return node_manifolds_.end(); }


template<size_t dim>
size_t  NodeManifoldManager<dim>::Manifolds() const
  { return node_manifolds_.size(); }




/**
    Creates double-node manifold with the classified meanings.
      
      @attention one can later add nodes to this manifold using  NodeManifold's methods
*/
template<size_t dim>
NodeManifold<dim>* const NodeManifoldManager<dim>::NewManifold( Node<dim>* const inside,
                                                                Node<dim>* const outside,
                                                                ManifoldType geom )
 {
    node_manifolds_.push_back( new NodeManifold<dim>( vector<Node<dim>*>({inside,outside}),
                                                      vector<INTERFACE_SIDE>({INSIDE,OUTSIDE}), geom ) );
    return node_manifolds_.back();
 }






    /// puts the nodes inside of the manifolds into the ascending order of values of the user specified  variable
template<size_t dim>
void NodeManifoldManager<dim>::SortManifoldsByVariableValue( std::string var_name, const csmp::Index& var_index )
 {
    for ( auto& nmf : node_manifolds_ ) {
         nmf->SortByVariableValue( var_index );
      }
    current_sort_variable_ = var_name;
   
 } // end SortManifoldsByVariableValue

 
 
 



template<size_t dim>
bool NodeManifoldManager<dim>::Delete( NodeManifold<dim>* md )
{
    sort( node_manifolds_.begin(), node_manifolds_.end() );
    if ( find(node_manifolds_.begin(), node_manifolds_.end(), md) != node_manifolds_.end() ) {
          remove(node_manifolds_.begin(),node_manifolds_.end(),md);
          return true;
      } else {
          ErrorHandler&  csmp_error( ErrorHandler::Instance() );
          csmp_error.notice( INFO, "NodeManifoldManager<dim>::Delete:",
                                   "NodeManifold does not exist. Nothing was done.");
      }
    return false;
}


/**
  If there are single-node manifolds, their nodes are deconnected and the manifolds are deleted.
  Later, the manager is searched for the remaining null pointers and these are removed.
*/
template<size_t dim>
size_t NodeManifoldManager<dim>::DeleteSingleNodeManifolds()
  {
     size_t n_single_node_manifolds(0U);
     for ( auto& nit : node_manifolds_ )
       if ( nit->Branches() < 2 ) {
            // disconnecting the remaining node
            nit->Remove( nit->N(0) );
            // deleting the manifold
            delete nit;
            nit = nullptr;
            n_single_node_manifolds++;
         }
     
      // compacting the manifold container
      sort( node_manifolds_.begin(), node_manifolds_.end() );
      node_manifolds_.erase( unique(node_manifolds_.begin(), node_manifolds_.end()), node_manifolds_.end() );
     
      return n_single_node_manifolds;
  }




/**
    writes manifolds to binary file; relies on unique node indexes
    
    @attention call DeleteSingleNodeManifolds() before to avoid storing any defunct manifolds
    @attention make sure that the nodes have unique numbers that match those in the other binary files
*/
template<size_t dim>
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
          topo_classifiers.push_back( static_cast<int8_t>(nmf->GeometricClassifier()) );
        binaryFileWrite( fp, topo_classifiers );
      }

      // 2. record of nodes per manifold (size, values)
      size_t n_manifold_node_entries(0U);
      {
        vector<int8_t>  nodes_per_manifold;
        nodes_per_manifold.reserve( records );
        for ( auto& nmf : node_manifolds_ ) {
             nodes_per_manifold.push_back( static_cast<int8_t>(nmf->Branches()) );
             n_manifold_node_entries += nmf->Branches();
          }
        binaryFileWrite( fp, nodes_per_manifold );
      }

      // 3. 'plist'-like record of node indices
      //   (the entries will have been sorted by variable values)
      {
        // name of variable that was used for sorting the nodes
        const string sort_variable(node_sorting_variable);
        binaryFileWrite( fp, sort_variable );
        // sorted nodes
        vector<size_t>  manifold_node_list;
        manifold_node_list.reserve( n_manifold_node_entries );
        for ( auto& nmf : node_manifolds_ ) {
             const size_t entries(nmf->Branches());
             for ( size_t i=0U; i<entries; ++i )
               manifold_node_list.push_back( nmf->N(i)->Idx() );
          }
        binaryFileWrite( fp, manifold_node_list );
      }
      // 4. like previous record but of INTERFACE_SIDE specifiers
      {
        vector<int8_t>  manifold_node_topo_list;
        manifold_node_topo_list.reserve( n_manifold_node_entries );
        for ( auto& nmf : node_manifolds_ ) {
             const size_t entries(nmf->Branches());
             for ( size_t i=0U; i<entries; ++i )
               manifold_node_topo_list.push_back( nmf->InterFaceSide(i) );
          }
        binaryFileWrite( fp, manifold_node_topo_list );
      }
  } // end writing the manifold records
  
} // end OutputNodeManifoldsToBinary
 
 
 
 
 

/** Reads manifolds from binary file, using indices to create pointer connections
 */
template<size_t dim>
string NodeManifoldManager<dim>::InputNodeManifoldsFromBinary( std::deque<Node<dim>*>& mesh_nodes, const char* file_name )
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
    vector<size_t>  nodes_of_manifolds;
    binaryFileRead( fp, nodes_of_manifolds );

    // 4. like previous record but of INTERFACE_SIDE specifiers
    vector<int8_t>  topo_of_nodes;
    binaryFileRead( fp, topo_of_nodes );


    // ------------------------------------
    // 3. reconstructing the node manifolds 
    // ------------------------------------
    vector<Node<dim>*>     nodes;
    vector<INTERFACE_SIDE> sides;
    size_t counter(0U);
    for ( size_t i=0U; i<manifolds; ++i ) {
         const size_t n_branches( nodes_per_manifold[i] );
         nodes.reserve( n_branches );
         sides.reserve( n_branches );
         for ( size_t j=0U; j<n_branches; ++j ) {
              assert( nodes_of_manifolds[counter] < mesh_nodes.size() );
              nodes.push_back( mesh_nodes[ nodes_of_manifolds[counter] ] );
              assert( topo_of_nodes[counter] <= OUTSIDE );
              sides.push_back( static_cast<INTERFACE_SIDE>(topo_of_nodes[counter]) );
              counter++;
           }
         assert( manifold_topology[i] < static_cast<int8_t>(ManifoldType::NOT_CLASSIFIED) );
         const ManifoldType topology = static_cast<ManifoldType>(manifold_topology[i]) ;
         node_manifolds_.push_back( new NodeManifold<dim>( nodes, sides, topology ) );
      }
    
    return current_sort_variable_;
 }
    



template class NodeManifoldManager<1U>;
template class NodeManifoldManager<2U>;
template class NodeManifoldManager<3U>;

} // csmp
