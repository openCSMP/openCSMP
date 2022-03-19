#include "MeshManager.h"
#include "MeshManagementUtilities.h"
#include "NodeManifoldManager.h"
#include "PropertyDatabase.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "ModelSubDomain.h"
#include "VSet.h"
#include "ErrorHandler.h"
#include "PropertyData.h"
#include "Box.h"
#include "ModelTopology.h"
// #include "MeshIterator.h"
#include "CSMP_highLevelUtilities.h"

#define MESH_MANAGER_DEBUG

using namespace std;

namespace csmp {

template<uint32_t dim>
MeshManager<dim>::MeshManager()
  : fem_manager_( dim, 1, true ), // linear interpolation functions, isoparametric elements
    fvm_manager_(fem_manager_),
    hybrid_element_mesh_( false )
{
}




/**
      input objects need to be fully constructed for this to work.
*/
template<uint32_t dim>
MeshManager<dim>::MeshManager( const PropertyDatabase<dim>& pref, const VSet<dim>& vset )
  : fem_manager_( dim, vset.OrderOfFiniteElementInterpolationFunctions(), vset.IsoparametricElementMesh() ),
    fvm_manager_(fem_manager_),
    hybrid_element_mesh_( vset.HybridElementTypeMesh() )
{
   assert( pref.VariableCount() > 0 );
   assert( vset.Vertices() > 0 );
   Initialize( pref, vset );
}




/**
     Deallocate all the dynamically allocated nodes, elements, faces, interfaces.
     
     @note Only the nodes of the volumetric elements are deleted because they are shared  with the surface mesh.
     @note Any pointer used to delete an object is set to null after the deletion.
          
 */
template<uint32_t dim>
MeshManager<dim>::~MeshManager()
  {
     delete node_manifold_manager_;
    
 } // end destructor








  /// returns number of nodes=vertices in the current mesh
template<uint32_t dim>
size_t MeshManager<dim>::Nodes() const
 { return nodes_.size(); }

  /// returns number of elements in the current mesh
template<uint32_t dim>
size_t MeshManager<dim>::Elements() const
 { return elements_.size(); }

  /// returns number of Faces=lower-dimensional elements in current mesh
template<uint32_t dim>
size_t MeshManager<dim>::Faces() const
 { return faces_.size(); }

  /// returns number of InterFaces=lower-dimensional elements in current mesh
template<uint32_t dim>
size_t MeshManager<dim>::InterFaces() const
 { return interfaces_.size(); }



template<uint32_t dim>
  typename plf::colony<Node<dim>>::iterator      MeshManager<dim>::NodesBegin()
  { return nodes_.begin(); }
  
template<uint32_t dim>
  typename plf::colony<Node<dim>>::iterator      MeshManager<dim>::NodesEnd()
  { return nodes_.end(); }

template<uint32_t dim>
  typename plf::colony<Element<dim>>::iterator   MeshManager<dim>::ElementsBegin()
  { return elements_.begin(); }
  
template<uint32_t dim>
  typename plf::colony<Element<dim>>::iterator   MeshManager<dim>::ElementsEnd()
  { return elements_.end(); }

template<uint32_t dim>
  typename plf::colony<Face<dim>>::iterator      MeshManager<dim>::FacesBegin()
  { return faces_.begin(); }

template<uint32_t dim>
  typename plf::colony<Face<dim>>::iterator      MeshManager<dim>::FacesEnd()
  { return faces_.end(); }

template<uint32_t dim>
  typename plf::colony<InterFace<dim>>::iterator MeshManager<dim>::InterFacesBegin()
  { return interfaces_.begin(); }
  
template<uint32_t dim>
  typename plf::colony<InterFace<dim>>::iterator MeshManager<dim>::InterFacesEnd()
  { return interfaces_.end(); }

  // const versions
template<uint32_t dim>
  typename plf::colony<Node<dim>>::const_iterator      MeshManager<dim>::NodesBegin() const
  { return nodes_.begin(); }

template<uint32_t dim>
  typename plf::colony<Node<dim>>::const_iterator      MeshManager<dim>::NodesEnd() const
  { return nodes_.end(); }

template<uint32_t dim>
  typename plf::colony<Element<dim>>::const_iterator   MeshManager<dim>::ElementsBegin() const
  { return elements_.begin(); }

template<uint32_t dim>
  typename plf::colony<Element<dim>>::const_iterator   MeshManager<dim>::ElementsEnd() const
  { return elements_.end(); }

template<uint32_t dim>
  typename plf::colony<Face<dim>>::const_iterator      MeshManager<dim>::FacesBegin() const
  { return faces_.begin(); }

template<uint32_t dim>
  typename plf::colony<Face<dim>>::const_iterator      MeshManager<dim>::FacesEnd() const
  { return faces_.end(); }

template<uint32_t dim>
  typename plf::colony<InterFace<dim>>::const_iterator MeshManager<dim>::InterFacesBegin() const
  { return interfaces_.begin(); }

template<uint32_t dim>
  typename plf::colony<InterFace<dim>>::const_iterator MeshManager<dim>::InterFacesEnd() const
  { return interfaces_.end(); }



// NODE MANIFOLD ITERATORS
  
template<uint32_t dim>
typename plf::colony<NodeManifold<dim>>::iterator MeshManager<dim>::NodeManifoldsBegin() {
     if ( node_manifold_manager_ == nullptr )
       throw csmp::Exception( ERROR, "MeshManager<dim>::NodeManifoldsBegin", "current model has no node manifolds");
     return node_manifold_manager_->ManifoldsBegin();
  }

template<uint32_t dim>
typename plf::colony<NodeManifold<dim>>::iterator MeshManager<dim>::NodeManifoldsEnd() {
     if ( node_manifold_manager_ == nullptr )
       throw csmp::Exception( ERROR, "MeshManager<dim>::NodeManifoldsEnd", "current model has no node manifolds");
     return node_manifold_manager_->ManifoldsEnd();
  }
  
template<uint32_t dim>
typename plf::colony<NodeManifold<dim>>::const_iterator MeshManager<dim>::NodeManifoldsBegin() const {
     if ( node_manifold_manager_ == nullptr )
       throw csmp::Exception( ERROR, "MeshManager<dim>::NodeManifoldsBegin", "current model has no node manifolds");
     return node_manifold_manager_->ManifoldsBegin();
  }

template<uint32_t dim>
typename plf::colony<NodeManifold<dim>>::const_iterator MeshManager<dim>::NodeManifoldsEnd() const {
     if ( node_manifold_manager_ == nullptr )
       throw csmp::Exception( ERROR, "MeshManager<dim>::NodeManifoldsEnd", "current model has no node manifolds");
     return node_manifold_manager_->ManifoldsEnd();
  }





/**
       Range checked access of entities by their place in the storage.
       throws out_of_range if abused.
*/
/*
template<uint32_t dim>
Node<dim>* const MeshManager<dim>::N( size_t idx ) const {
    return nodes_.at( idx );
 }
  
  
template<uint32_t dim>
Element<dim>* const MeshManager<dim>::E( size_t idx ) const {
    return elements_.at( idx );
 }
 
 
template<uint32_t dim>
Face<dim>* const MeshManager<dim>::F( size_t idx ) const {
    return faces_.at( idx );
 }
 
 
template<uint32_t dim>
InterFace<dim>* const MeshManager<dim>::I( size_t idx ) const {
    return interfaces_.at( idx );
}
*/




/**
Configures the distributed variable storage, inialises the finite element container
with the element types that are contained in the current mesh and builds the
mesh including the connectivity among its elements.

While the method detects insconsistencies / omissions in the element/face/interface neighbor arrays it will continue execution assigning nullptr where incorrect neighbors are encountered.
However it will report when and where this happened.

Reads the VSet in the format that is espoused by the method OutputMeshTo( VSet )

@note the neighbors of each element include only the elements of the same type, i.e.
a line element only has line neighbors, a surface element surface element neighbors
and so forth.

@attention lower dimensional elements may have multiple neighbors for each of their
faces. Yet only one of them will be assigned.

@todo SKM: create manifolds to deal with lower-dimensional elements that have
multiple neighbors per face.

@attention IMPORTANT: the order of entries in  plf::colony is not guaranteed once an erasure has occurred. This means
that all operations that build mesh with reference to the input VSet indexing must be complete before erasures occur.

*/
template<uint32_t dim>
bool  MeshManager<dim>::Initialize( const PropertyDatabase<dim>& phys_vars, const VSet<dim>& vset )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  hybrid_element_mesh_ = vset.HybridElementTypeMesh();
  cout <<"\nMeshManager<"<< dim <<">::Initialize: building mesh with "<< vset.Elements() <<" elements, "<< vset.Vertices() <<" nodes, ";
  cout << vset.Faces() <<" faces, and "<< vset.InterFaces() <<" interfaces.\n";
  if ( vset.HybridElementTypeMesh() ) cout <<"mesh consists of multiple element types.\n";
  if ( vset.Faces() > 0 ) cout <<"mesh contains 'Boundary' objects.\n";
  if ( vset.InterFaces() > 0 ) cout <<"mesh contains 'SplitBoundary' objects.\n";
  cout << endl;

  // ------------------------------------------------------------------------------------
  // 1. check availability of necessary finite element types, and valid model topology
  // ------------------------------------------------------------------------------------
  // initializing the finite-element manager true=isoparametric
  fem_manager_.InitializeElements( dim, vset.OrderOfFiniteElementInterpolationFunctions(), vset.IsoparametricElementMesh() );

  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: checking the availability of the necessary finite element types..." << endl;

  set<CSMP_FEM_TYPE> input_etypes;
  input_etypes.insert( parseFiniteElementTypeEnum( vset.ElementType( 0U ) ) );

  if ( vset.HybridElementTypeMesh() )
    for ( auto i = 0U; i<vset.TotalNumberOfCells(); ++i )
      input_etypes.insert( parseFiniteElementTypeEnum( vset.ElementType( i ) ) );

  cout << "\nMeshManager<" << dim << ">::Initialize: ";
  cout << "input VSet contains the following finite element types:\n\t";
  for ( typename set<CSMP_FEM_TYPE>::const_iterator
        iit = input_etypes.begin(); iit != input_etypes.end(); iit++ ) {
      cerr << parseFiniteElementType( (*iit) ) << "  ";
      if ( !fem_manager_.ContainsElementType( *iit ) ) {
        cerr << "\n\n\tFinite element type not available: " << parseFiniteElementType( *iit ) << endl;
        fem_manager_.Out();
        throw Exception( FATAL_ERROR,
                         "MeshManager<dim>::Initialize(VSet):",
                         "'FiniteElementManager' lacks finite-element type required by VSet." );
      }
    }
  cout << endl;

  // --------------------------------------------------------------------------
  // 2. construct nodes and elements using the VSet element type information
  // --------------------------------------------------------------------------
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: building storage and assigning nodes to elements..." << endl;

  // Storage for nodes
   {
      vector<double> coord( dim );
      const LocalVariables nvars( phys_vars.LocalVariablesAt( NODE ) );
      for ( size_t idx = 0U; idx < vset.Vertices(); ++idx ) {
          for ( auto j = 0U; j<dim; ++j ) coord[j] = vset.P( j, idx );
          nodes_.emplace( Node<dim>( idx, Point<dim>( coord ), nvars, static_cast<BOX_BOUNDARY>(vset.BFlag(idx)) ) );
        }
    }

  // storage for elements
   {
      const LocalVariables evars( phys_vars.LocalVariablesAt( ELEMENT ) );
      const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( ELEMENT ) );
      typename deque<vector<int64_t>>::const_iterator first( vset.PlistElmtsBegin() ), last( vset.PlistElmtsEnd() );

      size_t elmt_idx{0};
      // 2.1 If the MeshManager contains only one element type
      if ( !vset.HybridElementTypeMesh() ) {
          const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( 0U ));
          while ( first != last )
            {
              // create the element
              typename plf::colony<Element<dim>>::iterator
                eit = elements_.emplace( Element<dim>( elmt_idx, fem_manager_.E( csmpElementType ), fvm_manager_.Stencil( csmpElementType ),
                                                                                 evars, cvars, vset.Pmtrl(elmt_idx) ) );
              // assign the nodes
              const auto nodes( fem_manager_.E( csmpElementType )->Nodes() );
              for ( auto j = 0U; j < nodes; ++j )
                (*eit).Assign( j, &(*next(nodes_.begin(),vset.Plist( elmt_idx, j ))) );
                
              // assign the material
              (*eit).Material_ID( vset.Pmtrl( elmt_idx ) );
                
              elmt_idx++;
              first++;
            }
        }
      // 2.2 If there are multiple element types
      else {
          while ( first != last )
            {
              const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( elmt_idx ));
              typename plf::colony<Element<dim>>::iterator
                eit = elements_.emplace( Element<dim>( elmt_idx, fem_manager_.E( csmpElementType ), fvm_manager_.Stencil( csmpElementType ),
                                                                                 evars, cvars, vset.Pmtrl(elmt_idx) ) );
              const auto nodes( fem_manager_.E( csmpElementType )->Nodes() );
              for ( auto j = 0U; j < nodes; j++ ) (*eit).Assign( j, &(*next(nodes_.begin(),vset.Plist( elmt_idx, j ))) );
              elmt_idx++;
              first++;
            }
        }
      assert( elements_.size() == vset.Elements() );
   }
 
  // 2.3 Assign neighbor elements to elements
  const int64_t  n_elmts(elements_.size());
  if ( vset.WithNeighbourConnectivity() ) {
       if ( csmp_error.Verbose() )
          cout << "\nMeshManager<" << dim << ">::Initialize: assigning neighbors to elements..." << endl;

       for ( auto& e : elements_ ) {
            const int8_t csmpElementType = (!hybrid_element_mesh_) ? vset.ElementType( 0U ) : vset.ElementType( e.Idx() );
            const size_t n_neighbors( fem_manager_.E( csmpElementType )->Neighbors() );
//assert( n_neighbors == distance(vset.PfvertsBegin(e.Idx()),vset.PfvertsEnd(e.Idx())) );

            for ( auto j = 0U; j < n_neighbors; ++j ) {
                  const int64_t  index{ vset.Pfvert( e.Idx(), j ) };
                  if ( index >= n_elmts ) {
                       cerr <<"\n\t"<< index <<" vs. number of elements = "<< n_elmts << endl;
                       csmp_error.notice( ERROR, "MeshManager::Initialise: ", "element ID in 'pfverts' out of range.");
                    }
                  else if ( index >= 0 )
                    e.Assign( j, &(*next(elements_.begin(),index)) );
                  // negative numbers indicate boundaries where nothing needs to be done since neighbors are initialised to null pointers anyway
                  //  e.Assign( nidx++, static_cast<Element<dim>*>(nullptr) );
              }
         }
    }
  else
  csmp_error.notice( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity; nothing was done." );

  // ----------------------------------------------------------
  // 3. constructing the Faces using the VSet node information
  // ----------------------------------------------------------
  // (continuous running index 'idx' will be used so that face-IDs start with n-elements)
  // different face types are intrinsic to the VSet if so initialized
  const size_t n_nodes(nodes_.size());
  if ( vset.Faces() > 0 )
    {
       assert( vset.HybridElementTypeMesh() );
       if ( csmp_error.Verbose() )
         cout << "\nMeshManager<" << dim << ">::Initialize: assigning nodes to faces..." << endl;
      
       // since faces are lower-dimensional, the mesh must contain different element types
       const LocalVariables evars( phys_vars.LocalVariablesAt( FACE ) );
       const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( FACE ) );

       // the faces are numbered  elements to (elements + faces - 1), but they are stored in connector at Face 0..n-1
       int64_t   face_idx(vset.Elements());
       typename deque<vector<int64_t> >::const_iterator  first( vset.PlistFacesBegin() ), last( vset.PlistFacesEnd() );
       while ( first != last ) {
            const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( face_idx ));
            if ( csmpElementType == UNKNOWN ) {
                 cerr <<"\n\t"<< parseFiniteElementType(csmpElementType) <<" encountered for Face "<< face_idx <<"\n";
                 csmp_error.notice( FATAL_ERROR, "MeshManager::Initialise:", "encountered UNKNOWN Face element type." );
              }
            typename plf::colony<Face<dim>>::iterator
               fit = faces_.emplace( Face<dim>( face_idx, fem_manager_.E( csmpElementType ),
                                                          fvm_manager_.Stencil( csmpElementType ), evars, cvars ) );
            // assigning nodes to faces
            const auto nodes( (*fit).Nodes() );
            for ( auto j = 0U; j<nodes; ++j ) {
                const size_t node = vset.Plist( face_idx, j );
                assert( node < n_nodes );
                (*fit).Assign( j, &(*next(nodes_.begin(),node)) );
              }
            ++face_idx;
            ++first;
          }
       assert( faces_.size() == vset.Faces() );

       // connecting the faces to their equi- and higher-dimensional neighbors
       // --------------------------------------------------------------------
       // necessary info is stored in 'pfverts' record for Face
       if ( !vset.WithNeighbourConnectivity() )
         csmp_error.notice( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity for Face objects; nothing was done." );
       else
         {
            if ( csmp_error.Verbose() )
              cout << "\nMeshManager<" << dim << ">::Initialize: connecting faces to their equidimensional and higher-dimensional neighbors..." << endl;
            const int64_t  n_faces(faces_.size());
            for ( auto& e : faces_ )
              {
                 // Equidimensional Face neighbors first
                 // ------------------------------------
                 const auto neighbors( e.Neighbors() );
                 for ( auto j = 0U; j<neighbors; ++j )
                   {
                      // if there is a neighbor (as is the case if the stored index is greater than zero)
                      // (e->Idx() starts with elements=first face)
                      const int64_t  index( vset.Pfvert( e.Idx(), j ) );
                      if ( index >= n_elmts+n_faces ) {
                           cerr <<"\n\t"<< index <<" vs. number of elements+faces = "<< n_elmts + n_faces << endl;
                           csmp_error.notice( ERROR, "MeshManager::Initialise: ", "face ID in 'pfverts' out of range.");
                        }
                      // if the Face neighbor has an index smaller than n_elmts it must be a boundary indicator
                      if ( index >= n_elmts )
                        e.Assign( j, &(*next(faces_.begin(),index - n_elmts)) );
                      else {
                           assert( index < 0 );
                           e.Assign( j, static_cast<Face<dim>*>(nullptr) );
                        }
                   }
              
                 // Higher-dimensional Element neighbors (2) of Face
                 // ------------------------------------------------
                 // (are stored in VSet 'pfverts' record after the equidimensional neighbors)
                 // index of inner neighbor element i which is always there
                 // Inside neighbor 1
                 const int64_t  index1( vset.Pfvert( e.Idx(), neighbors ) );
                 if ( index1 >= n_elmts ) {
                      cerr <<"\n\tFace "<< e.Idx() <<": "<< index1 <<" vs. "<< n_elmts <<" elements.\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "Index of first higher-dimensional element of Face out of range.");
                   }
                 else if ( index1 < 0 ) {
                      cerr <<"\n\tFace "<< e.Idx() <<": "<< index1 <<"\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "First higher-dimensional element out of range.");
                   }

                 // Outside neighbor 2: outer neighbor element will only be there if Face on INTERNAL model boundary
                 const int64_t  index2( vset.Pfvert( e.Idx(), neighbors + 1U ) );
                 if ( index2 >= n_elmts ) {
                      cerr <<"\n\tFace "<< e.Idx() <<": "<< index2 <<" vs. "<< n_elmts <<" elements.\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "Index of second higher-dimensional element of Face out of range.");
                   }
                 // assignment of higher-dimensional neighbors
                 Element<dim>* const innerElement = (index1 < 0) ? nullptr : &(*next(elements_.begin(),index1));
                 Element<dim>* const outerElement = (index2 < 0) ? nullptr : &(*next(elements_.begin(),index2));
                 // assigning inner and outer higher-dimensional neighbors
                 e.Assign( innerElement, outerElement );
                 // and the corresponding face numbers of these higher dimensional elements
                 e.ParentFaceID( INSIDE, static_cast<uint32_t>(vset.Pfvert( e.Idx(), neighbors + 2U )) );
                 if ( outerElement )
                   e.ParentFaceID( OUTSIDE, static_cast<uint32_t>(vset.Pfvert( e.Idx(), neighbors + 3U )) );
                 
               } // end face loop
              
             } // end else
        
    } // end construction and initialisation of Face objects


  // ------------------------------------------------------------------
  // 4. constructing Interfaces objects using the VSet node information
  // ------------------------------------------------------------------
  // (continuous running index 'idx' will also be used for interfaces)
  if ( vset.InterFaces() > 0 )
    {
       assert( vset.HybridElementTypeMesh() );
       if ( csmp_error.Verbose() )
         cout << "\nMeshManager<" << dim << ">::Initialize: assigning nodes to interfaces..." << endl;
       
       const LocalVariables evars( phys_vars.LocalVariablesAt( INTER_FACE ) );
       const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( INTER_FACE ) );

       typename deque<vector<int64_t> >::const_iterator  first( vset.PlistInterFacesBegin() ),
                                                        last( vset.PlistInterFacesEnd() );

       size_t interface_idx(vset.Elements() + vset.Faces());
       while ( first != last )
         {
            const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( interface_idx ));
            typename plf::colony<InterFace<dim>>::iterator
              ifit = interfaces_.emplace( InterFace<dim>( interface_idx,
                                                          fem_manager_.E( csmpElementType ),
                                                          fvm_manager_.Stencil( csmpElementType ),
                                                          evars, cvars ) );
                                                       
            // number of nodes of the finite-element corresponding to the interface
            const auto nodes( (*ifit).FE()->Nodes() );
            // assigning nodes
            // inside
            for ( auto j = 0U; j<nodes; ++j ) {
                 const size_t node(vset.Plist( interface_idx, j ));
                 if ( node >= n_nodes ) {
                      cerr <<"\n\tInterFace "<< interface_idx <<": INSIDE node j "<< node <<" vs. "<< n_nodes <<" nodes.\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "Index of InterFace node out of range.");
                   }
                 (*ifit).Assign( j, &(*next(nodes_.begin(),node)), INSIDE );
              }
            // outside
            for ( auto j = 0U; j<nodes; ++j ) {
                 const size_t node(vset.Plist( interface_idx, j+nodes ));
                 if ( node >= n_nodes ) {
                      cerr <<"\n\tInterFace "<< interface_idx <<": OUTSIDE node j "<< node <<" vs. "<< n_nodes <<" nodes.\n";
                      csmp_error.notice( ERROR, "MeshManager::Initialise", "Index of InterFace node out of range.");
                   }
                 (*ifit).Assign( j, &(*next(nodes_.begin(),node)), OUTSIDE );
              }
            ++interface_idx;
            ++first;
          }
       assert( interfaces_.size() == vset.InterFaces() );

       // connecting the interfaces to their equi- and higher-dimensional neighbors
       // -------------------------------------------------------------------------
       // necessary info is stored in 'pfverts' record for InterFace: inner nbors first, then outer, then higher-dimensional ones
       if ( !vset.WithNeighbourConnectivity() )
         csmp_error.notice( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity for InterFace objects; nothing was done." );
       else {
         if ( csmp_error.Verbose() )
           cout << "\nMeshManager<" << dim << ">::Initialize: connecting interfaces to their higher-dimensional neighbors..." << endl;
         const int64_t  n_faces(faces_.size()), n_interfaces(interfaces_.size());
         const int64_t  cells(n_elmts+n_faces+n_interfaces);
         // connecting interfaces to their higher-dimensional neighbors
         for ( auto& e : interfaces_ )
           {
              // assigning equidimensional InterFace-type neighbors first
              // --------------------------------------------------------
              const auto neighbors( e.Neighbors() );
              for ( auto j = 0U; j<neighbors; ++j )
                {
                   // if there is a neighbor (as is the case if the stored index is greater than zero)
                   const int64_t  index( vset.Pfvert( e.Idx(), j ) );
                   if ( index >= cells ) {
                        cerr <<"\n\t"<< index <<" vs. number of elements+faces+interfaces = "<< cells << endl;
                        csmp_error.notice( ERROR, "MeshManager::Initialise: ", "interface ID in 'pfverts' out of range.");
                     }
                   assert( index >= vset.Elements() + vset.Faces() );
                   assert( index < vset.Elements() + vset.Faces() + vset.InterFaces() ); // (-) elements because face container is numbered from 0..n-1
                  
                   // the number of the interface in the container is the number from the VSet - elements and faces
                   // because the interface container is counts from 0..n-1
                   e.Assign( j, &(*next(interfaces_.begin(),index - n_elmts - n_faces)) );
                }

             // higher-dimensional Element-type neighbors
             // -----------------------------------------
             // (higher-dimensional neighbors are always present on both sides of the InterFace because interfaces exist only inside of a model)
             if ( vset.Pfvert( e.Idx(), neighbors ) < 0 || vset.Pfvert( e.Idx(), neighbors + 1U) < 0 )
               {
                  cerr <<"\n\tInterFace "<< e.Idx() <<": inner neighbor "<< vset.Pfvert( e.Idx(), neighbors );
                  cerr <<" and outer "<< vset.Pfvert( e.Idx(), neighbors+1U ) <<"\n";
                  csmp_error.notice( ERROR, "MeshManager::Initialise: ", "Higher dimensional neighbor of InterFace not defined in 'pfverts'.");
               }
             if ( vset.Pfvert( e.Idx(), neighbors ) >= cells || vset.Pfvert( e.Idx(), neighbors + 1U ) >= cells )
               {
                  cerr <<"\n\tInterFace "<< e.Idx() <<": inner neighbor "<< vset.Pfvert( e.Idx(), neighbors );
                  cerr <<" and outer "<< vset.Pfvert( e.Idx(), neighbors+1U ) <<"\n";
                  csmp_error.notice( ERROR, "MeshManager::Initialise: ", "Higher dimensional neighbor of InterFace out of range.");
               }
             // assignment: inner and outer Element objects
             const int64_t  index1 = vset.Pfvert( e.Idx(), neighbors );
             const int64_t  index2 = vset.Pfvert( e.Idx(), neighbors+1U );
             assert( index1 < cells );
             assert( index2 < cells );
             assert( index1 > MULTIPLE );
             assert( index2 > MULTIPLE );
             Element<dim>* const innerElement = (index1 < 0) ? nullptr : &(*next(elements_.begin(),index1));
             Element<dim>* const outerElement = (index2 < 0) ? nullptr : &(*next(elements_.begin(),index2));
             e.Assign( innerElement, outerElement );
             // assignment: Element face numbers adjacent to InterFace
             const auto inner_face_id = static_cast<uint32_t>(vset.Pfvert( e.Idx(), neighbors+2U ));
             const auto outer_face_id = static_cast<uint32_t>(vset.Pfvert( e.Idx(), neighbors+3U ));
             assert( inner_face_id < e.Faces() );
             assert( outer_face_id < e.Faces() );
             e.ParentFaceID( INSIDE,  inner_face_id );
             e.ParentFaceID( OUTSIDE, outer_face_id );
             // assignment: intervening Element or neighbor boundary flag
             const int64_t  index3 = vset.Pfvert( e.Idx(), neighbors+4U );
             assert( index3 < elements_.size() );
             assert( index3 > MULTIPLE );
             Element<dim>* const middleElement = (index3 < 0) ? nullptr : &(*next(elements_.begin(),index3));
             e.Assign( middleElement );
          }
      }
      
    } // end construction of Interface objects


  // ---------------------------------------------------------------------
  // 5. Flagging nodes at model boundary with BOX_BOUNDARY flags
  // ---------------------------------------------------------------------
  // NB: the node flags were already assigned further above where the nodes were created!
    
    
  // ------------------------------------------------------------------------------
  // 6. Assigning parent elements (these are the elements that share the node) and
  //    their respective internal node-id numbers to the nodes
  // ------------------------------------------------------------------------------
  // tested: OK
  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: assigning parent element information to nodes..." << endl;
  vector<uint32_t>  parent_elmts_per_node( vset.Vertices(), 0U );

  // counting how many parent elements each node has
  for ( const auto& e : elements_ ) {
      for ( auto nit = e.NodesBegin(); nit != e.NodesEnd(); ++nit ) {
           assert( (*nit)->Idx() < n_nodes );
           parent_elmts_per_node[ (*nit)->Idx() ]++;
        }
    }

  // reserving the memory for the parent storage and zeroing parent vector for next step
  auto i{0};
  for ( auto& n : nodes_ )
    n.ResizeParentStorage( parent_elmts_per_node[i++] );

   // assigning the parent element information to the nodes
   for ( auto& e : elements_ )
       for ( auto j = 0U; j<e.Nodes(); ++j )
         e.N( j )->Assign( j, &e );
       
   // sorting the parent element pointers stored by the nodes for searching
   for ( auto& n : nodes_ )
    n.SortParents();

   if ( csmp_error.Verbose() )
     cout << "\nMeshManager<" << dim << ">::Initialize: forming regions for contiguous subdomains..." << endl;
 

 
  // ------------------------------------------------------------------------------
  // 7. Connecting the nodes with one another through the neighbor pointers
  // ------------------------------------------------------------------------------
   // OK - consecutively numbered:  for ( auto nit : nodes_ ) cout <<" "<< nit.Idx();
   vector<set<size_t> > pnode;
   vset.EstablishNodeNeighborConnectivity( pnode );
   
   size_t node{0};
   for ( auto& n : nodes_ ) {
        // creating sorted vector of node points
        vector<Node<dim>*>  nptrs;
        nptrs.reserve( pnode[node].size() );
        for ( auto s : pnode[node] )
          nptrs.push_back( &(*next(nodes_.begin(),s)) );
        // assigning the neighbor pointers to the node
        n.Assign( nptrs, true );
        node++;
     }

     
  // ------------------------------------------------------------------------------
  // 8. reconstructing NodeManifolds if any
  // -------------------------------------------------------------------------------
   if ( !interfaces_.empty() ) {
       VData::vertexManifoldIndices  indexes;
       vset.ExtractNodeManifolds( indexes );
       node_manifold_manager_ = new NodeManifoldManager<dim>( indexes, nodes_ );
     }

   return true;
  
} // end Initialise

/*
#ifdef MESH_MANAGER_DEBUG
integrityCheck<dim,Element>( ElementsBegin(), ElementsEnd() );
if ( Faces() > 0 )
  integrityCheck<dim,Face>( FacesBegin(), FacesEnd() );
if ( InterFaces() > 0 ) {
      integrityCheck<dim,InterFace>( InterFacesBegin(), InterFacesEnd() );
     // add test for node manifolds
  }
#endif
*/



  // ==============================================================
  //
  // MESH MODIFICATION
  //
  // ==============================================================
  /*
      Think about use cases:
      - what are the individual steps
      - who in the hierarchy of nodes, elements is reponsible for what (Element is master?)
      - Do we need reference counting for nodes?


    // To erase elements the following steps are necessary
    // ---------------------------------------------------
    // 1. set the neighbor pointers to the element to zero
    // 2. remove the pointers from the connected nodes to this parent-element
    // 3. disconnect the nodes, deleting them when they have no parent elements anymore
    // 4. disconnect the neighbor elements
    // 5. delete element
 
  */






/**
   Recreates the node to parent connectivity in the entire mesh.

   @author SKM
   @date 13/10/21
*/
template<uint32_t dim>
void MeshManager<dim>::RebuildNodeParentElementRelationships()
{
  // renumbering the nodes
  size_t node_number{0};
  for ( auto& n : nodes_ )
    n.Idx( node_number++ );

  // counting how many parent elements each node has
  vector<uint32_t>  parent_elmts_per_node( Nodes(), 0U );
  for ( const auto& e : elements_ ) {
      const auto nodes_end{ e.NodesEnd() };
      for ( auto nit = e.NodesBegin(); nit != nodes_end; ++nit )
        parent_elmts_per_node[ (*nit)->Idx() ]++;
    }

  // reserving the memory for the parent storage and zeroing parent vector for next step
  node_number = 0;
  for ( auto& n : nodes_ ) {
       n.EraseParents();
       n.ResizeParentStorage( parent_elmts_per_node[node_number++] );
    }

   // assigning the parent element information to the nodes
   for ( auto& e : elements_ ) {
        const size_t n_nodes{ e.Nodes() };
        for ( auto j = 0U; j<n_nodes; ++j )
          e.N( j )->Assign( j, &e );
     }

} // end UpdateNodeParentElementRelationships











/**
    Inserts a new Node at the desired location.  No connection are made.
*/
template<uint32_t dim>
Node<dim>* const MeshManager<dim>::AddNodeAt( const Point<dim>& location,
                                              const LocalVariables& lvars,
                                              BOX_BOUNDARY bdry )
{
   typename plf::colony<Node <dim>>::iterator
     nit = nodes_.emplace( Node<dim>( nodes_.size(), location, lvars, bdry ) );
   return &(*nit);
}


/**
    Inserts  a new Node at the desired point, but only if there is not already a node there.
    
    @return if there is already a node at the point location, a pointer to that node is returned
    
        Search algorithm for the collocated node uses  "nearby" node as a starting point.
        
        Idea: start from nearby Node
        - loop over the neighbor nodes of the node ranking them in terms of their proximity from the target point
        - move to closest node and then repeat (remembering the shortest distance)
        - repeat until node is found while the distance decreases
        - if distance increases, the node does not exist and will be created
        - allow  to move across manifold member nodes in order to cross split boundaries
        
            /// the number of nodes that this Node is connected with
    size_t           Neighbors() const;
    /// access to any of these nodes
    Node<dim>*       Neighbor( size_t ) const;

*/
template<uint32_t dim>
Node<dim>* const MeshManager<dim>::AddNodeAtUniqueLocation( const Point<dim>& pt,
                                                            size_t nearby_node,
                                                            const LocalVariables& nvars,
                                                            BOX_BOUNDARY bflag )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
   // if the node location needs to be compared with existing ndes
   if ( nearby_node >=nodes_.size() ) {
        csmp_error.notice( WARNING, "MeshManager<dim>::AddNodeAt",
                          "nearby Node not contained in Mesh:", to_string(nearby_node) );
        // using the last node
        nearby_node = nodes_.size() - 1U;
     }

   // searching the mesh tree for a node with the same location (using the provided point as a start location)
   Node<dim>*          nptr( &(*next(nodes_.begin(),nearby_node)) );
   double              new_distance(pt.DistanceTo(nptr->Coordinate())), old_distance(1e30);
   map<double,size_t>  distances;
   // estimating a tolerance on the basis of the distance of the point to the node and the first node
   const double tolerance = 1.0e-7 * (new_distance + pt.DistanceTo((*nodes_.begin()).Coordinate())) / 2.;
   while ( old_distance > new_distance )
     {
        // tree travel: looping the neighbor nodes of the current node, finding the one that is the closest to the point
        const size_t n_nbors( nptr->Neighbors() );
        for ( auto i{0}; i<n_nbors; ++i )
          distances.insert( make_pair( pt.DistanceTo( nptr->Neighbor(i)->Coordinate() ), i ) );
        // since map defaults to less, its first entry is the node we want
        nptr = nptr->Neighbor( static_cast<uint32_t>((*distances.begin()).second) );
        old_distance = new_distance;
        new_distance = (*distances.begin()).first;
        assert( nptr != nullptr );
     }
   // TODO: deal with NodeManifolds - if IsManifold()...
   // if a node matching the point location was found, a pointer to it is returned
   if ( fabs(new_distance) < tolerance ) return nptr;
  
   // else a new node is created
   typename plf::colony<Node <dim>>::iterator
     nit = nodes_.emplace( Node<dim>( nodes_.size(), pt, nvars, bflag ) );
     
   return &(*nit);
}


/**
      Constructs element and connects it up with the supplied nodes and neighbor elements if any.
      
   if neighbors are not supplied, method tries to find neighbors through the parent connectivity of the nodes.
   Warning messages are issued if there are issues with the input data.
   
   @author SKM
   @date 17/9/21
*/
template<uint32_t dim>
Element<dim>*	const MeshManager<dim>::AddElement( CSMP_FEM_TYPE etype,
                                                  const LocalVariables& lvars,
                                                  const IntegrationPointVariables& ivars,
                                                  const std::vector<Node<dim>*>& nodes,
                                                  int32_t material_id )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // node vector
   if ( nodes.empty() )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddElement", "node vector is empty");

   // 1. constructing new element
   typename plf::colony<Element<dim>>::iterator
     eit = elements_.emplace( Element<dim>( elements_.size(),
                              fem_manager_.E(etype), fvm_manager_.Stencil(etype), lvars, ivars, material_id ) );

   // 2. assigning nodes
   const size_t n_nodes(nodes.size());
   for ( auto i{0}; i<n_nodes; ++i )
     (*eit).Assign( i, nodes[i] );


  // 3. trying to establish neighbor information from the nodes assuming that they have parent connectivity
  // ------------------------------------------------------------------------------------------------------ 
  //    checking whether the nodes have the necessary parent element information
  bool valid_parent_info(true);
  for ( auto i{0}; i<n_nodes; ++i )
    if ( (*eit).N(i)->Parents() == 0U ) {
         cerr <<"\n\tnode "<< i;
         valid_parent_info = false;
         csmp_error.notice( ERROR, "MeshManager<dim>::AddElement",
                          "neighbor information could not be created because node has no parent element info");
         return &(*eit);
      }
      
   // 4. if the nodes have parents, this method tries to find and connect the neighbors
   // ---------------------------------------------------------------------------------
   if ( connectNeighborsUsingNodeParents( &(*eit) ) < (*eit).FE()->Faces()-1 )
     csmp_error.notice( WARNING, "MeshManager<dim>::AddElement", "neighbor vector could not be used; found less neighbors than expected");
   
   return &(*eit);
  
} // AddElement







/**
     puts a lower-dimensional element inside of an InterFace, connecting it to its base pointer
     
     @attention the nodes need to be provided because they are shared among the lower-dimensional elements and collocated
          so that they cannot be told apart.
          
          @attention no neighbor connectivity is provided here because it is not known yet
*/
template<uint32_t dim>
Element<dim>*	const MeshManager<dim>::AddInterveningElement( csmp::InterFace<dim>* const ifptr,
                                                             const LocalVariables& lvars,
                                                             const IntegrationPointVariables& ivars,
                                                             const vector<Node<dim>*>& nodes,
                                                             int32_t material_id )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the InterFace
   if ( ifptr == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddInterveningElement", "finite element pointer not initialised");
   if ( ifptr->HasInterveningElement() ) {
        ifptr->InterveningElement()->Out();
        csmp_error.notice( ERROR, "MeshManager<dim>::AddInterveningElement", "InterFace already has intervening element");
     }
     
   // 1. checking the node vector
   if ( nodes.empty() )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddInterveningElement", "node vector is empty");
   if ( nodes.size() != ifptr->Nodes() )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddInterveningElement", "node vector has the wrong size");
     
   // 2. checking the validity of the node vector in debug mode
#ifdef DEBUG
   // node vector
   for ( auto i{0}; i<ifptr->FE()->Nodes(); ++i ) {
         if ( nodes[i] == nullptr ) {
              cerr <<"\n\tnode "<< i;
              csmp_error.notice( ERROR, "MeshManager<dim>::AddInterveningElement", "node vector contains a nullptr");
   break;
           }
         else if ( nodes[i]->Coordinate() != ifptr->N(i)->Coordinate() ) {
              cerr <<"\n\tnode "<< i;
              csmp_error.notice( ERROR, "MeshManager<dim>::AddInterveningElement", "node locations do not match");
   break;
           }
       }
#endif

   // 3. creating the new Element
   typename plf::colony<Element<dim>>::iterator
     eit = elements_.emplace( Element<dim>( elements_.size(), ifptr->FE(), ifptr->FV(), lvars, ivars, material_id ) );
   
   // 4. connecting the nodes to the element
   const size_t n_nodes( ifptr->FE()->Nodes() );
   for ( auto i{0}; i<n_nodes; ++i )
     (*eit).Assign( i, nodes[i] );
     
   // 5. Connecting the intervening element to interface
   ifptr->Assign( &(*eit) );
   
   return &(*eit);

} // end AddInterveningElement





/**
    Creates Face from lower-dimensional Element and assigns the higher dimensional neighbors on inside and outside.
    
    @attention the neighbor information is taken from supplied vector; nullptr entries are accepted so that the method may be used in an advancing front algorithm.
    
    @attention the original element is deleted and set to null. Later on, the MeshManager needs to be updated.
*/
template<uint32_t dim>
Face<dim>* const MeshManager<dim>::ReplaceElementByFace( csmp::Element<dim>* eptr,
                                                         csmp::Element<dim>* inner_eptr,
                                                         csmp::Element<dim>* outer_eptr,
                                                         uint32_t adjacent_face_of_inner_element,
                                                         uint32_t adjacent_face_of_outer_element,
                                                         const LocalVariables& lvars,
                                                         const IntegrationPointVariables& ivars )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( eptr == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::ReplaceElementByFace", "element pointer not initialised");

   // is the element indeed lower dimensional?
   if constexpr ( dim == 3 )
     if ( !eptr->IsSurfaceElement() )
     csmp_error.notice( ERROR, "MeshManager<3>::ReplaceElementByFace", "element to be replaced is not a lower-dimensional surface element");
   if constexpr ( dim == 2 )
     if ( !eptr->IsLineElement() )
     csmp_error.notice( ERROR, "MeshManager<2>::ReplaceElementByFace", "element to be replaced is not a lower-dimensional line element");

   if ( eptr->FV() == nullptr )
     csmp_error.notice( INFO, "MeshManager<dim>::ReplaceElementByFace", "finite volume stencil pointer not initialised");
   if ( inner_eptr == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::ReplaceElementByFace", "pointer to higher dimensional element on inside not initialised");
   if ( inner_eptr == outer_eptr ) {
        csmp_error.notice( ERROR, "MeshManager<dim>::ReplaceElementByFace", "cannot create Face"
                                  "pointer to higher dimensional elements are the same");
        return nullptr;
     }
       
   assert( adjacent_face_of_inner_element < inner_eptr->Faces() );
   if ( outer_eptr != nullptr ) assert( adjacent_face_of_outer_element < outer_eptr->Faces() );

   // constructing new face
   const size_t face_id = faces_.size(); // since the face will be added at the end of the colony
   typename plf::colony<Face<dim>>::iterator
     fit = faces_.emplace( Face<dim>( *eptr, inner_eptr, outer_eptr,
                                       adjacent_face_of_inner_element, adjacent_face_of_outer_element,
                                       lvars, ivars ) );
   (*fit).Idx( face_id );

   // 3. deleting original Element
   elements_.erase( elements_.get_iterator(eptr) );
   eptr = nullptr;
   
   return &(*fit);
   
 } // end ReplaceElementByFace
       


/**
   Like ReplaceElementByFace, but without deletion of the original element.
*/
template<uint32_t dim>
Face<dim>* const MeshManager<dim>::ConstructFaceFromElement( const csmp::Element<dim>* const eptr,
                                                             csmp::Element<dim>* inner_eptr,
                                                             csmp::Element<dim>* outer_eptr,
                                                             uint32_t adjacent_face_of_inner_element,
                                                             uint32_t adjacent_face_of_outer_element,
                                                             const LocalVariables& lvars,
                                                             const IntegrationPointVariables& ivars )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( eptr == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::ConstructFaceFromElement", "pointer to lower-dimensional input element is not initialised");
     
   // is the element indeed lower dimensional?
   if constexpr ( dim == 3 )
     if ( !eptr->IsSurfaceElement() )
     csmp_error.notice( ERROR, "MeshManager<3>::ConstructFaceFromElement", "element to be replaced is not a lower-dimensional surface element");
   if constexpr ( dim == 2 )
     if ( !eptr->IsLineElement() )
     csmp_error.notice( ERROR, "MeshManager<2>::ConstructFaceFromElement", "element to be replaced is not a lower-dimensional line element");

   if ( inner_eptr == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::ConstructFaceFromElement", "pointer to higher dimensional element on inside not initialised");
   if ( inner_eptr == outer_eptr ) {
        csmp_error.notice( ERROR, "MeshManager<dim>::ConstructFaceFromElement", "cannot create Face"
                                  "pointer to higher dimensional elements are the same");
        return nullptr;
     }
     
   assert( adjacent_face_of_inner_element < inner_eptr->Faces() );
   if ( outer_eptr != nullptr ) assert( adjacent_face_of_outer_element < outer_eptr->Faces() );
       
   // constructing new face
   const size_t face_id = faces_.size(); // since the face will be added at the end of the colony
   typename plf::colony<Face<dim>>::iterator
     fit = faces_.emplace( Face<dim>( *eptr, inner_eptr, outer_eptr,
                                       adjacent_face_of_inner_element, adjacent_face_of_outer_element,
                                       lvars, ivars ) );
   (*fit).Idx( face_id );

   return &(*fit);
   
 } // end ConstructFaceFromElement

       
       
       
       
       
  /// optionally, the neighbor element pointers might not be assigned; @note node pointers must be supplied in CCW order from outside looking in
template<uint32_t dim>
Face<dim>* const MeshManager<dim>::AddFace( Element<dim>* const inner_parent, uint32_t inner_parent_face_id,
                                            Element<dim>* const outer_parent, uint32_t outer_parent_face_id,
                                            const LocalVariables& lvars,
                                            const IntegrationPointVariables& ivars )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( inner_parent == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddFace", "pointer to higher dimensional element on inside not initialised");
   if ( outer_parent == nullptr )
     csmp_error.notice( INFO, "MeshManager<dim>::AddFace", "pointer to higher dimensional element on ouside not initialised");
   if ( outer_parent == inner_parent ) {
        csmp_error.notice( INFO, "MeshManager<dim>::AddFace", "cannot create Face",
                                 "pointers to higher dimensional elements are both the same.");
        return nullptr;
     }
   // 1. constructing new face
   const size_t face_id = faces_.size();
   typename plf::colony<Face<dim>>::iterator
     fit = faces_.emplace( Face<dim>( fem_manager_, fvm_manager_, inner_parent, outer_parent,
                                      inner_parent_face_id, outer_parent_face_id, lvars, ivars ) );
   (*fit).Idx( face_id );

   return &(*fit);

} // end AddFace





/*
    Creates lower-dimensional (line-element) face between supplied faces.
*/
template<uint32_t dim>
Face<dim>* const MeshManager<dim>::AddEdgeFace( Face<dim>* const adjacent_face1,
                                                uint32_t parent_elmt1_segm_id,
                                                Face<dim>* const adjacent_face2,
                                                uint32_t parent_elmt2_segm_id,
                                                const LocalVariables& lvars,
                                                const IntegrationPointVariables& ivars,
                                                const std::vector<Node<dim>*>& nodes )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( adjacent_face1 == nullptr ) {
       csmp_error.notice( ERROR, "MeshManager<dim>::AddFace", "pointer to higher dimensional Face1 on inside not initialised");
       return nullptr;
    }
   if ( adjacent_face1 == adjacent_face2 ) {
       csmp_error.notice( ERROR, "MeshManager<dim>::AddFace", "Face pointers point to same Face", "Nothing was done");
       return nullptr;
     }

   // 1. if there is only one higher-dimensional neighbor for line element face
   if ( adjacent_face1->InnerParent() == adjacent_face2->InnerParent() ) {
        // 1.1 if one of the two faces has yet a different parent
        
        // 1.2 if there is only one parent
     }

   // 2. if there are two neighbors so that the face can be constructed without issues
   const size_t face_id = faces_.size();
   const CSMP_FEM_TYPE etype = adjacent_face1->FE()->ElementTypeOfFace(adjacent_face1->InnerParentFaceID());
   typename plf::colony<Face<dim>>::iterator
     fit = faces_.emplace( Face<dim>( fem_manager_.E(etype), fvm_manager_,
                                      adjacent_face1->InnerParent(), adjacent_face2->InnerParent(),
                                      parent_elmt1_segm_id, parent_elmt2_segm_id,
                                      nodes, lvars, ivars ) );
   // assigning a face number
   (*fit).Idx( face_id );

   return &(*fit);

} // end AddFace(from Face objects_)





template<uint32_t dim>
Face<dim>* const MeshManager<dim>::AddBoundaryFace( csmp::Element<dim>* const eptr,
                                                    uint32_t local_face_id,
                                                    const LocalVariables& lvars,
                                                    const IntegrationPointVariables& ivars )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( eptr == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddBoundaryFace", "element pointer not initialised");
   if ( local_face_id >= eptr->Faces() )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddBoundaryFace", "face ID does not exist in element");
   if ( eptr->Neighbor(local_face_id) != nullptr )
     csmp_error.notice( WARNING, "MeshManager<dim>::AddBoundaryFace", "element face has a neighbor; is it really located at model boundary?");

   // 1. constructing new face, connecting it to its higher-dimensional neighbor on the inside, and assigning nodes
   const size_t face_number{faces_.size()};
   typename plf::colony<Face<dim>>::iterator
     fit = faces_.emplace( Face<dim>( *eptr, fem_manager_.E( eptr->FE()->ElementTypeOfFace(local_face_id) ), fvm_manager_, local_face_id, lvars, ivars ) );

   (*fit).Idx( face_number );

   return &(*fit);

} // end AddBoundaryFace






/**
    This method is for connecting the matching Element faces of a node-matched split mesh as created, for instance by ANSYS.
    Apart from creating the InterFace, Node manifolds are created and-or updated as necessary.
*/
template<uint32_t dim>
InterFace<dim>*	const	MeshManager<dim>::AddInterFace( Element<dim>* const inner_parent, uint32_t inner_element_face_id,
                                                      Element<dim>* const outer_parent, uint32_t outer_element_face_id,
                                                      const LocalVariables& lvars,
                                                      const IntegrationPointVariables& ivars )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( inner_parent == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddInterFace", "pointer to higher dimensional element on inside not initialised");
   if ( outer_parent == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::AddInterFace", "pointer to higher dimensional element on ouside not initialised");
   if ( inner_parent == outer_parent ) {
        csmp_error.notice( ERROR, "MeshManager<dim>::AddInterFace", "cannot create InterFace",
                                  "inner and outer parent pointers are the same");
        return nullptr;
     }

   assert( inner_element_face_id < inner_parent->Faces() );
   assert( outer_element_face_id < outer_parent->Faces() );

   // 1. establishing the finite element type if the interface
   const CSMP_FEM_TYPE etype = inner_parent->FE()->ElementTypeOfFace( inner_element_face_id );

   // 2. constructing new interface
   const size_t iface_id = interfaces_.size();
   typename plf::colony<InterFace<dim>>::iterator
     ifp = interfaces_.emplace( InterFace<dim>( fem_manager_.E(etype), fvm_manager_.Stencil(etype), lvars, ivars ) );
   const bool assign_nodes{true};
   (*ifp).Assign( inner_parent, inner_element_face_id, outer_parent, outer_element_face_id, assign_nodes );
   (*ifp).Idx( iface_id );

   assert( node_manifold_manager_ != nullptr );

   // 3. assigning node manifolds
   const size_t n_nodes{(*ifp).FE()->Nodes()};
   for ( auto i{0}; i<n_nodes; ++i )
     // if the inside node is different from the outside node so that there needs to be a manifold
     if ( (*ifp).N(i,INSIDE) != (*ifp).N(i,OUTSIDE) ) {
         // 1. if both nodes are not yet manifolds
         if ( !(*ifp).N(i,INSIDE)->IsManifold() && !(*ifp).N(i,OUTSIDE)->IsManifold() ) {
              node_manifold_manager_->NewManifold( (*ifp).N(i,INSIDE), (*ifp).N(i,OUTSIDE), ManifoldType::INTERFACE );
              continue;
           }
         // 2. if both nodes are already manifolds, they are merged into single one
         if ( (*ifp).N(i,INSIDE)->IsManifold() && (*ifp).N(i,OUTSIDE)->IsManifold() ) {
              // if they are different from one-another, they are merged
              if ( (*ifp).N(i,INSIDE)->Manifold() != (*ifp).N(i,OUTSIDE)->Manifold() )
                node_manifold_manager_->MergeManifolds( (*ifp).N(i,INSIDE)->Manifold(),
                                                        (*ifp).N(i,OUTSIDE)->Manifold() );
              continue;
           }
         // 3. if the inside node is already a manifold and does not contain the second one
         //    because the second is not a manifold
         if ( (*ifp).N(i,INSIDE)->IsManifold() ) {
              // the outside node is added to it
              (*ifp).N(i,INSIDE)->Manifold()->Add( (*ifp).N(i,OUTSIDE), OUTSIDE, ManifoldType::INTERSECTION );
           }
         // 4. if the outside node is already a manifold
         else if ( (*ifp).N(i,OUTSIDE)->IsManifold() ) {
              // the inside node is added to the outside nodes manifold
              (*ifp).N(i,OUTSIDE)->Manifold()->Add( (*ifp).N(i,INSIDE), INSIDE, ManifoldType::INTERSECTION );
           }
       }
   
   return &(*ifp);

} // end AddInterFace





/**
    Constructs a new Interface. The former nodes of the Face become those on the inside of the interface.
    Duplicates the nodes identified by the vector, adding them to the required manifolds.
    Other nodes will be shared across the sides of the interface.
    
    @attention assumes that nodes have already been duplicated as necessary and manifolds have been created and are uptodate
*/
template<uint32_t dim>
InterFace<dim>* const MeshManager<dim>::ReplaceFaceByInterFace( csmp::Face<dim>* fptr,
                                                                const LocalVariables& lvars,
                                                                const IntegrationPointVariables& ivars )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );

   // 0. verifying the input
   // pointers
   if ( fptr == nullptr )
     csmp_error.notice( ERROR, "MeshManager<dim>::ReplaceFaceByInterFace", "Face pointer not initialised");

   // 1. constructing new interface
   const size_t iface_id = interfaces_.size();
   typename plf::colony<InterFace<dim>>::iterator
     ifp = interfaces_.emplace( InterFace<dim>( fptr->FE(), fptr->FV(), lvars, ivars ) );
   // 2. also assigning nodes, expecting that the nodes on the opposite side of the interface are already there
   (*ifp).Assign( fptr->InnerParent(), fptr->OuterParent() );
   (*ifp).Idx( iface_id );

#ifdef DEBUG
   // verifying that the nodes on the inside matching those of the face
   for ( auto i{0}; i<fptr->Nodes(); ++i ) {
        assert( (*ifp).N(i) != nullptr );
        assert( (*ifp).N(i,INSIDE) != nullptr );
        assert( fptr->N(i) == (*ifp).N(i,INSIDE) );
        assert( (*ifp).N(i,OUTSIDE) != nullptr );
     }
#endif

   // 3. removing original face
   faces_.erase( faces_.get_iterator( fptr ) );
   fptr = nullptr;

   return &(*ifp);

} // end ReplaceFaceByInterFace







/**
    Duplicates existing node inside of the MeshManager and connects it to corresponding manifold.
         @attention the current node is assumed to be on the INSIDE of the Interface; when there is none yet, this distinction will be made automatically
*/
template<uint32_t dim>
Node<dim>* const MeshManager<dim>::Duplicate( Node<dim>* const nptr_inside,
                                              INTERFACE_SIDE new_node_side,
                                              ManifoldType geometry )
  {
    if ( nptr_inside == nullptr )
      throw csmp::Exception( ERROR, "MeshManager<dim>::Duplicate", "Node does not exist.");

    // copying the inside node
    typename plf::colony<Node<dim>>::iterator
      nit = nodes_.emplace( Node<dim>( *nptr_inside ) );
    
    // creating or updating the NodeManifold
    if ( nptr_inside->IsManifold() ) {
         // if we are already dealing with a manifold, the geomtric classifier is retained
         nptr_inside->Manifold()->Add( &(*nit), new_node_side, nptr_inside->Manifold()->GeometricClassifier() );
         (*nit).Assign( nptr_inside->Manifold() );
      }
    else // a new manifold is created with the provided geometric classifier
      node_manifold_manager_->NewManifold( nptr_inside, &(*nit), geometry );

    return &(*nit);
  }









/**
   replaces supplied lower-dimensional elements with Face objects, establishing their connectivity; the Elements are deleted afterwards, setting input pointers to NULL
   
      the Node flags of the Element are used to determine whether this is a boundary face
      
     @note RANGE ERASE DOES ONLY WORK FOR A CONSECUTIVE RANGE OF ITERATORS WHERE it1 < it2
     @code
     elements_.erase( (*elmt_iterators.begin()), (*elmt_iterators.end()) );
     @endcode

*/
template<uint32_t dim>
vector<Face<dim>*>  MeshManager<dim>::ReplaceElementsByFaces( const PropertyDatabase<dim>& pref,
                                                              typename vector<Element<dim>*>::iterator first,
                                                              typename vector<Element<dim>*>::iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    const size_t       n_original_elements{ elements_.size() };
    vector<Face<dim>*> face_ptrs;
    const long         n_faces_to_build{ distance(first,last) };

    if ( n_faces_to_build == 0U ) {
         csmp_error.notice( WARNING, "MeshManager<dim>::ReplaceElementsByFaces", "supplied iterator range is empty; nothing was done.");
         return face_ptrs;
      }
    else face_ptrs.reserve( n_faces_to_build );
    
    const LocalVariables             lvars(pref.LocalVariablesAt(FACE));
    const IntegrationPointVariables& ivars(pref.IntegrationPointVariablesAt(FACE));
    
    // backup copy used later for the deletion
    typename vector<Element<dim>*>::iterator  erase_it{ first };
    
    // 1. converting Elements into Faces
    // ---------------------------------
    size_t face_idx{0};
    
    while( first != last )
      {
         // 1.1 initial checks and labeling
         // (input range must not contain any nullptrs)
         assert( (*first) != nullptr );
         if constexpr ( dim == 3 )
           if ( !(*first)->IsSurfaceElement() ) {
                (*first)->Out();
                throw csmp::Exception( ERROR, "MeshManager<3>::ReplaceElementsByFaces", "supplied element is not a surface element and cannot be converted to Face.");
             }
         if constexpr ( dim == 2 )
           if ( !(*first)->IsLineElement() ) {
                (*first)->Out();
                throw csmp::Exception( ERROR, "MeshManager<2>::ReplaceElementsByFaces", "supplied element is not a line element and cannot be converted to Face.");
             }
        
         // 1.2 simplified construction of Face at model boundary
         bool boundary_face{true};
         for ( auto nit=(*first)->NodesBegin(); nit!=(*first)->NodesEnd(); ++nit )
           if ( (*nit)->AtBoundary() == NOT ) {
                boundary_face = false;
                break;
             }
         if ( boundary_face ) { // finding higher dimensional neighbor and its face idx
              pair<Element<dim>* const,uint32_t> pelmt = parentElement<dim>( (*first)->NodesBegin(), (*first)->NodesEnd() );
              // creating Face, storing a pointer to it
              face_ptrs.push_back( AddBoundaryFace( pelmt.first, pelmt.second, lvars, ivars ) );
              // numbering new Face consecutively
              face_ptrs.back()->Idx( face_idx++ );
           }
           
         // 1.3 more involved construction of Face object in the interior of a model
         //    (both neighbors are present)
         else { // finding higher-dimensional neighbors (2)
              pair<Element<dim>*,Element<dim>*>  pelmts = parentElementsSharedByFace<dim>( (*first)->NodesBegin(), (*first)->NodesEnd() );
              assert( pelmts.first  != nullptr );
              assert( pelmts.second != nullptr );
              // finding the face numbers of the parent elements
              pair<uint32_t,uint32_t> face_ids = findAdjacentElementFaces( pelmts.first, pelmts.second );
              // creating Face, storing a pointer to it
              face_ptrs.push_back( ConstructFaceFromElement( (*first), pelmts.first, pelmts.second,
                                                             face_ids.first, face_ids.second, lvars, ivars ) );
              // numbering new Face consecutively
              face_ptrs.back()->Idx( face_idx++ );
           }
       
         // NOTE: no Element erasure yet because this would invalidate node parent vector, corrupting this functionality
         first++;
      }
      
     // 2. remove elements replaced by Face objects
     // -------------------------------------------
     // (no attention needs to be paid to neighbor connectivity because the whole lower dimensional regions will be removed)
     cout <<"\nMeshManager::ReplaceElementsByFaces: (n_elements="<< elements_.size();
     cout <<") deleting "<< n_faces_to_build <<" elements...\n";
     while ( erase_it != last )
       {
// POINTER REMAINS VALID AFTER DELETION:          cerr <<"\n\t"<< (*erase_it)->Idx();
          // null the element parent pointers in the elements nodes
          for ( auto i{0}; i<(*erase_it)->Nodes(); ++i )
            (*erase_it)->N(i)->Unassign( (*erase_it) );
          // get element pointer for colony and delete the element TODO: inefficient because this is done multiple times for each node
          elements_.erase( elements_.get_iterator( *erase_it ) );
          // set the supplied element pointer to null TODO: this needs to be communicated to pointers of input regions?
          (*erase_it) = nullptr;
          // increment iterator
          erase_it++;
       }
     assert( n_original_elements - elements_.size() == n_faces_to_build );
     
     // RANGE ERASE DOES ONLY WORK FOR A CONSECUTIVE RANGE OF ITERATORS WHERE it1 < it2
     //elements_.erase( (*elmt_iterators.begin()), (*elmt_iterators.end()) );

     // 3. cleaning up inter-CELL and node to parent connectivity
     // ---------------------------------------------------------
     // TODO: these are global changes! - do this only for nodes that are affected
     UpdateConnectivity();
     
     return face_ptrs;
     
 } // end ReplaceElementsByFaces


/*
#ifdef MESH_MANAGER_DEBUG
integrityCheck<dim,Element>( ElementsBegin(), ElementsEnd() );
if ( Faces() > 0 )
  integrityCheck<dim,Face>( FacesBegin(), FacesEnd() );
if ( InterFaces() > 0 ) {
      integrityCheck<dim,InterFace>( InterFacesBegin(), InterFacesEnd() );
     // add test for node manifolds
  }
#endif
*/

/*
     // 3. cleaning up the node to parent connectivity
     // ----------------------------------------------
     // TODO: these are global changes! - do this only for nodes that are affected
     for ( auto& nit : nodes_ ) {
          nit.EraseNullPointerParents(); // element parents
          nit.UpdateNeighbors();         // node neighbors
       }

     // 4. connecting Faces to each other (Face neighbors)
     // --------------------------------------------------
     if constexpr( dim == 3 ) BuildSurfaceConnectivity<Face>( face_ptrs.begin(), face_ptrs.end() );
     if constexpr( dim == 2 ) BuildLineConnectivity<Face>( face_ptrs.begin(), face_ptrs.end() );
*/






/**
      Sets neighbor pointers of cells surrounding the domain to 'nullptr' if they were pointing to cells within the domain.
      
      The intention of this method is to avoid that these pointers wil accidentially be derefefenced causing crashes
      once the subdomain has been deleted.
      
      @todo: appears to have side effects.
      
      SKM 9/2/2022
*/
template<uint32_t dim>
template<template<uint32_t> class CELL>
size_t MeshManager<dim>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<dim,CELL>& subdomain )
 {
 
 throw csmp::Exception( ERROR, "MeshManager<dim>::DetachOutsideNeighborsAlongPerimeter:", "Not properly tested yet, may have side effects" );
 
    size_t n_detachments{0};
    const size_t n_cells{ subdomain.Elements() };
    for ( size_t i=subdomain.InteriorElements(); i < n_cells; ++i ) {
         const uint32_t n_perim_faces( subdomain.PerimeterFaces(i) );
         for ( auto j{0}; j < n_perim_faces; ++j ) {
              auto p_face = subdomain.PerimeterFace( i, j );
              // detaching outside neighbor, if any
              if ( subdomain.E(i)->Neighbor(p_face) != nullptr ) {
                   const auto n_nbor_nbors{ subdomain.E(i)->Neighbor(p_face)->Neighbors() };
                   for ( auto k{0}; k<n_nbor_nbors; ++k )
                     if ( subdomain.E(i)->Neighbor(p_face)->Neighbor(k) == subdomain.E(i) ) {
                          // detach subdomain cell
                          subdomain.E(i)->Neighbor(p_face)->Neighbor(k)->Unassign( subdomain.E(i) );
                          n_detachments++;
                       }
                }
           }
      }
      
   return n_detachments;
  
 } // end DetachOutsideNeighborsAlongPerimeter
  
template size_t MeshManager<1>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<1,Element>& );
template size_t MeshManager<2>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<2,Element>& );
template size_t MeshManager<3>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<3,Element>& );

template size_t MeshManager<1>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<1,Face>& );
template size_t MeshManager<2>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<2,Face>& );
template size_t MeshManager<3>::DetachOutsideNeighborsAlongPerimeter( ModelSubDomain<3,Face>& );



// ========================================================================================================================

// ERASURES

// ========================================================================================================================




/**
    Erases the supplied sequence of nodes in the MeshManager returning the number of erasures.
    The pointer to the deleted elements are set to 'nullptr'.
    Any potential node manifolds are updated.
*/
template<uint32_t dim>
size_t MeshManager<dim>::Delete( typename vector<Node<dim>*>::iterator first,
                                 typename vector<Node<dim>*>::iterator last )
 {
    size_t deleted_nodes( distance(first,last) );
 
     if ( deleted_nodes == 0U ) return 0U;
     
     // deleting the objects that are stored in the MeshManager's 'elements' container
     while( first != last ) {
          // only couting active nodes
          if ( (*first) == nullptr ) deleted_nodes--;
          // disconnecting the node from its manifold
          else if ( (*first)->IsManifold() ){
               (*first)->Manifold()->Remove( (*first) );
               // if no manifold remains it is removed from the manager
               if ( (*first)->Manifold()->Branches() == 1U )
               node_manifold_manager_->Delete( (*first)->Manifold() );
            }
          nodes_.erase( nodes_.get_iterator(*first) );
          (*first) = nullptr;
          first++;
       }
     
     return deleted_nodes;
    
 } // end Delete






/**
    Erases range of elements not counting  null-pointer cells in the supplied sequence,
    returning the number of genuine erasures.
    
    Also checks whether element deletion causes orphan nodes. If so, these are deleted as well.
    This check involves counting the nodes parent elements that are not null pointers to make
    sure that the true state of the node is captured.
    
    @attention The connectivity of the affected mesh neighborhood needs to get fixed separately.
    Inside all cell destructors the following steps are performed:

    1. set the neighbor pointers to the element to zero (= disconnect the neighbor elements)
    2. (remove the pointers from the connected nodes to this parent-element) - done later sweeping over the nodes
    3. delete orphanaged nodes
    4. disconnect nodes
    5. delete element
    
*/
template<uint32_t dim>
size_t MeshManager<dim>::Delete( typename vector<Element<dim>*>::iterator first,
                                 typename vector<Element<dim>*>::iterator last )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     size_t deleted_elements( distance(first,last) );
 
     if ( deleted_elements == 0U ) return 0U;
     
     // 1. deleting the objects that are stored in the MeshManager's 'elements' container
     vector<Node<dim>*>  node_ptrs;
     size_t              n_deleted_nodes{0};
     
     node_ptrs.reserve( Nodes() );
     
     while( first != last ) {
          // checking the pointer
          if ( (*first) == nullptr )  {
               csmp_error.notice( WARNING, "MeshManager::Delete",
                                           "supplied iterator dereferences to a nullptr" );
               deleted_elements--;
            }
          else {
               // remembering the nodes of the element
               const auto nodesEnd{ (*first)->NodesEnd() };
               for ( auto nit=(*first)->NodesBegin(); nit!=nodesEnd; ++nit )
                 node_ptrs.push_back( (*nit) );
               // deleting the element, setting the pointer in the input iterator to zero
               auto colony_it{ elements_.get_iterator( const_cast<Element<dim>*>(*first) ) };
               assert( colony_it!= elements_.end() );
               elements_.erase( colony_it );
               (*first) = nullptr;
               first++;
            }
       }

     // 2. deleting potential dangling nodes before the element
     sort( node_ptrs.begin(), node_ptrs.end() );
     node_ptrs.erase( unique( node_ptrs.begin(), node_ptrs.end() ), node_ptrs.end() );
  
     for ( auto& nit : node_ptrs )
       {
          // counting the parent elements of the nodes which have survived
          // the element deletion
          const auto n_parents{nit->Parents()};
          size_t n_active_parents{0};
          for ( auto j{0}; j<n_parents; ++j )
           if ( nit->Parent(j) != nullptr )
             n_active_parents++;

         // if this is indeed a node that has no parent element left
         // it is removed after removing it from potential manifolds
         if ( n_active_parents == 0 ) {
              if ( !interfaces_.empty() && nit->IsManifold() )
                nit->Manifold()->Remove( nit );
              nodes_.erase( nodes_.get_iterator(nit) );
              nit = nullptr;
              n_deleted_nodes++;
           }
         else if ( n_active_parents < n_parents )  {
              // 3. remove null pointers to parent elements that do not exist anymore
              nit->EraseNullPointerParents();
           }
       }

    if ( n_deleted_nodes > 0 )
      csmp_error.notice( INFO, "MeshManager::Delete", "Element deletion did create orphan nodes which were deleted." );

    return deleted_elements;
    
 } // end Delete




template<uint32_t dim>
size_t MeshManager<dim>::Delete( typename vector<Face<dim>*>::iterator first,
                                 typename vector<Face<dim>*>::iterator last )
 {
    size_t deleted_faces( distance(first,last) );
 
     if ( deleted_faces == 0U ) return 0U;
     
     // deleting the objects that are stored in the MeshManager's 'elements' container
     while( first != last ) {
          if ( (*first) == nullptr ) {
                deleted_faces--;
                faces_.erase( faces_.get_iterator(*first) );
               (*first) = nullptr;
            }
          first++;
       }
     
     return deleted_faces;
    
 } // end Delete




template<uint32_t dim>
size_t MeshManager<dim>::Delete( typename vector<InterFace<dim>*>::iterator first,
                                 typename vector<InterFace<dim>*>::iterator last )
 {
    size_t deleted_ifaces( distance(first,last) );
 
     if ( deleted_ifaces == 0U ) return 0U;
     
     // deleting the objects that are stored in the MeshManager's 'elements' container
     while( first != last ) {
          if ( (*first) == nullptr ) {
              deleted_ifaces--;
              interfaces_.erase( interfaces_.get_iterator(*first) );
              (*first) = nullptr;
            }
          first++;
       }

     return deleted_ifaces;
    
 } // end Delete











/**
    Connects elements, faces or interfaces with their equidimensional neighbors
    
    @todo deal with potential manifolds relating to lower-dimensional cells
    
*/
template<uint32_t dim>
template<template<uint32_t> class CELL>
void  MeshManager<dim>::BuildConnectivity( typename vector<CELL<dim>*>::const_iterator first,
                                           typename vector<CELL<dim>*>::const_iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( distance(first,last) == 0U ) {
         csmp_error.notice( WARNING, "MeshManager<dim>::BuildConnectivity:", "supplied cell vector is empty; nothing was done." );
         return;
      }
    //cout << "\nMeshManager<"<< dim <<">::BuildConnectivity: Establishing CSMP FE neighbor connectivity...\n";
    
    // making subranges for the volume, surface, and line elements
    vector<CELL<dim>*> volume_cells, surface_cells, line_cells;
    const size_t n_cells_max = distance(first,last);
    const auto   cellsEnd{last};
    
    // if we are dealing with element connectivity in 3D
    if constexpr ( dim == 3 && is_same< CELL<dim>,Element<dim> >::value ) {
        volume_cells.reserve( n_cells_max );
        surface_cells.reserve( n_cells_max/3 );
        line_cells.reserve( n_cells_max/6 );
        while ( first != cellsEnd ) {
             if      ( (*first)->FE()->IsVolumeElement() )  volume_cells.push_back(*first);
             else if ( (*first)->FE()->IsSurfaceElement() ) surface_cells.push_back(*first);
             else if ( (*first)->FE()->IsLineElement() )    line_cells.push_back(*first);
             first++;
          }
        // connecting the elements found
        if ( !volume_cells.empty() )  BuildVolumeConnectivity<Element>( volume_cells.begin(), volume_cells.end() );
        if ( !surface_cells.empty() ) BuildSurfaceConnectivity<Element>( surface_cells.begin(), surface_cells.end() );
        if ( !line_cells.empty() )    BuildLineConnectivity<Element>( line_cells.begin(), line_cells.end() );
      }
      
    // else
    if constexpr ( is_same< CELL<dim>,Face<dim> >::value || is_same< CELL<dim>,InterFace<dim> >::value  ) {
        // 2D case
        if constexpr ( dim == 3 || dim == 2 ) {
             surface_cells.reserve( n_cells_max );
             line_cells.reserve( n_cells_max/6 );
             while ( first != cellsEnd ) {
                  if ( (*first)->FE()->IsSurfaceElement() )   surface_cells.push_back(*first);
                  else if ( (*first)->FE()->IsLineElement() ) line_cells.push_back(*first);
                  first++;
               }
             // connecting the elements found
             if ( !surface_cells.empty() ) BuildSurfaceConnectivity<CELL>( surface_cells.begin(), surface_cells.end() );
             if ( !line_cells.empty() )    BuildLineConnectivity<CELL>( line_cells.begin(), line_cells.end() );
          }
        // 1D case
        if constexpr ( dim == 1 ) {
             line_cells.reserve( n_cells_max );
             while ( first != cellsEnd ) {
                  line_cells.push_back(*first);
                  first++;
               }
          }
      }
    
 } // end BuildConnectivity

template void MeshManager<3>::BuildConnectivity<Element>( typename vector<Element<3>*>::const_iterator, typename vector<Element<3>*>::const_iterator );
template void MeshManager<3>::BuildConnectivity<Face>( typename vector<Face<3>*>::const_iterator, typename vector<Face<3>*>::const_iterator );
template void MeshManager<3>::BuildConnectivity<InterFace>( typename vector<InterFace<3>*>::const_iterator, typename vector<InterFace<3>*>::const_iterator );

template void MeshManager<2>::BuildConnectivity<Element>( typename vector<Element<2>*>::const_iterator, typename vector<Element<2>*>::const_iterator );
template void MeshManager<2>::BuildConnectivity<Face>( typename vector<Face<2>*>::const_iterator, typename vector<Face<2>*>::const_iterator );
template void MeshManager<2>::BuildConnectivity<InterFace>( typename vector<InterFace<2>*>::const_iterator, typename vector<InterFace<2>*>::const_iterator );

template void MeshManager<1>::BuildConnectivity<Element>( typename vector<Element<1>*>::const_iterator, typename vector<Element<1>*>::const_iterator );
template void MeshManager<1>::BuildConnectivity<Face>( typename vector<Face<1>*>::const_iterator, typename vector<Face<1>*>::const_iterator );
template void MeshManager<1>::BuildConnectivity<InterFace>( typename vector<InterFace<1>*>::const_iterator, typename vector<InterFace<1>*>::const_iterator );



/**
    Element connectivity because volumetric elements never are faces or interfaces.
*/

template<>
template<template<uint32_t> class CELL>
void MeshManager<3>::BuildVolumeConnectivity( typename std::vector<CELL<3U>*>::const_iterator first,
                                              typename std::vector<CELL<3U>*>::const_iterator last )
   {
      // this method applies only to volumetric elements and 3D
      if constexpr ( is_same< CELL<3>,Element<3> >::value ) {
           // creating search keys from the corner nodes of the element faces
           // corner-nodes      elements that share face and their face id
           map<set<Node<3>*>,map<Element<3>*,uint32_t> >  elmt_pairs;
           const auto                                     elementsEnd{last};
           
           // pairing the elements up in the search map
           while ( first != elementsEnd ) {
                assert( (*first) != nullptr );
                const size_t n_faces{ (*first)->Faces() };
                for ( auto face{0}; face < n_faces; ++face ) {
                     // trying to insert it into the map
                     auto it = elmt_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<Element<3>*,uint32_t>{make_pair(*first,face)} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<Element<3>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the elements to one another
           for ( auto& it : elmt_pairs ) {
                size_t n_face_nbors{ it.second.size() };
                if ( n_face_nbors == 2 ) {
                     Element<3>* const eptr1 = (*it.second.begin()).first;
                     Element<3>* const eptr2 = (*it.second.rbegin()).first;
                     const uint32_t face_e1  = (*it.second.begin()).second;
                     const uint32_t face_e2  = (*it.second.rbegin()).second;
                     eptr1->Assign( face_e1, eptr2 );
                     eptr2->Assign( face_e2, eptr1 );
                  }
                // else no assignments have to be made as there is no neighbor
                assert( n_face_nbors <= 2 );
             }
        }
   
   } // end BuildVolumeElementConnectivity

template void MeshManager<3>::BuildVolumeConnectivity<Element>( typename vector<Element<3>*>::const_iterator,
                                                                typename vector<Element<3>*>::const_iterator );





/**
    Connects neighboring cells, disambiguating the potential manifolds by choosing co-planar elements.
    
    This is accomplished by finding the angle between 2 suface elements in 3D,  returning the acute angle in degrees (0..90o).
*/
template<uint32_t dim>
template<template<uint32_t> class CELL>
void MeshManager<dim>::BuildSurfaceConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                                                 typename std::vector<CELL<dim>*>::const_iterator last )
   {
      // this method applies only to surface elements in 3D
      if constexpr ( dim == 3 ) {
           // creating search keys from the corner nodes of the element faces
           // corner-nodes      elements that share face and their face id
           map<set<Node<3>*>,map<CELL<3>*,uint32_t> >  elmt_pairs;

           // pairing the cells up in the search map
           const auto cellsEnd{last};
           while ( first != cellsEnd ) {
                assert( (*first) != nullptr );
                const size_t n_faces{ (*first)->Faces() };
                for ( auto face{0}; face < n_faces; ++face ) {
                     // trying to insert cell into the map using a search key of node pointers
                     pair<typename map<set<Node<3>*>,map<CELL<3>*,uint32_t> >::iterator,bool>
                       it = elmt_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<CELL<3>*,uint32_t>{{*first,face}} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<CELL<3>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the cells to one another
           for ( auto& it : elmt_pairs ) {
                const size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2 ) {
                     CELL<3>* const ptr1  = (*it.second.begin()).first;
                     CELL<3>* const ptr2  = (*it.second.rbegin()).first;
                     assert( ptr1 != nullptr );
                     assert( ptr2 != nullptr );
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is a manifold and two most suitable neighbors must be found
                else if ( n_face_nbors > 2 ) { // TODO: test
                     // finding all possible combinations of surface elements
                     vector<int64_t> sequence( n_face_nbors );
                     iota( sequence.begin(), sequence.end(), 0 ); // fill 0..n-1
                     const size_t            n_samples{2};
                     deque<vector<int64_t> >  combinations;
                     const size_t n_combinations = createUniqueCombinations( sequence, n_samples, combinations );
                     // finding the combination of surfaces or line elements with the smallest acute angle between them
                     map<double,size_t>  ordered_combinations;
                     for ( auto i{0}; i < n_combinations; ++i ) {
                          CELL<3>* const ptr1 = (*next(it.second.begin(),combinations[i][0])).first;
                          CELL<3>* const ptr2 = (*next(it.second.begin(),combinations[i][1])).first;
                          assert( ptr1 != nullptr );
                          assert( ptr2 != nullptr );
                          const double angle = ( ptr1->IsSurfaceElement() && ptr2->IsSurfaceElement() ) ?
                                                   angleBetweenSurfaceCells( ptr1, ptr2 ) : angleBetweenLineCells( ptr1, ptr2 );
                          // using smallest angle
                          const double acute_angle = (angle > 90.) ? 180. -angle : angle;
                          // ordering
                          ordered_combinations.insert( make_pair(acute_angle,i) );
                       }
                     // the first element in the map has the smallest angle
                     const size_t combi   = (*ordered_combinations.begin()).second;
                     CELL<3>* const ptr1  = (*next(it.second.begin(),combinations[combi][0])).first;
                     CELL<3>* const ptr2  = (*next(it.second.begin(),combinations[combi][1])).first;
                     assert( ptr1 != nullptr );
                     assert( ptr2 != nullptr );
                     const uint32_t face_e1 = (*next(it.second.begin(),combinations[combi][0])).second;
                     const uint32_t face_e2 = (*next(it.second.begin(),combinations[combi][1])).second;
                     // uff! - finally.
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else no assignments have to be made as there is no neighbor
             }
        }
 
      // for surface elements, faces or interfaces in a 2D model
      if constexpr ( dim == 2 ) {
           map<set<Node<2>*>,map<CELL<2>*,uint32_t> > elmt_pairs;
           const auto                                 cellsEnd{last};
           
           // pairing the elements up in the search map
           while ( first != cellsEnd ) {
                assert( (*first) != nullptr );
                const size_t n_faces{ (*first)->Faces() };
                for ( auto face{0}; face < n_faces; ++face ) {
                     // trying to insert it into the map
                     auto it = elmt_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<CELL<2>*,uint32_t>{make_pair(*first,face)} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<CELL<2>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the elements to one another
           for ( auto& it : elmt_pairs ) {
                size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2 ) {
                     CELL<2>* const ptr1  = (*it.second.begin()).first;
                     CELL<2>* const ptr2  = (*it.second.rbegin()).first;
                     assert( ptr1 != nullptr );
                     assert( ptr2 != nullptr );
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is a manifold and two most suitable neighbors must be found
             }
        }
 
 } // end BuildSurfaceElementConnectivity
                                        
template void MeshManager<3>::BuildSurfaceConnectivity<Element>( typename vector<Element<3>*>::const_iterator,
                                                                 typename vector<Element<3>*>::const_iterator );
template void MeshManager<2>::BuildSurfaceConnectivity<Element>( typename vector<Element<2>*>::const_iterator,
                                                                 typename vector<Element<2>*>::const_iterator );
template void MeshManager<1>::BuildSurfaceConnectivity<Element>( typename vector<Element<1>*>::const_iterator,
                                                                typename vector<Element<1>*>::const_iterator );

template void MeshManager<3>::BuildSurfaceConnectivity<Face>( typename vector<Face<3>*>::const_iterator,
                                                              typename vector<Face<3>*>::const_iterator );
template void MeshManager<2>::BuildSurfaceConnectivity<Face>( typename vector<Face<2>*>::const_iterator,
                                                              typename vector<Face<2>*>::const_iterator );
template void MeshManager<1>::BuildSurfaceConnectivity<Face>( typename vector<Face<1>*>::const_iterator,
                                                              typename vector<Face<1>*>::const_iterator );

template void MeshManager<3>::BuildSurfaceConnectivity<InterFace>( typename vector<InterFace<3>*>::const_iterator,
                                                                   typename vector<InterFace<3>*>::const_iterator );
template void MeshManager<2>::BuildSurfaceConnectivity<InterFace>( typename vector<InterFace<2>*>::const_iterator,
                                                                   typename vector<InterFace<2>*>::const_iterator );
template void MeshManager<1>::BuildSurfaceConnectivity<InterFace>( typename vector<InterFace<1>*>::const_iterator,
                                                                   typename vector<InterFace<1>*>::const_iterator );


 // TESTING
 /*
cerr <<"\nFace connectivity at boundary:\n";
while ( first2 != cellsEnd ) {
     cerr <<"\n\t"<< (*first2)->Idx() <<": ";
     for ( size_t face{0}; face < (*first2)->Faces(); ++face ) {
          if ( (*first2)->Neighbor(face) == nullptr ) cerr << face <<":NO ";
          else cerr << face <<":"<< (*first2)->Neighbor(face)->Idx() <<" ";
       }
     first2++;
  }
*/




/**
   Line element manifolds exist in 3D and 2D.
*/
template<uint32_t dim>
template<template<uint32_t> class CELL>
void MeshManager<dim>::BuildLineConnectivity( typename std::vector<CELL<dim>*>::const_iterator first,
                                              typename std::vector<CELL<dim>*>::const_iterator last )
   {
      // this method applies only to line elements in 2 and 3D
      if constexpr ( dim != 1 ) {
           // creating search keys from the corner nodes of the element faces
           // corner-nodes      elements that share face and their face id
           map<set<Node<dim>*>,map<CELL<dim>*,uint32_t> >  elmt_pairs;
           const auto                                      elementsEnd{last};
           
           // pairing the elements up in the search map
           while ( first != elementsEnd ) {
                assert( (*first) != nullptr );
                const size_t n_faces{ (*first)->Faces() };
                for ( auto face{0}; face < n_faces; ++face ) {
                     // trying to insert it into the map
                     auto it = elmt_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<CELL<dim>*,uint32_t>{make_pair(*first,face)} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<CELL<dim>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the cells to one another
           for ( auto& it : elmt_pairs ) {
                size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2 ) {
                     CELL<dim>* const ptr1  = (*it.second.begin()).first;
                     CELL<dim>* const ptr2  = (*it.second.rbegin()).first;
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is a manifold and two most suitable neighbors must be found
                else if ( n_face_nbors > 2 ) { // TODO: test
                     // finding all possible combinations of surface elements
                     vector<int64_t> sequence( n_face_nbors );
                     iota( sequence.begin(), sequence.end(), 0 ); // fill 0..n-1
                     const size_t            n_samples{2};
                     deque<vector<int64_t> >  combinations;
                     const size_t n_combinations = createUniqueCombinations( sequence, n_samples, combinations );
                     // finding the combination of surfaces with the smallest acute angle between them
                     map<double,size_t>  ordered_combinations;
                     for ( auto i{0}; i < n_combinations; ++i ) {
                          CELL<dim>* const ptr1 = (*next(it.second.begin(),combinations[i][0])).first;
                          CELL<dim>* const ptr2 = (*next(it.second.begin(),combinations[i][1])).first;
                          const double angle = angleBetweenLineCells( ptr1, ptr2 );
                          const double acute_angle = (angle > 90.) ? 180. -angle : angle;
                          // ordering
                          ordered_combinations.insert( make_pair(acute_angle,i) );
                       }
                     // the first element in the map has the smallest angle
                     const size_t combi    = (*ordered_combinations.begin()).second;
                     CELL<dim>* const ptr1 = (*next(it.second.begin(),combinations[combi][0])).first;
                     CELL<dim>* const ptr2 = (*next(it.second.begin(),combinations[combi][1])).first;
                     const uint32_t face_e1  = (*next(it.second.begin(),combinations[combi][0])).second;
                     const uint32_t face_e2  = (*next(it.second.begin(),combinations[combi][1])).second;
                     // uff! - finally.
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else no assignments have to be made as there is no neighbor
             }
        }
 
      // for line elements, faces or interfaces in a 1D model
      if constexpr ( dim == 1 ) {
           map<set<Node<1>*>,map<CELL<1>*,uint32_t> >  elmt_pairs;
           const auto                                  elementsEnd{last};
           
           // pairing the elements up in the search map
           while ( first != elementsEnd ) {
                assert( (*first) != nullptr );
                const size_t n_faces{ (*first)->Faces() };
                for ( auto face{0}; face < n_faces; ++face ) {
                     // trying to insert it into the map
                     auto it = elmt_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<CELL<1>*,uint32_t>{make_pair(*first,face)} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<CELL<1>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the elements to one another
           for ( auto& it : elmt_pairs ) {
                size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2 ) {
                     CELL<1>* const ptr1  = (*it.second.begin()).first;
                     CELL<1>* const ptr2  = (*it.second.rbegin()).first;
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is a manifold and two most suitable neighbors must be found
             }
        }
        
   } // end BuildLineElementConnectivity


template void MeshManager<3>::BuildLineConnectivity<Element>( typename vector<Element<3>*>::const_iterator,
                                                              typename vector<Element<3>*>::const_iterator );
template void MeshManager<2>::BuildLineConnectivity<Element>( typename vector<Element<2>*>::const_iterator,
                                                              typename vector<Element<2>*>::const_iterator );

template void MeshManager<3>::BuildLineConnectivity<Face>( typename vector<Face<3>*>::const_iterator,
                                                           typename vector<Face<3>*>::const_iterator );
template void MeshManager<2>::BuildLineConnectivity<Face>( typename vector<Face<2>*>::const_iterator,
                                                           typename vector<Face<2>*>::const_iterator );

template void MeshManager<3>::BuildLineConnectivity<InterFace>( typename vector<InterFace<3>*>::const_iterator,
                                                                typename vector<InterFace<3>*>::const_iterator );
template void MeshManager<2>::BuildLineConnectivity<InterFace>( typename vector<InterFace<2>*>::const_iterator,
                                                                typename vector<InterFace<2>*>::const_iterator );







/**
    After some diagnostics that establish the extent of mesh modification, the cell neighbor connectivity is rebuilt.
    
    @attention calls EraseNullPointerCells()
    
    TODO: rewrite this in a form that makes use of existing connectivity.
*/
template<uint32_t dim>
void MeshManager<dim>::UpdateConnectivity()
 {
    // creating pointer vectors that are needed by BuildConnectivity() methods
    vector<Element<dim>*> element_ptrs;
    element_ptrs.reserve( elements_.size() );
    for ( auto& it : elements_ ) element_ptrs.push_back( &it );
    BuildConnectivity<csmp::Element>( element_ptrs.begin(), element_ptrs.end() );
    element_ptrs.clear();
    
    if ( !faces_.empty() ) {
         vector<Face<dim>*> face_ptrs;
         face_ptrs.reserve( faces_.size() );
         for ( auto& it : faces_ ) face_ptrs.push_back( &it );
         BuildConnectivity<csmp::Face>( face_ptrs.begin(), face_ptrs.end() );
      }
    if ( !interfaces_.empty() ) {
         vector<InterFace<dim>*> iface_ptrs;
         iface_ptrs.reserve( faces_.size() );
         for ( auto& it : interfaces_ ) iface_ptrs.push_back( &it );
         BuildConnectivity<csmp::InterFace>( iface_ptrs.begin(), iface_ptrs.end() );
      }
    
    // global node connectivity to parent elements
    // RebuildNodeParentElementRelationships( elements_.begin(), elements_.end() );
    // ----------------------------------------------------------------------------
    // TODO: restrict this to the neighborhood where changes occurred

    // counting the parent elements of each node
    map<Node<dim>*,set<Element<dim>*> >  parent_elmts_per_node;
    for ( auto& it : elements_ ) {
        const auto nodes_end{it.NodesEnd()};
        for ( auto nit = it.NodesBegin(); nit != nodes_end; ++nit ) {
             pair<typename map<Node<dim>*,set<Element<dim>*> >::iterator,bool>
               mit = parent_elmts_per_node.insert( make_pair( (*nit), set<Element<dim>*>{ &it } ) );
             if ( mit.second == false )
               (*mit.first).second.insert( &it );
          }
      }

    // reserving the memory for the parent storage and assigning the parent elements
    for ( auto& n : parent_elmts_per_node ) {
         const auto n_parents = static_cast<uint32_t>( n.second.size() );
         n.first->ResizeParentStorage( n_parents );
         // loopin over the future parents
         for ( const auto& it : n.second ) {
           size_t n_nodes{it->Nodes()};
           // assigning them to the node
           for ( auto j=0U; j<n_nodes; ++j )
             if ( n.first == it->N(j) ) {
                 it->N(j)->Assign( j, it );
                 break;
              }
           }
         assert( n.first->Parents() >= 1 );
      }
    // TODO: extend the repairs to include NodeManifolds
    cerr <<"\nMeshManager::UpdateConnectivity: WARNING: node manifolds are not reestablished here yet.\n";
    
 } // end UpdateConnectivity








/**
      Converts supplied range of faces into interfaces. Prior to InterFace creation, perimeter Faces of the Face region are identified, receiving a special treatment that
      avoids that the nodes along the tiplines are duplicated, see source.
      
      Supplied Face objects are replaced by InterFace objects (Face objects are deleted); returns pointer to first perimeter interface.
  Perimeter nodes are not duplicated, see CSMP User' s guide/

 @author SKM modified original function
 @date 16/05/2021
 
Duplicates nodes as necessary. 
@note Originally this was the Split() function from JC 1/7/2019

TODO: needs complete rewrite

*/







/** (Re)number all cells; either continuous (in_a_single_sequence=true): Elements then Faces then InterFaces) or in seperate 0..n-1 ranges for the different entity types.

@note this member function is constant because the idx_ is a mutable variable in the cell classes

*/
template<uint32_t dim>
void MeshManager<dim>::AssignUniqueNumbers( bool in_a_single_sequence )
{
  size_t n( 0U );
  for_each( nodes_.begin(), nodes_.end(), [&n]( Node<dim>& o ) { o.Idx( n++ ); return o; } );

  n = 0U; // resetting the counter
  for_each( elements_.begin(), elements_.end(), [&n]( Element<dim>& o ) { o.Idx( n++ ); return o; } );

  if ( !in_a_single_sequence ) n = 0U;
  for_each( faces_.begin(), faces_.end(), [&n]( Face<dim>& o ) { o.Idx( n++ ); return o; } );

  if ( !in_a_single_sequence ) n = 0U;
  for_each( interfaces_.begin(), interfaces_.end(), [&n]( InterFace<dim>& o ) { o.Idx( n++ ); return o; } );

} // end AssignUniqueNumbers




/**
    Puts nodes, elements, faces, and interfaces into an ascending order given by Idx() variables
    
    @attention const because it just changes the order but not the content of the container.
*/
template<uint32_t dim>
void MeshManager<dim>::ReorderObjectsByIndexes()
 {
    nodes_.sort( []( const auto& a, const auto& b ) -> bool
                    {
                        // (<) implies in ascending order
                        return a.Idx() < b.Idx();
                    } );

    elements_.sort( []( const auto& a, const auto& b ) -> bool
                      {
                          return a.Idx() < b.Idx();
                      } );
    
     if ( !faces_.empty() )
       faces_.sort( []( const auto& a, const auto& b ) -> bool
                      {
                         return a.Idx() < b.Idx();
                      } );

     if ( !interfaces_.empty() )
       interfaces_.sort( []( const auto& a, const auto& b ) -> bool
                          {
                            return a.Idx() < b.Idx();
                          } );

 } // end ReorderObjectsByIndexes





/** Writes VData information from MeshManager to VSet

Writes the elements, faces and interfaces
stored in the current MeshManager
to the supplied VSet.
This also includes the connectivity information,
i.e. the connections between these entities.

@param vset A VSet preferably empty, which is
resized first, if necessary and into which the
MeshManager connectivity information is input.

@param get_indices_from_stored_variables  renumber_uniquely = true is default,
else the current indexing of nodes, elements etc. is used.

@section conventions Conventions

Connections will exist between elements of the
topology, i.e., between:

- volume elements
- surface elements
- line elements

these connections are made via pointers and exist
across the faces, edges, and end end point of these entities.

Where there is no neighbor, like for a face that sits on a model
boundary, the neighbor pointer will be a 'NULL'.

Special rules apply to Face and Interface objects that
are essentially elements, but they are lower-dimensional
(for instance surfaces in a volumetric model)

Faces and Interfaces also have extra connections
to their higher-dimensional "parent" elements.

All these entities are stored in the VSet in a specific order:

1. highest dimensional elements (volumes in 3D)
2. dim-1 elements
3. dim-2 elements, if any
4. Face objects
5. InterFace objects

Since pointers cannot be written to file as such,
the first step in the output is to give all the entities
a unique numbering: 0..n-1 elements, then faces, last interfaces.

@section application Application

The data of the MeshManager is output to a VSet, for instance, when the
whole model is to be saved to a binary file.

@attention If the VSet is initialized successfully
the method will report this.

SKM 15/9/2014: fixed major bug in the output code blocks for
plist and pfverts.

@attention renumber the elements of the mesh first, using

@code
// renumbering the mesh; face numbers follow those of the elements
// last are the interface objects
// NB: we remember how many faces and interfaces there were to deduce the offsets at a later point
const bool in_a_single_sequence(true);
AssignUniqueNumbers( in_a_single_sequence );
@endcode

*/
template<uint32_t dim>
void MeshManager<dim>::OutputMeshTo( VSet<dim>& vset, bool get_indices_from_stored_variables )
{
  ErrorHandler& csmp_error( ErrorHandler::Instance() );
  
  if ( elements_.empty() ) {
       csmp_error.notice( WARNING, "MeshManager<dim>::OutputMeshTo", "Mesh Manager does not currently store a mesh; no output");
       return;
    }

  // 0. creating a mapping from pointers to integer values
  // -----------------------------------------------------
  if ( !get_indices_from_stored_variables ) {
      const bool in_a_single_sequence(true);
      AssignUniqueNumbers( in_a_single_sequence );
    }
  // else, MeshManager relies on current indexing of nodes, elements etc.
  // if this is non-unique or non-contiguous, this method call will fail.
  else ReorderObjectsByIndexes();

  // 1. resizing the VSet
  // --------------------
  if ( HybridElementMesh() || faces_.size() > 0 || interfaces_.size() > 0 )
    {
      const size_t higherDimParents( 2U );
      const size_t higherDimParentsFaceNum( 2U );
      const size_t interfaceMultiplier( 2U );
      const size_t interfaceExtras( 1U ); // 1 entry for potential high dim element

      deque<uint32_t>  nodes_per_element;
      deque<uint32_t>  neighbors_per_element;
      deque<int8_t>  csmp_fem_types;

      // 1.1 identifying how many nodes and neighbors there are per element
      for ( const auto& e : elements_ ) {
          nodes_per_element.push_back( e.Nodes() );
          neighbors_per_element.push_back( e.Neighbors() );
          csmp_fem_types.push_back( e.FE_Type() );
        }

      // 1.2 adding Face information after the elements
      if ( !faces_.empty() )
        for ( const auto& f : faces_ ) {
            nodes_per_element.push_back( f.Nodes() );
            neighbors_per_element.push_back( f.Neighbors() + higherDimParents + higherDimParentsFaceNum );
            csmp_fem_types.push_back( f.FE_Type() );
          }

      // 1.3 adding InterFace information after the faces
      if ( !interfaces_.empty() )
        for ( const auto& f : interfaces_ ) {
            // multiplier takes care of the multiplicated interface nodes that the InterFace will be connected to
            nodes_per_element.push_back( f.Nodes() * interfaceMultiplier );
            neighbors_per_element.push_back( f.Neighbors() + higherDimParents + higherDimParentsFaceNum + interfaceExtras );
            csmp_fem_types.push_back( f.FE_Type() );
          }

      // 1.4 resizing the VSet
      vset.Resize( csmp_fem_types, nodes_per_element, neighbors_per_element,
                   nodes_.size(), faces_.size(), interfaces_.size() );
                   
      // 1.5 setting the element types 'pelmt'
      vset.AddElementTypes( csmp_fem_types.begin(), csmp_fem_types.end() );
    }
  else { // if there is only a single element type
      const Element<dim>* const el = &(*elements_.begin());
      const FiniteElement*      fe = el->FE();
      assert( fe != nullptr );
      vset.Resize( fe->Nodes(),
                   fe->Neighbors(),
                   fe->ElementType(),
                   Nodes(), Elements() );

      // single element type mesh
      vset.ElementType( 0, el->FE_Type() );
    }

  // 2. adding node coordinates and boundary flags (BOX_BOUNDARY)
  // ------------------------------------------------------------
  size_t i( 0U );
  if constexpr ( dim == 1U ) {
      for ( const auto& n : nodes_ ) {
        vset.Px( i, n.x() );
        ++i;
      }
    }
  else if constexpr ( dim == 2U ) {
      for ( const auto& n : nodes_ ) {
        vset.Px( i, n.x() );
        vset.Py( i, n.y() );
        ++i;
      }
    }
  else {
    for ( const auto& n : nodes_ ) {
        vset.Px( i, n.x() );
        vset.Py( i, n.y() );
        vset.Pz( i, n.z() );
        ++i;
      }
   }
  // boundary flags
  vset.ResizeBFlags( /* nodes */ );
  for ( const auto& n : nodes_ )
    vset.AddBFlag( n.Idx(), n.AtBoundary() );
  
  
  // 'pelmt' was already set above
  
  // 3. adding 'plist' connectivity list and 'pmtrl'
  // -----------------------------------------------
  vector<int32_t>  pmtrl( elements_.size(), 0 );
  size_t           eidx{0};
  // elements
  for ( const auto& e : elements_ ) {
      const auto n_nodes{e.Nodes()};
      for ( auto j = 0U; j<n_nodes; ++j )
        vset.Plist( eidx, j, (e.N( j )->Idx()) );
      pmtrl[eidx] = e.Material_ID();
      ++eidx;
    }
  // adding the material identifiers
  vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );

  // faces
  if ( !faces_.empty() )
    for ( const auto& f : faces_ ) {
        const auto n_nodes{f.Nodes()};
        for ( auto j = 0U; j<n_nodes; ++j )
          vset.Plist( eidx, j, (f.N( j )->Idx()) );
        ++eidx;
      }

  // interfaces
  if ( !interfaces_.empty() )
    for ( const auto& f : interfaces_ ) {
        const auto n_nodes{f.Nodes()};
        for ( auto j = 0U; j<n_nodes; ++j )
          vset.Plist( eidx, j, (f.N( j )->Idx()) );
        ++eidx;
      }


  // 4. adding 'pfverts' neighbors per element list
  // ----------------------------------------------
   
  // 'pfverts' elements
  eidx = 0;
  for ( const auto& e : elements_ ) {
      const auto neighbors{e.Neighbors()};
      assert( neighbors <= 6 );
      for ( auto j{0}; j<neighbors; ++j ) {
            const Element<dim>* const ptr = e.Neighbor(j);
            if ( ptr ) {
                assert( ptr->Idx() < elements_.size() );
                vset.Pfvert( eidx, j, ptr->Idx() );
              }
            else vset.Pfvert( eidx, j, e.AtBoundary(j) );
         }
      ++eidx;
    }

  // 'pfverts' faces
  // ---------------
  // the faces are stored after the elements including connections to their higher-dimensional neighbors
  // add the end of the pfverts entries
  for ( const auto& f : faces_ ) {
      // equidimensional neighbors first
      const auto neighbors{ f.Neighbors() };
      for ( auto j = 0U; j<neighbors; ++j ) {
           const Face<dim>* const ptr = f.Neighbor(j);
           // if the neighbor exists (which it must on the inside of the Face)
           if ( ptr != nullptr )
             vset.Pfvert( eidx, j, ptr->Idx() );
           else vset.Pfvert( eidx, j, IRREGULAR );
        }
      // higher-dimensional neighbors second
      // inner neighbor
      assert( f.InnerParent()->IsEquidimensional() );
      assert( f.InnerParent()->Idx() < Elements() );
      vset.Pfvert( eidx, neighbors, f.InnerParent()->Idx() );
      // outer neighbor
      if ( f.OuterParent() != nullptr ) {
          assert( f.OuterParent()->IsEquidimensional() );
          assert( f.OuterParent()->Idx() < Elements() );
          vset.Pfvert( eidx, neighbors + 1U, f.OuterParent()->Idx() );
        }
      else {
          // getting the boundary placement of the inner element
          vset.Pfvert( eidx, neighbors + 1U, atBoundary( f.InnerParent(), f.InnerParentFaceID() ) );
        }
      // adding the local numbers of the faces that the Face is collocated with if any
      vset.Pfvert( eidx, neighbors + 2U, f.InnerParentFaceID() );
      // if there is no outer element, the face idx will initialised with NULL_IDX
      vset.Pfvert( eidx, neighbors + 3U, f.OuterParentFaceID() );
      ++eidx;
   }

  // 'pfverts' interfaces (which must always have two higher-dimensional neighbors)
  // ------------------------------------------------------------------------------
  for ( const auto& f : interfaces_ ) {
    // 1. equidimensional neighbors (=other interfaces) first
    //    they are written in the order in which they are stored in the interface
    //    To simplify things there is always an entry for the intervening element even if there is none.
    const auto neighbors{ f.Neighbors() };
    for ( auto j = 0U; j<neighbors; ++j ) {
         if ( f.Neighbor(j) != nullptr )
           vset.Pfvert( eidx, j, f.Idx() );
         else
           vset.Pfvert( eidx, j, INTERNAL );
      }
    // 2. inner and outer higher-dimensional neighbors (2 entries)
    //   (they must always exist because SplitBoundaries are internal model boundaries)
    assert( f.InnerParent() != nullptr );
    assert( f.OuterParent() != nullptr );
    assert( f.InnerParent()->Idx() < elements_.size() );
    assert( f.OuterParent()->Idx() < elements_.size() );
    vset.Pfvert( eidx, neighbors,      f.InnerParent()->Idx() );
    vset.Pfvert( eidx, neighbors + 1U, f.OuterParent()->Idx() );
      
    // 3. storing local number of face of the inner and outer elements that the InterFace is connected to (2 entries)
    assert( f.InnerParentFaceID() < f.InnerParent()->Faces() );
    assert( f.OuterParentFaceID() < f.OuterParent()->Faces() );
    vset.Pfvert( eidx, neighbors + 2U, f.InnerParentFaceID() );
    vset.Pfvert( eidx, neighbors + 3U, f.OuterParentFaceID() );
    
    // 4. storing number of intervening element or nullptr identifier (one entry)
    if ( f.HasInterveningElement() ) {
         assert( f.InterveningElement()->Idx() < elements_.size() );
         vset.Pfvert( eidx, neighbors + 4U, f.InterveningElement()->Idx() );
      }
    else
      vset.Pfvert( eidx, neighbors + 4U, INTERNAL );
      
    ++eidx;
  }
  
  cout << "\nMeshManager<" << dim << ">::OutputMeshTo: MeshManager successfully output to VSet..." << endl;

} // end OutputMeshTo( VSet )





/**
Writes VSet data from MeshManager to VSet, including 'pmtrl' material property information

Storing distributed variables associated with the mesh in the VSet

For all finite volumes and elements, their integration points and nodes,
but not for any of the variables stored on the model, region, boundaries or splitboundaries
this method stores the current values in the VSet.

To store the properties in the VSet, they are first written to PropertyData objects.
These are then added to the VSet property storage.

@attention  a continuous numbering of elements, faces, interfaces and nodes has to be
established, for instance with AssignUniqueNumbers() before this method is called.

@note region properties are stored together with the regions in respective binary files

@note this method is anything, but nice. Yet it will be quite a challenge to come up with a better
design; hopefully there will be no more extra variable placements or types in the future.

@attention in the VSet, the variables are identified only by their (unique) names. The property
database is therefore essential to retrieve all other variable related information.

@author SKM 5/5/2016

*/
template<uint32_t dim>
void MeshManager<dim>::OutputStoredVariablesTo( const PropertyDatabase<dim>& database, 
                                                VSet<dim>& vset ) const
{
  ErrorHandler&		csmp_error( ErrorHandler::Instance() );
  map<string, Index>  properties;
  
  // 'pmtrl' stored exclusively on elements
  // =========================================
  vector<int32_t> pmtrl;
  pmtrl.reserve( elements_.size() );
  for ( const auto& e : elements_ )
    pmtrl.push_back( e.Material_ID() );
  vset.AddPmtrl( pmtrl.begin(), pmtrl.end() );
  pmtrl.clear();
  
  
  // 'PropertyData'
  // =========================================

  database.ListProperties( NODE, properties );
  const size_t n_nodes( nodes_.size() );
  // for all node properties
  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    // setting the specifications for the property storage (no memory allocation yet)
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    // for the given property type
    const size_t flag_capacity( n_nodes * (*pit).second.flagDepth );
    const size_t data_capacity( n_nodes * (*pit).second.dataDepth );
    // allocating memory to store the property flags and values
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
      {
        case SCALAR: {
          // estimating the storage required
          ScalarVariable value;
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                pushBack( data, value );
              }
            }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                pushBack( data, value );
              }
            }
break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                pushBack( data, value );
              }
            }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                pushBack( data, value );
              }
            }
          break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                pushBack( data, value );
              }
            }
          break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of node variable not recognized." );
      }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }

  // -------------------------------------------------
  // -------------------------------------------------
  // element properties (including integration points)
  // -------------------------------------------------
  // -------------------------------------------------
  database.ListProperties( ELEMENT, properties );
  const size_t elements( elements_.size() );
  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    // creating the property storage
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    // for the given property type
    const size_t flag_capacity( elements * (*pit).second.flagDepth );
    const size_t data_capacity( elements * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( const auto& it : elements_ ) {
              it.Read( (*pit).second, value );
              pushBack( data, value );
            }
          }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( const auto& it : elements_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( const auto& it : elements_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : elements_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : elements_ ) {
              it.Read( (*pit).second, value );
              pushBack( data, value );
            }
          }
       break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of element variable not recognized." );
    }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }

  // element integration point properties
  // ------------------------------------
  if ( database.ListProperties( ELEMENT_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_ips( (*elements_.begin()).IntegrationPoints() ); // just an estimate

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      // creating the property storage
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( elements * elmt_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * elmt_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : elements_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
     break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of element integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // element sector integration point properties
  // -------------------------------------------
  if ( database.ListProperties( SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_sector_ips( (*elements_.begin()).IntegrationPointsPerSector() );
    const size_t sectors_per_element( (*elements_.begin()).Sectors() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : elements_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
                const auto sectors( it.Sectors() );
                for ( auto i = 0U; i<sectors; ++i ) {
                  const auto ips_per_sector( it.IntegrationPointsPerSector() );
                  for ( auto j = 0U; j<ips_per_sector; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    pushBack( data, value );
                  }
                }
              }
            }
         break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of element sector integraton point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // element facet integration point properties
  // ------------------------------------------
  if ( database.ListProperties( FACET_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_facet_ips( (*elements_.begin()).IntegrationPointsPerFacet() );
    const size_t facets_per_element( (*elements_.begin()).Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : elements_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // ----------------------------------------------
  // ----------------------------------------------
  // face properties (including integration points)
  // ----------------------------------------------
  // ----------------------------------------------
  database.ListProperties( FACE, properties );
  const size_t n_faces( Faces() );

  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    const size_t flag_capacity( n_faces * (*pit).second.flagDepth );
    const size_t data_capacity( n_faces  * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( const auto& it : faces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( const auto& it : faces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( const auto& it : faces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : faces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : faces_ ) {
              it.Read( (*pit).second, value );
              pushBack( data, value );
            }
          }
        break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face variable not recognized." );
    }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }


  // face integration point properties
  // ---------------------------------
  if ( database.ListProperties( FACE_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
    assert( Faces() > 0U );
    const size_t face_ips( (*faces_.begin()).IntegrationPoints() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( n_faces * face_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_faces * face_ips  * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : faces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
     break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of face integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // face sector integration point properties
  // ----------------------------------------
  if ( database.ListProperties( FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
    assert( Faces() > 0U );
    const size_t face_sector_ips( (*faces_.begin()).IntegrationPointsPerSector() );
    const size_t sectors_per_face( (*faces_.begin()).Sectors() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( n_faces * sectors_per_face * face_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_faces * sectors_per_face * face_sector_ips  * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : faces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of face-sector integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // face facet integration point properties
  // ---------------------------------------
  if ( database.ListProperties( FACE_FACET_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
    assert( Faces() > 0U );
    const size_t face_facet_ips( (*faces_.begin()).IntegrationPointsPerFacet() );
    const size_t facets_per_face( (*faces_.begin()).Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( n_faces * facets_per_face * face_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_faces * facets_per_face * face_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : faces_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto facets( it.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of face facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }


  // ---------------------------------------------------
  // ---------------------------------------------------
  // interface properties (including integration points)
  // ---------------------------------------------------
  // ---------------------------------------------------
  database.ListProperties( INTER_FACE, properties );
  const size_t n_interfaces( InterFaces() );

  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    const size_t flag_capacity( n_interfaces * (*pit).second.flagDepth );
    const size_t data_capacity( n_interfaces  * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          pushBack( data, value );
        }
      }
   break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface variable not recognized." );
    }
    // storing the data in the VSet
    vset.AddData( (*pit).first.c_str(), data );
  }

  // interface integration point properties
  // --------------------------------------
  if ( database.ListProperties( INTER_FACE_INTEGRATION_POINT, properties ) > 0 && InterFaces() > 0 ) {
    const size_t interface_ips( (*interfaces_.begin()).IntegrationPoints() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * interface_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * interface_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : interfaces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const size_t integration_points( it.IntegrationPoints() );
            for ( auto i = 0U; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              pushBack( data, value );
            }
          }
        }
     break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of interface integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }

  // interface sector integration point properties
  // ---------------------------------------------
  if ( database.ListProperties( INTER_FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 && InterFaces() > 0 ) {
    assert( InterFaces() > 0 );
    const size_t interface_sector_ips( (*interfaces_.begin()).IntegrationPointsPerSector() );
    const size_t sectors_per_interface( (*interfaces_.begin()).Sectors() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : interfaces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const auto sectors( it.Sectors() );
            for ( auto i = 0U; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( auto j = 0U; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of interface sector integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }

  // interface facet integration point properties
  // --------------------------------------------
  if ( database.ListProperties( INTER_FACE_FACET_INTEGRATION_POINT, properties ) > 0 && InterFaces() > 0 ) {
    assert( InterFaces() > 0 );
    const size_t interface_facet_ips( (*interfaces_.begin()).IntegrationPointsPerFacet() );
    const size_t facets_per_interface( (*interfaces_.begin()).Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * facets_per_interface * interface_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * facets_per_interface * interface_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );

      switch ( (*pit).second.type )
        {
          case SCALAR: {
                ScalarVariable value;
                for ( const auto& it : interfaces_ ) {
                  const auto facets( it.Facets() );
                  for ( auto i = 0U; i<facets; ++i ) {
                    const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                    for ( auto j = 0U; j<ips_per_facet; ++j ) {
                      it.Read( i, j, (*pit).second, value );
                      pushBack( data, value );
                    }
                  }
                }
              }
          break;
        case VECTOR: {
              VectorVariable<dim> value;
              for ( const auto& it : interfaces_ ) {
                const auto facets( it.Facets() );
                for ( auto i = 0U; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( auto j = 0U; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    pushBack( data, value );
                  }
                }
              }
            }
          break;
        case TENSOR: {
              TensorVariable<dim> value;
              for ( const auto& it : interfaces_ ) {
                const auto facets( it.Facets() );
                for ( auto i = 0U; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( auto j = 0U; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    pushBack( data, value );
                  }
                }
              }
            }
          break;
        case ARRAY: {
              ArrayVariable value( (*pit).second.dataDepth );
              for ( const auto& it : interfaces_ ) {
                const auto facets( it.Facets() );
                for ( auto i = 0U; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( auto j = 0U; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    pushBack( data, value );
                  }
                }
              }
            }
break;
        case FLAGGEDARRAY: {
              FlaggedArrayVariable value( (*pit).second.dataDepth );
              for ( const auto& it : interfaces_ ) {
                const auto facets( it.Facets() );
                for ( auto i = 0U; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( auto j = 0U; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    pushBack( data, value );
                  }
                }
              }
            }
break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of interface facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      vset.AddData( (*pit).first.c_str(), data );
    }
  }

} // end OutputStoredVariablesTo







  /**
  Read and assign property values written exactly by OutputStoredVariablesTo()
  
  @attention It is assumed that no erasures have occured in the MeshManager before the properties are assigned; else the order is no longer correct;

  @ test SKM 28/6/2016
  */
template<uint32_t dim>
void MeshManager<dim>::InputStoredVariablesFrom( const PropertyDatabase<dim>& database, const VSet<dim>& vset )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
  if ( vset.Vertices()   != nodes_.size() ||
       vset.Elements()   != elements_.size() ||
       vset.Faces()      != faces_.size() ||
       vset.InterFaces() != interfaces_.size() ) {
       csmp_error.notice( ERROR, "MeshManager<dim>::InputStoredVariablesFrom",
                         "mismatch between property data sizes and mesh stored in manager; no input." );
       return;
    }

  // 'pmtrl' flags for the elements
  // ==============================
  assert( elements_.size() == vset.Elements() );
  typename plf::colony<Element<dim>>::iterator it = elements_.begin();
  for ( auto pit=vset.PmtrlBegin(); pit!=vset.PmtrlEnd(); ++pit, ++it )
    (*it).Material_ID( (*pit) );

  // -------------------
  // assigns properties
  // -------------------
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
    {
      // apart from the name string key in the map, PropertyData contains the most important variable specifications
      if ( (*pit).second.Placement() != NODE ) continue;
      // some checks
      assert( database.IsDefined( (*pit).first.c_str() ) );
      const Index key( database.StorageKey( (*pit).first.c_str() ) );
      assert( key.place == NODE );
      assert( (*pit).second.Size() == nodes_.size() * key.dataDepth );

      // for the given property type
      switch ( key.type )
      {
        case SCALAR: {
          ScalarVariable value;
          size_t i( 0U );
          for ( auto& n : nodes_ ) {
                read( (*pit).second, i++, value );
                n.Store( key, value );
              }
            }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          size_t i( 0U );
          for ( auto& n : nodes_ ) {
                read( (*pit).second, i++, value );
                n.Store( key, value );
              }
            }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          size_t i( 0U );
          for ( auto& n : nodes_ ) {
                read( (*pit).second, i++, value );
                n.Store( key, value );
              }
            }
          break;
        case ARRAY: {
          ArrayVariable value( key.dataDepth );
          size_t i( 0U );
          for ( auto& n : nodes_ ) {
                read( (*pit).second, i++, value );
                n.Store( key, value );
              }
            }
          break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( key.dataDepth );
          size_t i( 0U );
          for ( auto& n : nodes_ ) {
                read( (*pit).second, i++, value );
                n.Store( key, value );
              }
            }
          break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of node variable not recognized." );
      }
    }

  // -------------------------------------------------
  // -------------------------------------------------
  // element properties (including integration points)
  // -------------------------------------------------
  // -------------------------------------------------
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
    {
      if ( (*pit).second.Placement() != ELEMENT ) continue;
      // some checks
      assert( database.IsDefined( (*pit).first.c_str() ) );
      const Index key( database.StorageKey( (*pit).first.c_str() ) );
      assert( key.place == ELEMENT );
      assert( (*pit).second.Size() / key.dataDepth == Elements() );

      switch ( key.type )
        {
          case SCALAR: {
            ScalarVariable value;
            for ( auto& e : elements_ ) {
                  read( (*pit).second, e.Idx(), value );
                  e.Store( key, value );
                }
              }
 break;
          case VECTOR: {
            VectorVariable<dim> value;
            size_t i( 0U );
                for ( auto& e : elements_ ) {
                  read( (*pit).second, i++, value );
                  e.Store( key, value );
                }
              }
 break;
          case TENSOR: {
            TensorVariable<dim> value;
            size_t i( 0U );
            for ( auto& e : elements_ ) {
                  read( (*pit).second, i++, value );
                  e.Store( key, value );
                }
              }
 break;
          case ARRAY: {
            ArrayVariable value( key.dataDepth );
            size_t i( 0U );
            for ( auto& e : elements_ ) {
                  read( (*pit).second, i++, value );
                  e.Store( key, value );
                }
              }
 break;
          case FLAGGEDARRAY: {
            FlaggedArrayVariable value( key.dataDepth );
            size_t i( 0U );
            for ( auto& e : elements_ ) {
                  read( (*pit).second, i++, value );
                  e.Store( key, value );
                }
              }
 break;
          default:
            csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                               (*pit).first, "type of element variable not recognized." );
        }
    }

  // element integration point properties
  // ------------------------------------
  // ELEMENT_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
    {
      if ( (*pit).second.Placement() != ELEMENT_INTEGRATION_POINT ) continue;
      // some checks
      assert( database.IsDefined( (*pit).first.c_str() ) );
      const Index key( database.StorageKey( (*pit).first.c_str() ) );
      assert( key.place == ELEMENT_INTEGRATION_POINT );

      switch ( key.type )
        {
          case SCALAR: {
            ScalarVariable value;
            size_t entry( 0U ); // running index
            for ( auto& e : elements_ ) {
              const size_t integration_points( e.IntegrationPoints() );
              for ( auto i = 0U; i<integration_points; ++i ) {
                    read( (*pit).second, entry, value );
                    e.Store( i, key, value );
                    entry++;
                  }
                }
              }
 break;
          case VECTOR: {
            VectorVariable<dim> value;
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              const size_t integration_points( e.IntegrationPoints() );
              for ( auto i = 0U; i<integration_points; ++i ) {
                    read( (*pit).second, entry, value );
                    e.Store( i, key, value );
                    entry++;
                  }
                }
              }
 break;
          case TENSOR: {
            TensorVariable<dim> value;
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              const size_t integration_points( e.IntegrationPoints() );
              for ( auto i = 0U; i<integration_points; ++i ) {
                    read( (*pit).second, entry, value );
                    e.Store( i, key, value );
                    entry++;
                  }
                }
              }
 break;
          case ARRAY: {
            ArrayVariable value( key.dataDepth );
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              const size_t integration_points( e.IntegrationPoints() );
              for ( auto i = 0U; i<integration_points; ++i ) {
                    read( (*pit).second, entry, value );
                    e.Store( i, key, value );
                    entry++;
                  }
                }
              }
 break;
          case FLAGGEDARRAY: {
            FlaggedArrayVariable value( key.dataDepth );
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              const size_t integration_points( e.IntegrationPoints() );
              for ( auto i = 0U; i<integration_points; ++i ) {
                    read( (*pit).second, entry, value );
                    e.Store( i, key, value );
                    entry++;
                  }
                }
              }
 break;
          default:
            csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                               (*pit).first, "type of element integration point variable not recognized." );
        }
    }

  // element sector integration point properties
  // -------------------------------------------
  // SECTOR_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
    {
      if ( (*pit).second.Placement() != SECTOR_INTEGRATION_POINT ) continue;
      // some checks
      assert( database.IsDefined( (*pit).first.c_str() ) );
      const Index key( database.StorageKey( (*pit).first.c_str() ) );
      assert( key.place == SECTOR_INTEGRATION_POINT );

      switch ( key.type )
        {
          case SCALAR: {
            ScalarVariable value;
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto sectors( e.Sectors() );
              for ( auto i = 0U; i<sectors; ++i ) {
                const auto ips_per_sector( e.IntegrationPointsPerSector() );
                for ( auto j = 0U; j<ips_per_sector; ++j ) {
                      read( (*pit).second, entry, value );
                      e.Store( i, j, key, value );
                      entry++;
                    }
                  }
                }
              }
 break;
          case VECTOR: {
            VectorVariable<dim> value;
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto sectors( e.Sectors() );
              for ( auto i = 0U; i<sectors; ++i ) {
                const auto ips_per_sector( e.IntegrationPointsPerSector() );
                for ( auto j = 0U; j<ips_per_sector; ++j ) {
                      read( (*pit).second, entry, value );
                      e.Store( i, j, key, value );
                      entry++;
                    }
                  }
                }
              }
 break;
          case TENSOR: {
            TensorVariable<dim> value;
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto sectors( e.Sectors() );
              for ( auto i = 0U; i<sectors; ++i ) {
                const auto ips_per_sector( e.IntegrationPointsPerSector() );
                for ( auto j = 0U; j<ips_per_sector; ++j ) {
                      read( (*pit).second, entry, value );
                      e.Store( i, j, key, value );
                      entry++;
                    }
                  }
                }
              }
 break;
          case ARRAY: {
            ArrayVariable value( key.dataDepth );
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto sectors( e.Sectors() );
              for ( auto i = 0U; i<sectors; ++i ) {
                const auto ips_per_sector( e.IntegrationPointsPerSector() );
                for ( auto j = 0U; j<ips_per_sector; ++j ) {
                      read( (*pit).second, entry, value );
                      e.Store( i, j, key, value );
                      entry++;
                    }
                  }
                }
              }
 break;
          case FLAGGEDARRAY: {
            FlaggedArrayVariable value( key.dataDepth );
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto sectors( e.Sectors() );
              for ( auto i = 0U; i<sectors; ++i ) {
                const auto ips_per_sector( e.IntegrationPointsPerSector() );
                for ( auto j = 0U; j<ips_per_sector; ++j ) {
                      read( (*pit).second, entry, value );
                      e.Store( i, j, key, value );
                      entry++;
                    }
                  }
                }
              }
 break;
          default:
            csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                               (*pit).first, "type of element sector integraton point variable not recognized." );
        }
    }

  // element facet integration point properties
  // ------------------------------------------
  // FACET_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACET_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACET_INTEGRATION_POINT );

    switch ( key.type )
      {
        case SCALAR: {
          ScalarVariable value;
          size_t entry( 0U );
          for ( auto& e : elements_ ) {
            if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
            const auto facets( e.Facets() );
            for ( auto i = 0U; i<facets; ++i ) {
              const auto ips_per_facet( e.IntegrationPointsPerFacet() );
              for ( auto j = 0U; j<ips_per_facet; ++j ) {
                read( (*pit).second, entry, value );
                e.Store( i, j, key, value );
                entry++;
              }
            }
          }
        }
        break;
      case VECTOR: {
            VectorVariable<dim> value;
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto facets( e.Facets() );
              for ( auto i = 0U; i<facets; ++i ) {
                const auto ips_per_facet( e.IntegrationPointsPerFacet() );
                for ( auto j = 0U; j<ips_per_facet; ++j ) {
                  read( (*pit).second, entry, value );
                  e.Store( i, j, key, value );
                  entry++;
                }
              }
            }
          }
        break;
      case TENSOR: {
            TensorVariable<dim> value;
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto facets( e.Facets() );
              for ( auto i = 0U; i<facets; ++i ) {
                const auto ips_per_facet( e.IntegrationPointsPerFacet() );
                for ( auto j = 0U; j<ips_per_facet; ++j ) {
                  read( (*pit).second, entry, value );
                  e.Store( i, j, key, value );
                  entry++;
                }
              }
            }
          }
        break;
      case ARRAY: {
            ArrayVariable value( key.dataDepth );
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto facets( e.Facets() );
              for ( auto i = 0U; i<facets; ++i ) {
                const auto ips_per_facet( e.IntegrationPointsPerFacet() );
                for ( auto j = 0U; j<ips_per_facet; ++j ) {
                  read( (*pit).second, entry, value );
                  e.Store( i, j, key, value );
                  entry++;
                }
              }
            }
          }
       break;
      case FLAGGEDARRAY: {
            FlaggedArrayVariable value( key.dataDepth );
            size_t entry( 0U );
            for ( auto& e : elements_ ) {
              if ( e.FE_Type() == ISOPARAMETRIC_LINEAR_BAR ) continue;
              const auto facets( e.Facets() );
              for ( auto i = 0U; i<facets; ++i ) {
                const auto ips_per_facet( e.IntegrationPointsPerFacet() );
                for ( auto j = 0U; j<ips_per_facet; ++j ) {
                  read( (*pit).second, entry, value );
                  e.Store( i, j, key, value );
                  entry++;
                }
              }
            }
          }
       break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of facet integration point variable not recognized." );
    }
  }



  // ----------------------------------------------
  // ----------------------------------------------
  // face properties (including integration points)
  // ----------------------------------------------
  // ----------------------------------------------
  // FACE
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE );
    assert( (*pit).second.Size() / key.dataDepth == Faces() );

    // for the given property type
    switch ( key.type )
      {
        case SCALAR: {
              ScalarVariable value;
              size_t i( 0U );
              for ( auto& f : faces_ ) {
                read( (*pit).second, i++, value );
                f.Store( key, value );
              }
            }
          break;
        case VECTOR: {
              VectorVariable<dim> value;
              size_t i( 0U );
              for ( auto& f : faces_ ) {
                read( (*pit).second, i++, value );
                f.Store( key, value );
              }
            }
          break;
        case TENSOR: {
              TensorVariable<dim> value;
              size_t i( 0U );
              for ( auto& f : faces_ ) {
                read( (*pit).second, i++, value );
                f.Store( key, value );
              }
            }
          break;
        case ARRAY: {
            ArrayVariable value( key.dataDepth );
            size_t i( 0U );
            for ( auto& f : faces_ ) {
              read( (*pit).second, i++, value );
              f.Store( key, value );
            }
          }
         break;
        case FLAGGEDARRAY: {
              FlaggedArrayVariable value( key.dataDepth );
              size_t i( 0U );
              for ( auto& f : faces_ ) {
                read( (*pit).second, i++, value );
                f.Store( key, value );
              }
            }
         break;
        default:
          csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                             (*pit).first, "type of face variable not recognized." );
      }
    }

  // face integration point properties
  // ---------------------------------
  // FACE_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face integration point variable not recognized." );
    }
  }

  // face sector integration point properties
  // ----------------------------------------
  // FACE_SECTOR_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE_SECTOR_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE_SECTOR_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face-sector integration point variable not recognized." );
    }
  }

  // face facet integration point properties
  // ---------------------------------------
  // FACE_FACET_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != FACE_FACET_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == FACE_FACET_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : faces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
                  read( (*pit).second, entry, value );
                  e.Store( i, j, key, value );
                  entry++;
                }
              }
            }
          }
        break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of face facet integration point variable not recognized." );
    }
  }



  // ---------------------------------------------------
  // ---------------------------------------------------
  // interface properties (including integration points)
  // ---------------------------------------------------
  // ---------------------------------------------------
  // INTER_FACE
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE );
    assert( (*pit).second.Size() / key.dataDepth == InterFaces() );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t i( 0U );
        for ( auto& e : interfaces_ ) {
          read( (*pit).second, i, value );
          e.Store( key, value );
          i++;
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t i( 0U );
        for ( auto& e : interfaces_ ) {
          read( (*pit).second, i, value );
          e.Store( key, value );
          i++;
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t i( 0U );
        for ( auto& e : interfaces_ ) {
          read( (*pit).second, i, value );
          e.Store( key, value );
          i++;
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( key.dataDepth );
        size_t i( 0U );
        for ( auto& e : interfaces_ ) {
          read( (*pit).second, i, value );
          e.Store( key, value );
          i++;
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
        size_t i( 0U );
        for ( auto& e : interfaces_ ) {
          read( (*pit).second, i, value );
          e.Store( key, value );
          i++;
        }
      }
   break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface variable not recognized." );
    }
  }

  // interface integration point properties
  // --------------------------------------
  // INTER_FACE_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const size_t integration_points( e.IntegrationPoints() );
          for ( auto i = 0U; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
   break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface integration point variable not recognized." );
    }
  }

  // interface sector integration point properties
  // ---------------------------------------------
  // INTER_FACE_SECTOR_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE_SECTOR_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE_SECTOR_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto sectors( e.Sectors() );
          for ( auto i = 0U; i<sectors; ++i ) {
            const auto ips_per_sector( e.IntegrationPointsPerSector() );
            for ( auto j = 0U; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
   break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface sector integration point variable not recognized." );
    }
  }

  // interface facet integration point properties
  // --------------------------------------------
  // INTER_FACE_FACET_INTEGRATION_POINT
  for ( auto pit = vset.PropertyValuesBegin(); pit != vset.PropertyValuesEnd(); ++pit )
  {
    if ( (*pit).second.Placement() != INTER_FACE_FACET_INTEGRATION_POINT ) continue;
    // some checks
    assert( database.IsDefined( (*pit).first.c_str() ) );
    const Index key( database.StorageKey( (*pit).first.c_str() ) );
    assert( key.place == INTER_FACE_FACET_INTEGRATION_POINT );

    switch ( key.type )
    {
      case SCALAR: {
        ScalarVariable value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( key.dataDepth );
        size_t entry( 0U );
        for ( auto& e : interfaces_ ) {
          const auto facets( e.Facets() );
          for ( auto i = 0U; i<facets; ++i ) {
            const auto ips_per_facet( e.IntegrationPointsPerFacet() );
            for ( auto j = 0U; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
   break;
      default:
        csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface facet integration point variable not recognized." );
    }
  }

} // end InputStoredVariablesFrom




// ============================================================================================
//
//           DIAGNOSTICS
//
// ============================================================================================


/**
    Performs a node-to-node, breadth-first traversal to identify whether the model consists  of disconnected mesh patches.
 */
template<uint32_t dim>
bool  MeshManager<dim>::IsContiguous() const
 {
    if ( elements_.empty() )
      throw csmp::Exception( ERROR, "MeshManager<dim>::IsContiguous", "'elements_' container is empty." );
 
   set<Node<dim>*> node_pointers;
   findInterconnectedNodeCluster<dim>( const_cast<Node<dim>*>(&(*nodes_.begin())), node_pointers );
   
   if ( node_pointers.size() < nodes_.size() ) return false;
   return true;

 } // end IsContiguous



/* ORIGINAL METHOD
     
    // finding the highest dimensional elements in the mesh
    const Element<dim>* eptr(nullptr);
    
    // looking for volume elements
    size_t elmt_count{0};
    if constexpr ( dim == 3 ) {
         for ( auto& it : elements_ )
           if ( it.IsVolumeElement() ) {
                eptr = &it;
                elmt_count++;
             }
      }
    else if constexpr ( dim == 2 ) {
         for ( auto& it : elements_ )
           if ( it.IsSurfaceElement() ) {
                eptr = &it;
                elmt_count++;
             }
      }
    else eptr = &(*elements_.begin());

    set<Element<dim>* const> contiguous_subset;
    floodFill( const_cast<Element<dim>* const>(eptr), contiguous_subset );

    // performing a floodfill on them
*/


template<uint32_t dim>
bool  MeshManager<dim>::HybridElementMesh() const
{
  return hybrid_element_mesh_;
}




/**
        @note SKM retained this method for the moment (6/9/2021)
        
       @author JCK
       @date 2018
*/
template<uint32_t dim>
size_t  MeshManager<dim>::CheckElementConnectivity() const
{
  cout <<"\nMeshManager::CheckElementConnectivity: checking mesh..."<< endl;
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  if ( !Elements() )
    throw csmp::Exception( ERROR, "MeshManager::CheckElementConnectivity", "the mesh contains no elements." );

  // --------------------------------------------------------------------
  // 3. detecting whether region contains lower-dimensional elements
  // --------------------------------------------------------------------
  size_t errors(0);

  bool with_volume_elements( false );
  bool with_surface_elements( false );
  bool with_line_elements( false );

  for ( const auto& it : elements_ ) {
      if      ( !with_line_elements    && it.IsLineElement() )		with_line_elements = true;
      else if ( !with_surface_elements && it.IsSurfaceElement() )	with_surface_elements = true;
      else if ( !with_volume_elements  && it.IsVolumeElement() )	with_volume_elements = true;
    }

  int32_t dimension_counter( 0 );
  int32_t highest_spatial_dim( 1 );
  if ( with_volume_elements )  dimension_counter++;
  if ( with_surface_elements ) dimension_counter++;
  if ( with_line_elements )    dimension_counter++;
  if ( with_volume_elements )  highest_spatial_dim = 3;
  else if ( with_surface_elements ) highest_spatial_dim = 2;

  // -----------------------------------------------------------------
  // 4. distinguishing boundary from interior elements, same for nodes
  //   (at this point the elements and nodes are already known)
  // -----------------------------------------------------------------
  size_t interior_elmts( 0 );
  size_t boundary_elmts( 0 );
  set<const Node<dim>*> boundary_nodes;
  vector<uint32_t>        fnids;

  // 4.1 If all elements have the same spatial dimension
  // ---------------------------------------------------
  if ( dimension_counter == 1 ) {
    for ( const auto& eit : elements_ ) {
      // identifying the boundary faces and their nodes
      // (each face potentially has a neighbor element)
      long  nbors_that_belong_to_group( eit.Neighbors() );
      for ( auto i = 0U; i<eit.Faces(); i++ )
        // if the face is at a model boundary
        if ( eit.Neighbor( i ) == nullptr )
          {
            // boundary nodes
            assert( eit.FE() != NULL );
            eit.FE()->NodesOfFace( i, fnids );
            for ( size_t j = 0U; j<fnids.size(); ++j )
              boundary_nodes.insert( eit.N( fnids[j] ) );
            // counting neighbors
            nbors_that_belong_to_group--;
          }

      // storing the distinguished elements in the respective vectors
      // ------------------------------------------------------------
      // interior elements
      if ( nbors_that_belong_to_group == eit.Neighbors() )
        ++interior_elmts;
      // elements with at least one face on the region boundary
      else
        ++boundary_elmts;
    }
  }
  // 4.2 If there are elements with different spatial dimensions
  // -----------------------------------------------------------
  //     the ones with highest dimensions are used to define perimeter
  //     all lower dimensional mesh that sticks out is flagged as perimeter as well.
  else {
    // a. identify the boundary elements among the highest dimensional elements,
    //    also collecting all their node pointers into a set.
    set<const Element<dim>*> lesser_dim_elmts;
    set<const Node<dim>*>    highest_dim_elmt_nodes;

    for ( const auto& eit : elements_ ) {
      // elements of the highest spatial dimension are used to define the boundary
      assert( parseFiniteElementDimension( eit.FE_Type() ) != 0 );
      if ( parseFiniteElementDimension( eit.FE_Type() ) == highest_spatial_dim )
      {
        // creating a subset with their nodes
        for ( auto i = 0U; i<eit.Nodes(); ++i ) {
          assert( eit.N( i ) != nullptr );
          highest_dim_elmt_nodes.insert( eit.N( i ) );
        }
        // if the element has faces that lie on the region boundary
        // it is considered a boudary element
        long  nbors_that_belong_to_group( eit.Neighbors() );
        for ( auto i = 0U; i<eit.Faces(); ++i )
          // 1) the element is on model boundary  or  2) one of its neighbors does not belong to its parent region
          if ( eit.Neighbor( i ) == nullptr )
            {
              nbors_that_belong_to_group--;
            }
        if ( nbors_that_belong_to_group == eit.Neighbors() ) ++interior_elmts;
        else ++boundary_elmts;
      }
      else lesser_dim_elmts.insert( const_cast<Element<dim>*>(&eit) );
    }
    assert( /* all elements are accounted for */ Elements() == interior_elmts + boundary_elmts + lesser_dim_elmts.size() );
    if ( Elements() != interior_elmts + boundary_elmts + lesser_dim_elmts.size() )
      errors++;

    set<const Element<dim>*> lesser_dim_elmts_detached; // to distinguish stand-alone lower dimensional mesh

    for ( const auto& e : lesser_dim_elmts )
      {
        // a) lower-dim elements that may be sticking out
        // ----------------------------------------------
        // lower-dimensional elements with nodes that do not belong to the node set of the
        // higher dimensional elements must be boundary elements
        size_t  exterior_nodes( 0U );
        for ( auto i = 0U; i<e->Nodes(); ++i )
          if ( !highest_dim_elmt_nodes.count( e->N( i ) ) )
            exterior_nodes++;

        // if individual nodes stick out the parent element sticks out as well.
        if ( exterior_nodes == e->Nodes() )
          lesser_dim_elmts_detached.insert( e );
      }

    if ( !lesser_dim_elmts_detached.empty() ) {
      csmp_error.notice( WARNING, "MeshManager::CheckElementConnectivity:",
                        "model contains lower-dimensional elements detached from higher dimensional elements" );

      // do some additional diagnostics on these elements
      // ------------------------------------------------
      cerr << "\n\tdetached elements: " << lesser_dim_elmts_detached.size() << ":";
      for ( const auto& e : lesser_dim_elmts_detached )
        cerr << " " << e->Idx();
      cerr << endl;
      errors++;
    }
  } // end multi-dim element region

  cout <<"\nMeshManager::CheckElementConnectivity: completed model check."<< endl;
  if ( errors > 0 ) cerr <<"\t"<< errors <<" major errors encountered."<< endl;
  
  return errors;
    
} // end CheckElementConnectivity






template<uint32_t dim>
void MeshManager<dim>::Out() const
{
  cout << "\nMeshManager<" << dim << ">::Out: " << endl;
  
  // nodes
  cout << "\nNODES: " << endl;
  size_t n_node{0};
  for ( const auto& n : nodes_ ) {
       string bound = parseBoundary( n.AtBoundary() );
       cout << "\nNode ID: " << n.Idx() << " ";
       cout << n.Coordinate();
       cout << " Boundary flag: " << bound << endl;
      n_node++;
    }

  // elements
  size_t n_elmt{0};
  cout << "\nELEMENTS: " << endl;
  for ( const auto& e : elements_ ) {
         string bound("NOT");
         for ( auto i{0}; i<e.Neighbors(); ++i )
           if ( e.Neighbor(i) == nullptr ) {
                bound = parseBoundary( atBoundary(&e,i) );
                break;
             }
         cout << "\nElement ID: " << e.Idx() <<" ("<< parseFiniteElementType(e.FE_Type());
         cout <<"), Boundary flag: " << bound << endl;
         cout << "Member Nodes: " << endl;
         for ( auto i = 0U; i < e.Nodes(); i++ )
           cout << e.N( i )->Idx() << "\t";
         cout << "\nNeighbor elements: " << endl;
         for ( auto i = 0U; i < e.Neighbors(); i++ )
           if ( e.Neighbor( i ) != nullptr )
             cout << e.Neighbor( i )->Idx() << "\t";
           else
             cout << "NO NEIGHBOR\t";
         cout << endl;
       n_elmt++;
    }

  // faces
  size_t n_face{0};
  cout << "\nFACES: " << endl;
  for ( const auto& f : faces_ ) {
        cout << "\nFace ID: " << f.Idx() <<" ("<< parseFiniteElementType(f.FE_Type()) <<")."<< endl;
        cout << "Member Nodes: " << endl;
        for ( auto i = 0U; i < f.Nodes(); i++ )
          cout << f.N( i )->Idx() << "\t";
        cout << "\nNeighbor faces: " << endl;
        for ( auto i = 0U; i < f.Neighbors(); i++ )
          if ( f.Neighbor( i ) != nullptr )
            cout << f.Neighbor( i )->Idx() << "\t";
          else
            cout << "NO NEIGHBOR\t";

        cout << endl;
      n_face++;
    }

  // inter faces
  size_t n_iface{0};
  cout << "\nINTERFACES: " << endl;
  for ( const auto& f : interfaces_ ) {
        cout << "\nInterFace ID: " << f.Idx() <<" ("<< parseFiniteElementType(f.FE_Type()) <<")."<< endl;
        cout << "Member Nodes: " << endl;
        for ( auto i = 0U; i < f.Nodes(); i++ )
          cout << f.N( i )->Idx() << "\t";
        cout << "\nNeighbor faces: " << endl;
        for ( auto i = 0U; i < f.Neighbors(); i++ )
          if ( f.Neighbor( i ) != NULL )
            cout << f.Neighbor( i )->Idx() << "\t";
          else
            cout << "NO NEIGHBOR\t";
        cout << endl;
      n_iface++;
    }

  // parent elements ID's for each node
  n_node = 0;
  cout << endl << endl;
  cout << "PARENT ELEMENT INFORMATION FOR ALL NODES: " << endl;
  for ( const auto& n : nodes_ ) {
      cout << "\nNode: " << n.Idx() << ", parent elements: " << endl;
      for ( auto i = 0u; i < n.Parents(); i++ )
        cout << n.Parent( i )->Idx() << " ";
      cout << endl;
      n_node++;
    }
 
  // finite element manager
  fem_manager_.Out();
  
  // finite volume stencils
  list<CSMP_FEM_TYPE> etypes;
  fem_manager_.CurrentElementTypes( etypes );
  cout <<"\nFinite volume stencils: ";
  for ( const auto& fit : etypes ) {
       cout <<"\n"<< parseFiniteElementType( fit );
       fvm_manager_.Stencil(fit)->Out();
    }
  cout << endl;
  
  // node manifolds
  if ( node_manifold_manager_ != nullptr ) {
      cout << endl << endl;
      cout << "NODE MANIFOLDS: " << endl;
      node_manifold_manager_->Out();
    }
    
} // end Out


// POTENTIALLY NEEDED METHODS


/*
  /// replace lower-dimensional element with Face object, deleting the Elements and establishing the neighbor connectivity of the new Faces
  size_t          ReplaceElementsByFaces( const PropertyDatabase<dim>&,
                                          typename std::vector<Element<dim>*>::iterator first,
                                          typename std::vector<Element<dim>*>::iterator last );
*/




template class MeshManager<1U>;
template class MeshManager<2U>;
template class MeshManager<3U>;


} // end namespace csmp 
