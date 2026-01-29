#include "MeshManager.h"
#include "MeshPatch.h"
#include "meshManagementUtilities.h"
#include "CSMP_highLevelUtilities.h"
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
#include "FaceConstructionData.h"
#include "UnionFind.h"


#define CSMP_MESH_MANAGER_DEBUG

using namespace std;

namespace csmp {

template<uint32_t dim>
MeshManager<dim>::MeshManager()
  : fem_manager_( dim, 1, true ), // linear interpolation functions, isoparametric elements
    fvm_manager_(nullptr),
    hybrid_element_mesh_( false )
{
}




/**
    Constructs a customised version of the MeshManager , appropriately initialising the FiniteElementManager.
    If the model uses linear, isoparametric elements, corresponding finite volume stencils are constructed and initialised.
*/
template<uint32_t dim>
MeshManager<dim>::MeshManager( const PropertyDatabase<dim>& pref, const VSet<dim>& vset )
  : fem_manager_( dim, vset.OrderOfFiniteElementInterpolationFunctions(), vset.IsoparametricElementMesh() ),
    fvm_manager_( (vset.OrderOfFiniteElementInterpolationFunctions()==1u && vset.IsoparametricElementMesh() ) ? new FiniteVolumeStencilManager<dim>(fem_manager_) : nullptr ),
    hybrid_element_mesh_( vset.HybridElementTypeMesh() )
{
   assert( pref.VariableCount() > 0 );
   assert( vset.Vertices() > 0 );
   const bool with_FV_variables = (finiteVolumeVariables(pref) > 0 ) ? true : false;
   Initialize( pref, vset, with_FV_variables );
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
     delete fvm_manager_;
    
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
size_t MeshManager<dim>::Interfaces() const
 { return interfaces_.size(); }

template<uint32_t dim>
size_t MeshManager<dim>::NodeManifolds() const
 {
    if ( node_manifold_manager_ == nullptr ) return 0U;
    return node_manifold_manager_->Manifolds();
 }

template<uint32_t dim>
bool MeshManager<dim>::HasNodeManifolds() const
 {
    if ( node_manifold_manager_ == nullptr ) return false;
    return true;
 }
 

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
  typename plf::colony<InterFace<dim>>::iterator MeshManager<dim>::InterfacesBegin()
  { return interfaces_.begin(); }
  
template<uint32_t dim>
  typename plf::colony<InterFace<dim>>::iterator MeshManager<dim>::InterfacesEnd()
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
  typename plf::colony<InterFace<dim>>::const_iterator MeshManager<dim>::InterfacesBegin() const
  { return interfaces_.begin(); }

template<uint32_t dim>
  typename plf::colony<InterFace<dim>>::const_iterator MeshManager<dim>::InterfacesEnd() const
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
bool  MeshManager<dim>::Initialize( const PropertyDatabase<dim>& phys_vars, const VSet<dim>& vset, bool initialise_FV_stencils )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );

  hybrid_element_mesh_ = vset.HybridElementTypeMesh();
  cout <<"\nMeshManager<"<< dim <<">::Initialize: building mesh with "<< vset.Elements() <<" elements, "<< vset.Vertices() <<" nodes, ";
  cout << vset.Faces() <<" faces, and "<< vset.Interfaces() <<" interfaces.\n";
  if ( vset.HybridElementTypeMesh() ) cout <<"mesh consists of multiple element types.\n";
  if ( vset.Faces() > 0 ) cout <<"mesh contains "<< vset.Faces() <<"'Boundary' objects.\n";
  if ( vset.Interfaces() > 0 ) cout <<"mesh contains "<< vset.Interfaces() <<"'SplitBoundary' objects.\n";
  cout << endl;

  // ------------------------------------------------------------------------------------
  // 1. check availability of necessary finite element & finite volume types
  // ------------------------------------------------------------------------------------
  // initializing the finite-element manager true=isoparametric
  fem_manager_.InitializeElements( dim, vset.OrderOfFiniteElementInterpolationFunctions(), vset.IsoparametricElementMesh() );

  if ( csmp_error.Verbose() )
    cout << "\nMeshManager<" << dim << ">::Initialize: checking the availability of the necessary finite element types..." << endl;

  set<CSMP_FEM_TYPE> input_etypes;
  input_etypes.insert( parseFiniteElementTypeEnum( vset.ElementType( 0U ) ) );

  if ( vset.HybridElementTypeMesh() )
    for ( size_t i{0U}; i<vset.Cells(); ++i )
      input_etypes.insert( parseFiniteElementTypeEnum( vset.ElementType( i ) ) );

  cout << "\nMeshManager<" << dim << ">::Initialize: ";
  cout << "input VSet contains the following finite element types:\n\t";
  for ( const auto& iit : input_etypes ) {
      cout << parseFiniteElementType( iit ) << "  ";
      if ( !fem_manager_.ContainsElementType( iit ) ) {
        cerr << "\n\n\tFinite element type not available: " << parseFiniteElementType( iit ) << endl;
        fem_manager_.Out();
        throw Exception( FATAL_ERROR,
                         "MeshManager<dim>::Initialize(VSet):",
                         "'FiniteElementManager' lacks finite-element type required by VSet." );
      }
    }
  cout << endl;
  
  // finite volume stencils
  if ( initialise_FV_stencils ) {
       const bool assign_stencils_to_elements{false}; // the elements have not been created yet
       InitializeFiniteVolumeStencils( phys_vars, assign_stencils_to_elements );
    }
    

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
            for ( uint32_t j{0U}; j<dim; ++j ) coord[j] = vset.P( j, idx );
            nodes_.emplace( Node<dim>( idx, Point<dim>( coord ), nvars,
                                       static_cast<BOX_BOUNDARY>(vset.BFlag(idx)),
                                       static_cast<TOPOTYPE>(vset.BREP_Flag(idx)) ) );
        }
    }

  // storage for elements
   {
      const LocalVariables evars( phys_vars.LocalVariablesAt( ELEMENT ) );
      const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( ELEMENT ) );
      auto first( vset.PlistElmtsBegin() ), last( vset.PlistElmtsEnd() );

      size_t elmt_idx{0ul};
      // 2.1 If the MeshManager contains only one element type
      if ( !vset.HybridElementTypeMesh() ) {
          const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( 0U ));
          const FiniteVolumeStencil<dim>* const stencil_ptr = (initialise_FV_stencils==true) ?
                                                               fvm_manager_->Stencil( csmpElementType ) :
                                                               static_cast<const FiniteVolumeStencil<dim>* const>(nullptr);
          while ( first != last )
            {
              // create the element
              typename plf::colony<Element<dim>>::iterator
                eit = elements_.emplace( Element<dim>( elmt_idx, fem_manager_.E( csmpElementType ),
                                                                                 stencil_ptr,
                                                                                 evars, cvars, vset.Pmtrl(elmt_idx) ) );
              // assign the nodes
              const auto nodes( fem_manager_.E( csmpElementType )->Nodes() );
              for ( uint32_t j{0U}; j < nodes; ++j )
                (*eit).Assign( j, &(*next(nodes_.begin(), static_cast<long>(vset.Plist( elmt_idx, j )))) );
                
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
              const FiniteVolumeStencil<dim>* const stencil_ptr = (initialise_FV_stencils==true) ?
                                                                   fvm_manager_->Stencil( csmpElementType ) :
                                                                   static_cast<const FiniteVolumeStencil<dim>* const>(nullptr);
              typename plf::colony<Element<dim>>::iterator
                eit = elements_.emplace( Element<dim>( elmt_idx, fem_manager_.E( csmpElementType ),
                                                       stencil_ptr,
                                                       evars, cvars, vset.Pmtrl(elmt_idx) ) );
                                                       
              const auto nodes( fem_manager_.E( csmpElementType )->Nodes() );
              for ( auto j{0U}; j < nodes; j++ ) (*eit).Assign( j, &(*next(nodes_.begin(), static_cast<long>(vset.Plist( elmt_idx, j )))) );
              elmt_idx++;
              first++;
            }
        }
      assert( elements_.size() == vset.Elements() );
   }
 
  // 2.3 Assign neighbor elements to elements
  const size_t  n_elmts = elements_.size();
  if ( vset.WithNeighbourConnectivity() ) {
       if ( csmp_error.Verbose() )
          cout << "\nMeshManager<" << dim << ">::Initialize: assigning neighbors to elements..." << endl;
       if ( !hybrid_element_mesh_ ) {
            const uint32_t n_neighbors = fem_manager_.E( vset.ElementType( 0U ) )->Neighbors();
            for ( auto& e : elements_ )
              for ( uint32_t j{0U}; j < n_neighbors; ++j ) {
                      const int64_t  index{ vset.Pfvert( e.Idx(), j ) };
                      if ( index >= static_cast<int64_t>(n_elmts) ) {
                           cerr <<"\n\t"<< index <<" vs. number of elements = "<< n_elmts << endl;
                           csmp_error.Note( ERROR, "MeshManager::Initialise: ", "element ID in 'pfverts' out of range.");
                        }
                      // only if there is a neighbor something needs to be done; else nullptr was aready assigned
                      else if ( index >= 0 )
                        e.Assign( j, &(*next(elements_.begin(),index)) );
                  }
         }
       // if this is a hybrid element mesh
       else {
           for ( auto& e : elements_ ) {
                const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( e.Idx() ));
                const uint32_t n_neighbors( fem_manager_.E( csmpElementType )->Neighbors() );
                for ( uint32_t j{0U}; j < n_neighbors; ++j ) {
                      const int64_t  index{ vset.Pfvert( e.Idx(), j ) };
                      if ( index >= static_cast<int64_t>(n_elmts) ) {
                           cerr <<"\n\t"<< index <<" vs. number of elements = "<< n_elmts << endl;
                           csmp_error.Note( ERROR, "MeshManager::Initialise: ", "element ID in 'pfverts' out of range.");
                        }
                      // only if there is a neighbor something needs to be done; else nullptr was aready assigned
                      else if ( index >= 0 )
                        e.Assign( j, &(*next(elements_.begin(),index)) );
                  }
             }
        }
    }
  else
    csmp_error.Note( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity; nothing was done." );


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
       size_t   face_idx(vset.Elements());
       auto  first( vset.PlistFacesBegin() ), last( vset.PlistFacesEnd() );
       while ( first != last ) {
            const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( face_idx ));
            const FiniteVolumeStencil<dim>* const stencil_ptr = (initialise_FV_stencils==true) ?
                                                                 fvm_manager_->Stencil( csmpElementType ) :
                                                                 static_cast<const FiniteVolumeStencil<dim>* const>(nullptr);
            if ( csmpElementType == UNKNOWN ) {
                 cerr <<"\n\t"<< parseFiniteElementType(csmpElementType) <<" encountered for Face "<< face_idx <<"\n";
                 csmp_error.Note( FATAL_ERROR, "MeshManager::Initialise:", "encountered UNKNOWN Face element type." );
              }
            typename plf::colony<Face<dim>>::iterator
               fit = faces_.emplace( Face<dim>( face_idx, fem_manager_.E( csmpElementType ),
                                                          stencil_ptr, evars, cvars ) );
            // assigning nodes to faces
            const auto nodes( (*fit).Nodes() );
            for ( uint32_t j{0U}; j<nodes; ++j ) {
                const size_t node = vset.Plist( face_idx, j );
                assert( node < n_nodes );
                (*fit).Assign( j, &(*next(nodes_.begin(),static_cast<long>(node))) );
              }
            ++face_idx;
            ++first;
          }
       assert( faces_.size() == vset.Faces() );

       // connecting the faces to their equi- and higher-dimensional neighbors
       // --------------------------------------------------------------------
       // necessary info is stored in 'pfverts' record for Face
       if ( !vset.WithNeighbourConnectivity() )
         csmp_error.Note( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity for Face objects; nothing was done." );
       else
         {
            if ( csmp_error.Verbose() )
              cout << "\nMeshManager<" << dim << ">::Initialize: connecting faces to their equidimensional and higher-dimensional neighbors..." << endl;
            const size_t  n_faces(faces_.size());
            for ( auto& e : faces_ )
              {
                 // Equidimensional Face neighbors first
                 // ------------------------------------
                 const auto neighbors( e.Neighbors() );
                 for ( uint32_t j{0U}; j<neighbors; ++j )
                   {
                      // if there is a neighbor (as is the case if the stored index is greater than zero)
                      // (e->Idx() starts with elements=first face)
                      const int64_t  index( vset.Pfvert( e.Idx(), j ) );
                      if ( index >= static_cast<int64_t>(n_elmts+n_faces) ) {
                           cerr <<"\n\t"<< index <<" vs. number of elements+faces = "<< n_elmts + n_faces << endl;
                           csmp_error.Note( ERROR, "MeshManager::Initialise: ", "face ID in 'pfverts' out of range.");
                        }
                      if ( index >= static_cast<int64_t>(n_elmts) )
                        e.Assign( j, &(*next(faces_.begin(),index - static_cast<int64_t>(n_elmts))) );
                      else {
                           // if the Face neighbor has an index smaller than n_elmts it must be a boundary indicator
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
                 if ( index1 >= static_cast<int64_t>(n_elmts) ) {
                      cerr <<"\n\tFace "<< e.Idx() <<": "<< index1 <<" vs. "<< n_elmts <<" elements.\n";
                      csmp_error.Note( ERROR, "MeshManager::Initialise", "Index of first higher-dimensional element of Face out of range.");
                   }
                 else if ( index1 < 0 ) {
                      cerr <<"\n\tFace "<< e.Idx() <<": "<< index1 <<"\n";
                      csmp_error.Note( ERROR, "MeshManager::Initialise", "First higher-dimensional element out of range.");
                   }

                 // Outside neighbor 2: outer neighbor element will only be there if Face on INTERNAL model boundary
                 const int64_t  index2( vset.Pfvert( e.Idx(), neighbors + 1U ) );
                 if ( index2 >= static_cast<int64_t>(n_elmts) ) {
                      cerr <<"\n\tFace "<< e.Idx() <<": "<< index2 <<" vs. "<< n_elmts <<" elements.\n";
                      csmp_error.Note( ERROR, "MeshManager::Initialise", "Index of second higher-dimensional element of Face out of range.");
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
                 // else, this would be a nullptr that is already initialised as such
                 
               } // end face loop
              
             } // end else
        
    } // end construction and initialisation of Face objects


  // ------------------------------------------------------------------
  // 4. constructing Interfaces objects using the VSet node information
  // ------------------------------------------------------------------
  // (continuous running index 'idx' will also be used for interfaces)
  if ( vset.Interfaces() > 0 )
    {
       assert( vset.HybridElementTypeMesh() );
       if ( csmp_error.Verbose() )
         cout << "\nMeshManager<" << dim << ">::Initialize: assigning nodes to interfaces..." << endl;
       
       const LocalVariables evars( phys_vars.LocalVariablesAt( INTER_FACE ) );
       const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt( INTER_FACE ) );

       auto  first( vset.PlistInterFacesBegin() ), last( vset.PlistInterFacesEnd() );

       size_t interface_idx(vset.Elements() + vset.Faces());
       while ( first != last )
         {
            const CSMP_FEM_TYPE csmpElementType = static_cast<CSMP_FEM_TYPE>(vset.ElementType( interface_idx ));
            const FiniteVolumeStencil<dim>* const stencil_ptr = (initialise_FV_stencils==true) ?
                                                                 fvm_manager_->Stencil( csmpElementType ) :
                                                                 static_cast<const FiniteVolumeStencil<dim>* const>(nullptr);
            typename plf::colony<InterFace<dim>>::iterator
              ifit = interfaces_.emplace( InterFace<dim>( interface_idx,
                                                          fem_manager_.E( csmpElementType ),
                                                          stencil_ptr,
                                                          evars, cvars ) );
                                                       
            // number of nodes of the finite-element corresponding to the interface
            const auto nodes( (*ifit).FE()->Nodes() );
            // assigning nodes
            // inside
            for ( uint32_t j{0U}; j<nodes; ++j ) {
                 const size_t node = vset.Plist( interface_idx, j );
                 if ( node >= n_nodes ) {
                      cerr <<"\n\tInterFace "<< interface_idx <<": INSIDE node j "<< node <<" vs. "<< n_nodes <<" nodes.\n";
                      csmp_error.Note( ERROR, "MeshManager::Initialise", "Index of InterFace node out of range.");
                   }
                 (*ifit).Assign( j, &(*next(nodes_.begin(),static_cast<long>(node))), INSIDE );
              }
            // outside
            for ( uint32_t j{0U}; j<nodes; ++j ) {
                 const size_t node = vset.Plist( interface_idx, j+nodes );
                 if ( node >= n_nodes ) {
                      cerr <<"\n\tInterFace "<< interface_idx <<": OUTSIDE node j "<< node <<" vs. "<< n_nodes <<" nodes.\n";
                      csmp_error.Note( ERROR, "MeshManager::Initialise", "Index of InterFace node out of range.");
                   }
                 (*ifit).Assign( j, &(*next(nodes_.begin(),static_cast<long>(node))), OUTSIDE );
              }
            ++interface_idx;
            ++first;
          }
       assert( interfaces_.size() == vset.Interfaces() );

       // connecting the interfaces to their equi- and higher-dimensional neighbors
       // -------------------------------------------------------------------------
       // necessary info is stored in 'pfverts' record for InterFace: inner nbors first, then outer, then higher-dimensional ones
       if ( !vset.WithNeighbourConnectivity() )
         csmp_error.Note( WARNING, "MeshManager::Initialize:", "Input VSet does not contain any neighbor connectivity for InterFace objects; nothing was done." );
       else {
         if ( csmp_error.Verbose() )
           cout << "\nMeshManager<" << dim << ">::Initialize: connecting interfaces to their higher-dimensional neighbors..." << endl;
         const size_t  n_faces(faces_.size()), n_interfaces(interfaces_.size());
         const size_t  n_cells(n_elmts+n_faces+n_interfaces);
         // connecting interfaces to their higher-dimensional neighbors
         for ( auto& itf : interfaces_ )
           {
              const int64_t iface_idx = static_cast<int64_t>(itf.Idx());
              assert(  iface_idx >= static_cast<int64_t>(n_elmts + n_faces) );
              
              // 1. Assigning equidimensional InterFace-type neighbors first
              // -----------------------------------------------------------
              const auto neighbors{ itf.Neighbors() };
              for ( uint32_t j{0U}; j<neighbors; ++j )
                {
                   // if there is a neighbor (as is the case if the stored index is greater than zero)
                   const int64_t  index = vset.Pfvert( static_cast<size_t>(iface_idx), j );
                   
                   // if there is no neighbor nothing needs to be done because all neighbor pointers
                   // are already set to 'null' per default
                   if ( index < 0 ) continue;
                   
                   // if the index is out of range
                   if ( index >= static_cast<int64_t>(n_cells) || index <= static_cast<int64_t>(n_elmts) ) {
                        cerr <<"\n\t"<< index <<" vs. number of elements+faces+interfaces = "<< n_cells << endl;
                        csmp_error.Note( ERROR, "MeshManager::Initialise: ", "interface ID in 'pfverts' out of range.");
                     }

                   // NB: the interface number in the container is the number from the VSet - elements - faces
                   // because the interface container indexes from 0..n-1 (int64_t ~long used by next() )
                   const int64_t neighbor_idx = index - static_cast<int64_t>(n_elmts + n_faces);
                   itf.Assign( j, &(*next(interfaces_.begin(),neighbor_idx)) );
                }

             // 2. Assigning the higher-dimensional neighbor Element objects
             // ------------------------------------------------------------
             // (both higher-dimensional neighbors must be defined because interfaces exist only on the inside of models)
             const int64_t  index1 = vset.Pfvert( static_cast<size_t>(iface_idx), neighbors );
             const int64_t  index2 = vset.Pfvert( static_cast<size_t>(iface_idx), neighbors+1U );
             
             if ( index1 < 0 || index2 < 0 ) {
                  cerr <<"\n\tInterFace "<< iface_idx <<": inner neighbor "<< index1 <<" and outer neighbor "<< index2 <<"\n";
                  csmp_error.Note( ERROR, "MeshManager::Initialise: ", "Higher dimensional neighbor of InterFace not defined in 'pfverts'.");
               }
             if ( index1 >= static_cast<int64_t>(n_elmts) || index2 >= static_cast<int64_t>(n_elmts) ) {
                  cerr <<"\n\tInterFace "<< iface_idx <<": inner neighbor "<< index1 <<" and outer "<< index2 <<"\n";
                  csmp_error.Note( ERROR, "MeshManager::Initialise: ", "Higher dimensional neighbor indices of InterFace out of range.");
               }

             // assignment: inner and outer Element objects
             assert( index1 > MULTIPLE );
             assert( index2 > MULTIPLE );
             Element<dim>* const innerElement = &(*next(elements_.begin(),index1));
             Element<dim>* const outerElement = &(*next(elements_.begin(),index2));
             // assignment: local number of faces adjacent to InterFace; these face numbers must always be defined
             const auto inner_face_id = static_cast<uint32_t>(vset.Pfvert( static_cast<size_t>(iface_idx), neighbors+2U ));
             const auto outer_face_id = static_cast<uint32_t>(vset.Pfvert( static_cast<size_t>(iface_idx), neighbors+3U ));
             assert( inner_face_id < innerElement->Faces() );
             assert( outer_face_id < outerElement->Faces() );
             itf.Assign( innerElement, inner_face_id, outerElement, outer_face_id );
             
             // assignment: intervening Element else boundary flag INTERNAL
             const int64_t  index3 = vset.Pfvert( static_cast<size_t>(iface_idx), neighbors+4U );
             assert( index3 < static_cast<int64_t>(n_elmts) );
             assert( index3 > MULTIPLE );
             Element<dim>* const middleElement = (index3 < 0) ? nullptr : &(*next(elements_.begin(),index3));
             if ( middleElement ) itf.Assign( middleElement );
          }
      }
      
    } // end construction of Interface objects


  // ---------------------------------------------------------------------
  // 5. Flagging nodes at model boundary with BOX_BOUNDARY flags
  // ---------------------------------------------------------------------
  // NB: the node flags were already assigned further above where the nodes were constructed!
    
    
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
  size_t i{0ul};
  for ( auto& n : nodes_ )
    n.ReserveParentStorage( parent_elmts_per_node[i++] );

   // assigning the parent element information to the nodes
   for ( auto& e : elements_ )
       for ( uint32_t j{0U}; j<e.Nodes(); ++j )
         e.N( j )->Assign( j, &e );
       
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
        for ( auto& s : pnode[node] )
          nptrs.push_back( &(*next(nodes_.begin(),static_cast<long>(s))) );
        // assigning the neighbor pointers to the node
        n.Assign( nptrs, true );
        node++;
     }

     
  // ------------------------------------------------------------------------------
  // 8. Constructing NodeManifolds if any
  // -------------------------------------------------------------------------------
   if ( !interfaces_.empty() )
     node_manifold_manager_ = new NodeManifoldManager<dim>( nodes_, vset.PmanifoldsBegin(), vset.PmanifoldsEnd() );
     
#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
if ( !interfaces_.empty() ) {
   cerr <<"\n\n"<<"\nMeshManager::Initialise: node manifolds initialised in NodeManifoldManager:\n";
   size_t counter{0U};
   for ( auto nit=node_manifold_manager_->ManifoldsBegin(); nit!=node_manifold_manager_->ManifoldsEnd(); ++nit ) {
        cerr <<"\n\t\t"<< counter++ <<": "<< parse( (*nit).GeometricClassifier() ) <<" ";
        for ( auto z{0U}; z<(*nit).Branches(); z++ )
          cerr << (*nit).N(z)->Idx() <<" ";
     }
    cerr << endl;
  }
#endif


// DEBUGGING VData - Element connectivity (ok for model FracBox
// ---------------------------------------------------------------------------------
/*
cout <<"\n"<<"MeshManager::Intialise(database,vset,FV): DEBUGGING neighbor connectivity: "<< endl;
// do the neighbor types match
for ( auto it=elements_.begin(); it!=elements_.end(); ++it ) {
     for ( auto eit=(*it).NeighborsBegin(); eit!=(*it).NeighborsEnd(); ++eit )
       // if the face has a neighbor
       if ( (*eit) )
         // if that neighbor is connected to the element
         for ( auto neit=(*eit)->NeighborsBegin(); neit!=(*eit)->NeighborsEnd(); ++neit )
           if ( (*neit) && (*neit)->Idx() == (*it).Idx() ) { //&& (*neit)->FE_Type() != (*eit)->FE_Type() ) {
                cout <<"\n"<<"nbor "<< (*eit)->Idx() <<": "<< parseAbbreviated_FE_Type((*eit)->FE_Type()) <<" is connected to ";
                cout << (*neit)->Idx() <<": "<< parseAbbreviated_FE_Type((*neit)->FE_Type()) <<" !";
                assert( (*neit)->FE_Type() == (*it).FE_Type() );
             }
  }
*/
   return true;
  
} // end Initialise








/**
  Initialises the finite volume policy of the elements, faces, and interfaces by assigning the finite volume stencil pointers of the elements to the
  corresponding finite volume stencils.
  
  @param assign_stencils_to_elements connects the stencils in the manager with the individual elements, faces and interfaces

@attention If a stencil is assigned already, noting is done. Remove stencil first (NULL ptr in elements)

*/
template<uint32_t dim>
void MeshManager<dim>::InitializeFiniteVolumeStencils( const PropertyDatabase<dim>& pref, bool assign_stencils_to_elements )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    // if the FV stencil manager has already been initialised
    if ( fvm_manager_ )  {
         csmp_error.Note( INFO, "MeshManager<dim>::InitializeFiniteVolumeStencils",
                         "The finite volume stencils were already initialised earlier.");
         return;
      }

    if ( !fem_manager_.UsesElementsWithLocalCoordinateSystem() ) {
         csmp_error.Note( FATAL_ERROR, "MeshManager<dim>::InitializeFiniteVolumeStencils",
                         "Currently FV stencils exist in parametric space, requiring FE with a local coordinate system.");
         return;
      }

    if ( fem_manager_.InterpolationOrder() != 1U ) {
         csmp_error.Note( FATAL_ERROR, "MeshManager<dim>::InitializeFiniteVolumeStencils",
                         "Currently FV stencils exist only for finite elements with linear interpolation functions.");
         return;
      }

    // 1. initialize the stencil manager (the stencils are build and assigned the correct properties
    fvm_manager_ = new FiniteVolumeStencilManager<dim>( fem_manager_ );

  if ( assign_stencils_to_elements )
   {
      if ( elements_.empty() ) {
           csmp_error.Note( ERROR, "MeshManager<dim>::InitializeFiniteVolumeStencils",
                           "Currently no Element objects exist to which stencils could be assigned.");
           return;
        }

      // 2. Now the stencil pointers in each finite element are connected to the correct corresponding stencils and update variable
      // storage for fv integration (sector/facet) point properties
      {
        const LocalVariables lvs( pref.LocalVariablesAt(ELEMENT) );
        const IntegrationPointVariables ipvs( pref.IntegrationPointVariablesAt(ELEMENT) );
        for ( auto& e : elements_ ) {
              if ( !e.FV() ) {
                  e.AssignFiniteVolume( fvm_manager_->Stencil( e.FE_Type() ) );
                  e.ResizePropertyStorage( lvs, ipvs );
                }
          }
      }
      // faces
      if ( !faces_.empty() )
        {
          const LocalVariables lvs( pref.LocalVariablesAt(ELEMENT) );
          const IntegrationPointVariables ipvs( pref.IntegrationPointVariablesAt(ELEMENT) );
          for ( auto& e : faces_ ) {
                if ( !e.FV() ) {
                    e.AssignFiniteVolume( fvm_manager_->Stencil( e.FE_Type() ) );
                    e.ResizePropertyStorage( lvs, ipvs );
                  }
            }
        }
      // interfaces
      if ( !interfaces_.empty() )
        {
          const LocalVariables lvs( pref.LocalVariablesAt(ELEMENT) );
          const IntegrationPointVariables ipvs( pref.IntegrationPointVariablesAt(ELEMENT) );
          for ( auto& e : interfaces_ ) {
                if ( !e.FV() ) {
                    e.AssignFiniteVolume( fvm_manager_->Stencil( e.FE_Type() ) );
                    e.ResizePropertyStorage( lvs, ipvs );
                  }
            }
        }
    }

 } // end InitializeFiniteVolumeStencils







  // ==============================================================
  //
  // MESH MODIFICATION (following methods)
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
    Inserts a new Node at the desired location.  No connection are made.
*/
template<uint32_t dim>
Node<dim>* const MeshManager<dim>::AddNodeAt( const Point<dim>& location,
                                              const LocalVariables& lvars,
                                              BOX_BOUNDARY bdry,
                                              TOPOTYPE topo )
{
   size_t node_idx{ nodes_.size() };
   typename plf::colony<Node <dim>>::iterator
     nit = nodes_.emplace( Node<dim>( nodes_.size(), location, lvars, bdry, topo ) );
     
   (*nit).Idx(node_idx);
   
   return &(*nit);
}







/**
   Constructs element and connects it up with the supplied nodes.
   Method also tries to find neighbor elements using the parent element connectivity of the nodes.
   A warning message is issued if there are issues with the input data or the neighbors cannot be identified.
   
   @author SKM
   @date 17/9/21
*/
template<uint32_t dim>
Element<dim>*	const MeshManager<dim>::AddElement( CSMP_FEM_TYPE etype,
                                                  const LocalVariables& lvars,
                                                  const IntegrationPointVariables& ivars,
                                                  const vector<Node<dim>*>& nodes,
                                                  int32_t material_id )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // verifying input node vector
   if ( nodes.empty() )
     csmp_error.Note( ERROR, "MeshManager<dim>::AddElement", "node vector is empty");
          
   // 1. constructing new element
   // ---------------------------
   size_t elmt_idx{ elements_.size() };
   typename plf::colony<Element<dim>>::iterator eit = ( fvm_manager_ ) ?
               elements_.emplace( Element<dim>( elements_.size(),
                                  fem_manager_.E(etype), fvm_manager_->Stencil(etype), lvars, ivars, material_id ) ) :
               elements_.emplace( Element<dim>( elements_.size(),
                                  fem_manager_.E(etype), static_cast<FiniteVolumeStencil<dim>*>(nullptr), lvars, ivars, material_id ) );
   (*eit).Idx( elmt_idx );
   
   // 2. assigning nodes
   // ------------------
   // verifying that the correct number of nodes is supplied
   if ( nodes.size() != (*eit).FE()->Nodes() )
     csmp_error.Note( ERROR, "MeshManager<dim>::AddElement", "supplied node vector has incorrect size");
   
   // node assignment
   const auto n_nodes{ nodes.size() };
   for ( uint32_t i{0U}; i<n_nodes; ++i ) {
        (*eit).Assign( i, nodes[i] );
        // assigning the element as a new parent element of the node
        nodes[i]->Assign( i, &(*eit) );
     }

  // 3. trying to establish neighbor information from the nodes assuming if they have parent connectivity
  // ----------------------------------------------------------------------------------------------------
  //    checking whether the nodes have the necessary parent element information
  bool valid_parent_info{true};
  for ( uint32_t i{0U}; i<n_nodes; ++i )
    if ( (*eit).N(i)->Parents() == 0U ) {
         cerr <<"\n\tnode "<< i;
         valid_parent_info = false;
         csmp_error.Note( WARNING, "MeshManager<dim>::AddElement",
                          "neighbor information could not be created because node(s) miss parent element info");
      }
      
   // 4. connect new element to its neighbors: this gets done only after a range of elements were created.
     
   return &(*eit);
  
} // AddElement







/**
     Puts a lower-dimensional element inside of an InterFace, connecting it to its base pointer
     
     @attention the nodes need to be provided because they are shared among the lower-dimensional elements and collocated
          so that they cannot be told apart.
          
          @attention no neighbor connectivity is provided here because it is not necessarily known when this metho gets called
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
   if ( ifptr == nullptr ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::AddInterveningElement", "finite element pointer not initialised");
        return nullptr;
     }
   if ( ifptr->HasInterveningElement() ) {
        ifptr->InterveningElement()->Out();
        csmp_error.Note( ERROR, "MeshManager<dim>::AddInterveningElement", "InterFace already has intervening element");
        return ifptr->InterveningElement();
     }

   // 1. checking the node vector
   if ( nodes.empty() )
     csmp_error.Note( ERROR, "MeshManager<dim>::AddInterveningElement", "node vector is empty");
   if ( nodes.size() != ifptr->FE()->Nodes() )
     csmp_error.Note( ERROR, "MeshManager<dim>::AddInterveningElement", "node vector has the wrong size");
     
   // 2. checking the validity of the node vector in debug mode
#ifdef DEBUG
   // node vector
   for ( uint32_t i{0U}; i<ifptr->FE()->Nodes(); ++i ) {
         if ( nodes[i] == nullptr ) {
              cerr <<"\n\tnode "<< i;
              csmp_error.Note( ERROR, "MeshManager<dim>::AddInterveningElement", "node vector contains a nullptr");
              break;
           }
         else if ( nodes[i]->Coordinate() != ifptr->N(i,INSIDE)->Coordinate() ) {
              cerr <<"\n\tnode "<< i;
              csmp_error.Note( ERROR, "MeshManager<dim>::AddInterveningElement", "node locations do not match");
              break;
           }
       }
#endif

   // 3. creating the new Element
   size_t elmt_idx{ elements_.size() };
   typename plf::colony<Element<dim>>::iterator
     eit = elements_.emplace( Element<dim>( elements_.size(), ifptr->FE(), ifptr->FV(), lvars, ivars, material_id ) );
     
   (*eit).Idx( elmt_idx );
   
   // 4. connecting the nodes to the element
   const uint32_t n_nodes{ ifptr->FE()->Nodes() };
   for ( uint32_t i{0U}; i<n_nodes; ++i )
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
Face<dim>* const MeshManager<dim>::ReplaceElementByFace( typename vector<Element<dim>*>::iterator eptr_it,
                                                         csmp::Element<dim>* inner_eptr,
                                                         csmp::Element<dim>* outer_eptr,
                                                         uint32_t adjacent_face_of_inner_element,
                                                         uint32_t adjacent_face_of_outer_element,
                                                         const LocalVariables& lvars,
                                                         const IntegrationPointVariables& ivars,
                                                         bool delete_original_face )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( (*eptr_it) == nullptr )
     csmp_error.Note( ERROR, "MeshManager<dim>::ReplaceElementByFace", "element pointer not initialised");

   // is the element indeed lower dimensional?
   if constexpr ( dim == 3 )
     if ( !(*eptr_it)->IsSurface() )
     csmp_error.Note( ERROR, "MeshManager<3>::ReplaceElementByFace", "element to be replaced is not a lower-dimensional surface element");
   if constexpr ( dim == 2 )
     if ( !(*eptr_it)->IsLine() )
     csmp_error.Note( ERROR, "MeshManager<2>::ReplaceElementByFace", "element to be replaced is not a lower-dimensional line element");

   if ( inner_eptr == nullptr )
     csmp_error.Note( ERROR, "MeshManager<dim>::ReplaceElementByFace", "pointer to higher dimensional element on inside not initialised");
   if ( inner_eptr == outer_eptr ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::ReplaceElementByFace", "cannot create Face"
                                  "pointer to higher dimensional elements are the same");
        return nullptr;
     }
       
   assert( adjacent_face_of_inner_element < inner_eptr->Faces() );
   if ( outer_eptr != nullptr ) assert( adjacent_face_of_outer_element < outer_eptr->Faces() );

   // 1. constructing new face
   const size_t face_id = faces_.size(); // since the face will be added at the end of the colony
   typename plf::colony<Face<dim>>::iterator
     fit = faces_.emplace( Face<dim>( *(*eptr_it), inner_eptr, outer_eptr,
                                      adjacent_face_of_inner_element, adjacent_face_of_outer_element,
                                      lvars, ivars ) );
   (*fit).Idx( face_id );
   
   // 2. removing the element from the parent-element vectors of its nodes before deleting it
   if ( delete_original_face ) Delete( eptr_it );
   return &(*fit);
   
 } // end ReplaceElementByFace
       



/**
  Replaces a lower dimensional element with an InterFace object and deletes lower dim element
  @attention Construction process assumes Nodes are ALREADY duplicated and assigned to the Inner and Outer Parent.

  @brief Performs input parameter checks. Constructs InterFace object using consgtructor and places in interfaces_ container.
  Deletes lower dimensional element and returns Interface object pointer
  Neighbor connectivity of parent elements are updated as they are unnasigned from each other (happens during interface construction).
*/
template<uint32_t dim>
InterFace<dim>* const MeshManager<dim>::ReplaceElementByInterFace( typename std::vector<Element<dim>*>::iterator eptr_it,
                                                                   csmp::Element<dim>* inner_eptr,
                                                                   csmp::Element<dim>* outer_eptr,
                                                                   uint32_t adjacent_face_of_inner_element,
                                                                   uint32_t adjacent_face_of_outer_element,
                                                                   const LocalVariables& lvars,
                                                                   const IntegrationPointVariables& ivars )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // ----------------------
   // pointers
   if ( (*eptr_it) == nullptr ) {
         csmp_error.Note( ERROR, "MeshManager<dim>::ReplaceElementByInterFace", "element pointer not initialised");
         return nullptr;
     }
   // is the element indeed lower dimensional?
   if constexpr ( dim == 3U ) if ( !(*eptr_it)->IsSurface() )
     csmp_error.Note( ERROR, "MeshManager<3>::ReplaceElementByInterFace", "element to be replaced is not a lower-dimensional surface element");
   if constexpr ( dim == 2U ) if ( !(*eptr_it)->IsLine() )
     csmp_error.Note( ERROR, "MeshManager<2>::ReplaceElementByInterFace", "element to be replaced is not a lower-dimensional line element");

   if ( inner_eptr == nullptr )
     csmp_error.Note( ERROR, "MeshManager<dim>::ReplaceElementByInterFace", "pointer to higher dimensional element on inside not initialised");

   if ( inner_eptr == outer_eptr ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::ReplaceElementByInterFace", "cannot create InterFace"
                                "pointers to higher dimensional elements are the same");
        return nullptr;
     }
   assert( adjacent_face_of_inner_element < inner_eptr->Faces() );
   if ( outer_eptr != nullptr ) assert( adjacent_face_of_outer_element < outer_eptr->Faces() );
   

   // 1. constructing the new interface
   // ---------------------------------
   const size_t iface_id = interfaces_.size(); // since the interface will be added at the end of the colony
   // creates interface, detaching InnerParent and OuterParent elements that share a face with the interface from each other
   typename plf::colony<InterFace<dim>>::iterator
     fit = interfaces_.emplace( InterFace<dim>( *(*eptr_it), inner_eptr, outer_eptr, // nodes from higher dim arent face should be already duplicated
                                                adjacent_face_of_inner_element, adjacent_face_of_outer_element,
                                                lvars, ivars ) );
   (*fit).Idx( iface_id );

   // 3. deleting the original Element
   // --------------------------------
   Delete( eptr_it );

   return &(*fit);
   
 } // end ReplaceElementByInterFace
       


       
       
       
       
       
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
     csmp_error.Note( ERROR, "MeshManager<dim>::AddFace", "pointer to higher dimensional element on inside not initialised");
   if ( outer_parent == nullptr )
     csmp_error.Note( INFO, "MeshManager<dim>::AddFace", "pointer to higher dimensional element on ouside not initialised");
   if ( outer_parent == inner_parent ) {
        csmp_error.Note( INFO, "MeshManager<dim>::AddFace", "cannot create Face",
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
                                                const vector<Node<dim>*>& nodes )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( adjacent_face1 == nullptr ) {
       csmp_error.Note( ERROR, "MeshManager<dim>::AddFace", "pointer to higher dimensional Face1 on inside not initialised");
       return nullptr;
    }
   if ( adjacent_face1 == adjacent_face2 ) {
       csmp_error.Note( ERROR, "MeshManager<dim>::AddFace", "Face pointers point to same Face", "Nothing was done");
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







/**
    After testing whether the provided element face lies on the model boundary (has no neigbhor),
    method caps it with a Face object that has the same normal orientation as the Element face.
*/
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
     csmp_error.Note( ERROR, "MeshManager<dim>::AddBoundaryFace", "element pointer not initialised");
   if ( local_face_id >= eptr->Faces() )
     csmp_error.Note( ERROR, "MeshManager<dim>::AddBoundaryFace", "face ID does not exist in element");
   if ( eptr->Neighbor(local_face_id) != nullptr )
     csmp_error.Note( WARNING, "MeshManager<dim>::AddBoundaryFace", "element face has a neighbor; is it really located at model boundary?");

   if ( eptr->IsEquidimensional() == false )
     csmp_error.Note( WARNING, "MeshManager<dim>::AddBoundaryFace", "element must have the same dimension as the model for this method to work");

   // is face indeed a boundary face?
   if ( eptr->Neighbor(local_face_id) != nullptr ) {
        cout <<"\n"<<"Element "<< eptr->Idx() <<" "<< parseAbbreviated_FE_Type(eptr->FE_Type()) <<": face "<< local_face_id << endl;
        cout <<"\t"<<"outside element: "<< eptr->Neighbor(local_face_id)->Idx() <<": ";
        cout << parseAbbreviated_FE_Type( eptr->Neighbor(local_face_id)->FE_Type() ) << endl;
        csmp_error.Note( WARNING, "MeshManager<dim>::AddBoundaryFace", "requested element face does not lie on model boundary");
        // advancing Face Idx until face without neighbor is found
        bool boundary_face_found{ false };
        local_face_id = 0u;
        while ( local_face_id < eptr->Faces() ) {
              if ( eptr->Neighbor(local_face_id) == nullptr ) {
                   boundary_face_found = true;
                   break;
                }
              local_face_id++;
          }
        if ( !boundary_face_found ) {
             csmp_error.Note( ERROR, "MeshManager<dim>::AddBoundaryFace", "requested element has no boundary faces");
             return nullptr;
          }
     }

   // 1. constructing new face, connecting it to its higher-dimensional neighbor on the inside, and assigning nodes
   const size_t global_face_number{faces_.size()};
   typename plf::colony<Face<dim>>::iterator
     fit = faces_.emplace( Face<dim>( *eptr, fem_manager_.E( eptr->FE()->ElementTypeOfFace(local_face_id) ),
                                      fvm_manager_, local_face_id, lvars, ivars ) );
   (*fit).Idx( global_face_number );

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
                                                      const IntegrationPointVariables& ivars,
                                                      vector<Node<dim>*> outside_nodes )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // pointers
   if ( inner_parent == nullptr )
     csmp_error.Note( ERROR, "MeshManager<dim>::AddInterFace", "pointer to higher dimensional element on inside not initialised");
   if ( outer_parent == nullptr )
     csmp_error.Note( ERROR, "MeshManager<dim>::AddInterFace", "pointer to higher dimensional element on ouside not initialised");
   if ( inner_parent == outer_parent ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::AddInterFace", "cannot create InterFace",
                                "inner and outer parent pointers are the same" );
        return nullptr;
     }

   assert( inner_element_face_id < inner_parent->Faces() );
   assert( outer_element_face_id < outer_parent->Faces() );

   // 1. establishing the finite element type if the interface
   const CSMP_FEM_TYPE etype = inner_parent->FE()->ElementTypeOfFace( inner_element_face_id );

   // 2. constructing new interface
   const size_t iface_id = interfaces_.size();
   const FiniteVolumeStencil<dim>* const stencil_ptr = (fvm_manager_) ? fvm_manager_->Stencil(etype) : nullptr;
   typename plf::colony<InterFace<dim>>::iterator
     ifp = interfaces_.emplace( InterFace<dim>( fem_manager_.E(etype), stencil_ptr, lvars, ivars ) );
     
   // 3. assigning higher dimensional elements and faces
   (*ifp).Assign( inner_parent, inner_element_face_id, outer_parent, outer_element_face_id );
   (*ifp).Idx( iface_id );
   
   // 4. assigning the inside nodes to the new InterFace (which are those of the face of the inside element)
   uint32_t n_count{0U};
   for ( const auto& n : inner_parent->FE()->NodesOfFace( inner_element_face_id ) )
     (*ifp).Assign( n_count++, inner_parent->N(n), INSIDE );
   
   // 5. assigning the outside nodes to the InterFace
   const auto n_nodes_per_face{ outside_nodes.size() };
   for ( auto n{0U}; n<n_nodes_per_face; ++n )
     (*ifp).Assign( n, outside_nodes[n], OUTSIDE );

   // 6. replacing the original nodes of the face of the outside element with the new nodes
   n_count = 0U;
   for ( const auto& n : outer_parent->FE()->NodesOfFace( outer_element_face_id ) )
     outer_parent->Assign( n, outside_nodes[n_count++] );
     
   // 7. detaching the input Elements from one-another
   inner_parent->Unassign( outer_parent );
   outer_parent->Unassign( inner_parent );

   return &(*ifp);

} // end AddInterFace






/**
  Replaces a lower dimensional element with an InterFace object and deletes lower dim element
  @attention Construction process assumes Nodes are ALREADY duplicated and assigned to the Inner and Outer Parent.

  @brief Performs input parameter checks. Constructs InterFace object using consgtructor and places in interfaces_ container.
  Deletes lower dimensional element and returns Interface object pointer
  Neighbor connectivity of parent elements are updated as they are unnasigned from each other (happens during interface construction).
*/
template<uint32_t dim>
InterFace<dim>* const MeshManager<dim>::WrapInterFaceAroundElement( csmp::Element<dim>* const eptr,
                                                                    csmp::Element<dim>* inner_eptr,
                                                                    csmp::Element<dim>* outer_eptr,
                                                                    uint32_t adjacent_face_of_inner_element,
                                                                    uint32_t adjacent_face_of_outer_element,
                                                                    const LocalVariables& lvars,
                                                                    const IntegrationPointVariables& ivars )
 {
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
   // 0. verifying the input
   // ----------------------
   // pointers
   if ( eptr == nullptr ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::WrapInterFaceAroundElement", "element pointer not initialised");
        return nullptr;
     }
   // is the element indeed lower dimensional?
   if constexpr ( dim == 3U ) if ( !eptr->IsSurface() )
     csmp_error.Note( ERROR, "MeshManager<3>::WrapInterFaceAroundElement", "element to be replaced is not a lower-dimensional surface element");
   if constexpr ( dim == 2U ) if ( !eptr->IsLine() )
     csmp_error.Note( ERROR, "MeshManager<2>::WrapInterFaceAroundElement", "element to be replaced is not a lower-dimensional line element");

   if ( inner_eptr == nullptr )
     csmp_error.Note( ERROR, "MeshManager<dim>::WrapInterFaceAroundElement", "pointer to higher dimensional element on inside not initialised");

   if ( inner_eptr == outer_eptr ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::WrapInterFaceAroundElement", "cannot create InterFace"
                                "pointers to higher dimensional elements are the same");
        return nullptr;
     }
   assert( adjacent_face_of_inner_element < inner_eptr->Faces() );
   if ( outer_eptr != nullptr ) assert( adjacent_face_of_outer_element < outer_eptr->Faces() );
   

   // 1. constructing the new interface
   // ---------------------------------
   const size_t iface_id = interfaces_.size(); // since the interface will be added at the end of the colony
   // creates interface, detaching InnerParent and OuterParent elements that share a face with the interface from each other
   typename plf::colony<InterFace<dim>>::iterator
     fit = interfaces_.emplace( InterFace<dim>( *eptr, inner_eptr, outer_eptr, // nodes from higher dim arent face should be already duplicated
                                                adjacent_face_of_inner_element, adjacent_face_of_outer_element,
                                                lvars, ivars ) );
   (*fit).Idx( iface_id );
   (*fit).Assign( eptr ); // intervening element

   return &(*fit);
   
 } // end WrapInterFaceAroundElement









/**
   Relying on the parent element information from its corner nodes, method tries finds  neighbor elements for each element Face (or boundary).
   If a neighbor can be found, the function connects the element to it.
 
   @return the number of neighbors that were identified
 */
template<uint32_t dim>
uint32_t MeshManager<dim>::ConnectNeighborsUsingNodeParents( Element<dim>* const eptr )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     if ( eptr == nullptr ) {
             csmp_error.Note( ERROR, "MeshManager<dim>::ConnectNeighborsUsingNodeParents:",
                               "invalid element pointer" );
            return 0U;
        }
    
     vector<Node<dim>*> face_nodes;
     uint32_t nbors_found{0u};

     const auto n_faces(eptr->Faces());
     for ( uint32_t i{0U}; i<n_faces; ++i )
       {
          // finding the corner nodes of the face
          const auto crn_nodes = eptr->CornerNodesOfFace(i);
          face_nodes.assign( crn_nodes.begin(), crn_nodes.end() );
          
          // finding same-dimensional parent elements that share the elements face nodes (there should only be one!)
          vector<Element<dim>*>  elmts_with_all_nodes = parentElementsContaining<dim>( face_nodes.begin(), face_nodes.end() );
          eraseDifferentDimensionalOrInvalidCells( eptr, elmts_with_all_nodes );
          
          // if there are none nothing can be done, but there should not be more than one
          if ( elmts_with_all_nodes.size() <= 1 ) {
               if ( eptr != elmts_with_all_nodes[0] )
                 csmp_error.Note( ERROR, "MeshManager<dim>::ConnectNeighborsUsingNodeParents:", "nodes identified do not match element found" );
               continue;
            }
          else if ( elmts_with_all_nodes.size() > 2 ) {
               cerr <<"\n\t"<< eptr->Idx() <<": face: "<< i;
               csmp_error.Note( ERROR, "MeshManager<dim>::ConnectNeighborsUsingNodeParents:", "faces can only have one equidimensional neighbor" );
            }
            
          // assigning the neighbor element
          uint32_t nbor_elmt = ( eptr != elmts_with_all_nodes[0] ) ? 0u : 1u;
          eptr->Assign( i, elmts_with_all_nodes[nbor_elmt] );

          // finding which face of the identified neighbor element contains the shared nodes
          uint32_t face_of_nbor = faceWithCornerNodes( elmts_with_all_nodes[nbor_elmt], face_nodes.begin(), face_nodes.end() );
          
          // assigning the current element to the neighbor
          elmts_with_all_nodes[nbor_elmt]->Assign( face_of_nbor, eptr );
          
          nbors_found++;
       }
       
    return nbors_found;
        
 } // end ConnectNeighborsUsingNodeParents









/**
     Loops over  valid cell neighbors  and sets their neighbor pointers to  this cell to nullptr.
     
     @attention Never call this method in a Cell deletion loop where the next cell might be the one you just deleted!
 */
template<uint32_t dim>
template<template<uint32_t> class CELL>
void MeshManager<dim>::DetachNeighborsFrom( CELL<dim>* const eptr )
 {
    assert( eptr != nullptr );
    // nulling the connections of neighbor neighbor elements to this element
    // (neighbor pointer to this element is nulled)
    const uint32_t n_nbors{ eptr->Neighbors() };
    for ( uint32_t i{0U}; i<n_nbors; ++i ) {
          if ( eptr->Neighbor(i) ) {
              const uint32_t n_nbor_nbors{ eptr->Neighbor(i)->Neighbors() };
              for ( uint32_t j{0U}; j<n_nbor_nbors; ++j )
                if ( eptr->Neighbor(i)->Neighbor(j) == eptr ) {
                     // detaching neighbor from cell
                     eptr->Neighbor(i)->UnassignNeighbor( j );
                     break;
                  }
              // detaching itself from neighbor (only necessary if there is one) : works!
// NOT NECESSARY  eptr->UnassignNeighbor( i );
           }
       }
 }

template void MeshManager<1>::DetachNeighborsFrom( Element<1>* const );
template void MeshManager<2>::DetachNeighborsFrom( Element<2>* const );
template void MeshManager<3>::DetachNeighborsFrom( Element<3>* const );

template void MeshManager<1>::DetachNeighborsFrom( Face<1>* const );
template void MeshManager<2>::DetachNeighborsFrom( Face<2>* const );
template void MeshManager<3>::DetachNeighborsFrom( Face<3>* const );

template void MeshManager<1>::DetachNeighborsFrom( InterFace<1>* const );
template void MeshManager<2>::DetachNeighborsFrom( InterFace<2>* const );
template void MeshManager<3>::DetachNeighborsFrom( InterFace<3>* const );




/**
     Deletes Node without disconnecting its neighbors from it and
     without altering its parent element vector.
     
     @note any potential connections of neighbors to the node need be removed before potentially neighboring nodes start to get deleted,
     else already deleted nodes might get accessed again causing undefined behaviour.
*/
template<uint32_t dim>
auto	MeshManager<dim>::Delete( typename vector<Node<dim>*>::iterator nptr_ref ) -> typename plf::colony< Node<dim> >::iterator
 {
    // TODO: check that the node is not connected to elements anymore
    //     throw csmp::Exception( ERROR, "MeshManager<dim>::Delete(Node-ptr-ref)",
    //                           "Node  still belongs to Element(s). Therefore it cannot be deleted");
 
    auto pfl_it = nodes_.get_iterator( (*nptr_ref) );
    if ( pfl_it != nodes_.end() ) {
        (*nptr_ref) = nullptr;
         return nodes_.erase( pfl_it );
      }
    return nodes_.end();
 }



/**
      Remove selement from the parent-element vectors of its nodes
      (this way this gets done only once for each node concerned).
      No other connectivity changes are made because these would have side effects.
      
      @param eptr_ref reference to a pointer that is passed by reference so that it can be nulled
      @return iterator to next element in colony or end() if element could not be found.
      
      @attention Do not detach neighbour cells within this method (=sawing-off the branch you are sitting on).
*/
template<uint32_t dim>
auto	MeshManager<dim>::Delete( typename vector<Element<dim>*>::iterator eptr_ref ) -> typename plf::colony< Element<dim> >::iterator
 {
    assert( (*eptr_ref) != nullptr );
   // remove element from the parent element list of its connected nodes
   for ( auto nit=(*eptr_ref)->NodesBegin(); nit!=(*eptr_ref)->NodesEnd(); ++nit )
     (*nit)->Unassign( (*eptr_ref) );

    // geting iterator to element ('get_iterator' returns end, if element cannot be found)
    auto pfl_it = elements_.get_iterator( (*eptr_ref) );
    // if the element exists in the colony
    if ( pfl_it != elements_.end() ) {
         (*eptr_ref) = nullptr;
         // deleting element and returning colony iterator to next element in colony
         return elements_.erase( pfl_it );
      }
    // does an element after to the one erased exist in the colony? - if so, an iterator to it is returned
    return elements_.end();
 }

/* CUT OUT: ALL THESE HAVE SIDE EFFECTS when involved regions touch eachother
   // remove element from the parent element list of its connected nodes
   for ( auto nit=(*eptr_ref)->NodesBegin(); nit!=(*eptr_ref)->NodesEnd(); ++nit )
     (*nit)->Unassign( (*eptr_ref) );
*/



/**
     Deletion of a Face only affects its face neighbors

      @attention Do not detach neighbour cells within this method (=sawing-off the branch you are sitting on).
*/
template<uint32_t dim>
auto	MeshManager<dim>::Delete( typename vector<Face<dim>*>::iterator fptr_ref ) -> typename plf::colony< Face<dim> >::iterator
 {
    assert( (*fptr_ref) != nullptr );
    auto pfl_it = faces_.get_iterator( (*fptr_ref) );
    if ( pfl_it != faces_.end() ) {
         (*fptr_ref) = nullptr;
         return faces_.erase( pfl_it );
      }
    return faces_.end();
 }


/**
     Deletion of an InterFace only affects its face neighbors
      
      @attention Do not detach neighbour cells within this method (=sawing-off the branch you are sitting on).
*/
template<uint32_t dim>
auto	MeshManager<dim>::Delete( typename vector<InterFace<dim>*>::iterator fptr_ref ) -> typename plf::colony< InterFace<dim> >::iterator
 {
    assert( (*fptr_ref) != nullptr );
    auto pfl_it = interfaces_.get_iterator( (*fptr_ref) );
    if ( pfl_it != interfaces_.end() ) {
         (*fptr_ref) = nullptr;
         return interfaces_.erase( pfl_it );
      }
    return interfaces_.end();
 }



/**
   Unassigns any non-nullpointer neighbors the element type of which is unknown (global operation on all cells)
   @return number of such neighbors that were found and unassigned.
   
   @attention Method cannot be used because if such elements exist there already is undefined behaviour!
*/
/*
template<uint32_t dim>
template<template<uint32_t> class CELL>
size_t MeshManager<dim>::RemoveDegenerateNeighbors()
 {
    size_t false_nbors{0ul};

    if constexpr( is_same<CELL<dim>,Element<dim>>::value ) {
      for ( auto& cell : elements_ )
        for ( uint32_t i{0u}; i<cell.Neighbors(); ++i )
          if ( cell.Neighbor(i) )
            if ( !cell.Neighbor(i)->FE() || cell.Neighbor(i)->FE_Type() == UNKNOWN ) {
                 cell.UnassignNeighbor( i );
                 false_nbors++;
              }
      }
    else if constexpr( is_same<CELL<dim>,Face<dim>>::value ) {
      for ( auto& cell : faces_ )
        for ( uint32_t i{0u}; i<cell.Neighbors(); ++i )
          if ( cell.Neighbor(i) )
            if ( !cell.Neighbor(i)->FE() || cell.Neighbor(i)->FE_Type() == UNKNOWN ) {
                 cell.UnassignNeighbor( i );
                 false_nbors++;
              }
      }
    else if constexpr( is_same<CELL<dim>,InterFace<dim>>::value ) {
      for ( auto& cell : interfaces_ )
        for ( uint32_t i{0u}; i<cell.Neighbors(); ++i )
          if ( cell.Neighbor(i) )
            if ( !cell.Neighbor(i)->FE() || cell.Neighbor(i)->FE_Type() == UNKNOWN ) {
                 cell.UnassignNeighbor( i );
                 false_nbors++;
              }
      }
      
   return false_nbors;
   
 } // end

template size_t MeshManager<3>::RemoveDegenerateNeighbors<Element>();
template size_t MeshManager<2>::RemoveDegenerateNeighbors<Element>();
template size_t MeshManager<1>::RemoveDegenerateNeighbors<Element>();

template size_t MeshManager<3>::RemoveDegenerateNeighbors<Face>();
template size_t MeshManager<2>::RemoveDegenerateNeighbors<Face>();
template size_t MeshManager<1>::RemoveDegenerateNeighbors<Face>();

template size_t MeshManager<3>::RemoveDegenerateNeighbors<InterFace>();
template size_t MeshManager<2>::RemoveDegenerateNeighbors<InterFace>();
template size_t MeshManager<1>::RemoveDegenerateNeighbors<InterFace>();
*/


/**
     As method above but assuming that the Elements on either side of the interface are already detached, having separate nodes.
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
   if ( inner_parent == nullptr ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::AddInterFace", "pointer to higher dimensional element on inside not initialised");
        return nullptr;
     }
   if ( outer_parent == nullptr ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::AddInterFace", "pointer to higher dimensional element on ouside not initialised");
        return nullptr;
     }
   if ( inner_parent == outer_parent ) {
        csmp_error.Note( ERROR, "MeshManager<dim>::AddInterFace", "cannot create InterFace",
                                  "inner and outer parent pointers are the same");
        return nullptr;
     }

   assert( inner_element_face_id < inner_parent->Faces() );
   assert( outer_element_face_id < outer_parent->Faces() );

   // 1. establishing the finite element type if the interface
   const CSMP_FEM_TYPE etype = inner_parent->FE()->ElementTypeOfFace( inner_element_face_id );

   // 2. constructing new interface
   const size_t iface_id = interfaces_.size();
   const FiniteVolumeStencil<dim>* const stencil_ptr = (fvm_manager_) ? fvm_manager_->Stencil(etype) : nullptr;
   typename plf::colony<InterFace<dim>>::iterator
     ifp = interfaces_.emplace( InterFace<dim>( fem_manager_.E(etype), stencil_ptr, lvars, ivars ) );
     
   // 3. assigning higher dimensional elements and faces
   (*ifp).Assign( inner_parent, inner_element_face_id, outer_parent, outer_element_face_id );
   (*ifp).Idx( iface_id );

   // 4. assigning nodes
   (*ifp).InitialiseNodeVector();
   
   // 5. detaching the input Elements from one-another
   inner_parent->Unassign( outer_parent );
   outer_parent->Unassign( inner_parent );

   return &(*ifp);

} // end AddInterFace










/**
    Constructs a new Interface. The former nodes of the Face become those on the inside of the interface.
    Duplicates the nodes identified by the vector, adding them to the required manifolds.
    Other nodes will be shared across the sides of the interface.
    
    @attention assumes that nodes have already been duplicated as necessary and manifolds have been created and are uptodate
*/
template<uint32_t dim>
InterFace<dim>* const MeshManager<dim>::ReplaceFaceByInterFace( typename vector<Face<dim>*>::iterator fptr_it,
                                                                const LocalVariables& lvars,
                                                                const IntegrationPointVariables& ivars,
                                                                vector<Node<dim>*> outside_nodes )
{
   ErrorHandler&  csmp_error( ErrorHandler::Instance() );

   // 0. verifying the input
   // pointers
   if ( (*fptr_it) == nullptr )
     csmp_error.Note( ERROR, "MeshManager<dim>::ReplaceFaceByInterFace", "Face pointer not initialised");

   // 1. constructing new interface
   // -----------------------------
   // - constructor takes care of assigning the higher dimensional elements
   // - constructor assigns outside nodes also changing the nodes of the higher-dimensional element
   // - constructor detaches the neighbor connection between the higher dimensional elements
   const size_t iface_id = interfaces_.size();
   typename plf::colony<InterFace<dim>>::iterator
     ifp = interfaces_.emplace( InterFace<dim>( (*fptr_it), lvars, ivars, outside_nodes ) );
     
   // 2. assigning idx
   (*ifp).Idx( iface_id );

#ifdef DEBUG
   // verifying that the nodes on the inside matching those of the face
   for ( uint32_t i{0U}; i<(*fptr_it)->Nodes(); ++i ) {
        assert( (*ifp).N(i) != nullptr );
        assert( (*ifp).N(i,INSIDE) != nullptr );
        assert( (*fptr_it)->N(i) == (*ifp).N(i,INSIDE) );
        assert( (*ifp).N(i,OUTSIDE) != nullptr );
     }
#endif

   // 3. removing original face
   Delete( fptr_it );

   return &(*ifp);

} // end ReplaceFaceByInterFace







/**
    Duplicates existing node and connects it to corresponding manifold, else, the existing node is returned.
    
    @attention the current node is assumed to be on the INSIDE of the Interface; when there is no manifold yet.
    
    @param nptr_inside pointer to the node that will be on the inside of the InterFace that gets created if any.
    @return pointer to the new node now stored by the MeshManager.
*/
template<uint32_t dim>
Node<dim>* const MeshManager<dim>::Duplicate( Node<dim>* const nptr_inside,
                                              const LocalVariables& lvars )
  {
    if ( nptr_inside == nullptr )
      throw csmp::Exception( ERROR, "MeshManager<dim>::Duplicate", "Node pointer is a 'nullptr'.");

    // copying the inside node to create a new node
    auto nit = AddNodeAt( nptr_inside->Coordinate(), lvars, nptr_inside->AtBoundary() );
    
    // copying the properties over
    (*nit).CopyPropertyValuesFrom( *nptr_inside );

    // creating or updating the NodeManifold
    if ( nptr_inside->IsManifold() ) {
         // if we are already dealing with a manifold, the new node is added to it
         nptr_inside->Manifold()->Add( &(*nit) );
         // (*nit).Assign( (*nptr_inside->Manifold()) ); is already done by Add()
      }
    else {
         // checking that the NodeManifoldManager has been initialised
        if ( !node_manifold_manager_ )
           node_manifold_manager_ = new NodeManifoldManager<dim>();
           
         // a new manifold from the old and the new node using the provided default geometric classifier
         auto nmf = node_manifold_manager_->AddManifold( nodes_, nptr_inside, &(*nit), ManifoldType::SPLIT_BOUNDARY );
         // and its nodes are connected to it
         nptr_inside->Assign( (*nmf) );
         (*nit).Assign( (*nmf) );
      }

    // working out whether the original classification as an interface was correct
#ifdef DEBUG
    const bool verbose{false};
    consistencyCheck( (*(*nit).Manifold()), verbose );
#endif
    return &(*nit);
    
  } // end Duplicate
    



/**
    Rotates node numbers while keeping the separation of nodes into corner nodes, midside nodes, further nodes (the latter no being touched)

    Turning the nodes from counterclockwise looking from the outside in to clockwise.
*/
template<uint32_t dim>
void reverseOrderOfFaceNodes( uint32_t crn_nodes,
                              typename vector<Node<dim>*>::iterator first,
                              typename vector<Node<dim>*>::iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    const long n_nodes = distance(first,last);
    assert( crn_nodes >= n_nodes );
         
    // for linear elements which only have corner nodes
    if ( crn_nodes == n_nodes ) {
         reverse( first, last );
         return;
      }

    // for quadratic elements which only have corner nodes
    if ( (crn_nodes * 2) % n_nodes == 0 ) {
         reverse( first, next(first,crn_nodes) ); // reverse corner nodes
         reverse( next(first,crn_nodes), last ); // reverse midside nodes
         return;
      }
   
    csmp_error.Note( ERROR, "reverseOrderOfFaceNodes", "cases where n-nodes > 2 * crn_nodes not handled yet");

 } // end reverseNodeOrder





/**
   replaces supplied lower-dimensional elements with Face objects, establishing their connectivity; the Elements are deleted afterwards, setting input pointers to NULL
   
      the Node flags of the Element are used to determine whether this is a boundary face
      
     @note RANGE ERASE DOES ONLY WORK FOR A CONSECUTIVE RANGE OF ITERATORS WHERE it1 < it2
     @code
     elements_.erase( (*elmt_iterators.begin()), (*elmt_iterators.end()) );
     @endcode

*/
template<uint32_t dim>
vector<Face<dim>*>  MeshManager<dim>::ReplaceInteriorElementsByFaces( const PropertyDatabase<dim>& pref,
                                                                      typename vector<Element<dim>*>::iterator first,
                                                                      typename vector<Element<dim>*>::iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    vector<Face<dim>*> face_ptrs;
    const auto         n_faces_to_build{ distance(first,last) };
#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
    auto               n_elements = distance(first,last);
#endif
    if ( n_faces_to_build == 0U ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::ReplaceInteriorElementsByFaces", "supplied iterator range is empty; nothing was done.");
         return face_ptrs;
      }
    else face_ptrs.reserve( static_cast<size_t>(n_faces_to_build) );
    
    const LocalVariables             lvars(pref.LocalVariablesAt(FACE));
    const IntegrationPointVariables& ivars(pref.IntegrationPointVariablesAt(FACE));
    
    // 1. converting Elements into Faces
    // ---------------------------------
    auto first_copy = first;
    size_t face_idx{0};
    
    while( first_copy != last )
      {
         // 1.1 initial checks and labeling
         // (input range must not contain any nullptrs)
         assert( (*first_copy) != nullptr );
         if constexpr ( dim == 3U )
           if ( !(*first_copy)->IsSurface() ) {
                (*first_copy)->Out();
                throw csmp::Exception( ERROR, "MeshManager<3>::ReplaceInteriorElementsByFaces", "supplied element is not a surface element and cannot be converted to Face.");
             }
         if constexpr ( dim == 2U )
           if ( !(*first_copy)->IsLine() ) {
                (*first_copy)->Out();
                throw csmp::Exception( ERROR, "MeshManager<2>::ReplaceInteriorElementsByFaces", "supplied element is not a line element and cannot be converted to Face.");
             }
        
         // 1.2 construction of Face object in the interior of a model where both neighbors are present
         // -------------------------------------------------------------------------------------------
         // 1.2.1 finding higher dimensional elements to connect the Face to
         auto  pelmts = parentElements<dim>( (*first)->CornerNodesBegin(), (*first)->CornerNodesEnd() );
         assert( pelmts.first.first  != nullptr );
         assert( pelmts.second.first != nullptr );
         // 1.2.2 create Face, storing pointer to it; deletes lower-dimensional element from which it was constructed
         face_ptrs.push_back( ReplaceElementByFace( first_copy, pelmts.first.first, pelmts.second.first,
                                                    pelmts.first.second, pelmts.second.second, lvars, ivars, false ) );
         // numbering new Face consecutively
         face_ptrs.back()->Idx( face_idx++ );
       
         // making sure that the elements on the perimeter no longer connect to deleted elements
         DetachNeighborsFrom( (*first_copy) );
         // NOTE: no Element erasure yet because this would invalidate node parent vector, corrupting this functionality
         first++;
      }
      
     // 2. Deleting the replaced elements
     // ---------------------------------
     while( first != last )
      {
          Delete( first );
          first++;
      }
     // RANGE ERASE DOES ONLY WORK FOR A CONSECUTIVE RANGE OF ITERATORS WHERE it1 < it2
     //elements_.erase( (*elmt_iterators.begin()), (*elmt_iterators.end()) );
#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
     cout <<"\n\nMeshManager: ReplaceInteriorElementsByFaces: created "<< face_ptrs.size();
     cout <<" faces and deleted "<< n_elements <<" elements."<< endl;
#endif

     // 2. cleaning up inter-CELL and node to parent connectivity
     // ---------------------------------------------------------
     BuildConnectivity<csmp::Face>( face_ptrs.begin(), face_ptrs.end() ); // between the faces
     
     return face_ptrs;
     
 } // end ReplaceInteriorElementsByFaces







/**
     Replaces supplied lower-dimensional elements with Face objects, establishing their connectivity; the Elements are deleted afterwards, setting input pointers to NULL
   
      the Node flags of the Element are used to determine whether this is a boundary face.
      
      @attention If any of the lower-dimensional elements that are to be converted do not lie on the model boundary, they are not converted into faces,
      but still deleted. This avoids that faces are created away from the outer surface of the model, but it may leave gaps in the Face cover of the boundary.
      In that case the input regions need to be repaired before this method is applied.
      
     @note PLF::COLONY RANGE ERASE DOES ONLY WORK FOR A CONSECUTIVE RANGE OF ITERATORS WHERE it1 < it2
     @code
     elements_.erase( (*elmt_iterators.begin()), (*elmt_iterators.end()) );
     @endcode

*/
template<uint32_t dim>
vector<Face<dim>*>  MeshManager<dim>::ReplaceBoundaryElementsByFaces( const PropertyDatabase<dim>& pref,
                                                                      typename vector<Element<dim>*>::iterator first,
                                                                      typename vector<Element<dim>*>::iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    vector<Face<dim>*> face_ptrs;
    const auto         n_faces_to_build{ distance(first,last) };
#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
    auto               n_elements   = distance(first,last);
    size_t             n_duplicates = duplicatesCheck<dim,Element>( first, last );
    assert( n_duplicates == 0 );
#endif
    if ( n_faces_to_build == 0U ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::ReplaceBoundaryElementsByFaces", "supplied iterator range is empty; nothing was done.");
         return face_ptrs;
      }
    else face_ptrs.reserve( static_cast<size_t>(n_faces_to_build) );
    
    const LocalVariables             lvars(pref.LocalVariablesAt(FACE));
    const IntegrationPointVariables& ivars(pref.IntegrationPointVariablesAt(FACE));
    
    // 1. converting Elements into Faces
    // ---------------------------------
    size_t face_idx{0ul};
    // remembering the first iterator
    auto first_copy{ first };
    
    while( first_copy != last )
      {
         // 1.1 initial checks and labeling
         // -------------------------------
         // (input range must not contain any nullptrs)
         assert( (*first_copy) != nullptr );
         if constexpr ( dim == 3U )
           if ( !(*first_copy)->IsSurface() ) {
                (*first_copy)->Out();
                throw csmp::Exception( ERROR, "MeshManager<3>::ReplaceBoundaryElementsByFaces",
                                      "supplied element is not a surface element and cannot be converted to Face.");
             }
         if constexpr ( dim == 2U )
           if ( !(*first_copy)->IsLine() ) {
                (*first_copy)->Out();
                throw csmp::Exception( ERROR, "MeshManager<2>::ReplaceBoundaryElementsByFaces",
                                      "supplied element is not a line element and cannot be converted to Face.");
             }
        
         // 1.2 constructing Faces at model boundary, verifying that the supplied elemnent actually is located on the boundary
         // ------------------------------------------------------------------------------------------------------------------
         // trying to find higher dimensional neighbor of the face-element and its face idx
         pair<Element<dim>* const,uint32_t> pelmt = parentElement<dim>( (*first_copy)->CornerNodesBegin(), (*first_copy)->CornerNodesEnd() );
         
         // if this is unsuccesful, boundary face creation is stopped, but element will be deleted
         if ( pelmt.first == nullptr ) {
              cerr <<"\n\t"<< (*first_copy)->Idx();
              csmp_error.Note( WARNING, "MeshManager<dim>::ReplaceBoundaryElementsByFaces", "unable to identify higher-dim element neighbor of Element.");
#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
              n_elements--;
#endif
              first_copy++;
              continue;
           }
           
         // creating Face, storing a pointer to it
         //                                    element ptr  local element ID in face
         face_ptrs.push_back( AddBoundaryFace( pelmt.first, pelmt.second, lvars, ivars ) );
         //                   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
         face_ptrs.back()->Idx( face_idx++ );
         
         // verifying the creation of th face and indexing it
         if ( face_ptrs.back() == nullptr ) {
              cout <<"\n\t"<<"Element<"<< dim <<">"<< pelmt.first->Idx() <<" "<< parseFiniteElementType( pelmt.first->FE_Type() );
              throw csmp::Exception( ERROR, "MeshManager<dim>::ReplaceBoundaryElementsByFaces",
                                    "could not create suitable face matching element");
           }
          // NOTE: no Element erasure yet because this would invalidate node parent vectors, corrupting this functionality
         first_copy++;
      }

     // 2. Deleting the Elements
     // ------------------------------------------------------------------------------------
     bool first_call{true};
     
     while( first != last )
       {
#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
          auto pfl_iter  = elements_.get_iterator(*first);
#endif
          auto next_elmt = Delete( first );
          //               ^^^^^^^^^^^^^^^
#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
          if ( (*first) != nullptr and
               next_elmt != elements_.end() and
               pfl_iter  != elements_.end() and
               next_elmt != next(pfl_iter,1) ) {
               if ( first_call ) {
                   csmp_error.Note( WARNING, "MeshManager<dim>::ReplaceBoundaryElementsByFaces", "did not delete Element");
                   cout <<"\n\t"<<"Hexadecimal addresses of Element(s) not deleted: "<< (*first);
                   first_call = false;
                 }
               else cout <<" "<< (*first);
            }
#endif
          first++;
        }
      if ( !first_call ) cout << endl;

#if defined(DEBUG) && defined(CSMP_MESH_MANAGER_DEBUG)
     cout <<"\n\nMeshManager: ReplaceBoundaryElementsByFaces: created "<< face_ptrs.size();
     cout <<" faces and deleted "<< n_elements <<" elements."<< endl;
     size_t surface_elmts{0};
     for ( const auto& it : elements_ ) if ( it.IsSurface() ) surface_elmts++;
     cout <<"\t\t"<<"there are "<< surface_elmts <<" surface elements left in the model."<< endl;
#endif

     // 3. cleaning up inter-CELL and node-to-parent connectivity
     // ---------------------------------------------------------
     BuildConnectivity<csmp::Face>( face_ptrs.begin(), face_ptrs.end() ); // between the faces

     return face_ptrs;
     
 } // end ReplaceBoundaryElementsByFaces








/**
   Replaces supplied dim-1 Element objects with InterFace objects, adding the necessary multiplicated nodes, identifying outside elements and flagging regions
   for rebuild
   and establishing their connectivity. Since the Face objects are deleted the supplied pointer ranges to them (interior and perimeter) are invalidated
   by this method.

   @param first iterator to the first FaceConstructionData of the supplied dim-1 region

   @param last iterator to the last FaceConstructionData of the dim-1 region which also points behind the last perimeter face

   @param perim_first iterator to perimeter node vector of dim-1 region

   @param perim_last iterator to perimeter node vector of dim-1 region

   @param region_material_ids   for material identifiers placed on Elements, Regions etc.
   
   A set which will be overwritten by method with region_identifier ids which should be scheduled for rebuilding since the elements (outside) have new nodes

   @assumption Assumes element property "region identifier" is defined.

   @attention: This method sets TOPO flags PERIMETER_LINE and PERIMETER_POINT


   Steps - Creation of interfaces from FaceConstructionData:

   1. Loops over FaceConstructionData and duplicates nodes of InnerParent if they are found to Not be within the Perimeter Nodes supplied.
      If nodes are at a perimeter, no duplication occurs. If nodes are already a manifold (intersection), then duplication occurs (even at perimeter).

   2. INTERNAL flags are added to nodes on both INSIDE and OUTSIDE unless at a model boundary.

   3. Outer Parent is assigned the new nodes created

   4. Interface objects are constructed based on Former higher-dimensional Element neighbors that are now connected to different nodes, but which share
      the same coordinates. Interface construction uses the coordinates to distinguish matching nodes. Higher dim elements unnassign each other as neighbors in the process.

   5. A neighbor search algorithm is performed on the outside of the interface objects which replaces the old inside node with the new duplicate node for all parent elements
      (perimeter nodes must be well defined for this to exclude the inside nodes from the assignment)
      All the outside elements are queried for their unique Region they belong to using "region identifier" property, this is stored for SplitBoundaryInterface

   6. If a node on an existing interface was split (even if on perimeter), then the interface object is updated to contain (potentially) the duplicate node.

   7. Finally the node-to-parent element connectivity for the whole mesh is rebuilt.
   8. Likewise, the node-to-parent interface connectivity is also rebuilt, since all interfaces (new and old), now have correct nodes within them.

   Cases:
   3D: An X intersection, where the perimeter of split boundary 1 (not split), should be split in the dimension of splitboundary 2.
   Therefore, splitboundary node on INSIDE,OUTSIDE of SB1 is not split, and not manifold.
   However, same node, is classified a manifold by SB2, and is split, with different nodes on INSIDE OUTSIDE of SB2. How to disambiguate?

   If split if manifold --> error when we split SB1.


   @author E.P
   @date 18/8/22

   // TODO: Take away dependency on PerimeterNodes iterators when the TOPO flags can be relied upon

*/
template<uint32_t dim>
vector<InterFace<dim>*>  MeshManager<dim>::ReplaceElementsByInterFaces( const PropertyDatabase<dim>& dbase,
                                                                        typename vector<FaceConstructionData<dim>>::iterator first,
                                                                        typename vector<FaceConstructionData<dim>>::iterator last,
                                                                        typename vector<Node<dim>*>::const_iterator perim_first,
                                                                        typename vector<Node<dim>*>::const_iterator perim_last,
                                                                        set<Node<dim>*>& split_perimeter_nodes,
                                                                        set<size_t>& region_material_ids,
                                                                        bool convert_original_elements_to_intervening_elements )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    // vector of interfaces which will be returned
    vector<InterFace<dim>*>  iface_ptrs;
    const size_t  n_original_elmts{ static_cast<size_t>(distance(first,last))};

    const csmp::Index region_key = dbase.StorageKey( "region identifier");

    if (n_original_elmts == 0U ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::ReplaceElementsByInterFaces", "supplied iterator range is empty; nothing was done.");
         return iface_ptrs;
      }
    else iface_ptrs.reserve( n_original_elmts );

    // Constructing the node manifold manager if necessary
    const bool no_previous_manifolds = ( node_manifold_manager_ == nullptr ) ? true : false;
    if ( no_previous_manifolds )
      node_manifold_manager_ = new NodeManifoldManager<dim>();

    // establish the storage requirements for face variables
    const LocalVariables             lvsNode( dbase.LocalVariablesAt(NODE) );
    const LocalVariables             lvsInterfaces(dbase.LocalVariablesAt(INTER_FACE) );
    const IntegrationPointVariables  lvsIntegrationPoints( dbase.IntegrationPointVariablesAt(INTER_FACE) );

    set<Element<dim>*>    outside_neighbors_to_search;
    set<Element<dim>*>    inside_parents;
    // tracking already duplicated nodes to avoid duplicates
    //  original,  duplicate
    map<Node<dim>*,Node<dim>*> new_nodes;
    map<Node<dim>*,Node<dim>*> new_nodes_with_existing_manifold;
    map<Node<dim>*,Node<dim>*> new_perimeter_nodes_without_manifold;


    // 1. converting interior Face objects into InterFace ones, duplicating their nodes
    // --------------------------------------------------------------------------------
    typename vector<FaceConstructionData<dim>>::iterator first_copy = first;
    
    // Note: WE do NOT split nodes on the perimeter of this subdomain
    while( first_copy != last )
      {
         // collecting neighbors of outer parent that dont include inner parent
         //Getting Inside Element
         Element<dim>* inside_elmt    = first_copy->InnerElement();
         Element<dim>* outside_elmt   = first_copy->OuterElement();
         Element<dim>* lower_dim_elmt = first_copy->LowerDimElement();   //used to attribute TOPO Flag

         outside_neighbors_to_search.insert( outside_elmt );
         inside_parents.insert( inside_elmt );

         // 1.2 duplicating the nodes creating manifolds as necessary and adding boundary flags
         // -----------------------------------------------------------------------------------
         //getting vector of inside nodes to iterater
         map<Node<dim>*, Node<dim>*> in_out_nodes;
         auto            nit{ new_nodes.end() };
         // creating the node vector and reverting its order so that it matches the face of the higher dimensional outside element
         for ( const auto& i : inside_elmt->FE()->NodesOfFace( first_copy->InnerElementFace() ) ) {
           Node<dim>* inside_node = inside_elmt->N(i);
           //If we are at not at perimeter, or if we are at intersection (manifold)
           bool inside_was_manifold  = inside_node->IsManifold();
           bool inside_was_perimeter = (inside_node->Attribute() == PERIMETER_LINE ||
                                        inside_node->Attribute() == PERIMETER_POINT );  //if we hit another regions perimeter (not the region we are splitting)
           if ( find( perim_first, perim_last, inside_node ) == perim_last && split_perimeter_nodes.find(inside_node) == split_perimeter_nodes.end() ){
             if ( (nit=new_nodes.find( inside_node )) == new_nodes.end() ) {                // if a matching outside node has not been created yet
                  //duplicating outside node if not already duplicated
                  Node<dim>* out_node = Duplicate( inside_node, lvsNode );
                  in_out_nodes.insert( make_pair( inside_node, out_node ));
                  new_nodes.insert( make_pair( inside_node, out_node ) );
                  if (inside_was_manifold){
                    assert(!inside_was_perimeter);
                    new_nodes_with_existing_manifold.insert(make_pair(inside_node, out_node));
                  }
                  if (inside_was_perimeter){
                    assert(!inside_was_manifold);
                    new_perimeter_nodes_without_manifold.insert(make_pair(inside_node,out_node));
                  }
               }
             // if the necessary new node was already created earlier it was retrieved and is assigned here
             else in_out_nodes.insert( make_pair( inside_node, (*nit).second ));
           } else {
             //we are on the perimeter
             in_out_nodes.insert(make_pair( inside_node, inside_node ));            //take inside node when node is at perimeter of model
             //assign perimeter topo flag (needed if we have to split a perimeter later)
             if (lower_dim_elmt->IsSurface()) inside_node->Attribute( PERIMETER_LINE  );
             if (lower_dim_elmt->IsLine() )   inside_node->Attribute( PERIMETER_POINT );

           }

           // add Inside and outside node INTERNAL flag if not at boundary
           if ( inside_node->AtBoundary() == NOT ) {
             inside_node->AtBoundary(INTERNAL);
             in_out_nodes[inside_node]->AtBoundary(INTERNAL);
           }

         }//end of node loop and duplication

         // assign outside nodes to OuterParent element
         for ( const uint32_t& n : outside_elmt->FE()->NodesOfFace( first_copy->OuterElementFace() ) ) {
           // TODO: test with a meaningful number that takes model size into account rather than this is arbitrary tolerance!
           assert( distance( outside_elmt->N(n)->Coordinate(), in_out_nodes[outside_elmt->N(n)]->Coordinate()) < 0.001 );
           outside_elmt->Assign(n, in_out_nodes[ outside_elmt->N(n) ] );
         }

         // 1.4 construction of InterFace from parent elements
         // --------------------------------------------------
         //     - higher-dimensional nbors are already known (INSIDE , OUTSIDE)
         //     - faces of higher dimensional neighbors are also known
         //     - nodes on outside are known and ASSIGNED to the OuterParent (old nodes are on inside, new nodes on outside)
         //     - nodes on inside and outside are the same for perimeter interfaces away from boundaries
         iface_ptrs.push_back( WrapInterFaceAroundElement( first_copy->LowerDimElement(),
                                                           first_copy->InnerElement(),
                                                           first_copy->OuterElement(),
                                                           first_copy->InnerElementFace(),
                                                           first_copy->OuterElementFace(),
                                                           lvsInterfaces, lvsIntegrationPoints ) );
                                                           
         // STILL NEEDED? - DetachNeighborsFrom( (*first_copy).LowerDimElement() );
         first_copy++;
      }


    // 2. Assigning new nodes to outside elements - Neighbor search being performed
    // ----------------------------------------------------------------------------
     // Loop over outer parents and search neighbors with old node
     while ( outside_neighbors_to_search.empty() == false ){
       typename set<Element<dim>*>::iterator eit = outside_neighbors_to_search.begin(); //take start of set

       //search for neighbor
       for (uint32_t nbor{0U}; nbor< (*eit)->Neighbors(); nbor++){
         Element<dim>* e_nbr = (*eit)->Neighbor(nbor);
         if ( e_nbr != nullptr ){                                              //if neighbor exists
           //search the nodes of neighbor for inside node
           uint32_t n_nodes = e_nbr->Nodes();
           for (uint32_t n{0U}; n < n_nodes; n++ ){
             typename map<Node<dim>*,Node<dim>*>::iterator found_it = new_nodes.find( e_nbr->N(n) ); //searching for inside Node in element
             if ( found_it != new_nodes.end() ){
               e_nbr->Assign(n, found_it->second ); //replacing inside node of neighbor with outside node
               assert( inside_parents.find(e_nbr) == inside_parents.end() );
               outside_neighbors_to_search.insert( e_nbr ); //add to search, so that neighbors of this neighbor are searched

               //Extract outside neighbor and his Reg ID (which is defined on the elements within SplitBoundaryInterface
               ScalarVariable reg_id;
               e_nbr->Read( region_key, reg_id );
               region_material_ids.insert( static_cast<uint32_t>(reg_id()) );
             }
           } //looped over all nodes
         } //valid neighbor
       }//end of neighbor search - all relevant neighbors have been added to the future neighbor search - and node numbers updated

       // Current element has completed its neighbor search and node assignment
       outside_neighbors_to_search.erase( eit );

     }//search continues until all inside nodes are updated


     // 3.0 Update existing Interfaces which were intersected
     // ---------------------------------------------------------
     // 3.1 For each duplicated node that was already a manifold (intersection)
     for ( auto& nit : new_nodes_with_existing_manifold ) {
       Node<dim>* old_node = nit.first;
       const auto interfaces = old_node->Manifold()->InterFaces(old_node);
       //iterate over all interfaces of the node
       uint32_t found{0U};
       for ( uint32_t i{0U}; i < interfaces; i++ ){
         //getting interface information
         pair<InterFace<dim>*, pair<uint32_t, INTERFACE_SIDE>> interface_index = old_node->Manifold()->InterFaceIndex(old_node, i);
         INTERFACE_SIDE side = interface_index.second.second;
         //getting nodes face of higher-dim parent element
         set<Node<dim>*> nds_of_parent = interface_index.first->Parent(side)->CornerNodesOfFace(interface_index.first->ParentFaceID(side));
         // if Parent was updated with new node
         if (nds_of_parent.find( old_node ) == nds_of_parent.end() ){
           assert(nds_of_parent.find(nit.second) != nds_of_parent.end()); //check the parent element has the new duplicate node
           //update interface with new node:         index             , new node,           side
           interface_index.first->Assign(interface_index.second.first, nit.second, interface_index.second.second );
           found++;
         }
       }
       //should only update one interface per node
       if constexpr (dim==2) assert(found == 1); //for a 1d splitboundary in 2d model
     }


     //3.2 For each duplicated node that was on a splitboundaries perimeter (which is not a manifold)
     for (auto& nit : new_perimeter_nodes_without_manifold){
       Node<dim>* perimeter_node = nit.first;

       //i) search for potential interfaces
       set<InterFace<dim>*> potential_interfaces;
       const uint32_t neighbors = perimeter_node->Neighbors();
       //search for manifold nodes that are neighbors, and grab their interfaces
       for (uint32_t n{0U}; n<neighbors; ++n){
         Node<dim>* neighbor_node = perimeter_node->Neighbor(n);
         //We need a manifold object that has been configured (that hasnt just been created above)
         if (neighbor_node->IsManifold() && neighbor_node->Manifold()->NodeMapSize() != 0 ){
           const uint32_t interfaces = neighbor_node->Manifold()->InterFaces(neighbor_node);
           assert(interfaces > 0);
           for ( uint32_t i{0U}; i<interfaces; ++i ){
             potential_interfaces.insert(neighbor_node->Manifold()->I(neighbor_node, i)); //inserting interface as potentially having perimeter node
           }
         }//found potential interface
       }//end interface search
       assert(!potential_interfaces.empty()); //This can be the case if all nodes on interface are on the perimeter!! (I.e they are not split)

       //ii) Search all interfaces for perimeter node
//       bool found_perimter_interface = false;
       for ( auto& ifit : potential_interfaces){
         //search for perimeter node (inside outside nodes match)
         for ( uint32_t n{0U}; n<ifit->FE()->Nodes(); ++n ){
            if (ifit->N(n,INSIDE) == perimeter_node && ifit->MatchingN(n,INSIDE) == ifit->MatchingN(n,OUTSIDE)){
              //We should also have perimeter node on the outside
              assert(ifit->MatchingN(n,OUTSIDE) == perimeter_node );
              //Update both sides of interface with new node IF higher dim parents say so!
              set<Node<dim>*> nds_of_in_face  = ifit->InnerParent()->CornerNodesOfFace( ifit->InnerParentFaceID() );
              set<Node<dim>*> nds_of_out_face = ifit->OuterParent()->CornerNodesOfFace( ifit->OuterParentFaceID() );
              if (nds_of_in_face.find(perimeter_node) == nds_of_in_face.end()){
                assert(nds_of_out_face.find(perimeter_node) == nds_of_out_face.end());  //Check also outside face doesnt have perimeter node (not true perimeter)
                assert(nds_of_in_face.find(nit.second) != nds_of_in_face.end());        //check we have the duplicate node instead
                assert(nds_of_out_face.find(nit.second) != nds_of_out_face.end());      //check we have the duplicate node instead

                //assign interface the new node instead
                ifit->Assign(n, nit.second, INSIDE);
                assert(n < ifit->FE()->CornerNodes()); //We must have a corner node
                assert( ifit->MatchingN(n,OUTSIDE) == ifit->N(ifit->FE()->CornerNodes() - 1 - n, OUTSIDE)); //Only applied if corner node, and assumes outside nodes rotate the other way and start from last node
                ifit->Assign(ifit->FE()->CornerNodes() - 1 - n , nit.second, OUTSIDE );
              }// end of interface correction
            }//end of perimeter node found
         }//end of perimeter node search
       }//end of interfaces corrections
     }//end of perimeter nodes fixes

    
    // 4. Deleting intervening elements if requested
    // ---------------------------------------------
    if ( convert_original_elements_to_intervening_elements ) {
        // 4.1 deleting the intervening elements
        while( first != last ) {
            // geting iterator to element ('get_iterator' returns end, if element cannot be found)
            auto pfl_it = elements_.get_iterator( (*first).LowerDimElement() );
            // if the element exists in the colony
            if ( pfl_it != elements_.end() ) {
                 // remove element from the parent element list of its connected nodes
                 for ( auto nit=(*first).LowerDimElement()->NodesBegin(); nit!=(*first).LowerDimElement()->NodesEnd(); ++nit )
                   (*nit)->Unassign( (*first).LowerDimElement() );
                 // here one should null the element pointer, but nothing can be done to FaceConstructionData
                 //(*first) = nullptr;
                 // deleting element and returning colony iterator to next element in colony
                 elements_.erase( pfl_it );
              }
            first++;
          }
        // 4.2 setting the intervening element pointers in the InterFace objects to zero
        for ( auto& interface_ptr : iface_ptrs )
          interface_ptr->Assign( nullptr );
      }


    BuildConnectivity<csmp::InterFace>( iface_ptrs.begin(), iface_ptrs.end() ); // between the interfaces

    cout <<"\n"<<"MeshManager<"<< dim <<">::ReplaceElementsByInterFaces: created "<< iface_ptrs.size() <<" new interfaces and ";
    cout << new_nodes.size() <<" new nodes."<< endl;

    return iface_ptrs;

 } // end ReplaceFacesByInterFaces







/**
   Replaces supplied Face objects with InterFace objects, adding the necessary multiplicated nodes
   and establishing their connectivity. Since the Face objects are deleted the supplied pointer ranges to them (interior and perimeter) are invalidated
   by this method.
   
   @param dbase is needed for the inialisation of the LocalVariableStorage associated with the Face objects
   
   @param first iterator to the first Face of the supplied boundary
   
   @param bfirst  is an iterator that simultaneously is the end of the interior faces and the beginning of the boundary faces
   
   @param last iterator to the last Face of the boundary which also points behind the last perimeter face
   
   Steps - starting with the processing of interior faces:
   
   1. The creation of InterFace objects in the interior and the necessary duplication of nodes occur simultaneously.
   
   2. Interface objects are constructed from the second range, also duplicating nodes at the model- or domain boundaries. No new nodes are inserted if the boundary terminates inside of a higher dimensional region. In this case, interior and exterior perimeter nodes on the SplitBoundary perimeter are assigned the same, pre-existing perimeter node inherited from the converted boundary.

   During the construction of all InterFace objects, the higher dimensional neighbors are assigned and the corresponding Element faces are remembered.

   3. Former higher-dimensional Element neighbors that are now separated by the SplitBoundary are disconnected from one another assigning their neighbor  pointers to 'nullptr'.
     (this step is not necessary if the whole connectivty is rebuilt anyway)

   4. The new InterFace objects are connected with one-another so that a SplitBoundary constructor has all the necessary information to distinguish interior from perimeter.
   
   5. Finally the node-to-parent element connectivity of the InterFace nodes needs to be rebuilt restricting parent element access to the side of the interface that the node forms part of
   
   @author SKM
   @date 6/4/22


   //TODO: Take away dependency on PerimeterNodes iterators when the TOPO flags can be relied upon!!

*/
template<uint32_t dim>
vector<InterFace<dim>*>  MeshManager<dim>::ReplaceFacesByInterFaces( const PropertyDatabase<dim>& dbase,
                                                                     typename vector<Face<dim>*>::iterator first,
                                                                     typename vector<Face<dim>*>::iterator bfirst,
                                                                     typename vector<Face<dim>*>::iterator last,
                                                                     typename vector<Node<dim>*>::const_iterator perim_first,
                                                                     typename vector<Node<dim>*>::const_iterator perim_last  )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    vector<InterFace<dim>*>  iface_ptrs;
    const size_t  n_original_faces{ faces_.size() };

    if ( distance(first,last) == 0U ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::ReplaceFacesByInterFaces", "supplied iterator range is empty; nothing was done.");
         return iface_ptrs;
      }
    else iface_ptrs.reserve( n_original_faces );
    
    // Constructing the node manifold manager if necessary
    const bool no_previous_manifolds = ( node_manifold_manager_ == nullptr ) ? true : false;
    if ( no_previous_manifolds )
      node_manifold_manager_ = new NodeManifoldManager<dim>();
    
    const LocalVariables             nvars(dbase.LocalVariablesAt(NODE));
    const LocalVariables             lvars(dbase.LocalVariablesAt(INTER_FACE));
    const IntegrationPointVariables& ivars(dbase.IntegrationPointVariablesAt(INTER_FACE));
    
    // vector of interfaces which will be returned
    vector<InterFace<dim>*>  interface_ptrs;
    set<Element<dim>*>    outside_neighbors_to_search;
    set<Element<dim>*>    inside_parents;
    interface_ptrs.reserve( static_cast<size_t>(distance(first,last)) );
    // tracking already duplicated nodes to avoid duplicates
    //  original,  duplicate
    map<Node<dim>*,Node<dim>*>  new_nodes;

    // 1. converting interior Face objects into InterFace ones, duplicating their nodes
    // --------------------------------------------------------------------------------
    // (original Face objects are removed)
    while( first != bfirst )
      {
         // 1.1 initial checks
         // (input range must not contain any nullptrs)
         assert( (*first) != nullptr );
         
         // collecting neighbors of outer parent that dont include inner parent
         outside_neighbors_to_search.insert( (*first)->OuterParent() );
         inside_parents.insert( (*first)->InnerParent() );

         // 1.2 duplicating the nodes creating manifolds as necessary
         // ---------------------------------------------------------
         const auto         n_nodes{ (*first)->Nodes() };
         vector<Node<dim>*> outside_nodes( n_nodes, nullptr );
         auto               nit{ new_nodes.end() };
         // creating the node vector and reverting its order so that it matches the face of the higher dimensional outside element
         for ( uint32_t i{0U}; i<n_nodes; i++ ) {
           if ( find( perim_first, perim_last, (*first)->N(i) ) == perim_last ){             //using the perimeter nodes (Costly find operation)
           //if ( (*first)->N(i)->Attribute() != PERIMETER_POINT && (*first)->N(i)->Attribute() != PERIMETER_LINE ){        //WHEN WE CAN RELY ON TOPO FLAGS
             if ( (nit=new_nodes.find((*first)->N(i))) == new_nodes.end() ) {                // if a matching outside node has not been created yet
                  outside_nodes[i] = Duplicate( (*first)->N(i), nvars );
                  new_nodes.insert( make_pair( (*first)->N(i), outside_nodes[i] ) );
               }
             // if the necessary new node was already created earlier it was retrieved and is assigned here
             else outside_nodes[i] = (*nit).second;
           } else
             outside_nodes[i] = (*first)->N(i);            //take inside node when node is at perimeter of model
         }
         // reversing the sequence once outside nodes are calibrated
         reverseOrderOfFaceNodes<dim>( (*first)->FE()->CornerNodes(),
                                       outside_nodes.begin(), outside_nodes.end() );
  
           
         // 1.3 construction of InterFace away from boundaries
         // --------------------------------------------------
         //     - higher-dimensional nbors are already known
         //     - faces of higher dimensional neighbors are also known
         //     - nodes on outside are not known (old nodes are on inside, new nodes on outside)
         //     - nodes on inside and outside are the same for perimeter interfaces away from boundaries
         interface_ptrs.push_back( ReplaceFaceByInterFace( first, lvars, ivars, outside_nodes ) );
         first++;
      }


    // 2. Converting perimeter faces into interfaces, dealing with boundaries
    // ----------------------------------------------------------------------
    // (original Face objects are removed)
    first = bfirst;
    while( first != last )
      {
         // 2.1 initial checks
         // (input range must not contain any nullptrs)
         assert( (*first) != nullptr );

         // 2.2 duplicating nodes but only if we are at a model boundary or the node already is a manifold
         // ----------------------------------------------------------------------------------------------
         // (if the node is not duplicated, the original node is inserted into the InterFace outside node vector)
         const auto         n_nodes{ (*first)->Nodes() };
         vector<Node<dim>*> outside_nodes( n_nodes, nullptr );
         auto               nit{ new_nodes.end() };
         for ( uint32_t i{0U}; i<n_nodes; i++ )
           if ( ((*first)->N(i)->AtBoundary() != NOT && (*first)->N(i)->AtBoundary() != INTERNAL) || (*first)->N(i)->IsManifold() ) {
                if ( (nit=new_nodes.find((*first)->N(i))) == new_nodes.end() ) {
                     outside_nodes[i] = Duplicate( (*first)->N(i), nvars );
                     new_nodes.insert( make_pair( (*first)->N(i), outside_nodes[i] ) );
                  }
                else outside_nodes[i] = (*nit).second;
             }
           else outside_nodes[i] = (*first)->N(i);

         //reverse( outside_nodes.begin(), outside_nodes.end() );
         reverseOrderOfFaceNodes<dim>( (*first)->FE()->CornerNodes(),
                                       outside_nodes.begin(), outside_nodes.end() );

         // 2.3 construction of InterFace on the model perimeter
         // ----------------------------------------------------
         //     - higher-dimensional nbors are already known
         //     - faces of higher dimensional neighbors are also known
         //     - nodes on outside are not known (old nodes are on inside, new nodes on outside)
         //     - nodes on inside and outside are the same for perimeter interfaces away from boundaries
         //     - boundaries are inferred, when:
         //       - BOX_BOUNDARY flag is !NOT
         interface_ptrs.push_back( ReplaceFaceByInterFace( first, lvars, ivars, outside_nodes ) );
         first++;
      }
 

     // 3. cleaning up inter-CELL and node to parent connectivity
     // ---------------------------------------------------------
      // Reassigning outside parent elements with new node
      // 3.1 Loop over outer parents and search neighbors with old node
      while ( outside_neighbors_to_search.empty() == false ){
        typename set<Element<dim>*>::iterator eit = outside_neighbors_to_search.begin(); //take start of set

        //search for neighbor
        for (uint32_t nbor{0U}; nbor< (*eit)->Neighbors(); nbor++){
          Element<dim>* e_nbr = (*eit)->Neighbor(nbor);
          if ( e_nbr != nullptr ){                                              //if neighbor exists
            //search the nodes of neighbor for inside node
            uint32_t n_nodes = e_nbr->Nodes();
            for (uint32_t n{0U}; n < n_nodes; n++ ){
              typename map<Node<dim>*,Node<dim>*>::iterator found_it = new_nodes.find( e_nbr->N(n) ); //searching for inside Node in element
              if ( found_it != new_nodes.end() ){
                e_nbr->Assign(n, found_it->second ); //replacing inside node of neighbor with outside node
                assert( inside_parents.find(e_nbr) == inside_parents.end() );
                outside_neighbors_to_search.insert( e_nbr ); //add to search, so that neighbors of this neighbor are searched
              }
            } //looped over all nodes
          } //valid neighbor
        }//end of neighbor search - all relevant neighbors have been added to the future neighbor search - and node numbers updated

        // Current element has completed its neighbor search and node assignment
        outside_neighbors_to_search.erase( eit );

      }//search continues untill all inside nodes are updated

     BuildConnectivity<csmp::InterFace>( interface_ptrs.begin(), interface_ptrs.end() );
     
     cout <<"\n"<<"MeshManager<"<< dim <<">::ReplaceFacesByInterFaces: created "<< interface_ptrs.size() <<" new interfaces and ";
     cout << new_nodes.size() <<" new nodes."<< endl;
     
     return interface_ptrs;
     
 } // end ReplaceFacesByInterFaces






template<uint32_t dim>
vector<Face<dim>*>  MeshManager<dim>::CreateFacesBetweenNodeSharingElements( const PropertyDatabase<dim>& dbase,
                                                                             const vector<pair<pair<Element<dim>*,uint32_t>,
                                                                                               pair<Element<dim>*,uint32_t> > >& face_nbor_elmts )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( face_nbor_elmts.empty() ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::CreateFacesBetweenNodeSharingElements", "supplied range of element pairs is empty; nothing was done.");
         return vector<Face<dim>*>{}; // empty vec
      }
    
    // Property storage
    const LocalVariables             lvars(dbase.LocalVariablesAt(FACE));
    const IntegrationPointVariables& ivars(dbase.IntegrationPointVariablesAt(FACE));

    // 1. Creating Face objects
    // ------------------------
    vector<Face<dim>*> face_ptrs;
    face_ptrs.reserve( face_nbor_elmts.size() );
    
    for ( const auto& it : face_nbor_elmts )
      {
         // construction of Face
         // --------------------
         //     - higher-dimensional nbors are already known
         //     - faces of higher dimensional neighbors are also known
         //     - nodes on outside are not known (old nodes are on inside, new nodes on outside)
         face_ptrs.push_back( AddFace( it.first.first, it.first.second,
                                       it.second.first, it.second.second,
                                       lvars, ivars ) );
      }

     // 2. cleaning up inter-CELL and node to parent connectivity
     // ---------------------------------------------------------
    BuildConnectivity<csmp::Face>( face_ptrs.begin(), face_ptrs.end() ); // between the faces
     
     cout <<"\n"<<"MeshManager<"<< dim <<">::CreateFacesBetweenNodeSharingElements: created "<< face_ptrs.size() <<" new faces."<< endl;
     
     return face_ptrs;
     
  } // end CreateInterFacesBetweenNodeSharingElements






/**
    Creates InterFace objects between face and node sharing Elements adding the necessary extra nodes on the outside,
    and node manifolds where multiple nodes end up collocated.
    Finally the method creates the neighbor connectivity among the interfaces so that they can be organised into interior and perimeter.
    Then the connectivity of the overall mesh is updated, taking into account that elements and nodes loose their connections
    across the new interface.

    @param dbase a reference to property database needed for the initialisation of the LocalVariableStorage of potential interface variables
    @param interface_nbor_elmts vector of higher-dimensional elements that share a face where the new interface will be created.
    @param multiplicate_perimeter_nodes whether the nodes on the perimeter should be copied so that they are different in adjacent regions
    @param outside_elmts output of a set of element pointers to the volumetric elements neighboring the faces of interest
    
    The iinside elements are first in Element-element face pairs.
    
    @return returns vector of pointers to the newly created InterFace objects.
 */
template<uint32_t dim>
vector<InterFace<dim>*>  MeshManager<dim>::CreateInterfacesBetweenNodeSharingElements( const PropertyDatabase<dim>& dbase,
                                                                                       const vector<pair<pair<Element<dim>*,uint32_t>,
                                                                                                         pair<Element<dim>*,uint32_t> > >& interface_nbor_elmts,
                                                                                       bool multiplicate_perimeter_nodes,
                                                                                       unordered_set<Element<dim>*>& outside_elmts )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( interface_nbor_elmts.empty() ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::CreateInterfacesBetweenNodeSharingElements",
                                   "supplied range of element pairs is empty; nothing was done.");
         return vector<InterFace<dim>*>{}; // empty vec
      }
    
    // Constructing the node manifold manager if necessary
    const bool no_previous_manifolds = ( node_manifold_manager_ == nullptr ) ? true : false;
    if ( no_previous_manifolds )
      node_manifold_manager_ = new NodeManifoldManager<dim>();
    
    const LocalVariables             nvars(dbase.LocalVariablesAt(NODE));
    const LocalVariables             lvars(dbase.LocalVariablesAt(INTER_FACE));
    const IntegrationPointVariables& ivars(dbase.IntegrationPointVariablesAt(INTER_FACE));
    vector<Node<dim>*>               perimeter_node_ptrs;
    
    if ( multiplicate_perimeter_nodes == false )
      {
        // 1. finding the perimeter nodes of the interface patch that will be created
        // --------------------------------------------------------------------------
        perimeter_node_ptrs.reserve( interface_nbor_elmts.size() * dim ); // just a guess
        
        for ( const auto& it : interface_nbor_elmts ) {
             for ( const auto& i : it.first.first->FE()->NodesOfFace( it.first.second ) )
               perimeter_node_ptrs.push_back( it.first.first->N(i) );
          }
          
        // sorting and removing duplicates from node vector, making it searchable
        sort( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end() );
        perimeter_node_ptrs.erase( unique( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end() ), perimeter_node_ptrs.end() );
        // printNodeCoordinates<dim>( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end() );
        // NB: tested: at this point we only have perimeter nodes left


        // 2. Finding the subset of these nodes which will also be on the perimeter of the SplitBoundary
        // ---------------------------------------------------------------------------------------------
        // IMPORTANT: this needs to be done before creating the InterFace objects because it influences which of its nodes will have to be manifolds
        // How? - in 2D, these are the nodes that are only contained in one of the faces
        if constexpr (dim == 2U ) {
             vector<size_t> face_count( perimeter_node_ptrs.size(), 0U );
             // again
             for ( const auto& it : interface_nbor_elmts ) {
                  for ( const auto& i : it.first.first->FE()->NodesOfFace( it.first.second ) ) {
                       // finding the vector index corresponding to the perimeter node
                       auto lb = perimeter_node_ptrs.end();
                       if ( (lb=lower_bound( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end(), it.first.first->N(i) )) !=  perimeter_node_ptrs.end() )
                       face_count[ static_cast<size_t>(distance(perimeter_node_ptrs.begin(),lb)) ]++;
                    }
               }
             // eliminating those pointers from 'perimeter_node_ptrs' that are shared by multiple elements
             for ( size_t i{0}; i<face_count.size(); i++ )
               if ( face_count[i] > 1 ) perimeter_node_ptrs[i] = nullptr;
             perimeter_node_ptrs.erase( remove( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end(), nullptr ), perimeter_node_ptrs.end() );
          }
        // In 3D, does one need a triangulation
        if constexpr (dim == 3U ) {
             throw csmp::Exception( ERROR, "CreateInterfacesBetweenNodeSharingElements", "interface patch perimeter identification not implemented yet");
             MeshPatch<3U> patch( SURFACE );
             patch.BuildInterveningPatch( interface_nbor_elmts, false, perimeter_node_ptrs );
          }
    
      } // end multiplicate perimeter nodes
      
      
    // 3. Creating InterFace objects, duplicating nodes and dealing with boundaries
    // ----------------------------------------------------------------------------
    // tracking already duplicated nodes to avoid further duplication
    //  original,  duplicate
    map<Node<dim>*,Node<dim>*>  new_nodes;
    // vector of interfaces which will be returned
    vector<InterFace<dim>*>     interface_ptrs;
    interface_ptrs.reserve( interface_nbor_elmts.size() );
    
    for ( const auto& it : interface_nbor_elmts )
      {
         // 2.1 duplicating nodes but only if we are at a model boundary or the node already is a manifold
         // ----------------------------------------------------------------------------------------------
         // (if the node is not duplicated, the original node is inserted into the InterFace outside node vector)
         // - from the corner nodes of the element faces that will be at the interface, segment keys are made
         // - the faces are numbered
         // - segments that have element faces on either side, are inside the patch, the other ones are at the perimeter
         vector<Node<dim>*> outside_nodes( it.first.first->FE()->NodesPerFace( it.first.second ), nullptr );
         auto               nit{ new_nodes.end() };
         uint32_t           nd_count{0u};
         
         for ( const auto& i : it.first.first->FE()->NodesOfFace( it.first.second ) ) {
             // nodes are multiplicated, always if 'multiplicate_perimeter_nodes=true'
             // or if they do not lie on the perimeter of the new interface patch
             if ( perimeter_node_ptrs.empty() ) {
                  // nodes are duplicated unless they were already duplicated
                  if ( (nit=new_nodes.find(it.first.first->N(i))) == new_nodes.end() ) {
                       outside_nodes[nd_count] = Duplicate( it.first.first->N(i), nvars );
                       new_nodes.insert( make_pair( it.first.first->N(i), outside_nodes[nd_count] ) );
                    }
                  else outside_nodes[nd_count] = (*nit).second;
                }
              else {
                  // if perimeter nodes are excluded, nodes are duplicated if they do not lie on the perimeter
                  if ( !binary_search( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end(), it.first.first->N(i) ) ) {
                        // and only if these nodes have not already been duplicated
                        if ( (nit=new_nodes.find(it.first.first->N(i))) == new_nodes.end() ) {
                             outside_nodes[nd_count] = Duplicate( it.first.first->N(i), nvars );
                             new_nodes.insert( make_pair( it.first.first->N(i), outside_nodes[nd_count] ) );
                          }
                        else outside_nodes[nd_count] = (*nit).second;
                     }
                   // or if they are located on the model boundary or if the are already manifolds
                   else if ( it.first.first->N(i)->AtBoundary() != NOT || it.first.first->N(i)->IsManifold() ) {
                        if ( (nit=new_nodes.find(it.first.first->N(i))) == new_nodes.end() ) {
                             // NB: Duplicate adds the duplicated manifold nodes to the respective manifolds
                             outside_nodes[nd_count] = Duplicate( it.first.first->N(i), nvars );
                             new_nodes.insert( make_pair( it.first.first->N(i), outside_nodes[nd_count] ) );
                          }
                        else outside_nodes[nd_count] = (*nit).second;
                     }
                   // if the perimeter is considered perimeter nodes are just copied to the opposite side
                   else outside_nodes[nd_count] = it.first.first->N(i);
               }
             nd_count++;
           }

         // turning the nodes from counterclockwise to clockwise because they will go on the opposite side of the interface
         // reverse( outside_nodes.begin(), outside_nodes.end() );
         reverseOrderOfFaceNodes<dim>( static_cast<uint32_t>(it.first.first->CornerNodesOfFace( it.first.second ).size()),
                                       outside_nodes.begin(), outside_nodes.end() );

         // 2.3 construction of InterFace on the model perimeter
         // ----------------------------------------------------
         //     - higher-dimensional nbors are already known
         //     - faces of higher dimensional neighbors are also known
         //     - nodes on outside are not known (old nodes are on inside, new nodes on outside)
         //     - nodes on inside and outside are the same for perimeter interfaces away from boundaries
         //     - boundaries are inferred, when:
         //       - BOX_BOUNDARY flag is !NOT
         // NB: the nbor-elmt pointers corresponding to the shared face get nulled here 
         interface_ptrs.push_back( AddInterFace( it.first.first, it.first.second,
                                                 it.second.first, it.second.second,
                                                 lvars, ivars, outside_nodes ) );
      }

     // 3. cleaning up inter-CELL and node to parent connectivity
     // ---------------------------------------------------------
     // 3.1 find parent elements of old nodes that are in outside region, unassign them from old nodes and assign them to corresponding new nodes
// TODO: turn into function and test (it seems to miss the elements that do not have a face at the boundary)
     unordered_map<Node<dim>*, Element<dim>*> old_node_unassigned_element_map;
     for( auto& nd_pair : new_nodes ) {
       auto old_node = nd_pair.first;
       for (uint32_t i{0u}; i < old_node->Parents(); ++i) {
         auto parent_elmt = old_node->Parent(i);
         // for all outside elements (not just the equidimensional ones)
         if ( outside_elmts.find(parent_elmt) != outside_elmts.end() ) {
           auto local_nd_index = old_node->ParentNodeNumber(i);
           old_node_unassigned_element_map.insert(make_pair(old_node,parent_elmt));
           auto new_node = nd_pair.second;
           // assigning the new node to parent element
           parent_elmt->Assign(local_nd_index, new_node );
         }
       }
     }
     // 3.2 unassigning the old nodes from the parent element maps of the outside elements
     for( auto& old_node_unassigned_element : old_node_unassigned_element_map ) {
       auto old_node = old_node_unassigned_element.first;
       auto parent_elmt = old_node_unassigned_element.second;
       // removing the old node from the parent element vector
       old_node->Unassign(parent_elmt);
     }
     
    // connecting the new interfaces with one another
    BuildConnectivity<InterFace>( interface_ptrs.begin(), interface_ptrs.end() );

     // no entities are removed by this function
     
     cout <<"\n"<<"MeshManager<"<< dim <<">::CreateInterfacesBetweenNodeSharingElements: created "<< interface_ptrs.size() <<" new interfaces and ";
     cout << new_nodes.size() <<" new nodes."<< endl;
     
     return interface_ptrs;
     
  } // end CreateInterfacesBetweenNodeSharingElements








/**
    Reorganises 'outside' node pointer vector so that the locations of the nodes pointed to match those nodes
    in the 'inside' vector but in reverse order.
    
    TODO: method fails in DykePartiallySplit testcase, perhaps because of a tolerance issue
*/
template<uint32_t dim>
bool matchNodesByPosition( const vector<Node<dim>*>& inside_nodes,
                           vector<Node<dim>*>& outside_nodes )
 {
    assert( !inside_nodes.empty() );
    assert( inside_nodes.size() ==  outside_nodes.size() );
    
    vector<Node<dim>*>  temp;
    temp.reserve( inside_nodes.size() );
    
    for ( auto& iit : inside_nodes )
      for ( auto& oit : outside_nodes ) {
           // if less than fails twice the points are assumed to be equal
           if ( !(oit->Coordinate() < iit->Coordinate()) &&
                !(oit->Coordinate() > iit->Coordinate()) ) {
                temp.push_back( oit );
                break;
             }
        }
    
    reverse( temp.begin(), temp.end() );
    outside_nodes = temp;
    
    return ( temp.size() == inside_nodes.size() );
    
 } // end matchNodesByPosition

template bool matchNodesByPosition( const vector<Node<3U>*>&, vector<Node<3U>*>& );



/**
    Creates InterFace objects between face/node sharing Elements adding the necessary node manifolds and InterFace connectivity; inside elements are first in pair
*/
template<uint32_t dim>
vector<InterFace<dim>*>  MeshManager<dim>::CreateInterfacesBetweenNodeMatchingElements( const PropertyDatabase<dim>& dbase,
                                           const vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > >& interface_nbor_elmts )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    if ( interface_nbor_elmts.empty() ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::CreateInterfacesBetweenNodeMatchingElements",
                                   "supplied range of element pairs is empty; nothing was done.");
         return vector<InterFace<dim>*>{}; // empty vec
      }
    
    // 1. Constructing the node manifold manager if necessary
    // ------------------------------------------------------
    const bool no_previous_manifolds = ( node_manifold_manager_ == nullptr ) ? true : false;
    if ( no_previous_manifolds )
      node_manifold_manager_ = new NodeManifoldManager<dim>();
    
    const LocalVariables             nvars(dbase.LocalVariablesAt(NODE));
    const LocalVariables             lvars(dbase.LocalVariablesAt(INTER_FACE));
    const IntegrationPointVariables& ivars(dbase.IntegrationPointVariablesAt(INTER_FACE));
        
    // 2. Creating InterFace objects
    // -----------------------------
    // (all nodes are already duplicated)
    // vector of interfaces which will be returned
    vector<InterFace<dim>*>  interface_ptrs;
    interface_ptrs.reserve( interface_nbor_elmts.size() );
    // unique set of node pairs needed to create the manifolds later on
    // (inside,outside)
    set<pair<Node<dim>*,Node<dim>*> > node_ptr_pairs;

    for ( const auto& it : interface_nbor_elmts )
      {
         // 2.1 Collecting node pairs to form manifolds and outside nodes to construct interface
         // ------------------------------------------------------------------------------------
         const auto n_face_nodes{ it.first.first->FE()->NodesPerFace( it.first.second ) };
         vector<Node<dim>*> inside_nodes, outside_nodes;
         inside_nodes.reserve( n_face_nodes );
         outside_nodes.reserve( n_face_nodes );
         
         for ( const auto& i : it.first.first->FE()->NodesOfFace( it.first.second ) )
           inside_nodes.push_back( it.first.first->N(i) );
         for ( const auto& i : it.second.first->FE()->NodesOfFace( it.second.second ) )
           outside_nodes.push_back( it.second.first->N(i) );
         
         // organising the interface nodes in the outside vector such that they match the inside ones by position
         const bool all_nodes_matched = matchNodesByPosition( inside_nodes, outside_nodes );
         if ( all_nodes_matched == false ) {
              cerr <<"\n\n"<<"interface between elements "<< it.first.first->Idx() <<" and "<< it.second.first->Idx();
              cerr <<" with node coordinates (inside vs. outside):\n";
              for ( uint32_t i{0U}; i<n_face_nodes; i++ ) {
                   cerr <<"\t"<< i <<": "<< inside_nodes[i]->Coordinate() <<" vs. ";
                   if ( (n_face_nodes-i-1U) < outside_nodes.size() ) cerr << outside_nodes[n_face_nodes-i-1U]->Coordinate();
                   else cerr <<"no matching node found.";
                   cerr << endl;
                }
              cerr << endl;
              csmp_error.Note( WARNING, "MeshManager<dim>::CreateInterfacesBetweenNodeMatchingElements",
                                        "nodes at interface between region could not be matched.");
           }
         else {
             // and storing them as pairs that will become manifolds
             for ( uint32_t i{0U}; i<n_face_nodes; i++ )
               node_ptr_pairs.insert( make_pair( inside_nodes[i],
                                                 outside_nodes[n_face_nodes-i-1U] ) );
                
             // 2.2 constructing InterFace objects
             // ----------------------------------
             //     - higher-dimensional nbors are already known
             //     - faces of higher dimensional neighbors are also known
             //     - nodes on inside are deduced by constructor, outside nodes are supplied as 'outside_nodes'
             interface_ptrs.push_back( AddInterFace( it.first.first, it.first.second,
                                                     it.second.first, it.second.second,
                                                     lvars, ivars, outside_nodes ) );
          }
      }

    // 3. Creating the node manifolds
    // ------------------------------
    for ( const auto& nit : node_ptr_pairs )
      {
//cerr << nit.first->Idx() <<"--"<< nit.second->Idx() <<" ";
          // if inside or outside nodes already are manifolds, the non-manifold nodes are added to them
          // (Note: Add() also assigns the argument node to this manifold)
          if ( nit.first->IsManifold() && !nit.second->IsManifold() )
             nit.first->Manifold()->Add( nit.second );
          else if ( !nit.first->IsManifold() && nit.second->IsManifold() )
             nit.second->Manifold()->Add( nit.first );
          else {
               // a new manifold is created using the provided default geometric classifier
               auto nmf = node_manifold_manager_->AddManifold( nodes_,
                                                               nit.first,
                                                               nit.second,
                                                               ManifoldType::SPLIT_BOUNDARY );
               // and its nodes are connected to it
               nit.first->Assign( (*nmf) );
               nit.second->Assign( (*nmf) );
            }

           // working out whether the original classification as an interface was correct
#ifdef DEBUG
           const bool verbose{false};
           consistencyCheck( (*nit.first->Manifold()), verbose );
#endif
        }

     // 4. cleaning up inter-CELL and node to parent connectivity
     // ---------------------------------------------------------
     BuildConnectivity<InterFace>( interface_ptrs.begin(), interface_ptrs.end() );
     
     cout <<"\n"<<"MeshManager<"<< dim <<">::CreateInterfacesBetweenNodeMatchingElements: created "<< interface_ptrs.size();
     cout <<" new interfaces."<< endl;
     
     return interface_ptrs;
     
 } // end CreateInterfacesBetweenNodeMatchingElements


// DEBUGGING - OK
/*
cerr <<"\n"<<"inside-outside matching node points:\n";
cerr << it.first.first->N(inside_fnids[i])->Coordinate() <<" ";
cerr << it.second.first->N(outside_fnids[n_face_nodes-i-1U])->Coordinate();

*/



// ========================================================================================================================

// ERASURES

// ========================================================================================================================



// TODO: by comparison with a traversal of a node graph, this includes many redundant binary searches
/**
    Erases the nodes in the MeshManager pointed to by the supplied range of pointers.
    Only nodes that do not have parent elements are deleted. If there are nodes that are still interconnected, they are reported.
    Method returns the number of erasures made.
    Any potential node manifolds are updated.
    
    @attention checks whether the nodes to delete are still connected to any elements, reporting them if this is the case.
    Whatch out! - this check is not done for the nodes that belong to Face or InterFace objects because nodes do not keep any parent info fom them.
*/
template<uint32_t dim>
size_t MeshManager<dim>::DeleteNodesAndRepairNodeConnnectivity( typename vector<Node<dim>*>::iterator first,
                                                                typename vector<Node<dim>*>::iterator last )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

     auto n_nodes_to_delete( static_cast<size_t>(distance(first,last)) );
 
     if ( n_nodes_to_delete == 0 ) return 0U;

     // 1. making sure that the nodes do not belong to any cells anymore
     // ----------------------------------------------------------------
     vector<Node<dim>*> nodes_to_delete; nodes_to_delete.reserve( n_nodes_to_delete );
     set<Node<dim>*>    nodes_still_connected;
     
     while ( first != last ) {
           assert( (*first) != nullptr );
           bool still_connected_node_found{ false };
           const uint32_t n_parents{ (*first)->Parents() };
           for ( uint32_t i{0u}; i<n_parents; i++ ) {
                 assert( (*first)->Parent(i) );
                 if ( (*first)->Parent(i)->N( (*first)->ParentNodeNumber(i) ) != nullptr ) {
                      nodes_still_connected.insert( (*first) );
                      still_connected_node_found = true;
                      n_nodes_to_delete--;
                      break;
                   }
               }
           if ( !still_connected_node_found ) nodes_to_delete.push_back( (*first) );
           first++;
        }
     // reporting
     if ( !nodes_still_connected.empty() ) {
          cout <<"\n"<<"MeshManager<dim>::DeleteNodesAndRepairNodeConnnectivity: ";
          cout <<"you are trying to delete nodes still connected to elements:"<< endl;
          for ( const auto& nit : nodes_still_connected ) {
               cout <<"\t\t"<< nit->Idx() <<": "<< nit->Coordinate() << endl;
            }
          csmp_error.Note( WARNING, "MeshManager<dim>::DeleteNodesAndRepairNodeConnnectivity",
                          "detected non-eligible nodes; continuing deletion without these." );
            
          // are there any nodes left to delete?
          if ( n_nodes_to_delete <= 0 ) return 0ul;
       }


     // 2. disconnecting neighbor nodes from the nodes that will be removed
     // -------------------------------------------------------------------
     // but only for neighbor nodes which are not within the range of the nodes to be deleted
     sort( nodes_to_delete.begin(), nodes_to_delete.end() );
     for ( const auto& nit : nodes_to_delete ) {
        const uint32_t n_nbors{ nit->Neighbors() };
        for ( uint32_t i{0U}; i<n_nbors; i++ ) {
              const uint32_t n_nbor_nbors{ nit->Neighbor(i)->Neighbors() };
              for ( uint32_t j{0U}; j<n_nbor_nbors; j++ )
                if ( binary_search( nodes_to_delete.begin(), nodes_to_delete.end(), nit->Neighbor(i)->Neighbor(j) ) ) {
                     // neighbor node only moves current node beyond the end of the neighbors container
                     nit->RemoveNeighbor( nit );
                     break;
                  }
          }
      }

     // 3. removing nodes from potential manifolds before removing themselves
     // ---------------------------------------------------------------------
     for ( auto& nit : nodes_to_delete ) {
          // disconnecting the node from its manifold
          if ( nit->IsManifold() ) {
               nit->Manifold()->Remove( nit );
               // if the node no longer is a manifold, it is removed from the node-manifold manager
               if ( nit->Manifold()->Branches() == 1U )
               node_manifold_manager_->Delete( nit->Manifold() );
            }
          nodes_.erase( nodes_.get_iterator(nit) );
       }
     
     return static_cast<size_t>(n_nodes_to_delete);
    
 } // end DeleteNodesAndRepairNodeConnnectivity






/**
    Erases range of elements not counting  null-pointer cells in the supplied sequence,
    returning the number of genuine erasures.
    
    @param first iterator to first element of a range of either of line, surface, or volume  elements
    @param last end of range of same celltype elements
    
    Also checks whether element deletion causes orphan nodes. If so, these are deleted as well.
    This check involves counting the node's parent elements that are not null pointers to make
    sure that the true state of the node is captured.
    
    @attention method assumes that the neighbor connectivity of the elements is valid and Region is intact (no corrupt pointers)

    The following steps are performed:

    1. set the neighbor pointers to the elements that will be deleted to zero (= disconnect the neighbor elements)
    2. (remove the pointers from the connected nodes to this parent-element) - done later sweeping over the nodes (Delete())
    3. delete elements
    
    @attention us this method only when the iterator ranges spans all the elements that are going to be deleted, else undefined behaviour may result.
*/
template<uint32_t dim>
size_t MeshManager<dim>::DeleteElementsAfterDisconnectingRemainingOnes( typename vector<Element<dim>*>::iterator first,
                                                                        typename vector<Element<dim>*>::iterator last )
 {
     long elmts_to_delete( distance(first,last) );
 
     if ( elmts_to_delete == 0 ) return 0U;
     
     // 1. distinguishing 2 cases: 1) equidimensional elements, and 2) lower-dim elements that share their nodes equidim ones
     // ---------------------------------------------------------------------------------------------------------------------
     updateHaloCellConnectivity<dim,Element>(first,last);

     // 2. deleting the elements, setting pointers in input range to zero
     // -----------------------------------------------------------------
     size_t deleted_elements{ 0U };
     while ( first != last ) {
          // updates node-parent element connectivity
          Delete( first );
          // checking whether deletion was successful
          if ( (*first) == nullptr ) deleted_elements++;
          first++;
       }

    return deleted_elements;
    
 } // end DeleteElementsAfterDisconnectingRemainingOnes



/* SUPERSEDED VERSION

template<uint32_t dim>
size_t MeshManager<dim>::DeleteElementsAndRepairConnnectivity( typename vector<Element<dim>*>::iterator first,
                                                               typename vector<Element<dim>*>::iterator last )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     
     throw csmp::Exception( ERROR, "MeshManager<dim>::DeleteElementsAndRepairConnnectivity", "Method deletes Elements and is not tested yet");

     long elmts_to_delete( distance(first,last) );
 
     if ( elmts_to_delete == 0 ) return 0U;
     
     // 0. distinguishing 2 cases: 1) equidimensional elements, and 2) lower-dim elements that share their nodes equidim ones
     // ---------------------------------------------------------------------------------------------------------------------
     const CELL_SHAPE cell_shape = parseFiniteElementDimension( (*first)->FE_Type() );
     if ( parseFiniteElementDimension( (*prev(last,1))->FE_Type() ) != cell_shape ) {
          csmp_error.Note( ERROR, "MeshManager<dim>::DeleteElementsAndRepairConnnectivity",
                            "range of supplied elements appear to be of different cell shape; nothing was done" );
          return 0U;
       }
     const bool elmts_are_equidimensional = ( (cell_shape == VOLUME  && dim == 3U) ||
                                              (cell_shape == SURFACE && dim == 2U) ) ? true : false;
     // for finding perimeter
     vector<Element<dim>*> element_ptrs( first, last );
     sort( element_ptrs.begin(), element_ptrs.end() );


     // 1. if the elements have the same dimensions as the model
     // --------------------------------------------------------
     // - nodes in the inside loose their parents and will be deleted
     // - nodes on the perimeter will be updated
     // - the elements will be deleted
     if ( elmts_are_equidimensional )
       {
         // 1.1  finding perimeter nodes for updating and disconnecting element neighbors across perimeter faces
         vector<Element<dim>*> adjacent_elmts;
         vector<Node<dim>*>    interior_nodes;
         adjacent_elmts.reserve( static_cast<size_t>(elmts_to_delete/2) ); // elements outside of range but adjacent to the interior ones
         interior_nodes.reserve( static_cast<size_t>(elmts_to_delete) );   // nodes that are not on such boundaries
         auto first1{ first };
         while ( first != last ) {
              assert( (*first) != nullptr );
              bool inside_elmt{ true };
              // an element face with no neighbor is a perimeter face
              for ( uint32_t i{0U}; i<(*first)->Neighbors(); i++ )
                // if the element neighbor exists and is not part of the elements that shall be deleted
                if ( (*first)->Neighbor(i) != nullptr &&
                     !binary_search( element_ptrs.begin(), element_ptrs.end(), (*first)->Neighbor(i) ) )
                  {
                      // - recording abbuting elements on the outside needed later for updating node-parent connectivity
                      adjacent_elmts.push_back( (*first)->Neighbor(i) );
                      // - disconnecting abutting neighbor elements from the ones that will be deleted
                      for ( uint32_t k{0U}; k<(*first)->Neighbor(i)->Neighbors(); k++ )
                        if ( (*first)->Neighbor(i)->Neighbor(k) == (*first) )
                          (*first)->Neighbor(i)->Assign( k, static_cast<Element<dim>*>(nullptr) );
                      
                  }
                // - noting the element place
                else inside_elmt = false;
              // collecting the inside nodes for deletion
              if ( inside_elmt ) {
                   for ( uint32_t i{0U}; i<(*first)->Nodes(); i++ )
                      interior_nodes.push_back( (*first)->N(i) );
                }
              first++;
           }
           
          // deleting potential duplicates from the element and node vectors
          sort( interior_nodes.begin(), interior_nodes.end() );
          interior_nodes.erase( unique( interior_nodes.begin(), interior_nodes.end() ), interior_nodes.end() );
          sort( adjacent_elmts.begin(), adjacent_elmts.end() );
          adjacent_elmts.erase( unique( adjacent_elmts.begin(), adjacent_elmts.end() ), adjacent_elmts.end() );
             
          // 1.2 deleting the elements, setting pointers to zero
          size_t deleted_elements{ 0U };
          while ( first1 != last ) {
               Delete( first1 );
               deleted_elements++;
               first1++;
            }
          
          // 1.3 deleting the interior nodes, updating neighbor connectivity with perimeter ones
          DeleteNodesAndRepairNodeConnnectivity( interior_nodes.begin(), interior_nodes.end() );
          
          // 1.4 updating the perimeter nodes
          ConnectNodesToParentsAndNeighbors( adjacent_elmts.begin(), adjacent_elmts.end() );
          
          return deleted_elements;
         
       } // end equidimensional
       
       
     // 2. if the elements are lower-dimensional all their nodes are on the perimeter
     // ------------------------------------------------------------------------------
     // - no nodes have to be deleted
     // - the element parents of all nodes are updated to account for the deleted elements
     // - the elements will be deleted
     vector<Element<dim>*> adjacent_elmts;
     adjacent_elmts.reserve( static_cast<size_t>(elmts_to_delete/2) );
     auto first1{ first };
     while ( first != last )
       {
          // 2.1 disconnecting abbuting elements, but remembering them for updating
          for ( uint32_t i{0U}; i<(*first)->Neighbors(); i++ )
            if ( (*first)->Neighbor(i) != nullptr &&
                 !binary_search( element_ptrs.begin(), element_ptrs.end(), (*first)->Neighbor(i) ) )
              {
                  // - recording abbuting elements on the outside needed later for updating node-parent connectivity
                  adjacent_elmts.push_back( (*first)->Neighbor(i) );
                  // - disconnecting these neighbor elements from the ones that will be deleted
                  for ( uint32_t k{0U}; k<(*first)->Neighbor(i)->Neighbors(); k++ )
                    if ( (*first)->Neighbor(i)->Neighbor(k) == (*first) )
                      (*first)->Neighbor(i)->Assign( k, static_cast<Element<dim>*>(nullptr) );
                  
              }
          first++;
       }

     // 2.2 deleting the elements, setting pointers to zero
     size_t deleted_elements{ 0U };
     while ( first1 != last ) {
          Delete( first1 );
          deleted_elements++;
          first1++;
       }

    // 2.3 updating the parent connectivity of nodes as well as the neighbor connectivity
    ConnectNodesToParentsAndNeighbors( adjacent_elmts.begin(), adjacent_elmts.end() );

    return deleted_elements;
    
 } // end DeleteElementsAndRepairConnnectivity

END SUPERSEDED VERSION */







/**
      Deletes range of Faces after detecting and disconnecting potential neighbor faces around the perimeter of the face patch.
*/
template<uint32_t dim>
size_t MeshManager<dim>::DeleteFacesAndRepairConnnectivity( typename vector<Face<dim>*>::iterator first,
                                                            typename vector<Face<dim>*>::iterator last )
 {
//     throw csmp::Exception( ERROR, "MeshManager<dim>::DeleteFacesAndRepairConnnectivity", "method not tested yet");

     auto faces_to_delete( static_cast<size_t>(distance(first,last)) );
 
     if ( faces_to_delete == 0 ) return 0U;
     
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     
     // disconnecting the Face objects that surround the Face patch that will be deleted
     updateHaloCellConnectivity<dim,Face>( first, last );

     size_t deleted_faces{0ul};
     size_t nullptr_faces{0ul};
     
     while ( first != last )
       {
          Delete( first );
       // ^^^^^^^^^^^^^^^^
          if ( (*first) == nullptr ) deleted_faces++;
          else nullptr_faces++;
          first++;
       }
     
     if ( nullptr_faces > 0 ) {
          cout <<"\n\t"<<"detected "<< nullptr_faces <<" nullptr faces in the input interator range.";
          csmp_error.Note( ERROR, "MeshManager<dim>::DeleteFacesAndRepairConnnectivity",
                           "deleted less faces than in input range");
       }
     else if ( deleted_faces != faces_to_delete )
       csmp_error.Note( ERROR, "MeshManager<dim>::DeleteFacesAndRepairConnnectivity",
                       "deleted less faces than in input range");
     
     
     return deleted_faces;
    
 } // end DeleteAndRepairConnnectivity(Face)







/**
    Deletes range of InterFaces after detecting and disconnecting potential neighbor faces around the perimeter of the face patch.
    This is accimplished by calling 'updateHaloCellConnectivity'.
    Then cells and duplicated nodes are removed and the mesh is reconnected across the fomer SplitBoundary.
    
    @Note IMPORTANT: any potential intervening Element regions should be deleted because they will no longer be connected to any of the elements in the rest of the mesh.
    This MeshManager method cannot do this because it has no access to the RegionInterFace which would be required for deleting the associated region as well.
    
    @attention IMPORTANT - after this operation, any (unique or non-unique) outside regions must be updated because their node vector will have changed!
    
    @attention method does not reconnect the mesh where interfaces are removed.
*/
template<uint32_t dim>
size_t MeshManager<dim>::DeleteInterfacesAndRepairConnnectivity( typename vector<InterFace<dim>*>::iterator first,
                                                                 typename vector<InterFace<dim>*>::iterator last )
 {
//     throw csmp::Exception( ERROR, "MeshManager<dim>::DeleteInterfacesAndRepairConnnectivity", "method not tested yet");
  
     auto interfaces_to_delete( static_cast<size_t>(distance(first,last)) );
 
     if ( interfaces_to_delete == 0 ) return 0U;
     
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );

#ifndef NDEBUG
     // checking for nullptr interfaces
     auto n_nullptrs_in_input = count_if( first, last, [](InterFace<dim>* const ptr) { return ptr == nullptr; } );
     if ( n_nullptrs_in_input > 0 )
       csmp_error.Note( ERROR, "MeshManager<dim>::DeleteInterfacesAndRepairConnnectivity", to_string(n_nullptrs_in_input),
                       "'nullptr' InterFace objects counted in input range; have these InterFaces already been deleted?");
#endif

     // 0. disconnecting the InterFace objects that surround the InterFace patch that will be deleted
     // ---------------------------------------------------------------------------------------------
     updateHaloCellConnectivity<dim,InterFace>( first, last );
     

     // 1. processing potential lower-dimensional elements sandwiched between interface sides
     //    reconnecting higher-dimensional neighbor elements with one-another and dealing with duplicated nodes
     // -------------------------------------------------------------------------------------------------------
     vector<Node<dim>*> nodes_to_delete;
     nodes_to_delete.reserve( interfaces_to_delete );
     for_each( first, last, [&nodes_to_delete, &csmp_error, this]( InterFace<dim>* const ptr) {
                                // 1.1 detaching intervening elements
                                // ----------------------------------
                                if ( ptr->HasInterveningElement() )
                                  ptr->UnAssignInterveningElement();
                                // 1.2 reconnecting higher-dimensional neighbors
                                // ---------------------------------------------
                                ptr->InnerParent()->Assign( ptr->InnerParentFaceID(), ptr->OuterParent() );
                                ptr->OuterParent()->Assign( ptr->OuterParentFaceID(), ptr->InnerParent() );
                                // 1.3 collecting a unique set of outer element nodes that will be deleted if they are manifold nodes
                                // --------------------------------------------------------------------------------------------------
                                const uint32_t n_nodes{ ptr->FE()->Nodes() };
                                for ( uint32_t i{0u}; i<n_nodes; ++i )
                                  // if the node is a manifold with at least 2 nodes within it
                                  if ( ptr->N(i,OUTSIDE)->IsManifold() && ptr->N(i,OUTSIDE)->Manifold()->Branches() > 1 ) {
                                       nodes_to_delete.push_back(ptr->N(i,OUTSIDE));
                                       // 1.4 removing outside node from manifold to prepare it for deletion
                                       // ------------------------------------------------------------------
                                       if ( ptr->N(i,OUTSIDE)->Manifold()->Branches() > 1 ) {
                                            ptr->N(i,OUTSIDE)->Manifold()->Remove( ptr->N(i,OUTSIDE) );
                                         }
                                       // the node no longer is a manifold and therefore is removed from node-manifold manager
                                       else node_manifold_manager_->Delete( ptr->N(i,OUTSIDE)->Manifold() );
                                       // 1.5 Disconnecting the node from its neighbors
                                       // ---------------------------------------------
                                       for ( auto nit=ptr->N(i,OUTSIDE)->NeighborsBegin(); nit!=ptr->N(i,OUTSIDE)->NeighborsEnd(); ++nit )
                                         if ( (*nit) )
                                            (*nit)->RemoveNeighbor( ptr->N(i,OUTSIDE) );
                                    }
                                // 1.6 replacing the face nodes of the outside element with those from the inside
                                // ------------------------------------------------------------------------------
                                uint32_t node_counter{0u};
                                for ( auto& fnid : ptr->OuterParent()->FE()->NodesOfFace( ptr->OuterParentFaceID() ) ) {
                                     if ( ptr->OuterParent()->N(fnid) != ptr->N(node_counter,OUTSIDE) ) {
                      #ifndef NDEBUG
                                          cout <<"\n\t"<<"interface node: "<< ptr->N(node_counter,OUTSIDE)->Idx() <<" and ";
                                          cout <<"outer parent face node: "<< ptr->OuterParent()->N(fnid)->Idx() <<" do not match.";
                                          csmp_error.Note( WARNING, "MeshManager<dim>::DeleteInterfacesAndRepairConnnectivity",
                                                          "detected node-numbering inconsistency between InterFace and outside Element face; nodes could not be fused");
                      #endif
                                          continue;
                                       }
                                     // inverting the node sequence to get the correct matching node
                                     ptr->OuterParent()->Assign( fnid, ptr->N( ptr->FE()->Nodes() - 1 - node_counter, INSIDE) );
                                     // 1.7 adding the newly connected outside element to the node parents
                                     // ------------------------------------------------------------------
                                     ptr->OuterParent()->N(fnid)->Assign( fnid, ptr->OuterParent() );
                                     // testing
                                     // cout <<"\n"<<"outside node: "<< node_counter;
                                     // cout <<": "<< ptr->N(node_counter,OUTSIDE)->Coordinate() <<" vs inside node: ";
                                     //cout << ptr->N(node_counter,INSIDE)->Coordinate();
                                     node_counter++;
                                  }
                                  
                              } );
       
     // 2. deleting disconnected outside nodes that are no longer needed after the InterFace nodes were fused
     // -----------------------------------------------------------------------------------------------------
     // (eliminating potential duplicate nodes from 'nodes_to_delete' vector)
     sort( nodes_to_delete.begin(), nodes_to_delete.end() );
     nodes_to_delete.erase( unique( nodes_to_delete.begin(), nodes_to_delete.end() ), nodes_to_delete.end() );
#ifndef NDEBUG
     vector<size_t> deleted_nodes;
     deleted_nodes.reserve( nodes_to_delete.size() );
#endif
     for ( auto nit=nodes_to_delete.begin(); nit!=nodes_to_delete.end(); ++nit ) {
           // tracking the IDs of the deleted nodes
           size_t node_idx = (*nit)->Idx();
           this->Delete( nit );
#ifndef NDEBUG
           if ( (*nit) == nullptr ) deleted_nodes.push_back( node_idx );
#endif
       }

     // 3. deleting the interfaces and nulling the pointers to them
     // -----------------------------------------------------------
     size_t deleted_interfaces{ static_cast<size_t>(interfaces_to_delete) };
     while( first != last ) {
          Delete( first );
//        ^^^^^^^^^^^^^^^
          if ( (*first) != nullptr ) deleted_interfaces--;
          first++;
       }

#ifndef NDEBUG
     if ( !deleted_nodes.empty() ) {
          cout <<"\nMeshManager::DeleteInterfacesAndRepairConnnectivity: deleted nodes:\n\t";
          sort( deleted_nodes.begin(), deleted_nodes.end() );
          for ( const auto& nit : deleted_nodes ) cout << nit <<" ";
          cout << endl;
       }
     if ( deleted_interfaces < interfaces_to_delete ) {
          cout <<"\nMeshManager::DeleteInterfacesAndRepairConnnectivity: deleted less interfaces than intended:\n\t";
          cout << deleted_interfaces <<" vs "<< interfaces_to_delete << endl;
       }
#endif

     return deleted_interfaces;
    
 } // end DeleteInterfacesAndRepairConnnectivity

/* TESTING
               if ( (*first)->MatchingN( node_counter, OUTSIDE ) != (*first)->OuterParent()->N( fnid ) ) {
                    cout <<"\n"<<"interface node pointers: ";
                    for ( uint32_t i{0U}; i<(*first)->FE()->Nodes(); i++ ) cout <<" "<< (*first)->N(i,INSIDE);
                    for ( uint32_t i{0U}; i<(*first)->FE()->Nodes(); i++ ) cout <<" "<< (*first)->N(i,OUTSIDE);
                    cout << endl;
                    cout <<"\noutside element face node: "<< (*first)->OuterParent()->N( fnid );
                    cout <<"\ninside element face node:  "<< (*first)->N(node_counter);
                    cout <<"\noutside matching node:     "<< (*first)->MatchingN( node_counter, OUTSIDE );
                    cout <<"\ninside matching node:      "<< (*first)->MatchingN( node_counter, INSIDE );
                    cout << endl;
                 }
*/
/* REPLACED CODE
          // 1. detaching potential lower-dimensional elements sandwiched between interface sides
          // ------------------------------------------------------------------------------------
          // (these elements and their nodes are recorded for deletion further below)
          if ( (*first)->HasInterveningElement() ) {
#ifndef NDEBUG
                if ( first_call ) {
                     cout <<"\n\t"<<"detaching intervening Element: "<< (*first)->InterveningElement()->Idx() <<": ";
                     cout << parseFiniteElementType( (*first)->InterveningElement()->FE_Type() );
                     csmp_error.Note( WARNING, "MeshManager<dim>::DeleteInterfacesAndRepairConnnectivity",
                                     "detaching intervening lower-dimensional Element");
                     first_call = false;
                  }
#endif
                for ( auto nit=(*first)->InterveningElement()->NodesBegin(); nit!=(*first)->InterveningElement()->NodesEnd(); ++nit )
                  nodes_to_delete.push_back( (*nit) );

                // detaching the intervening element
                // =================================
                (*first)->UnAssignInterveningElement();
            }

*/






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
    const auto     n_cells_max = distance(first,last);
    if ( n_cells_max == 0 ) {
         csmp_error.Note( WARNING, "MeshManager<dim>::BuildConnectivity:", "supplied cell vector is empty; nothing was done." );
         return;
      }
    //cout << "\nMeshManager<"<< dim <<">::BuildConnectivity: Establishing CSMP FE neighbor connectivity...\n";
    
    // making subranges for the volume, surface, and line elements
    vector<CELL<dim>*> volume_cells, surface_cells, line_cells;
    const auto cellsEnd{last};
    
    // if we are dealing with element connectivity in 3D
    if constexpr ( is_same< CELL<dim>,Element<dim> >::value ) {
        if constexpr( dim == 3U ) volume_cells.reserve( static_cast<size_t>(n_cells_max) );
        if constexpr( dim != 1U ) surface_cells.reserve( static_cast<size_t>(n_cells_max/3) );
        line_cells.reserve( static_cast<size_t>(n_cells_max/6) );
        while ( first != cellsEnd ) {
             if constexpr( dim == 3U ) if ( (*first)->FE()->IsVolume() )  volume_cells.push_back(*first);
             if constexpr( dim != 1U ) if ( (*first)->FE()->IsSurface() ) surface_cells.push_back(*first);
             if ( (*first)->FE()->IsLine() ) line_cells.push_back(*first);
             first++;
          }
        // connecting the elements found
        if constexpr( dim == 3U )
          if ( !volume_cells.empty() )
            BuildVolumeConnectivity<Element>( volume_cells.begin(), volume_cells.end() );
            
        if ( !surface_cells.empty() ) BuildSurfaceConnectivity<Element>( surface_cells.begin(), surface_cells.end() );
        if ( !line_cells.empty() )    BuildLineConnectivity<Element>( line_cells.begin(), line_cells.end() );
      }
      
    // else
    else { // ( is_same< CELL<dim>,Face<dim> >::value || is_same< CELL<dim>,InterFace<dim> >::value  )
        // dim-1 case
        if constexpr ( dim != 1U ) {
             surface_cells.reserve( static_cast<size_t>(n_cells_max) );
             line_cells.reserve( static_cast<size_t>(n_cells_max/6) );
             while ( first != cellsEnd ) {
                  if ( (*first)->FE()->IsSurface() )   surface_cells.push_back(*first);
                  else if ( (*first)->FE()->IsLine() ) line_cells.push_back(*first);
                  first++;
               }
             // connecting the elements found
             if ( !surface_cells.empty() ) BuildSurfaceConnectivity<CELL>( surface_cells.begin(), surface_cells.end() );
             if ( !line_cells.empty() )    BuildLineConnectivity<CELL>( line_cells.begin(), line_cells.end() );
          }
        // 1D case
        else { // ( dim == 1U )
             line_cells.reserve( static_cast<size_t>(n_cells_max) );
             while ( first != cellsEnd ) {
                  line_cells.push_back(*first);
                  first++;
               }
             if ( !line_cells.empty() ) BuildLineConnectivity<CELL>( line_cells.begin(), line_cells.end() );
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
void MeshManager<3>::BuildVolumeConnectivity( typename vector<CELL<3U>*>::const_iterator first,
                                              typename vector<CELL<3U>*>::const_iterator last )
   {
      ErrorHandler& csmp_error( ErrorHandler::Instance() );
  
#ifndef NDEBUG
      if ( (detectDuplicateCells<3U,CELL>( first, last, true )) > 0 )
        csmp_error.Note( WARNING, "MeshManager<dim>::BuildVolumeConnectivity", "detected duplicate cells in iterator range");
#endif

      // this method applies only to volumetric elements and 3D
      if constexpr ( is_same< CELL<3U>,Element<3U> >::value ) {
           // creating search keys from the corner nodes of the element faces
           // corner-nodes      elements that share face and their face id
           map<set<Node<3>*>,map<Element<3U>*,uint32_t> >  elmt_pairs;
           const auto                                      elementsEnd{last};

           // pairing the elements up in the search map
           while ( first != elementsEnd ) {
                assert( (*first) != nullptr );
                const uint32_t n_faces{ (*first)->Faces() };
                for ( uint32_t face{0U}; face < n_faces; ++face ) {
                     // trying to insert it into the map
                     auto it = elmt_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face),
                                                             map<Element<3U>*,uint32_t>{make_pair(*first,face)} )
                                                );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // ATTENTION: we must to null current neighbor connectivity out of reach because
                     // this creates potential disconnections within the mesh
                     //(*first)->Assign( face, static_cast<CELL<3>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the elements to one another
           for ( auto& it : elmt_pairs ) {
                size_t n_face_nbors{ it.second.size() };
                if ( n_face_nbors == 2U ) {
                     Element<3U>* const eptr1 = (*it.second.begin()).first;
                     Element<3U>* const eptr2 = (*it.second.rbegin()).first;
                     const uint32_t face_e1  = (*it.second.begin()).second;
                     const uint32_t face_e2  = (*it.second.rbegin()).second;
                     eptr1->Assign( face_e1, eptr2 );
                     eptr2->Assign( face_e2, eptr1 );
                  }
                // else no assignments have to be made as there is no neighbor
                assert( n_face_nbors <= 2U );
             }
        }
   
   } // end BuildVolumeElementConnectivity

template void MeshManager<3U>::BuildVolumeConnectivity<Element>( typename vector<Element<3U>*>::const_iterator,
                                                                 typename vector<Element<3U>*>::const_iterator );







/**
    Connects neighboring cells, disambiguating the potential manifolds by choosing co-planar elements.
    
    This is accomplished by finding the angle between 2 suface elements in 3D,  returning the acute angle in degrees (0..90o).
*/
template<uint32_t dim>
template<template<uint32_t> class CELL>
void MeshManager<dim>::BuildSurfaceConnectivity( typename vector<CELL<dim>*>::const_iterator first,
                                                 typename vector<CELL<dim>*>::const_iterator last )
   {
      ErrorHandler& csmp_error( ErrorHandler::Instance() );
  
#ifndef NDEBUG
      if ( (detectDuplicateCells<dim,CELL>( first, last, true )) > 0 )
        csmp_error.Note( WARNING, "MeshManager<dim>::BuildSurfaceConnectivity", "detected duplicate cells in iterator range");
#endif
      // this method applies only to surface elements in 3D
      if constexpr ( dim == 3U ) {
           // creating search keys from the corner nodes of the element faces
           // corner-nodes      elements that share face and their face id
           map<set<Node<3>*>,map<CELL<3>*,uint32_t> >  cell_pairs;

           // pairing the cells up in the search map
           const auto cellsEnd{last};
           while ( first != cellsEnd ) {
                assert( (*first) != nullptr );
                const uint32_t n_faces{ (*first)->Faces() };
                for ( uint32_t face{0u}; face < n_faces; ++face ) {
                     // trying to insert cell into the map using a search key of node pointers
                     pair<typename map<set<Node<3>*>,map<CELL<3>*,uint32_t> >::iterator,bool>
                       it = cell_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<CELL<3>*,uint32_t>{{*first,face}} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // ATTENTION: we must to null current neighbor connectivity out of reach because
                     // this creates potential disconnections within the mesh
                     //(*first)->Assign( face, static_cast<CELL<3>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the cells to one another
           for ( auto& it : cell_pairs ) {
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
                else if ( n_face_nbors > 2 ) {
                     // finding all possible combinations of surface elements
                     vector<uint32_t> sequence( n_face_nbors );
                     iota( sequence.begin(), sequence.end(), 0 ); // fill 0..n-1
                     const uint32_t            n_samples{2U};
                     deque<vector<uint32_t> >  combinations;
                     const size_t n_combinations = createUniqueCombinations( sequence, n_samples, combinations );
                     // finding the combination of surfaces or line elements with the smallest acute angle between them
                     map<double,uint32_t>  ordered_combinations;
                     for ( uint32_t i{0U}; i < n_combinations; ++i ) {
                          CELL<3U>* const ptr1 = (*next(it.second.begin(),combinations[i][0])).first;
                          CELL<3U>* const ptr2 = (*next(it.second.begin(),combinations[i][1])).first;
                          assert( ptr1 != nullptr );
                          assert( ptr2 != nullptr );
                          const double acute_angle = ( ptr1->IsSurface() && ptr2->IsSurface() ) ?
                                                       angleBetweenSurfaceCells( ptr1, ptr2 ) : angleBetweenLineCells( ptr1, ptr2 );
                          // ordering
                          ordered_combinations.insert( make_pair(acute_angle,i) );
                       }
                     // processing the combinations until there are no more pairs of cells left to process
                     for ( const auto& n : ordered_combinations ) {
                           // starting with the first combination with the smallest angle between the line elements
                           CELL<3U>* const ptr1 = (*next(it.second.begin(),combinations[n.second][0])).first;
                           CELL<3U>* const ptr2 = (*next(it.second.begin(),combinations[n.second][1])).first;
                           const uint32_t face_e1  = (*next(it.second.begin(),combinations[n.second][0])).second;
                           const uint32_t face_e2  = (*next(it.second.begin(),combinations[n.second][1])).second;
                           // checking whether cells have already been assigned a neighbor
                           if ( ptr1->Neighbor(face_e1) != nullptr || ptr2->Neighbor(face_e2) != nullptr )
                             continue;
                           ptr1->Assign( face_e1, ptr2 );
                           ptr2->Assign( face_e2, ptr1 );
                        }
                   }
                // else no assignments have to be made as there is no neighbor
             }
        }
 
      // for surface elements, faces or interfaces in a 2D model
      if constexpr ( dim == 2U ) {
           map<set<Node<2>*>,map<CELL<2>*,uint32_t> > cell_pairs;
           const auto                                 cellsEnd{last};
           
           // pairing the elements up in the search map
           while ( first != cellsEnd ) {
                assert( (*first) != nullptr );
                const uint32_t n_faces{ (*first)->Faces() };
                for ( uint32_t face{0u}; face < n_faces; ++face ) {
                     // trying to insert it into the map
                     auto it = cell_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<CELL<2>*,uint32_t>{make_pair(*first,face)} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<CELL<2>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the elements to one another
           for ( auto& it : cell_pairs ) {
                size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2U ) {
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



/**
   Line element manifolds exist in 3D and 2D.
*/
template<uint32_t dim>
template<template<uint32_t> class CELL>
void MeshManager<dim>::BuildLineConnectivity( typename vector<CELL<dim>*>::const_iterator first,
                                              typename vector<CELL<dim>*>::const_iterator last )
   {
      ErrorHandler& csmp_error( ErrorHandler::Instance() );
  
#ifndef NDEBUG
      if ( (detectDuplicateCells<dim,CELL>( first, last, true )) > 0 )
        csmp_error.Note( WARNING, "MeshManager<dim>::BuildLineConnectivity", "detected duplicate cells in iterator range");
#endif
  
      // this method applies only to line elements in 2 and 3D
      if constexpr ( dim != 1U ) {
           // creating search keys from the corner nodes of the element faces
           // corner-nodes      elements that share face and their face id
           map<Node<dim>*,map<CELL<dim>*,uint32_t> >  cell_pairs;
           const auto                                 cellsEnd{last};
           
           // pairing the elements up in the search map
           while ( first != cellsEnd ) {
                assert( (*first) != nullptr );
                const uint32_t n_faces{ (*first)->Faces() };
                for ( uint32_t face{0U}; face < n_faces; ++face ) {
                     // now there is only a single corner node corresponding to the opposite face of the line element
                     // (node 1 is at Face 0 and node 0 at Face 1 as for all simplex elements)
                     auto it = cell_pairs.insert( make_pair( (*first)->N( n_faces - face - 1U ),
                                                              map<CELL<dim>*,uint32_t>{make_pair((*first),face)} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // ATTENTION: we must to null current neighbor connectivity out of reach because
                     // this creates potential disconnections within the mesh
                     //(*first)->Assign( face, static_cast<CELL<3>*>(nullptr) );
                  }
                first++;
             }
 
           // processing the results, connecting the cells to one another
           for ( auto& it : cell_pairs ) {
                auto n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2U ) {
                     CELL<dim>* const ptr1  = (*it.second.begin()).first;
                     CELL<dim>* const ptr2  = (*it.second.rbegin()).first;
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is a manifold and two most suitable neighbors must be found
                else if ( n_face_nbors > 2U ) {
                     // finding all possible combinations of line elements
                     vector<uint32_t> sequence( n_face_nbors );
                     iota( sequence.begin(), sequence.end(), 0 ); // fill 0..n-1
                     const uint32_t           n_samples{2U};
                     deque<vector<uint32_t>>  combinations;
                     const size_t n_combinations = createUniqueCombinations( sequence, n_samples, combinations );
                     // finding the combination of lines with the smallest acute angle between them
                     map<double,uint32_t>  ordered_combinations;
                     for ( uint32_t i{0U}; i < n_combinations; ++i ) {
                          CELL<dim>* const ptr1 = (*next(it.second.begin(),combinations[i][0])).first;
                          CELL<dim>* const ptr2 = (*next(it.second.begin(),combinations[i][1])).first;
                          const double acute_angle = angleBetweenLineCells( ptr1, ptr2 );
                          // ordering
                          ordered_combinations.insert( make_pair(acute_angle,i) );
                       }
                     // processing the combinations until there are no more pairs of cells left to process
                     for ( const auto& n : ordered_combinations ) {
                           // starting with the first combination with the smallest angle between the line elements
                           CELL<dim>* const ptr1 = (*next(it.second.begin(),combinations[n.second][0])).first;
                           CELL<dim>* const ptr2 = (*next(it.second.begin(),combinations[n.second][1])).first;
                           const uint32_t face_e1  = (*next(it.second.begin(),combinations[n.second][0])).second;
                           const uint32_t face_e2  = (*next(it.second.begin(),combinations[n.second][1])).second;
                           // checking whether cells have already been assigned a neighbor
                           if ( ptr1->Neighbor(face_e1) != nullptr || ptr2->Neighbor(face_e2) != nullptr )
                             continue;
                           ptr1->Assign( face_e1, ptr2 );
                           ptr2->Assign( face_e2, ptr1 );
                        }
                  }
                // else no assignments have to be made as there is no neighbor
             }
        }
 
      // for line elements, faces or interfaces in a 1D model
      if constexpr ( dim == 1U ) {
           map<set<Node<1>*>,map<CELL<1>*,uint32_t> >  cell_pairs;
           const auto                                  cellsEnd{last};
           
           // pairing the elements up in the search map
           while ( first != cellsEnd ) {
                assert( (*first) != nullptr );
                const uint32_t n_faces{ (*first)->Faces() };
                for ( uint32_t face{0u}; face < n_faces; ++face ) {
                     // trying to insert it into the map
                     auto it = cell_pairs.insert( make_pair( (*first)->CornerNodesOfFace(face), map<CELL<1>*,uint32_t>{make_pair(*first,face)} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<CELL<1>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the elements to one another
           for ( auto& it : cell_pairs ) {
                auto n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2U ) {
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
template void MeshManager<1>::BuildLineConnectivity<Element>( typename vector<Element<1>*>::const_iterator,
                                                              typename vector<Element<1>*>::const_iterator );

template void MeshManager<3>::BuildLineConnectivity<Face>( typename vector<Face<3>*>::const_iterator,
                                                           typename vector<Face<3>*>::const_iterator );
template void MeshManager<2>::BuildLineConnectivity<Face>( typename vector<Face<2>*>::const_iterator,
                                                           typename vector<Face<2>*>::const_iterator );
template void MeshManager<1>::BuildLineConnectivity<Face>( typename vector<Face<1>*>::const_iterator,
                                                           typename vector<Face<1>*>::const_iterator );

template void MeshManager<3>::BuildLineConnectivity<InterFace>( typename vector<InterFace<3>*>::const_iterator,
                                                                typename vector<InterFace<3>*>::const_iterator );
template void MeshManager<2>::BuildLineConnectivity<InterFace>( typename vector<InterFace<2>*>::const_iterator,
                                                                typename vector<InterFace<2>*>::const_iterator );
template void MeshManager<1>::BuildLineConnectivity<InterFace>( typename vector<InterFace<1>*>::const_iterator,
                                                                typename vector<InterFace<1>*>::const_iterator );





/**
       Special method to obtain the connectivity of InterFace objects, which is tricker because InterFace nodes are multiplicated. Yet these node clusters are grouped into NodeManifold objects.
       Thus, rather than creating keys from Node pointers to match the faces of the InterFaces with one-another, this method uses the NodeManifold pointers.
       
       @note When InterFace nodes are the same on either side of the InterFace and there is no manifold object, the method casts the node pointers into Manfifolds to still permit a comparison.
       
       @note because InterFaces are lines in 2D and surfaces in 3D, the disambiguation in the case where there are multiple neighbors is dimension dependent.
       Correspondingly, there are constexpr compile time pathways that distinguish the proceduures for them.
       
       @attention method does not touch the connectivity of interfaces to their higher dimensional neighbors
       
       @author SKM
       @date 10/4/22
*/
template<uint32_t dim>
void MeshManager<dim>::BuildInterFaceConnectivity( typename vector<InterFace<dim>*>::const_iterator first,
                                                   typename vector<InterFace<dim>*>::const_iterator last )
   {
      // for InterFaces = surface elements in s 3D model
      // -----------------------------------------------
      if constexpr ( dim == 3U ) {
           // creating search keys from the corner nodes of the interface faces
           // corner-nodes      interfaces that share face and their face id
           map<set<NodeManifold<3U>*>,map<InterFace<3U>*,uint32_t> >  iface_pairs;

           // pairing the cells up in the search map
           const auto cellsEnd{last};
           while ( first != cellsEnd ) {
                assert( (*first) != nullptr );
                const size_t n_faces{ (*first)->Faces() };
                for ( uint32_t face{0u}; face < n_faces; ++face ) {
                     // getting node-manifold pointers from set of pointers to face corner nodes
                     set<NodeManifold<3U>*> manifold_ptrs;
                     for ( const auto& nit : (*first)->CornerNodesOfFace(face) ) {
                           if ( nit->IsManifold() ) manifold_ptrs.insert( nit->Manifold() );
                           else manifold_ptrs.insert( reinterpret_cast<NodeManifold<3U>*>(nit) );
                        }
                     auto it = iface_pairs.insert( make_pair( manifold_ptrs, map<InterFace<3U>*,uint32_t>{{*first,face}} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<InterFace<3U>*>(nullptr) );
                  }
                first++;
             }
             
           // processing the results, connecting the cells to one another
           for ( auto& it : iface_pairs ) {
                const size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2U ) {
                     InterFace<3U>* const ptr1  = (*it.second.begin()).first;
                     InterFace<3U>* const ptr2  = (*it.second.rbegin()).first;
                     assert( ptr1 != nullptr );
                     assert( ptr2 != nullptr );
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is an interface manifold and two most suitable surface neighbors must be found
                else if ( n_face_nbors > 2U ) { // TODO: test
                     // finding all possible combinations of surface elements
                     vector<uint32_t> sequence( n_face_nbors );
                     iota( sequence.begin(), sequence.end(), 0 ); // fill 0..n-1
                     const uint32_t n_samples{2U};
                     deque<vector<uint32_t> >  combinations;
                     const size_t n_combinations = createUniqueCombinations( sequence, n_samples, combinations );
                     // finding the combination of surfaces or line elements with the smallest acute angle between them
                     map<double,size_t>  ordered_combinations;
                     for ( auto i{0U}; i < n_combinations; ++i ) {
                          InterFace<3U>* const ptr1 = (*next(it.second.begin(),combinations[i][0])).first;
                          InterFace<3U>* const ptr2 = (*next(it.second.begin(),combinations[i][1])).first;
                          assert( ptr1 != nullptr );
                          assert( ptr2 != nullptr );
                          const double acute_angle = ( ptr1->IsSurface() && ptr2->IsSurface() ) ?
                                                       angleBetweenSurfaceCells( ptr1, ptr2 ) : angleBetweenLineCells( ptr1, ptr2 );
                          // ordering
                          ordered_combinations.insert( make_pair(acute_angle,i) );
                       }
                     // processing the combinations until there are no more pairs of cells left to process
                     for ( const auto& n : ordered_combinations ) {
                           // starting with the first combination with the smallest angle between the line elements
                           InterFace<3U>* const ptr1 = (*next(it.second.begin(),combinations[n.second][0])).first;
                           InterFace<3U>* const ptr2 = (*next(it.second.begin(),combinations[n.second][1])).first;
                           const uint32_t face_e1  = (*next(it.second.begin(),combinations[n.second][0])).second;
                           const uint32_t face_e2  = (*next(it.second.begin(),combinations[n.second][1])).second;
                           // checking whether cells have already been assigned a neighbor
                           if ( ptr1->Neighbor(face_e1) != nullptr || ptr2->Neighbor(face_e2) != nullptr )
                             continue;
                           ptr1->Assign( face_e1, ptr2 );
                           ptr2->Assign( face_e2, ptr1 );
                        }
                  }
                // else no assignments have to be made as there is no neighbor
             }
        }
 
      // for line elements, interfaces in a 2D model
      // -------------------------------------------
      if constexpr ( dim == 2U ) {
           map<NodeManifold<2U>*,map<InterFace<2U>*,uint32_t> > iface_pairs;
           const auto ifacesEnd{last};
           
           // pairing the elements up in the search map
           while ( first != ifacesEnd ) {
                assert( (*first) != nullptr );
                const auto n_faces{ (*first)->Faces() };
                for ( auto face{0U}; face < n_faces; ++face ) {
                     // getting node-manifold pointers: the first 2 nodes must correspond to the corner nodes of the Face
                     NodeManifold<2U>* mptr = ( (*first)->N(face)->IsManifold() ) ? (*first)->N(face)->Manifold() : reinterpret_cast<NodeManifold<2U>*>((*first)->N(face));
                     assert( mptr );
                    // creating the map entry
                     auto it = iface_pairs.insert( make_pair( mptr, map<InterFace<2U>*,uint32_t>{{*first,face}} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( it.second == false )
                       (*it.first).second.insert( make_pair( (*first), face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     (*first)->Assign( face, static_cast<InterFace<2U>*>(nullptr) );
                  }
                first++;
             }
 
           // processing the results, connecting the elements to one another
           for ( auto& it : iface_pairs ) {
                size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2U ) {
                     InterFace<2U>* const ptr1 = (*it.second.begin()).first;
                     InterFace<2U>* const ptr2 = (*it.second.rbegin()).first;
                     assert( ptr1 != nullptr );
                     assert( ptr2 != nullptr );
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is a manifold and two most suitable neighbors must be found
                else if ( n_face_nbors > 2U ) {
                     // finding all possible combinations of surface elements
                     vector<uint32_t> sequence( n_face_nbors );
                     iota( sequence.begin(), sequence.end(), 0 ); // fill 0..n-1
                     const uint32_t n_samples{2U};
                     deque<vector<uint32_t> >  combinations;
                     const size_t n_combinations = createUniqueCombinations( sequence, n_samples, combinations );
                     // finding the combination of surfaces with the smallest acute angle between them
                     map<double,uint32_t>  ordered_combinations;
                     for ( auto i{0U}; i < n_combinations; ++i ) {
                          InterFace<2U>* const ptr1 = (*next(it.second.begin(),combinations[i][0])).first;
                          InterFace<2U>* const ptr2 = (*next(it.second.begin(),combinations[i][1])).first;
                          const double acute_angle = angleBetweenLineCells( ptr1, ptr2 );
                          // ordering
                          ordered_combinations.insert( make_pair(acute_angle,i) );
                       }
                     // processing the combinations until there are no more pairs of cells left to process
                     for ( const auto& n : ordered_combinations ) {
                           // starting with the first combination with the smallest angle between the line elements
                           InterFace<2U>* const ptr1 = (*next(it.second.begin(),combinations[n.second][0])).first;
                           InterFace<2U>* const ptr2 = (*next(it.second.begin(),combinations[n.second][1])).first;
                           const uint32_t face_e1  = (*next(it.second.begin(),combinations[n.second][0])).second;
                           const uint32_t face_e2  = (*next(it.second.begin(),combinations[n.second][1])).second;
                           // checking whether cells have already been assigned a neighbor
                           if ( ptr1->Neighbor(face_e1) != nullptr || ptr2->Neighbor(face_e2) != nullptr )
                             continue;
                           ptr1->Assign( face_e1, ptr2 );
                           ptr2->Assign( face_e2, ptr1 );
                        }
                  }
                // else no assignments have to be made as there is no neighbor
             }
        }
 
 } // end BuildSurfaceElementConnectivity




/**
      Wholesale re-establishment of the connections between Elements, Faces and Interfaces, however, restricted to those of the same dimenisonality
      (volumetric Elements to volumetric Elements, surface to surface and line to line). Since lower-dimensional cells can be manifolds
      (multiple neighbors per face), these are disambiguated by chosing those cells as neigbors that lie in the same plane (or closest to)
      or define the same direction.
      
      Once the the connections between cells are re-established, the method reconnectes the nodes to parent elements.
      Using this information, the nodes are connected to their neighbor elements (avoiding connections between the nodes that form part the same manifold).
      
      @note The connection between cells is not affected by Regions, Boundaries or Splitboundaries. The Face objects on the outside
      of a model, for instance, are all interconnected with one another.
      
      @note Elements are not connected across SplitBoundaries, but these connections can be retrieved from corresponding InterFaces.
      
      @note The parent connectivity of nodes includes only Element objects, i.e., not Face or InterFace objects.
    
    @attention calls EraseNullPointerCells()
    
    @note After some diagnostics that establish the extent of mesh modification, the cell neighbor connectivity is rebuilt.
    
    @author SKM
    @date 12/10/21
*/
template<uint32_t dim>
void MeshManager<dim>::UpdateConnectivity()
 {
    cout <<"\n"<<"MeshManager<dim>::UpdateConnectivity: rebuilding cell connectivity for entire model..."<< endl;
    
    // 1. creating cell vectors needed by BuildConnectivity() methods and reconnecting equidimensional cells
    // -----------------------------------------------------------------------------------------------------
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
         BuildInterFaceConnectivity( iface_ptrs.begin(), iface_ptrs.end() );
      }
    
    // 2. (Re)-creating node connectivity to parent elements
    // -----------------------------------------------------
    // 2.1 counting the parent elements of each node
    unordered_map<Node<dim>*,set<Element<dim>*> >  parent_elmts_per_node;
    for ( auto& it : elements_ ) {
        assert( it.FE() );
        const auto nodes_end{ it.NodesEnd() };
        for ( auto nit = it.NodesBegin(); nit != nodes_end; ++nit ) {
             auto mit = parent_elmts_per_node.insert( make_pair( (*nit), set<Element<dim>*>{ &it } ) );
             if ( mit.second == false )
               (*mit.first).second.insert( &it );
          }
      }

    // 2.2 re-assigning the parent elements to nodes
    for ( auto& n : parent_elmts_per_node ) {
         const auto n_parents = static_cast<uint32_t>( n.second.size() );
         n.first->ReserveParentStorage( n_parents );
         // looping over the future parents
         for ( const auto& it : n.second ) {
           const uint32_t n_nodes{it->Nodes()};
           // assigning them to the node
           for ( uint32_t j{0U}; j<n_nodes; ++j )
             if ( n.first == it->N(j) ) {
                 it->N(j)->Assign( j, it );
                 break;
              }
           }
         assert( n.first->Parents() >= 1 );
      }
    
    // 2.3 Rebuilding the node connectivity
    for ( auto& nit : nodes_ ) nit.AssignNodeNeighbors();
    

    // 4. Update node manifolds
    // ------------------------
    // 4.1 (Re)-creating node connectivity to parent interfaces
    // --------------------------------------------------------
    unordered_map<Node<dim>*,set<pair<InterFace<dim>*,pair<uint32_t,INTERFACE_SIDE>>> >  parent_ifaces_per_node;
    for ( auto& it : interfaces_ ) {
      assert( it.FE() );
      const uint32_t nodes{ it.FE()->Nodes() };
      for ( uint32_t n{0U}; n < nodes; ++n ) {
        vector<INTERFACE_SIDE> sides{INSIDE,OUTSIDE};
        for ( INTERFACE_SIDE side : sides ){
          if( it.N(n,side)->IsManifold() ){
            set<pair<InterFace<dim>*,pair<uint32_t,INTERFACE_SIDE>>> trial_set{make_pair( &it, make_pair(n,side))};
            auto mit = parent_ifaces_per_node.insert( make_pair( it.N(n,side), trial_set )) ; //inserting NodeManifold - InterFace pair if it doesnt exist
            if (mit.second == false) //if insertion didnt happen because manifold already exists
              mit.first->second.insert( *trial_set.begin() ); //add new interface to already existing set of InterFaces for existing NodeManifold
          }//end of found manifold
        }//end of side iteration
      }//end of node loop
    }//end of interface loop


    // 5. re-assigning parent interface index pairs to the node manifold
    // -----------------------------------------------------------------
    for ( auto& n : parent_ifaces_per_node ) {
         assert(n.first->IsManifold());
         n.first->Manifold()->Assign(n.first, n.second); //assigning set of interfaces and index pairs to the manifold
      }

    cout <<"\t"<<"established cell-to-cell connectivity."<< endl;

 } // end UpdateConnectivity




/**
       Same as above but only for the Elements and Nodes in the given range.
*/
template<uint32_t dim>
void MeshManager<dim>::UpdateConnectivity(  typename vector<InterFace<dim>*>::const_iterator first,
                                            typename vector<InterFace<dim>*>::const_iterator last )
 {
    cout <<"\n"<<"MeshManager<dim>::UpdateConnectivity: rebuilding connectivity in the vicinity of supplied InterFaces..."<< endl;
    
    // 0. Finding all Elements sharing a node with the interface
    // ---------------------------------------------------------
    vector<Element<dim>*> element_ptrs;
    element_ptrs.reserve( static_cast<size_t>(distance(first,last) * 3) );
    vector<Node<dim>*> core_node_ptrs;
    core_node_ptrs.reserve( static_cast<size_t>(distance(first,last) * 4) );
    
    for ( auto it=first; it!=last; ++it ) {
         // disconnects inner from outer parent element
         assert( (*it)->InnerParent() != nullptr );
         if ( (*it)->InnerParent()->Neighbor( (*it)->InnerParentFaceID() ) )
           (*it)->InnerParent()->Neighbor( (*it)->InnerParentFaceID() )->UnassignNeighbor( (*it)->InnerParentFaceID() );
         assert( (*it)->OuterParent() != nullptr );
         if ( (*it)->OuterParent()->Neighbor( (*it)->OuterParentFaceID() ) )
           (*it)->OuterParent()->Neighbor( (*it)->OuterParentFaceID() )->UnassignNeighbor( (*it)->OuterParentFaceID() );
         // higher-dimensional neighbors and their neighbors
         // inside
         element_ptrs.push_back( (*it)->InnerParent() );
         for ( auto eit=(*it)->InnerParent()->NeighborsBegin(); eit!=(*it)->InnerParent()->NeighborsEnd(); ++eit )
           if ( (*eit) != nullptr )
             element_ptrs.push_back( (*eit) );
         // outside
         element_ptrs.push_back( (*it)->OuterParent() );
         for ( auto eit=(*it)->OuterParent()->NeighborsBegin(); eit!=(*it)->OuterParent()->NeighborsEnd(); ++eit )
           if ( (*eit) != nullptr )
             element_ptrs.push_back( (*eit) );
         // nodes
         for ( auto nit=(*it)->NodesBegin(); nit!=(*it)->NodesEnd(); ++nit ) {
              assert( (*nit) != nullptr );
              core_node_ptrs.push_back( (*nit) );
           }
      }
    
    // eliminating duplicates
    // elements
    sort( element_ptrs.begin(), element_ptrs.end() );
    element_ptrs.erase( unique( element_ptrs.begin(), element_ptrs.end() ), element_ptrs.end() );
    // core nodes
    sort( core_node_ptrs.begin(), core_node_ptrs.end() );
    core_node_ptrs.erase( unique( core_node_ptrs.begin(), core_node_ptrs.end() ), core_node_ptrs.end() );
    
    // distinguishing halo nodes from core nodes
    /* not needed yet
    vector<Node<dim>*> node_ptrs; // halo_node_ptrs.reserve( distance(first,last) );
    for ( auto& eit : element_ptrs )
      for ( auto nit=eit->NodesBegin(); nit!=eit->NodesEnd(); ++nit )
        node_ptrs.push_back( (*nit) );
    // eliminating duplicates
    sort( node_ptrs.begin(), node_ptrs.end() );
    node_ptrs.erase( unique( node_ptrs.begin(), node_ptrs.end() ), node_ptrs.end() );
    // halo nodes
    vector<Node<dim>*> halo_node_ptrs; // halo_node_ptrs.reserve( distance(first,last) );
    set_difference( node_ptrs.begin(), node_ptrs.end(),
                   core_node_ptrs.begin(), core_node_ptrs.end(),
                   back_inserter(halo_node_ptrs) );
    */
    
    // 1. creating Element vector needed by BuildConnectivity() method, reconnecting equidimensional cells
    // ---------------------------------------------------------------------------------------------------
    // (this method preserves connections to cells outside the current range of cells)
    BuildInterFaceConnectivity( first, last );
    BuildConnectivity<csmp::Element>( element_ptrs.begin(), element_ptrs.end() );

    // 2. (Re)-creating node connectivity to parent elements
    // -----------------------------------------------------
    // 2.1 counting the parent elements of each node (considering all elements connected to core nodes guarantees that none are missed
    map<Node<dim>*,set<Element<dim>*> >  parent_elmts_per_node;
    for ( auto& it : element_ptrs ) {
        assert( it->FE() );
        const auto nodes_end{ it->NodesEnd() };
        for ( auto nit = it->NodesBegin(); nit != nodes_end; ++nit ) {
             auto mit = parent_elmts_per_node.insert( make_pair( (*nit), set<Element<dim>*>{ it } ) );
             if ( mit.second == false )
               (*mit.first).second.insert( it );
          }
      }

    // 2.2 re-assigning the parent elements to nodes
    for ( auto& n : parent_elmts_per_node ) {
         const auto n_parents = static_cast<uint32_t>( n.second.size() );
         n.first->ReserveParentStorage( n_parents );
         // looping over the future parents
         for ( const auto& it : n.second ) {
           const uint32_t n_nodes{it->Nodes()};
           // assigning them to the node
           for ( uint32_t j{0U}; j<n_nodes; ++j )
             if ( n.first == it->N(j) ) {
                 it->N(j)->Assign( j, it );
                 break;
              }
           }
         assert( n.first->Parents() >= 1 );
      }
    
    // 2.3 Rebuilding the node connectivity
    for ( auto& nit : core_node_ptrs ) nit->AssignNodeNeighbors();
//    for ( auto& nit : nodes_ ) nit.AssignNodeNeighbors();
    

    // 4. Update node manifolds
    // ------------------------
    // 4.1 (Re)-creating node connectivity to parent interfaces
    // --------------------------------------------------------
    map<Node<dim>*,set<pair<InterFace<dim>*,pair<uint32_t,INTERFACE_SIDE>>> >  parent_ifaces_per_node;
    for ( auto it=first; it!=last; ++it ) {
      assert( (*it)->FE() );
      const uint32_t nodes{ (*it)->FE()->Nodes() };
      for ( uint32_t n{0U}; n < nodes; ++n ) {
        vector<INTERFACE_SIDE> sides{INSIDE,OUTSIDE};
        for ( INTERFACE_SIDE side : sides ){
          if( (*it)->N(n,side)->IsManifold() ){
            set<pair<InterFace<dim>*,pair<uint32_t,INTERFACE_SIDE>>> trial_set{make_pair( (*it), make_pair(n,side))};
            auto mit = parent_ifaces_per_node.insert( make_pair( (*it)->N(n,side), trial_set )) ; //inserting NodeManifold - InterFace pair if it doesnt exist
            if (mit.second == false) //if insertion didnt happen because manifold already exists
              mit.first->second.insert( *trial_set.begin() ); //add new interface to already existing set of InterFaces for existing NodeManifold
          }//end of found manifold
        }//end of side iteration
      }//end of node loop
    }//end of interface loop


    // 5. re-assigning parent interface index pairs to the node manifold
    // -----------------------------------------------------------------
    for ( auto& n : parent_ifaces_per_node ) {
         assert(n.first->IsManifold());
         n.first->Manifold()->Assign(n.first, n.second); //assigning set of interfaces and index pairs to the manifold
      }

    cout <<"\t"<<"(re)established connectivity."<< endl;

 } // end UpdateConnectivity( near interfaces )




/**
       Reconnects Faces and adjacent Elements,  rebuilds node-parents for Face nodes.
       TODO: part 2 seems to break element connectivity
*/
template<uint32_t dim>
void MeshManager<dim>::UpdateConnectivity(  typename vector<Face<dim>*>::const_iterator first,
                                            typename vector<Face<dim>*>::const_iterator last )
 {
    cout <<"\n"<<"MeshManager<dim>::UpdateConnectivity: rebuilding the connectivity in the vicinity of supplied Faces..."<< endl;
    
    // 0. Finding all Elements sharing a node with the interface
    // ---------------------------------------------------------
    vector<Element<dim>*> element_ptrs;
    element_ptrs.reserve( static_cast<size_t>(distance(first,last) * 3) );
    vector<Node<dim>*> core_node_ptrs;
    core_node_ptrs.reserve( static_cast<size_t>(distance(first,last) * 4) );
    
    for ( auto it=first; it!=last; ++it ) {
         // (re)connects inner from outer parent element
         assert( (*it)->InnerParent() != nullptr );
         assert( (*it)->OuterParent() != nullptr );
         (*it)->InnerParent()->Neighbor( (*it)->InnerParentFaceID() )->Assign( (*it)->InnerParentFaceID(), (*it)->OuterParent() );
         (*it)->OuterParent()->Neighbor( (*it)->OuterParentFaceID() )->Assign( (*it)->OuterParentFaceID(), (*it)->InnerParent() );
         // higher-dimensional neighbors and their neighbors
         // inside
         element_ptrs.push_back( (*it)->InnerParent() );
         for ( auto eit=(*it)->InnerParent()->NeighborsBegin(); eit!=(*it)->InnerParent()->NeighborsEnd(); ++eit )
           if ( (*eit) != nullptr )
             element_ptrs.push_back( (*eit) );
         // outside
         element_ptrs.push_back( (*it)->OuterParent() );
         for ( auto eit=(*it)->OuterParent()->NeighborsBegin(); eit!=(*it)->OuterParent()->NeighborsEnd(); ++eit )
           if ( (*eit) != nullptr )
             element_ptrs.push_back( (*eit) );
         // nodes
         for ( auto nit=(*it)->NodesBegin(); nit!=(*it)->NodesEnd(); ++nit ) {
              assert( (*nit) != nullptr );
              core_node_ptrs.push_back( (*nit) );
           }
      }
    
    // eliminating duplicates
    // elements
    sort( element_ptrs.begin(), element_ptrs.end() );
    element_ptrs.erase( unique( element_ptrs.begin(), element_ptrs.end() ), element_ptrs.end() );
    // core nodes
    sort( core_node_ptrs.begin(), core_node_ptrs.end() );
    core_node_ptrs.erase( unique( core_node_ptrs.begin(), core_node_ptrs.end() ), core_node_ptrs.end() );


    // 1. creating Element vector needed by BuildConnectivity() method, reconnecting equidimensional cells
    // ---------------------------------------------------------------------------------------------------
    BuildConnectivity<Face>( first, last );
    BuildConnectivity<csmp::Element>( element_ptrs.begin(), element_ptrs.end() );

    // 2. (Re)-creating node connectivity to parent elements
    // -----------------------------------------------------
    // 2.1 counting the parent elements of each node
    map<Node<dim>*,set<Element<dim>*> >  parent_elmts_per_node;
    for ( auto& it : element_ptrs ) {
        assert( it->FE() );
        const auto nodes_end{ it->NodesEnd() };
        for ( auto nit = it->NodesBegin(); nit != nodes_end; ++nit ) {
             auto mit = parent_elmts_per_node.insert( make_pair( (*nit), set<Element<dim>*>{ it } ) );
             if ( mit.second == false )
               (*mit.first).second.insert( it );
          }
      }

    // 2.2 re-assigning the parent elements to nodes
    for ( auto& n : parent_elmts_per_node ) {
         const auto n_parents = static_cast<uint32_t>( n.second.size() );
         n.first->ReserveParentStorage( n_parents );
         // looping over the future parents
         for ( const auto& it : n.second ) {
           const uint32_t n_nodes{it->Nodes()};
           // assigning them to the node
           for ( uint32_t j{0U}; j<n_nodes; ++j )
             if ( n.first == it->N(j) ) {
                 it->N(j)->Assign( j, it );
                 break;
              }
           }
         assert( n.first->Parents() >= 1 );
      }

    // 2.3 Rebuilding the node connectivity
    for ( auto& nit : core_node_ptrs ) nit->AssignNodeNeighbors();

    cout <<"\t"<<"(re)established connectivity."<< endl;

 } // end UpdateConnectivity( near interfaces )





///  for nodes attached to elements in the supplied element range, the parent and the neighbor connectivity is reconstructed from scratch
// TODO: this forgets the elements that only have a single node on the deleted region
template<uint32_t dim>
void MeshManager<dim>::ConnectNodesToParentsAndNeighbors( typename vector<Element<dim>*>::iterator first,
                                                          typename vector<Element<dim>*>::iterator last )
 {
    assert( distance(first,last) >= 1U );
 
    set<Node<dim>*>  nodes_to_update;
        
    // 1. Connecting nodes to their parent elements
    // --------------------------------------------
    // counting the parent elements of each node
    map<Node<dim>*,set<Element<dim>*> >  parent_elmts_per_node;
    while ( first != last ) {
        if ( (*first) == nullptr ) continue;
        assert( (*first)->FE() );
        const auto nodes_end{ (*first)->NodesEnd() };
        for ( auto nit = (*first)->NodesBegin(); nit != nodes_end; ++nit ) {
             nodes_to_update.insert( (*nit) );
             auto mit = parent_elmts_per_node.insert( make_pair( (*nit), set<Element<dim>*>{ (*first) } ) );
             if ( mit.second == false )
               (*mit.first).second.insert( (*first) );
          }
        first++;
      }

    // 2. Assigning parent elements to the nodes
    // -----------------------------------------
    for ( auto& n : parent_elmts_per_node ) {
         const auto n_parents = static_cast<uint32_t>( n.second.size() );
         n.first->ReserveParentStorage( n_parents );
         // looping over the future parents
         for ( const auto& it : n.second ) {
           const auto n_nodes{ it->Nodes() };
           // assigning them to the node
           for ( uint32_t j{0U}; j<n_nodes; ++j )
             if ( n.first == it->N(j) ) {
                 it->N(j)->Assign( j, it );
                 break;
              }
           }
         assert( n.first->Parents() >= 1 );
      }
    
    // 3. Connecting the nodes to their node neighbors
    // -----------------------------------------------
    for ( auto& nit : nodes_to_update ) nit->AssignNodeNeighbors();
    
    // TODO: update potential node manifolds here as well
    
 } // end ConnectNodesToParentsAndNeighbors




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
  size_t n{ 0ul };
  for_each( nodes_.begin(), nodes_.end(), [&n]( Node<dim>& nd ) { nd.Idx( n++ ); } );

  n = 0U; // resetting the counter
  for_each( elements_.begin(), elements_.end(), [&n]( Element<dim>& e ) { e.Idx( n++ ); } );

  if ( !in_a_single_sequence ) n = 0U;
  for_each( faces_.begin(), faces_.end(), [&n]( Face<dim>& f ) { f.Idx( n++ ); } );

  if ( !in_a_single_sequence ) n = 0U;
  for_each( interfaces_.begin(), interfaces_.end(), [&n]( InterFace<dim>& i ) { i.Idx( n++ ); } );

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
       csmp_error.Note( WARNING, "MeshManager<dim>::OutputMeshTo", "Mesh Manager does not currently store a mesh; no output");
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
      const size_t interfaceExtras( 1U ); // 1 entry for potential high dim element

      deque<uint32_t>  nodes_per_element;
      deque<uint32_t>  neighbors_per_element;
      deque<int8_t>    csmp_fem_types;

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
            nodes_per_element.push_back( f.Nodes() );
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
    vset.BFlag( n.Idx(), n.AtBoundary() );
    
  // geometry flags
  vset.ResizeBREP_Flags( /* nodes */ );
  for ( const auto& n : nodes_ )
    vset.BREP_Flag( n.Idx(), n.Attribute() );
  
  // 'pelmt' was already set above


  // 3. adding 'plist' connectivity list and 'pmtrl'
  // -----------------------------------------------
  vector<int32_t>  pmtrl( elements_.size(), 0 );
  size_t           eidx{0};
  // elements
  for ( const auto& e : elements_ ) {
      const auto n_nodes{e.Nodes()};
      for ( uint32_t j{0U}; j<n_nodes; ++j )
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
        for ( uint32_t j{0U}; j<n_nodes; ++j )
          vset.Plist( eidx, j, (f.N( j )->Idx()) );
        ++eidx;
      }

  // interfaces
  if ( !interfaces_.empty() )
    for ( const auto& f : interfaces_ ) {
        const auto n_nodes{f.FE()->Nodes()};
        for ( uint32_t j{0U}; j<n_nodes; ++j )
          vset.Plist( eidx, j, (f.N( j, INSIDE )->Idx()) );
        for ( uint32_t j{0U}; j<n_nodes; ++j )
          vset.Plist( eidx, n_nodes + j, (f.N( j, OUTSIDE )->Idx()) );
        ++eidx;
      }


  // 4. adding 'pfverts' neighbors per element list
  // ----------------------------------------------   
  // 'pfverts' elements
  eidx = 0;
  for ( const auto& e : elements_ ) {
      const auto neighbors{e.Neighbors()};
      assert( neighbors <= 6 );
      for ( uint32_t j{0U}; j<neighbors; ++j ) {
            const Element<dim>* const ptr = e.Neighbor(j);
            if ( ptr ) {
                assert( e.IsVolume() == ptr->IsVolume() && e.IsSurface() == ptr->IsSurface() && e.IsLine() == ptr->IsLine() );
                assert( ptr->Idx() < elements_.size() );
                vset.Pfvert( eidx, j, static_cast<int64_t>(ptr->Idx()) );
              }
            else {
                int64_t bflag = e.AtBoundary(j);
                if ( bflag == 0 ) {
                   e.Out();
                   csmp_error.Note( ERROR, "MeshManager<dim>::OutputMeshTo", "'pfvert' atBoundary evaluated to zero");
                   e.AtBoundary(j);
                }
                vset.Pfvert( eidx, j, bflag );
              }
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
      for ( uint32_t j{0U}; j<neighbors; ++j ) {
           // if the neighbor exists
           if ( f.Neighbor(j) != nullptr )
             vset.Pfvert( eidx, j, static_cast<int64_t>(f.Neighbor(j)->Idx()) );
           else {
                // the face is either located on a model boundary or an internal boundary
                bool all_nodes_at_external_boundary{true};
                for ( const auto& fnit : f.CornerNodesOfFace(j) )
                  if ( fnit->AtBoundary() == NOT || fnit->AtBoundary() == INTERNAL ) {
                       all_nodes_at_external_boundary = false;
                       break;
                    }
                if ( all_nodes_at_external_boundary ) vset.Pfvert( eidx, j, IRREGULAR );
                else vset.Pfvert( eidx, j, INTERNAL );
             }
        }
      // higher-dimensional neighbors second
      // inner neighbor
      assert( f.InnerParent()->IsEquidimensional() );
      assert( f.InnerParent()->Idx() < Elements() );
      vset.Pfvert( eidx, neighbors, static_cast<int64_t>(f.InnerParent()->Idx()) );
      // outer neighbor
      if ( f.OuterParent() != nullptr ) {
          assert( f.OuterParent()->IsEquidimensional() );
          assert( f.OuterParent()->Idx() < Elements() );
          vset.Pfvert( eidx, neighbors + 1U, static_cast<int64_t>(f.OuterParent()->Idx()) );
        }
      else {
          // getting the placement of the inner element at the boundary (negative number), if it cannot be determined, INTERNAL boundary is set
          // TODO: perhaps this requires an extra corner case in atBoundary()
          vset.Pfvert( eidx, neighbors + 1U, std::min(atBoundary( f.InnerParent(), f.InnerParentFaceID() ),INTERNAL) );
        }
      // face id's converted to
      // adding the local numbers of the faces that the Face is collocated with if any
      vset.Pfvert( eidx, neighbors + 2U, f.InnerParentFaceID() );
      // if there is no outer element, the face idx will initialised with UNSPECIFIED
      if ( f.OuterParentFaceID() == numeric_limits<uint32_t>::max() )
        vset.Pfvert( eidx, neighbors + 3U, UNSPECIFIED );
      else vset.Pfvert( eidx, neighbors + 3U, f.OuterParentFaceID() );
      ++eidx;
   }

  // 'pfverts' interfaces (which must always have two higher-dimensional neighbors)
  // ------------------------------------------------------------------------------
  for ( const auto& f : interfaces_ ) {
    // 1. equidimensional neighbors (=other interfaces) first
    //    they are written in the order in which they are stored in the interface
    //    To simplify things there is always an entry for the intervening element even if there is none.
    const auto neighbors{ f.Neighbors() };
    for ( uint32_t j{0U}; j<neighbors; ++j ) {
         if ( f.Neighbor(j) != nullptr )
           vset.Pfvert( eidx, j, static_cast<int64_t>(f.Neighbor(j)->Idx()) );
         else
           vset.Pfvert( eidx, j, INTERNAL );
      }
    // 2. inner and outer higher-dimensional neighbors (2 entries)
    //   (they must always exist because SplitBoundaries are internal model boundaries)
    assert( f.InnerParent() != nullptr );
    assert( f.OuterParent() != nullptr );
    assert( f.InnerParent()->Idx() < elements_.size() );
    assert( f.OuterParent()->Idx() < elements_.size() );
    vset.Pfvert( eidx, neighbors,      static_cast<int64_t>(f.InnerParent()->Idx()) );
    vset.Pfvert( eidx, neighbors + 1U, static_cast<int64_t>(f.OuterParent()->Idx()) );
      
    // 3. storing local number of face of the inner and outer elements that the InterFace is connected to (2 entries)
    assert( f.InnerParentFaceID() < f.InnerParent()->Faces() );
    assert( f.OuterParentFaceID() < f.OuterParent()->Faces() );
    vset.Pfvert( eidx, neighbors + 2U, f.InnerParentFaceID() );
    vset.Pfvert( eidx, neighbors + 3U, f.OuterParentFaceID() );
    
    // 4. storing number of intervening element or nullptr identifier (one entry)
    if ( f.HasInterveningElement() ) {
         assert( f.InterveningElement()->Idx() < elements_.size() );
         vset.Pfvert( eidx, neighbors + 4U, static_cast<int64_t>(f.InterveningElement()->Idx()) );
      }
    else
      vset.Pfvert( eidx, neighbors + 4U, INTERNAL );
      
    ++eidx;
  }

  // 'pmanifolds' supporting interfaces: initialising VData::manifold_container
  // --------------------------------------------------------------------------
  if ( !interfaces_.empty() ) {
       assert( node_manifold_manager_ );
       VData::manifoldContainer  node_manifolds;
       node_manifolds.reserve( node_manifold_manager_->Manifolds() );
       for ( auto nmf=node_manifold_manager_->ManifoldsBegin(); nmf!=node_manifold_manager_->ManifoldsEnd(); ++nmf )
         node_manifolds.push_back( (*nmf).Data() );
       vset.AddNodeManifolds( node_manifolds.begin(), node_manifolds.end() );
    }

  cout << "\nMeshManager<" << dim << ">::OutputMeshTo: MeshManager successfully output to VSet..." << endl;

} // end OutputMeshTo( VSet )


/* DEBUGGING
cout <<"\n\n"<<"testing node numbering:\n";
for ( const auto& n : nodes_ ) cout << n.Idx() <<" ";
cout <<"\n\n"<<"testing element numbering:\n";
for ( const auto& n : elements_ ) cout << n.Idx() <<" ";
cout <<"\n\n"<<"testing face numbering:\n";
for ( const auto& n : faces_ ) cout << n.Idx() <<" ";
cout <<"\n\n"<<"testing interface numbering:\n";
for ( const auto& n : interfaces_ ) cout << n.Idx() <<" ";
cout << endl;
*/




/**

Writes scalar, vector, tensor, array and flagged array property data from MeshManager to VSet, including 'pmtrl' material property information.

Process is applied to all finite volumes and elements, their integration points and nodes,
but not for any of the variables stored on the model, region, boundaries or splitboundaries.
A method for these entities is provided separately.

To store the properties in the VSet, they are first written to PropertyData objects.
These are then added to the VSet property storage.

@attention  a continuous numbering of elements, faces, interfaces and nodes has to be
established before this method is called, for instance, using AssignUniqueNumbers().

@attention property values are output if the corresponding dataset does not contain any NAN values.

@note the values of region, boundary or splitboundary properties are stored together with the these in respective binary files

@attention in the VSet,  variables are identified only by their (unique) names. Very litlle extra information is contained in the
PropertyData class. The property database is therefore essential to retrieve all other variable related information.

@author SKM 5/5/2016

*/
template<uint32_t dim>
void MeshManager<dim>::OutputStoredVariablesTo( const PropertyDatabase<dim>& database, 
                                                VSet<dim>& vset ) const
{
  ErrorHandler&		    csmp_error( ErrorHandler::Instance() );
  map<string, Index>  properties;
  bool                detected_first_NaN_value{ false };
  
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

  database.ListVariables( NODE, properties );
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
    bool variable_contains_NaN_values{ false };

    switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : nodes_ ) {
                it.Read( (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of node variable not recognized." );
      }
    // storing the data in the VSet
    if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
    else {
        if ( !detected_first_NaN_value ) {
              csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                               (*pit).first, "contained NaN values and was therefore not stored in VSet." );
             detected_first_NaN_value = true;
          }
        else cout <<" '"<< (*pit).first <<"'";
      }
  }

  // -------------------------------------------------
  // -------------------------------------------------
  // element properties (including integration points)
  // -------------------------------------------------
  // -------------------------------------------------
  database.ListVariables( ELEMENT, properties );
  const size_t elements( elements_.size() );
  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    // creating the property storage
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    // for the given property type
    const size_t flag_capacity( elements * (*pit).second.flagDepth );
    const size_t data_capacity( elements * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );
    bool variable_contains_NaN_values{ false };

// DEBUGGING
#ifndef NDEBUG
//set<string> NaN_scalars;
#endif

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( const auto& it : elements_ ) {
              it.Read( (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
//                     NaN_scalars.insert( (*pit).first + to_string(it.Idx()) );
                     break;
                  }
              pushBack( data, value );
            }
          }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( const auto& it : elements_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( const auto& it : elements_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : elements_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : elements_ ) {
              it.Read( (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
       break;
      default:
        csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                           (*pit).first, "type of element variable not recognized." );
    }
    // storing the data in the VSet
    if ( !variable_contains_NaN_values ) {
          vset.AddData( (*pit).first.c_str(), data );
      }
    else {
        if ( !detected_first_NaN_value ) {
              csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                               (*pit).first, "contained NaN values and was therefore not stored in VSet." );
             detected_first_NaN_value = true;
          }
        else cout <<" '"<< (*pit).first <<"'";
      }
  }

  // element integration point properties
  // ------------------------------------
  if ( database.ListVariables( ELEMENT_INTEGRATION_POINT, properties ) > 0 ) {
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
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : elements_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
     break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                          (*pit).first, "type of element integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }


  // element sector integration point properties
  // -------------------------------------------
  if ( database.ListVariables( SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_sector_ips( (*elements_.begin()).IntegrationPointsPerSector() );
    const size_t sectors_per_element( (*elements_.begin()).Sectors() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : elements_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
                const auto sectors{ it.Sectors() };
                for ( uint32_t i{0U}; i<sectors; ++i ) {
                  const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
                  for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    if ( value.Has_NaN_Values() ) {
                         variable_contains_NaN_values = true;
                         break;
                      }
                    pushBack( data, value );
                  }
                }
              }
            }
         break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of element sector integraton point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }


  // element facet integration point properties
  // ------------------------------------------
  if ( database.ListVariables( FACET_INTEGRATION_POINT, properties ) > 0 ) {
    assert( Elements() > 0 );
    const size_t elmt_facet_ips( (*elements_.begin()).IntegrationPointsPerFacet() );
    const size_t facets_per_element( (*elements_.begin()).Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : elements_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : elements_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : elements_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }


  // ----------------------------------------------
  // ----------------------------------------------
  // face properties (including integration points)
  // ----------------------------------------------
  // ----------------------------------------------
  database.ListVariables( FACE, properties );
  const size_t n_faces( Faces() );

  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    const size_t flag_capacity( n_faces * (*pit).second.flagDepth );
    const size_t data_capacity( n_faces  * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );
    bool variable_contains_NaN_values{ false };

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( const auto& it : faces_ ) {
              it.Read( (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( const auto& it : faces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( const auto& it : faces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : faces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : faces_ ) {
              it.Read( (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        break;
      default:
        csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                           (*pit).first, "type of face variable not recognized." );
    }
    // storing the data in the VSet
    if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
    else {
        if ( !detected_first_NaN_value ) {
              csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                               (*pit).first, "contained NaN values and was therefore not stored in VSet." );
             detected_first_NaN_value = true;
          }
        else cout <<" '"<< (*pit).first <<"'";
      }
  }


  // face integration point properties
  // ---------------------------------
  if ( database.ListVariables( FACE_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
    assert( Faces() > 0U );
    const size_t face_ips( (*faces_.begin()).IntegrationPoints() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      // for the given property type
      const size_t flag_capacity( n_faces * face_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_faces * face_ips  * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : faces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
     break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of face integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }


  // face sector integration point properties
  // ----------------------------------------
  if ( database.ListVariables( FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
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
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : faces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const uint32_t ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of face-sector integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }


  // face facet integration point properties
  // ---------------------------------------
  if ( database.ListVariables( FACE_FACET_INTEGRATION_POINT, properties ) > 0 && Faces() > 0 ) {
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
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : faces_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : faces_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : faces_ ) {
            const auto facets{ it.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const uint32_t ips_per_facet( it.IntegrationPointsPerFacet() );
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of face facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }


  // ---------------------------------------------------
  // ---------------------------------------------------
  // interface properties (including integration points)
  // ---------------------------------------------------
  // ---------------------------------------------------
  database.ListVariables( INTER_FACE, properties );
  const size_t n_interfaces( Interfaces() );

  for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
  {
    PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
    const size_t flag_capacity( n_interfaces * (*pit).second.flagDepth );
    const size_t data_capacity( n_interfaces  * (*pit).second.dataDepth );
    data.Reserve( flag_capacity, data_capacity );
    bool variable_contains_NaN_values{ false };

    switch ( (*pit).second.type )
    {
      case SCALAR: {
        ScalarVariable value;
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
        break;
      case VECTOR: {
        VectorVariable<dim> value;
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
        break;
      case TENSOR: {
        TensorVariable<dim> value;
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
        break;
      case ARRAY: {
        ArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
       break;
      case FLAGGEDARRAY: {
        FlaggedArrayVariable value( (*pit).second.dataDepth );
        for ( const auto& it : interfaces_ ) {
          it.Read( (*pit).second, value );
          if ( value.Has_NaN_Values() ) {
               variable_contains_NaN_values = true;
               break;
            }
          pushBack( data, value );
        }
      }
   break;
      default:
        csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                           (*pit).first, "type of interface variable not recognized." );
    }
    // storing the data in the VSet
    if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
    else {
        if ( !detected_first_NaN_value ) {
              csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                               (*pit).first, "contained NaN values and was therefore not stored in VSet." );
             detected_first_NaN_value = true;
          }
        else cout <<" '"<< (*pit).first <<"'";
      }
  }

  // interface integration point properties
  // --------------------------------------
  if ( database.ListVariables( INTER_FACE_INTEGRATION_POINT, properties ) > 0 && Interfaces() > 0 ) {
    const size_t interface_ips{ (*interfaces_.begin()).IntegrationPoints() };

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * interface_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * interface_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : interfaces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const uint32_t integration_points( it.IntegrationPoints() );
            for ( uint32_t i{0U}; i<integration_points; ++i ) {
              it.Read( i, (*pit).second, value );
              if ( value.Has_NaN_Values() ) {
                   variable_contains_NaN_values = true;
                   break;
                }
              pushBack( data, value );
            }
          }
        }
     break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of interface integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }

  // interface sector integration point properties
  // ---------------------------------------------
  if ( database.ListVariables( INTER_FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 && Interfaces() > 0 ) {
    assert( Interfaces() > 0 );
    const size_t interface_sector_ips{ (*interfaces_.begin()).IntegrationPointsPerSector() };
    const size_t sectors_per_interface{ (*interfaces_.begin()).Sectors() };

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
      {
        case SCALAR: {
          ScalarVariable value;
          for ( const auto& it : interfaces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case VECTOR: {
          VectorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case TENSOR: {
          TensorVariable<dim> value;
          for ( const auto& it : interfaces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
          break;
        case ARRAY: {
          ArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
         break;
        case FLAGGEDARRAY: {
          FlaggedArrayVariable value( (*pit).second.dataDepth );
          for ( const auto& it : interfaces_ ) {
            const auto sectors{ it.Sectors() };
            for ( uint32_t i{0U}; i<sectors; ++i ) {
              const auto ips_per_sector( it.IntegrationPointsPerSector() );
              for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                it.Read( i, j, (*pit).second, value );
                if ( value.Has_NaN_Values() ) {
                     variable_contains_NaN_values = true;
                     break;
                  }
                pushBack( data, value );
              }
            }
          }
        }
     break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of interface sector integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }

  // interface facet integration point properties
  // --------------------------------------------
  if ( database.ListVariables( INTER_FACE_FACET_INTEGRATION_POINT, properties ) > 0 && Interfaces() > 0 ) {
    assert( Interfaces() > 0 );
    const size_t interface_facet_ips( (*interfaces_.begin()).IntegrationPointsPerFacet() );
    const size_t facets_per_interface( (*interfaces_.begin()).Facets() );

    for ( auto pit = properties.begin(); pit != properties.end(); ++pit )
    {
      PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
      const size_t flag_capacity( n_interfaces * facets_per_interface * interface_facet_ips * (*pit).second.flagDepth );
      const size_t data_capacity( n_interfaces * facets_per_interface * interface_facet_ips * (*pit).second.dataDepth );
      data.Reserve( flag_capacity, data_capacity );
      bool variable_contains_NaN_values{ false };

      switch ( (*pit).second.type )
        {
          case SCALAR: {
                ScalarVariable value;
                for ( const auto& it : interfaces_ ) {
                  const auto facets{ it.Facets() };
                  for ( uint32_t i{0U}; i<facets; ++i ) {
                    const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                    for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                      it.Read( i, j, (*pit).second, value );
                      if ( value.Has_NaN_Values() ) {
                           variable_contains_NaN_values = true;
                           break;
                        }
                      pushBack( data, value );
                    }
                  }
                }
              }
          break;
        case VECTOR: {
              VectorVariable<dim> value;
              for ( const auto& it : interfaces_ ) {
                const auto facets{ it.Facets() };
                for ( uint32_t i{0U}; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    if ( value.Has_NaN_Values() ) {
                         variable_contains_NaN_values = true;
                         break;
                      }
                    pushBack( data, value );
                  }
                }
              }
            }
          break;
        case TENSOR: {
              TensorVariable<dim> value;
              for ( const auto& it : interfaces_ ) {
                const auto facets{ it.Facets() };
                for ( uint32_t i{0U}; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    if ( value.Has_NaN_Values() ) {
                         variable_contains_NaN_values = true;
                         break;
                      }
                    pushBack( data, value );
                  }
                }
              }
            }
          break;
        case ARRAY: {
              ArrayVariable value( (*pit).second.dataDepth );
              for ( const auto& it : interfaces_ ) {
                const auto facets{ it.Facets() };
                for ( uint32_t i{0U}; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    if ( value.Has_NaN_Values() ) {
                         variable_contains_NaN_values = true;
                         break;
                      }
                    pushBack( data, value );
                  }
                }
              }
            }
          break;
        case FLAGGEDARRAY: {
              FlaggedArrayVariable value( (*pit).second.dataDepth );
              for ( const auto& it : interfaces_ ) {
                const auto facets{ it.Facets() };
                for ( uint32_t i{0U}; i<facets; ++i ) {
                  const auto ips_per_facet( it.IntegrationPointsPerFacet() );
                  for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                    it.Read( i, j, (*pit).second, value );
                    if ( value.Has_NaN_Values() ) {
                         variable_contains_NaN_values = true;
                         break;
                      }
                    pushBack( data, value );
                  }
                }
              }
            }
          break;
        default:
          csmp_error.Note( ERROR, "MeshManager<dim>::OutputStoredVariablesTo:",
                             (*pit).first, "type of interface facet integration point variable not recognized." );
      }
      // storing the data in the VSet
      if ( !variable_contains_NaN_values ) vset.AddData( (*pit).first.c_str(), data );
      else {
          if ( !detected_first_NaN_value ) {
                csmp_error.Note( WARNING, "MeshManager<dim>::OutputStoredVariablesTo:",
                                 (*pit).first, "contained NaN values and was therefore not stored in VSet." );
               detected_first_NaN_value = true;
            }
          else cout <<" '"<< (*pit).first <<"'";
        }
    }
  }
  if ( detected_first_NaN_value ) cout <<"  are other variables that were not stored because they contained NaN values."<< endl;

} // end OutputStoredVariablesTo



/// helper function that searches property sets for variables with specific placement
static bool containsPropertiesWithPlacement( PLACEMENT place,
                                             map<string,PropertyData>::const_iterator pit,
                                             const map<string,PropertyData>::const_iterator end ) {
   while ( pit != end  ) {
        if ( (*pit).second.Placement() == place ) return true;
        pit++;
     }
   return false;
}


  /**
  Read and assign property values written exactly by OutputStoredVariablesTo()
  
  @attention It is assumed that no erasures have occured in the MeshManager before the properties are assigned; else the order is no longer correct;

  @ test SKM 28/6/2016
  */
template<uint32_t dim>
void MeshManager<dim>::InputStoredVariablesFrom( const PropertyDatabase<dim>& database, const VSet<dim>& vset )
{
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
  
  if ( vset.Vertices() != nodes_.size() &&
       containsPropertiesWithPlacement( NODE, vset.PropertyValuesBegin(), vset.PropertyValuesEnd() ) ) {
       csmp_error.Note( WARNING, "MeshManager<dim>::InputStoredVariablesFrom",
                               "mismatch between number of vset and model nodes; NODE property input aborted." );
       return;
    }
  if ( vset.Elements() != elements_.size() &&
       containsPropertiesWithPlacement( ELEMENT, vset.PropertyValuesBegin(), vset.PropertyValuesEnd() ) ) {
       csmp_error.Note( WARNING, "MeshManager<dim>::InputStoredVariablesFrom",
                               "mismatch between number of vset and model elements; ELEMENT property input aborted." );
       return;
    }
  if ( vset.Faces() != faces_.size() &&
       containsPropertiesWithPlacement( FACE, vset.PropertyValuesBegin(), vset.PropertyValuesEnd() ) ) {
       csmp_error.Note( WARNING, "MeshManager<dim>::InputStoredVariablesFrom",
                               "mismatch between number of vset and model faces; FACE property input aborted." );
       return;
    }
  if ( vset.Interfaces() != interfaces_.size() &&
       containsPropertiesWithPlacement( INTER_FACE, vset.PropertyValuesBegin(), vset.PropertyValuesEnd() ) ) {
       csmp_error.Note( WARNING, "MeshManager<dim>::InputStoredVariablesFrom",
                               "mismatch between number of vset and model interfaces; INTERFACE property input aborted." );
       return;
    }

  // 'pmtrl' flags for the elements
  // ==============================
  assert( elements_.size() == vset.Elements() );
  typename plf::colony<Element<dim>>::iterator it = elements_.begin();
  for ( auto pit=vset.PmtrlBegin(); pit!=vset.PmtrlEnd(); ++pit, ++it )
    (*it).Material_ID( (*pit) );

  // -----------------------------------------
  // assigns properties, NODE properties first
  // -----------------------------------------
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
          csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
            csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
              for ( uint32_t i{0U}; i<integration_points; ++i ) {
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
              for ( uint32_t i{0U}; i<integration_points; ++i ) {
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
              for ( uint32_t i{0U}; i<integration_points; ++i ) {
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
              for ( uint32_t i{0U}; i<integration_points; ++i ) {
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
              for ( uint32_t i{0U}; i<integration_points; ++i ) {
                    read( (*pit).second, entry, value );
                    e.Store( i, key, value );
                    entry++;
                  }
                }
              }
            break;
          default:
            csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
              const auto sectors{ e.Sectors() };
              for ( uint32_t i{0U}; i<sectors; ++i ) {
                const auto ips_per_sector{ e.IntegrationPointsPerSector() };
                for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
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
              const auto sectors{ e.Sectors() };
              for ( uint32_t i{0U}; i<sectors; ++i ) {
                const auto ips_per_sector{ e.IntegrationPointsPerSector() };
                for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
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
              const auto sectors{ e.Sectors() };
              for ( uint32_t i{0U}; i<sectors; ++i ) {
                const auto ips_per_sector{ e.IntegrationPointsPerSector() };
                for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
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
              const auto sectors{ e.Sectors() };
              for ( uint32_t i{0U}; i<sectors; ++i ) {
                const auto ips_per_sector{ e.IntegrationPointsPerSector() };
                for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
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
              const auto sectors{ e.Sectors() };
              for ( uint32_t i{0U}; i<sectors; ++i ) {
                const auto ips_per_sector{ e.IntegrationPointsPerSector() };
                for ( uint32_t j{0U}; j<ips_per_sector; ++j ) {
                      read( (*pit).second, entry, value );
                      e.Store( i, j, key, value );
                      entry++;
                    }
                  }
                }
              }
            break;
          default:
            csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
            const auto facets{ e.Facets() };
            for ( uint32_t i{0U}; i<facets; ++i ) {
              const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
              for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
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
              const auto facets{ e.Facets() };
              for ( uint32_t i{0U}; i<facets; ++i ) {
                const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
                for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
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
              const auto facets{ e.Facets() };
              for ( uint32_t i{0U}; i<facets; ++i ) {
                const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
                for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
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
              const auto facets{ e.Facets() };
              for ( uint32_t i{0U}; i<facets; ++i ) {
                const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
                for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
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
              const auto facets{ e.Facets() };
              for ( uint32_t i{0U}; i<facets; ++i ) {
                const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
                for ( uint32_t j{0U}; j<ips_per_facet; ++j ) {
                  read( (*pit).second, entry, value );
                  e.Store( i, j, key, value );
                  entry++;
                }
              }
            }
          }
       break;
      default:
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
          csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
        break;
      default:
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
        break;
      default:
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
                  read( (*pit).second, entry, value );
                  e.Store( i, j, key, value );
                  entry++;
                }
              }
            }
          }
        break;
      default:
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
    assert( (*pit).second.Size() / key.dataDepth == Interfaces() );

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
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
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
          for ( auto i{0U}; i<integration_points; ++i ) {
            read( (*pit).second, entry, value );
            e.Store( i, key, value );
            entry++;
          }
        }
      }
   break;
      default:
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
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
          const auto sectors{ e.Sectors() };
          for ( auto i{0U}; i<sectors; ++i ) {
            const auto ips_per_sector{ e.IntegrationPointsPerSector() };
            for ( auto j{0U}; j<ips_per_sector; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
   break;
      default:
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
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
          const auto facets{ e.Facets() };
          for ( auto i{0U}; i<facets; ++i ) {
            const auto ips_per_facet{ e.IntegrationPointsPerFacet() };
            for ( auto j{0U}; j<ips_per_facet; ++j ) {
              read( (*pit).second, entry, value );
              e.Store( i, j, key, value );
              entry++;
            }
          }
        }
      }
   break;
      default:
        csmp_error.Note( ERROR, "Region<dim>::OutputVariableTo:",
                           (*pit).first, "type of interface facet integration point variable not recognized." );
    }
  }

} // end InputStoredVariablesFrom




// ============================================================================================
//
//           DIAGNOSTICS
//
// ============================================================================================


template<uint32_t dim>
set<uint32_t> MeshManager<dim>::OrderOfShapeFunctions() const
 {
    if ( elements_.empty() )
      throw csmp::Exception( ERROR, "MeshManager<dim>::OrderOfShapeFunctions", "'elements_' container is empty." );

    // order of the interpolation functions of the elements, faces and interfaces used in the mesh
    set<uint32_t>  shape_function_orders;
    
    for ( const auto& it : elements_ )
      shape_function_orders.insert( it.FE()->Interpolation() );

    if ( !faces_.empty() )
      for ( const auto& it : faces_ )
        shape_function_orders.insert( it.FE()->Interpolation() );

    if ( !interfaces_.empty() )
      for ( const auto& it : interfaces_ )
        shape_function_orders.insert( it.FE()->Interpolation() );

    return shape_function_orders;
 }
 




/**
    Performs a node-to-node, breadth-first traversal to identify whether the model consists  of disconnected mesh patches.
    
    @note use sparingly! - this test can be quite memory intensive and slow for large models.
    
 */
template<uint32_t dim>
bool  MeshManager<dim>::IsContiguous() const
 {
    if ( elements_.empty() )
      throw csmp::Exception( ERROR, "MeshManager<dim>::IsContiguous", "'elements_' container is empty." );
      
    // using a floodfill on the model-dimensional elements to check for contiguity
    vector<const Element<dim>*>  equi_dim_elmts;
    equi_dim_elmts.reserve( Elements() );
    for ( const auto& it : elements_ )
     if ( it.IsEquidimensional() )
       equi_dim_elmts.push_back( &it );
       
    // performing the floodfill
    set<Element<dim>*> elements_contiguous_subset;
    floodFill( const_cast<Element<dim>*>(*equi_dim_elmts.begin()), elements_contiguous_subset );
 
    if ( elements_contiguous_subset.size() < equi_dim_elmts.size() ) return false;
    return true;

 } // end IsContiguous


/* COLLECTION OF FAILED TRIALS

// produces the same results
   set<Node<dim>*> node_pointers;
   findInterconnectedNodeCluster<dim>( const_cast<Node<dim>*>(&(*nodes_.begin())), node_pointers );

    // creating a set of corner node-pointer pairs representing the element edges
    fail = always true!
    set<pair<const Node<dim>*,const Node<dim>*>> validEdges;
    vector<uint32_t>  snids;
    
    for ( const auto& element : elements_ ) {
         for ( uint32_t segm_id{0u}; segm_id<element.Segments(); ++segm_id )
           element.FE()->NodesOfSegment( segm_id, snids );
         // ignoring any mid-side nodes
         validEdges.insert( make_pair( element.N(snids[0]), element.N(snids[1]) ) );
      }
 
    auto interconnected_nodes = findInterconnectedNodes( nodes_, validEdges );

   if ( interconnected_nodes < Nodes() ) return false;
   return true;
   
    // Perform union-find
  UnionFind<const Element<dim>*> unionFind;
  for ( const auto& eit : elements_ ) {
      const Element<dim>* e = &eit;
      const uint32_t neighbors( e->Neighbors() );
      for ( uint32_t i = 0U; i<neighbors; i++ )
        if ( e->Neighbor( i ) != nullptr )
          unionFind.SameComponent( e, e->Neighbor( i ) );
    }
  deque<pair<size_t,const Element<dim>*>> components;
  unionFind.Components( components );
  
  if ( components.size() <= 1 ) return true;

*/



template<uint32_t dim>
bool  MeshManager<dim>::HybridElementMesh() const
{
  return hybrid_element_mesh_;
}



  ///  Reports nodes that do not belong to any parent elements, faces or interfaces
template<uint32_t dim>
size_t MeshManager<dim>::OrphanNodes() const {
     return MeshManager<dim>::OrphanNodeVector().size();
  }


template<uint32_t dim>
vector<const Node<dim>*> MeshManager<dim>::OrphanNodeVector() const
  {
     // renumbering the Node objects consecutively and creating a corresponding integer vector
     size_t counter{0};
     for ( const auto& node : nodes_ )  node.Idx( counter++ );
     vector<long> all_nodes( nodes_.size(), UNSPECIFIED );
     iota( all_nodes.begin(), all_nodes.end(), 0 );
 
     plf::colony<size_t>  connected_nodes;
 
     for ( const auto& elmt : elements_ )
       for ( auto nit=elmt.NodesBegin(); nit!=elmt.NodesEnd(); ++nit )
         connected_nodes.insert( (*nit)->Idx() );
     
     // better than for vector, these are efficient member functions
     connected_nodes.sort();
     connected_nodes.unique();
 
     if ( !faces_.empty() ) {
          for ( const auto& face : faces_ )
           for ( auto nit=face.NodesBegin(); nit!=face.NodesEnd(); ++nit )
             connected_nodes.insert( (*nit)->Idx() );
         
          connected_nodes.sort();
          connected_nodes.unique();
       }
     if ( !interfaces_.empty() ) {
          for ( const auto& iface : interfaces_ )
           for ( auto nit=iface.NodesBegin(); nit!=iface.NodesEnd(); ++nit )
             connected_nodes.insert( (*nit)->Idx() );
         
          connected_nodes.sort();
          connected_nodes.unique();
       }

     // assumes that there are no duplicates
     if ( connected_nodes.size() == nodes_.size() ) return vector<const Node<dim>*>( /* empty vec */ );
  
     vector<long> nodes_missing;
     set_difference( all_nodes.begin(), all_nodes.end(),
                     connected_nodes.begin(), connected_nodes.end(),
                     back_inserter(nodes_missing) );
     
     vector<const Node<dim>*> orphan_nodes;
     
     if ( !nodes_missing.empty() ) {
         for ( const auto& node : nodes_missing )
           orphan_nodes.push_back( &(*next(nodes_.begin(),node)) );
       }
     
     return orphan_nodes;
  
  } // end OrphanNodeVector





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
      if      ( !with_line_elements    && it.IsLine() )		with_line_elements = true;
      else if ( !with_surface_elements && it.IsSurface() )	with_surface_elements = true;
      else if ( !with_volume_elements  && it.IsVolume() )	with_volume_elements = true;
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
      for ( auto i{0U}; i<eit.Faces(); i++ )
        // if the face is at a model boundary
        if ( eit.Neighbor( i ) == nullptr )
          {
            // boundary nodes
            assert( eit.FE() != NULL );
            for ( const auto& j : eit.FE()->NodesOfFace(i) )
              boundary_nodes.insert( eit.N(j) );
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
        for ( uint32_t i{0U}; i<eit.Nodes(); ++i ) {
          assert( eit.N( i ) != nullptr );
          highest_dim_elmt_nodes.insert( eit.N( i ) );
        }
        // if the element has faces that lie on the region boundary
        // it is considered a boudary element
        long  nbors_that_belong_to_group( eit.Neighbors() );
        for ( auto i{0U}; i<eit.Faces(); ++i )
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
        for ( uint32_t i{0U}; i<e->Nodes(); ++i )
          if ( !highest_dim_elmt_nodes.count( e->N( i ) ) )
            exterior_nodes++;

        // if individual nodes stick out the parent element sticks out as well.
        if ( exterior_nodes == e->Nodes() )
          lesser_dim_elmts_detached.insert( e );
      }

    if ( !lesser_dim_elmts_detached.empty() ) {
      csmp_error.Note( WARNING, "MeshManager::CheckElementConnectivity:",
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
         for ( auto i{0U}; i<e.Neighbors(); ++i )
           if ( e.Neighbor(i) == nullptr ) {
                bound = parseBoundary( atBoundary(&e,i) );
                break;
             }
         cout << "\nElement ID: " << e.Idx() <<" ("<< parseFiniteElementType(e.FE_Type());
         cout <<"), Boundary flag: " << bound << endl;
         cout << "Member Nodes: " << endl;
         for ( uint32_t i{0U}; i < e.Nodes(); i++ )
           cout << e.N( i )->Idx() << "\t";
         cout << "\nNeighbor elements: " << endl;
         for ( uint32_t i{0U}; i < e.Neighbors(); i++ )
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
        for ( uint32_t i{0U}; i < f.Nodes(); i++ )
          cout << f.N( i )->Idx() << "\t";
        cout << "\nNeighbor faces: " << endl;
        for ( uint32_t i{0U}; i < f.Neighbors(); i++ )
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
        for ( uint32_t i{0U}; i < f.FE()->Nodes(); i++ )
          cout << f.N( i, INSIDE )->Idx() << "\t";
        for ( uint32_t i{0U}; i < f.FE()->Nodes(); i++ )
          cout << f.N( i, OUTSIDE )->Idx() << "\t";
        cout << "\nNeighbor faces: " << endl;
        for ( uint32_t i{0U}; i < f.Neighbors(); i++ )
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
      for ( uint32_t i = 0u; i < n.Parents(); i++ )
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
       if ( fvm_manager_ ) fvm_manager_->Stencil(fit)->Out();
    }
  cout << endl;
  
  // node manifolds
  if ( node_manifold_manager_ != nullptr ) {
      cout << endl << endl;
      cout << "NODE MANIFOLDS: " << endl;
      node_manifold_manager_->Out();
    }
    
} // end Out



template class MeshManager<1U>;
template class MeshManager<2U>;
template class MeshManager<3U>;


// NON-MEMBER FUNCTIONS

/// finds cells that have the same nodes and reports their numbers
template<uint32_t dim, template<uint32_t> class CELL>
size_t detectDuplicateCells( typename plf::colony<CELL<dim>>::const_iterator first,
                             typename plf::colony<CELL<dim>>::const_iterator last,
                             bool verbose )
 {
    map<set<Node<dim>*>,set<const CELL<dim>*> > potential_duplicates;
    size_t  n_duplicates{0ul};
    size_t  cell_counter{0ul};
    
    while( first != last ) {
         // creating cell keys from their node pointers
         set<Node<dim>*> node_set;
         const auto n_nodes{ first->Nodes() };
         if constexpr( is_same<CELL<dim>,InterFace<dim>>::value ) {
              for ( uint32_t i{0U}; i<n_nodes/2u; i++ ) {
                   node_set.insert( first->N(i) );
                }
           }
         else if constexpr( is_same<CELL<dim>,Element<dim>>::value ||
                            is_same<CELL<dim>,Face<dim>>::value )  {
              for ( uint32_t i{0U}; i<n_nodes; i++ ) {
                   node_set.insert( first->N(i) );
                }
           }
         // recording the cells
         auto it = potential_duplicates.insert( make_pair( node_set, set<const CELL<dim>*>{ &(*first) } ) );
         // if there is a cell with the same nodes but a different pointer, it is recorded
         if ( it.second == false ) {
              auto cit =(*it.first).second.insert( &(*first) );
              if ( verbose && cit.second == false ) {
                  cout <<"\n"<<"detectDuplicateCells: input range contains multiple copies of:";
                  first->Out();
                }
              n_duplicates++;
           }
         first++;
         cell_counter++;
      }
      
    // printing the duplicate cells if any
    if ( n_duplicates  > 0U && verbose ) {
         cout <<"\n\n"<<"detectDuplicateCells: found "<< n_duplicates <<" cells sharing all nodes in input range:";
         for ( auto pd : potential_duplicates )
           if ( pd.second.size() > 1U )
             (*pd.second.begin())->Out();
      }
  
    return n_duplicates;
    
 } // end detectDuplicateCells

template size_t detectDuplicateCells<3,Element>( plf::colony<Element<3>>::const_iterator, plf::colony<Element<3>>::const_iterator, bool );
template size_t detectDuplicateCells<2,Element>( plf::colony<Element<2>>::const_iterator, plf::colony<Element<2>>::const_iterator, bool );
template size_t detectDuplicateCells<1,Element>( plf::colony<Element<1>>::const_iterator, plf::colony<Element<1>>::const_iterator, bool );

template size_t detectDuplicateCells<3,Face>( plf::colony<Face<3>>::const_iterator, plf::colony<Face<3>>::const_iterator, bool );
template size_t detectDuplicateCells<2,Face>( plf::colony<Face<2>>::const_iterator, plf::colony<Face<2>>::const_iterator, bool );
template size_t detectDuplicateCells<1,Face>( plf::colony<Face<1>>::const_iterator, plf::colony<Face<1>>::const_iterator, bool );

template size_t detectDuplicateCells<3,InterFace>( plf::colony<InterFace<3>>::const_iterator, plf::colony<InterFace<3>>::const_iterator, bool );
template size_t detectDuplicateCells<2,InterFace>( plf::colony<InterFace<2>>::const_iterator, plf::colony<InterFace<2>>::const_iterator, bool );
template size_t detectDuplicateCells<1,InterFace>( plf::colony<InterFace<1>>::const_iterator, plf::colony<InterFace<1>>::const_iterator, bool );



/**

Method supports the deletion a set of finite elements from your mesh, locally updating the connectivity of the remaining elements,
while it tries to minimize the number of operations by focusing only on the immediate neighbors of the elements being deleted.

Step-by-Step Algorithm

1.	Mark the Elements for Deletion: A unique negative region ID of -999  is assigned to each element that will be deleted. This attribute allows to identify deleted elements without affecting other mesh operations.
	
2.	Identify the Halo Neighbors of the Deleted Elements:
		For each element in the set to be deleted, iterate over its neighbors (from the std::vector of neighbor pointers).
		For each face, if the neighboring element on that face is not marked for deletion, mark it as a halo element for updating or collect a pointer to it in a separate container as it needs to have its connectivity adjusted.
  
3.	Update the Connectivity of Halo Elements:
		For each halo neighbor identified in the previous step, iterate over its list of neighbors (the neighbor pointers for each face).
		For each face with a neighbor that is marked for deletion, set the corresponding neighbor pointer to nullptr because this connection will become invalid.
  
@todo If the face with the deleted neighbor belongs to a lower dimensional element, it may be a manfold that is potentially  connected to other remaining element(s).
In this case, replace the deleted neighbor that will be deleted with another adjacent lower-dimensional element sharing its face nodes.
This requires searching because there is no knowledge stored about face-node sharing alternative neighbors.

Now that these tasks have been completed one can:

4.	Delete the Elements from the Mesh:
		Now that the border neighbors have been updated, you can safely delete the elements in the set marked for deletion.
  
5.	Clean-Up Neighbor Pointers:
		Remove any redundant or invalid pointers from the std::vector of neighbor pointers in each element to ensure the integrity of your mesh data structure.

@attention function assumes that all elements in the iterator range have the same dimension, i.e. are lines, surfaces or volumes

@section Application

call method with template arguments as follows

@code
          updateHaloCellConnectivity<3,Element>(first,last);
@endcode

@author SKM
@date 30/10/2024

*/
template<uint32_t dim>
void updateHaloElementConnectivity( typename vector<Element<dim>*>::iterator first,
                                    typename vector<Element<dim>*>::iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( distance(first,last) == 0 ) {
         csmp_error.Note( ERROR, "updateHaloElementConnectivity", "input element iterator range is empty");
         return;
      }
    
    const int32_t delete_mark{-999};
    size_t        issues{0ul};

    // ======================================================================
    // PROCESSING ELEMENTS THAT HAVED THE SAME SPATIAL DIMENSION AS THE MODEL
    // ======================================================================
    if ( (*first)->IsEquidimensional() ) {
        // 1. mark elements (using region ID), and remove them from the parent-element vectors of its nodes
        // (NB: changing the region ID does not matter because these elements will be deleted anyway)
        for_each( first, last, [ &csmp_error, &issues ](Element<dim>* cptr) {
                                     if ( cptr ) {
                                          // verifying that iterator points to an equidimensional element
                                          if ( !cptr->IsEquidimensional() ) {
                                               cout <<"\n\t"<<"Element<"<< dim <<">"<< cptr->Idx() <<": "<< parseAbbreviated_FE_Type(cptr->FE_Type());
                                               csmp_error.Note( WARNING, "updateHaloElementConnectivity", "element pointer points to lower-dimensional element");
                                               issues++;
                                            }
                                          // marking the elements via their region ID
                                          cptr->Region_ID(delete_mark);
                                       }
                                     else  {
                                         csmp_error.Note( WARNING, "updateHaloElementConnectivity", "input iterator contains nullptr");
                                         issues++;
                                       }
                                 } );

        // 2. collect the halo cells so that their neighbor connectivity can later be updated efficiently
        unordered_set<csmp::Element<dim>*> halo_cells;
        for_each( first, last, [ &halo_cells ](Element<dim>* it) {
                                    if (it) {
                                      const auto n_nbors = it->Neighbors();
                                      for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
                                        if ( it->Neighbor(eidx) && it->Neighbor(eidx)->Region_ID() != delete_mark )
                                          halo_cells.insert( it->Neighbor(eidx) );
                                    }
                                 } );
        
        if ( halo_cells.empty() ) {
             csmp_error.Note( INFO, "updateHaloElementConnectivity",
                             "equidimensional input mesh has no halo elements; no element-neighbor updates were made.");
             return;
          }
        
        // 3. update the connectivity of the halo elements
        // (here extra checks are performed in DEBUG mode to determine whether they have other nullptr neighbors like at model boundary)
        for ( auto& it : halo_cells ) {
             assert( it );
             const auto n_nbors = it->Neighbors();
             for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
               // if the (face) neighbor of this object is going to be deleted, set it to nullptr
               if ( it->Neighbor(eidx) && it->Neighbor(eidx)->Region_ID() == delete_mark )
                 {
                    // removing the element from the parent-element vectors of face nodes on the domain perimeter
                    for ( const auto& fnid : it->FE()->NodesOfFace(eidx) )
                      it->N(fnid)->Unassign( it->Neighbor(eidx) );
                    // remove neighbor connection
                    it->UnassignNeighbor(eidx);
                 }
          }
         if ( issues > 0 ) {
             csmp_error.Note( ERROR, "updateHaloElementConnectivity", "update failed; expect follow on errors.");
           }
         return;
         
      } // end equidimensional elements

    
    // ========================================================
    // PROCESSING LOWER_DIMENSIONAL ELEMENTS
    // ========================================================
    // (in this case the parent-element vector of all nodes need updating because these nodes are shared with equidimensional elements)
    const bool domain_of_surf_elmts = (*first)->IsSurface();
    const bool domain_of_line_elmts = (*first)->IsLine();
    cout <<"\n"<<"updateHaloElementConnectivity: provided input iterator range contains (lower-dimensional) ";
    if ( domain_of_surf_elmts ) cout <<"surface";
    else if ( domain_of_line_elmts ) cout <<"line";
    cout <<" elements.";
    cout <<" This method does not reconnect manifold neighbors of these if there are any."<< endl;
    
    // 1. mark elements (using region ID), and remove them from the parent-element vectors of its nodes
    // (NB: changing the region ID does not matter because these elements will be deleted anyway)
    for_each( first, last, [ &csmp_error, &issues, &domain_of_surf_elmts, &domain_of_line_elmts ](Element<dim>* cptr) {
                                 if ( cptr ) {
                                      // checking element dimensionality
                                      if constexpr( dim == 3 ) {
                                            if ( domain_of_surf_elmts && !cptr->IsSurface() ) {
                                                 cout <<"\n\t"<<"Element<"<< dim <<">"<< cptr->Idx() <<": "<< parseAbbreviated_FE_Type(cptr->FE_Type());
                                                 csmp_error.Note( WARNING, "updateHaloElementConnectivity", "wrong element dimension");
                                                 issues++;
                                               }
                                            else if ( domain_of_line_elmts  && !cptr->IsLine() ) {
                                                 cout <<"\n\t"<<"Element<"<< dim <<">"<< cptr->Idx() <<": "<< parseAbbreviated_FE_Type(cptr->FE_Type());
                                                 csmp_error.Note( WARNING, "updateHaloElementConnectivity", "wrong element dimension");
                                                 issues++;
                                               }
                                         }
                                      else if constexpr( dim == 2 ) {
                                            if ( !cptr->IsLine() ) {
                                                 cout <<"\n\t"<<"Element<"<< dim <<">"<< cptr->Idx() <<": "<< parseAbbreviated_FE_Type(cptr->FE_Type());
                                                 csmp_error.Note( WARNING, "updateHaloElementConnectivity", "wrong element dimension");
                                                 issues++;
                                               }
                                         }
                                      // removing Element from the parent-element vectors of its nodes
                                      for ( auto nit=cptr->NodesBegin(); nit!=cptr->NodesEnd(); ++nit )
                                         (*nit)->Unassign( cptr );
                                      // marking the elements via their region ID
                                      cptr->Region_ID(delete_mark);
                                   }
                                 else  {
                                     csmp_error.Note( WARNING, "updateHaloElementConnectivity", "input iterator contains nullptr");
                                     issues++;
                                   }
                             } );

    // 2. collect the halo cells so that their neighbor connectivity can later be updated efficiently
    unordered_set<csmp::Element<dim>*> halo_cells2;
    for_each( first, last, [ &halo_cells2 ](Element<dim>* it) {
                                if (it) {
                                  const auto n_nbors = it->Neighbors();
                                  for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
                                    if ( it->Neighbor(eidx) && it->Neighbor(eidx)->Region_ID() != delete_mark )
                                      halo_cells2.insert( it->Neighbor(eidx) );
                                }
                             } );
    
        if ( halo_cells2.empty() ) {
             csmp_error.Note( INFO, "updateHaloElementConnectivity",
                             "lower-dimensional input mesh has no halo elements; no element-neighbor updates were made.");
             return;
          }
    
    // 3. update the connectivity of the halo elements
    // (here extra checks are performed in DEBUG mode to determine whether they have other nullptr neighbors like at model boundary)
    for ( auto& it : halo_cells2 ) {
         assert( it );
         const auto n_nbors = it->Neighbors();
         for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
           // if the (face) neighbor of this object is going to be deleted, set it to nullptr
           if ( it->Neighbor(eidx) && it->Neighbor(eidx)->Region_ID() == delete_mark )
             {
                // remove neighbor connection
                it->UnassignNeighbor(eidx);
             }
      }

 } // end updateHaloElementConnectivity

template void updateHaloElementConnectivity<3>( vector<Element<3>*>::iterator, vector<Element<3>*>::iterator );
template void updateHaloElementConnectivity<2>( vector<Element<2>*>::iterator, vector<Element<2>*>::iterator );
template void updateHaloElementConnectivity<1>( vector<Element<1>*>::iterator, vector<Element<1>*>::iterator );






/**
        For Face and InterFace objects. Therefore does not deal with intervening elements, but only disconnects objects to be deleted
        from equal ones forming their halo, i.e., sharing faces with them on the Boundary or SplitBoundary perimeter.
*/
template<uint32_t dim, template<uint32_t> class CELL>
void updateHaloCellConnectivity( typename vector<CELL<dim>*>::iterator first,
                                 typename vector<CELL<dim>*>::iterator last )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( distance(first,last) == 0 ) {
         csmp_error.Note( ERROR, "updateHaloCellConnectivity", "input cell range is empty");
         return;
      }
    
    const size_t delete_mark{ numeric_limits<size_t>::max() };
 
    // 1. mark elements (using region ID), change does not matter because these elements will be deleted
    for_each( first, last, [](CELL<dim>* cptr) {
                                 if ( cptr ) {
                                      // marking the elements via their region ID
                                      cptr->Idx(delete_mark);
                                   }
                                 else cerr <<"\n"<<"updateHaloElementConnectivity: found 'nullptr' in input iterator range.";
                             } );

    // 2. collect the halo cells so that their neighbor connectivity can later be updated efficiently
    set<CELL<dim>*> halo_cells;
    for_each( first, last, [ &halo_cells ](CELL<dim>* cptr) {
                                const auto n_nbors = cptr->Neighbors();
                                for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
                                  if ( cptr->Neighbor(eidx) && cptr->Neighbor(eidx)->Idx() != delete_mark )
                                    halo_cells.insert( cptr->Neighbor(eidx) );
                             } );
    
    if ( halo_cells.empty() ) return;
    
    // 3. update the connectivity of the halo elements
    // (here extra checks are performed in DEBUG mode to determine whether they have other nullptr neighbors like at model boundary)
    for ( auto& it : halo_cells ) {
         assert( it );
         const auto n_nbors = it->Neighbors();
         for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
           // if the (face) neighbor of this object is going to be deleted, set it to nullptr
           if ( it->Neighbor(eidx) && it->Neighbor(eidx)->Idx() == delete_mark )
             it->UnassignNeighbor(eidx);
      }

 } // end updateHaloCellConnectivity

template void updateHaloCellConnectivity<3,Face>( vector<Face<3>*>::iterator, vector<Face<3>*>::iterator );
template void updateHaloCellConnectivity<2,Face>( vector<Face<2>*>::iterator, vector<Face<2>*>::iterator );
template void updateHaloCellConnectivity<1,Face>( vector<Face<1>*>::iterator, vector<Face<1>*>::iterator );

template void updateHaloCellConnectivity<3,InterFace>( vector<InterFace<3>*>::iterator, vector<InterFace<3>*>::iterator );
template void updateHaloCellConnectivity<2,InterFace>( vector<InterFace<2>*>::iterator, vector<InterFace<2>*>::iterator );
template void updateHaloCellConnectivity<1,InterFace>( vector<InterFace<1>*>::iterator, vector<InterFace<1>*>::iterator );




// COLONY ITERATOR VERSION
/*
template<uint32_t dim, template<uint32_t> class CELL>
void updateHaloCellConnectivity( typename plf::colony<CELL<dim>>::iterator first,
                                 typename plf::colony<CELL<dim>>::iterator last )
 {
    // 1. mark elements (using region ID), change does not matter because these elements will be deleted
    for_each( first, last, [](CELL<dim>& obj) { obj.Region_ID(-999); } );

    // 2. collect the halo cells so that their neighbor connectivity can later be updated efficiently
    set<CELL<dim>&> halo_cells;
    for_each( first, last, [ &halo_cells ](CELL<dim>& obj) {
                                const auto n_nbors = obj.Neighbors();
                                for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
                                  if ( obj.Neighbor(eidx) && obj.Neighbor(eidx)->Region_ID() != -999 )
                                    halo_cells.insert( obj.Neighbor(eidx) );
                             } );

    // 3. update the connectivity of the halo elements
    // (here extra checks are performed in DEBUG mode to determine whether they have other nullptr neighbors like at model boundary)
    for ( auto& it : halo_cells ) {
         assert( it );
         const auto n_nbors = it->Neighbors();
         for ( uint32_t eidx{0u}; eidx<n_nbors; ++eidx )
           // if the (face) neighbor of this object is going to be deleted, set it to nullptr
           if ( it.Neighbor(eidx) && it.Neighbor(eidx)->Region_ID() == -999 )
             it.UnassignNeighbor(eidx);
      }

 } // end updateHaloCellConnectivity

template void updateHaloCellConnectivity<3,Element>( plf::colony<Element<3>>::iterator, plf::colony<Element<3>>::iterator );
*/



} // end namespace csmp
