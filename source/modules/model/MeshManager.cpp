#include "MeshManager.h"
#include "PropertyDatabase.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Region.h"
#include "VSet.h"
#include "Exception.h"
#include "AP_BoolVector.h"
#include "CSMP_highLevelUtilities.h"
#include "ErrorHandler.h"
#include "PropertyData.h"
#include "Box.h"

using namespace std;

namespace csmp {

// refactored
template<size_t dim>
MeshManager<dim>::MeshManager()
 : hybrid_element_mesh_(false)
 {
 }


// refactored
template<size_t dim>
MeshManager<dim>::MeshManager( const PropertyDatabase<dim>& pref, const FiniteElementManager& fem_manager, VSet<dim>& vset )
 : hybrid_element_mesh_(vset.HybridElementTypeMesh())
 {
    BuildElementsAndVariableStorage( pref, fem_manager, vset );
    InitializeConnectivity( vset );
 }



// NB: here we have to deallocate all the dynamic storage, NOT DONE YET! 
template<size_t dim>
MeshManager<dim>::~MeshManager()
 {
 }
 
// TO FIX PROPERLY
template<size_t dim>
MeshManager<dim>::MeshManager( const MeshManager<dim>& mmgr )
 {
    *this = mmgr;
 }

 

/**
 
Assigns a copy of the supplied MeshManager to this MeshManager manager and
re-creates the connectivity for the new storage locations. 

@param mmgr A MeshManager object.


@section implementation Implementation

Apart from copying the storage vectors for Nodes, IntegrationPoints and
Elements, the Elements must be initialized by updating their connections
to the new storage locations of their member Nodes, IntegrationPoints and
neighbor elements. This extra step is required since the 
connections are pointers that are not re-addressed when the vectors are
copied. 
*/
template<size_t dim>
MeshManager<dim>&  MeshManager<dim>::operator=( const MeshManager<dim>& mmgr )
 {
    if ( this == &mmgr ) return *this;
    
    throw csmp::Exception( WARNING, "MeshManager<dim>::operator=", 
                                    "operator has not been tested yet");

    hybrid_element_mesh_ = mmgr.hybrid_element_mesh_;
    node_collection_.Assign(mmgr.node_collection_);
    elmt_collection_.Assign(mmgr.elmt_collection_);
    face_collection_.Assign(mmgr.face_collection_);
    interface_collection_.Assign(mmgr.interface_collection_);
   
    return *this;
   
 } // end operator=





/**
    Configures the distributed variable storage, inialises the finite element container 
    with the element types that are contained in the current mesh and builds the 
    mesh including the connectivity among its elements.
    
    @note the neighbors of each element include only the elements of the same type, i.e.
    a line element only has line neighbors, a surface element surfaec element neighbors 
    and so forth.
    
    @attention lower dimensional elements may have multiple neighbors for each of their
    faces. Yet only one of them will be assigned.
    
    @todo SKM: create manifolds to deal with lower-dimensional elements that have 
    multiple neighbors per face.
*/
template<size_t dim>
bool  MeshManager<dim>::Initialize( const PropertyDatabase<dim>& phys_vars,
                                    const FiniteElementManager& fem_manager,
                                    const VSet<dim>& vset )
 {
    hybrid_element_mesh_ = vset.HybridElementTypeMesh();
    BuildElementsAndVariableStorage( phys_vars, fem_manager, vset );
    return InitializeConnectivity( vset );
 }



/**
     Rebuild the mesh from an input dataset that has already been verified 
     including the set-up of the property storage.
*/
template<size_t dim>
bool  MeshManager<dim>::Reconstruct( const PropertyDatabase<dim>& phys_vars,
                                     const FiniteElementManager& fem_manager,
                                     const VSet<dim>& vset )
 {
    hybrid_element_mesh_ = vset.HybridElementTypeMesh();
    ReconstructMeshAndVariableStorage( phys_vars, fem_manager, vset );
    return InitializeVerifiedConnectivity( vset );
 }


 
 

/**
 
Builds arbitrary (potentially mixed element type) MeshManageres from the VSet
supplied. Nodes and elements are numbered 0..n-1U.  

@section implementation Implementation

First, the validity of the supplied finite element types is checked.

The mesh is stored into elmt_collection_  and node_collection_ deques.

@attention these deques are not needed or need to be removed if tree-only structure is wanted
to facilitate dynamic mesh adaptation.
By default, adaptive_remeshing_ is set to false.

@section application Application

This is the one that is mostly used todate.

 */
template<size_t dim>
bool MeshManager<dim>::BuildElementsAndVariableStorage( const PropertyDatabase<dim>& phys_vars,
                                                        const FiniteElementManager& fem_manager,
                                                        const VSet<dim>& vset )
 {
    // 0. if the MeshManager has already been initialized all that storage is deallocated and things 
    //    are build from scratch
    assert( elmt_collection_.empty() );    

    // 1. Checking the availability of the necessary finite element types
    // ------------------------------------------------------------------
    set<CSMP_FEM_TYPE> input_etypes;
    input_etypes.insert( parseFiniteElementTypeEnum(vset.ElementType(0U)) );
    
    if ( vset.HybridElementTypeMesh() )
      for ( size_t i=0U; i<vset.Simplices(); ++i )
        input_etypes.insert( parseFiniteElementTypeEnum(vset.ElementType(i)) );

    cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage(VSet): ";
    cout <<"input VSet contains the following finite element types:\n\t";
    for ( auto iit=input_etypes.begin(); iit!=input_etypes.end(); iit++ ) {
         cerr << parseFiniteElementType((*iit)) <<"  ";
         if ( !fem_manager.ContainsElementType(*iit) ) {
              cerr <<"\n\n\tFinite element type not available: "<< parseFiniteElementType( *iit ) << endl;
              fem_manager.Out();
              throw csmp::Exception( FATAL_ERROR, 
                             "MeshManager<dim>::BuildElementsAndVariableStorage(VSet):",
                             "'FiniteElementManager' lacks finite-element type required by VSet." );
          }
      }
    cout << endl;

    // 2. constructing the elements using the VSet nodes information
    // -------------------------------------------------------------
    // different element types are intrinsic to the VSet if so initialized
    typename deque<vector<size_t> >::const_iterator  first(vset.PlistBegin()), last(vset.PlistEnd());
    size_t  idx(0U);
    
    try {    
        /// Variable storage for elements and integration points
        const LocalVariables evars( phys_vars.LocalVariablesAt(ELEMENT) );
        const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt(ELEMENT) );

        // 2.1 If the MeshManager contains only one element type
        if ( !vset.HybridElementTypeMesh() ) {
             const int32 csmpElementType = vset.ElementType(0U);
             while( first != last ) {
                  auto e = elmt_collection_.Emplace( fem_manager.E(csmpElementType) );
                  e->ResizePropertyStorage( evars, cvars );
                  e->Idx(idx);
                  idx++;
                  first++;
             }
          }
        // 2.2 If there are multiple element types
        else
            while( first != last ) { 
                  // checking whether the element type matches the number of nodes for the given element
                  if ( (*first).size() != fem_manager.E(vset.ElementType(idx))->Nodes() ) {
                       cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
                       cout <<"Element: "<< parseFiniteElementType(vset.ElementType(idx));
                       cout <<" has a different number of nodes than are stored in the VSet: "<< (*first).size() << endl;
                    } 
                  int32 csmpElementType = vset.ElementType(idx);
                  auto e = elmt_collection_.Emplace( fem_manager.E(csmpElementType) );
                  e->ResizePropertyStorage( evars, cvars );
                  e->Idx(idx);
                  idx++;
                  first++;
               }

        // 3. construction the objects within the new memory
        // -------------------------------------------------
        const LocalVariables nvars( phys_vars.LocalVariablesAt(NODE) );

        Node<dim>  default_node;
        for ( size_t i=0U; i<vset.Vertices(); i++ ) {
             auto n = node_collection_.Emplace( default_node );
             n->ResizePropertyStorage(nvars);
             n->Idx(i);
          } 
     }
     
    // catching all possible standard exceptions & csmp::Exceptions
    catch( bad_alloc& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"bad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
      }
    catch( bad_exception& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"bad_exception: Exception error caused by: "<< ba.what() << endl;
      }
    catch( bad_typeid& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"bad_typeid: Type ID error caused by: "<< ba.what() << endl;
      }
    catch( ios_base::failure& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"ios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
      }
    // standard logic errors
    catch( domain_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"domain_error: Logic error caused by: "<< ba.what() << endl;
      }
    catch( invalid_argument& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"invalid_argument: Logic error caused by: "<< ba.what() << endl;
      }
    catch( length_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"length_error: Logic error caused by: "<< ba.what() << endl;
       }
    catch( out_of_range& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"out_of_range: Logic error caused by: "<< ba.what() << endl;
      }
    // runtime errors
    catch( overflow_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"overflow_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch( range_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"range_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch( underflow_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"underflow_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch(...) {
         cout <<"\nMeshManager<"<< dim <<">::BuildElementsAndVariableStorage: ";
         cout <<"\nunspecified exception occurred."<< endl;
      }

    cout <<"\nMeshManager<"<< dim;
    cout <<">::BuildElementsAndVariableStorage(VSet): Collection sizes:" << endl;
    cout <<"Node collection:            " << node_collection_.size() << endl;
    cout <<"Element collection:         " << elmt_collection_.size() << endl;
    cout << endl;
 
    return true;     
 } // end BuildElementsAndVariableStorage
 



/**
    Rebuilds the storage for the mesh from the supplied VSet.
    
    @author SKM 4/4/2016
    
    @todo TODO: SKM do not seperate storage construction from Neighbor assigmnents, but build complete elements, faces, and interfaces
*/
template<size_t dim>
bool MeshManager<dim>::ReconstructMeshAndVariableStorage( const PropertyDatabase<dim>& phys_vars,
                                                          const FiniteElementManager& fem_manager,
                                                          const VSet<dim>& vset )
 {
    // 0. if the MeshManager has already been initialized all that storage is deallocated and things 
    //    are build from scratch
    assert( elmt_collection_.empty() );    
    assert( face_collection_.empty() );
    assert( interface_collection_.empty() );

    // 1. Checking the availability of the necessary finite element types
    // ------------------------------------------------------------------
    set<CSMP_FEM_TYPE> input_etypes;
    input_etypes.insert( parseFiniteElementTypeEnum(vset.ElementType(0U)) );
    
    if ( vset.HybridElementTypeMesh() )
      for ( size_t i=0U; i<vset.Simplices(); ++i )
        input_etypes.insert( parseFiniteElementTypeEnum(vset.ElementType(i)) );

    cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage(VSet): ";
    cout <<"input VSet contains the following finite element types:\n\t";
    for ( typename set<CSMP_FEM_TYPE>::const_iterator
          iit=input_etypes.begin(); iit!=input_etypes.end(); iit++ ) {
         cerr << parseFiniteElementType((*iit)) <<"  ";
         if ( !fem_manager.ContainsElementType(*iit) ) {
              cerr <<"\n\n\tFinite element type not available: "<< parseFiniteElementType( *iit ) << endl;
              fem_manager.Out();
              throw csmp::Exception( FATAL_ERROR, 
                             "MeshManager<dim>::ReconstructMeshAndVariableStorage(VSet):",
                             "'FiniteElementManager' lacks finite-element type required by VSet." );
          }
      }
    cout << endl;

    try {
        // 2. constructing the elements using the VSet element type information
        // --------------------------------------------------------------------
        typename deque<vector<size_t> >::const_iterator  first(vset.PlistElmtsBegin()), last(vset.PlistElmtsEnd());
        // running index that is associated with the assumption that the elements are numbered consecutively
        // followed by the faces and the interfaces
        size_t  idx(0U);
    
        // variable storage for elements and integration points
        const LocalVariables evars( phys_vars.LocalVariablesAt(ELEMENT) );
        const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt(ELEMENT) );

        // 2.1 If the MeshManager contains only one element type (and therefore no faces nor interfaces)
        if ( !vset.HybridElementTypeMesh() ) {
             assert( vset.ElementTypes() == 1 );
             const int32 csmpElementType = vset.ElementType(0U);
             while( first != last ) {
                  elmt_collection_.Emplace( idx++, fem_manager.E(csmpElementType), evars, cvars, NOT );
                  ++first;
              }
          }
        // 2.2 If there are multiple element types
        else {
            while( first != last ) { 
#ifndef NDEBUG
                  assert( vset.ElementTypes() == vset.Elements() + vset.Faces() + vset.InterFaces() );
                  // checking whether the element type matches the number of nodes for the given element
                  if ( (*first).size() != fem_manager.E(vset.ElementType(idx))->Nodes() ) {
                       cerr <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
                       cerr <<"Element: "<< parseFiniteElementType(vset.ElementType(idx));
                       cerr <<" has a different number of nodes than are stored in the VSet: "<< (*first).size() << endl;
                    } 
#endif
                  const int32 csmpElementType = vset.ElementType(idx);
                  elmt_collection_.Emplace( idx, fem_manager.E(csmpElementType), evars, cvars, NOT );
                  ++idx;
                  ++first;
               }
          std::cerr << "Number of elements: " << idx << '\n';
          std::cerr << "Number of elements: " << elmt_collection_.size() << '\n';
        }

        // 3. constructing the faces using the VSet nodes information
        // ----------------------------------------------------------
        // (continuous running index 'idx' will be used so that face-IDs start with n-elements)
        // different face types are intrinsic to the VSet if so initialized
        if ( vset.Faces() > 0 ) {
            // since faces are lower-dimensional, the mesh must contain different element types
            assert( vset.HybridElementTypeMesh() );
            const LocalVariables evars( phys_vars.LocalVariablesAt(FACE) );
            const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt(FACE) );
          
            typename deque<vector<size_t> >::const_iterator  first(vset.PlistFacesBegin()), last(vset.PlistFacesEnd());
            while( first != last ) {
#ifndef NDEBUG
                  // checking whether the element type matches the number of nodes for the given element
                  if ( (*first).size() != fem_manager.E(vset.ElementType(idx))->Nodes() ) {
                       cerr <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
                       cerr <<"Face: "<< parseFiniteElementType(vset.ElementType(idx));
                       cerr <<" has a different number of nodes than are stored in the VSet: "<< (*first).size() << endl;
                    } 
#endif
                  const int32 csmpElementType = vset.ElementType(idx);
                  face_collection_.Emplace( idx, fem_manager.E(csmpElementType), evars, cvars );
                  ++first;
                  ++idx;
               }
          } // end faces
          
        // 4. constructing the interfaces using the VSet nodes information
        // ---------------------------------------------------------------
        // (continuous running index 'idx' will be used)
        if ( vset.InterFaces() > 0 ) {
            assert( vset.HybridElementTypeMesh() );
            const LocalVariables evars( phys_vars.LocalVariablesAt(INTER_FACE) );
            const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt(INTER_FACE) );

            typename deque<vector<size_t> >::const_iterator  first(vset.PlistInterFacesBegin()), last(vset.PlistInterFacesEnd());
            while( first != last ) {
#ifndef NDEBUG
                      // checking whether the element type matches the number of nodes for the given element
                      if ( (*first).size() != fem_manager.E(vset.ElementType(idx))->Nodes() ) {
                           cerr <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
                           cerr <<"Element: "<< parseFiniteElementType(vset.ElementType(idx));
                           cerr <<" has a different number of nodes than are stored in the VSet: "<< (*first).size() << endl;
                        } 
#endif
                      const int32 csmpElementType = vset.ElementType(idx);
                      interface_collection_.Emplace( idx++, fem_manager.E(csmpElementType), evars, cvars );
                      ++idx;
                      ++first;
                   }

            } // end interfaces


        // 5. construction of Node objects including coordinates and property storage allocation
        // -------------------------------------------------------------------------------------
        const LocalVariables nvars( phys_vars.LocalVariablesAt(NODE) );
        vector<double64>     coord(dim);

        for ( size_t i=0U; i<vset.Vertices(); ++i ) {
             for ( size_t j=0U; j<dim; ++j ) coord[j] = vset.P( j, i );
             // emplaced construction: Node( size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY=NOT );
             node_collection_.Emplace( i, Point<dim>(coord), nvars, NOT );
          }
     }
    // catching all possible standard exceptions & csmp::Exceptions
    catch( bad_alloc& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"bad_alloc: Memory allocation error caused by: "<< ba.what() << endl;
      }
    catch( bad_exception& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"bad_exception: Exception error caused by: "<< ba.what() << endl;
      }
    catch( bad_typeid& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"bad_typeid: Type ID error caused by: "<< ba.what() << endl;
      }
    catch( ios_base::failure& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"ios_base::failure: Probable I/O error caused by: "<< ba.what() << endl;
      }
    // standard logic errors
    catch( domain_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"domain_error: Logic error caused by: "<< ba.what() << endl;
      }
    catch( invalid_argument& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"invalid_argument: Logic error caused by: "<< ba.what() << endl;
      }
    catch( length_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"length_error: Logic error caused by: "<< ba.what() << endl;
       }
    catch( out_of_range& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"out_of_range: Logic error caused by: "<< ba.what() << endl;
      }
    // runtime errors
    catch( overflow_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"overflow_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch( range_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"range_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch( underflow_error& ba ) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"underflow_error: Runtime error caused by: "<< ba.what() << endl;
      }
    catch(...) {
         cout <<"\nMeshManager<"<< dim <<">::ReconstructMeshAndVariableStorage: ";
         cout <<"\nunspecified exception occurred."<< endl;
      }

    cout <<"\nMeshManager<"<< dim;
    cout <<">::ReconstructMeshAndVariableStorage(VSet): Collection sizes:" << endl;
    cout <<"\t\tNode collection:            " << node_collection_.size() << endl;
    cout <<"\t\tElement collection:         " << elmt_collection_.size() << endl;
    cout <<"\t\tFace collection:            " << face_collection_.size() << endl;
    if ( !interface_collection_.empty() )
      cout <<"\t\tInterFace collection:       " << interface_collection_.size() << endl;
    cout << endl;
 
    return true;
   
 } // end ReconstructMeshAndVariableStorage
 





/**

  Initializes the MeshManager from Vdata using following 'VSet' format:
  
  pelement:  CSMP element type of each element from which information needed to 
             read the plist of a multi-element type mesh can be deduced.

  plist:     node ID's per element

  pfverts:   list of elements opposite to respective nodes in each triangle
             (pfverts is written such that the face edges are stored 
             corresponding to the node entries in the pelements array) 

  px,py,pz:  coordinates of nodes [nodes]

 the supplied pointers to specific Finite Elements can be used if the
 VSet contains flags that identify which element types the specific elements
 are (this could be done using a number). The default is fes[0] and is
 used if no other information is supplied.

     - the property collection must be build independently:
       property_collection.SetStorageSpecifications( phys_vars, MeshManager_specifications );
       property_collection.BuildPropertyStorage(); 
*/
template<size_t dim>
bool MeshManager<dim>::InitializeConnectivity( const VSet<dim>& vset )
 {

    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: building storage..."<< endl;

    // --------------------------------
    // 1. assigning node coordinates
    // --------------------------------
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: assigning node coordinates..."<< endl;
   if      ( dim == 1U ) {
     size_t i(0);
     for ( auto& n : node_collection_ ) {
       n.x( vset.Px( i ) );
       ++i;
     }
   }
   else if ( dim == 2U ) {
     size_t i(0);
     for ( auto& n : node_collection_ ) {
           n.x( vset.Px( i ) );
           n.y( vset.Py( i ) );
           ++i;
        }
   }
   else {
     size_t i(0);
     for ( auto& n : node_collection_ ) {
           n.x( vset.Px( i ) );
           n.y( vset.Py( i ) );
           n.z( vset.Pz( i ) );
           ++i;
        }
   }


    // --------------------------------------------------
    // 2. assigning nodes to elements using 'plist' array
    // --------------------------------------------------
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: assigning nodes to elements..."<< endl;
    for ( auto& e : elmt_collection_ ) {
       for ( size_t j=0U; j<e.Nodes(); j++ )
          e.Assign( j, &node_collection_.Index( vset.Plist(e.Idx(),j) ) );
    }


    // --------------------------------------------------------------
    // 3. Assigning neighbor elements to elements (face-verts) if any
    // --------------------------------------------------------------
    if ( vset.WithNeighbourConnectivity() ) {
         for ( auto& e : elmt_collection_ ) {
               const size_t neighbors(e.Neighbors());
               for ( size_t j=0U, nidx=0u; j<neighbors; ++j ) {
                 // if there is a neighbor (as is the case if the stored index is greater than zero)
                 const int32 index(static_cast<int32>(vset.Pfvert( e.Idx(), j )) );
                 if ( index >= 0 )
                   e.Assign( nidx++, &elmt_collection_.Index( static_cast<size_t>(index) ) );
                 else
                   e.Assign( nidx++, static_cast<Element<dim>*>(nullptr) );
             }
         }
      }
    else
      csmp_error.notice( WARNING, "MeshManager::InitializeConnectivity:", "Input VSet does not contain any neighbor connectivity; nothing was done.");

    // ---------------------------------------------------------------------
    // 4. Flagging nodes located at the model boundary
    // ---------------------------------------------------------------------
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: flagging boundary objects..."<< endl;
    if ( vset.BFlags() > 0 )
      {
        // nodes were initially constructed as not located at the model boundary
        for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); ++bit )
          node_collection_.Index( (*bit).first ).AtBoundary( intToBOX_BOUNDARY( (*bit).second ) );
      }
    

    // ---------------------------------------------------------------------
    // 5. Flagging the elements using boundary flags from the nodes
    // ---------------------------------------------------------------------
    flagElementUsingNodal_BOX_BOUNDARY_Flags<dim>( ElementsBegin(), ElementsEnd() );
    /* TODO: clean after verification that code works as expected
    const bool vset_has_connectivity_info( vset.ElementNeighbors() > 0 );
    if( vset_has_connectivity_info )
    {
        size_t        n, bnodes;
        BOX_BOUNDARY  bflag;
        bool          assigned;

        for ( typename deque<Element<dim> >::iterator
              eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
          {
             // element is assumed to have been flagged as NOT
             assert( (*eit).AtBoundary() == NOT );

             // getting those elements which lie at an internal boundary, i.e., one of their
             // neighbors has been flagged as internal or none of the neighbor flags can be matched
             for ( size_t j=0U; j<(*eit).Neighbors(); j++ )
                 if ( (*eit).Neighbor(j) == NULL  &&  vset.Pfvert((*eit).Idx(),j) <= REGION_BOUNDARY ) {
                     (*eit).AtBoundary( INTERNAL );
                     break;
                 }

             if ( (*eit).AtBoundary() == NOT )
               for ( assigned=false, bnodes=n=0U; n<(*eit).Nodes(); n++ )
                 if ( (bflag=(*eit).N(n)->AtBoundary()) != NOT ) {
                      // dealing with the element corners
                      if      ( bflag == CNR1 ) { (*eit).AtBoundary( CNR1 ); break; }
                      else if ( bflag == CNR2 ) { (*eit).AtBoundary( CNR2 ); break; }
                      else if ( bflag == CNR3 ) { (*eit).AtBoundary( CNR3 ); break; }
                      else if ( bflag == CNR4 ) { (*eit).AtBoundary( CNR4 ); break; }
                      else if ( bflag == CNR5 ) { (*eit).AtBoundary( CNR5 ); break; }
                      else if ( bflag == CNR6 ) { (*eit).AtBoundary( CNR6 ); break; }
                      else if ( bflag == CNR7 ) { (*eit).AtBoundary( CNR7 ); break; }
                      else if ( bflag == CNR8 ) { (*eit).AtBoundary( CNR8 ); break; }
                      // now dealing with the other cases (if more than 2 bnodes are equal,
                      // their flag will be used to define the boundary)
                      if ( ++bnodes >= 2U  or  (bnodes >= 1U  and  dim == 1U) ) {
                           (*eit).AtBoundary( ((*eit).N(n))->AtBoundary() );
                           assigned = true;
                           break;
                        }
                   }
          }
    }
   */

   // ------------------------------------------------------------------------------
   // 6. Assigning parent elements (these are the elements that share the node) and
   // their respective internal node-id numbers to the nodes
   // -------------------------------------------------------------------------------
   if( csmp_error.Verbose() )
       cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: assigning parent element information to nodes..."<< endl;
   vector<size_t>  parent_elmts_per_node( node_collection_.size(), 0U );

   // counting how many parent elements each node has
   for ( auto& e : elmt_collection_ )
     for ( auto nit=e.NodesBegin(); nit!=e.NodesEnd(); nit++ )
       parent_elmts_per_node[ (*nit)->Idx() ]++;

   // reserving the memory for the parent storage and zeroing parent vector for next step
   size_t i(0);
   for ( auto& n : node_collection_ )
     n.ResizeParentStorage( parent_elmts_per_node[i++] );

   // assigning the parent element information to the nodes
    for ( auto& e : elmt_collection_ )
      for ( size_t j=0U; j<e.Nodes(); j++ )
        e.N(j)->Assign( j, &e );

   return true;
 
 } // end InitializeConnectivity(VSet)


/**
   Rebuild parent relationships.

TODO: Better docs
*/
template<size_t dim>
void MeshManager<dim>::RebuildParentRelationships(
        typename std::vector<csmp::Node<dim>*>::iterator begin, typename std::vector<csmp::Node<dim>*>::iterator end)
{
   std::deque< std::pair<size_t,Element<dim>*>> parents;
   for ( auto it = begin; it != end; ++it ) {
      auto n = *it;
      for ( size_t i = 0; i < n->Parents(); ++i ) {
          auto e = n->Parent(i);
          if (Verify(e)) {
              parents.push_back(std::make_pair(n->ParentNodeNumber(i),e));
          }
      }
      n->EraseParents();
      n->ResizeParentStorage(parents.size());
      for ( auto& p : parents ) {
          n->Assign(p.first, p.second);
      }
      parents.clear();
   }
}


/**
    Efficient initialisation without checks, assuming that these were done
    earlier before a model was saved to binary file.
 
    @note This method not recreates the element connectivity but also the box flags.
    
    @note use this method to bring a pre-existing CSMP model back from
    binary file.
    
    @attention methods assumes that elements - faces - interfaces were written in
    consecutive sequence to the VSet.
*/
template<size_t dim>
bool MeshManager<dim>::InitializeVerifiedConnectivity( const VSet<dim>& vset )
 {
    // prior initialisation expected
    assert( !elmt_collection_.empty() );
    assert( !node_collection_.empty() );
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: building storage..."<< endl;

    // -------------------------------------
    // 1. assigning coordinates to the nodes
    // -------------------------------------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning node coordinates..."<< endl;
   if      ( dim == 1U ) {
     size_t i(0);
     for ( auto& n : node_collection_ )
        n.x( vset.Px( i++ ) );
   }
   else if ( dim == 2U ) {
     size_t i(0);
     for ( auto& n : node_collection_ ) {
       n.x( vset.Px( i ) );
       n.y( vset.Py( i ) );
       ++i;
     }
   }
   else {
     size_t i(0);
     for ( auto& n : node_collection_ ) {
       n.x( vset.Px( i ) );
       n.y( vset.Py( i ) );
       n.z( vset.Pz( i ) );
       ++i;
      }
   }

    // ------------------------------------------------------------------------
    // 2. assigning nodes to elements, faces and interfaces using 'plist' array
    // ------------------------------------------------------------------------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning nodes to elements..."<< endl;
    for ( auto& e : elmt_collection_ ) {
           const size_t nodes(e.Nodes());
           for ( size_t j=0U; j<nodes; ++j )
             e.Assign( j, &node_collection_.Index( vset.Plist(e.Idx(),j) ) );
       }
    // faces
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning nodes to faces..."<< endl;
    for ( auto& e : face_collection_ ) {
           const size_t nodes(e.Nodes());
           for ( size_t j=0U; j<nodes; ++j )
             // assigning node indices
             e.Assign( j, &node_collection_.Index( vset.Plist(e.Idx(),j) ) );
       }
    // interfaces
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning nodes to interfaces..."<< endl;
    for ( auto& e : interface_collection_ ) {
           const size_t nodes(e.Nodes());
           for ( size_t j=0U; j<nodes; ++j ) {
                // assigning node indices
                e.Assign( j, &node_collection_.Index( vset.Plist(e.Idx(),j) ), INSIDE );
             }
       }

    // ----------------------------------------------------------------
    // 3. Assigning neighbor elements to elements, faces and interfaces
    // ----------------------------------------------------------------
    // 3.1 elements
    // ------------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning neighbors to elements..."<< endl;
    for ( auto& e : elmt_collection_ ) {
          const size_t neighbors(e.Neighbors());
          for ( size_t j=0U, nidx=0u; j<neighbors; ++j ) {
            // if there is a neighbor (as is the case if the stored index is greater than zero)
            const int32 index(static_cast<int32>(vset.Pfvert( e.Idx(), j )) );
            if (index >= 0) {
              e.Assign( nidx++, &elmt_collection_.Index( static_cast<size_t>(index) ) );
            }
            else {
                e.Assign( nidx++, static_cast<Element<dim>*>(nullptr) );
            }
        }
      }
    // 3.2 faces
    // ---------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: connecting faces to their higher-dimensional neighbors..."<< endl;
    const size_t elements(elmt_collection_.size());
    const size_t faces(face_collection_.size());
    for ( auto& e : face_collection_ ) {
         // equidimensional neighbors first
         // -------------------------------
         const size_t neighbors(e.Neighbors());
         for ( size_t j=0U; j<neighbors; ++j ) {
              // if there is a neighbor (as is the case if the stored index is greater than zero)
              // (e.Idx() starts with elements=first face)
              const long64 index(vset.Pfvert( e.Idx(), j ));
              if (  index >= 0 ) {
                   assert( index >= elements );
                   assert( index < elements + faces ); // (-) elements because face container is numbered from 0..n-1
                   e.Assign( j, &face_collection_.Index( static_cast<size_t>(index)-elements ) );
                }
              // nullptr assignment is not needed since this is the default initialisation
              //else e.Assign( j, static_cast<Face<dim>*>(nullptr) );
           }
        // higher-dimensional neighbors
        // ----------------------------
        // (are stored in VSet 'pfverts' record after the equidimensional neighbors)
        // index of inner neighbor element i which is always there
        const long64 index1(vset.Pfvert( e.Idx(), neighbors ));
        // index of outer neighbor element which may be there
        const long64 index2(vset.Pfvert( e.Idx(), neighbors+1U ));
        Element<dim>* const outerElement = (index2 < 0) ? nullptr : &elmt_collection_.Index(index2);
        // assigning inner and outer higher-dimensional neighbors
        //              inner element             outer element
        e.Assign( &elmt_collection_.Index(index1), outerElement );
      }
    // 3.3 interfaces
    // --------------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: connecting interfaces to their higher-dimensional neighbors..."<< endl;
    const size_t interfaces(interface_collection_.size());
    // connecting interfaces to their higher-dimensional neighbors
    for ( auto& e : interface_collection_ ) {
         // equidimensional neighbors first
         const size_t neighbors(e.Neighbors());
         for ( size_t j=0U; j<neighbors; ++j ) {
            // if there is a neighbor (as is the case if the stored index is greater than zero)
            const long64 index(vset.Pfvert( e.Idx(), j ));
            if (  index >= 0 ) {
                   assert( index >= elements + faces );
                   assert( index < elements + faces + interfaces ); // (-) because interface container is numbered from 0..n-1
                   e.Assign( j, &interface_collection_.Index( static_cast<size_t>(index)-elements-faces ) );
              }
            // nullptr assignment is not needed since this is the default initialisation
            //else e.Assign( j, static_cast<InterFace<dim>*>(nullptr) );
         }
        // higher-dimensional neighbors
        // ----------------------------
        // (the 2 sides will always be present because interfaces exist only on internal boundaries)
        e.Assign( &elmt_collection_.Index( vset.Pfvert( e.Idx(), neighbors ) ),
                  &elmt_collection_.Index( vset.Pfvert( e.Idx(), neighbors+1U ) ) );
      }

    // ---------------------------------------------------------------------
    // 4. Flagging nodes located at the model boundary
    // ---------------------------------------------------------------------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: flagging boundary objects..."<< endl;
    if ( vset.BFlags() > 0 )
      {
        // nodes were initially constructed as not located at the model boundary
        for ( auto bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++ )
          node_collection_.Index( (*bit).first ).AtBoundary( intToBOX_BOUNDARY( (*bit).second ) );
      }
    
    // ---------------------------------------------------------------------
    // 5. Flagging the elements using boundary flags from the nodes
    // ---------------------------------------------------------------------
    flagElementUsingNodal_BOX_BOUNDARY_Flags<dim>( ElementsBegin(), ElementsEnd() );

    // ------------------------------------------------------------------------------
    // 6. Assigning parent elements (these are the elements that share the node) and
    // their respective internal node-id numbers to the nodes
    // -------------------------------------------------------------------------------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning parent element information to nodes..."<< endl;
    vector<size_t>  parent_elmts_per_node( node_collection_.size(), 0U );

    // counting how many parent elements each node has
    for ( auto& e : elmt_collection_ ) {
      for ( auto nit=e.NodesBegin(); nit!=e.NodesEnd(); nit++ )
        parent_elmts_per_node[ (*nit)->Idx() ]++;
    }

    // reserving the memory for the parent storage and zeroing parent vector for next step
    {
        size_t i(0);
        for ( auto& n : node_collection_ ) {
          n.ResizeParentStorage( parent_elmts_per_node[i++] );
        }
    }

    // assigning the parent element information to the nodes
     for ( auto& e : elmt_collection_ ) {
       for ( size_t j=0U; j<e.Nodes(); j++ )
         e.N(j)->Assign( j, &e );
     }

    return true;
 
 } // end InitializeWithVerifiedConnectivity(VSet)








/**
    Counts elements the nodes of which are all located on the model boundary.
    
    @attention such elements typically give rise to problems with the assignment of 
    boundary conditions and should be eliminated.
    
    @attention method will work only if the elements are stored in elmt_collection_ deque.
    
    @author SKM
*/
template<size_t dim>
size_t MeshManager<dim>::DetectElementsWithAllNodesOnBoundary( set<size_t>& belmts ) const
 {
   belmts.clear();
   
   assert( !elmt_collection_.empty() );
   if ( elmt_collection_.empty() ) return 0U;
   
   size_t boundary_only_elements(0U);
   
   // counting how many parent elements each node has
   for ( auto& e : elmt_collection_ ) {
         size_t counter(0U);
         for ( auto nit=e.NodesBegin(); nit!=e.NodesEnd(); nit++ )
           if ( (*nit)->AtBoundary() != NOT )
             counter++;
         if ( counter == e.Nodes() ) {
               // cerr <<"\nMeshManager<dim>::DetectElementsWithAllNodesOnBoundary: found all-node-on-boundary element:";
               // (*eit).Out();
               belmts.insert( e.Idx() );
               boundary_only_elements++;
           }
     }
   
    return boundary_only_elements;
 }




/**
Assigns fem_manager to fvs_manager. Assigns the finite volume stencil pointers of the elements to the
corresponding finite volume stencils.

@attention If a stencil is assigned already, noting is done. Remove stencil first (NULL ptr in elements)

@todo (2-C) This should not also assign the femManager to the stencilManager - method does too much
*/
template<size_t dim>
void MeshManager<dim>::InitializeFiniteVolumeStencils( const PropertyDatabase<dim>& pref,
                                                       const FiniteElementManager& fem_manager,
                                                       FiniteVolumeStencilManager<dim>& fvs_manager )
 {
    if ( elmt_collection_.size() <= 1U )
       throw csmp::Exception( FATAL_ERROR, "MeshManager<dim>::InitializeFiniteVolumeStencils",
                                    "Currently no model exists to which stencils could be assigned.");

    // 1. initialize the stencil manager (the stencils are build and assigned the correct properties
    fvs_manager.Initialize( fem_manager );

    // 2. Now the stencil pointers in each finite element are connected to the correct corresponding stencils and update variable
    // storage for fv integration (sector/facet) point properties
    const LocalVariables lvs( pref.LocalVariablesAt(ELEMENT) );
    const IntegrationPointVariables ipvs( pref.IntegrationPointVariablesAt(ELEMENT) );
    for ( auto& e : elmt_collection_ ) {
          if ( !e.FV() )
            {
              e.AssignFiniteVolume( fvs_manager.Stencil( e.FE_Type() ) );
              e.ResizePropertyStorage( lvs, ipvs );
            }
    }

 } // end InitializeFiniteVolumeStencils



/// Unique Insert of corresponding entity

/// Insert of corresponding entity

/**
   avoids the call to the copy constructor of face
*/
template<size_t dim>
Node<dim>* MeshManager<dim>::PushBack( Node<dim>&& node )
{
    return node_collection_.Emplace( node );
}

template<size_t dim>
Element<dim>* MeshManager<dim>::PushBack( Element<dim>&& elmt )
{
    return elmt_collection_.Emplace( elmt );
}

template<size_t dim>
Face<dim>* const MeshManager<dim>::PushBack( Face<dim>&& face )
{
    return face_collection_.Emplace( face );
}

template<size_t dim>
InterFace<dim>* MeshManager<dim>::PushBack( InterFace<dim>&& interface )
{
    return interface_collection_.Emplace( interface );
}



template<size_t dim>
Node<dim>* MeshManager<dim>::PushBackIfUnique( Node<dim>&& node )
{
    auto nit = find( node_collection_.begin(), node_collection_.end(), node );
    if ( nit != node_collection_.end() )
       return ( &(*nit) );

    return node_collection_.Emplace( node );
}

template<size_t dim>
Element<dim>* MeshManager<dim>::PushBackIfUnique( Element<dim>&& elmt )
{
    auto eit = find( elmt_collection_.begin(), elmt_collection_.end(), elmt );
    if( eit != elmt_collection_.end() )
        return ( &(*eit) );

    return elmt_collection_.Emplace( elmt );
}

template<size_t dim>
Face<dim>* MeshManager<dim>::PushBackIfUnique( Face<dim>&& face )
{
    auto fit = find( face_collection_.begin(), face_collection_.end(), face );
    if( fit != face_collection_.end() )
        return ( &(*fit) );

    return face_collection_.Emplace( face );
}

template<size_t dim>
InterFace<dim>* MeshManager<dim>::PushBackIfUnique( InterFace<dim>&& interface )
{
    auto ifit = find( interface_collection_.begin(), interface_collection_.end(), interface );
    if( ifit != interface_collection_.end() )
        return ( &(*ifit) );
  
    return interface_collection_.Emplace( interface );
}


template<size_t dim>
void MeshManager<dim>::Erase( Node<dim>& node )
{
    node_collection_.Free(&node);
}

template<size_t dim>
void MeshManager<dim>::Erase( Element<dim>& elmt )
{
    elmt_collection_.Free(&elmt);
}

template<size_t dim>
void MeshManager<dim>::Erase( Face<dim>& face )
{
    face_collection_.Free(&face);
}

template<size_t dim>
void MeshManager<dim>::Erase( InterFace<dim>& interface )
{
    interface_collection_.Free(&interface);
}

/**
   erases all nodes, discerning those that do not have any parent element connections.
*/
template<size_t dim>
bool MeshManager<dim>::EraseNodes()
{
    bool emptyNodes(false);
    auto nit = node_collection_.begin();
    for( ; nit != node_collection_.end() ; )
    {
        if( (*nit).Parents() )
        {
            ++nit;
            //emptyNodes = true;
            //nit = node_collection_.erase( nit );
        }
        else
            ++nit;
    }
    return emptyNodes;
}

template<size_t dim>
bool MeshManager<dim>::EraseElements()
{
    bool emptyElements(false);
    auto eit = elmt_collection_.begin();
    for( ; eit != elmt_collection_.end() ; )
    {
        if( (*eit).NeighborElementVector().empty() && (*eit).NodeVector().empty() )
        {
            ++eit;
            //emptyElements = true;
            //eit = elmt_collection_.erase( eit );
        }
        else if( (*eit).NeighborElementVector().empty() )
        {
            (*eit).UnassignNodes();
            ++eit;
            //emptyElements = true;
            //eit = elmt_collection_.erase( eit );
        }
        else if( (*eit).NodeVector().empty() )
        {
            (*eit).UnassignNeighbors();
            ++eit;
            //emptyElements = true;
            //eit = elmt_collection_.erase( eit );
        }
        else
        {
            bool null_neighbors( true );
            for( size_t i = 0; i < (*eit).Neighbors(); i++ )
                if( (*eit).Neighbor( i ) != NULL )
                {
                    null_neighbors = false;
                    break;
                }
            if( null_neighbors )
            {
                (*eit).UnassignNodes();
                (*eit).UnassignNeighbors();
                ++eit;
                //emptyElements = true;
                //eit = elmt_collection_.erase( eit );
            }
            else
                ++eit;
        }
    }
    return emptyElements;
}

template<size_t dim>
bool MeshManager<dim>::EraseFaces( )
{
    bool emptyFaces(false);
    auto fit = face_collection_.begin();
    for( ; fit != face_collection_.end() ; )
    {
        if( (*fit).NeighborElementVector().empty() )
        {
            ++fit;
            //emptyFaces = true;
            //fit = face_collection_.erase( fit );
        }
        else
        {
            bool null_neighbors( true );
            for( size_t i = 0; i < (*fit).Neighbors(); i++ )
                if( (*fit).Neighbor( i ) != NULL )
                {
                    null_neighbors = false;
                    break;
                }
            if( null_neighbors )
            {
                ++fit;
                //emptyFaces = true;
                //fit = face_collection_.erase( fit );
            }
            else
                ++fit;
        }
    }
    return emptyFaces;
}


template<size_t dim>
bool MeshManager<dim>::EraseInterFaces( )
{
    bool emptyInterFaces(false);
    auto ifit = interface_collection_.begin();
    for( ; ifit != interface_collection_.end() ; )
    {
        if( (*ifit).NeighborElementVector().empty() )
        {
            ++ifit;
            //emptyInterFaces = true;
            //ifit = interface_collection_.erase( ifit );
        }
        else
        {
            bool null_neighbors( true );
            for( size_t i = 0; i < (*ifit).Neighbors(); i++ )
                if( (*ifit).Neighbor( i ) != NULL )
                {
                    null_neighbors = false;
                    break;
                }
            if( null_neighbors )
            {
                ++ifit;
                //emptyInterFaces = true;
                //ifit = interface_collection_.erase( ifit );
            }
            else
                ++ifit;
        }
    }
    return emptyInterFaces;
}





/** (Re)number all cells; either continuous for all cells or seperate ranges for all entity types

    @note this member function is constant because the idx_ is a mutable variable in the cell classes
*/
template<size_t dim>
void MeshManager<dim>::AssignUniqueNumbers( bool in_a_single_sequence ) const 
 {
    size_t n(0U);

    for_each( node_collection_.begin(), node_collection_.end(), [&n]( const Node<dim>& o ){ o.Idx( n++ ); return o; } );

    n = 0U; // resetting the counter
    for_each( elmt_collection_.begin(), elmt_collection_.end(), [&n]( const Element<dim>& o ){ o.Idx( n++ ); return o; } );

    if ( !in_a_single_sequence ) n = 0U;
    for_each( face_collection_.begin(), face_collection_.end(), [&n]( const Face<dim>& o ){ o.Idx( n++ ); return o; } );

    if ( !in_a_single_sequence ) n = 0U;
    for_each( interface_collection_.begin(), interface_collection_.end(), [&n]( const InterFace<dim>& o ){ o.Idx( n++ ); return o; } );
   
 } // end AssignUniqueNumbers





/** Writes MeshManager to VSet

Writes the elements, faces and interfaces 
stored in the current MeshManager 
to the supplied VSet. 
This also includes the connectivity information,
i.e. the connections between these entities.

@param vset A VSet preferably empty, which is
resized first, if necessary and into which the
MeshManager connectivity information is input.

@section conventions Conventions

Connections will exist between elements of the 
topology, i.e., between:

- volume elements
- surface elements
- line elements

these connections are made via pointers and exist
across the faces, edges, and end end point of these entities.

Where there is no neighbor, like for a face that sits on a model
boundary, the neighbor pointer will be a 'nullptr'.

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
template<size_t dim>
void MeshManager<dim>::OutputMeshTo( VSet<dim>& vset ) const
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    // 1. resizing the VSet
    // --------------------
    const size_t higherDimParents(2U);
    const size_t interfaceMultiplier(2U);
   
    if ( HybridElementMesh() or  Faces() > 0 or InterFaces() > 0 ) {
       deque<size_t>  nodes_per_element;
       deque<size_t>  neighbors_per_element;
       deque<int32>   csp_fem_types;

       // 1.1 identifying how many nodes and neighbors there are per element
      for ( auto& e : elmt_collection_ ) {
            nodes_per_element.push_back(e.Nodes());
            neighbors_per_element.push_back(e.Neighbors());
            csp_fem_types.push_back(e.FE_Type());
         }
      
       // 1.2 adding Face information after the elements
      for ( auto& e : face_collection_ ) {
            nodes_per_element.push_back(e.Nodes());
            neighbors_per_element.push_back( e.Neighbors() + higherDimParents );
            csp_fem_types.push_back(e.FE_Type());
         }
      
       // 1.3 adding InterFace information after the faces
      for ( auto& e : interface_collection_ ) {
            // multiplier takes care of the multiplicated interface nodes that the InterFace will be connected to
            nodes_per_element.push_back( e.Nodes() * interfaceMultiplier );
            neighbors_per_element.push_back( e.Neighbors() * interfaceMultiplier + higherDimParents );
            csp_fem_types.push_back(e.FE_Type());
         }
       // 1.4 resizing the VSet
       vset.Resize( csp_fem_types, nodes_per_element, neighbors_per_element,
                    node_collection_.size(), Faces(), InterFaces() );
      }
    else { // if there is only a single element type
        auto& el = elmt_collection_.Root();
        auto fe = el.FE();
        assert( fe != nullptr );
        vset.Resize( fe->Nodes(),
                     fe->Neighbors(),
                     fe->ElementType(),
                     node_collection_.size(), elmt_collection_.size() );
                     
        vset.ElementType( 0, el.FE_Type() );
      }
      
    // 2. adding node coordinates and boundary flags (BOX_BOUNDARY)
    // ------------------------------------------------------------
    size_t i(0U);
    if ( dim == 1U )
      for( auto& n : node_collection_ )
        vset.Px(i++,n.x());
    else if ( dim == 2U )
      for( auto& n : node_collection_ ) {
           vset.Px(i,n.x());
           vset.Py(i,n.y());
           ++i;
        }
    else
      for( auto& n : node_collection_ ) {
           vset.Px(i,n.x());
           vset.Py(i,n.y());
           vset.Pz(i,n.z());
           ++i;
        }

    // 3. adding 'plist' connectivity list
    // -----------------------------------
    size_t eidx = 0;
    // elements
    for ( auto& e : elmt_collection_ ) {
      for ( size_t j=0U; j<e.Nodes(); ++j )
        vset.Plist( eidx, j, (e.N(j)->Idx()) );
      ++eidx;
    }

    // faces
    for ( auto& e : face_collection_ ) {
      for ( size_t j=0U; j<e.Nodes(); ++j )
        vset.Plist( eidx, j, (e.N(j)->Idx()) );
      ++eidx;
    }
 
    // interfaces
    for ( auto& e : interface_collection_ ) {
      for ( size_t j=0U; j<e.Nodes(); ++j )
        vset.Plist( eidx, j, (e.N(j)->Idx()) );
      ++eidx;
    }
   
   
    // 4. adding 'pfverts' neighbors per element list
    // ----------------------------------------------
    eidx = 0;

    // 'pfverts' elements
    for ( auto& e : elmt_collection_ ) {
      for ( size_t j=0U; j<e.Neighbors(); ++j ) {
             Element<dim>* const ptr(e.Neighbor(j));
             if ( ptr != nullptr ) vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
             else                  vset.Pfvert( eidx, j, e.AtBoundary() );
        }
      ++eidx;
    }
   
    // 'pfverts' faces
    // ---------------
    // the faces are stored after the elements including connections to their higher-dimensional neighbors
    // add the end of the pfverts entries
    if ( !face_collection_.empty() )
      {
        for ( auto& f : face_collection_ ) {
             // equidimensional neighbors first
             const size_t neighbors(f.Neighbors());
             for ( size_t j=0U; j<neighbors; ++j ) {
                    Face<dim>* const ptr(f.Neighbor(j));
                    // if the neighbor exists (which it must on the inside of the Face)
                    if ( ptr != nullptr ) {
                         vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
                      }
                    else vset.Pfvert( eidx, j, f.InnerParent()->AtBoundary() );
               }
             // higher-dimensional neighbors second
             // inner neighbor
             if      ( dim == 3U ) assert( f.InnerParent()->IsVolumeElement() );
             else if ( dim == 2U ) assert( f.InnerParent()->IsSurfaceElement() );
             assert( f.InnerParent()->Idx() < elmt_collection_.size() );
             vset.Pfvert( eidx, neighbors, static_cast<int32>(f.InnerParent()->Idx()) );
             // outer neighbor
             if      ( f.OuterParent() != nullptr && dim == 3U ) assert( f.OuterParent()->IsVolumeElement() );
             else if ( f.OuterParent() != nullptr && dim == 2U ) assert( f.OuterParent()->IsSurfaceElement() );
             if ( f.OuterParent() != nullptr ) {
                  assert( f.OuterParent()->Idx() < elmt_collection_.size() );
                  vset.Pfvert( eidx, neighbors + 1U, static_cast<int32>(f.OuterParent()->Idx()) );
               }
             else {
                  // if there is no neighbor, the inner element parent should be at the model boundary
                  if ( f.InnerParent()->AtBoundary() == NOT ) {
                       f.Out();
                       csmp_error.notice( WARNING, "MeshManager<dim>::OutputMeshTo (face neighbors):",
                                         "inner dim+1 element should be at model boundary because Face has is no outer element.");
                       f.InnerParent()->AtBoundary( IRREGULAR );
                    }
                  vset.Pfvert( eidx, neighbors + 1U, f.InnerParent()->AtBoundary() );
               }
            ++eidx;
          }
      }
    // 'pfverts' interfaces
    // --------------------
    if ( !interface_collection_.empty() )
      {
        for ( auto& f : interface_collection_ ) {
             // equidimensional neighbors (=other interfaces) first
// TODO: each side will have neighbors on the separated sides of the interface; track!
             const size_t neighbors(f.Neighbors());
             for ( size_t j=0U; j<neighbors; ++j ) {
                    InterFace<dim>* const ptr(f.Neighbor(j));
                    if ( ptr != nullptr ) {
                         vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
                      }
                    else vset.Pfvert( eidx, j, REGION_BOUNDARY );
               }
             // inner and outer higher-dimensional neighbors
             // ( they must always exist because SplitBoundaries are internal model boundaries)
             assert( f.InnerParent() != nullptr );
             assert( f.OuterParent() != nullptr );
             vset.Pfvert( eidx, neighbors, static_cast<int32>(f.InnerParent()->Idx()) );
             vset.Pfvert( eidx, neighbors + 1U, static_cast<int32>(f.OuterParent()->Idx()) );
             ++eidx;
          }
      }
   
    // 5. adding boundary flags
    // ------------------------
   for( auto& n : node_collection_ ) {
      if ( n.AtBoundary() != NOT ) vset.AddBFlag( n.Idx(), n.AtBoundary() );
   }
         
    cout <<"\nMeshManager<"<< dim <<">::OutputMeshTo: MeshManager successfully output to VSet..."<< endl;

 } // end OutputMeshTo







/**
     Builds MeshManager using the VSet constructed exactly in the way as done by OutputMeshTo()
*/
/*
template<size_t dim>
void MeshManager<dim>::InputMeshFrom(  const PropertyDatabase<dim>& database, FiniteElementManager& fem_mgr,
                                       const VSet<dim>& vset, bool rebuild_mesh_from_scratch  )
 {
    ReconstructMeshAndVariableStorage( database, fem_mgr, vset );

    InitializeVerifiedConnectivity( vset );
   
    cout <<"\nMeshManager<"<< dim <<">::InputMeshFrom: MeshManager successfully initialised from VSet."<< endl;

 } // end InputMeshFrom
*/






/**  
    Storing distributed variables associated with the mesh in the VSet
    
    For all finite volumes and elements, their integration points and nodes, 
    but not for any of the variables stored on the model, region, boundaries or splitboundaries
    this method stores the current values in the VSet.
    
    To store the properties in the VSet, they are first written to PropertyData objects.
    These are then added to the VSet property storage.
    
    @attention  a continuous numbering of elements, faces, interfaces and nodes has to be 
    created with AssignUniqueNumbers() before this method is called.
    
    @note region properties are stored together with the regions in respective binary files
    
    @note this method is anything, but nice. Yet it will be quite a challenge to come up with a better
    design; hopefully there will be no more extra variable placements or types in the future.
    
    @attention in the VSet, the variables are identified only by their (unique) names. The property
    database is therefore essential to retrieve all other variable related information.
    
    @author SKM 5/5/2016
 
*/
template<size_t dim>
void MeshManager<dim>::OutputStoredVariablesTo( const PropertyDatabase<dim>& database, VSet<dim>& vset ) const
 {
    ErrorHandler&      csmp_error( ErrorHandler::Instance() );
    map<string,Index>  properties;
   
    // ---------------
    // ---------------
    // node properties
    // ---------------
    // ---------------
    database.ListProperties( NODE, properties );
    const size_t nodes(node_collection_.size());
    // for all node properties
    for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
      {
         // setting the specifications for the property storage (no memory allocation yet)
         PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
         const auto nodesEnd(node_collection_.end());
         // for the given property type
         const size_t flag_capacity( nodes * (*pit).second.flagDepth );
         const size_t data_capacity( nodes  * (*pit).second.dataDepth );
         // allocating memory to store the property flags and values
         data.Reserve( flag_capacity, data_capacity );

         switch( (*pit).second.type )
           {
             case SCALAR: {
                    // estimating the storage required
                    ScalarVariable value;
                    for ( auto it=node_collection_.begin(); it!=nodesEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    for ( auto it=node_collection_.begin(); it!=nodesEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    for ( auto it=node_collection_.begin(); it!=nodesEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    for ( auto it=node_collection_.begin(); it!=nodesEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    for ( auto it=node_collection_.begin(); it!=nodesEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of node variable not recognized.");
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
    const size_t elements(elmt_collection_.size());
    for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
      {
         // creating the property storage
         PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
         // estimating the storage required
         const auto elementsEnd(elmt_collection_.end());
         // for the given property type
         const size_t flag_capacity( elements * (*pit).second.flagDepth );
         const size_t data_capacity( elements  * (*pit).second.dataDepth );
         data.Reserve( flag_capacity, data_capacity );
        
         switch( (*pit).second.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                         (*it).Read( (*pit).second, value );
                         pushBack( data, value );
                      }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of element variable not recognized.");
           }
         // storing the data in the VSet
         vset.AddData( (*pit).first.c_str(), data );
      }

    // element integration point properties
    // ------------------------------------
    if ( database.ListProperties( ELEMENT_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !elmt_collection_.empty() );
        const size_t elmt_ips(elmt_collection_.Root().IntegrationPoints()); // just an estimate
       
        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             // creating the property storage
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto elementsEnd(elmt_collection_.end());
             // for the given property type
             const size_t flag_capacity( elements * elmt_ips * (*pit).second.flagDepth );
             const size_t data_capacity( elements  * elmt_ips * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of element integration point variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
          }
      }
   
   
    // element sector integration point properties
    // -------------------------------------------
    if ( database.ListProperties( SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !elmt_collection_.empty() );
        const size_t elmt_sector_ips(elmt_collection_.Root().IntegrationPointsPerSector());
        const size_t sectors_per_element(elmt_collection_.Root().Sectors());
       
        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto elementsEnd(elmt_collection_.end());

             const size_t flag_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.flagDepth );
             const size_t data_capacity( elements * sectors_per_element * elmt_sector_ips * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of element sector integraton point variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
         }
      }


    // element facet integration point properties
    // ------------------------------------------
    if ( database.ListProperties( FACET_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !elmt_collection_.empty() );
        const size_t elmt_facet_ips(elmt_collection_.Root().IntegrationPointsPerFacet());
        const size_t facets_per_element(elmt_collection_.Root().Facets());

        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto   elementsEnd(elmt_collection_.end());

             const size_t flag_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.flagDepth );
             const size_t data_capacity( elements * facets_per_element * elmt_facet_ips * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=elmt_collection_.begin(); it!=elementsEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of facet integration point variable not recognized.");
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
    const size_t faces(face_collection_.size());
   
    if ( !face_collection_.empty() )
      {
        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto facesEnd(face_collection_.end());

             const size_t flag_capacity( faces * (*pit).second.flagDepth );
             const size_t data_capacity( faces  * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of face variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
          }
      }


    // face integration point properties
    // ---------------------------------
    if ( !face_collection_.empty() && database.ListProperties( FACE_INTEGRATION_POINT, properties ) > 0 ) {
        const size_t face_ips(face_collection_.Root().IntegrationPoints());
       
        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto facesEnd(face_collection_.end());
             // for the given property type
             const size_t flag_capacity( faces * face_ips * (*pit).second.flagDepth );
             const size_t data_capacity( faces * face_ips  * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of face integration point variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
         }
      }
   

    // face sector integration point properties
    // ----------------------------------------
    if ( !face_collection_.empty() && database.ListProperties( FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
        const size_t face_sector_ips(face_collection_.Root().IntegrationPointsPerSector());
        const size_t sectors_per_face(face_collection_.Root().Sectors());

        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto facesEnd(face_collection_.end());
             // for the given property type
             const size_t flag_capacity( faces * sectors_per_face * face_sector_ips * (*pit).second.flagDepth );
             const size_t data_capacity( faces * sectors_per_face * face_sector_ips  * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of face-sector integration point variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
         }
      }
   
   
    // face facet integration point properties
    // ---------------------------------------
    if ( !face_collection_.empty() && database.ListProperties( FACE_FACET_INTEGRATION_POINT, properties ) > 0 ) {
        const size_t face_facet_ips(face_collection_.Root().IntegrationPointsPerFacet());
        const size_t facets_per_face(face_collection_.Root().Facets());

        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto facesEnd(face_collection_.end());
             // for the given property type
             const size_t flag_capacity( faces * facets_per_face * face_facet_ips * (*pit).second.flagDepth );
             const size_t data_capacity( faces * facets_per_face * face_facet_ips * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );
            
             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=face_collection_.begin(); it!=facesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of face facet integration point variable not recognized.");
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
    const size_t interfaces(interface_collection_.size());
    // only if there actually are interface in the model
    if ( !interface_collection_.empty() )
      {
        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto interfacesEnd(interface_collection_.end());
             const size_t flag_capacity( interfaces * (*pit).second.flagDepth );
             const size_t data_capacity( interfaces  * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             (*it).Read( (*pit).second, value );
                             pushBack( data, value );
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of interface variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
          }
      }

    // interface integration point properties
    // --------------------------------------
    if ( !interface_collection_.empty() && database.ListProperties( INTER_FACE_INTEGRATION_POINT, properties ) > 0 ) {
        const size_t interface_ips(interface_collection_.Root().IntegrationPoints());
       
        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto interfacesEnd(interface_collection_.end());
             const size_t flag_capacity( interfaces * interface_ips * (*pit).second.flagDepth );
             const size_t data_capacity( interfaces * interface_ips * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t integration_points((*it).IntegrationPoints());
                             for ( size_t i=0U; i<integration_points; ++i ) {
                                  (*it).Read( i, (*pit).second, value );
                                  pushBack( data, value );
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of interface integration point variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
          }
      }
   
    // interface sector integration point properties
    // ---------------------------------------------
    if ( !interface_collection_.empty() && database.ListProperties( INTER_FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
        const size_t interface_sector_ips(interface_collection_.Root().IntegrationPointsPerSector());
        const size_t sectors_per_interface(interface_collection_.Root().Sectors());

        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto interfacesEnd(interface_collection_.end());
             const size_t flag_capacity( interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.flagDepth );
             const size_t data_capacity( interfaces * sectors_per_interface * interface_sector_ips * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t sectors((*it).Sectors());
                             for ( size_t i=0U; i<sectors; ++i ) {
                                  const size_t ips_per_sector((*it).IntegrationPointsPerSector());
                                  for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of interface sector integration point variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
          }
      }
   
    // interface facet integration point properties
    // --------------------------------------------
    if ( !interface_collection_.empty() && database.ListProperties( INTER_FACE_FACET_INTEGRATION_POINT, properties ) > 0 ) {
        const size_t interface_facet_ips(face_collection_.Root().IntegrationPointsPerFacet());
        const size_t facets_per_interface(face_collection_.Root().Facets());

        for ( auto pit=properties.begin(); pit!=properties.end(); ++pit )
          {
             PropertyData  data( (*pit).second.place, (*pit).second.type, dim, (*pit).second.dataDepth );
             const auto interfacesEnd(interface_collection_.end());
             const size_t flag_capacity( interfaces * facets_per_interface * interface_facet_ips * (*pit).second.flagDepth );
             const size_t data_capacity( interfaces * facets_per_interface * interface_facet_ips * (*pit).second.dataDepth );
             data.Reserve( flag_capacity, data_capacity );

             switch( (*pit).second.type )
               {
                 case SCALAR: {
                        ScalarVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case VECTOR: {
                        VectorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case TENSOR: {
                        TensorVariable<dim> value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case ARRAY: {
                        ArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 case FLAGGEDARRAY: {
                        FlaggedArrayVariable value;
                        for ( auto it=interface_collection_.begin(); it!=interfacesEnd; ++it ) {
                             const size_t facets((*it).Facets());
                             for ( size_t i=0U; i<facets; ++i ) {
                                  const size_t ips_per_facet((*it).IntegrationPointsPerFacet());
                                  for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                       (*it).Read( i, j, (*pit).second, value );
                                       pushBack( data, value );
                                    }
                               }
                          }
                     }
                   break;
                 default:
                   csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                      (*pit).first, "type of interface facet integration point variable not recognized.");
               }
             // storing the data in the VSet
             vset.AddData( (*pit).first.c_str(), data );
          }
      }
 
 } // end OutputStoredVariablesTo







/**
    Read and assign property values written exactly by OutputStoredVariablesTo()
    
    @ test SKM 28/6/2016
*/
template<size_t dim>
void MeshManager<dim>::InputStoredVariablesFrom( const PropertyDatabase<dim>& database, const VSet<dim>& vset )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    // ---------------
    // ---------------
    // node properties
    // ---------------
    // ---------------
    // for all node properties
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         // apart from the name string key in the map, PropertyData contains the most important variable specifications
         if ( (*pit).second.Placement() != NODE ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == NODE );
         assert( (*pit).second.Size() == node_collection_.size() * key.dataDepth );
        
         // for the given property type
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                     size_t i = 0;
                     for ( auto& n : node_collection_) {
                         read( (*pit).second, i, value );
                         n.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t i = 0;
                    for ( auto& n : node_collection_) {
                      read( (*pit).second, i, value );
                      n.Store( key, value );
                      ++i;
                    }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t i = 0;
                    for ( auto& n : node_collection_) {
                      read( (*pit).second, i, value );
                      n.Store( key, value );
                      ++i;
                    }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t i = 0;
                    for ( auto& n : node_collection_) {
                      read( (*pit).second, i, value );
                      n.Store( key, value );
                      ++i;
                    }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t i = 0;
                    for ( auto& n : node_collection_) {
                      read( (*pit).second, i, value );
                      n.Store( key, value );
                      ++i;
                    }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of node variable not recognized.");
           }
      }
   
    // -------------------------------------------------
    // -------------------------------------------------
    // element properties (including integration points)
    // -------------------------------------------------
    // -------------------------------------------------
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != ELEMENT ) continue;
         // some checks
          assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == ELEMENT );
         assert( (*pit).second.Size() / key.dataDepth == elmt_collection_.size() );

         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t i(0);
                    for ( auto& e : elmt_collection_ ) {
                          read( (*pit).second, i, value );
                          e.Store( key, value );
                          ++i;
                       }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t i(0);
                    for ( auto& e : elmt_collection_ ) {
                          read( (*pit).second, i, value );
                          e.Store( key, value );
                          ++i;
                       }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t i(0);
                    for ( auto& e : elmt_collection_ ) {
                          read( (*pit).second, i, value );
                          e.Store( key, value );
                          ++i;
                       }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t i(0);
                    for ( auto& e : elmt_collection_ ) {
                          read( (*pit).second, i, value );
                          e.Store( key, value );
                          ++i;
                       }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t i(0);
                    for ( auto& e : elmt_collection_ ) {
                          read( (*pit).second, i, value );
                          e.Store( key, value );
                          ++i;
                       }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of element variable not recognized.");
           }
      }

    // element integration point properties
    // ------------------------------------
    // ELEMENT_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != ELEMENT_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == ELEMENT_INTEGRATION_POINT );

         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U); // running index
                    for ( auto& e : elmt_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                      const size_t integration_points(e.IntegrationPoints());
                      for ( size_t i=0U; i<integration_points; ++i ) {
                        read( (*pit).second, entry, value );
                        e.Store( i, key, value );
                        entry++;
                      }
                    }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U); 
                    for ( auto& e : elmt_collection_ ) {
                      const size_t integration_points(e.IntegrationPoints());
                      for ( size_t i=0U; i<integration_points; ++i ) {
                        read( (*pit).second, entry, value );
                        e.Store( i, key, value );
                        entry++;
                      }
                    }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                      const size_t integration_points(e.IntegrationPoints());
                      for ( size_t i=0U; i<integration_points; ++i ) {
                        read( (*pit).second, entry, value );
                        e.Store( i, key, value );
                        entry++;
                      }
                    }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U); 
                    for ( auto& e : elmt_collection_ ) {
                      const size_t integration_points(e.IntegrationPoints());
                      for ( size_t i=0U; i<integration_points; ++i ) {
                        read( (*pit).second, entry, value );
                        e.Store( i, key, value );
                        entry++;
                      }
                    }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of element integration point variable not recognized.");
           }
      }

    // element sector integration point properties
    // -------------------------------------------
    // SECTOR_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != SECTOR_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == SECTOR_INTEGRATION_POINT );

         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                                  (*pit).first, "type of element sector integraton point variable not recognized.");
           }
      }

    // element facet integration point properties
    // ------------------------------------------
    // FACET_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != FACET_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == FACET_INTEGRATION_POINT );

         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : elmt_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                                  (*pit).first, "type of facet integration point variable not recognized.");
           }
      }


   
    // ----------------------------------------------
    // ----------------------------------------------
    // face properties (including integration points)
    // ----------------------------------------------
    // ----------------------------------------------
    // FACE
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != FACE ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == FACE );
         assert( (*pit).second.Size() / key.dataDepth == face_collection_.size() );

         // for the given property type
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t i(0);
                    for ( auto& f : face_collection_ ) {
                         read( (*pit).second, i, value );
                         f.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t i(0);
                    for ( auto& f : face_collection_ ) {
                         read( (*pit).second, i, value );
                         f.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t i(0);
                    for ( auto& f : face_collection_ ) {
                         read( (*pit).second, i, value );
                         f.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t i(0);
                    for ( auto& f : face_collection_ ) {
                         read( (*pit).second, i, value );
                         f.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t i(0);
                    for ( auto& f : face_collection_ ) {
                         read( (*pit).second, i, value );
                         f.Store( key, value );
                         ++i;
                      }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of face variable not recognized.");
           }
      }

    // face integration point properties
    // ---------------------------------
    // FACE_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != FACE_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == FACE_INTEGRATION_POINT );

         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of face integration point variable not recognized.");
           }
      }

    // face sector integration point properties
    // ----------------------------------------
    // FACE_SECTOR_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != FACE_SECTOR_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == FACE_SECTOR_INTEGRATION_POINT );

         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                                  (*pit).first, "type of face-sector integration point variable not recognized.");
           }
      }

    // face facet integration point properties
    // ---------------------------------------
    // FACE_FACET_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != FACE_FACET_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == FACE_FACET_INTEGRATION_POINT );

         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : face_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                                  (*pit).first, "type of face facet integration point variable not recognized.");
           }
      }


   
    // ---------------------------------------------------
    // ---------------------------------------------------
    // interface properties (including integration points)
    // ---------------------------------------------------
    // ---------------------------------------------------
    // INTER_FACE
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != INTER_FACE ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == INTER_FACE );
         assert( (*pit).second.Size() / key.dataDepth == interface_collection_.size() );
        
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t i(0);
                    for ( auto& e : interface_collection_ ) {
                         read( (*pit).second, i, value );
                         e.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t i(0);
                    for ( auto& e : interface_collection_ ) {
                         read( (*pit).second, i, value );
                         e.Store( key, value );
                         ++i;
                      }
                }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t i(0);
                    for ( auto& e : interface_collection_ ) {
                         read( (*pit).second, i, value );
                         e.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t i(0);
                    for ( auto& e : interface_collection_ ) {
                         read( (*pit).second, i, value );
                         e.Store( key, value );
                         ++i;
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t i(0);
                    for ( auto& e : interface_collection_ ) {
                         read( (*pit).second, i, value );
                         e.Store( key, value );
                         ++i;
                      }
                }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of interface variable not recognized.");
           }
      }

    // interface integration point properties
    // --------------------------------------
    // INTER_FACE_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != INTER_FACE_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == INTER_FACE_INTEGRATION_POINT );
        
         switch( key.type )
          {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t integration_points(e.IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              e.Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             default:
               csmp_error.notice( ERROR, "Region<dim>::OutputVariableTo:",
                                  (*pit).first, "type of interface integration point variable not recognized.");
           }
      }

    // interface sector integration point properties
    // ---------------------------------------------
    // INTER_FACE_SECTOR_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != INTER_FACE_SECTOR_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == INTER_FACE_SECTOR_INTEGRATION_POINT );
        
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t sectors(e.Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(e.IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
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
                                  (*pit).first, "type of interface sector integration point variable not recognized.");
           }
      }

    // interface facet integration point properties
    // --------------------------------------------
    // INTER_FACE_FACET_INTEGRATION_POINT
    for ( auto pit=vset.PropertyValuesBegin(); pit!=vset.PropertyValuesEnd(); ++pit )
      {
         if ( (*pit).second.Placement() != INTER_FACE_FACET_INTEGRATION_POINT ) continue;
         // some checks
         assert( database.IsDefined( (*pit).first.c_str() ) );
         const csmp::Index key(database.StorageKey( (*pit).first.c_str() ));
         assert( key.place == INTER_FACE_FACET_INTEGRATION_POINT );
        
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   e.Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( auto& e : interface_collection_ ) {
                         const size_t facets(e.Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(e.IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
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
                                  (*pit).first, "type of interface facet integration point variable not recognized.");
           }
      }

 } // end InputStoredVariablesFrom




template<size_t dim>
bool  MeshManager<dim>::HybridElementMesh() const 
 { return hybrid_element_mesh_; }

// size of containers

template<size_t dim>
size_t  MeshManager<dim>::Nodes() const
 { return node_collection_.size(); }

template<size_t dim>
size_t  MeshManager<dim>::Elements() const
 { return elmt_collection_.size(); }

template<size_t dim>
size_t  MeshManager<dim>::Faces() const
 { return face_collection_.size(); }

template<size_t dim>
size_t  MeshManager<dim>::InterFaces() const
 { return interface_collection_.size(); }


// accessors

template<size_t dim>
csmp::Node<dim>&    MeshManager<dim>::RootNode()
 { return node_collection_.Root(); }

template<size_t dim>
csmp::Element<dim>&  MeshManager<dim>::RootElement()
 { return elmt_collection_.Root(); }

template<size_t dim>
csmp::Face<dim>&  MeshManager<dim>::RootFace()
 { return face_collection_.Root(); }

template<size_t dim>
csmp::InterFace<dim>&  MeshManager<dim>::RootInterFace()
 { return interface_collection_.Root(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Node<dim> >::iterator  MeshManager<dim>::NodesBegin()
 { return node_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Node<dim> >::iterator  MeshManager<dim>::NodesEnd()
 { return node_collection_.end(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Element<dim> >::iterator  MeshManager<dim>::ElementsBegin()
 { return elmt_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Element<dim> >::iterator  MeshManager<dim>::ElementsEnd()
 { return elmt_collection_.end(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Face<dim> >::iterator  MeshManager<dim>::FacesBegin()
 { return face_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Face<dim> >::iterator  MeshManager<dim>::FacesEnd()
 { return face_collection_.end(); }

template<size_t dim>
typename PrimitiveContainer<csmp::InterFace<dim> >::iterator  MeshManager<dim>::InterFacesBegin()
 { return interface_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::InterFace<dim> >::iterator  MeshManager<dim>::InterFacesEnd()
 { return interface_collection_.end(); }

// const accessors

template<size_t dim>
const csmp::Node<dim>&    MeshManager<dim>::RootNode() const
 { return node_collection_.Root(); }

template<size_t dim>
const csmp::Element<dim>&  MeshManager<dim>::RootElement() const
 { return elmt_collection_.Root(); }

template<size_t dim>
const csmp::Face<dim>&  MeshManager<dim>::RootFace() const
 { return face_collection_.Root(); }

template<size_t dim>
const csmp::InterFace<dim>&  MeshManager<dim>::RootInterFace() const
 { return interface_collection_.Root(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Node<dim> >::const_iterator  MeshManager<dim>::NodesBegin() const
 { assert( !node_collection_.empty() ); return node_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Node<dim> >::const_iterator  MeshManager<dim>::NodesEnd() const
 { assert( !node_collection_.empty() ); return node_collection_.end(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Element<dim> >::const_iterator  MeshManager<dim>::ElementsBegin() const
 { assert( !elmt_collection_.empty() ); return elmt_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Element<dim> >::const_iterator  MeshManager<dim>::ElementsEnd() const
 { assert( !elmt_collection_.empty() ); return elmt_collection_.end(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Face<dim> >::const_iterator  MeshManager<dim>::FacesBegin() const
 { assert( !face_collection_.empty() ); return face_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::Face<dim> >::const_iterator  MeshManager<dim>::FacesEnd() const
 { assert( !face_collection_.empty() ); return face_collection_.end(); }

template<size_t dim>
typename PrimitiveContainer<csmp::InterFace<dim> >::const_iterator  MeshManager<dim>::InterFacesBegin() const
 { assert( !interface_collection_.empty() ); return interface_collection_.begin(); }

template<size_t dim>
typename PrimitiveContainer<csmp::InterFace<dim> >::const_iterator  MeshManager<dim>::InterFacesEnd() const
 { assert( !interface_collection_.empty() ); return interface_collection_.end(); }


// NODES

template<size_t dim>
bool MeshManager<dim>::NodeAtIndexIsSafe() const
  { return node_collection_.IndexOperationIsSafe(); }
  
template<size_t dim>
csmp::Node<dim>& MeshManager<dim>::NodeAtIndex(size_t i)
  { assert( i < node_collection_.size() ); return node_collection_.Index(i); }
  
template<size_t dim>
const csmp::Node<dim>& MeshManager<dim>::NodeAtIndex(size_t i) const
  { assert( i < node_collection_.size() ); return node_collection_.Index(i); }

template<size_t dim>
csmp::Node<dim>* const MeshManager<dim>::PointerToNodeAtIndex( size_t i ) const
 { return const_cast<csmp::Node<dim>* const>(&NodeAtIndex(i)); }


// ELEMENTS

template<size_t dim>
bool MeshManager<dim>::ElementAtIndexIsSafe() const
  { return elmt_collection_.IndexOperationIsSafe(); }
  
template<size_t dim>
csmp::Element<dim>& MeshManager<dim>::ElementAtIndex(size_t i)
  { assert( i < elmt_collection_.size() ); return elmt_collection_.Index(i); }
  
template<size_t dim>
const csmp::Element<dim>& MeshManager<dim>::ElementAtIndex(size_t i) const
  { assert( i < elmt_collection_.size() ); return elmt_collection_.Index(i); }

template<size_t dim>
csmp::Element<dim>* const MeshManager<dim>::PointerToElementAtIndex( size_t i ) const
  { return const_cast<csmp::Element<dim>* const>(&ElementAtIndex(i)); }


// FACE

template<size_t dim>
bool MeshManager<dim>::FaceAtIndexIsSafe() const
  { return face_collection_.IndexOperationIsSafe(); }
  
template<size_t dim>
csmp::Face<dim>& MeshManager<dim>::FaceAtIndex(size_t i)
  { assert( i < face_collection_.size() ); return face_collection_.Index(i); }
  
template<size_t dim>
const csmp::Face<dim>& MeshManager<dim>::FaceAtIndex(size_t i) const
  { assert( i < face_collection_.size() ); return face_collection_.Index(i); }

template<size_t dim>
csmp::Face<dim>* const MeshManager<dim>::PointerToFaceAtIndex( size_t i ) const
 { return const_cast<csmp::Face<dim>* const>(&FaceAtIndex(i)); }


// INTERFACE

template<size_t dim>
bool MeshManager<dim>::InterFaceAtIndexIsSafe() const
  { return interface_collection_.IndexOperationIsSafe(); }
  
template<size_t dim>
csmp::InterFace<dim>& MeshManager<dim>::InterFaceAtIndex(size_t i)
  { assert( i < interface_collection_.size() ); return interface_collection_.Index(i); }
  
template<size_t dim>
const csmp::InterFace<dim>& MeshManager<dim>::InterFaceAtIndex(size_t i) const
  { assert( i < interface_collection_.size() ); return interface_collection_.Index(i); }

template<size_t dim>
csmp::InterFace<dim>* const MeshManager<dim>::PointerToInterFaceAtIndex( size_t i ) const
 { return const_cast<csmp::InterFace<dim>* const>(&InterFaceAtIndex(i)); }





template<size_t dim>
void MeshManager<dim>::Out() const
 {
    cout <<"\nMeshManager<"<< dim <<">::Out: "<< endl;
    
    // nodes
    cout <<"\nNODES: "<< endl;
    for ( auto& n : node_collection_ ) {
         string bound = parseBoundary(n.AtBoundary());
         cout <<"\nNode ID: "<< n.Idx() <<" ";
         cout << n.Coordinate();
         cout <<" Boundary flag: "<< bound << endl;
      }

    // elements
    cout <<"\nELEMENTS: "<< endl;
    for ( auto& e : elmt_collection_ ) {
         string bound = parseBoundary(e.AtBoundary());
         cout <<"\nElement ID: "<< e.Idx() <<" Boundary flag: "<< bound << endl;
         cout <<"Member Nodes: "<< endl;
         for ( size_t i=0U; i<e.Nodes(); i++ )
           cout << e.N(i)->Idx() <<"\t";
          cout <<"\nNeighbor elements: "<< endl;
         for ( size_t i=0U; i<e.Neighbors(); i++ )
           if ( e.Neighbor(i) != NULL )
             cout << e.Neighbor(i)->Idx() <<"\t";
           else
           cout <<"NO NEIGHBOR\t";        
 
         cout << endl;
     }

   // faces
    cout <<"\nFACES: "<< endl;
    for ( auto& e : face_collection_ ) {
         cout <<"\nFace ID: "<< e.Idx() <<" No boundary flags."<< endl;
         cout <<"Member Nodes: "<< endl;
         for ( size_t i=0U; i<e.Nodes(); i++ )
           cout << e.N(i)->Idx() <<"\t";
          cout <<"\nNeighbor elements: "<< endl;
         for ( size_t i=0U; i<e.Neighbors(); i++ )
           if ( e.Neighbor(i) != NULL )
             cout << e.Neighbor(i)->Idx() <<"\t";
           else
           cout <<"NO NEIGHBOR\t";
         // higher dimensional neighbors
         cout <<"\nhigher-dimensional neighbor elements:\n";
         assert( e.InnerParent() != nullptr );
         cout <<"\tinner: "<< e.InnerParent()->Idx() <<"\t";
         if ( e.OuterParent() != nullptr )
            cout <<"\touter: "<< e.OuterParent()->Idx() <<"\t";
          else
            cout <<"NO NEIGHBOR\t";
         cout << endl;
     }

   // inter faces
    cout <<"\nINTERFACES: "<< endl;
    for ( auto& e : interface_collection_ ) {
         cout <<"\nInterFace ID: "<< e.Idx() <<" No boundary flags."<< endl;
         cout <<"Member Nodes: "<< endl;
         for ( size_t i=0U; i<e.Nodes(); i++ )
           cout << e.N(i)->Idx() <<"\t";
          cout <<"\nNeighbor elements: "<< endl;
         for ( size_t i=0U; i<e.Neighbors(); i++ )
           if ( e.Neighbor(i) != NULL )
             cout << e.Neighbor(i)->Idx() <<"\t";
           else
           cout <<"NO NEIGHBOR\t";
         // higher dimensional neighbors
         cout <<"\nhigher-dimensional neighbor elements:\n";
         assert( e.InnerParent() != nullptr );
         cout <<"\tinner: "<< e.InnerParent()->Idx() <<"\t";
         if ( e.OuterParent() != nullptr )
            cout <<"\touter: "<< e.OuterParent()->Idx() <<"\t";
          else
            cout <<"NO NEIGHBOR\t";
         cout << endl;
     }

   // parent elements ID's for each node
   cout << endl << endl;   
    cout <<"PARENT ELEMENT INFORMATION FOR ALL NODES: "<< endl;
   for ( auto& n : node_collection_ )
     {
        cout <<"\nNode: "<< n.Idx() <<", parent elements: "<< endl;
        for ( size_t i=0u; i<n.Parents(); i++ )
          cout << n.Parent( i )->Idx() <<" ";
        cout << endl;
     } 

 } // end Out
 
 
 
template class MeshManager<1U>;
template class MeshManager<2U>;
template class MeshManager<3U>;



} // end namespace csmp 















