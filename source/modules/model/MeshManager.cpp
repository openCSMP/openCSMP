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
 : adaptive_remeshing_(false),  
   hybrid_element_mesh_(false)
 {
 }


// refactored
template<size_t dim>
MeshManager<dim>::MeshManager( const PropertyDatabase<dim>& pref, const FiniteElementManager& fem_manager, VSet<dim>& vset )
 : adaptive_remeshing_(false), 
   hybrid_element_mesh_(vset.HybridElementTypeMesh())
 {
    BuildElementsAndVariableStorage( pref, fem_manager, vset );
    InitializeConnectivity( vset );
 }



// NB: here we have to deallocate all the dynamic storage, NOT DONE YET! 
template<size_t dim>
MeshManager<dim>::~MeshManager()
 {
    if ( adaptive_remeshing_ ) {
          // doing the nodes first
//          MeshManagerAdapter<dim>::const_iterator  elmt_iterator(ElementsBegin());
/*
          while ( elmt_iterator != ElementsEnd() ) {
               for ( typename vector<Element<dim>*>::iterator
                     nit=(*elmt_iterator).NodesBegin(); nit!=(*elmt_iterator).NodesEnd(); nit++ ) {
                    delete (*nit);
                    (*nit) = 0;
                 }
               elmt_iterator++;
            }
*/            
         // now doing the elements
         throw csmp::Exception( ERROR, "MeshManager<dim>::~MeshManager:",
                                "adaptively refined mesh not implemented yet");
      }
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

    adaptive_remeshing_     = mmgr.adaptive_remeshing_;
    hybrid_element_mesh_    = mmgr.hybrid_element_mesh_;
    node_collection_        = mmgr.node_collection_;
    elmt_collection_        = mmgr.elmt_collection_;
    face_collection_        = mmgr.face_collection_;
    interface_collection_   = mmgr.interface_collection_;

    // now all pointers inside the elements and nodes must be 
    // correctly assigned to the new locations
    if ( !adaptive_remeshing_ ) {
        // element nodes and neighbor pointers
        typename deque<csmp::Element<dim> >::iterator  ite(elmt_collection_.begin());
        
        for ( typename deque<Element<dim> >::const_iterator
              eit=mmgr.elmt_collection_.begin(); eit!=mmgr.elmt_collection_.end(); eit++, ite++ ) 
          {
             // node pointers
             for ( size_t i=0U; i<(*eit).Nodes(); i++ )
               (*ite).Assign( i, &node_collection_[ (*eit).N(i)->Idx() ] );

             // neighbor pointers
             for ( size_t i=0U; i<(*eit).Neighbors(); i++ )
               if ( (*eit).Neighbor(i) != NULL )
                (*ite).Assign( i, &elmt_collection_[ (*eit).Neighbor(i)->Idx() ] );
          }

        // parent elements to nodes
        typename deque<csmp::Node<dim> >::iterator  itn(node_collection_.begin());

        for ( typename deque<csmp::Node<dim> >::const_iterator
              nit=mmgr.node_collection_.begin(); nit!=mmgr.node_collection_.end(); nit++, itn++ ) 
          for ( size_t i=0U; i<(*nit).Parents(); i++ )
            (*itn).Assign( (*nit).ParentNodeNumber(i), &elmt_collection_[ (*nit).Parent(i)->Idx() ] );
      }
    else { // the model has been build as a tree structure (adaptive_remeshing_=true)
        throw csmp::Exception( FATAL_ERROR, "MeshManager<dim>::operator=", 
                       "for graph-style MeshManageres, copy construction has not been implemented yet");
      }   
       
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
    for ( typename set<CSMP_FEM_TYPE>::const_iterator
          iit=input_etypes.begin(); iit!=input_etypes.end(); iit++ ) {
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
                  elmt_collection_.push_back( (Element<dim>(fem_manager.E(csmpElementType))) );
                  elmt_collection_[idx].ResizePropertyStorage( evars, cvars ); 
                  elmt_collection_[idx].Idx(idx);
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
                  elmt_collection_.push_back( Element<dim>( fem_manager.E(csmpElementType) ) ); 
                  elmt_collection_[idx].ResizePropertyStorage( evars, cvars );
                  elmt_collection_[idx].Idx(idx);
                  idx++;
                  first++;
               }

        // 3. construction the objects within the new memory
        // -------------------------------------------------
        const LocalVariables nvars( phys_vars.LocalVariablesAt(NODE) );

        Node<dim>  default_node;
        for ( size_t i=0U; i<vset.Vertices(); i++ ) {
             node_collection_.push_back( default_node );   
             node_collection_[i].ResizePropertyStorage(nvars);
             node_collection_[i].Idx(i);
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

    // since deques were used to store nodes and elements there is no flexibility
    adaptive_remeshing_ = false;

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
        typename deque<vector<size_t> >::const_iterator  first(vset.PlistBegin()), last(vset.PlistFacesBegin());
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
                  elmt_collection_.push_back( Element<dim>( idx++, fem_manager.E(csmpElementType),
                                                            evars, cvars, NOT ) );
                  first++;
              }
          }
        // 2.2 If there are multiple element types
        else
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
                  elmt_collection_.push_back( Element<dim>( idx++, fem_manager.E(csmpElementType),
                                                            evars, cvars, NOT ) );
                  first++;
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
          
            typename deque<vector<size_t> >::const_iterator  first(vset.PlistFacesBegin()), last(vset.PlistInterFacesBegin());
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
                  face_collection_.push_back( Face<dim>( idx++, fem_manager.E(csmpElementType),
                                                         evars, cvars ) );
                  first++;
               }
          } // end faces
          
        // 4. constructing the interfaces using the VSet nodes information
        // ---------------------------------------------------------------
        // (continuous running index 'idx' will be used)
        if ( vset.InterFaces() > 0 ) {
            assert( vset.HybridElementTypeMesh() );
            const LocalVariables evars( phys_vars.LocalVariablesAt(INTER_FACE) );
            const IntegrationPointVariables cvars( phys_vars.IntegrationPointVariablesAt(INTER_FACE) );

            typename deque<vector<size_t> >::const_iterator  first(vset.PlistInterFacesBegin()), last(vset.PlistEnd());
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
                     interface_collection_.push_back( InterFace<dim>( idx++, fem_manager.E(csmpElementType),
                                                                      evars, cvars ) );
                      first++;
                   }

            } // end interfaces


        // 5. construction of Node objects including coordinates and property storage allocation
        // -------------------------------------------------------------------------------------
        const LocalVariables nvars( phys_vars.LocalVariablesAt(NODE) );
        vector<double64>     coord(dim);

        for ( size_t i=0U; i<vset.Vertices(); ++i ) {
             for ( size_t j=0U; j<dim; ++j ) coord[j] = vset.P( j, i );
             // emplaced construction: Node( size_t idx, const Point<dim>&, const LocalVariables&, BOX_BOUNDARY=NOT );
             node_collection_.push_back( Node<dim>( i, Point<dim>(coord), nvars ) );
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

    // since deques were used to store nodes and elements there is no flexibility
    adaptive_remeshing_ = false;

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
    if      ( dim == 1U )
      for ( size_t i=0U; i<node_collection_.size(); i++ )
        node_collection_[i].x( vset.Px( i ) );
    else if ( dim == 2U )
      for ( size_t i=0U; i<node_collection_.size(); i++ ) {
           node_collection_[i].x( vset.Px( i ) );
           node_collection_[i].y( vset.Py( i ) );
        }
    else
      for ( size_t i=0U; i<node_collection_.size(); i++ ) {
           node_collection_[i].x( vset.Px( i ) );
           node_collection_[i].y( vset.Py( i ) );
           node_collection_[i].z( vset.Pz( i ) );
        }


    // --------------------------------------------------
    // 2. assigning nodes to elements using 'plist' array
    // --------------------------------------------------
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: assigning nodes to elements..."<< endl;
    for ( typename deque<Element<dim> >::iterator
          eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
      for ( size_t j=0U; j<(*eit).Nodes(); j++ )
         (*eit).Assign( j, &node_collection_[ vset.Plist((*eit).Idx(),j) ] );


    // --------------------------------------------------------------
    // 3. Assigning neighbor elements to elements (face-verts)
    // --------------------------------------------------------------
    /* SKM: THIS CODE SEGMENT BREAKS establishNeighborConnectivity() ?
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: assigning neighbors to elements..."<< endl;
      for ( typename deque<Element<dim> >::iterator
        eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
        for ( size_t j=0U; j<(*eit).Neighbors(); j++ )
          // if there is a neighbor
          if ( vset.Pfvert((*eit).Idx(),j) >= 0 )
              (*eit).Assign( j, &elmt_collection_[ static_cast<size_t>(vset.Pfvert((*eit).Idx(),j)) ] );
          else
              (*eit).Assign( j, static_cast<Element<dim>*>(NULL) );
    */

    // ---------------------------------------------------------------------
    // 4. Flagging nodes located at the model boundary
    // ---------------------------------------------------------------------
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeConnectivity: flagging boundary objects..."<< endl;
    if ( vset.BFlags() > 0 )
      {
        // nodes were initially constructed as not located at the model boundary
        for ( typename map<size_t,long64>::const_iterator
              bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++ )
          node_collection_[ (*bit).first ].AtBoundary( intToBOX_BOUNDARY( (*bit).second ) );
      }
    

    // ---------------------------------------------------------------------
    // 5. Flagging the elements using boundary flags from the nodes
    // ---------------------------------------------------------------------
    flagElementUsingNodalAtBoundaryFlags<dim>( ElementsBegin(), ElementsEnd() );
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
   for ( typename deque<Element<dim> >::const_iterator
         eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
     for ( typename vector<csmp::Node<dim>*>::const_iterator
           nit=(*eit).NodesBegin(); nit!=(*eit).NodesEnd(); nit++ )
       parent_elmts_per_node[ (*nit)->Idx() ]++;

   // reserving the memory for the parent storage and zeroing parent vector for next step
   for ( size_t i=0U; i<node_collection_.size(); i++ )
     node_collection_[i].ResizeParentStorage( parent_elmts_per_node[i] );

   // assigning the parent element information to the nodes
    for ( typename deque<Element<dim> >::iterator
          eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
      for ( size_t j=0U; j<(*eit).Nodes(); j++ )
        (*eit).N(j)->Assign( j, &(*eit) );

   return true;
 
 } // end InitializeConnectivity(VSet)



 



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

    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: building storage..."<< endl;

    // -------------------------------------
    // 1. assigning coordinates to the nodes
    // -------------------------------------
    const size_t  nodes(node_collection_.size());
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning node coordinates..."<< endl;
    if      ( dim == 1U )
      for ( size_t i=0U; i<nodes; ++i )
        node_collection_[i].x( vset.Px( i ) );
    else if ( dim == 2U )
      for ( size_t i=0U; i<nodes; ++i ) {
           node_collection_[i].x( vset.Px( i ) );
           node_collection_[i].y( vset.Py( i ) );
        }
    else
      for ( size_t i=0U; i<nodes; ++i ) {
           node_collection_[i].x( vset.Px( i ) );
           node_collection_[i].y( vset.Py( i ) );
           node_collection_[i].z( vset.Pz( i ) );
        }

    // ------------------------------------------------------------------------
    // 2. assigning nodes to elements, faces and interfaces using 'plist' array
    // ------------------------------------------------------------------------
    const typename deque<Element<dim> >::iterator elementsEnd(elmt_collection_.end());
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning nodes to elements..."<< endl;
    for ( typename deque<Element<dim> >::iterator
          eit=elmt_collection_.begin(); eit!=elementsEnd; ++eit ) {
           const size_t nodes((*eit).Nodes());
           for ( size_t j=0U; j<nodes; ++j )
             (*eit).Assign( j, &node_collection_[ vset.Plist((*eit).Idx(),j) ] );
       }
    // faces
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning nodes to faces..."<< endl;
    const typename deque<Face<dim> >::iterator facesEnd(face_collection_.end());
    for ( typename deque<Face<dim> >::iterator
          fit=face_collection_.begin(); fit!=facesEnd; ++fit ) {
           const size_t nodes((*fit).Nodes());
           for ( size_t j=0U; j<nodes; ++j )
             // assigning node indices
             (*fit).Assign( j, &node_collection_[ vset.Plist((*fit).Idx(),j) ] );
       }
    // interfaces
    if ( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning nodes to interfaces..."<< endl;
    const typename deque<InterFace<dim> >::iterator interfacesEnd(interface_collection_.end());
    for ( typename deque<InterFace<dim> >::iterator
          fit=interface_collection_.begin(); fit!=interfacesEnd; ++fit ) {
           const size_t nodes((*fit).Nodes());
           for ( size_t j=0U; j<nodes; ++j ) {
                // assigning node indices
                (*fit).Assign( j, &node_collection_[ vset.Plist((*fit).Idx(),j) ], INSIDE );
             }
       }

    // ----------------------------------------------------------------
    // 3. Assigning neighbor elements to elements, faces and interfaces
    // ----------------------------------------------------------------
    // 3.1 elements
    // ------------
    cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning neighbors to elements..."<< endl;
    for ( typename deque<Element<dim> >::iterator
          eit=elmt_collection_.begin(); eit!=elementsEnd; ++eit ) {
          const size_t neighbors((*eit).Neighbors());
          for ( size_t j=0U; j<neighbors; ++j ) {
            // if there is a neighbor (as is the case if the stored index is greater than zero)
            const int32 index(static_cast<int32>(vset.Pfvert( (*eit).Idx(), j )) );
            if (  index >= 0 )
                (*eit).Assign( j, &elmt_collection_[ static_cast<size_t>(index) ] );
            else
                (*eit).Assign( j, static_cast<Element<dim>*>(nullptr) );
        }
      }
    // 3.2 faces
    // ---------
    // connecting faces to their higher-dimensional neighbors
    const size_t elements(elmt_collection_.size());
    const size_t faces(face_collection_.size());
    for ( typename deque<Face<dim> >::iterator
          eit=face_collection_.begin(); eit!=facesEnd; ++eit ) {
         // equidimensional neighbors first
         // -------------------------------
         const size_t neighbors((*eit).Neighbors());
         for ( size_t j=0U; j<neighbors; ++j ) {
              // if there is a neighbor (as is the case if the stored index is greater than zero)
              // (e.Idx() starts with elements=first face)
              const long64 index(vset.Pfvert( (*eit).Idx(), j ));
              if (  index >= 0 ) {
                   assert( index >= elements );
                   assert( index < elements + faces ); // (-) elements because face container is numbered from 0..n-1
                   (*eit).Assign( j, &face_collection_[ static_cast<size_t>(index)-elements ] );
                }
              else (*eit).Assign( j, static_cast<Face<dim>*>(nullptr) );
           }
        // higher-dimensional neighbors
        // ----------------------------
        // (are stored in VSet 'pfverts' record after the equidimensional neighbors)
        // index of inner neighbor element i which is always there
        const long64 index1(vset.Pfvert( (*eit).Idx(), neighbors ));
        // index of outer neighbor element which may be there
        const long64 index2(vset.Pfvert( (*eit).Idx(), neighbors+1U ));
        Element<dim>* const outerElement = (index2 < 0) ? nullptr : &elmt_collection_[index2];
        // assigning inner and outer higher-dimensional neighbors
        //              inner element             outer element
        (*eit).Assign( &elmt_collection_[index1], outerElement );
      }
    // 3.3 interfaces
    // --------------
    const size_t interfaces(interface_collection_.size());
    // connecting interfaces to their higher-dimensional neighbors
    for ( typename deque<InterFace<dim> >::iterator
          eit=interface_collection_.begin(); eit!=interfacesEnd; ++eit ) {
         // equidimensional neighbors first
         const size_t neighbors((*eit).Neighbors());
         for ( size_t j=0U; j<neighbors; ++j ) {
            // if there is a neighbor (as is the case if the stored index is greater than zero)
            const long64 index(vset.Pfvert( (*eit).Idx(), j ));
            if (  index >= 0 ) {
                   assert( index >= elements + faces );
                   assert( index < elements + faces + interfaces ); // (-) because interface container is numbered from 0..n-1
                   (*eit).Assign( j, &interface_collection_[ static_cast<size_t>(index)-elements-faces ] );
              }
            else (*eit).Assign( j, static_cast<InterFace<dim>*>(nullptr) );
         }
        // higher-dimensional neighbors
        // ----------------------------
        // (the 2 sides will always be present because interfaces exist only on internal boundaries)
        (*eit).Assign( &elmt_collection_[ vset.Pfvert( (*eit).Idx(), neighbors ) ],
                       &elmt_collection_[ vset.Pfvert( (*eit).Idx(), neighbors+1U ) ] );
      }

    // ---------------------------------------------------------------------
    // 4. Flagging nodes located at the model boundary
    // ---------------------------------------------------------------------
    if( csmp_error.Verbose() )
        cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: flagging boundary objects..."<< endl;
    if ( vset.BFlags() > 0 )
      {
        // nodes were initially constructed as not located at the model boundary
        for ( typename map<size_t,long64>::const_iterator
              bit=vset.BFlagsBegin(); bit!=vset.BFlagsEnd(); bit++ )
          node_collection_[ (*bit).first ].AtBoundary( intToBOX_BOUNDARY( (*bit).second ) );
      }
    
    // ---------------------------------------------------------------------
    // 5. Flagging the elements using boundary flags from the nodes
    // ---------------------------------------------------------------------
    flagElementUsingNodalAtBoundaryFlags<dim>( ElementsBegin(), ElementsEnd() );

    // ------------------------------------------------------------------------------
    // 6. Assigning parent elements (these are the elements that share the node) and
    // their respective internal node-id numbers to the nodes
    // -------------------------------------------------------------------------------
    if ( csmp_error.Verbose() )
      cout <<"\nMeshManager<"<< dim <<">::InitializeVerifiedConnectivity: assigning parent element information to nodes..."<< endl;
    vector<size_t>  parent_elmts_per_node( node_collection_.size(), 0U );

    // counting how many parent elements each node has
    for ( typename deque<Element<dim> >::const_iterator
          eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
      for ( typename vector<csmp::Node<dim>*>::const_iterator
            nit=(*eit).NodesBegin(); nit!=(*eit).NodesEnd(); nit++ )
        parent_elmts_per_node[ (*nit)->Idx() ]++;

    // reserving the memory for the parent storage and zeroing parent vector for next step
    for ( size_t i=0U; i<node_collection_.size(); i++ )
      node_collection_[i].ResizeParentStorage( parent_elmts_per_node[i] );

    // assigning the parent element information to the nodes
     for ( typename deque<Element<dim> >::iterator
           eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
       for ( size_t j=0U; j<(*eit).Nodes(); j++ )
         (*eit).N(j)->Assign( j, &(*eit) );

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
   for ( typename deque<Element<dim> >::const_iterator
         eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ ) {
         size_t counter(0U);
         for ( typename vector<csmp::Node<dim>*>::const_iterator
               nit=(*eit).NodesBegin(); nit!=(*eit).NodesEnd(); nit++ )
           if ( (*nit)->AtBoundary() != NOT )
             counter++;
         if ( counter == (*eit).Nodes() ) {
               // cerr <<"\nMeshManager<dim>::DetectElementsWithAllNodesOnBoundary: found all-node-on-boundary element:";
               // (*eit).Out();
               belmts.insert( (*eit).Idx() );
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
    for ( typename deque<Element<dim> >::iterator
          it=elmt_collection_.begin(); it!=elmt_collection_.end(); it++ )
          if( !(*it).FV_Stencil() )
            {
              (*it).Assign( fvs_manager.Stencil( (*it).FE_Type() ) );
              (*it).ResizePropertyStorage( lvs, ipvs );
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
    node_collection_.emplace_back( node );
    return &node_collection_.back();
}

template<size_t dim>
Element<dim>* MeshManager<dim>::PushBack( Element<dim>&& elmt )
{
    elmt_collection_.emplace_back( elmt );
    return &elmt_collection_.back();
}

template<size_t dim>
Face<dim>* const MeshManager<dim>::PushBack( Face<dim>&& face )
{
    face_collection_.emplace_back( face );
    return &face_collection_.back();
}

template<size_t dim>
InterFace<dim>* MeshManager<dim>::PushBack( InterFace<dim>&& interface )
{
    interface_collection_.push_back( interface );
    return &interface_collection_.back();
}



template<size_t dim>
Node<dim>* MeshManager<dim>::PushBackIfUnique( Node<dim>&& node )
{
    typename std::deque< Node<dim> >::iterator nit = find( node_collection_.begin(), node_collection_.end(), node );
    if ( nit != node_collection_.end() )
       return ( &(*nit) );

    node_collection_.emplace_back( node );
    return &node_collection_.back();
}

template<size_t dim>
Element<dim>* MeshManager<dim>::PushBackIfUnique( Element<dim>&& elmt )
{
    typename std::deque< Element<dim> >::iterator eit = find( elmt_collection_.begin(), elmt_collection_.end(), elmt );
    if( eit != elmt_collection_.end() )
        return ( &(*eit) );

    elmt_collection_.emplace_back( elmt );
    return &elmt_collection_.back();
}

template<size_t dim>
Face<dim>* MeshManager<dim>::PushBackIfUnique( Face<dim>&& face )
{
    typename std::deque< Face<dim> >::iterator fit = find( face_collection_.begin(), face_collection_.end(), face );
    if( fit != face_collection_.end() )
        return ( &(*fit) );

    face_collection_.emplace_back( face );
    return &face_collection_.back();
}

template<size_t dim>
InterFace<dim>* MeshManager<dim>::PushBackIfUnique( InterFace<dim>&& interface )
{
    typename std::deque< InterFace<dim> >::iterator ifit = find( interface_collection_.begin(), interface_collection_.end(), interface );
    if( ifit != interface_collection_.end() )
        return ( &(*ifit) );
  
    interface_collection_.push_back( interface );
    return &interface_collection_.back();
}


/**
    Erases nodes (as long as they can be found) returning true as long as nodes remain, else false is returned
    
    @note method uses the erase(remove) paradigm to avoid that any iterators, pointers 
    and references related to the container are invalidated.
*/
template<size_t dim>
bool MeshManager<dim>::Erase( const Node<dim>& node )
{
   // making sure that the node that is erased is moved to the end of the queue first so that all other iterators stay intact
    if ( node_collection_.erase( remove( node_collection_.begin(), node_collection_.end(), node ), node_collection_.end() )
                                         == node_collection_.end() ) return false;
    return true;
}

template<size_t dim>
bool MeshManager<dim>::Erase( const Element<dim>& elmt )
{
    if ( elmt_collection_.erase( remove( elmt_collection_.begin(), elmt_collection_.end(), elmt ), elmt_collection_.end() )
                                         == elmt_collection_.end() ) return false;
    return true;
}

template<size_t dim>
bool MeshManager<dim>::Erase( const Face<dim>& face )
{
    if ( face_collection_.erase( remove( face_collection_.begin(), face_collection_.end(), face ), face_collection_.end() )
                                         == face_collection_.end() ) return false;
    return true;
}

template<size_t dim>
bool MeshManager<dim>::Erase( const InterFace<dim>& interface )
{
    if ( interface_collection_.erase( remove( interface_collection_.begin(),
                                              interface_collection_.end(), interface ),
                                      interface_collection_.end() ) == interface_collection_.end() ) return false;
    return true;
}

/**
   erases all nodes, discerning those that do not have any parent element connections.
*/
template<size_t dim>
bool MeshManager<dim>::EraseNodes()
{
    bool emptyNodes(false);
    typename std::deque<Node<dim> >::iterator nit = node_collection_.begin();
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
    typename std::deque<Element<dim> >::iterator eit = elmt_collection_.begin();
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
    typename std::deque<Face<dim> >::iterator fit = face_collection_.begin();
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
    typename std::deque<InterFace<dim> >::iterator ifit = interface_collection_.begin();
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
    assert( !adaptive_remeshing_ );

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
       for ( typename deque<Element<dim> >::const_iterator
             eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); ++eit ) {
            nodes_per_element.push_back((*eit).Nodes());
            neighbors_per_element.push_back((*eit).Neighbors());
            csp_fem_types.push_back((*eit).FE_Type());
         }
      
       // 1.2 adding Face information after the elements
       for ( typename deque<Face<dim> >::const_iterator
             eit=face_collection_.begin(); eit!=face_collection_.end(); ++eit ) {
            nodes_per_element.push_back((*eit).Nodes());
            neighbors_per_element.push_back( (*eit).Neighbors() + higherDimParents );
            csp_fem_types.push_back((*eit).FE_Type());
         }
      
       // 1.3 adding InterFace information after the faces
       for ( typename deque<InterFace<dim> >::const_iterator
             eit=interface_collection_.begin(); eit!=interface_collection_.end(); ++eit ) {
            // multiplier takes care of the multiplicated interface nodes that the InterFace will be connected to
            nodes_per_element.push_back( (*eit).Nodes() * interfaceMultiplier );
            neighbors_per_element.push_back( (*eit).Neighbors() * interfaceMultiplier + higherDimParents );
            csp_fem_types.push_back((*eit).FE_Type());
         }
       // 1.4 resizing the VSet
       vset.Resize( csp_fem_types, nodes_per_element, neighbors_per_element,
                    node_collection_.size(), Faces(), InterFaces() );
      }
    else { // if there is only a single element type
        assert( elmt_collection_[0].FE() != nullptr );
        vset.Resize( elmt_collection_[0].FE()->Nodes(), 
                     elmt_collection_[0].FE()->Neighbors(),
                     elmt_collection_[0].FE()->ElementType(),
                     node_collection_.size(), elmt_collection_.size() );
                     
        vset.ElementType( 0, elmt_collection_[0].FE_Type() );
      }
      
    // 2. adding node coordinates and boundary flags (BOX_BOUNDARY)
    // ------------------------------------------------------------
    size_t i(0U);
    const typename deque<Node<dim> >::const_iterator nodesEnd(node_collection_.end());
    if ( dim == 1U )
      for( typename deque<Node<dim> >::const_iterator
           nit=node_collection_.begin(); nit!=nodesEnd; ++nit++, ++i )
        vset.Px(i,(*nit).x()); 
    else if ( dim == 2U )
      for( typename deque<Node<dim> >::const_iterator
           nit=node_collection_.begin(); nit!=nodesEnd; ++nit, ++i ) {
           vset.Px(i,(*nit).x()); 
           vset.Py(i,(*nit).y());
        }
    else
      for( typename deque<Node<dim> >::const_iterator
           nit=node_collection_.begin(); nit!=nodesEnd; ++nit, ++i ) {
           vset.Px(i,(*nit).x()); 
           vset.Py(i,(*nit).y());
           vset.Pz(i,(*nit).z());
        }

    // 3. adding 'plist' connectivity list
    // -----------------------------------
    // elements
    const size_t  elements(elmt_collection_.size());
    for ( size_t eidx=0U; eidx<elements; ++eidx )
      for ( size_t j=0U; j<elmt_collection_[eidx].Nodes(); ++j )
        vset.Plist( eidx, j, (elmt_collection_[eidx].N(j)->Idx()) );
    // faces
    const size_t faces(face_collection_.size());
    for ( size_t eidx=0U; eidx<faces; ++eidx )
      for ( size_t j=0U; j<face_collection_[eidx].Nodes(); ++j )
        vset.Plist( elements + eidx, j, (face_collection_[eidx].N(j)->Idx()) );
    // interfaces
    const size_t interfaces(interface_collection_.size());
    for ( size_t eidx=0U; eidx<interfaces; ++eidx )
      for ( size_t j=0U; j<interface_collection_[eidx].Nodes(); ++j )
        vset.Plist( elements + faces + eidx, j, (interface_collection_[eidx].N(j)->Idx()) );
   
   
    // 4. adding 'pfverts' neighbors per element list
    // ----------------------------------------------
    // 'pfverts' elements
    for ( size_t eidx=0U; eidx<elements; ++eidx )
      for ( size_t j=0U; j<elmt_collection_[eidx].Neighbors(); ++j ) {
             Element<dim>* const ptr(elmt_collection_[eidx].Neighbor(j));
             if ( ptr != nullptr ) vset.Pfvert( eidx, j, static_cast<int32>(ptr->Idx()) );
             else                  vset.Pfvert( eidx, j, elmt_collection_[eidx].AtBoundary() );
        }
    // 'pfverts' faces
    // ---------------
    // the faces are stored after the elements including connections to their higher-dimensional neighbors
    // add the end of the pfverts entries
    if ( faces > 0 )
      {
         for ( size_t fidx=0U; fidx<faces; ++fidx ) {
             // equidimensional neighbors first
             const size_t neighbors(face_collection_[fidx].Neighbors());
             for ( size_t j=0U; j<neighbors; ++j ) {
                    Face<dim>* const ptr(face_collection_[fidx].Neighbor(j));
                    // if the neighbor exists (which it must on the inside of the Face)
                    if ( ptr != nullptr ) {
                         vset.Pfvert( fidx + elements, j, static_cast<int32>(ptr->Idx()) );
                      }
                    else vset.Pfvert( fidx + elements, j, face_collection_[fidx].InnerParent()->AtBoundary() );
               }
             // higher-dimensional neighbors second
             // inner neighbor
             if      ( dim == 3U ) assert( face_collection_[fidx].InnerParent()->IsVolumeElement() );
             else if ( dim == 2U ) assert( face_collection_[fidx].InnerParent()->IsSurfaceElement() );
             assert( face_collection_[fidx].InnerParent()->Idx() < elements );
             vset.Pfvert( fidx + elements, neighbors, static_cast<int32>(face_collection_[fidx].InnerParent()->Idx()) );
             // outer neighbor
             if      ( face_collection_[fidx].OuterParent() != nullptr && dim == 3U ) assert( face_collection_[fidx].OuterParent()->IsVolumeElement() );
             else if ( face_collection_[fidx].OuterParent() != nullptr && dim == 2U ) assert( face_collection_[fidx].OuterParent()->IsSurfaceElement() );
             if ( face_collection_[fidx].OuterParent() != nullptr ) {
                  assert( face_collection_[fidx].OuterParent()->Idx() < elements );
                  vset.Pfvert( fidx + elements, neighbors + 1U, static_cast<int32>(face_collection_[fidx].OuterParent()->Idx()) );
               }
             else {
                  // if there is no neighbor, the inner element parent should be at the model boundary
                  if ( face_collection_[fidx].InnerParent()->AtBoundary() == NOT ) {
                       face_collection_[fidx].Out();
                       csmp_error.notice( WARNING, "MeshManager<dim>::OutputMeshTo (face neighbors):",
                                         "inner dim+1 element should be at model boundary because Face has is no outer element.");
                       face_collection_[fidx].InnerParent()->AtBoundary( IRREGULAR );
                    }
                  vset.Pfvert( fidx + elements, neighbors + 1U, face_collection_[fidx].InnerParent()->AtBoundary() );
               }
          }
      }
    // 'pfverts' interfaces
    // --------------------
    if ( interfaces > 0 )
      {
         const size_t interfaces(interface_collection_.size());
         for ( size_t fidx=0U; fidx<interfaces; ++fidx ) {
             // equidimensional neighbors (=other interfaces) first
// TODO: each side will have neighbors on the separated sides of the interface; track!
             const size_t neighbors(interface_collection_[fidx].Neighbors());
             for ( size_t j=0U; j<neighbors; ++j ) {
                    InterFace<dim>* const ptr(interface_collection_[fidx].Neighbor(j));
                    if ( ptr != nullptr ) {
                         vset.Pfvert( fidx + elements + faces, j, static_cast<int32>(ptr->Idx()) );
                      }
                    else vset.Pfvert( fidx + elements + faces, j, REGION_BOUNDARY );
               }
             // inner and outer higher-dimensional neighbors
             // ( they must always exist because SplitBoundaries are internal model boundaries)
             assert( interface_collection_[fidx].InnerParent() != nullptr );
             assert( interface_collection_[fidx].OuterParent() != nullptr );
             vset.Pfvert( fidx + elements + faces, neighbors, static_cast<int32>(interface_collection_[fidx].InnerParent()->Idx()) );
             vset.Pfvert( fidx  + elements + faces, neighbors + 1U, static_cast<int32>(interface_collection_[fidx].OuterParent()->Idx()) );
          }
      }
   
    // 5. adding boundary flags
    // ------------------------
    for( typename deque<Node<dim> >::const_iterator
         nit=node_collection_.begin(); nit!=nodesEnd; ++nit )
      if ( (*nit).AtBoundary() != NOT ) vset.AddBFlag( (*nit).Idx(), (*nit).AtBoundary() );
         
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
        const size_t elmt_ips(elmt_collection_[0].IntegrationPoints()); // just an estimate
       
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
        const size_t elmt_sector_ips(elmt_collection_[0].IntegrationPointsPerSector());
        const size_t sectors_per_element(elmt_collection_[0].Sectors());
       
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
        const size_t elmt_facet_ips(elmt_collection_[0].IntegrationPointsPerFacet());
        const size_t facets_per_element(elmt_collection_[0].Facets());

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


    // face integration point properties
    // ---------------------------------
    if ( database.ListProperties( FACE_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !face_collection_.empty() );
        const size_t face_ips(face_collection_[0].IntegrationPoints());
       
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
    if ( database.ListProperties( FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !face_collection_.empty() );
        const size_t face_sector_ips(face_collection_[0].IntegrationPointsPerSector());
        const size_t sectors_per_face(face_collection_[0].Sectors());

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
    if ( database.ListProperties( FACE_FACET_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !face_collection_.empty() );
        const size_t face_facet_ips(face_collection_[0].IntegrationPointsPerFacet());
        const size_t facets_per_face(face_collection_[0].Facets());

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

    // interface integration point properties
    // --------------------------------------
    if ( database.ListProperties( INTER_FACE_INTEGRATION_POINT, properties ) > 0 ) {
        const size_t interface_ips(interface_collection_[0].IntegrationPoints());
       
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
    if ( database.ListProperties( INTER_FACE_SECTOR_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !interface_collection_.empty() );
        const size_t interface_sector_ips(interface_collection_[0].IntegrationPointsPerSector());
        const size_t sectors_per_interface(interface_collection_[0].Sectors());

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
    if ( database.ListProperties( INTER_FACE_FACET_INTEGRATION_POINT, properties ) > 0 ) {
        assert( !interface_collection_.empty() );
        const size_t interface_facet_ips(face_collection_[0].IntegrationPointsPerFacet());
        const size_t facets_per_interface(face_collection_[0].Facets());

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
         assert( (*pit).second.Size() / key.dataDepth == node_collection_.size() );
        
         const size_t nodes(node_collection_.size());
         // for the given property type
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    for ( size_t i=0U; i<nodes; ++i ) {
                         read( (*pit).second, i, value );
                         node_collection_[i].Store( key, value );
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    for ( size_t i=0U; i<nodes; ++i ) {
                         read( (*pit).second, i, value );
                         node_collection_[i].Store( key, value );
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    for ( size_t i=0U; i<nodes; ++i ) {
                         read( (*pit).second, i, value );
                         node_collection_[i].Store( key, value );
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    for ( size_t i=0U; i<nodes; ++i ) {
                         read( (*pit).second, i, value );
                         node_collection_[i].Store( key, value );
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    for ( size_t i=0U; i<nodes; ++i ) {
                         read( (*pit).second, i, value );
                         node_collection_[i].Store( key, value );
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

         const size_t elements(elmt_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    for ( size_t i=0U; i<elements; ++i ) {
                          read( (*pit).second, i, value );
                          elmt_collection_[i].Store( key, value );
                       }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    for ( size_t i=0U; i<elements; ++i ) {
                          read( (*pit).second, i, value );
                          elmt_collection_[i].Store( key, value );
                       }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    for ( size_t i=0U; i<elements; ++i ) {
                          read( (*pit).second, i, value );
                          elmt_collection_[i].Store( key, value );
                       }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    for ( size_t i=0U; i<elements; ++i ) {
                          read( (*pit).second, i, value );
                          elmt_collection_[i].Store( key, value );
                       }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    for ( size_t i=0U; i<elements; ++i ) {
                          read( (*pit).second, i, value );
                          elmt_collection_[i].Store( key, value );
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

         const size_t elements(elmt_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U); // running index
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t integration_points(elmt_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              elmt_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t integration_points(elmt_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              elmt_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U); 
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t integration_points(elmt_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              elmt_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U); 
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t integration_points(elmt_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              elmt_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U); 
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t integration_points(elmt_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              elmt_collection_[k].Store( i, key, value );
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

         const size_t elements(elmt_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t sectors(elmt_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(elmt_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t sectors(elmt_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(elmt_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t sectors(elmt_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(elmt_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t sectors(elmt_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(elmt_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t sectors(elmt_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(elmt_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
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

         const size_t elements(elmt_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t facets(elmt_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(elmt_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t facets(elmt_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(elmt_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t facets(elmt_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(elmt_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t facets(elmt_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(elmt_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<elements; ++k ) {
                         const size_t facets(elmt_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(elmt_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   elmt_collection_[k].Store( i, j, key, value );
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

         const size_t faces(face_collection_.size());
         // for the given property type
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    for ( size_t i=0U; i<faces; ++i ) {
                         read( (*pit).second, i, value );
                         face_collection_[i].Store( key, value );
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    for ( size_t i=0U; i<faces; ++i ) {
                         read( (*pit).second, i, value );
                         face_collection_[i].Store( key, value );
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    for ( size_t i=0U; i<faces; ++i ) {
                         read( (*pit).second, i, value );
                         face_collection_[i].Store( key, value );
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    for ( size_t i=0U; i<faces; ++i ) {
                         read( (*pit).second, i, value );
                         face_collection_[i].Store( key, value );
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    for ( size_t i=0U; i<faces; ++i ) {
                         read( (*pit).second, i, value );
                         face_collection_[i].Store( key, value );
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

         const size_t faces(face_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t integration_points(face_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              face_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t integration_points(face_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              face_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t integration_points(face_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              face_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t integration_points(face_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              face_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t integration_points(face_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              face_collection_[k].Store( i, key, value );
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

         const size_t faces(face_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t sectors(face_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(face_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t sectors(face_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(face_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t sectors(face_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(face_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t sectors(face_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(face_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t sectors(face_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(face_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
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

         const size_t faces(face_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t facets(face_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(face_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t facets(face_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(face_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t facets(face_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(face_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t facets(face_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(face_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<faces; ++k ) {
                         const size_t facets(face_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(face_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   face_collection_[k].Store( i, j, key, value );
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
        
         const size_t interfaces(interface_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    for ( size_t i=0U; i<interfaces; ++i ) {
                         read( (*pit).second, i, value );
                         interface_collection_[i].Store( key, value );
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    for ( size_t i=0U; i<interfaces; ++i ) {
                         read( (*pit).second, i, value );
                         interface_collection_[i].Store( key, value );
                      }
                }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    for ( size_t i=0U; i<interfaces; ++i ) {
                         read( (*pit).second, i, value );
                         interface_collection_[i].Store( key, value );
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    for ( size_t i=0U; i<interfaces; ++i ) {
                         read( (*pit).second, i, value );
                         interface_collection_[i].Store( key, value );
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    for ( size_t i=0U; i<interfaces; ++i ) {
                         read( (*pit).second, i, value );
                         interface_collection_[i].Store( key, value );
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
        
         const size_t interfaces(interface_collection_.size());
         switch( key.type )
          {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t integration_points(interface_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              interface_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t integration_points(interface_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              interface_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t integration_points(interface_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              interface_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t integration_points(interface_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              interface_collection_[k].Store( i, key, value );
                              entry++;
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t integration_points(interface_collection_[k].IntegrationPoints());
                         for ( size_t i=0U; i<integration_points; ++i ) {
                              read( (*pit).second, entry, value );
                              interface_collection_[k].Store( i, key, value );
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
        
         const size_t interfaces(interface_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t sectors(interface_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(interface_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t sectors(interface_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(interface_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t sectors(interface_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(interface_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t sectors(interface_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(interface_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t sectors(interface_collection_[k].Sectors());
                         for ( size_t i=0U; i<sectors; ++i ) {
                              const size_t ips_per_sector(interface_collection_[k].IntegrationPointsPerSector());
                              for ( size_t j=0U; j<ips_per_sector; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
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
        
         const size_t interfaces(interface_collection_.size());
         switch( key.type )
           {
             case SCALAR: {
                    ScalarVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t facets(interface_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(interface_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case VECTOR: {
                    VectorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t facets(interface_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(interface_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case TENSOR: {
                    TensorVariable<dim> value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t facets(interface_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(interface_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case ARRAY: {
                    ArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t facets(interface_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(interface_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
                                   entry++;
                                }
                           }
                      }
                 }
               break;
             case FLAGGEDARRAY: {
                    FlaggedArrayVariable value;
                    size_t entry(0U);
                    for ( size_t k=0U; k<interfaces; ++k ) {
                         const size_t facets(interface_collection_[k].Facets());
                         for ( size_t i=0U; i<facets; ++i ) {
                              const size_t ips_per_facet(interface_collection_[k].IntegrationPointsPerFacet());
                              for ( size_t j=0U; j<ips_per_facet; ++j ) {
                                   read( (*pit).second, entry, value );
                                   interface_collection_[k].Store( i, j, key, value );
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
void MeshManager<dim>::Out() const
 {
    cout <<"\nMeshManager<"<< dim <<">::Out: "<< endl;
    
    // nodes
    cout <<"\nNODES: "<< endl;
    for ( typename deque<Node<dim> >::const_iterator
          nit=node_collection_.begin(); nit!=node_collection_.end(); nit++ )
      { 
         string bound = parseBoundary((*nit).AtBoundary());
         cout <<"\nNode ID: "<< (*nit).Idx() <<" ";
         cout << (*nit).Coordinate();
         cout <<" Boundary flag: "<< bound << endl;
      }

    // elements
    cout <<"\nELEMENTS: "<< endl;
    for ( typename deque<Element<dim> >::const_iterator
          eit=elmt_collection_.begin(); eit!=elmt_collection_.end(); eit++ )
      { 
         string bound = parseBoundary((*eit).AtBoundary());
         cout <<"\nElement ID: "<< (*eit).Idx() <<" Boundary flag: "<< bound << endl;
         cout <<"Member Nodes: "<< endl;
         for ( size_t i=0U; i<(*eit).Nodes(); i++ ) 
           cout << (*eit).N(i)->Idx() <<"\t";
          cout <<"\nNeighbor elements: "<< endl;
         for ( size_t i=0U; i<(*eit).Neighbors(); i++ )
           if ( (*eit).Neighbor(i) != NULL )
             cout << (*eit).Neighbor(i)->Idx() <<"\t";
           else
           cout <<"NO NEIGHBOR\t";        
 
         cout << endl;
     }

   // faces
   // TODO:

   // inter faces
   // TODO:

   // parent elements ID's for each node
   cout << endl << endl;   
    cout <<"PARENT ELEMENT INFORMATION FOR ALL NODES: "<< endl;
   for ( typename deque<Node<dim> >::const_iterator
         nit=node_collection_.begin(); nit!=node_collection_.end(); nit++ )
     {
        cout <<"\nNode: "<< (*nit).Idx() <<", parent elements: "<< endl;
        for ( size_t i=0u; i<(*nit).Parents(); i++ )
          cout << (*nit).Parent( i )->Idx() <<" ";
        cout << endl;
     } 

 } // end Out
 
 
 
template class MeshManager<1U>;
template class MeshManager<2U>;
template class MeshManager<3U>;



} // end namespace csmp 















