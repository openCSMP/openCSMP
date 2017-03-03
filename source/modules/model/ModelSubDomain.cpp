/*
 *  ModelSubDomain.cpp
 *
 *  Created by Stephan Matthai on 7/18/10.
 *  Copyright 2010 SKM private. All rights reserved.
 *
 */
#include <type_traits>
#include "ModelSubDomain.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Box.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "FiniteElementManager.h"
#include "MeshManager.h"
#include "FiniteVolumeStencilManager.h"

#include "Visitor.h"
#include "PropertyDatabase.h"
#include "Interrelation.h"
#include "CopyReplaceVisitor.h"

#include "Exception.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"
#include "binaryReadWrite.h"

//#define MODEL_SUBDOMAIN_DEBUG

using namespace std;

namespace csmp {

template<size_t dim, template<size_t> class SIMPLEX>
ModelSubDomain<dim,SIMPLEX>::ModelSubDomain( const string& subdomain_name,
                                             const PropertyDatabase<dim>& pref )
 :  pref_(pref),
    subdomain_name_(subdomain_name),
    verbose_(true)
 {
 }


template<size_t dim, template<size_t> class SIMPLEX>
ModelSubDomain<dim,SIMPLEX>::ModelSubDomain( const ModelSubDomain& ed )
 : pref_(ed.pref_),
   elmt_vec_(ed.elmt_vec_),
   node_vec_(ed.node_vec_),
   first_bd_node_(ed.first_bd_node_),
   bd_face_vec_(ed.bd_face_vec_),
   subdomain_name_(ed.subdomain_name_),
   verbose_(true)
 {
   this->LVS( ed.LVS() );
 }


/// move constructor; @attention remove verbose output after testing
template<size_t dim, template<size_t> class SIMPLEX>
ModelSubDomain<dim,SIMPLEX>::ModelSubDomain( ModelSubDomain&& ed )
 : pref_{ed.pref_},
   elmt_vec_{ed.elmt_vec_},
   node_vec_{ed.node_vec_},
   first_bd_node_{ed.first_bd_node_},
   bd_face_vec_{ed.bd_face_vec_},
   subdomain_name_{ed.subdomain_name_},
   verbose_{true}
 {
   this->LVS( move(ed.LVS()) );
   cout <<"\nModelSubDomain<dim,SIMPLEX>::ModelSubDomain: called MOVE constructor.\n";
 }





/* DOES NOT WORK BECAUSE "All Elements" region gets sorted before

template<size_t dim, template<size_t> class SIMPLEX>
ModelSubDomain<dim,SIMPLEX>::ModelSubDomain( const PropertyDatabase<dim>& pref,
                                             const ModelSubDomain<dim,SIMPLEX>& mesh, "All Elements"
                                             const SubDomainInfo& info )
 : pref_(pref),
   first_bd_node_(info.interior_nodes.size()),
   subdomain_name_(info.name),
   verbose_(true)
 {
    // building the element vector
    // ---------------------------
    elmt_vec_.reserve( info.interior_elmts.size() + info.perimeter_elmts.size() );
    // assigning pointers to the interior elements
    for ( auto it=info.interior_elmts.begin(); it!=info.interior_elmts.end(); ++it )
      elmt_vec_.push_back( mesh.elmt_vec_[ (*it) ] );
    // assigning pointers to the perimeter elements
    for ( auto it=info.perimeter_elmts.begin(); it!=info.perimeter_elmts.end(); ++it )
      elmt_vec_.push_back( mesh.elmt_vec_[ (*it) ] );
    // sorting the subvectors for future searching
    const auto perimeterElementsBegin( next(elmt_vec_.begin(), info.interior_elmts.size()) );
    sort( elmt_vec_.begin(), perimeterElementsBegin );
    sort( perimeterElementsBegin, elmt_vec_.end() );

    // building the vector of vectors of those faces of the simplices that lie on the subdomain perimeter
    // --------------------------------------------------------------------------------------------------
    bd_face_vec_.reserve( info.perimeter_faces.size() );
    //std::vector<std::vector<int8 > > perimeter_faces
    for ( auto it=info.perimeter_faces.begin(); it!=info.perimeter_faces.end(); ++it ) {
          const size_t perimeter_faces((*it).size());
          std::vector<ONE_BYTE_NUMBER> face_vec;
          face_vec.reserve(perimeter_faces);
          for ( auto fit=(*it).begin(); fit!=(*it).end(); ++fit )
            face_vec.push_back( static_cast<ONE_BYTE_NUMBER>( (*fit) ) );
          this->bd_face_vec_.emplace_back( face_vec );
      }

    // building the node vector
    // ------------------------
    node_vec_.reserve( info.interior_nodes.size() + info.perimeter_nodes.size() );
    // assigning pointers to the interior nodes
    for ( auto it=info.interior_nodes.begin(); it!=info.interior_nodes.end(); ++it )
      node_vec_.push_back( mesh.node_vec_[ (*it) ] );

    // assigning pointers to the perimeter nodes
    this->first_bd_node_ = info.interior_nodes.size();
    for ( auto it=info.perimeter_nodes.begin(); it!=info.perimeter_nodes.end(); ++it )
      node_vec_.push_back( mesh.node_vec_[ (*it) ] );

    // sorting the subvectors for future searching
    const auto perimeterNodesBegin( next(node_vec_.begin(), info.interior_nodes.size()) );
    sort( node_vec_.begin(), perimeterNodesBegin );
    sort( perimeterNodesBegin, node_vec_.end() );

    // allocating the storage for subdomain properties
    // -----------------------------------------------
    this->ResizePropertyStorage( pref.LocalVariablesAt( parsePlacement<dim,SIMPLEX>() ) );
 
 } // end constructor
*/



template<size_t dim, template<size_t> class SIMPLEX>
ModelSubDomain<dim,SIMPLEX>::~ModelSubDomain()
 {
 }


template<size_t dim, template<size_t> class SIMPLEX>
ModelSubDomain<dim,SIMPLEX>&  ModelSubDomain<dim,SIMPLEX>::operator=( const ModelSubDomain& ed )
 {
     if ( &ed != this ) {
          elmt_vec_       = ed.elmt_vec_;
          node_vec_       = ed.node_vec_;
          first_bd_node_  = ed.first_bd_node_;
          bd_face_vec_    = ed.bd_face_vec_;
          verbose_        = ed.verbose_;
          subdomain_name_ = ed.subdomain_name_;
          this->LVS( ed.LVS() );
       }
     return *this;
 }

template<size_t dim, template<size_t> class SIMPLEX>
string  ModelSubDomain<dim,SIMPLEX>::Name() const
 {
    return subdomain_name_;
 }


template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::Name( const string& name )
 {
    subdomain_name_ = name;
 }


template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::Verbose( bool verbose )
{
    this->verbose_ = verbose;
}

template<size_t dim, template<size_t> class SIMPLEX>
bool ModelSubDomain<dim,SIMPLEX>::Verbose()
{
    return this->verbose_;
}


/**
    Helper functions that checks complex boundary flags of a variable that shall be assigned.
    They functions only transfer values to it if the variable has not got the specified flag.
    
    Generic version for variables that are placed on the Node, Element, Face etc.
    
    @TODO: make this part of new PDE_Intetgrator class
    
    @author SKM 10/9/2014
*/
template<size_t dim, template<size_t> class SIMPLEX, class Var>
void writeVariableIf( SIMPLEX<dim>*, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
 } // end generic specification


/// write guard for scalar variables
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             const csmp::Index& idx,
                             const ScalarVariable& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(idx) != dont_overwrite )
      ptr->Store( idx, var );
 } // end version for scalars
 

/// write guard for vector variables
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             const csmp::Index& idx,
                             const VectorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( idx, vc );
    for ( size_t i=0U; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             const csmp::Index& idx,
                             const TensorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( idx, ts );
    for ( size_t i=0U; i<dim; i++ )
      if ( ptr->Status(idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( size_t j=0U; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( idx, ts );
 } // end version for tensors



/**
    Helper functions that checks a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
    
    Generic version for variables that are placed on Element/Face/Interface integration points.
*/
template<size_t dim, template<size_t> class SIMPLEX, class Var>
void writeVariableIf( SIMPLEX<dim>*, size_t ip, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
 } // end generic specification

/// write guard for scalar variables
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t ip,
                             const csmp::Index& idx,
                             const ScalarVariable& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for scalars
 

/// write guard for vector variables
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t ip,
                             const csmp::Index& idx,
                             const VectorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( ip, idx, vc );
    for ( size_t i=0U; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( ip, idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t ip,
                             const csmp::Index& idx,
                             const TensorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( ip, idx, ts );
    for ( size_t i=0U; i<dim; i++ )
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( size_t j=0U; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( ip, idx, ts );
 } // end version for tensors



/**
    Helper functions that checks a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
    
    Generic version for finite volume-related integration points.
*/
template<size_t dim, template<size_t> class SIMPLEX, class Var>
void writeVariableIf( SIMPLEX<dim>*, size_t sector_or_facet,
                      size_t ip, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
 } // end generic specification

/// write guard for scalar variables
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t sector_or_facet,
                             size_t ip,
                             const csmp::Index& idx,
                             const ScalarVariable& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for scalars
 

/// write guard for vector variables
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t sector_or_facet,
                             size_t ip,
                             const csmp::Index& idx,
                             const VectorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( sector_or_facet, ip, idx, vc );
    for ( size_t i=0U; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( sector_or_facet, ip, idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<size_t dim, template<size_t> class SIMPLEX>
void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t sector_or_facet,
                             size_t ip,
                             const csmp::Index& idx,
                             const TensorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( sector_or_facet, ip, idx, ts );
    for ( size_t i=0U; i<dim; i++ )
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( size_t j=0U; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( sector_or_facet, ip, idx, ts );
 } // end version for tensors



// inlined methods

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::Nodes() const
  {
     return node_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::InteriorNodes() const
  {
     return first_bd_node_;
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::PerimeterNodes() const
  {
     return node_vec_.size() - InteriorNodes();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::Elements() const
  {
     return elmt_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::InteriorElements() const
  {
     return elmt_vec_.size() - bd_face_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::PerimeterElements() const
  {
     return bd_face_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
bool ModelSubDomain<dim,SIMPLEX>::Empty() const
  {
     return elmt_vec_.empty();
  }


/// returns how many faces of the target element lie on the subdomain boundary
template<size_t dim, template<size_t> class SIMPLEX>
size_t  ModelSubDomain<dim,SIMPLEX>::PerimeterFaces( size_t e ) const
 {
    assert( e >= InteriorElements() );
    assert( e < elmt_vec_.size() );
    return bd_face_vec_[e-InteriorElements()].size();
 }

/// returns the elements local face number of the n'th face that is on the subdomain boundary
template<size_t dim, template<size_t> class SIMPLEX>
size_t  ModelSubDomain<dim,SIMPLEX>::PerimeterFace( size_t e, size_t face ) const
 {
    assert( e >= InteriorElements() );
    assert( e < elmt_vec_.size() );
    assert( face < PerimeterFaces(e) );
    return static_cast<size_t>(bd_face_vec_[e-InteriorElements()][face]);
 }


template<size_t dim, template<size_t> class SIMPLEX>
csmp::Node<dim>*  ModelSubDomain<dim,SIMPLEX>::N( size_t nd ) const
 { assert( nd < node_vec_.size() ); return node_vec_[nd]; }


template<size_t dim, template<size_t> class SIMPLEX>
SIMPLEX<dim>*  ModelSubDomain<dim,SIMPLEX>::E( size_t e ) const
 { assert( e < elmt_vec_.size() ); return elmt_vec_[e]; }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::SimplexVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::ElementVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::FaceVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::InterFaceVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<Node<dim>*>&  ModelSubDomain<dim,SIMPLEX>::NodeVector()
  { return node_vec_; }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::NodesBegin() const
 { return node_vec_.begin(); }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterNodesBegin() const
 { return std::next( node_vec_.begin(), InteriorNodes() ); }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::NodesEnd() const
 { return node_vec_.end(); }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::ElementsBegin() const
 { return elmt_vec_.begin(); }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterElementsBegin() const
 { return std::next( elmt_vec_.begin(), InteriorElements() ); }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::ElementsEnd() const
 { return elmt_vec_.end(); }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<csmp::Node<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::NodesBegin()
 { return node_vec_.begin(); }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<csmp::Node<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterNodesBegin()
 { return std::next( node_vec_.begin(), InteriorNodes() ); }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<csmp::Node<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::NodesEnd()
 { return node_vec_.end(); }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::ElementsBegin()
 { return elmt_vec_.begin(); }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterElementsBegin()
 { return std::next( elmt_vec_.begin(), InteriorElements() ); }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::ElementsEnd()
 { return elmt_vec_.end(); }




/**
    Detects of how many spatial dimensions element types are contained in model.
    It returns a pair: first value gives number of different spatial dimensions contained,
    second value returns the highest spatial dimension contained.

    @author SKM 1/11/2013

*/
template<size_t dim, template<size_t> class SIMPLEX>
pair<int32,int32>  ModelSubDomain<dim,SIMPLEX>::SpatialDimensions() const
 {
    bool with_volume_elements(false);
    bool with_surface_elements(false);
    bool with_line_elements(false);

    for( typename vector<SIMPLEX<dim>*>::const_iterator
         it=this->elmt_vec_.begin(); it!=this->elmt_vec_.end(); it++ ) {
         if      ( (*it)->IsLineElement() )    with_line_elements = true;
         else if ( (*it)->IsSurfaceElement() ) with_surface_elements = true;
         else if ( (*it)->IsVolumeElement() )  with_volume_elements = true;
      }

    int32 counter(0);
    if ( with_volume_elements )  counter++;
    if ( with_surface_elements ) counter++;
    if ( with_line_elements )    counter++;
    int32 highest_spatial_dim(1);
    if      ( with_volume_elements )  highest_spatial_dim = 3;
    else if ( with_surface_elements ) highest_spatial_dim = 2;

    return make_pair( counter, highest_spatial_dim );

 } // end ElementSpatialDimensions








/**
    Connects the simplices with their equidimensional neighbors
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::EstablishNeighborConnectivity()
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( elmt_vec_.empty() ) {
         csmp_error.notice( WARNING, "ModelSubDomain<dim,SIMPLEX>::EstablishNeighborConnectivity:", "supplied cell vector is empty; nothing was done." );
         return;
      }
    cout << "\nModelSubDomain<"<< dim <<",SIMPLEX>::EstablishNeighborConnectivity: Establishing CSMP FE neighbor connectivity...\n";
 
    // 1. making separate search vectors of face keys for surface and line elements
    // ----------------------------------------------------------------------------
    cout << "  Building a list of the faces of the cells...\n";
    //       key             face number,neighbor
    multimap<set<Node<dim>*>,pair<size_t,SIMPLEX<dim>*> >  volume_neighbor_keys,
                                                           surface_neighbor_keys, 
                                                           line_neighbor_keys;
    vector<size_t>                 fnids;
    typename std::set<Node<dim>*>  key; // region, boundary and split boundary all use nodes

    for ( typename vector<SIMPLEX<dim>*>::const_iterator it=elmt_vec_.begin(); it!=elmt_vec_.end(); ++it ) {
          const size_t faces((*it)->Faces());
          for ( size_t face=0U; face<faces; ++face )
            {
               if ( (*it) == nullptr ) {
                    csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::EstablishNeighborConnectivity:",
                                      "supplied element contains NULL pointer to elements; nothing was done." );
                    return;
                 }
               // creating face key of node pointers from indices of face nodes
               (*it)->FE()->NodesOfFace( face, fnids );
               const size_t nodes(fnids.size());
               for ( size_t j=0U; j<nodes; ++j )
                 key.insert( (*it)->N( fnids[j] ) );
                 
               // inserting newly generated keys into multimap
               // TODO: use move constructor here
               if ( (*it)->IsVolumeElement() )
                 volume_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
               else if ( (*it)->IsSurfaceElement() )
                 surface_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
               else // for all line elements
                 line_neighbor_keys.insert( make_pair( key, make_pair( face, (*it) ) ) );
               key.clear();
            }
        }

    // 2. (re)building element neigborhoods
    // ------------------------------------
    // (the assumption here is that adjacent neighbors are arranged consecutively in the multimap)
    cout << "  (Re)building neighbor connectivity...";

    // 2.1 line elements
    // -----------------
    if ( !line_neighbor_keys.empty() )
      {
        SIMPLEX<dim>* e1Ptr(nullptr);
        SIMPLEX<dim>* e2Ptr(nullptr);

        cout << "\n\t\tline elements...";
        //                key                  n-face, neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,SIMPLEX<dim>*> >::iterator it1(line_neighbor_keys.begin()),
                                                                                 it2(line_neighbor_keys.begin());
        it2++;

        while ( it2 != line_neighbor_keys.end() )
          { 
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first )
                {
                   assert( (*it1).second.second != nullptr );
                   assert( (*it2).second.second != nullptr );
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr ); // avoid self-assignment

                   // assigning the two cells face neighbors to one another
                   //                        face pointer  nbor face idx   neighbor pointer
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   
                   // both iterators are advanced (so that with the second increment a new pair of faces is reached)
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == line_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // line elements
    
    // 2.2 surface elements
    // --------------------
    if ( dim >= 2U and !surface_neighbor_keys.empty() )
      {
        cout << "\n\t\tsurface elements...";
        //                key              n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,SIMPLEX<dim>*> >::iterator it1(surface_neighbor_keys.begin()),
                                                                                 it2(surface_neighbor_keys.begin());
        it2++;

        while ( it2 != surface_neighbor_keys.end() )
          { 
              // if there is a pair of valid neighbor elements, neighbor assignments are made
              if ( (*it1).first == (*it2).first )
                {
                    assert( (*it1).second.second != nullptr );
                    assert( (*it2).second.second != nullptr );
                    SIMPLEX<dim>* e1Ptr((*it1).second.second);
                    SIMPLEX<dim>* e2Ptr((*it2).second.second);
                    if ( e1Ptr != e2Ptr ) {
                         ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                         ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                      }
                    else csmp_error.notice( WARNING, "ModelSubDomain<dim,SIMPLEX>::EstablishNeighborConnectivity:",
                                            "discovered potentially duplicate surface simplex.");
                    ++it1;
                    ++it2;
                }
              if ( it2 == surface_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // surface elements
      
    // 2.3 volume elements
    // -------------------
    if ( dim == 3U and !volume_neighbor_keys.empty() ) {

        SIMPLEX<dim>* e1Ptr(nullptr);
        SIMPLEX<dim>* e2Ptr(nullptr);

        cout << "\n\t\tvolume elements...\n";
        //                key             n-face neighbor
        typename multimap<set<Node<dim>*>,pair<size_t,SIMPLEX<dim>*> >::iterator it1(volume_neighbor_keys.begin()),
                                                                                 it2(volume_neighbor_keys.begin());
        it2++;

        while ( it2 != volume_neighbor_keys.end() )
          { 
              // if there is a pair of valid neighbor elements, neighbor assignments are made
//              if ( (*it1).first == (*it2).first and ( (*it1).second.second != nullptr and (*it2).second.second != nullptr ) )
              if ( (*it1).first == (*it2).first )
                {
                   assert( (*it1).second.second != nullptr );
                   assert( (*it2).second.second != nullptr );
                   e1Ptr = (*it1).second.second;
                   e2Ptr = (*it2).second.second;
                   assert( e1Ptr != e2Ptr );
                   ((*it1).second.second)->Assign( (*it1).second.first, e2Ptr );
                   ((*it2).second.second)->Assign( (*it2).second.first, e1Ptr );
                   ++it1;
                   ++it2;
                }

              // both iterators are advanced
              if ( it2 == volume_neighbor_keys.end() ) break;
              ++it1;
              ++it2;
          } 
      } // dim=3
    
 } // end EstablishNeighborConnectivity



/* TESTING - EstablishNeighborConnectivity

/ printing the multimap
RenumberElements();
cerr <<"\nline element face key map:\n";
cerr <<"\n\tnode-id, face-id, nbor elmt id, neighbor 1 and 2";
for ( auto it=line_neighbor_keys.begin(); it!=line_neighbor_keys.end(); ++it ) {
      cerr <<"\n\t"<< (*(*it).first.begin())->Idx() <<", "<< (*it).second.first <<", ";
      if ( (*it).second.second != nullptr ) {
           cerr << (*it).second.second->Idx() <<", ";
           if ( (*it).second.second->Neighbor(0) != nullptr )
             cerr << (*it).second.second->Neighbor(0)->Idx() <<", ";
           else cerr <<"nullptr" <<", ";
           if ( (*it).second.second->Neighbor(1) != nullptr )
             cerr << (*it).second.second->Neighbor(1)->Idx() <<", ";
           else cerr <<"nullptr" <<", ";
        }
      else cerr <<"nullptr.";
   }
cerr << endl;
cerr <<"\nprinting the elements:";
for ( auto it=elmt_vec_.begin(); it!=elmt_vec_.end(); ++it )
  (*it)->Out();
cerr << endl;
      
*/






/**
    In order to use the region node/element flags 'INTERIOR' or 'PERIMETER',
    boundary nodes and elements must be identified first.
    This is done every time when a new region is formed.

    @attention convention: elements are considered 'PERIMETER' elements
    of a Region only if they have at least one edge or face at the
    region boundary. This is the case irrespective of the dimensionality
    of such elements. In the extreme case of line elements inside a volume,
    they are flagged boundary, if they have one node on the perimeter of
    the volume region.
    Another way of looking at this is to consider all elements boundary
    elements, that share a face with the boundary.

    @attention convention: If a Region consists of elements of different
    dimensionality, the highest dimensional ones define the position of the boundary, i.e.
    all lower-dimensional elements that stick out of the region (and including all their
    nodes) are flagged 'PERIMETER'.
    Lower-dimensional elements inside a higher dimensional region are considered
    'INTERIOR'.

    @section application Application

    IdentifyPerimeter() is called as part of the region-forming process.
    The perimeter flagging is used to access boundaries selectively.

    @todo (1) Does not work for quadratic regions

   @section implementation Implementation

   1. Perimeter elements and nodes are found looking for the absence
      of element neighbors, but only the highest dimensionality elements
      inside the Region

      - all the nodes of the highest-dimensionality elements are stored
        in a set.

   2. Lower-dimensional elements are found that do not share all of their nodes
      with the highest dimensional elements.

      - These elements nodes must be outside of the highest dimensional region.
      - These elements and their nodes are flagged 'PERIMETER'.

   3. All lower dimensional elements that share all of their nodes with the
      highest dimensional elements AND have at least nodes flagged 'PERIMETER'
      on one edge (a single node for a line element)
      will be flagged 'PERIMETER' elements.

   It follows that all lower-dimensional mesh that 'sticks out' of a higher
   dimensional region has elements and nodes flagged as PERIMETER.

*/
template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::IdentifyPerimeter()
 {
    // 1. putting the perimeter elements at the end of the elmt_vec and sorting
    //    interior and perimeter cell ranges subsequently
    // ----------------------------------------------------
    const size_t first_bd_elmt = this->PartitionCellVector();
    assert( first_bd_elmt <= this->elmt_vec_.size() );

 } // end IdentifyPerimeter










/**

  Remarks on binary search:

  - The sorted source range referenced must be valid; all pointers must be dereferenceable and,
    within the sequence, the last position must be reachable from the first by incrementation.

  - The sorted range must each be arranged as a precondition to the application of the binary_search
    algorithm in accordance with the same ordering as is to be used by the algorithm
    to sort the combined ranges.

    @attention any lower-dimensional elements and their nodes that stick outside of a higher
    dimensional region will be flagged as boundary.
*/
template<size_t dim, template<size_t> class SIMPLEX>
size_t  ModelSubDomain<dim,SIMPLEX>::PartitionCellVector()
 {
    ErrorHandler&  csmp_error(ErrorHandler::Instance());

    if ( this->elmt_vec_.empty() )
      throw logic_error( (string("ModelSubDomain<dim>::PartitionCellVector: method called on empty region: ") + Name()).c_str() );

    sort( this->elmt_vec_.begin(), this->elmt_vec_.end() );

    // --------------------------------------------------------------------
    // 0. detecting whether this region contains lower-dimensional elements
    // --------------------------------------------------------------------
    // pair contains: 1) number of spatial element dimensions in region, 2)  highest contained dimension
    const pair<int32,int32>  elmt_dim = SpatialDimensions();

    // -----------------------------------------------------------------
    // 1. distinguishing boundary from interior elements, same for nodes
    //   (at this point the elements and nodes are already known)
    // -----------------------------------------------------------------
    set<SIMPLEX<dim>*> interior_elmts, boundary_elmts;
    set<Node<dim>*>                   boundary_nodes;
    set<pair<SIMPLEX<dim>*,size_t> >  boundary_faces;
    vector<size_t>  fnids;

    // 1.1 If all elements have the same spatial dimension
    // ---------------------------------------------------
    if ( elmt_dim.first == 1 )
      {
        for ( typename vector<SIMPLEX<dim>*>::const_iterator
              eit=this->elmt_vec_.begin(); eit!=this->elmt_vec_.end(); eit++ )
           {
             // identifying the boundary faces and their nodes
             // (each face potentially has a neighbor element)
             long  nbors_that_belong_to_group((*eit)->Neighbors());
             for ( size_t i=0U; i<(*eit)->Faces(); i++ )
               // if the face is at a model boundary or has a neighbor that does not belong to the region
               if ( (*eit)->Neighbor(i) == NULL  or  !binary_search( this->elmt_vec_.begin(), this->elmt_vec_.end(), (*eit)->Neighbor(i) ) )
                 {
                    // boundary faces
                    boundary_faces.insert( make_pair( (*eit), i ) );
                    // boundary nodes
                    assert( (*eit)->FE() != NULL );
                    (*eit)->FE()->NodesOfFace( i, fnids );
                    for ( size_t j=0U; j<fnids.size(); ++j )
                      boundary_nodes.insert( (*eit)->N(fnids[j]) );
                    // counting neighbors
                    nbors_that_belong_to_group--;
                 }

             // storing the distinguished elements in the respective vectors
             // ------------------------------------------------------------
             // interior elements
             if ( nbors_that_belong_to_group == (*eit)->Neighbors() )
                 interior_elmts.insert( (*eit) );
             // elements with at least one face on the region boundary
             else
                 boundary_elmts.insert( (*eit) );
          }
       }

    // 1.2 If there are elements with different spatial dimensions
    // -----------------------------------------------------------
    //     the ones with highest dimensions are used to define perimeter
    //     all lower dimensional mesh that sticks out is flagged as perimeter as well.
    else
      {
        // a. identify the boundary elements among the highest dimensional elements,
        //    also collecting all their node pointers into a set.
        set<SIMPLEX<dim>*> lesser_dim_elmts;
        set<Node<dim>*>    highest_dim_elmt_nodes;

        for ( typename vector<SIMPLEX<dim>*>::const_iterator
              eit=this->elmt_vec_.begin(); eit!=this->elmt_vec_.end(); eit++ )
          {
              // elements of the highest spatial dimension are used to define the boundary
              assert( parseFiniteElementDimension( (*eit)->FE_Type() ) != 0 );
              if ( parseFiniteElementDimension( (*eit)->FE_Type() ) == elmt_dim.second )
                {
                   // creating a subset with their nodes
                   for ( size_t i=0U; i<(*eit)->Nodes(); ++i ) {
                        assert( (*eit)->N(i) != NULL );
                        highest_dim_elmt_nodes.insert( (*eit)->N(i) );
                     }
                   // if the element has faces that lie on the region boundary
                   // it is considered a boudary element
                   long  nbors_that_belong_to_group((*eit)->Neighbors());
                   for ( size_t i=0U; i<(*eit)->Faces(); ++i )
                     // 1) the element is on model boundary  or  2) one of its neighbors does not belong to its parent region
                     if ( (*eit)->Neighbor(i) == NULL        or  !binary_search( this->elmt_vec_.begin(), this->elmt_vec_.end(), (*eit)->Neighbor(i) ) )
                       {
                          // the element pointer and the face number are used to create a unique key for the discovered boundary face
                          boundary_faces.insert( make_pair( (*eit), i ) );
                          // recording the nodes of the boundary face as boundary nodes
                          (*eit)->FE()->NodesOfFace( i, fnids );
                          for ( size_t j=0U; j<fnids.size(); ++j )
                            boundary_nodes.insert( (*eit)->N(fnids[j]) );

                          nbors_that_belong_to_group--;
                       }
                   if ( nbors_that_belong_to_group == (*eit)->Neighbors() ) interior_elmts.insert( (*eit) );
                   else boundary_elmts.insert( (*eit) );
                }
               else lesser_dim_elmts.insert( (*eit) );
            }
        assert( /* all elements are accounted for */ this->elmt_vec_.size() == interior_elmts.size() + boundary_elmts.size() + lesser_dim_elmts.size() );

#ifdef MODEL_SUBDOMAIN_DEBUG
cout <<"\nModelSubDomain<" << dim <<">::PartitionCellVector: '"<< this->Name() <<"' highest subdomain element dim: "<< elmt_dim.second <<"\n";
cout <<"\n\tlesser-dim elements: "<< lesser_dim_elmts.size() <<", boundary elements: "<< boundary_elmts.size();
cout <<"\n\ttotal nodes: "<< this->node_vec_.size() <<", highest-dim element nodes: "<< highest_dim_elmt_nodes.size() << endl;
cout.flush();
#endif

         // ---------------------------------------------------------------------------
         // 1.3  processing lower dimensional elements and their nodes in the subdomain
         // ---------------------------------------------------------------------------
         // 1.3.1 nodes that are not contained in the higher-dimensional element subset are identified as extra boundary node
         for ( typename set<SIMPLEX<dim>*>::const_iterator it=lesser_dim_elmts.begin(); it!=lesser_dim_elmts.end(); ++it )
           for ( size_t i=0U; i<(*it)->Nodes(); ++i )
             if ( highest_dim_elmt_nodes.find( (*it)->N(i) ) == highest_dim_elmt_nodes.end() )
               boundary_nodes.insert( (*it)->N(i) );

         // 1.3.2 finding the lesser dimensional elements on the region boundary
         set<SIMPLEX<dim>*> lesser_dim_elmts_detached; // to distinguish stand-alone lower dimensional mesh

         for ( typename set<SIMPLEX<dim>*>::const_iterator it=lesser_dim_elmts.begin(); it!=lesser_dim_elmts.end(); ++it )
           {
              // a) lower-dim elements sticking out
              // ----------------------------------
              // lower-dimensional elements with nodes that do not belong to the node set of the
              // higher dimensional elements must be boundary elements
              size_t  exterior_nodes(0U);
              for ( size_t i=0U; i<(*it)->Nodes(); ++i )
                if ( highest_dim_elmt_nodes.find( (*it)->N(i) ) == highest_dim_elmt_nodes.end() )
                  exterior_nodes++;

              // if individual nodes stick out the parent element sticks out as well.
              if ( exterior_nodes >= 1U ) {
                   // adding boundary elements and boundary faces
                   for ( size_t i=0U; i<(*it)->Faces(); ++i )
                     if ( (*it)->Neighbor(i) == nullptr  or
                          lesser_dim_elmts.find( static_cast<SIMPLEX<dim>*>((*it)->Neighbor(i)) ) == lesser_dim_elmts.end() ) {
                          // the element is a boundary element that sticks out of the region
                          boundary_elmts.insert( (*it) );
                          boundary_faces.insert( make_pair( (*it), i ) );
                       }
                   // if the entire element sticks out of the higher dimensional domain, it is saved for a warning issued later
                   if ( exterior_nodes == (*it)->Nodes() )
                     lesser_dim_elmts_detached.insert( (*it) );
                }
              // b) lower-dim element that shares all their nodes with the higher dimensional ones
              // ---------------------------------------------------------------------------------
              // are boundary elements if they have at least one face on the subdomain boundary:
              // - for line elements this means at least one node
              // - for surface elements this means at least one edge
              else {
                   // line elements (assuming that the faces correspond to the nodes)
                   if ( (*it)->IsLineElement() ) {
                        for ( size_t i=0U; i<(*it)->Nodes(); ++i )
                          // if the node is a boundary noode
                          if ( boundary_nodes.find( (*it)->N(i) ) != boundary_nodes.end() ) {
                               boundary_elmts.insert( (*it) );
                               // boundary faces for line elements
                               boundary_faces.insert( make_pair( (*it), i ) );
                            }
                     }
                   // surface elements (this will only be possible in a 3D model) are on the boundary
                   // if they share at least one face with it
                   else {
                        for ( size_t i=0U; i<(*it)->Faces(); ++i ) {
                             (*it)->FE()->NodesOfFace( i, fnids );
                            size_t  bnodes(0U);
                            for ( size_t j=0U; j<fnids.size(); ++j )
                              if ( boundary_nodes.find( (*it)->N(fnids[j]) ) != boundary_nodes.end() )
                                bnodes++;
                            // if all the nodes of at least one face lie at the boundary, so does the element
                            if ( bnodes == fnids.size() ) {
                                 boundary_elmts.insert( (*it) );
                                 // boundary faces
                                 boundary_faces.insert( make_pair( (*it), i ) );
                              }
                          }
                     }
                }
           }

         // adding lesser-dimensional elements that are not at the boundary to the interior domain
         for ( typename set<SIMPLEX<dim>*>::const_iterator it=lesser_dim_elmts.begin(); it!=lesser_dim_elmts.end(); ++it )
           if ( boundary_elmts.find( (*it) ) == boundary_elmts.end() )
             interior_elmts.insert( (*it) );

         if ( !lesser_dim_elmts_detached.empty() ) {             
              csmp_error.notice( WARNING, "ModelSubdomain<dim,SIMPLEX>::PartitionCellVector:", Name().c_str(),
                                "subdomain contains lower-dimensional elements detached from higher dimensional domain; these will be treated as boundary.");
              // do some additional diagnostics on these elements
              // ------------------------------------------------
              cerr <<"\n\tdetached elements: "<< lesser_dim_elmts_detached.size() <<":";
              for ( typename set<SIMPLEX<dim>*>::const_iterator
                    it=lesser_dim_elmts_detached.begin(); it!=lesser_dim_elmts_detached.end(); ++it ) cerr <<" "<< (*it)->Idx();
              cerr << endl;
           }

#ifdef MODEL_SUBDOMAIN_DEBUG
cout <<"\n\ttotal elements: "<< this->elmt_vec_.size() <<", interior ones: "<< interior_elmts.size() <<", boundary elements: "<< boundary_elmts.size();
cout <<", stand-alone lower-dim elements: "<< lesser_dim_elmts_detached.size();
cout <<" (sum="<< interior_elmts.size() + boundary_elmts.size() <<").";
cout <<"\n\ttotal nodes: "<< this->node_vec_.size() <<", boundary nodes: "<< boundary_nodes.size() << endl << endl;
cout.flush();
#endif

       } // end multi-dim element region


    // --------------------------------------------------
    // 2. rebuilding the element vector
    // --------------------------------------------------
//    assert( interior_elmts.size() + boundary_elmts.size() == this->elmt_vec_.size() );
    // appending the boundary element vector<double64> to the interior element vector
    this->elmt_vec_.clear();
    this->elmt_vec_.assign( interior_elmts.begin(), interior_elmts.end() );
    back_insert_iterator<vector<SIMPLEX<dim>*> >  back_it(this->elmt_vec_);
    copy( boundary_elmts.begin(), boundary_elmts.end(), back_it );

    // swap trick to trim excess memory from end of vector
    vector<SIMPLEX<dim>*>( this->elmt_vec_ ).swap( this->elmt_vec_ );


#ifdef MODEL_SUBDOMAIN_DEBUG
// TESTING - is there an element with a boundary face that is not in the boundary element vector and vice versa
bool no_error_yet(true);
set<SIMPLEX<dim>*> elmts_with_bfaces;
for ( typename set<pair<SIMPLEX<dim>*,size_t> >::const_iterator it=boundary_faces.begin(); it!=boundary_faces.end(); ++it ) {
      if ( boundary_elmts.find( (*it).first ) == boundary_elmts.end() ) {
           if ( no_error_yet ) {
                cerr <<"\nboundary face parent elements vx. boundary elements:\n";
                no_error_yet=false;
             }
           cerr <<" "<< (*it).first->Idx() <<": "<< parseFiniteElementType( (*it).first->FE_Type() );
        }
      elmts_with_bfaces.insert( (*it).first );
    }
if ( elmts_with_bfaces.size() != boundary_elmts.size() )
  cerr <<"\n\tmismatch between boundary face-  and boundary element set: "<< elmts_with_bfaces.size() <<" vs. "<< boundary_elmts.size() << endl;
assert( elmts_with_bfaces.size() == boundary_elmts.size() );
#endif

    // --------------------------------------------------
    // 3. creating the boundary face vector
    // --------------------------------------------------
    if ( !this->bd_face_vec_.empty() ) this->bd_face_vec_.clear();
    this->bd_face_vec_.reserve( this->elmt_vec_.size() - boundary_elmts.size() );
    //       parent element of face, face
    typename set<pair<SIMPLEX<dim>*,size_t> >::const_iterator  bfit( boundary_faces.begin() );
    typename set<pair<SIMPLEX<dim>*,size_t> >::const_iterator  ffit( boundary_faces.begin() );
    vector<ONE_BYTE_NUMBER>  bface_data;
    bface_data.reserve(3);
    size_t counter(0U);

    while ( bfit != boundary_faces.end() ) {
         assert( (*ffit).first == this->elmt_vec_[counter + interior_elmts.size()] );
         while ( (*bfit).first == (*ffit).first ) {
              bface_data.push_back( static_cast<ONE_BYTE_NUMBER>( (*bfit).second ) );
              bfit++;
              if ( bfit == boundary_faces.end() )
                  break;
           }
         ffit = bfit;
         counter++;
         this->bd_face_vec_.push_back( bface_data );
         bface_data.clear();
      }
    vector<vector<ONE_BYTE_NUMBER> >( this->bd_face_vec_ ).swap( this->bd_face_vec_ );

    assert( this->bd_face_vec_.size() == this->Elements() - this->InteriorElements() );


    // --------------------------------------------------
    // 4. partitioning the node vector
    // --------------------------------------------------
    this->first_bd_node_ = this->node_vec_.size() - boundary_nodes.size();
    // rebuilding and sorting the node vector (set nodes are already sorted)
    vector<csmp::Node<dim>*>  temp;
    temp.reserve(this->node_vec_.size());
    // first, the interior elements are inserted
    for ( typename vector<csmp::Node<dim>*>::const_iterator
          nit=this->node_vec_.begin(); nit!=this->node_vec_.end(); nit++ )
      if ( boundary_nodes.find(*nit) == boundary_nodes.end() )
        temp.push_back( *nit );
    // now the vector is sorted
    sort( temp.begin(), temp.end() );
    // second, the already sorted boundary nodes are appended
    for ( typename set<csmp::Node<dim>*>::const_iterator
          nit=boundary_nodes.begin(); nit!=boundary_nodes.end(); nit++ )
        temp.push_back( *nit );
    // now the temporary vector is assigned to the permanent one
    this->node_vec_ = temp;

#ifdef MODEL_SUBDOMAIN_DEBUG
cout <<"\nModelSubDomain<dim,SIMPLEX>::EstablishNeighborConnectivity: '"<< this->Name() <<"': of the ";
cout << this->elmt_vec_.size() <<" elements, "<< boundary_elmts.size() <<" lie at the model boundary."<< endl;
cout.flush();
#endif

    return this->elmt_vec_.size() - boundary_elmts.size();

 } // end PartitionElementVector









// PLACEMENT OF PROPERTY ASSIGNMENT

SUBDOMAIN_PART parseSubdomainPart( const char* subdomain )
{
    std::string ssubdomain(subdomain);

    std::transform(ssubdomain.begin(),ssubdomain.end(),ssubdomain.begin(),::toupper);

    if (  ssubdomain == "COMPLETE"  )  return COMPLETE;
    if (  ssubdomain == "INTERIOR"  )  return INTERIOR;
    if (  ssubdomain == "PERIMETER" )  return PERIMETER;
    return COMPLETE;
}


std::string parseSubdomainPart( SUBDOMAIN_PART ssubdomain )
 {
    if (  ssubdomain == COMPLETE  )  return std::string("COMPLETE");
    if (  ssubdomain == INTERIOR  )  return std::string("INTERIOR");
    if (  ssubdomain == PERIMETER )  return std::string("PERIMETER");
    return std::string("COMPLETE");
}



// INTERRELATIONS INTERFACE

/** applies Interrelation to ModelSubDomain
*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::Apply( Interrelation<dim>& relation )
 {
    relation.Apply( *this );

 } // end Apply


// VISITORS INTERFACE


/// virtual function stub that will be overwritten by base classes
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::Accept( csmp::Visitor<dim>& )
 {
    throw csmp::Exception( ERROR,
                           "ModelSubDomain<dim,SIMPLEX>::Accept",
                           "Subclass method should be called; nothing was done");

 } // end Accept


// INTEGRATION POINTS

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::IntegrationPoints() const
 {
    if ( !(*elmt_vec_.begin())->FE()->UsesLocalCoordinates() ) return 0U;

    size_t  cpoints(0U);
    for( typename vector<SIMPLEX<dim>*>::const_iterator
         it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
      cpoints += (*it)->IntegrationPoints();

    return cpoints;
 }


template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::SectorIntegrationPoints() const
  {
  if ( !(*elmt_vec_.begin())->FE()->UsesLocalCoordinates() ) return 0U;

  size_t  cpoints(0U);
  for( typename vector<SIMPLEX<dim>*>::const_iterator
    it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
    cpoints += (*it)->Sectors()*(*it)->IntegrationPointsPerSector();

  return cpoints;
  }

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::FacetIntegrationPoints() const
  {
  if ( !(*elmt_vec_.begin())->FE()->UsesLocalCoordinates() ) return 0U;

  size_t  cpoints(0U);
  for( typename vector<SIMPLEX<dim>*>::const_iterator
    it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
    cpoints += (*it)->Facets()*(*it)->IntegrationPointsPerFacet();

  return cpoints;
  }




// INDEXES



/**

Returns a vector<double64> with the ID numbers of the Elements which belong
to the Region.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::MemberElementIndexes( vector<size_t>& ids ) const
 {
    ids.clear();
    ids.reserve( elmt_vec_.size() );

    for ( typename vector<SIMPLEX<dim>*>::const_iterator
          it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ ) ids.push_back( (*it)->Idx() );
 }

/** Renumbers nodes from 0 to n-1.
*/

template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::RenumberNodes() const
 {
    size_t  counter(0U);

    for ( const auto it : node_vec_ ) it->Idx( counter++ );

    return counter;
 }

/** Renumbers elements from 0 to n-1.
*/
template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::RenumberElements() const
 {
    size_t counter(0U);

    for( const auto it : elmt_vec_ ) it->Idx(counter++);

    return counter;
 } // end RenumberElements


/** Renumbers elements and nodes from 0 to n-1.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::UpdateMemberIndexes() const
 {
    RenumberNodes();
    RenumberElements();
   
 } // end UpdateRegionMemberIndexes







// ACCESSORY

///// looks for simplex with provided index
//template<size_t dim, template<size_t> class SIMPLEX>
//bool ModelSubDomain<dim,SIMPLEX>::Contains( size_t eidx ) const
//  {
//    const typename vector<SIMPLEX<dim>*>::const_iterator simplexEnd( this->elmt_vec_.end() );
//    for( typename vector<SIMPLEX<dim>*>::const_iterator simplex( this->elmt_vec_.begin() ); simplex != simplexEnd; ++simplex )
//      if( (*simplex)->Idx() == eidx )
//        return true;
//    return false;
//  } // end

template<size_t dim, template<size_t> class SIMPLEX>
bool ModelSubDomain<dim,SIMPLEX>::Contains( const SIMPLEX<dim>* e ) const
 {
    if ( binary_search( elmt_vec_.begin(), elmt_vec_.begin() + static_cast<long>(InteriorElements()), e ) )
       return true;

    if ( binary_search( elmt_vec_.begin() + static_cast<long>(InteriorElements()), elmt_vec_.end(), e ) )
       return true;

    return false;

 } // end


template<size_t dim, template<size_t> class SIMPLEX>
bool ModelSubDomain<dim,SIMPLEX>::Contains( const Node<dim>* nptr ) const
 {
    assert( nptr != NULL );

    if ( binary_search( node_vec_.begin(), node_vec_.begin() + static_cast<long>(InteriorNodes()), nptr ) )
       return true;

    if ( binary_search( node_vec_.begin() + static_cast<long>(InteriorNodes()), node_vec_.end(), nptr ) )
       return true;

    return false;

 } // end

template<size_t dim, template<size_t> class SIMPLEX>
bool ModelSubDomain<dim,SIMPLEX>::IsPerimeterNode( const csmp::Node<dim>* nd_ptr ) const
 {
    assert( nd_ptr != NULL );
    return std::binary_search( PerimeterNodesBegin(), NodesEnd(), nd_ptr );
 }


template<size_t dim, template<size_t> class SIMPLEX>
bool  ModelSubDomain<dim,SIMPLEX>::IsPerimeterElement( const SIMPLEX<dim>* e_ptr ) const
 {
    assert( e_ptr != NULL );
    return std::binary_search( PerimeterElementsBegin(), ElementsEnd(), e_ptr );
 }


template<size_t dim, template<size_t> class SIMPLEX>
bool  ModelSubDomain<dim,SIMPLEX>::IsPerimeterNode( size_t i ) const
 {
    if ( i >= Nodes() ) return false;
    return ( i >=  first_bd_node_ );
 }

/**
    Returns true when the queried element index (0..n-1)
    is among those of the elements located on the boundary of the region;
    else false.
*/
template<size_t dim, template<size_t> class SIMPLEX>
bool  ModelSubDomain<dim,SIMPLEX>::IsPerimeterElement( size_t e ) const
 {
    if ( e >= elmt_vec_.size() ) return false;
    return (e < InteriorElements()) ? false : true;
 }



// GEOMETRY


template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const
 {
    xyz_min = xyz_max = (*node_vec_.begin())->Coordinate();

    // only nodes at the group boundary have to be checked
    for ( typename vector<csmp::Node<dim>*>::const_iterator
          bit=PerimeterNodesBegin(); bit!=NodesEnd(); bit++ )
      {
         csmp::Point<dim> p = (*bit)->Coordinate();
         xyz_min[0] = std::min( p[0], xyz_min[0] );
         xyz_max[0] = std::max( p[0], xyz_max[0] );
         xyz_min[1] = std::min( p[1], xyz_min[1] );
         xyz_max[1] = std::max( p[1], xyz_max[1] );
         xyz_min[2] = std::min( p[2], xyz_min[2] );
         xyz_max[2] = std::max( p[2], xyz_max[2] );
      }

 } // end MinMaxCoordinates

/**

AssignElementCharacteristicsTo() allows to assign a number of Element
characteristics as identified by strings (second argument) to scalar physical
variables. These characteristics are:

"volume" (3D)
"inner radius"
"aspect ratio" (longest / shortest segment)

@section arguments Input Arguments

The element characteristic indicated by the first string argument is
assigned to a property of the user's choice specified by the second method
argument and as identified by its name in the variable database.

@section implementation Implementation

The method was implemented mainly for two-dimensional calculations.
Therefore, height and width only work in two-dimensional calculations.
The characteristics are obtained from the FiniteElement which is
bridged to the current Element that is queried.

@section application Application

The characteristic "inner radius" is commonly used to find the
appropriate resolution for an advection or visualization grid on which a
variable is to be mapped on.

AssignElementCharacteristicsTo() can be used to:
- test the shape of elements in a mesh for their suitability for a
  computation (aspect ratio).
- testing how skewed the elements got by deformation
- integrate element properties in specific calculations, for instance,
  if a fracture is just one element wide, the minimum height of
  these fracture elements may be equivalent to the fracture aperture.

@section messages Messages

Firstly, the method will report an error and return, if the target
property to which the element characteristic shall be assigned to is not
an element variable.

Further errors may be reported if the user tries to use the characteristic
area in a 3D computation, or volume in a 2D computation. Also, element
height and width are thus far only available in 2D.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::AssignElementCharacteristicsTo( const char* characteristic,
                                                                  const char* var )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
     csmp::Index     prop_key = pref_.StorageKey(var);
     ScalarVariable  sc;

     if ( prop_key.place != ELEMENT and prop_key.type != SCALAR ) {
          csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignElementCharacteristicsTo",
                                    "method assigns only to scalar element properties");
          return;
       }
     if ( dim != 1U and !strcmp( characteristic, "length" ) )
       csmp_error.notice( WARNING, "ModelSubDomain<dim,SIMPLEX>::AssignElementCharacteristicsTo",
                                   "'length' will only be assigned to line elements, nothing done to others");

     if ( dim != 2U and !strcmp( characteristic, "area" ) )
       csmp_error.notice( WARNING, "ModelSubDomain<dim,SIMPLEX>::AssignElementCharacteristicsTo",
                                   "'area' will only be assigned to surface elements, nothing done to others");

     if ( dim != 3U and !strcmp( characteristic, "volume" ) ) {
          csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignElementCharacteristicsTo",
                                    "'volume' can only be assigned to volume elements. Nothing was done");
          return;
       }

     if ( !strcmp( characteristic, "length" ) )
       for ( typename vector<SIMPLEX<dim>*>::iterator
             eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
         {
             if ( (*eit)->FE()->IsLineElement() ) {
                  sc = fabs( (*eit)->Volume() );
                  (*eit)->Store( prop_key, sc );
               }
         }
     else if ( !strcmp( characteristic, "area" ) )
       for ( typename vector<SIMPLEX<dim>*>::iterator
             eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
         {
             if ( (*eit)->FE()->IsSurfaceElement() ) {
                  sc = fabs( (*eit)->Volume() );
                  (*eit)->Store( prop_key, sc );
               }
         }
     else if ( !strcmp( characteristic, "volume" ) )
       for ( typename vector<SIMPLEX<dim>*>::iterator
             eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
         {
             if ( (*eit)->FE()->IsVolumeElement() ) {
                  sc = fabs( (*eit)->Volume() );
                  (*eit)->Store( prop_key, sc );
               }
         }
     else if ( !strcmp( characteristic, "aspect ratio" ) )
       {
          for ( typename vector<SIMPLEX<dim>*>::iterator
                eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
            {
               sc = (*eit)->AspectRatio();
               (*eit)->Store( prop_key, sc );
            }
      }
     else if ( !strcmp( characteristic, "inner radius" ) )
       {
          for ( typename vector<SIMPLEX<dim>*>::iterator
                eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
            {
               sc = (*eit)->InnerRadius();
               (*eit)->Store( prop_key, sc );
            }
      }
     else {
          csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignElementCharacteristicsTo",
                             characteristic,  " entered as characteristic was not identified; nothing done");
          cout <<"\nYour options are: "<< endl;
          cout <<"\n\t inner radius"   << endl;
          cout <<"\t aspect ratio"     << endl;
          cout <<"\t length" << endl;
          cout <<"\t area (>=2D models)" << endl;
          cout <<"\t volume (3D models only)" << endl;
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignElementCharacteristicsTo", "Now exciting");
      }

 } // end





/**

Assigns node coordinates to the supplied node variable which may be either
of scalar or vector<double64> type. If either the X, Y, or Z coordinate are to be
assigned to a scalar variable an index from 0...dim-1 must be supplied
to indicate which coordinate axis shall be mapped.

@section arguments Input Arguments

The name of the physical variable to which the node coordinate shall
be assigned to, and, if a scalar property shall be initialized with
coordinate values, an int index to indicate which coordinate axis is
desired (0=X, 1=Y, 2=Z).

@section implementation Implementation

Method uses the Coordinates() and X(), Y(), and Z() interfaces of the
Node class.

@section application Application

To arrive at coordinate values for a model where topography is
important etc.

@section messages Messages

The method will report a fatal error and terminate the computation, if
the target variable is not placed on the node or if it is not of the
correct type.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo(  const char* vector_variable )
 {
     csmp::Index          prop_key = pref_.StorageKey(vector_variable);
     VectorVariable<dim>  vc;

      if ( prop_key.type != VECTOR )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo",
                                    "Node coordinates can only be assigned to a vector<double64> variable");

     if ( prop_key.place != NODE )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo",
                                    "Node coordinates can only be assigned to a node variable");

     for ( typename vector<csmp::Node<dim>*>::iterator
           nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
          vc = (*nit)->Coordinate();
          (*nit)->Store( prop_key, vc );
       }

 } // end AssignNodeCoordinatesTo



template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo( const char* scalar_variable, char c )
 {
     csmp::Index     prop_key = pref_.StorageKey(scalar_variable);
     ScalarVariable  sc;

     if ( c != 'x' && c != 'y' && c != 'z' &&
          c != 'X' && c != 'Y' && c != 'Z')
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo",
                                    "coordinate index letter is not valid",
                                    "should be either of 'x,y,z,X,Y,Z'");

     if ( prop_key.type != SCALAR )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo",
                                    "Individual node coordinates can only be assigned to a scalars");

     if ( prop_key.place != NODE )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo",
                                    "Node coordinates can only be assigned to a node variable");

     if ( dim == 2U and c == 'z' )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo",
                                    "2D models have no z coordinates");

     if ( dim == 1U and c != 'x' )
       throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::AssignNodeCoordinatesTo",
                                    "1D models only have x coordinates");

     for ( typename vector<csmp::Node<dim>*>::iterator
           nit=NodesBegin(); nit!=NodesEnd(); nit++ )
       {
         if      ( c == 'x' || c == 'X' ) sc = (*nit)->x();
         else if ( c == 'y' || c == 'Y' ) sc = (*nit)->y();
         else if ( c == 'z' || c == 'Z' ) sc = (*nit)->z();
         (*nit)->Store( prop_key, sc );
       }

 } // end AssignNodeCoordinatesTo










// MIN-MAX PROPERTIES


/** Returns min/max of property values inside a region into its arguments.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::MinMaxOf( const char* property, double64& gmin, double64& gmax ) const
 {
    ErrorHandler&  csmp_error(ErrorHandler::Instance());

    if ( !pref_.IsDefined(property) ) {
          csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::MinMaxOf", property, "is undefined; nothing could be done" );
          return;
      }
    csmp::Index  gprop_key(pref_.StorageKey(property));
    MinMaxOf( gprop_key, gmin, gmax );

 } // end MinMaxOf




/**
    Recovering and sorting the minimum and maximum Eigen values of the tensor variable.

    @todo SKM this wants to be a lambda function in the next method.
*/
template<size_t dim>
void minMaxEigenValues( const TensorVariable<dim>& ts, double64& tmin, double64& tmax )
 {
    VectorVariable<dim>  evals;
    ts.EigenValues( evals );
    std::set<double64> min_max;
    for ( size_t i=0U; i<dim; i++ ) min_max.insert( evals[i] );
    tmin = (*min_max.begin());
    tmax = (*min_max.rbegin());
 }




/**
    Recovers ranges of variable values from regions, boundaries or split boundaries.
    
    @attention for vector variables the length range is returned.
    @attention for tensor variables the maximum Eigenvalue is recovered.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::MinMaxOf( const csmp::Index& prop_key, double64& vmin, double64& vmax ) const
 {
    // properties / variables placed on the model
    if ( prop_key.place == REGION || prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY ) {
          if ( prop_key.type == SCALAR ) vmin = vmax = this->Read( prop_key );
          else if ( prop_key.type == VECTOR ) {
               VectorVariable<dim>  vc;
               this->Read( prop_key, vc );
               vmin = vmax = vc.Length();
            }
          else if ( prop_key.type == TENSOR ) {
               TensorVariable<dim>  ts;
               this->Read( prop_key, ts );
               minMaxEigenValues( ts, vmin, vmax );
            }
          else if ( prop_key.type == ARRAY ) {
               ArrayVariable  a(prop_key.dataDepth);
               this->Read( prop_key, a );
               a.MinMax( vmin, vmax );
            }
          else if ( prop_key.type == FLAGGEDARRAY ) {
               FlaggedArrayVariable  a(prop_key.dataDepth);
               this->Read( prop_key, a );
               a.MinMax( vmin, vmax );
            }
          return;
      }


   // node properties of any kind
   if ( prop_key.place == NODE ) {
       typename vector<csmp::Node<dim>*>::const_iterator  nit(node_vec_.begin());
       assert( (*nit) != NULL );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*nit)->Read(prop_key); nit++;
            while ( nit!=node_vec_.end() ) {
                 vmin = std::min( vmin, (*nit)->Read(prop_key) );
                 vmax = std::max( vmax, (*nit)->Read(prop_key) );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*nit)->Read(prop_key,vc); nit++;
            vmin = vmax = vc.Length();
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, vc );
                 const double64  vlength(vc.Length());
                 vmin = std::min( vmin, vlength );
                 vmax = std::max( vmax, vlength );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*nit)->Read(prop_key,ts); nit++;
            minMaxEigenValues( ts, vmin, vmax );
            double64 tmin, tmax;
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, ts );
                 minMaxEigenValues( ts, tmin, tmax );
                 vmin = std::min( vmin, tmin );
                 vmax = std::max( vmax, tmax );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*nit)->Read(prop_key,a); nit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 nit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*nit)->Read(prop_key,a); nit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( nit!=node_vec_.end() ) {
                 (*nit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 nit++;
              }
            return;
         }
    } // end node properties


   // element-based properties of any kind
   if ( prop_key.place == ELEMENT || prop_key.place == FACE || prop_key.place == INTER_FACE ) {
       typename vector<SIMPLEX<dim>*>::const_iterator  eit(elmt_vec_.begin());
       assert( (*eit) != NULL );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(prop_key); eit++;
            while ( eit!=elmt_vec_.end() ) {
                 vmin = std::min( vmin, (*eit)->Read(prop_key) );
                 vmax = std::max( vmax, (*eit)->Read(prop_key) );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=elmt_vec_.end() ) {
                 (*eit)->Read( prop_key, vc );
                 const double64  vlength(vc.Length());
                 vmin = std::min( vmin, vlength );
                 vmax = std::max( vmax, vlength );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double64 tmin, tmax;
            while ( eit!=elmt_vec_.end() ) {
                 (*eit)->Read( prop_key, ts );
                 minMaxEigenValues( ts, tmin, tmax );
                 vmin = std::min( vmin, tmin );
                 vmax = std::max( vmax, tmax );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                 (*eit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                 (*eit)->Read( prop_key, a );
                 a.MinMax( amin, amax );
                 vmin = std::min( vmin, amin );
                 vmax = std::max( vmax, amax );
                 eit++;
              }
            return;
         }
    } // end element, face, interface properties


   // element-based integration-point properties of any kind
   if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
       typename vector<SIMPLEX<dim>*>::const_iterator  eit(elmt_vec_.begin());
       assert( (*eit) != NULL );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(0U,prop_key); eit++;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     vmin = std::min( vmin, (*eit)->Read(i,prop_key) );
                     vmax = std::max( vmax, (*eit)->Read(i,prop_key) );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(0U,prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, vc );
                     const double64  vlength(vc.Length());
                     vmin = std::min( vmin, vlength );
                     vmax = std::max( vmax, vlength );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(0U,prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double64 tmin, tmax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, ts );
                     minMaxEigenValues( ts, tmin, tmax );
                     vmin = std::min( vmin, tmin );
                     vmax = std::max( vmax, tmax );
                 }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=1U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
    } // end element integration point properties on element, face or interface


   // FV facet integration-point properties of any kind placed on element, face, or interface
   if ( prop_key.place == FACET_INTEGRATION_POINT ||
        prop_key.place == FACE_FACET_INTEGRATION_POINT ||
        prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT )
     {
       typename vector<SIMPLEX<dim>*>::const_iterator  eit(elmt_vec_.begin());
       assert( (*eit) != NULL );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(0U,0U,prop_key); eit++;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Facets(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     vmin = std::min( vmin, (*eit)->Read(i,j,prop_key) );
                     vmax = std::max( vmax, (*eit)->Read(i,j,prop_key) );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(0U,0U,prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Facets(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, vc );
                     const double64  vlength(vc.Length());
                     vmin = std::min( vmin, vlength );
                     vmax = std::max( vmax, vlength );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(0U,0U,prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double64 tmin, tmax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Facets(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, ts );
                     minMaxEigenValues( ts, tmin, tmax );
                     vmin = std::min( vmin, tmin );
                     vmax = std::max( vmax, tmax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Facets(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Facets(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerFacet(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
    } // end FV facet integration point properties on element, face or interface


   // FV sector integration-point properties of any kind placed on element, face, or interface
   if ( prop_key.place == SECTOR_INTEGRATION_POINT ||
        prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
        prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT  )
     {
       typename vector<SIMPLEX<dim>*>::const_iterator  eit(elmt_vec_.begin());
       assert( (*eit) != NULL );
       if ( prop_key.type == SCALAR ) {
            vmin = vmax = (*eit)->Read(0U,0U,prop_key); eit++;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Sectors(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     vmin = std::min( vmin, (*eit)->Read(i,j,prop_key) );
                     vmax = std::max( vmax, (*eit)->Read(i,j,prop_key) );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == VECTOR ) {
            VectorVariable<dim>  vc;
            (*eit)->Read(0U,0U,prop_key,vc); eit++;
            vmin = vmax = vc.Length();
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Sectors(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, vc );
                     const double64  vlength(vc.Length());
                     vmin = std::min( vmin, vlength );
                     vmax = std::max( vmax, vlength );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == TENSOR ) {
            TensorVariable<dim>  ts;
            (*eit)->Read(0U,0U,prop_key,ts); eit++;
            minMaxEigenValues( ts, vmin, vmax );
            double64 tmin, tmax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Sectors(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, ts );
                     minMaxEigenValues( ts, tmin, tmax );
                     vmin = std::min( vmin, tmin );
                     vmax = std::max( vmax, tmax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == ARRAY ) {
            ArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Sectors(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
       if ( prop_key.type == FLAGGEDARRAY ) {
            FlaggedArrayVariable a(prop_key.dataDepth);
            (*eit)->Read(0U,0U,prop_key,a); eit++;
            a.MinMax( vmin, vmax );
            double64  amin, amax;
            while ( eit!=elmt_vec_.end() ) {
                for ( size_t i=0U; i<(*eit)->Sectors(); i++ )
                  for ( size_t j=0U; j<(*eit)->IntegrationPointsPerSector(); j++ ) {
                     (*eit)->Read( i, j, prop_key, a );
                     a.MinMax( amin, amax );
                     vmin = std::min( vmin, amin );
                     vmax = std::max( vmax, amax );
                  }
                 eit++;
              }
            return;
         }
     } // end sector integration points for any element, face or interface

    vmin = vmax = std::numeric_limits<double64>::quiet_NaN();

    ErrorHandler&  csmp_error(ErrorHandler::Instance());
    csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::MinMaxOf(index)",
                      "placement of property coould not be indentified");

 } // end MinMaxOf(index)













// MANIPULATIONS WITH PROPERTIES




/// @todo (3-D) Refactor!! This is a potential bug if wrong integration properties are used!
template<size_t dim, template<size_t> class SIMPLEX>
template<typename Var>
void ModelSubDomain<dim,SIMPLEX>::InputPropertyValue( const char* input_prop,
                                                      const Var& var,
                                                      SUBDOMAIN_PART sdp )
 {
     csmp::Index prop_key = pref_.StorageKey(input_prop);

     if ( sdp == PERIMETER  and (prop_key.place == REGION  or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY) ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                            input_prop, "placed on Region cannot be assigned just on perimeter");
          return;
       }

     if( prop_key.place == REGION || prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY )
       {
        this->Store( prop_key, var );
        return;
       }

     if ( sdp == COMPLETE ) {
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
                  (*it)->Store( prop_key, var );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
                  for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                    (*it)->Store( i, prop_key, var );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Facets(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerFacet(); j++ )
                  (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Sectors(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerSector(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( typename vector<csmp::Node<dim>*>::iterator
                      nit=node_vec_.begin(); nit!=node_vec_.end(); nit++ )
                  (*nit)->Store( prop_key, var );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                                                "Property placement not recognized");
       }
     else if ( sdp == PERIMETER ) { // property is assigned only to perimeter of boundary
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
                  (*it)->Store( prop_key, var );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
                  for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                    (*it)->Store( i, prop_key, var );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Facets(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerFacet(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Sectors(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerSector(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( typename vector<csmp::Node<dim>*>::iterator
                      nit=PerimeterNodesBegin(); nit!=node_vec_.end(); nit++ )
                  (*nit)->Store( prop_key, var );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                                                "Property placement not recognized");
       }
     else { // property is assigned only to interior of subdomain
           if ( prop_key.place == ELEMENT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
                     (*it)->Store( prop_key, var );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT || prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
                  for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                    (*it)->Store( i, prop_key, var );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT || prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
               for ( size_t i=0U; i<(*it)->Facets(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerFacet(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT || prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
               for ( size_t i=0U; i<(*it)->Sectors(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerSector(); j++ )
                   (*it)->Store( i, j, prop_key, var );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( typename vector<csmp::Node<dim>*>::iterator
                      nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                  (*nit)->Store( prop_key, var );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                                                    "Property placement not supported");
       }

 } // end InputPropertyValue

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const VectorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const VectorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const VectorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const TensorVariable<1>&, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const TensorVariable<2>&, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const TensorVariable<3>&, SUBDOMAIN_PART );



/**
    Helper function that checks a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
*/


/**
    Assigns variable values to model subdomains, but only if the components of the variable do not have
    the flag identified as 'do_not_overwrite'.

    @attention Since tensor variables only have one flag per row, the overwrite protection is applied
    on a row-by-row basis.

    @attention This special treatment is not applied to FV integration points.

*/
template<size_t dim, template<size_t> class SIMPLEX>
template<typename Var>
void ModelSubDomain<dim,SIMPLEX>::InputPropertyValue( const char* input_prop,
                                                      const Var& var,
                                                      VARIABLE_FLAG do_not_overwrite,
                                                      SUBDOMAIN_PART sdp )
 {
     csmp::Index prop_key = pref_.StorageKey(input_prop);

     if ( sdp == PERIMETER  and (prop_key.place == REGION  or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY) ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                            input_prop, "placed on Region cannot be assigned just on perimeter");
          return;
       }

     if( prop_key.place == REGION || prop_key.place == BOUNDARY || prop_key.place == SPLIT_BOUNDARY )
       {
        this->Store( prop_key, var );
        return;
       }

     if ( sdp == COMPLETE ) {
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
                  writeVariableIf( (*it), prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
                  for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                    writeVariableIf( (*it), i, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Facets(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerFacet(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Sectors(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerSector(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( typename vector<csmp::Node<dim>*>::iterator
                      nit=node_vec_.begin(); nit!=node_vec_.end(); nit++ )
                  writeVariableIf( (*nit), prop_key, var, DIRICH );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                                                    "Property placement not recognized");
       }
     else if ( sdp == PERIMETER ) { // property is assigned only to perimeter of boundary
           if ( prop_key.place == ELEMENT or prop_key.place == FACE or prop_key.place == INTER_FACE ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
                  writeVariableIf( (*it), prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
                  for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                    writeVariableIf( (*it), i, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Facets(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerFacet(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=PerimeterElementsBegin(); it!=elmt_vec_.end(); it++ )
               for ( size_t i=0U; i<(*it)->Sectors(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerSector(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( typename vector<csmp::Node<dim>*>::iterator
                      nit=PerimeterNodesBegin(); nit!=node_vec_.end(); nit++ )
                  writeVariableIf( (*nit), prop_key, var, do_not_overwrite );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                                                "Property placement not recognized");
       }
     else { // property is assigned only to interior of subdomain
           if ( prop_key.place == ELEMENT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
                  writeVariableIf( (*it), prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == ELEMENT_INTEGRATION_POINT || prop_key.place == FACE_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_INTEGRATION_POINT ) {
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
                  for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                    writeVariableIf( (*it), i, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == FACET_INTEGRATION_POINT || prop_key.place == FACE_FACET_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
               for ( size_t i=0U; i<(*it)->Facets(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerFacet(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == SECTOR_INTEGRATION_POINT || prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
                     prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
             for ( typename vector<SIMPLEX<dim>*>::iterator
               it=ElementsBegin(); it!=PerimeterElementsBegin(); it++ )
               for ( size_t i=0U; i<(*it)->Sectors(); i++ )
                 for ( size_t j=0U; j<(*it)->IntegrationPointsPerSector(); j++ )
                    writeVariableIf( (*it), i, j, prop_key, var, do_not_overwrite );
             }
           else if ( prop_key.place == NODE ) { // for nodes on first side of interface
                for ( typename vector<csmp::Node<dim>*>::iterator
                      nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                  writeVariableIf( (*nit), prop_key, var, do_not_overwrite );
             }
           else throw csmp::Exception( FATAL_ERROR, "ModelSubDomain<dim,SIMPLEX>::InputPropertyValue",
                                                    "Property placement not supported");
       }

 } // end InputPropertyValue

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ScalarVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const ArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const FlaggedArrayVariable&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const VectorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const VectorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const VectorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const VectorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const VectorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const VectorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const VectorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const VectorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const VectorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Element>::InputPropertyValue( const char*, const TensorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Element>::InputPropertyValue( const char*, const TensorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Element>::InputPropertyValue( const char*, const TensorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,Face>::InputPropertyValue( const char*, const TensorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,Face>::InputPropertyValue( const char*, const TensorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,Face>::InputPropertyValue( const char*, const TensorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );

template void ModelSubDomain<1U,InterFace>::InputPropertyValue( const char*, const TensorVariable<1>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<2U,InterFace>::InputPropertyValue( const char*, const TensorVariable<2>&, VARIABLE_FLAG, SUBDOMAIN_PART );
template void ModelSubDomain<3U,InterFace>::InputPropertyValue( const char*, const TensorVariable<3>&, VARIABLE_FLAG, SUBDOMAIN_PART );



/**
    Returns status of whole or parts of model subdomain if there is a single flag value; else ANY is returned.
    For scalars, a only a single variable flag can be returned, while for vectors and tensors the flag from component 0
    is provided by default, unless any other component number is asked for.
    Note: no matter which component is asked for in a vector or tensor, if any component has flag ANY,
    then the ANY flag is returned by default. (this is checked for in loops from 0 < dim for vectors and 0 <dim*dim for tensors)
*/
template<size_t dim, template<size_t> class SIMPLEX>
VARIABLE_FLAG  ModelSubDomain<dim,SIMPLEX>::PropertyStatus( const char* property, SUBDOMAIN_PART group_flag , size_t component) const
{
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( !pref_.IsDefined(property) ) {
        csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::PropertyStatus", property, "is undefined; nothing could be done" );
        return ANY;
    }
    csmp::Index  prop_key = pref_.StorageKey(property);
    string src("ModelSubDomain<");
    string cache(to_string(dim));
    src += cache;
    src += ",";
    src += typeid(SIMPLEX<dim>).name();
    src +=">::PropertyStatus";

    if ( prop_key.place == FACET_INTEGRATION_POINT ||
         prop_key.place == SECTOR_INTEGRATION_POINT ||
         prop_key.place == FACE_SECTOR_INTEGRATION_POINT ||
         prop_key.place == FACE_FACET_INTEGRATION_POINT ||
         prop_key.place == INTER_FACE_INTEGRATION_POINT ||
         prop_key.place == INTER_FACE_SECTOR_INTEGRATION_POINT ||
         prop_key.place == INTER_FACE_FACET_INTEGRATION_POINT )
    {
        csmp_error.notice( ERROR, src.c_str(), "Method is not implemented for facet and sector integration points");
        return ANY;
    }

    if ( elmt_vec_.empty() ) {
        csmp_error.notice( ERROR, src.c_str(), "ModelSubDomain is empty");
        return ANY;
    }

    if ( prop_key.type == SCALAR || prop_key.type == ARRAY )
    {
        // if all the flags in the region shall be changed
        if ( group_flag == COMPLETE ) {
            if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                return this->Status( prop_key );
            }
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                    if ( status0 != (*eit)->Status(prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                    for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        if ( status0 != (*eit)->Status(i,prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key));
                for ( typename vector<csmp::Node<dim>*>::const_iterator
                      nit=NodesBegin(); nit!=NodesEnd(); nit++ )
                    if ( status0 != (*nit)->Status(prop_key) ) return ANY;
                return status0;
            }
        }
        // if the property shall only be changed on the Subdomain boundary
        else if ( group_flag == PERIMETER ) {
            if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                return this->Status( prop_key );
            }
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(InteriorElements())->Status(prop_key));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                    if ( status0 != (*eit)->Status(prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(InteriorElements())->Status(0,prop_key));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                    for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        if ( status0 != (*eit)->Status(i,prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(InteriorNodes())->Status(prop_key));
                for ( typename vector<csmp::Node<dim>*>::const_iterator
                      nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                    if ( status0 != (*nit)->Status(prop_key) ) return ANY;
                return status0;
            }
        }
        else if ( group_flag == INTERIOR ) {
            if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                return this->Status( prop_key );
            }
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                    if ( status0 != (*eit)->Status(prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                    for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        if ( status0 != (*eit)->Status(i,prop_key) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key));
                for ( typename vector<csmp::Node<dim>*>::const_iterator
                      nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                    if ( status0 != (*nit)->Status(prop_key) ) return ANY;
                return status0;
            }
        }

    } // END SCALARS, ARRAYS


    if ( prop_key.type == VECTOR || prop_key.type == TENSOR || prop_key.type == FLAGGEDARRAY )
    {
        size_t length( prop_key.flagDepth );

        // if all the flags in the region shall be changed
        if ( group_flag == COMPLETE ) {
            if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                VARIABLE_FLAG  status0(this->Status(prop_key,component));
                for ( size_t n=1U; n<length; n++ )
                    if ( status0 != this->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key,component));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                    for ( size_t n=0U; n<length; n++ )
                        if ( status0 != (*eit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key,component));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                    for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        for ( size_t n=0U; n<length; n++ )
                            if ( status0 != (*eit)->Status(i,prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key,component));
                for ( typename vector<csmp::Node<dim>*>::const_iterator
                      nit=NodesBegin(); nit!=NodesEnd(); nit++ )
                    for ( size_t n=0U; n<length; n++ )
                        if ( status0 != (*nit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
        }
        // if the property shall only be changed on the Subdomain boundary
        else if ( group_flag == PERIMETER ) {
            if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                VARIABLE_FLAG  status0(this->Status(prop_key,component));
                for ( size_t n=1U; n<length; n++ )
                    if ( status0 != this->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(InteriorElements())->Status(prop_key,component));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                    for ( size_t n=0U; n<length; n++ )
                        if ( status0 != (*eit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(InteriorElements())->Status(0,prop_key,component));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                    for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        for ( size_t n=0U; n<length; n++ )
                            if ( status0 != (*eit)->Status(i,prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(InteriorNodes())->Status(prop_key,component));
                for ( typename vector<csmp::Node<dim>*>::const_iterator
                      nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                    for ( size_t n=0U; n<length; n++ )
                        if ( status0 != (*nit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
        }
        else if ( group_flag == INTERIOR ) {
            if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                VARIABLE_FLAG  status0(this->Status(prop_key,component));
                for ( size_t n=1U; n<length; n++ )
                    if ( status0 != this->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT or
                 prop_key.place == FACE or
                 prop_key.place == INTER_FACE ) {
                VARIABLE_FLAG  status0(E(0)->Status(prop_key,component));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                    for ( size_t n=0U; n<length; n++ )
                        if ( status0 != (*eit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                VARIABLE_FLAG  status0(E(0)->Status(0,prop_key,component));
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                    for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        for ( size_t n=0U; n<length; n++ )
                            if ( status0 != (*eit)->Status(i,prop_key,n) ) return ANY;
                return status0;
            }
            if ( prop_key.place == NODE ) {
                VARIABLE_FLAG  status0(N(0)->Status(prop_key,component));
                for ( typename vector<csmp::Node<dim>*>::const_iterator
                      nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                    for ( size_t n=0U; n<length; n++ )
                        if ( status0 != (*nit)->Status(prop_key,n) ) return ANY;
                return status0;
            }
        }

    } // END VECTORS, TENSORS, FLAGGEDARRAYS

    return ANY;

} // end PropertyStatus




template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::ChangePropertyStatus( const char* property,
                                                        VARIABLE_FLAG new_status_of_scalar,
                                                        SUBDOMAIN_PART sdpart )
 {
    vector<VARIABLE_FLAG>  status(1U,new_status_of_scalar);
    ChangePropertyStatus( property, status, sdpart );
 }



template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::ChangePropertyStatusWhere( const char* property,
                                                             VARIABLE_FLAG new_status_of_scalar,
                                                             double64 min_value_to_change,
                                                             double64 max_value_to_change )
 {
    vector<VARIABLE_FLAG>  status(1U,new_status_of_scalar);
    ChangePropertyStatusWhere( property, status, min_value_to_change, max_value_to_change );
 }



/// changes the variable flag to status for those group members which carry the group_flag.
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::ChangePropertyStatus( const char* property,
                                                        const std::vector<VARIABLE_FLAG>& status,
                                                        SUBDOMAIN_PART group_flag )
 {
    /// @todo (2-P) Rm rtti
    csmp::Index  prop_key = pref_.StorageKey(property);
    string src("ModelSubDomain<");
    src += to_string(dim);
    src += ",";
    src += "Element/Face/Interface";
    src +=">::ChangePropertyStatus:";

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, src.c_str(), "Method not implemented for tensor properties yet");

    if ( status.empty() )
      throw csmp::Exception( ERROR, src.c_str(), "Status vector has not been initialized");

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( elmt_vec_.empty() ) {
         csmp_error.notice( ERROR, src.c_str(), "Region is empty.");
         return;
      }

    if ( prop_key.type == SCALAR || prop_key.type == ARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                    this->Status( prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); ++eit )
                      (*eit)->Status( prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); ++eit )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=NodesEnd(); ++nit )
                      (*nit)->Status( prop_key, status[0U] );
                    return;
                 }
            }
          // if the flags shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); ++eit )
                      (*eit)->Status( prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); ++eit )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=PerimeterNodesBegin(); nit!=NodesEnd(); ++nit )
                      (*nit)->Status( prop_key, status[0U] );
                    return;
                 }
            }
          // if the flags for all entities in the interior of the region shall be changed
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); ++eit )
                      (*eit)->Status( prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); ++eit )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status[0U] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=PerimeterNodesBegin(); ++nit )
                      (*nit)->Status( prop_key, status[0U] );
                    return;
                 }
            }

      } // END SCALARS, ARRAYS

    if ( prop_key.type == VECTOR || prop_key.type == TENSOR || prop_key.type == FLAGGEDARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                    for ( size_t n=0U; n<status.size(); n++ )
                        this->Status( prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t n=0U; n<status.size(); n++ )
                        (*eit)->Status( prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        for ( size_t n=0U; n<status.size(); n++ )
                          (*eit)->Status( i, prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=NodesEnd(); nit++ )
                      for ( size_t n=0U; n<status.size(); n++ )
                        (*nit)->Status( prop_key, n, status[n] );
                    return;
                 }
            }
          // if the property shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t n=0U; n<status.size(); n++ )
                        (*eit)->Status( prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        for ( size_t n=0U; n<status.size(); n++ )
                          (*eit)->Status( i, prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                      for ( size_t n=0U; n<status.size(); n++ )
                        (*nit)->Status( prop_key, n, status[n] );
                    return;
                 }
            }
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                      for ( size_t n=0U; n<status.size(); n++ )
                        (*eit)->Status( prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        for ( size_t n=0U; n<status.size(); n++ )
                          (*eit)->Status( i, prop_key, n, status[n] );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                      for ( size_t n=0U; n<status.size(); n++ )
                        (*nit)->Status( prop_key, n, status[n] );
                    return;
                 }
            }

      } // END VECTORS,TENSORS,FLAGGEDARRAYS

    cout <<"\n'"<< property <<"' ";
    throw csmp::Exception( ERROR, src.c_str(), "Property placement not recognized");

 } // end ChangePropertyStatus






/// changes the variable flag to status for those group members which carry the group_flag.
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::ChangePropertyStatus( const char* property,
                                                        size_t position,
                                                        VARIABLE_FLAG status,
                                                        SUBDOMAIN_PART group_flag )
 {
    /// @todo (2-P) Rm rtti
    csmp::Index  prop_key = pref_.StorageKey(property);
    string src("ModelSubDomain<");
    string cache(to_string(dim));
    src += cache;
    src += ",";
    src += typeid(SIMPLEX<dim>).name();
    src +=">::ChangePropertyStatus";

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, src.c_str(), "Method not implemented for tensor properties yet");

    ErrorHandler&  csmp_error( ErrorHandler::Instance() );

    if ( elmt_vec_.empty() ) {
         csmp_error.notice( ERROR, src.c_str(), "Region is empty");
         return;
      }

    if ( prop_key.type == SCALAR || prop_key.type == ARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                    this->Status( prop_key, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                      (*eit)->Status( prop_key, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=NodesEnd(); nit++ )
                      (*nit)->Status( prop_key, status );
                    return;
                 }
            }
          // if the property shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                      (*eit)->Status( prop_key, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                      (*nit)->Status( prop_key, status );
                    return;
                 }
            }
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                      (*eit)->Status( prop_key, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                        (*eit)->Status( i, prop_key, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                      (*nit)->Status( prop_key, status );
                    return;
                 }
            }

      } // END SCALARS, ARRAYS

    if ( prop_key.type == VECTOR || prop_key.type == TENSOR || prop_key.type == FLAGGEDARRAY )
      {
          // if all the flags in the region shall be changed
          if ( group_flag == COMPLETE ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY ) {
                    this->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                        (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                          (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=NodesEnd(); nit++ )
                        (*nit)->Status( prop_key, position, status );
                    return;
                 }
            }
          // if the property shall only be changed on the Subdomain boundary
          else if ( group_flag == PERIMETER ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                        (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=PerimeterElementsBegin(); eit!=ElementsEnd(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                          (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=PerimeterNodesBegin(); nit!=NodesEnd(); nit++ )
                        (*nit)->Status( prop_key, position, status );
                    return;
                 }
            }
          else if ( group_flag == INTERIOR ) {
               if ( prop_key.place == REGION or prop_key.place == BOUNDARY or prop_key.place == SPLIT_BOUNDARY )
                throw csmp::Exception( ERROR, src.c_str(),
                   "Flags of subdomain properties can only be assigned to complete subdomains");
               if ( prop_key.place == ELEMENT or
                    prop_key.place == FACE or
                    prop_key.place == INTER_FACE ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                        (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
                    for ( typename vector<SIMPLEX<dim>*>::iterator
                          eit=ElementsBegin(); eit!=PerimeterElementsBegin(); eit++ )
                      for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                          (*eit)->Status( prop_key, position, status );
                    return;
                 }
               if ( prop_key.place == NODE ) {
                    for ( typename vector<csmp::Node<dim>*>::iterator
                          nit=NodesBegin(); nit!=PerimeterNodesBegin(); nit++ )
                        (*nit)->Status( prop_key, position, status );
                    return;
                 }
            }

      } // END VECTORS,TENSORS, FLAGGEDARRAYS

    cout <<"\n'"<< property <<"' ";
    throw csmp::Exception( ERROR, src.c_str(), "Property placement not recognized");

 } // end ChangePropertyStatus



template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::ChangePropertyStatusWhere( const char* property,
                                                             const std::vector<VARIABLE_FLAG>& status,
                                                             double64 pmin, double64 pmax )
 {
    csmp::Index  prop_key = pref_.StorageKey(property);

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere",
                            "Method not implemented for tensor properties yet");
    if ( status.empty() )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere",
                            "Status vector has not been initialized");
    if ( elmt_vec_.empty() )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere", "Region is empty");

    if ( prop_key.type == SCALAR )
      {
        ScalarVariable  sc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     (*eit)->Status( prop_key, status[0U] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, sc );
                     if ( sc.IsWithinRange( pmin, pmax ) )
                       (*eit)->Status( i, prop_key, status[0U] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     (*nit)->Status( prop_key, status[0U] );
                }
           }
       } // END SCALARS

    if ( prop_key.type == ARRAY )
      {
         ArrayVariable  av;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     (*eit)->Status( prop_key, status[0U] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, av );
                     if ( av.IsWithinRange( pmin, pmax ) )
                       (*eit)->Status( i, prop_key, status[0U] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     (*nit)->Status( prop_key, status[0U] );
                }
           }
       } // END ARRAYS

    if ( prop_key.type == VECTOR )
      {
         VectorVariable<dim>  vc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                     for ( size_t n=0U; n<status.size(); n++ )
                       (*eit)->Status( prop_key, n, status[n] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, vc );
                     if ( vc.IsWithinRange( pmin, pmax ) )
                       for ( size_t n=0U; n<status.size(); n++ )
                         (*eit)->Status( i, prop_key, n, status[n] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                     for ( size_t n=0U; n<status.size(); n++ )
                       (*nit)->Status( prop_key, n, status[n] );
                }
           }
       } // END VECTORS

    if ( prop_key.type == TENSOR )
      {
         TensorVariable<dim>  ts;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                     for ( size_t n=0U; n<status.size(); n++ )
                       (*eit)->Status( prop_key, n, status[n] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, ts );
                     if ( ts.IsWithinRange( pmin, pmax ) )
                       for ( size_t n=0U; n<status.size(); n++ )
                         (*eit)->Status( i, prop_key, n, status[n] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                     for ( size_t n=0U; n<status.size(); n++ )
                       (*nit)->Status( prop_key, n, status[n] );
                }
           }
       } // END TENSORS

    if ( prop_key.type == FLAGGEDARRAY )
      {
         FlaggedArrayVariable fv;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                     for ( size_t n=0U; n<status.size(); n++ )
                       (*eit)->Status( prop_key, n, status[n] );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, fv );
                     if ( fv.IsWithinRange( pmin, pmax ) )
                       for ( size_t n=0U; n<status.size(); n++ )
                         (*eit)->Status( i, prop_key, n, status[n] );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                     for ( size_t n=0U; n<status.size(); n++ )
                       (*nit)->Status( prop_key, n, status[n] );
                }
           }
       } // END FLAGGEDARRAYS

 } // end ChangePropertyStatusWhere (if in value range)









template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::ChangePropertyStatusWhere( const char* property,
                                                             size_t position,
                                                             VARIABLE_FLAG status,
                                                             double64 pmin, double64 pmax )
 {
    csmp::Index  prop_key = pref_.StorageKey(property);

    if ( prop_key.type == TENSOR )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere",
                            "Method not implemented for tensor properties yet");

    if ( elmt_vec_.empty() )
      throw csmp::Exception( ERROR, "ModelSubDomain<dim>::ChangePropertyStatusWhere", "Region is empty");

    if ( prop_key.type == SCALAR )
      {
        ScalarVariable  sc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     (*eit)->Status( prop_key, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, sc );
                     if ( sc.IsWithinRange( pmin, pmax ) )
                       (*eit)->Status( i, prop_key, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, sc );
                   if ( sc.IsWithinRange( pmin, pmax ) )
                     (*nit)->Status( prop_key, status );
                }
           }
       } // END SCALARS

    if ( prop_key.type == ARRAY )
      {
         ArrayVariable  av;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     (*eit)->Status( prop_key, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, av );
                     if ( av.IsWithinRange( pmin, pmax ) )
                       (*eit)->Status( i, prop_key, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, av );
                   if ( av.IsWithinRange( pmin, pmax ) )
                     (*nit)->Status( prop_key, status );
                }
           }
       } // END ARRAYS

    if ( prop_key.type == VECTOR )
      {
         VectorVariable<dim>  vc;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                       (*eit)->Status( prop_key, position, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, vc );
                     if ( vc.IsWithinRange( pmin, pmax ) )
                         (*eit)->Status( prop_key, position, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, vc );
                   if ( vc.IsWithinRange( pmin, pmax ) )
                       (*nit)->Status( prop_key, position, status );
                }
           }
       } // END VECTORS

    if ( prop_key.type == TENSOR )
      {
         TensorVariable<dim>  ts;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                       (*eit)->Status( prop_key, position, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, ts );
                     if ( ts.IsWithinRange( pmin, pmax ) )
                         (*eit)->Status( prop_key, position, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, ts );
                   if ( ts.IsWithinRange( pmin, pmax ) )
                       (*nit)->Status( prop_key, position, status );
                }
           }
       } // END TENSORS

    if ( prop_key.type == FLAGGEDARRAY )
      {
         FlaggedArrayVariable fv;

         if ( prop_key.place == ELEMENT or
              prop_key.place == FACE or
              prop_key.place == INTER_FACE ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ ) {
                   (*eit)->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                       (*eit)->Status( prop_key, position, status );
                }
           }
         else if ( prop_key.place == ELEMENT_INTEGRATION_POINT ) {
              for ( typename vector<SIMPLEX<dim>*>::iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                     (*eit)->Read( i, prop_key, fv );
                     if ( fv.IsWithinRange( pmin, pmax ) )
                         (*eit)->Status( prop_key, position, status );
                  }
           }
         else if ( prop_key.place == NODE ) {
              for ( typename vector<csmp::Node<dim>*>::iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ ) {
                   (*nit)->Read( prop_key, fv );
                   if ( fv.IsWithinRange( pmin, pmax ) )
                       (*nit)->Status( prop_key, position, status );
                }
           }
       } // END VECTORS

 } // end ChangePropertyStatusWhere (if in value range)




















// INTERPOLATIONS AND EXTRAPOLATIONS


/**

InterpolateNodePropertyToElementProperty() interpolates node property to the
'barycenter' of the triangle. The resulting
value is different from the result obtained by applying the interrelation
subclass NodeToElementProperty. The Calculate() method of the latter assigns
the average value of the 3 element nodes to the element property 'eprop'.

@section arguments Input Arguments

The two string arguments define the name of the node property which is
interpolated and the name of the element property to which the
interpolated value is written.

@section application Application

This method was developed for two-dimensional convection calculations
where the fluid density as element property must be exactly the same for in
the two adjacent triangles which make up a square in regular-gridded meshes.

@section messages Messages

If the variable placement or type of the specified properties fails to
match the specifications outlined above,
InterpolateNodePropertyToElementProperty() will report an error and
return without completing its task.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToElementProperty( const char* nprop, const char* eprop )
 {
     csmp::Index  e_key = pref_.StorageKey(eprop),
                  n_key = pref_.StorageKey(nprop);

     // 1. check whether conditions for operation are O.K.
     if ( e_key.place != ELEMENT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToElementProperty",
                                 "Property arg2 is not an element property, nothing was done...");
          return;
       }
     if ( n_key.place != NODE ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToElementProperty",
                                 "Property arg1 is not a node property, nothing was done...");
          return;
       }
     if ( e_key.type != n_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToElementProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

     switch( e_key.type )
       {
           case SCALAR: {
                ScalarVariable  sce;
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  {
                     (*eit)->PropertyValueAtBaryCenter( n_key, sce );
                     (*eit)->Store( e_key, sce );
                  }
                }
             break;
           case VECTOR: {
                VectorVariable<dim>  vce;
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  {
                     (*eit)->PropertyValueAtBaryCenter( n_key, vce );
                     (*eit)->Store( e_key, vce );
                  }
               }
            break;
          case TENSOR: {
               TensorVariable<dim>  tse;
               for ( typename vector<SIMPLEX<dim>*>::iterator
                     eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                 {
                     (*eit)->PropertyValueAtBaryCenter( n_key, tse );
                     (*eit)->Store( e_key, tse );
                 }
              }
            break;
          default:
            throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToElementProperty",
                                   "Property type could not be identified, nothing was done...");

     } // end switch

 } // end InterpolateNodeToElementProperty





template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* cprop )
 {
     csmp::Index  c_key = pref_.StorageKey(cprop),
                  n_key = pref_.StorageKey(nprop);

     // 1. check whether conditions for operation are O.K.
     if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToIntegrationPointProperty",
                                 "Property arg2 is not an constraint point property, nothing was done...");
          return;
       }
     if ( n_key.place != NODE ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToIntegrationPointProperty",
                                 "Property arg1 is not a node property, nothing was done...");
          return;
       }
     if ( c_key.type != n_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToIntegrationPointProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

     vector<double64>  IPOL;

     switch( c_key.type )
       {
           case SCALAR: {
                ScalarVariable  sce;
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                    {
                       assert( (*eit)->IntegrationPoints() == (*eit)->FE()->IntegrationPoints() );
                       (*(*eit)).N_AtIntegrationPoint( i, IPOL );
                       sce=0.;
                       for ( size_t j=0U; j<(*eit)->Nodes(); j++ ) sce += IPOL[j] * (*eit)->N(j)->Read( n_key );
                       (*eit)->Store( i, c_key, sce );
                    }
                }
             break;
           case VECTOR: {
                VectorVariable<dim>  vce, vce2;
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                    {
                       assert( (*eit)->IntegrationPoints() == (*eit)->FE()->IntegrationPoints() );
                       (*(*eit)).N_AtIntegrationPoint( i, IPOL );
                       vce=0.;
                       for ( size_t j=0U; j<(*eit)->Nodes(); j++ )
                         {
                            (*eit)->N(j)->Read( n_key, vce2 );
                            vce += (vce2 * IPOL[j]);
                         }
                       (*eit)->Store( i, c_key, vce );
                    }
                }
            break;
          case TENSOR: {
               TensorVariable<dim>  tse, tse2;
                for ( typename vector<SIMPLEX<dim>*>::iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                    {
                       assert( (*eit)->IntegrationPoints() == (*eit)->FE()->IntegrationPoints() );
                       (*(*eit)).N_AtIntegrationPoint( i, IPOL );
                       tse=0.;
                       for ( size_t j=0U; j<(*eit)->Nodes(); j++ )
                         {
                            (*eit)->N(j)->Read( n_key, tse2 );
                            tse += (tse2 * IPOL[j]);
                         }
                       (*eit)->Store( i, c_key, tse );
                    }
                }
            break;
          default:
            throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateNodeToIntegrationPointProperty",
                                   "Property type could not be identified, nothing was done...");

     } // end switch

 } // end InterpolateNodeToIntegrationPointProperty





/**

This method interpolates the specified constraint point property
to the target element property.

@param cprop The name of the constraint point property which shall be interpolated to the
@param eprop element property (arg2).

The results are returned to the Model.

@section implementation Implementation

Currently the method averages the IntegrationPoint variable values
to find the element property value.

@section application Application

To visualise constraint point properties as CELL_CENTERED variables in
VTK they have to have a unique value in the element = CELL.

Alternatively, one could integrate the property over the element and
divide the result value by the element area. However, this would work
only for scalar properties.

@section messages Messages

Consistency checks are performed on the placement and type of the
input variables.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::InterpolateIntegrationPointToElementProperty( const char* cprop, const char* eprop )
 {
     csmp::Index  e_key = pref_.StorageKey(eprop);
     csmp::Index  c_key = pref_.StorageKey(cprop);

     if ( e_key.place != ELEMENT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateIntegrationPointToElementProperty",
                                 "Property arg2 is not an element property, nothing was done...");
          return;
       }
     if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateIntegrationPointToElementProperty",
                                 "Property arg1 is not a constraint point property, nothing was done...");
          return;
       }
     if ( e_key.type != c_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::InterpolateIntegrationPointToElementProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

    if ( e_key.type == SCALAR ) {
        ScalarVariable           sc;
          vector<ScalarVariable >  sc_vec;
          for ( typename vector<SIMPLEX<dim>*>::iterator
                it=ElementsBegin(); it!=ElementsEnd(); it++ ) {
              (*it)->IntegrationPointPropertyVector( c_key, sc_vec );
              average( sc_vec, sc );
              (*it)->Store( e_key, sc );
            }
      }
    else if ( e_key.type == VECTOR ) {
        VectorVariable<dim>           vc;
        vector<VectorVariable<dim> >  vc_vec;
          for ( typename vector<SIMPLEX<dim>*>::iterator
                it=ElementsBegin(); it!=ElementsEnd(); it++ ) {
              (*it)->IntegrationPointPropertyVector( c_key, vc_vec );
              average( vc_vec, vc );
              (*it)->Store( e_key, vc );
            }
      }
    if ( e_key.type == TENSOR ) {
        TensorVariable<dim>           ts;
        vector<TensorVariable<dim> >  ts_vec;
          for ( typename vector<SIMPLEX<dim>*>::iterator
                it=ElementsBegin(); it!=ElementsEnd(); it++ ) {
              (*it)->IntegrationPointPropertyVector( c_key, ts_vec );
              average( ts_vec, ts );
              (*it)->Store( e_key, ts );
            }
      }

 } // end InterpolateIntegrationPointToElementProperty

template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToIntegrationPointProperty( const char* eprop,
                                                                                const char* cprop )
{
    assert( !elmt_vec_.empty() );

    csmp::Index  e_key = pref_.StorageKey(eprop);
    csmp::Index  c_key = pref_.StorageKey(cprop);

    // 1. check whether conditions for operation are O.K.
    if ( e_key.place != ELEMENT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToIntegrationPointProperty",
                                "Property arg1 is not an element property, nothing was done...");
         return;
      }
    if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToIntegrationPointProperty",
                                "Property arg2 is not an integration point property, nothing was done...");
         return;
      }
    if ( e_key.type != c_key.type ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToIntegrationPointProperty",
                                "Properties are not of the same type, nothing was done...");
         return;
      }

    // 2. For each element the property value is assigned for each integration point inside that element

    // SCALAR VARIABLES
    // ----------------
    if ( e_key.type == SCALAR )
    {
         ScalarVariable sc;
         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                        eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
             {
                 // reading element property
                 sc = (*eit)->Read( e_key );

                 for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                     (*eit)->Store( i, c_key, sc );
              }
    }

    // VECTOR VARIABLES
    // ----------------
    if ( e_key.type == VECTOR )
    {
         VectorVariable<dim>  vc;

         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                        eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
             {
                 // reading element property
                 (*eit)->Read( e_key, vc );

                 for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                     (*eit)->Store( i, c_key, vc );
             }
    }

    // TENSOR VARIABLES
    // ----------------
    if ( e_key.type == TENSOR )
    {
         TensorVariable<dim>  ts;

         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                        eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
             {
                 (*eit)->Read( e_key, ts );

                 for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                     (*eit)->Store( i, c_key, ts );
             }
    }

} // end ExtrapolateElementToIntegrationPointProperty


template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToFacetIntegrationPointProperty( const char* eprop,
                                                                                      const char* fipprop )
{
    assert( !elmt_vec_.empty() );

    csmp::Index  e_key = pref_.StorageKey(eprop);
    csmp::Index  fip_key = pref_.StorageKey(fipprop);

    // 1. check whether conditions for operation are O.K.
    if ( e_key.place != ELEMENT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToFacetIntegrationPointProperty",
                                "Property arg1 is not an element property, nothing was done...");
         return;
      }
    if ( fip_key.place != FACET_INTEGRATION_POINT ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToFacetIntegrationPointProperty",
                                "Property arg2 is not an integration point property, nothing was done...");
         return;
      }
    if ( e_key.type != fip_key.type ) {
         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToFacetIntegrationPointProperty",
                                "Properties are not of the same type, nothing was done...");
         return;
      }

    // 2. For each element the property value is assigned for each facet integration point inside that element

    // SCALAR VARIABLES
    // ----------------
    if ( e_key.type == SCALAR )
    {
         ScalarVariable sc;
         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                        eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
             {
                 // reading element property
                 sc = (*eit)->Read( e_key );

                 for ( size_t i=0U; i<(*eit)->FV()->Facets(); i++ )
                    for ( size_t j=0U; j<(*eit)->FV()->IntegrationPointsPerFacet(); j++ )
                     (*eit)->Store( i, j, fip_key, sc );
              }
    }

    // VECTOR VARIABLES
    // ----------------
    if ( e_key.type == VECTOR )
    {
         VectorVariable<dim>  vc;

         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                        eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
             {
                 // reading element property
                 (*eit)->Read( e_key, vc );

                 for ( size_t i=0U; i<(*eit)->FV()->Facets(); i++ )
                    for ( size_t j=0U; j<(*eit)->IntegrationPointsPerFacet(); j++ )
                     (*eit)->Store( i, j, fip_key, vc );
             }
    }

    // TENSOR VARIABLES
    // ----------------
    if ( e_key.type == TENSOR )
    {
         TensorVariable<dim>  ts;

         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                        eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
             {
                 (*eit)->Read( e_key, ts );

                 for ( size_t i=0U; i<(*eit)->FV()->Facets(); i++ )
                    for ( size_t j=0U; j<(*eit)->FV()->IntegrationPointsPerFacet(); j++ )
                     (*eit)->Store( i, j, fip_key, ts );
             }
    }

} // end ExtrapolateElementToFacetIntegrationPointProperty




/**
    Interpolates piecewise constant element properties to the nodes.
    Since each node belongs to multiple elements their contributions to this node have to be weighted.
    This method offers two ways to do this.
    
    'by-distance' is the default. In this approach, contributions are weighted by the inverse of the 
    distance of each element's barycentre from the node. Thus, small elements have a bigger role than
    big ones.
    
    'default' the contributions are weighted by element length/area/volume. This weighs bigger elements
    more than small ones. In many cases this is not what one wants, but it is still offered for
    some specific applications.
*/
template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToNodeProperty( const char* eprop,
                                                                     const char* nprop,
                                                                     bool  by_distance )
 {
     assert( !node_vec_.empty() );
     assert( !elmt_vec_.empty() );

     csmp::Index  e_key = pref_.StorageKey(eprop),
                  n_key = pref_.StorageKey(nprop);

     // 1. check whether conditions for operation are O.K.
     if ( e_key.place != ELEMENT ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToNodeProperty",
                                 "Property arg1 is not an element property, nothing was done...");
          return;
       }
     if ( n_key.place != NODE ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToNodeProperty",
                                 "Property arg2 is not a node property, nothing was done...");
          return;
       }
     if ( e_key.type != n_key.type ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateElementToNodeProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

     // 2. The element property values are weighted by the area of the element and then
     //    added for each node to a property array vector, the number of additions to each
     //    node is counted
     vector<double64>  sum( Nodes(), 0. );

     RenumberNodes();

     // SCALAR VARIABLES
     // ----------------
     if ( e_key.type == SCALAR ) {
          vector<double64>  sc_data( Nodes(), 0. );
          if ( !by_distance )
            {
               for ( typename vector<SIMPLEX<dim>*>::const_iterator
                     eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                 {
                    // getting the properties
                    double64 sc = (*eit)->Read( e_key );

                    // weighting property value by element volume
                    double64 volume = (*eit)->Volume();
                    sc *= volume;

                    for ( size_t i=0U; i<(*eit)->Nodes(); i++ ) {
                         sc_data[ (*eit)->N(i)->Idx() ] += sc;
                         sum[ (*eit)->N(i)->Idx() ]     += volume;
                      }
                 }
            }
          else
            {
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  {
                     // reading element property
                     double64 sc = (*eit)->Read( e_key );
                     // measuring distances from element to nodes to
                     // weight property values
                     Point<dim>  ctr((*eit)->BaryCenter());

                     for ( size_t i=0U; i<(*eit)->Nodes(); i++ ) {
                          double64 dist((ctr - (*eit)->N(i)->Coordinate()).Length());
                          // memorizing distances and weighted values
                          sc_data[ (*eit)->N(i)->Idx() ] += sc / dist;
                          sum[ (*eit)->N(i)->Idx() ]     += 1. / dist;
                       }
                  }
            }

          // calculating the nodal averages and mapping them back to the nodes
          for ( typename vector<csmp::Node<dim>*>::iterator
                nit=NodesBegin(); nit!=NodesEnd(); nit++ )
            {
               ScalarVariable res( (*nit)->Status( n_key ), sc_data[ (*nit)->Idx() ] / sum[ (*nit)->Idx() ] );
               (*nit)->Store( n_key, res );
            }
       }



     // VECTOR VARIABLES
     // ----------------
     if ( e_key.type == VECTOR ) {
          VectorVariable<dim>  vc; vc = 0.;
          vector<VectorVariable<dim> >  vc_data(Nodes(), vc );
          if ( !by_distance )
            {
               for ( typename vector<SIMPLEX<dim>*>::const_iterator
                     eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                 {
                    (*eit)->Read( e_key, vc );
                    double64  volume = (*eit)->Volume();
                    vc *= volume;

                    for ( size_t i=0; i<(*eit)->Nodes(); i++ )
                      {
                         vc_data[ (*eit)->N(i)->Idx() ] += vc;
                         sum[ (*eit)->N(i)->Idx() ] += volume;
                      }
                 }
            }
          else
            {
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  {
                     // reading element property
                     (*eit)->Read( e_key, vc );
                     // measuring distances from element to nodes to
                     // weight property values
                     Point<dim>  ctr((*eit)->BaryCenter());

                     for ( size_t i=0U; i<(*eit)->Nodes(); i++ ) {
                          double64 dist((ctr - (*eit)->N(i)->Coordinate()).Length());
                          // memorizing distances and weighted values
                          vc_data[ (*eit)->N(i)->Idx() ] += vc / dist;
                          sum[ (*eit)->N(i)->Idx() ]     += 1. / dist;
                       }
                  }
            }
          for ( typename vector<csmp::Node<dim>*>::iterator
                nit=NodesBegin(); nit!=NodesEnd(); nit++ )
            {
               vc = vc_data[ (*nit)->Idx() ] / sum[ (*nit)->Idx() ];
               (*nit)->Store( n_key, vc );
            }
       }



     // TENSOR VARIABLES
     // ----------------
     if ( e_key.type == TENSOR ) {
          TensorVariable<dim>  ts; ts = 0.;
          vector<TensorVariable<dim> >  ts_data(Nodes(), ts );
          if ( !by_distance )
            {
               for ( typename vector<SIMPLEX<dim>*>::const_iterator
                     eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                 {
                    (*eit)->Read( e_key, ts );
                    double64 volume = (*eit)->Volume();
                    ts *= volume;

                    for ( size_t i=0U; i<(*eit)->Nodes(); i++ )
                      {
                         ts_data[ (*eit)->N(i)->Idx() ] += ts;
                         sum[ (*eit)->N(i)->Idx() ] += volume;
                      }
                 }
            }
          else
            {
                for ( typename vector<SIMPLEX<dim>*>::const_iterator
                      eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                  {
                     // reading element property
                     (*eit)->Read( e_key, ts );
                     // measuring distances from element to nodes to
                     // weight property values
                     Point<dim>  ctr((*eit)->BaryCenter());

                     for ( size_t i=0U; i<(*eit)->Nodes(); i++ ) {
                          double64 dist((ctr - (*eit)->N(i)->Coordinate()).Length());
                          // memorizing distances and weighted values
                          ts_data[ (*eit)->N(i)->Idx() ] += ts / dist;
                          sum[ (*eit)->N(i)->Idx() ]     += 1. / dist;
                       }
                  }
            }
          for ( typename vector<csmp::Node<dim>*>::iterator
                nit=NodesBegin(); nit!=NodesEnd(); nit++ )
            {
               ts = ts_data[ (*nit)->Idx() ] / sum[ (*nit)->Idx() ];
               (*nit)->Store( n_key, ts );
            }
       }
    if (this->Verbose()){
        cout <<"\nModel<"<<dim<<">::ExtrapolateElementToNodeProperty: ";
        cout <<"'" << eprop <<"' has been successfully extrapolated to '"<< nprop <<"'." << endl;
    }

 } // end ExtrapolateElementToNodeProperty



 

/**

Using the element interpolation functions, this method extrapolates the
desired integration point property to the nodes. Constributions of adjacent
elements are averaged but no weighting by element size or proximity of
barycentre to the node is applied.

@section arguments Input Arguments

The names of the targeted integration point and node variables.

The result of the extrapolation is returned into the Model property
storage.

@section implementation Implementation

The method depends on a corresponding function of the FiniteElement
which is for all numerically integrated element types.
For elements with quadratic inpterpolation functions, the method
uses a linear extrapolation which is justified for solution variable derivatives
as they define a trilinear field on the IntegrationPoints.

@section application Application

Use this for numerically integrated finite elements, in particular ones with 
an interpolation order greater than 1.

@section messages Messages

A consistency check on variable type and placement is performed.
 */
template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop )
 {
     csmp::Index  n_key = pref_.StorageKey(nprop);
     csmp::Index c_key = pref_.StorageKey(cprop);

     ErrorHandler& csmp_error( ErrorHandler::Instance() );

     if ( n_key.place != NODE ) {
          csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateIntegrationPointToNodeProperty",
                                 "Property arg2 is not a node property, nothing was done...");
          return;
       }
     if ( c_key.place != ELEMENT_INTEGRATION_POINT ) {
          csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateIntegrationPointToNodeProperty",
                                 "Property arg1 is not a constraint point property, nothing was done...");
          return;
       }
     if ( n_key.type != c_key.type ) {
          csmp_error.notice( ERROR, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateIntegrationPointToNodeProperty",
                                 "Properties are not of the same type, nothing was done...");
          return;
       }

    RenumberNodes();

    // required vectors
    vector<double64>  cp_var, n_var;

    // SCALAR PROPERTIES
    if ( n_key.type == SCALAR ) {
         vector<double64> temp( Nodes(), 0. ), sum( Nodes(), 0. );
         const size_t  variable_components(1U);
         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                it=ElementsBegin(); it!=ElementsEnd(); it++ ) {
              // getting the integration-point property
              cp_var.resize( (*it)->IntegrationPoints() );
              for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ ) cp_var[i] = (*it)->Read( i, c_key );
              // extrapolating it to nodes
              n_var.resize( (*it)->Nodes() );                 // components = 1
              (*it)->ExtrapolateIntegrationPointVariableToNodes( variable_components, cp_var, n_var );
              // weighting and accumulating it into vector<double64> for later averaging
              const Point<dim>  barycenter((*it)->BaryCenter());
              for ( size_t i=0U; i<(*it)->Nodes(); i++ )
                {
                   // finding the distance of the node from the barycentre
                   const double64 distance = barycenter.DistanceTo((*it)->N(i)->Coordinate());
                   temp[ (*it)->N(i)->Idx() ] += n_var[i] / distance;
                   sum[  (*it)->N(i)->Idx() ] += 1. / distance;
                }
           }
          // distance weighting and averaging the extrapolated node values
         for ( typename vector<csmp::Node<dim>*>::iterator nit=NodesBegin(); nit!=NodesEnd(); nit++ )
           (*nit)->Store( n_key, makeScalar((*nit)->Status(n_key), temp[ (*nit)->Idx() ] / sum[ (*nit)->Idx() ]) );
      }

    // VECTOR PROPERTIES
    else if ( n_key.type == VECTOR ) {
         // creating a zero-initialized temporary vector<double64>
         VectorVariable<dim>  zero_vec; zero_vec=0.;
         vector<VectorVariable<dim> > cpvec, temp( Nodes(), zero_vec ), sum( Nodes(), zero_vec );
         const size_t  vcomponents(dim);
         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                it=ElementsBegin(); it!=ElementsEnd(); it++ ) {
              // getting the constraint point property
              cpvec.resize( (*it)->IntegrationPoints() );
              for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                (*it)->Read( i, c_key, cpvec[i] );
              // rolling the vector<double64> variables out into linear vector<double64> cp_var
              cp_var.resize( (*it)->IntegrationPoints() * vcomponents );
              for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                for ( size_t j=0U; j<vcomponents; j++ ) cp_var[ i * vcomponents + j ] = cpvec[i][j];
              // extrapolating constraint point property to nodes
              n_var.resize( (*it)->Nodes() * vcomponents ); // components = 1
              (*it)->ExtrapolateIntegrationPointVariableToNodes( vcomponents, cp_var, n_var );

              // weighting and accumulating it into vector<double64> for later averaging
              const Point<dim>  barycenter((*it)->BaryCenter());
              // accumulating result into vector<double64> for later averaging
              for ( size_t i=0U; i<(*it)->Nodes(); i++ )
                {
                   // finding the distance of the node from the barycentre
                   const double64 distance = barycenter.DistanceTo((*it)->N(i)->Coordinate());
                   for ( size_t j=0U; j<vcomponents; j++ ) {
                         temp[ (*it)->N(i)->Idx() ](j) += n_var[ i * vcomponents + j ] / distance;
                         sum[  (*it)->N(i)->Idx() ](j) += 1. / distance;
                     }
                }
           }
          // distance weighting and averaging the extrapolated node values
         for ( typename vector<csmp::Node<dim>*>::iterator it=NodesBegin(); it!=NodesEnd(); it++ )
           {
              // averaging the vector<double64> variable
              for ( size_t j=0U; j<vcomponents; j++ ) {
                   temp[ (*it)->Idx() ](j) /= sum[ (*it)->Idx() ](j);
                   temp[ (*it)->Idx() ].Flag(j) = (*it)->Status( n_key, j );
                }
              // storing it
              (*it)->Store( n_key, temp[ (*it)->Idx() ] );
           }
      }

    // TENSOR PROPERTIES
    if ( n_key.type == TENSOR ) {
          csmp_error.notice( WARNING, "ModelSubDomain<dim,SIMPLEX>::ExtrapolateIntegrationPointToNodeProperty",
                            "distance weighting is not applied; tensor values are simply averaged at the nodes.");
         // creating a zero-initialized temporary vector<double64>
         TensorVariable<dim>  zero_ts; zero_ts=0.;
         vector<TensorVariable<dim> > cpts, temp( Nodes(), zero_ts );
         const size_t  tcomponents(dim * dim);
         for ( typename vector<SIMPLEX<dim>*>::const_iterator
                it=ElementsBegin(); it!=ElementsEnd(); it++ ) {
              // getting the constraint point property
              cpts.resize( (*it)->IntegrationPoints() );
              for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                (*it)->Read( i, c_key, cpts[i] );
              // rolling the tensor variable rows out sequentially into the linear vector<double64> cp_var
              cp_var.resize( (*it)->IntegrationPoints() * tcomponents );
              for ( size_t i=0U; i<(*it)->IntegrationPoints(); i++ )
                for ( size_t j=0U; j<dim; j++ )
                  for ( size_t k=0U; k<dim; k++ ) cp_var[ i * tcomponents + j * dim + k ] = cpts[i](j,k);
              // extrapolating constraint point property to nodes
              n_var.resize( (*it)->Nodes() * tcomponents ); // components = 1
              (*it)->ExtrapolateIntegrationPointVariableToNodes( tcomponents, cp_var, n_var );
              // accumulating result into tensor variable vector<double64> for later averaging
              for ( size_t i=0U; i<(*it)->Nodes(); i++ )
                for ( size_t j=0U; j<dim; j++ )
                  for ( size_t k=0U; k<dim; k++ )
                    temp[ (*it)->N(i)->Idx() ](j,k) += n_var[ i * tcomponents + j * dim + k ];
           }
         // averaging the resulting node property and storing it
         for ( typename vector<csmp::Node<dim>*>::iterator
               it=NodesBegin(); it!=NodesEnd(); it++ ) {
              // averaging the tensor variable
              for ( size_t j=0U; j<dim; j++ )
                for ( size_t k=0U; k<dim; k++ )
                  temp[ (*it)->Idx() ](j,k) /= static_cast<double64>( (*it)->Parents() );
              // storing it
              (*it)->Store( n_key, temp[ (*it)->Idx() ] );
           }
      }

 } // end ExtrapolateIntegrationPointToNodeProperty








// CALCULATIONS


/**

Takes the arithmetic (not element size weighted) average of all property values of nodes, integration points,
or elements that belong to the region, depending on where the target property is placed.

@param prop The name of the property which shall be averaged.
*/
template<size_t dim, template<size_t> class SIMPLEX>
double64  ModelSubDomain<dim,SIMPLEX>::Average( const char* prop ) const
 {
    csmp::Index  idx = pref_.StorageKey(prop);
    size_t       counter(0U);

     if ( (elmt_vec_.empty()) ) {
          throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::RegionPropertyAverage",
                                        "Region is empty");
       }

    switch( idx.place )
      {
        case REGION:
        case BOUNDARY:
        case SPLIT_BOUNDARY:
            if ( idx.type == SCALAR ) {
                 return this->Read( idx );
              }
            if ( idx.type == VECTOR ) {
                 VectorVariable<dim>  vc;
                 this->Read( idx, vc );
                 return vc.Length();
              }
            if ( idx.type == TENSOR ) {
                 VectorVariable<dim>  evals;
                 TensorVariable<dim>  ts;
                 this->Read( idx, ts );
                 ts.EigenValues( evals );
                 double64  avg(0.);
                 for ( size_t j=0U; j<dim; j++ ) avg += evals[j];
                 return avg / static_cast<double64>(dim);
              }
            break;

        case ELEMENT:
            if ( idx.type == SCALAR ) {
                 ScalarVariable  sc;
                 double64  avg(0.);
                 for ( typename vector<SIMPLEX<dim>*>::const_iterator
                       eit=elmt_vec_.begin(); eit!=elmt_vec_.end(); eit++ ) {
                      (*eit)->Read( idx, sc );
                      avg += sc();
                   }
                 return avg / static_cast<double64>(elmt_vec_.size());
              }
            if ( idx.type == VECTOR ) {
                 VectorVariable<dim>  vc;
                 double64  avg(0.);
                 for ( typename vector<SIMPLEX<dim>*>::const_iterator
                       eit=elmt_vec_.begin(); eit!=elmt_vec_.end(); eit++ ) {
                      (*eit)->Read( idx, vc );
                      avg += vc.Length();
                   }
                 return avg /= static_cast<double64>(elmt_vec_.size());
              }
            if ( idx.type == TENSOR ) {
                 VectorVariable<dim>  evals;
                 TensorVariable<dim>  ts;
                 double64  avg(0.);
                 for ( typename vector<SIMPLEX<dim>*>::const_iterator
                       eit=elmt_vec_.begin(); eit!=elmt_vec_.end(); eit++ ) {
                      (*eit)->Read( idx, ts );
                      ts.EigenValues( evals );
                      double64 ts_avg(0.);
                      for ( size_t j=0U; j<dim; j++ ) ts_avg += evals[j];
                      avg += ts_avg / static_cast<double64>(dim);
                   }
                 return avg / static_cast<double64>(elmt_vec_.size());
              }
          break;
        case ELEMENT_INTEGRATION_POINT:
            assert( pref_.VariableCount(ELEMENT_INTEGRATION_POINT) > 0U );
            counter = 0U;
            if ( idx.type == SCALAR ) {
                 ScalarVariable  sc;
                 double64  avg(0.);
                 for ( typename vector<SIMPLEX<dim>*>::const_iterator
                       eit=elmt_vec_.begin(); eit!=elmt_vec_.end(); eit++ )
                   for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                        (*eit)->Read( i, idx, sc );
                        counter++;
                        avg += sc();
                     }
                 return avg / static_cast<double64>(elmt_vec_.size());
              }
            if ( idx.type == VECTOR ) {
                 VectorVariable<dim>  vc;
                 double64  avg(0.);
                 for ( typename vector<SIMPLEX<dim>*>::const_iterator
                       eit=elmt_vec_.begin(); eit!=elmt_vec_.end(); eit++ )
                   for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                        (*eit)->Read( i, idx, vc );
                        counter++;
                        avg += vc.Length();
                     }
                 return avg /= static_cast<double64>(counter);
              }
            if ( idx.type == TENSOR ) {
                 VectorVariable<dim>  evals;
                 TensorVariable<dim>  ts;
                 double64  avg(0.);
                 for ( typename vector<SIMPLEX<dim>*>::const_iterator
                       eit=elmt_vec_.begin(); eit!=elmt_vec_.end(); eit++ )
                   for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                        (*eit)->Read( i, idx, ts );
                        counter++;
                        ts.EigenValues( evals );
                        double64 ts_avg(0.);
                        for ( size_t j=0U; j<dim; j++ ) ts_avg += evals[j];
                        avg += ts_avg / static_cast<double64>(dim);
                     }
                 return avg / static_cast<double64>(counter);
              }
          break;
        case NODE:
            if ( idx.type == SCALAR ) {
                 ScalarVariable  sc;
                 double64  avg(0.);
                 for ( typename vector<csmp::Node<dim>*>::const_iterator
                       it=node_vec_.begin(); it!=node_vec_.end(); it++ ) {
                      (*it)->Read( idx, sc );
                      avg += sc();
                   }
                 return avg / static_cast<double64>(node_vec_.size());
              }
            if ( idx.type == VECTOR ) {
                 VectorVariable<dim>  vc;
                 double64  avg(0.);
                 for ( typename vector<csmp::Node<dim>*>::const_iterator
                       it=node_vec_.begin(); it!=node_vec_.end(); it++ ) {
                      (*it)->Read( idx, vc );
                      avg += vc.Length();
                   }
                 return avg /= static_cast<double64>(node_vec_.size());
              }
            if ( idx.type == TENSOR ) {
                 VectorVariable<dim>  evals;
                 TensorVariable<dim>  ts;
                 double64  avg(0.);
                 for ( typename vector<csmp::Node<dim>*>::const_iterator
                       it=node_vec_.begin(); it!=node_vec_.end(); it++ ) {
                      (*it)->Read( idx, ts );
                      ts.EigenValues( evals );
                      double64 ts_avg(0.);
                      for ( size_t j=0U; j<dim; j++ ) ts_avg += evals[j];
                      avg += ts_avg / static_cast<double64>(dim);
                   }
                 return avg / static_cast<double64>(node_vec_.size());
              }
         break;
         default:
            throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::RegionAverage:", "property placement not handled yet.");
       }

    return std::numeric_limits<double64>::quiet_NaN();

 } // end Average




/**

Provided that property 'a' is a scalar node variable and property 'b' is a
vector<double64> variable placed on the element, CopyGradientOfProperty_A_To_B() will
calculate the gradient of property 'a' for each element and assign the
result to the element variable 'b'.

@section arguments Input Arguments

The method takes two name strings as arguments. The first string specifies
the scalar node variable of which the gradient will be calculated and the
second string identifies the vector<double64> variable placed on the element which
will store the calculated gradient of property 'a'.

@section application Application

The definition of many geophysical or geochemical interrelations requires
a knowledge of property gradients which can be calculated with this method,
using the element interpolation functions of type of finite element which
was used to build the mesh. The method can also be used for post-processing
following the application of algorithms.

@section messages Messages

CopyGradientOfProperty_A_To_B() will return without completing any
calculations, if the input variables do not comply with the specifications
outlined above.

*/
template<size_t dim, template<size_t> class SIMPLEX>
bool  ModelSubDomain<dim,SIMPLEX>::CopyGradientOfProperty_A_To_B( const char* a, const char* b )
 {
    csmp::Index  a_key = pref_.StorageKey(a),
                 b_key = pref_.StorageKey(b);

    // 1. testing variable A for suitability
    // -------------------------------------
    if ( !pref_.IsDefined(a) ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property A does not exist. Now exciting..." << endl;
         return false;
      }
    if ( a_key.place != NODE ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property A is not a NODE variable. Gradient can't be calculated..." << endl;
         return false;
      }
    if ( a_key.type == TENSOR ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property A is a TENSOR variable. Nothing is done..." << endl;
         return false;
      }

    // 2. testing variable B for suitability
    // -------------------------------------
    if ( !pref_.IsDefined(b) ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property B does not exist. Now exciting..." << endl;
         return false;
      }
    if ( b_key.place != ELEMENT and b_key.place != ELEMENT_INTEGRATION_POINT ) {
         cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cerr << "Property B is not an ELEMENT or ELEMENT_INTEGRATION_POINT variable. Gradient can't be calculated..." << endl;
         return false;
      }

    DenseMatrix<DM_MIN>  DN;

    // 3. SCALAR to VECTOR calculations
    // --------------------------------
    if ( a_key.type == SCALAR )
      {
         if ( b_key.type != VECTOR ) {
               cerr << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
               cerr << "Property B is not a VECTOR variable. Nothing is done..." << endl;
               return false;
            }

         vector<ScalarVariable >  SC;
         VectorVariable<dim>      vc;
 
        if ( b_key.place == ELEMENT )
          {
             for ( typename vector<SIMPLEX<dim>*>::iterator
                  eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
               {
                  if ( (*eit)->FE()->UsesLocalCoordinates() ) (*eit)->dN_AtBaryCenter( DN );
                  else (*eit)->dN( DN );
                  (*eit)->NodePropertyVector( a_key, SC );

                  vc = 0.;
                  const size_t nodes((*eit)->Nodes());
                  for ( size_t i=0U; i<nodes; i++ )
                    for ( size_t j=0U; j<dim; j++ )
                      vc(j) += DN(j,i) * SC[i]();
     
                  (*eit)->Store( b_key, vc );
               }
             return true;
           }
        if ( b_key.place == ELEMENT_INTEGRATION_POINT )
          {
             for ( typename vector<SIMPLEX<dim>*>::iterator
                  eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
               {
                  (*eit)->NodePropertyVector( a_key, SC );

                  const size_t ipoints((*eit)->IntegrationPoints());
                  const size_t nodes((*eit)->Nodes());
                  for ( size_t i=0U; i<ipoints; ++i )
                    {
                       (*eit)->dN_AtIntegrationPoint( DN, i );
                       vc = 0.;
                       for ( size_t i=0U; i<nodes; i++ )
                         for ( size_t j=0U; j<dim; j++ )
                           vc(j) += DN(j,i) * SC[i]();

                       (*eit)->Store( i, b_key, vc );
                    }
               }
             return true;
           }
       } // end a.type=scalar
   
   
    // 4. VECTOR to TENSOR calculations
    // --------------------------------
    if ( a_key.type == VECTOR )
      {
         if ( b_key.type != TENSOR ) {
              cerr << "ModelSubDomain<"<<dim<<">::CopyGradientOfProperty_A_To_B: "<< endl;
              cerr << "Property A is a VECTOR so Property B should be a TENSOR variable. Nothing is done..." << endl;
              return false;
           }

         vector<VectorVariable<dim> > VC;
         TensorVariable<dim>          ts;

         if ( b_key.place == ELEMENT )
           {
             for ( typename vector<SIMPLEX<dim>*>::iterator
                   eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
               {
                  if ( (*eit)->FE()->UsesLocalCoordinates() ) (*eit)->dN_AtBaryCenter( DN );
                  else (*eit)->dN( DN );
                  (*eit)->NodePropertyVector( a_key, VC );

                  ts = 0.;

                  // the gradients become rows of the tensor
                  const size_t nodes((*eit)->Nodes());
                  for ( size_t n=0U; n<nodes; n++ )
                    for ( size_t i=0U; i<dim; i++ )
                      for ( size_t j=0U; j<dim; j++ )
                        ts(i,j) += DN(i,n) * VC[n][j];

                  // saving the resulting vector<double64>
                  (*eit)->Store( b_key, ts );
               }
             return true;
           }
         if ( b_key.place == ELEMENT_INTEGRATION_POINT )
           {
             for ( typename vector<SIMPLEX<dim>*>::iterator
                   eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
               {
                  (*eit)->NodePropertyVector( a_key, VC );
                  // the gradients become rows of the tensor
                  const size_t ipoints((*eit)->IntegrationPoints());
                  const size_t nodes((*eit)->Nodes());
                  for ( size_t i=0U; i<ipoints; ++i )
                    {
                      (*eit)->dN_AtIntegrationPoint( DN, i );
                      ts = 0.;
                      for ( size_t n=0U; n<nodes; n++ )
                        for ( size_t i=0U; i<dim; i++ )
                          for ( size_t j=0U; j<dim; j++ )
                            ts(i,j) += DN(i,n) * VC[n][j];

                      // saving the resulting vector<double64>
                      (*eit)->Store( i, b_key, ts );
                    }
               }
             return true;
           }


      }

   return true;

 }  // end CopyGradientOfProperty_A_To_B
 
/* ORIGINAL

template<size_t dim, template<size_t> class SIMPLEX>
bool  ModelSubDomain<dim,SIMPLEX>::CopyGradientOfProperty_A_To_B( const char* a, const char* b )
 {
    csmp::Index  a_key = pref_.StorageKey(a),
                 b_key = pref_.StorageKey(b);

    // 1. testing variable A for suitability
    // -------------------------------------
    if ( !pref_.IsDefined(a) )
      {
         cout << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cout << "Property A does not exist. Now exciting..." << endl;
         return false;
      }
    if ( a_key.place != NODE )
      {
         cout << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cout << "Property A is not a NODE variable. Gradient can't be calculated..." << endl;
         return false;
      }
    if ( a_key.type == TENSOR )
      {
         cout << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cout << "Property A is a TENSOR variable. Nothing is done..." << endl;
         return false;
      }

    // 2. testing variable B for suitability
    // -------------------------------------
    if ( !pref_.IsDefined(b) )
      {
         cout << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cout << "Property B does not exist. Now exciting..." << endl;
         return false;
      }
    if ( b_key.place != ELEMENT )
      {
         cout << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
         cout << "Property B is not an ELEMENT variable. Gradient can't be calculated..." << endl;
         return false;
      }

    DenseMatrix<DM_MIN>  DN;

    if ( a_key.type == SCALAR ) {
        if ( b_key.type != VECTOR )
          {
             cout << "ModelSubDomain<"<< dim << ">::CopyGradientOfProperty_A_To_B: "<< endl;
             cout << "Property B is not a VECTOR variable. Nothing is done..." << endl;
             return false;
          }

      vector<ScalarVariable >  SC;
      VectorVariable<dim>      vc;

        for ( typename vector<SIMPLEX<dim>*>::iterator
              eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
          {
             if ( (*eit)->FE()->UsesLocalCoordinates() ) (*eit)->dN_AtBaryCenter( DN );
             else (*eit)->dN( DN );
             (*eit)->NodePropertyVector( a_key, SC );

             vc = 0.;

             for ( size_t i=0U; i<(*eit)->Nodes(); i++ )
               for ( size_t j=0U; j<dim; j++ )
                 vc(j) += DN(j,i) * SC[i]();

             // saving the resulting vector<double64>
             (*eit)->Store( b_key, vc );
          }
      }
    else if ( a_key.type == VECTOR ) {
        if ( b_key.type != TENSOR ) {
             cout << "ModelSubDomain<"<<dim<<">::CopyGradientOfProperty_A_To_B: "<< endl;
             cout << "Property A is a VECTOR so Property B should be a TENSOR variable. Nothing is done..." << endl;
             return false;
          }

      vector<VectorVariable<dim> > VC;
      TensorVariable<dim>          ts;

        for ( typename vector<SIMPLEX<dim>*>::iterator
              eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
          {
             if ( (*eit)->FE()->UsesLocalCoordinates() ) (*eit)->dN_AtBaryCenter( DN );
             else (*eit)->dN( DN );
             (*eit)->NodePropertyVector( a_key, VC );

             ts = 0.;

             // the gradients become rows of the tensor
             for ( size_t n=0U; n<(*eit)->Nodes(); n++ )
               for ( size_t i=0U; i<dim; i++ )
                 for ( size_t j=0U; j<dim; j++ )
                   ts(i,j) += DN(i,n) * VC[n][j];

             // saving the resulting vector<double64>
             (*eit)->Store( b_key, ts );
          }
      }

   return true;

 }  // end CopyGradientOfProperty_A_To_B

*/
 
 
 

template<size_t dim, template<size_t> class SIMPLEX>
void  ModelSubDomain<dim,SIMPLEX>::CopyReplace( const char* from, const char* to )
 {
     if ( pref_.Type(from) == SCALAR ) {
          CopyReplaceVisitor<ScalarVariable,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == VECTOR ) {
          CopyReplaceVisitor<VectorVariable<dim>,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == TENSOR ) {
          CopyReplaceVisitor<TensorVariable<dim>,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == ARRAY ) {
          CopyReplaceVisitor<ArrayVariable,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else if ( pref_.Type(from) == FLAGGEDARRAY ) {
          CopyReplaceVisitor<FlaggedArrayVariable,dim>
            cpvisitor( pref_, from, to );
          Accept( cpvisitor );
       }
     else throw csmp::Exception( ERROR,
                                 "ModelSubDomain<dim,SIMPLEX>::CopyReplace",
                                 "Property type not supported yet");

 } // end CopyReplace




















// SCREEN OUTPUT



/**

Prints the values and flags of the distributed variable of interest on
stdout.

@param prop The name of the physical variable of interest.

@section application Application

To check variable values in small test problems.

@section messages Messages

OutputVariableToScreen() will report the variable type, its placement,
and property index. Then it will print the host object IDs followed
by the variable flags and values.

*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::OutputVariableToScreen( const char* prop ) const
 {
     csmp::Index          prop_key = pref_.StorageKey(prop);
     ScalarVariable       sc;
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;

     string place_name = parsePlacement(prop_key.place);

     cout <<"\n\nModelSubDomain<"<< dim << ">::OutputVariableToScreen: output of ";
     switch(prop_key.type)
      {
         case SCALAR: cout <<"SCALAR variable: '";
           break;
         case VECTOR: cout <<"VECTOR variable: '";
           break;
         case TENSOR: cout <<"TENSOR variable: '";
           break;
         case ARRAY:  cout <<"ARRAY variable: '";
           break;
         case FLAGGEDARRAY:  cout <<"FLAGGEDARRAY variable: '";
           break;
         default:
           cout <<"\nModelSubDomain<"<< dim <<">::OutputVariableToScreen: type of '"<< prop <<"' not recognized."<< endl;
           return;
      }
     cout << prop <<"' placed on the: "<< place_name << endl;

     switch( prop_key.place )
       {
          case NODE:
              for ( typename vector<csmp::Node<dim>*>::const_iterator
                    nit=NodesBegin(); nit!=NodesEnd(); nit++ )
                {
                  cout <<"\nIdx: " << (*nit)->Idx() <<"\t\t";
                  cout << string(parseBoundary((*nit)->AtBoundary())) <<" ";
                  switch (prop_key.type)
                    {
                       case SCALAR:
                           (*nit)->Read( prop_key, sc );
                           cout << sc;
                         break;
                       case VECTOR:
                           (*nit)->Read( prop_key, vc );
                           cout << vc;
                         break;
                       case TENSOR:
                           (*nit)->Read( prop_key, ts );
                           cout << ts;
                         break;
                       case ARRAY: {
                            ArrayVariable ary;
                            (*nit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            (*nit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::OutputVariableToScreen:",
                                               "type of Node variable not recognized." );
                    }
                }
              break;
          case ELEMENT_INTEGRATION_POINT:
              for ( typename vector<SIMPLEX<dim>*>::const_iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                {
                  cout <<"\nParent element Idx: " << (*eit)->Idx() <<"\t\t";
                  for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                    switch (prop_key.type)
                      {
                         case SCALAR:
                             (*eit)->Read( i, prop_key, sc ); cout << i <<":"<< sc;
                           break;
                         case VECTOR:
                             (*eit)->Read( i, prop_key, vc ); cout << i <<":"<< vc;
                           break;
                         case TENSOR:
                             (*eit)->Read( i, prop_key, ts ); cout << i <<":"<< ts;
                           break;
                       case ARRAY: {
                            ArrayVariable ary;
                            (*eit)->Read( i, prop_key, ary );
                            ary.Out();
                         }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            (*eit)->Read( i, prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::OutputVariableToScreen:",
                                               "type of Element integration point variable not recognized." );
                      }
                }
              break;
           case ELEMENT:
              for ( typename vector<SIMPLEX<dim>*>::const_iterator
                    eit=ElementsBegin(); eit!=ElementsEnd(); eit++ )
                {
                  cout <<"\nIdx: " << (*eit)->Idx() <<"\t\t";
                  /// Roman, 2014 (Face&InterFace): Should Face contain AtBoundary flag?
                  //cout << string(parseBoundary((*eit)->AtBoundary())) <<" ";
                  switch (prop_key.type)
                    {
                       case SCALAR: (*eit)->Read( prop_key, sc ); cout << sc;
                         break;
                       case VECTOR: (*eit)->Read( prop_key, vc ); cout << vc;
                         break;
                       case TENSOR: (*eit)->Read( prop_key, ts ); cout << ts;
                         break;
                       case ARRAY: {
                            ArrayVariable ary;
                            (*eit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            (*eit)->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::OutputVariableToScreen:",
                                               "type of Element variable not recognized." );
                    }
                }
              break;
           case SPLIT_BOUNDARY:
                  cout <<"\nSplitBoundary: ";
                  switch (prop_key.type)
                    {
                       case SCALAR: this->Read( prop_key, sc ); cout << sc;
                         break;
                       case VECTOR: this->Read( prop_key, vc ); cout << vc;
                         break;
                       case TENSOR: this->Read( prop_key, ts ); cout << ts;
                         break;
                       case ARRAY: {
                            ArrayVariable ary;
                            this->Read( prop_key, ary );
                            ary.Out();
                          }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            this->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::OutputVariableToScreen:",
                                               "type of split boundary variable not recognized." );
                    }
                  break;
           case BOUNDARY:
                {
                  cout <<"\nRegion: ";
                  switch (prop_key.type)
                    {
                       case SCALAR: this->Read( prop_key, sc ); cout << sc;
                         break;
                       case VECTOR: this->Read( prop_key, vc ); cout << vc;
                         break;
                       case TENSOR: this->Read( prop_key, ts ); cout << ts;
                         break;
                       case ARRAY: {
                            ArrayVariable ary;
                            this->Read( prop_key, ary );
                            ary.Out();
                          }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            this->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::OutputVariableToScreen:",
                                               "type of boundary variable not recognized." );
                    }
                }
              break;
           case REGION:
                {
                  cout <<"\nRegion: ";
                  switch (prop_key.type)
                    {
                       case SCALAR: this->Read( prop_key, sc ); cout << sc;
                         break;
                       case VECTOR: this->Read( prop_key, vc ); cout << vc;
                         break;
                       case TENSOR: this->Read( prop_key, ts ); cout << ts;
                         break;
                       case ARRAY: {
                            ArrayVariable ary;
                            this->Read( prop_key, ary );
                            ary.Out();
                          }
                         break;
                       case FLAGGEDARRAY: {
                            FlaggedArrayVariable ary;
                            this->Read( prop_key, ary );
                            ary.Out();
                         }
                         break;
                       default:
                         throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::OutputVariableToScreen:",
                                               "type of region variable not recognized." );
                    }
                }
              break;
          default:
               cerr << "\nModelSubDomain<"<< dim <<">::OutputVariableToScreen: The property ";
               cerr << prop <<"  "<< place_name <<" could not be retrieved from the Region";
               break;
       }
    cout << endl;
    cout.flush();

 } // end OutputVariableToScreen







/**
    prints state of the object manifest in:
    
    subdomain_name_ - passed down when region is created so that it can be referred to

    elmt_vec_ - doubly sorted, interior elements first
    
    bd_face_vec_ - face-IDs of elements at the domain boundary as in second segment of elmt_vec_

    node_vec_ - doubly sorted, interior nodes first
    
    first_bd_node_ - first node in the boundary range

    verbose_ - or not

*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::Out() const
 {
    cout <<"\nModelSubDomain<dim,SIMPLEX>::Out(): name: '"<< subdomain_name_ <<"'";
    cout <<" member elements: interior="<< InteriorElements();
    cout <<", boundary="<< elmt_vec_.size()-InteriorElements() <<": "<< endl;

    for ( typename vector<SIMPLEX<dim>*>::const_iterator
          it=elmt_vec_.begin(); it!=elmt_vec_.end(); it++ ) {
         if ( (*it) == NULL )
           throw csmp::Exception( ERROR, "ModelSubDomain<dim,SIMPLEX>::Out",
                                 "member element pointer not initialised");
//         else (*it)->Out();
      }

    cout <<"\n\n edge elements and their edge faces (current local numbering): "<< endl;
    vector<vector<ONE_BYTE_NUMBER> >::const_iterator  bit(bd_face_vec_.begin());
    for ( size_t i=InteriorElements(); i<elmt_vec_.size(); i++, bit++ ) {
         cout <<"\nelement "<< i <<": edge face numbers: ";
         for ( vector<ONE_BYTE_NUMBER>::const_iterator
               ft=(*bit).begin(); ft!=(*bit).end(); ft++ ) cout << (*ft) <<" ";
      }

    cout <<"\n\n edge nodes: "<< node_vec_.size() - first_bd_node_ <<" (current local numbering):"<< endl;
    for ( size_t i=first_bd_node_; i<node_vec_.size(); i++ ) {
         if ( node_vec_[i] == NULL )
           throw csmp::Exception( ERROR, "ModelSubDomain<dim>::Out", "member node pointer not initialised.");
         else cout << node_vec_[i]->Idx() <<" ";
      }

    cout << endl;
 }




/**
    Writes data block containing all information that is needed to reconstruct a ModelSubDomain
    
    @note template template parameters for functions are not allowed
*/
template<size_t dim, template<size_t> class SIMPLEX>
void ModelSubDomain<dim,SIMPLEX>::WriteDomainIndexesToBinaryFile( FILE* fp ) const
 {
    assert( fp != nullptr );
   
    // 1. writing name of the region
    skm_C_fwrite( fp, Name().c_str() );
   
    // 2. writing the interior element records of the region
    std::vector<size_t> IDs( distance(ElementsBegin(), PerimeterElementsBegin()) );
    transform( ElementsBegin(), PerimeterElementsBegin(),
               IDs.begin(), []( const SIMPLEX<dim>* const ptr ){ return ptr->Idx(); } );
    skm_C_fwrite( fp, IDs );

    // 3. writing the perimeter element records of the region
    IDs.resize( distance(PerimeterElementsBegin(), ElementsEnd()) );
    transform( PerimeterElementsBegin(), ElementsEnd(),
               IDs.begin(), []( const SIMPLEX<dim>* const ptr ){ return ptr->Idx(); } );
    skm_C_fwrite( fp, IDs );
   
    // 4. writing the boundary faces
    // -----------------------------
    /* 
        since this is vector of vectors predominated by single value entries,
        it is collapsed into a flat vector in which all entries that refer to 
        multiple values per element are prefaced by a negative numer that indicates
        how many multiple faces per element follow, for example
        1 5  5 3 -2 6 2 3 3 5 6
                    ^^^          marking the 2 local face indices that relate to an element that has
        2 faces on the model boundary.
    */
    std::vector<int8> faceIDs; // signed byte -127..128: small because only the local face IDs are needed
    faceIDs.reserve( PerimeterElements() );
    for ( size_t eid(InteriorElements()); eid<Elements(); ++eid ) {
         const size_t perimeter_faces(PerimeterFaces(eid));
         // writing the number of perimeter faces as negative number, but only if there are more than 1 perimeter faces
         if ( perimeter_faces > 1 ) faceIDs.push_back( static_cast<int8>(-perimeter_faces) );
         // writing the perimeter face ids
         for ( size_t j=0U; j<perimeter_faces; ++j )
           faceIDs.push_back( static_cast<int8>(PerimeterFace(eid,j)) );
      }
    skm_C_fwrite( fp, faceIDs );

    // 5. writing the interior nodes
    IDs.resize( distance(NodesBegin(), PerimeterNodesBegin()) );
    transform( NodesBegin(), PerimeterNodesBegin(),
               IDs.begin(), []( const Node<dim>* const ptr ){ return ptr->Idx(); } );
    skm_C_fwrite( fp, IDs );
   
    // 6. writing the perimeter nodes
    IDs.resize( distance(PerimeterNodesBegin(), NodesEnd()) );
    transform( PerimeterNodesBegin(), NodesEnd(),
               IDs.begin(), []( const Node<dim>* const ptr ){ return ptr->Idx(); } );
    skm_C_fwrite( fp, IDs );
   
 } // end WriteDomainIndexesToBinaryFile


/*


cerr <<"\n\n\nregion: "<< Name() <<"\n";
cerr <<"interior nodes";
out( IDs );
cerr <<"\n perimeter nodes";
out( IDs );
*/




/**
    Reads all the data required to fully reconstruct a ModelSubDomain (without search operations)
*/
void readDomainIndexesFromBinaryFile( FILE* fp, SubDomainInfo& info )
 {
    assert( fp != nullptr );
   
    // 1. reading name of the region
    char name[INFO_STRING];
    skm_C_fread( fp, name );
    info.name = name;
    assert( !info.name.empty() );
   
    // 2. reading the interior element records of the region
    skm_C_fread( fp, info.interior_elmts );
    assert( !info.interior_elmts.empty() );
   
    // 3. reading the perimeter element records of the region
    skm_C_fread( fp, info.perimeter_elmts );
    assert( !info.perimeter_elmts.empty() );
   
    // 4. reading the boundary faces
    // -----------------------------
    /* 
        expects flat vector in which all entries that refer to
        multiple values per element are prefaced by a negative numer that indicates
        how many multiple faces per element follow, for example
        1 5  5 3 -2 6 2 3 3 5 6
                    ^^^          marking the 2 local face indices that relate to an element that has
        2 faces on the model boundary.
        where there is no negative number, a single entry is assumed
    */
    std::vector<int8> faceIDs; // signed byte -127..128: small because only the local face IDs are needed
    skm_C_fread( fp, faceIDs );
    assert( !faceIDs.empty() );
 
    if ( !info.perimeter_faces.empty() ) info.perimeter_faces.clear();
    info.perimeter_faces.reserve( faceIDs.size() );
    for ( std::vector<int8>::const_iterator it=faceIDs.begin(); it!=faceIDs.end(); ++it ) {
         const size_t perimeter_faces = ((*it) < 0) ? abs( (*it) ) : 1;
         std::vector<int8> face_ids;
         face_ids.reserve(3);
         for ( size_t i=0U; i<perimeter_faces; ++i ) {
              if ( perimeter_faces > 1 ) it++;
              assert( it != faceIDs.end() );
              face_ids.push_back( (*it) );
           }
         info.perimeter_faces.push_back( move(face_ids) );
      }

    // 5. reading the interior nodes
    skm_C_fread( fp, info.interior_nodes );
    assert( !info.interior_nodes.empty() );
  
    // 6. reading the perimeter nodes
    skm_C_fread( fp, info.perimeter_nodes );
    assert( !info.perimeter_nodes.empty() );
   
 } // end readRegionIndexesFromBinaryFile





/**
    @return returns the number of nodes on the subdomain perimeter which are shared by the subdomain and a given model boundary
*/
template<size_t dim, template<size_t> class SIMPLEX>
size_t ModelSubDomain<dim,SIMPLEX>::SharedPerimeterNodes( typename std::vector<csmp::Node<dim>*>::const_iterator start,
                                                          typename std::vector<csmp::Node<dim>*>::const_iterator end ) const
 {
    if ( start == end ) return 0U;
 
    typename vector<Node<dim>*>::const_iterator  first1(this->PerimeterNodesBegin());
    size_t shared_nodes(0U);

    // comparing the boundary nodes
    while ( first1 != this->NodesEnd() and start != end )
      {
        if ( *first1 < *start ) ++first1;
        else if ( *start < *first1 ) ++start;
        else {
             shared_nodes++;
             first1++;
             start++;
          }
      }

    return shared_nodes;
 }





/**
     Counts and returns the number of nodes shared between the two regions.
     Special attention is paid to the fact the nodes are sorted in two ranges.
     
     @attention function assumes that both subdomains are valid, containing multiple nodes.
 */
template<size_t dim, template<size_t> class SIMPLEX>
size_t  sharedNodes( const ModelSubDomain<dim,SIMPLEX>& g1, const ModelSubDomain<dim,SIMPLEX>& g2 )
 {
    typename vector<Node<dim>*>::const_iterator  first1(g1.NodesBegin());
    typename vector<Node<dim>*>::const_iterator  first2(g2.NodesBegin());
    size_t shared_nodes(0U);

    // comparing the interior nodes
    while ( first1 != g1.PerimeterNodesBegin() and first2 != g2.PerimeterNodesBegin() )
      {
        if ( *first1 < *first2 ) ++first1;
        else if ( *first2 < *first1 ) ++first2;
        else {
             shared_nodes++;
             first1++;
             first2++;
          }
      }

    // comparing the boundary nodes
    first1 = g1.PerimeterNodesBegin();
    first2 = g2.PerimeterNodesBegin();

    while ( first1 != g1.NodesEnd() and first2 != g2.NodesEnd() )
      {
        if ( *first1 < *first2 ) ++first1;
        else if ( *first2 < *first1 ) ++first2;
        else {
             shared_nodes++;
             first1++;
             first2++;
          }
      }

    return shared_nodes;

 } // end sharedNodes

template size_t sharedNodes( const ModelSubDomain<1U,Element>& g1, const ModelSubDomain<1U,Element>& g2 );
template size_t sharedNodes( const ModelSubDomain<2U,Element>& g1, const ModelSubDomain<2U,Element>& g2 );
template size_t sharedNodes( const ModelSubDomain<3U,Element>& g1, const ModelSubDomain<3U,Element>& g2 );


/** 
    As for sharedNodes() - but application is restricted to nodes that sit
    on the perimeter / surface of the region.
*/
template<size_t dim, template<size_t> class SIMPLEX>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,SIMPLEX>& g1, const ModelSubDomain<dim,SIMPLEX>& g2 )
 {
    typename vector<Node<dim>*>::const_iterator  first1(g1.PerimeterNodesBegin());
    typename vector<Node<dim>*>::const_iterator  first2(g2.PerimeterNodesBegin());
    size_t shared_nodes(0U);

    // comparing the boundary nodes
    while ( first1 != g1.NodesEnd() and first2 != g2.NodesEnd() )
      {
        if ( *first1 < *first2 ) ++first1;
        else if ( *first2 < *first1 ) ++first2;
        else {
             shared_nodes++;
             first1++;
             first2++;
          }
      }

    return shared_nodes;

 } // end sharedPerimeterNodes








 template class ModelSubDomain<1U,Element>;
 template class ModelSubDomain<2U,Element>;
 template class ModelSubDomain<3U,Element>;

 template class ModelSubDomain<1U,Face>;
 template class ModelSubDomain<2U,Face>;
 template class ModelSubDomain<3U,Face>;

 template class ModelSubDomain<1U,InterFace>;
 template class ModelSubDomain<2U,InterFace>;
 template class ModelSubDomain<3U,InterFace>;

} // end namespace

