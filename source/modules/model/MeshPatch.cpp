//
//  MeshPatch.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 23/5/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "MeshPatch.h"
#include "MeshManagementUtilities.h"
#include "ErrorHandler.h"
#include "Element.h"
#include "Node.h"
#include "IsoparametricLinearLineElement.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearQuadrilateral.h"


using namespace std;

namespace csmp {
 
 
template<uint32_t dim>
MeshPatch<dim>::~MeshPatch()
 {
    for ( auto ptr : fe_ptrs_ )
      delete ptr.second;
 }

 
template<uint32_t dim>
size_t MeshPatch<dim>::Cells() const
 {
    return elements_.size();
 }


template<uint32_t dim>
size_t MeshPatch<dim>::Nodes() const
 {
    return nodes_.size();
 }


template<uint32_t dim>
size_t MeshPatch<dim>::CellTypes() const
 {
   return fe_ptrs_.size();
 }


template<uint32_t dim>
CELL_SHAPE MeshPatch<dim>::CellGeometry() const
 {
    return cell_dimension_;
 }


template<uint32_t dim>
typename plf::colony<Element<dim>>::iterator  MeshPatch<dim>::ElementsBegin()
  { return elements_.begin(); }
  
template<uint32_t dim>
typename plf::colony<Element<dim>>::iterator MeshPatch<dim>::ElementsEnd()
 { return elements_.end(); }

template<uint32_t dim>
typename plf::colony<Node<dim>>::iterator MeshPatch<dim>::NodesBegin()
 { return nodes_.begin(); }

template<uint32_t dim>
typename plf::colony<Node<dim>>::iterator MeshPatch<dim>::NodesEnd()
 { return nodes_.end(); }


/**
    Creates an interconnected element patch from the shared faces of the supplied elements.
    There is no property storage associated with these elements from the start.
    The elements are rather barebone.
    
    @param input_elmts pairs of highest-dim elements that share a face that will be turned into a new lower-dimensional element
    @param copy_nodes  the Elements already share nodes, but light copies of these are made when true
    @param perimeter_node_ptrs whether copies are made or not, the method returns vector of pointers of the perimeter nodes of the new element patch
    @return the number of new elements in the patch
*/
template<uint32_t dim> ///<                                                                        inside element and its outside face number           outside element and inside face number
size_t MeshPatch<dim>::BuildInterveningPatch( const vector<pair<pair<Element<dim>*,uint32_t>,pair<Element<dim>*,uint32_t> > >& input_elmts,
                                              bool copy_nodes, vector<Node<dim>*>  perimeter_node_ptrs )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
    const size_t  n_elmts_to_create{ input_elmts.size() };
    if ( n_elmts_to_create == 0U ) {
         csmp_error.notice( ERROR, "MeshPatch<dim>::BuildInterveningPatch",
                                   "no interface elements were provided as input; nothing was done");
         return 0U;
      }
    
    for ( auto& it : input_elmts )
      {
         if constexpr ( dim == 3U ) assert( it.first.first->IsVolume() );
         if constexpr ( dim == 2U ) assert( it.first.first->IsSurface() );

         // 1. constructing the most basic Element objects from the information for the intervening faces
         // ---------------------------------------------------------------------------------------------
         const CSMP_FEM_TYPE etype = it.first.first->FE_Type();
         auto                fit   = fe_ptrs_.end();
         auto                eit   = elements_.end();
         
         switch( etype ) {
             case ISOPARAMETRIC_LINEAR_TRIANGLE:
                  if ( (fit=fe_ptrs_.find(ISOPARAMETRIC_LINEAR_TRIANGLE)) != fe_ptrs_.end() )
                    eit = elements_.insert( Element<dim>( (*fit).second ) );
                  else {
                       auto fem_it = fe_ptrs_.insert( make_pair(ISOPARAMETRIC_LINEAR_TRIANGLE,new IsoparametricLinearTriangle(dim)) );
                       eit = elements_.insert( Element<dim>( (*fem_it.first).second ) );
                    }
               break;
             case ISOPARAMETRIC_LINEAR_QUADRILATERAL:
                  if ( (fit=fe_ptrs_.find(ISOPARAMETRIC_LINEAR_QUADRILATERAL)) != fe_ptrs_.end() )
                    eit = elements_.insert( Element<dim>( (*fit).second ) );
                  else {
                       auto fem_it = fe_ptrs_.insert( make_pair(ISOPARAMETRIC_LINEAR_QUADRILATERAL,new IsoparametricLinearQuadrilateral(dim)) );
                       eit = elements_.insert( Element<dim>( (*fem_it.first).second ) );
                    }
               break;
             case ISOPARAMETRIC_LINEAR_BAR:
                  if ( (fit=fe_ptrs_.find(ISOPARAMETRIC_LINEAR_BAR)) != fe_ptrs_.end() )
                    eit = elements_.insert( Element<dim>( (*fit).second ) );
                  else {
                       auto fem_it = fe_ptrs_.insert( make_pair(ISOPARAMETRIC_LINEAR_BAR,new IsoparametricLinearLineElement(dim)) );
                       eit = elements_.insert( Element<dim>( (*fem_it.first).second ) );
                    }
               break;
             default: {
                  cerr <<"\n\n"<< parseFiniteElementType( etype );
                  csmp_error.notice( ERROR, "MeshPatch<dim>::BuildInterveningPatch", "element type not recognised; nothing could be done.");
                  return 0U;
               }
           }
         
         // 2. connecting nodes
         // -------------------
         Element<dim>& elmt = (*eit);
         vector<uint32_t> fnids;
         elmt.FE()->NodesOfFace( it.first.second, fnids );
         uint32_t nd_count{0U};
         for ( auto i : fnids )
           elmt.Assign( nd_count++, it.first.first->N(i) );
      }
      
      
    // 3. connecting neighbors
    // -----------------------
    if constexpr ( dim == 3U ) {
           // creating search keys from the corner nodes of the element faces
           // corner-nodes      elements that share face and their face id
           map<set<Node<3>*>,map<Element<3>*,uint32_t> >  elmt_pairs;
           // pairing the cells up in the search map
           for ( auto& it : elements_ )
             {
                const size_t n_faces{ it.Faces() };
                for ( auto face{0}; face < n_faces; ++face ) {
                     // trying to insert cell into the map using a search key of node pointers
                     auto pit = elmt_pairs.insert( make_pair( it.CornerNodesOfFace(face), map<Element<3>*,uint32_t>{{&it,face}} ) );
                     // if the face record already exists, the new element pointer - face is added to it
                     if ( pit.second == false )
                       (*pit.first).second.insert( make_pair( &it, face ) );
                       
                     // nulling the current neighbor connectivity of the elements if any
                     it.Assign( face, static_cast<Element<3>*>(nullptr) );
                  }
             }
           // processing the results, connecting the cells to one another
           for ( auto& it : elmt_pairs ) {
                const size_t n_face_nbors{ it.second.size() };
                // if there is just a single matching neighbor
                if ( n_face_nbors == 2 ) {
                     Element<3>* const ptr1  = (*it.second.begin()).first;
                     Element<3>* const ptr2  = (*it.second.rbegin()).first;
                     assert( ptr1 != nullptr );
                     assert( ptr2 != nullptr );
                     const uint32_t face_e1 = (*it.second.begin()).second;
                     const uint32_t face_e2 = (*it.second.rbegin()).second;
                     ptr1->Assign( face_e1, ptr2 );
                     ptr2->Assign( face_e2, ptr1 );
                  }
                // else this is a manifold and two most suitable neighbors must be found
                else if ( n_face_nbors > 2 ) {
                     (*it.second.begin()).first->Out();
                     (*it.second.rbegin()).first->Out();
                     csmp_error.notice( ERROR, "creatInterveningElementPatchFrom", "Elements on either side of the face cannot be matched");
                  }
       } // end 3D
      
    }
 
  // 4. identifying perimeter nodes
  // ------------------------------
  if ( !perimeter_node_ptrs.empty() ) perimeter_node_ptrs.clear();
  perimeter_node_ptrs.reserve( elements_.size() ); // rough guess
  
  for ( const auto& it : elements_ )
    for ( auto i{0U}; i<it.Neighbors(); i++ )
      if ( it.Neighbor(i) == nullptr ) {
           vector<uint32_t> fnids;
           it.FE()->NodesOfFace( i, fnids );
           for ( auto j : fnids )
             perimeter_node_ptrs.push_back( it.N(j) );
        }
        
  sort( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end() );
  perimeter_node_ptrs.erase( unique( perimeter_node_ptrs.begin(), perimeter_node_ptrs.end()), perimeter_node_ptrs.end() );
  perimeter_node_ptrs.shrink_to_fit();
  

  // 2. copying pre-existing nodes if so requested
  // ---------------------------------------------
  if ( copy_nodes ) {
       // 2.1 making a unique map of pointers to old and new nodes
       map<Node<dim>*,Node<dim>*>  node_mapping;
       for ( auto& it : elements_ )
         for ( auto i{0U}; i<it.Nodes(); i++ ) {
               auto inserted = node_mapping.insert( make_pair( it.N(i), nullptr ) );
               // if an extra node is required
               if ( inserted.second == true )
                   (*inserted.first).second = &(*nodes_.insert( Node<dim>( nodes_.size(),
                                                                           it.N(i)->Coordinate(),
                                                                           LocalVariables(),
                                                                           it.N(i)->AtBoundary() ) ));
          }
      // 2.2 replacing the existing nodes with the new locally stored ones
      // (can't be done in one loop as it would confuse old and new nodes)
       for ( auto& it : elements_ )
         for ( auto i{0U}; i<it.Nodes(); i++ ) {
              // finding the replacement nodes in the map
              auto nit = node_mapping.find( it.N(i) );
              assert( nit != node_mapping.end() ); // assert that the node has been found
              it.Assign( i,( *nit).second );
           }
   }

  return elements_.size();
    
} // end creatInterveningElementPatchFrom

template class MeshPatch<3U>;

} // end csmp
