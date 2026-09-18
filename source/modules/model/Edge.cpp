// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  Edge.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 19/10/2022.
//

#include "Edge.h"
#include "Node.h"
#include "Face.h"
#include "Visitor.h"
#include "FiniteElement.h"
#include "FiniteElementManager.h"
#include "FiniteVolumeStencilManager.h"

using namespace std;

namespace csmp {

template<uint32_t DIM>
Edge<DIM>::Edge( FiniteElement* feptr,
                 const FiniteVolumeStencil<DIM>* fvptr,
                 Face<DIM>* const inner_parent,
                 Face<DIM>* const outer_parent, 
                 const LocalVariables& lv,
                 const IntegrationPointVariables& iv )
  : FiniteElementPolicy<DIM,Edge>(feptr),
    FiniteVolumePolicy<DIM,Edge>(fvptr),
    node_connector_( feptr->Nodes(), nullptr ),
    idx_(NULL_IDX),
    innerParent_(inner_parent),
    outerParent_(outer_parent)
 {
    assert( feptr );
    assert( innerParent_ != nullptr );
    if ( innerParent_ == outerParent_ ) {
         cerr <<"\n\Edge<3U>(ctor: of Edge between parents): supplied pointers point to the same element: ";
         cerr << inner_parent->Idx();
         inner_parent->Out();
         innerParent_ = nullptr;
         outerParent_ = nullptr;
         assert( innerParent_ != outerParent_ );
         return;
      }
    
    // finding the face which is shared and assigning the nodes
    // --------------------------------------------------------
    bool  matching_face_found{false};
    const size_t n_faces_inner{innerParent_->Faces()};
    
    // inner parent is always present
    for ( auto i{0U}; i<n_faces_inner; ++i )
      // if the faces match
      if ( inner_parent->Neighbor(i) == innerParent_ )
        {
           inner_parent_face_id_ = i;
           // assigning the nodes
           uint32_t n_count{0U};
           for ( const auto& k : inner_parent->FE()->NodesOfFace(i) )
             node_connector_[n_count++] = inner_parent->N(k);

           matching_face_found = true;
           break;
        }
         
    // outer parent: exists only if Edge is on inside of model
    if ( outerParent_ ) {
         // finding the face number of outer element
         const auto n_faces_outer{outerParent_->Faces()};
         for ( auto j{0U}; j<n_faces_outer; ++j )
           if ( outer_parent->Neighbor(j) == innerParent_ ) {
                outer_parent_face_id_ = j;
                break;
             }
      }
 
    // reporting the failed construction
    if ( !matching_face_found ) {
         cerr <<"\n\nEdge<3U>(ctor: Edge between parents): parent faces "<< inner_parent->Idx();
         cerr <<" and "<< outer_parent->Idx() <<" do not seem to share a face:\n";
         inner_parent->Out();
         outer_parent->Out();
         innerParent_ = nullptr;
         outerParent_ = nullptr;
         return;
      }

    // pointers should point to the adjacent element
    assert( innerParent_->Neighbor(inner_parent_face_id_) == outerParent_ );
    if ( outerParent_ )
      assert( outerParent_->Neighbor(outer_parent_face_id_) == innerParent_ );
    
    // creating local storage for face and face integration point variables
    if ( this->UsesLocalCoordinates() )
        this->ResizePropertyStorage( lv, iv );
    else
        this->ResizePropertyStorage( lv );
      
 } // end (constructor that infers nodes faces from higher-dimensional neighbor faces)



template<uint32_t DIM>
uint32_t  Edge<DIM>::ConnectedNeighbors() const {
    uint32_t interconnected_edges{ 0U };
    if ( innerParent_ ) interconnected_edges++;
    if ( outerParent_ ) interconnected_edges++;
    return interconnected_edges;
 }



/// only constant iterators are provided because the user is not supposed to change the node pr neighbor connectivity (done by MeshManager)
template<uint32_t DIM>
typename vector<Node<DIM>*>::const_iterator Edge<DIM>::NodesBegin() const {
     return node_connector_.begin();
  }
  


template<uint32_t DIM>
typename vector<Node<DIM>*>::const_iterator Edge<DIM>::NodesEnd() const {
     return node_connector_.end();
  }


/// access the nodes that are connected to the Edge
template<uint32_t DIM>
Node<DIM>* const Edge<DIM>::N( uint32_t n_local ) const
 {
    assert( n_local < node_connector_.size() );
    return node_connector_[n_local];
 }
 

/// access the neighbor faces of this face
template<uint32_t DIM>
Edge<DIM>* const Edge<DIM>::Neighbor( uint32_t n ) const
 {
    if ( n == 1U ) return nbor_edge1_;
    return nbor_edge0_;
 }


/// to apply visitors whose application level is Boundary and target is Edge
template<uint32_t DIM>
void Edge<DIM>::Accept( Visitor<DIM>& vis )
 {
    if ( vis.ApplicationTarget() == EDGE ) {
         vis.Visit(this);
         return;
      }
    if ( vis.ApplicationTarget() == NODE ) {
         uint32_t n_nodes{ this->Nodes() };
         for ( auto i{0U}; i< n_nodes; i++ )
             this->N(i)->Accept( vis );
         return;
      }
    throw logic_error("Edge<3U>::Accept: target of visitation unresolved.");

 } // end Accept






/// on-the-fly 0..n-1 numbering stored in a mutable local variable (therefore const)
template<uint32_t DIM>
void   Edge<DIM>::Idx( size_t id ) const { idx_ = id; }

template<uint32_t DIM>
size_t Edge<DIM>::Idx() const { return idx_; }



/// access the higher dimensional elements on either side of face; @attention returns nullptr if outside is not present
template<uint32_t DIM>
Face<DIM>* const Edge<DIM>::Parent( INTERFACE_SIDE side ) const {
     assert( side != MIDDLE );
     if ( side == INSIDE ) return innerParent_;
     return outerParent_;
  }


/// higher-dimensional element located on the side of the face to which the unit normal points; @attention does not exist on model boundary
template<uint32_t DIM>
Face<DIM>* const Edge<DIM>::InnerParent() const { return innerParent_; }

template<uint32_t DIM>
Face<DIM>* const Edge<DIM>::OuterParent() const { return outerParent_; }

/// returns which Face of the higher dimensional inner neighbor Face this Edge shares its nodes with
template<uint32_t DIM>
uint32_t Edge<DIM>::InnerParentFaceID() const { return inner_parent_face_id_; }

template<uint32_t DIM>
uint32_t Edge<DIM>::OuterParentFaceID() const { return outer_parent_face_id_; }



/// returns length of the Edge; method assumes same role as Volume() for the element
template<uint32_t DIM>
double Edge<DIM>::Length() const {
    // done by policy already: this->CoordinateMatrix();
    return this->Volume();
 }


/// the centre of gravity of the element as determined from the corner nodes
template<uint32_t DIM>
Point<DIM> Edge<DIM>::BaryCenter() const {
     Point<DIM> temp( node_connector_[0]->Coordinate() + node_connector_[1]->Coordinate() );
     temp /= 2.;
     return temp;
  }
  
 
template<uint32_t DIM>
void  Edge<DIM>::NodeCoordinateMatrix( DenseMatrix<DM_MIN>& XY ) const
{
    const auto n_nodes{ Nodes() };
    XY.Resize( n_nodes, DIM );
    for ( auto i{0U}; i<n_nodes; ++i )
        XY.AssignRow( i, N(i)->Coordinate() );

} // end NodeCoordinateMatrix
 
 
template<uint32_t DIM>
void Edge<DIM>::Out() const
 {
    cout <<"\n"<<"Edge<"<< DIM <<">::Out: Idx: <<"<< idx_ << endl;
    cout <<"\t"<<" with the nodes: ";
    for ( const auto& nit : node_connector_ ) {
         cout << nit->Idx() <<" ";
      }
    if ( !node_connector_.empty() ) cout << endl;
    if ( innerParent_ )
      cout <<"\t"<<"connected inner parent Face: "<< inner_parent_face_id_ <<": "<< innerParent_->Idx() << endl;
    if ( outerParent_ )
      cout <<"\t"<<"connected outer parent Face: "<< outer_parent_face_id_ <<": "<< outerParent_->Idx() << endl;
    cout <<"\t"<<"Edge neighbors: ";
    if ( nbor_edge0_ )
      cout <<"\t"<<"connected neighbor 0: "<< nbor_edge0_->Idx() << endl;
    if ( nbor_edge1_ )
      cout <<"\t"<<"connected neighbor 1: "<< nbor_edge1_->Idx() << endl;

 } // end Out
  

template class Edge<1U>;
template class Edge<2U>;
template class Edge<3U>;

} // end csmp
