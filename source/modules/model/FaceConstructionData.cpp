//
//  FaceConstructionData.cpp
//
//  Created by Stephan Matthai on 28/03/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "FaceConstructionData.h"
#include "Element.h"
#include "Index.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
FaceConstructionData<dim>::FaceConstructionData( Element<dim>& parent_element,
                                                 pair<Element<dim>*,Element<dim>*>& neighbors,
                                                 pair<uint32_t,uint32_t>& nbor_faces,
                                                 pair<long,long>& materials,
                                                 long material )
 : dim_m1_element_(parent_element),
   neighbors_(neighbors),
   nbor_faces_(nbor_faces),
   materials_(materials),
   material_(material),
   patch_number_(numeric_limits<uint32_t>::max())
 {
 }


template<uint32_t dim>
FaceConstructionData<dim>::FaceConstructionData( Element<dim>&  parent_element,
                                                 pair<Element<dim>*,Element<dim>*>& neighbors,
                                                 pair<uint32_t,uint32_t>& nbor_faces,
                                                 pair<long,long>& materials,
                                                 long material,
                                                 uint32_t patch_number )
 : dim_m1_element_(parent_element),
   neighbors_(neighbors),
   nbor_faces_(nbor_faces),
   materials_(materials),
   material_(material),
   patch_number_(patch_number)
 {
 }
 
 
  
template<uint32_t dim>
FaceConstructionData<dim>::FaceConstructionData( const FaceConstructionData& fcd )
 : dim_m1_element_(fcd.dim_m1_element_),
   neighbors_(fcd.neighbors_),
   nbor_faces_(fcd.nbor_faces_),
   materials_(fcd.materials_),
   material_(fcd.material_),
   patch_number_(fcd.patch_number_)
 {
 }



/// inner (first) and out (second) higher-order parent element
template<uint32_t dim>
pair<Element<dim>*,Element<dim>*> FaceConstructionData<dim>::NeighborElements() const
 {
     return neighbors_;
 }

/**
    lower-dimensional parent element from which face data were recovered in part
    
    @attention this pointer is non-const because it has to be null-assignable
*/
template<uint32_t dim>
Element<dim>* FaceConstructionData<dim>::LowerDimElement()
 {
    return &dim_m1_element_;
 }


/// on the opposite site of the outward pointing normal
template<uint32_t dim>
Element<dim>* const FaceConstructionData<dim>::InnerElement() const
 {
    return neighbors_.first;
 }
 
 
template<uint32_t dim>
uint32_t FaceConstructionData<dim>::InnerElementFace() const
 {
    return nbor_faces_.first;
 }
 

/// on the side to which the normal points to
template<uint32_t dim>
Element<dim>* const FaceConstructionData<dim>::OuterElement() const
 {
    return neighbors_.second;
 }
 
 
template<uint32_t dim>
uint32_t FaceConstructionData<dim>::OuterElementFace() const
 {
    return nbor_faces_.second;
 }


/// the material out of which the element consists from which the boundary face shall be constructed
template<uint32_t dim>
long FaceConstructionData<dim>::ElementMaterial() const
 {
     return material_;
 }



/// in the order: 1) parent, 2) inner, 3) outer
template<uint32_t dim>
void FaceConstructionData<dim>::Materials( vector<long>& mtrls ) const
 {
     mtrls.resize(3U);
     mtrls[0] = material_;
     mtrls[1] = materials_.first;
     mtrls[2] = materials_.second;
 }
 
 
 /// those of the neighboring higher-dimensional elements
template<uint32_t dim>
pair<long,long> FaceConstructionData<dim>::Materials() const
 {
     return materials_;
 }

 
template<uint32_t dim>
void FaceConstructionData<dim>::PatchNumber( uint32_t number )
 {
     patch_number_ = number;
 }


 
template<uint32_t dim>
uint32_t FaceConstructionData<dim>::PatchNumber() const
 {
     return patch_number_;
 }



template<uint32_t dim>
void FaceConstructionData<dim>::Out() const
 {
    // parent_element_
    cout <<"\nFaceConstructionData::Out: index of parent element: " <<" belonging to patch: "<< patch_number_;
    cout <<"\n\tindices of higher-dimensional elements on inside (first) and outside (second): ";
    cout << neighbors_.first <<" and "<< neighbors_.second;
    cout <<"\n\tmatching faces on inside (first) and outside (second):                         ";
    cout << nbor_faces_.first <<" and "<< nbor_faces_.second;
    cout <<"\n\tinteger codified juxtaposed regions: "<< materials_.first <<" and "<< materials_.second;
    cout <<" as well as parent region: "<< material_;
    cout << endl;
 }

template class FaceConstructionData<3U>;
template class FaceConstructionData<2U>;
template class FaceConstructionData<1U>;




/**
     higherDimensionalNeighbors() is used by CreateInternalBoundaryFrom( dim-minus1-region...)  and for a similar method which creates a SplitBoundary.

     higherDimensionalNeighbors() - finds the higher-dim neighbor elements of
     a dim-1 element embedded within the higher-dim mesh.
     
     - finds the IDs of the higher dimensional neighbors of the current lower-dimensional element
 
     - identifies which of the neighbors is on the inside and which on the outside
       as indicated by the normal direction of the lower dimensional element
 
     - identifies the materials on either side
     
     @return two inside-outside pairs of element idx numbers and corresponding materials on either side
     all the data are stored in the returnd FaceConstructionData object.
     
     assumptions
     - assumes that the nodes and elements in the entire model domain are numbered continuously
     
     application
     - use this function for finding neighbors of a surface element that sits on the inside of another region
 
      @attention both neighbors have to be present for this to work.
      
     TODO: if speed becomes critical, perhaps refactor using set_intersection() on nodes of elements
 
     @author SKM
     @date 2016
  
*/
template<uint32_t dim>
FaceConstructionData<dim>  higherDimensionalNeighbors( Element<dim>& e, const csmp::Index& mtrl_key )
 {
     if constexpr ( dim == 2U ) assert( e.IsLine() );
     if constexpr ( dim == 3U ) assert( e.IsSurface() );
     assert( mtrl_key.place != UNDEFINED );

     // 1. looping over the parent elements of the nodes searching for the faces which are shared with the lower dimensional element
     // ----------------------------------------------------------------------------------------------------------------------------
     auto nbors = parentElements<dim>( e.CornerNodesBegin(), e.CornerNodesEnd() );
     pair<uint32_t,uint32_t> faces{ nbors.first.second, nbors.second.second };
     pair<long,long> materials{ (*nbors.first.first).Read( mtrl_key ), (*nbors.second.first).Read( mtrl_key ) };

#ifndef NDEBUG
    assert( faces.first  != numeric_limits<uint32_t>::max() );
    assert( faces.second != numeric_limits<uint32_t>::max() );

    // VERIFICATION
    // ------------
    // verifying the inside neighbor element by projecting the face normal onto lower-dim-element normal

    Point<dim> enrml = e.UnitNormal();

    // first element
    Point<dim> fnrml = (*nbors.first.first).UnitNormalToFace( faces.first );
    // if the projection is negative, the first element lies on the outside
    enum POSITION { INNER_ELMT, OUTER_ELMT };
    POSITION  epos_elmt1 = ( dotProduct( enrml, fnrml ) < 0. ) ? OUTER_ELMT : INNER_ELMT;

    // swapping sides if necessary
    if ( epos_elmt1 != INNER_ELMT ) {
         cerr <<"\nhigherDimensionalNeighbors: normal inconsistency found and fixed; element:"<< e.Idx();
         cerr <<" "<< parseFiniteElementType( e.FE_Type() ) << endl;
         swap( nbors.first.first, nbors.second.first );
         swap( faces.first, faces.second );
         swap( materials.first, materials.second );
      }
#endif

    // initialise with nbors, their faces, adjacent materials, and patch numbers
    auto neighbors = make_pair( nbors.first.first, nbors.second.first );
    if ( mtrl_key.place != UNDEFINED )
      return FaceConstructionData( e, neighbors, faces, materials,
                                   static_cast<long>(e.Read( mtrl_key)) );
    else
      return FaceConstructionData( e, neighbors, faces, materials, UNSPECIFIED );
      
 } // end higherDimensionalNeighbors

/* DEPRECATED CODE

     // making a set of element nodes to serve as a key for later identification of faces
     set<Node<dim>*> node_set;
     const uint32_t nodes(e.Nodes());
     for ( uint32_t i{0U}; i<nodes; ++i ) node_set.insert( e.N(i) );
     // neighbor elements and their faces
     map<Element<dim>*,uint32_t>  nbor_elmts;
     for ( uint32_t i{0U}; i<nodes; i++ ) {
          const uint32_t parents{ e.N(i)->Parents() };
          for ( uint32_t j{0U}; j<parents; ++j ) {
               assert( e.N(i)->Parent(j) );
               const uint32_t faces{ e.N(i)->Parent(j)->Faces() };
               for ( uint32_t k{0U}; k<faces; ++k ) {
                     set<Node<dim>*> test_set = e.N(i)->Parent(j)->CornerNodesOfFace(k);
                     // if the face is shared with the element its face and face number are recorded
                     if ( node_set == test_set ) {
                          // storing a pointer to this element and its local face number
                          // making sure that no duplicate is received
                          nbor_elmts.insert( make_pair( e.N(i)->Parent(j), k ) );
                       }
                     test_set.clear();
                 }
            }
       }
*/


template FaceConstructionData<3U>  higherDimensionalNeighbors( Element<3U>&, const csmp::Index& );
template FaceConstructionData<2U>  higherDimensionalNeighbors( Element<2U>&, const csmp::Index& );
template FaceConstructionData<1U>  higherDimensionalNeighbors( Element<1U>&, const csmp::Index& );

} // end csmp
