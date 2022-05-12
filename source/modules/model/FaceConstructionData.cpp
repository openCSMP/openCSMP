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

FaceConstructionData::FaceConstructionData( size_t  parent_element,
                                            std::pair<size_t,size_t>& neighbors,
                                            std::pair<uint32_t,uint32_t>& nbor_faces,
                                            std::pair<long,long>&     materials,
                                            long material )
 : parent_element_(parent_element),
   neighbors_(neighbors),
   nbor_faces_(nbor_faces),
   materials_(materials),
   material_(material),
   patch_number_(numeric_limits<uint32_t>::max())
 {
 }


FaceConstructionData::FaceConstructionData( size_t  parent_element,
                                            std::pair<size_t,size_t>& neighbors,
                                            std::pair<uint32_t,uint32_t>& nbor_faces,
                                            std::pair<long,long>&     materials,
                                            long material,
                                            uint32_t patch_number )
 : parent_element_(parent_element),
   neighbors_(neighbors),
   nbor_faces_(nbor_faces),
   materials_(materials),
   material_(material),
   patch_number_(patch_number)
 {
 }
 
 
  
FaceConstructionData::FaceConstructionData( const FaceConstructionData& fcd )
 : parent_element_(fcd.parent_element_),
   neighbors_(fcd.neighbors_),
   nbor_faces_(fcd.nbor_faces_),
   materials_(fcd.materials_),
   material_(fcd.material_),
   patch_number_(fcd.patch_number_)
 {
 }



/// inner (first) and out (second) higher-order parent element
pair<size_t,size_t> FaceConstructionData::NeighborElements() const
 {
     return neighbors_;
 }

/**
    lower-dimensional parent element from which face data were recovered in part
    
    @attention if this element does not exist, NULL_IDX is returned
*/
size_t FaceConstructionData::Element() const
 {
    return parent_element_;
 }

/// on the opposite site of the outward pointing normal
size_t FaceConstructionData::InnerElement() const
 {
    return neighbors_.first;
 }
 
 
uint32_t FaceConstructionData::InnerElementFace() const
 {
    return nbor_faces_.first;
 }
 

/// on the side to which the normal points to
size_t FaceConstructionData::OuterElement() const
 {
    return neighbors_.second;
 }
 
 
uint32_t FaceConstructionData::OuterElementFace() const
 {
    return nbor_faces_.second;
 }


/// the material out of which the element consists from which the boundary face shall be constructed
long FaceConstructionData::ElementMaterial() const
 {
     return material_;
 }



/// in the order: 1) parent, 2) inner, 3) outer
void FaceConstructionData::Materials( vector<long>& mtrls ) const
 {
     mtrls.resize(3U);
     mtrls[0] = material_;
     mtrls[1] = materials_.first;
     mtrls[2] = materials_.second;
 }
 
 
 /// those of the neighboring higher-dimensional elements
pair<long,long> FaceConstructionData::Materials() const
 {
     return materials_;
 }

 
void FaceConstructionData::PatchNumber( uint32_t number )
 {
     patch_number_ = number;
 }


 
uint32_t FaceConstructionData::PatchNumber() const
 {
     return patch_number_;
 }



void FaceConstructionData::Out() const
 {
    cout <<"\nFaceConstructionData::Out: index of parent element: "<< parent_element_ <<" belonging to patch: "<< patch_number_;
    cout <<"\n\tindices of higher-dimensional elements on inside (first) and outside (second): ";
    cout << neighbors_.first <<" and "<< neighbors_.second;
    cout <<"\n\tmatching faces on inside (first) and outside (second):                         ";
    cout << nbor_faces_.first <<" and "<< nbor_faces_.second;
    cout <<"\n\tinteger codified juxtaposed regions: "<< materials_.first <<" and "<< materials_.second;
    cout <<" as well as parent region: "<< material_;
    cout << endl;
 }



/**
     Method used by CreateInternalBoundaryFrom( dim_minus1_region...

     higherDimensionalNeighbors() - finds the higher-dim neighbor elements of
     a dim-1 element embedded within the higher-dim mesh.
     
     @attention both neighbors have to be present for this to work.
     
     = LOCAL METHOD ONLY KNOWN TO THIS COMPILATION UNIT
 
     - finds the IDs of the higher dimensional neghibors of the current element
 
     - identifies which of the neighbors is on the inside and which on the outside
       as indicated by the normal direction of the lower dimensional element
 
     - identifies the materials on either side
     
     @return two inside-outside pairs of element idx numbers and corresponding materials on either side
     all the data are stored in the returnd FaceConstructionData object.
     
     assumptions
     - assumes that the nodes and elements in the entire model domain are numbered continuously
     
     application
     - use this function for finding neighbors of a surface element that sits on the inside of another region
     
     @author SKM
     @date 2016
  
*/
template<uint32_t dim>
FaceConstructionData  higherDimensionalNeighbors( const Element<dim>& e, const csmp::Index& mtrl_key )
 {
     if constexpr ( dim == 2U ) assert( e.IsLine() );
     if constexpr ( dim == 3U ) assert( e.IsSurface() );
     assert( mtrl_key.place != UNDEFINED );

     // 1. looping over the parent elements of the nodes searching for the faces which are shared with the lower dimensional element
     // ----------------------------------------------------------------------------------------------------------------------------
     
     // making a set of element nodes to later identify faces by comparison
     set<size_t> node_set, test_set;
     const auto  nodes(e.Nodes());
     for ( auto i{0U}; i<nodes; ++i )   node_set.insert(e.N(i)->Idx());
     // neighbor elements and their faces
     map<const Element<dim>*,uint32_t>  nbor_elmts;
     vector<uint32_t> fnids;
     for ( auto i{0U}; i<nodes; i++ ) {
          const auto parents(e.N(i)->Parents());
          for ( auto j{0U}; j<parents; ++j ) {
               const Element<dim>* const eptr(e.N(i)->Parent(j));
               const auto faces(eptr->Faces());
               for ( auto k{0U}; k<faces; ++k ) {
                     eptr->FE()->NodesOfFace( k, fnids );
                     size_t fnodes(fnids.size());
                     for ( auto l{0U}; l<fnodes; ++l )
                       test_set.insert( eptr->N( fnids[l])->Idx() );
                     // if the face is shared the element and its face are recorded
                     if ( node_set == test_set ) {
                          // storing a pointer to this element and its local face number
                          // making sure that no duplicate is received
                          nbor_elmts.insert( make_pair(eptr,k) );
                       }
                     test_set.clear();
                 }
            }
       }

#ifdef DEBUG
    // VERIFICATION
    assert( nbor_elmts.size() == 2U );
    const auto e1{ nbor_elmts.begin() };
    const auto e2{ nbor_elmts.rbegin() };
    assert( (*e1).first->Neighbor( (*e1).second ) == (*e2).first );
    assert( (*e2).first->Neighbor( (*e2).second ) == (*e1).first );
#endif

    // 2. finding the inside neighbors by projecting face normals onto lower dim element normal
    // ----------------------------------------------------------------------------------------
    vector<double>  enrml, fnrml;
    e.UnitNormal( enrml );
    const auto nbor1{ nbor_elmts.begin() };
    const auto nbor2{ nbor_elmts.rbegin() };
    enum POSITION { INNER_ELMT, OUTER_ELMT };
    //                             inside first           outside second
    pair<size_t,size_t>     nbors{ (*nbor1).first->Idx(), (*nbor2).first->Idx() };
    pair<uint32_t,uint32_t> faces{ (*nbor1).second, (*nbor2).second };
    pair<long,long>         materials{ (*nbor1).first->Read( mtrl_key ), (*nbor2).first->Read( mtrl_key ) };

    // first element
    // -------------
    (*nbor1).first->UnitNormalToFace( faces.first, fnrml );
    double dotproduct(0.);
    for ( auto k{0U}; k<dim; ++k ) dotproduct += enrml[k] * fnrml[k];
   
    // if the projection is negative, the first element lies on the outside
    POSITION  epos_elmt1 = ( dotproduct < 0. ) ? OUTER_ELMT : INNER_ELMT;

    // second element
    // --------------
    (*nbor2).first->UnitNormalToFace( faces.second, fnrml );
    dotproduct = 0.;
    for ( auto k{0U}; k<dim; ++k ) dotproduct += enrml[k] * fnrml[k];

    POSITION  epos_elmt2 = ( dotproduct < 0. ) ? OUTER_ELMT : INNER_ELMT;

    // checking that we have no duplication here
    assert( epos_elmt1 != epos_elmt2 );

    // swapping sides if necessary
    if ( epos_elmt1 != INNER_ELMT ) {
         swap( nbors.first, nbors.second );
         swap( faces.first, faces.second );
         swap( materials.first, materials.second );
      }

    // initialise with nbors, their faces, adjacent materials, and patch numbers
    if ( mtrl_key.place != UNDEFINED )
      return FaceConstructionData( e.Idx(), nbors,  faces, materials,
                                   static_cast<long>(e.Read( mtrl_key)) );
    else
      return FaceConstructionData( e.Idx(), nbors, faces, materials, -1 );
      
 } // end higherDimensionalNeighbors

template FaceConstructionData  higherDimensionalNeighbors( const Element<3U>&, const csmp::Index& );
template FaceConstructionData  higherDimensionalNeighbors( const Element<2U>&, const csmp::Index& );
template FaceConstructionData  higherDimensionalNeighbors( const Element<1U>&, const csmp::Index& );



} // end csmp
