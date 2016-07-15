//
//  FaceConstructionData.cpp
//
//  Created by Stephan Matthai on 28/03/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include "FaceConstructionData.h"

using namespace std;

namespace csmp {

FaceConstructionData::FaceConstructionData( size_t  parent_element,
                                            std::pair<size_t,size_t>& neighbors,
                                            std::pair<size_t,size_t>& nbor_faces,
                                            std::pair<long,long>&     materials,
                                            long material )
 : parent_element_(parent_element),
   neighbors_(neighbors),
   nbor_faces_(nbor_faces),
   materials_(materials),
   material_(material),
   patch_number_(numeric_limits<size_t>::max())
 {
 }


FaceConstructionData::FaceConstructionData( size_t  parent_element,
                                            std::pair<size_t,size_t>& neighbors,
                                            std::pair<size_t,size_t>& nbor_faces,
                                            std::pair<long,long>&     materials,
                                            long material,
                                            size_t patch_number )
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
 
 
size_t FaceConstructionData::InnerElementFace() const
 {
    return nbor_faces_.first;
 }
 

/// on the side to which the normal points to
size_t FaceConstructionData::OuterElement() const
 {
    return neighbors_.second;
 }
 
 
size_t FaceConstructionData::OuterElementFace() const
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

 
void FaceConstructionData::PatchNumber( size_t number )
 {
     patch_number_ = number;
 }


 
size_t FaceConstructionData::PatchNumber() const
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
 }



} // end csmp
