//
//  FaceConstructionData.hpp
//
//  Created by Stephan Matthai on 28/03/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#ifndef FACE_CONSTRUCTION_DATA_H
#define FACE_CONSTRUCTION_DATA_H

#include <iostream>
#include <vector>
#include <limits>

namespace csmp {

/// all it takes to build a face later
class FaceConstructionData {
  public:
    FaceConstructionData( size_t  parent_element,
                          std::pair<size_t,size_t>& neighbor_elements,
                          std::pair<uint32_t,uint32_t>& neighbor_element_faces,
                          std::pair<long,long>& neighbor_element_materials,
                          long   this_material );

    FaceConstructionData( size_t  parent_element,
                          std::pair<size_t,size_t>& neighbor_elements,
                          std::pair<uint32_t,uint32_t>& neighbor_element_faces,
                          std::pair<long,long>& neighbor_element_materials,
                          long   this_material,
                          uint32_t number_of_patch_this_data_belongs_to );
  
    FaceConstructionData( const FaceConstructionData& );
     // NO ASSIGMENT OPERATOR BECAUSE ALL CLASS MEMBERS ARE CONSTANT
  
    /// inner (first) and out (second) higher-order parent element
    std::pair<size_t,size_t> NeighborElements() const;
  
    /// idx of lower-dimensional element from which face data were derived if it exists
    size_t Element() const;
  
    /// idx of element on the opposite site of the outward pointing normal
    size_t   InnerElement() const;
    /// the local id of the face located at the boundary
    uint32_t InnerElementFace() const;
  
    /// idx of element on the side to which the normal points to
    size_t   OuterElement() const;
    uint32_t OuterElementFace() const;
  
    /// the material out of which the element consists from which the boundary face shall be constructed
    long ElementMaterial() const;
  
    /// in the order: 1) parent, 2) inner, 3) outer
    void Materials( std::vector<long>& ) const;
  
    /// those of the neighboring higher-dimensional elements
    std::pair<long,long> Materials() const;
  
    /// identifier of the boundary that the face will belong to
    void     PatchNumber( uint32_t number );
    uint32_t PatchNumber() const;
  
    void Out() const;
  
  private:
    FaceConstructionData();
    FaceConstructionData& operator=( const FaceConstructionData& );
  
  private:
    const size_t  parent_element_;                  ///< idx of lower-dimensional parent element
    const std::pair<size_t,size_t>     neighbors_;  ///< indices of higher-dimensional elements on inside (first) and outside (second)
    const std::pair<uint32_t,uint32_t> nbor_faces_; ///< matching faces on inside (first) and outside (second)
    const std::pair<long,long>         materials_;  ///< integer codified juxtaposed regions
    const long                         material_;   ///< of the element from which the Face shall be constructed
    uint32_t patch_number_;                           ///< unique identifier for the internal boundary patch this data relates to
};


} // end csmp

#endif /* FACE_CONSTRUCTION_DATA_H */
