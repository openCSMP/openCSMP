// Copyright © 2016 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FACE_CONSTRUCTION_DATA_H
#define FACE_CONSTRUCTION_DATA_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Element;

/**
    @brief Container of pairs of pointers to Element objects juxtaposed along a shared face identified by its local face number.
    Used in the construction of Boundary and SplitBoundary objects.
    
    The juxtaposed elements must have the same dimension as the model.
    The first element in the pair is the inside one and the second one the outside one.
    Material values identify which regions they belong to.
    There is also a Face or InterFace patch number used to identify the new modelsubdomain that will be created with the help
    of FaceConstructionData.
    
    @attention FaceConstructionData is populated by the function
    
       higherDimensionalNeighbors( csmp::Element<dim>&, const csmp::Index& );
*/
template<uint32_t dim>
class FaceConstructionData {
  public:
    FaceConstructionData( Element<dim>& parent_element,
                          std::pair<Element<dim>*,Element<dim>*>& neighbor_elements,
                          std::pair<uint32_t,uint32_t>& neighbor_element_faces,
                          std::pair<long,long>& neighbor_element_materials,
                          long   this_material );

    FaceConstructionData( Element<dim>& parent_element,
                          std::pair<Element<dim>*,Element<dim>*>& neighbor_elements,
                          std::pair<uint32_t,uint32_t>& neighbor_element_faces,
                          std::pair<long,long>& neighbor_element_materials,
                          long   this_material,
                          uint32_t number_of_patch_this_data_belongs_to );
    
     // NO ASSIGMENT OPERATOR BECAUSE ALL CLASS MEMBERS ARE CONSTANT
  
    /// inner (first) and out (second) higher-order parent element
    std::pair<Element<dim>*,Element<dim>*> NeighborElements() const;
  
    /// idx of lower-dimensional element from which face data were derived if it exists
    Element<dim>* LowerDimElement();
  
    /// idx of element on the opposite site of the outward pointing normal
    Element<dim>* const InnerElement() const;
    
    /// the local id of the face located at the boundary
    uint32_t InnerElementFace() const;
  
    /// idx of element on the side to which the normal points to
    Element<dim>* const OuterElement() const;
    uint32_t OuterElementFace() const;
  
    /// the material out of which the element consists from which the boundary face shall be constructed
    long ElementMaterial() const;
  
    /// in the order: 1) parent, 2) inner, 3) outer
    void Materials( std::vector<long>& ) const;
  
    /// Returns the materials of the neighboring higher-dimensional elements, inside followed by outside
    std::pair<long,long> Materials() const;
  
    /// identifier of the boundary that the face will belong to
    void     PatchNumber( uint32_t number );
    uint32_t PatchNumber() const;
  
    void Out() const;
  
  private:
    FaceConstructionData();
  
  private:
    Element<dim>&                                dim_m1_element_; ///< idx of lower-dimensional parent element
    const std::pair<Element<dim>*,Element<dim>*> neighbors_;    ///< indices of higher-dimensional elements on inside (first) and outside (second)
    const std::pair<uint32_t,uint32_t>           nbor_faces_;   ///< matching faces on inside (first) and outside (second)
    const std::pair<long,long>                   materials_;    ///< integer codified juxtaposed regions
    const long                                   material_;     ///< of the element from which the Face shall be constructed
    uint32_t                                     patch_number_; ///< unique identifier for the internal boundary patch this data relates to
};


/// finds the neighbors of dim-1 element, and their faces that connect to it; index records neighbor materials
template<uint32_t dim>
FaceConstructionData<dim>  higherDimensionalNeighbors( csmp::Element<dim>&, const csmp::Index& );

} // end csmp

#endif /* FACE_CONSTRUCTION_DATA_H */
