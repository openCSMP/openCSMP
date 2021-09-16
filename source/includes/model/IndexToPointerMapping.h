//
//  MeshPointers.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 16/5/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_INDEX_TO_POINTER_MAPPING_H
#define CSMP_INDEX_TO_POINTER_MAPPING_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;


/**
     IndexToPointerMapping  - is used to build the dynamic tree structure that connects elements with nodes and defines regions
     when a model is build from a VSet.

     This is dum data struct, keeping a record of the index-pointer relationship from when the model was created.
     Class is used for the swift construnction of regions, boundaries and splitboundaries without a need to traverse the mesh.
     
     @author SKM
     @date 16/5/2021
*/
template<size_t dim>
struct IndexToPointerMapping {
  std::vector<Node<dim>*>      nodeConnector_;
  std::vector<Element<dim>*>   elementConnector_;
  std::vector<Face<dim>*>      faceConnector_;
  std::vector<InterFace<dim>*> interFaceConnector_;
  bool WithFaces() const { return !faceConnector_.empty(); }
  bool WithInterFaces() const { return !interFaceConnector_.empty(); }
  void Clear();
};

} // end csmp

#endif /* CSMP_INDEX_TO_POINTER_MAPPING_H */
