//
//  MeshPatchAttributes.h
//  CSMP_unit_tests
//
//  Created by Stephan Matthai on 23/5/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#ifndef CSMP_MESH_PATCH_ATTRIBUTES_H
#define CSMP_MESH_PATCH_ATTRIBUTES_H

#include "FiniteElement.h"

namespace csmp {

/**
      Information on contiguous cell-patch representing a submesh of the model.
      Used by MeshManager.
      Size info is used also to keep number information uptodate.
          
         @author SKM
         @date 23/5/21
*/
class MeshPatchAttributes {
  public:
    MeshPatchAttributes( size_t cells, CELL_SHAPE dim );
    
    void       Cells( size_t );
    size_t     Cells() const;
    
    /// whether this is a LINE, SURFACE or VOLUME mesh
    void       Geometry( CELL_SHAPE );
    CELL_SHAPE Geometry() const;
    
  private:
    size_t      cells_;          ///< number of cells in the patch
    CELL_SHAPE  cell_dimension_; ///<  types of elements in the patch
};

} // end csmp

#endif /* CSMP_MESH_PATCH_ATTRIBUTES_H */
