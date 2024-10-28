//
//  vsetMakers_Test.hpp
//  Open CSMP++SAMG
//
//  Created by Stephan Matthai on 29/5/2024.
//  Copyright © 2024 The University of Melbourne. All rights reserved.
//

#ifndef CSMP_VSET_MAKERS_TEST_H
#define CSMP_VSET_MAKERS_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class VSet;

/**
    Tests the suite of  meshes/boundary flags, data etc. that serve as the input for many of the unit tests created for the CSMP library.
    
    TODO: test prism_test for 'invalid' elements using geometric calculations
    TODO: test whether the faces are in the expected place
    
*/
class vsetMakers_Test : public Test {
  public:
    virtual void run();

    /// test numbering/faces/neighbors/segments compliance of single element meshes with element numbering etc. for CSMP_FEM_conventions.pdf
    bool TestConsistencyWithCSMP_Conventions( const VSet<2U>& );
    bool TestConsistencyWithCSMP_Conventions( const VSet<3U>& );
  
  private: // auxiliary functions
  
    /// returns node-sets (node idx using global numbering) of element faces, in the form of a vector
    template<uint32_t dim>
    std::vector<std::vector<uint32_t>> NodesOfElementFaces( const VSet<dim>& vset, size_t elmt );
    
    /// as determined using the nodes that make up the faces; tests also whether nodes are correctly assigned to faces
    bool FaceNumbersMatch_CSMP_Conventions();
    bool doNodeNumbersMatchCSMP_Conventions();

    // TODO: to test for consistency with 2024 CSMP conventions
    /*
    /// quadratic quadrilateral single-element mesh
    VSet<2U>  test_CreateVSet();
    
    void create_TrianglePatch_VSet( VSet<2U>& );
    
    void create_One_Square_VSet(VSet<2U>&, double length_of_sides, bool bSkewed=false );
    
    void create_Square_VSet( VSet<2U>&, size_t size_sides, double dimension, bool skewed=false );
    
    void create_SlitRectangle_VSet(VSet<2U>&, size_t x_dimension, size_t y_dimension, double x_length, 
                                        double y_length, size_t depth_of_slit, bool bSkewed=false );
 
    */
    
    // 2D poly-element sets and topology info
    // --------------------------------------
    /*
    ///  Rectangle-shaped MODEL_TINY, consisting of 1 line element two triangles, 1 quadrilateral and 6 face object marking the box boundary.
    ModelTopology  test_CreateSimplestPolyElement2DModel( VSet<2U>& );

    /// Rectangle shaped mixed model with 2 intersecting line element regions (USED IN UNIT TESTS)
    ModelTopology  create_MeshPatchWithLineElements_VSet( VSet<2U>& );

    /// model SPLIT22_BASIC with box boundaries (Faces) and one through-going and one internal crossing split boundary (USED IN UNIT TESTS)
    ModelTopology  create_BoundarySplitBoundaryPatch( VSet<2U>& );
    */
    
    // 3D meshes of a single element type
    // ----------------------------------
    /*
    /// 6 tets from a cube, no midside nodes
    void testCreateTetra_VSet( VSet<3U>& );
    
    void create_Pyramid_VSet( VSet<3U>&, bool bSkewed=false );
    
    void create_One_Hexahedra_VSet( VSet<3U>&, bool bSkewed = false );
    
    void create_Hexahedra_VSet( VSet<3U>&, bool bSkewed=false );
    
    void create_One_Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    void create_Prism_VSet( VSet<3U>&, bool bSkewed=false );
    */
    
    // 3D element sets with different types of elements
    // ------------------------------------------------
    /* USED IN 3D UNIT TESTS
    void create_Pyramid_Hexa_VSet( VSet<3U>&, bool bSkewed=false );

    void create_Prism_Hexa_VSet( VSet<3U>&, bool bSkewed=false );
    
    // more complex models created externally but available as VSets
    // -------------------------------------------------------------
    // USED IN UNIT TEST
    void create_FracBox( ModelTopology&, VSet<3U>& );
    */
};


} // end csmp

#endif /* CSMP_VSET_MAKERS_TEST_H */
