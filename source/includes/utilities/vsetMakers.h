#ifndef VSET_MAKERS_H
#define VSET_MAKERS_H

#include "VSet.h"
#include "ModelTopology.h"

namespace csmp {

    /// quadratic quadrilateral single-element mesh
    VSet<2U>  test_CreateVSet();
    
    // 2D single element sets
    // ----------------------
    void test_Create_TrianglePatch_VSet( VSet<2U>& );
    
    void test_Create_One_Square_VSet(VSet<2U>&, double length_of_sides, bool bSkewed=false );
    
    void test_Create_Square_VSet( VSet<2U>&, size_t size_sides, double dimension, bool skewed=false );
    
    void test_Create_SlitRectangle_VSet(VSet<2U>&, size_t x_dimension, size_t y_dimension, double x_length, 
                                        double y_length, size_t depth_of_slit, bool bSkewed=false );
                                        
    // 2D poly-element sets and topology info
    // --------------------------------------
    
    ///  Rectangle-shaped MODEL_TINY, consisting of 1 line element two triangles, 1 quadrilateral and 6 face object marking the box boundary.
    ModelTopology  test_CreateSimplestPolyElement2DModel( VSet<2U>& );

    /// Rectangle shaped mixed model with 2 intersecting line element regions
    ModelTopology  test_Create_MeshPatchWithLineElements_VSet( VSet<2U>& );

    /// model SPLIT22_BASIC with box boundaries (Faces) and one through-going and one internal crossing split boundary
    ModelTopology  test_Create_BoundarySplitBoundaryPatch( VSet<2U>& );

    // 3D single element sets
    // ----------------------
    
    /// 6 tets from a cube, no midside nodes
    void testCreateTetra_VSet( VSet<3U>& );
    
    void test_Create_Pyramid_VSet( VSet<3U>&, bool bSkewed=false );
    
    void test_Create_One_Hexahedra_VSet( VSet<3U>&, bool bSkewed = false );
    
    void test_Create_Hexahedra_VSet( VSet<3U>&, bool bSkewed=false );
    
    void test_Create_One_Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    void test_Create_Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    // 3D element sets with different types of elements
    // ------------------------------------------------
    void test_Create_Pyramid_Hexa_VSet( VSet<3U>&, bool bSkewed=false );

    void test_Create_Prism_Hexa_VSet( VSet<3U>&, bool bSkewed=false );
    
    // more complex models created externally but available as VSets
    // -------------------------------------------------------------
    void test_Create_FracBox( ModelTopology&, VSet<3U>& );
}

#endif
