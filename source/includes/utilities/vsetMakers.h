#ifndef VSET_MAKERS_H
#define VSET_MAKERS_H

#include "VSet.h"
#include "ModelTopology.h"

namespace csmp {

    /// all methods create VSets numbered 0..n-1
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
    ModelTopology  test_Create_MeshPatchWithLineElements_VSet( VSet<2U>& );


    // 3D single element sets
    // ----------------------
    void test_Create_Pyramid_VSet( VSet<3U>&, bool bSkewed=false );
    
    void test_Create_One_Hexahedra_VSet( VSet<3U>&, bool bSkewed = false );
    
    void test_Create_Hexahedra_VSet( VSet<3U>&, bool bSkewed=false );
    
    void test_Create_One_Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    void test_Create_Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    // 3D element sets with different types of elements
    // ------------------------------------------------
    void test_Create_Pyramid_Hexa_VSet( VSet<3U>&, bool bSkewed=false );

    void test_Create_Prism_Hexa_VSet( VSet<3U>&, bool bSkewed=false );
}

#endif
