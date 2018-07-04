#ifndef VSET_MAKERS_H
#define VSET_MAKERS_H

#include "Model.h"

namespace csmp {

    /// all methods create VSets numbered 0..n-1
    VSet<2U>  test_CreateVSet();
    
    // 2D single element sets
    // ----------------------
    void test_Create_TrianglePatch_VSet( VSet<2U>& vset );
    
    void test_Create_One_Square_VSet(VSet<2U>& vset, double64 length_of_sides, bool bSkewed=false );
    
    void test_Create_Square_VSet( VSet<2U>& vset, size_t size_sides, double64 dimension, bool skewed=false );
    
    void test_Create_SlitRectangle_VSet(VSet<2U>& vset, size_t x_dimension, size_t y_dimension, double64 x_length, 
                                        double64 y_length, size_t depth_of_slit, bool bSkewed=false );
                                        
    // 3D single element sets
    // ----------------------
    void test_Create_Pyramid_VSet(VSet<3U>& vset, bool bSkewed=false );
    
    void test_Create_One_Hexahedra_VSet(VSet<3U>& vset, bool bSkewed = false );
    
    void test_Create_Hexahedra_VSet( VSet<3U>& vset, bool bSkewed=false );
    
    void test_Create_One_Prism_VSet(VSet<3U>& vset, bool bSkewed=false );
    
    void test_Create_Prism_VSet(VSet<3U>& vset, bool bSkewed=false );
    
    // 3D element sets with different types of elements
    // ------------------------------------------------
    void test_Create_Pyramid_Hexa_VSet( VSet<3U>& vset, bool bSkewed=false );

    void test_Create_Prism_Hexa_VSet(VSet<3U>& vset, bool bSkewed=false );
}

#endif
