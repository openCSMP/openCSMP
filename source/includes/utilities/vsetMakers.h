#ifndef VSET_MAKERS_H
#define VSET_MAKERS_H

#include "VSet.h"
#include "ModelTopology.h"

namespace csmp {

    /// mesh of a 4 linear quadrilateral elements
    VSet<2U>  create_Quadrilateral_VSet();
    
    // 2D single element-type meshes
    // -----------------------------
    
    /// mesh patch with 10 triangles
    void create_TrianglePatch_VSet( VSet<2U>& );
    
    void create_1Square_VSet(VSet<2U>&, double length_of_sides, bool bSkewed=false );
    
    /// creates a  mesh of rectangles the dimensions of which are defined by the user. The model has a zero width slit in the middle where nodes are duplicated
    void create_SlitRectangle_VSet( VSet<2U>&, int x_dimension, int y_dimension, double x_length,
                                    double y_length, int depth_of_slit, bool bSkewed=false );
                                        
    /// creates either 1 square or 1 split rectangle mesh using the methods above
    void create_Square_VSet( VSet<2U>&, int size_sides, double dimension, bool skewed=false );

    // 2D poly-element type meshes wirh topology info
    // ----------------------------------------------
    ///  Rectangle-shaped MODEL_TINY, consisting of 1 line element two triangles, 1 quadrilateral and 6 face object marking the box boundary.
    ModelTopology  create_SimplePolyElement2DModel( VSet<2U>& );

    /// Rectangle shaped mixed model with 2 intersecting line element regions
    ModelTopology  create_MeshPatchWithLineElements_VSet( VSet<2U>& );

    /// model SPLIT22_BASIC with box boundaries (Faces) and one through-going and one internal crossing split boundary
    ModelTopology  create_BoundarySplitBoundaryPatch( VSet<2U>& );
    
    /// 2D rectangular model the two halfs of each are offset from one another
    void create_Disconnected2D_VSet( VSet<2U>& );
    

    // 3D single-element-type meshes
    // -----------------------------
    /// Rubik cube of 27 hexahedra, three XY planes of 3x3 elements, numbered from the left to the right (x=0..3), from the bottom to the top (y=0..3)
    void create_RubikCube( VSet<3U>& );

    /// 6 four-noded tetrahedra filling a cube, no midside nodes
    void create_Tetra_VSet( VSet<3U>& );
    
    /// 27 pyramids packed into cube
    void create_Pyramids_VSet( VSet<3U>&, bool bSkewed=false );
    
    /// single 8-noded hexahedron
    void create_1Hexahedron_VSet( VSet<3U>&, bool bSkewed = false );
    
    /// 27 hexahedra forming a  cube
    void create_Hexahedra_VSet( VSet<3U>&, bool bSkewed=false );
    
    /// single 6-noded prism element
    void create_1Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    /// mesh with 54 prism elements
    void create_Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    // 3D element sets with different types of elements
    // ------------------------------------------------
    void create_Pyramid_Hexa_VSet( VSet<3U>&, bool bSkewed=false );

    /// three layers of 3 x 3 hexahedral elements (in the XY plane, but with a Z-axis aligned stack of prism elements in the middle
    void create_Prism_Hexa_VSet( VSet<3U>&, bool bSkewed=false );
    
    // more complex models created externally but available as VSets
    // -------------------------------------------------------------
    /// simplex mesh with tetra, triangles and line elements
    void create_FracBox( ModelTopology&, VSet<3U>& );
}

#endif
