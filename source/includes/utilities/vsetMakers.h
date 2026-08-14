#ifndef VSET_MAKERS_H
#define VSET_MAKERS_H

#include "VSet.h"
#include "ModelTopology.h"

/**
# CSMP++ Mesh Conventions Summary
see also doc/CSMP_FEM_conventions.pdf

CSMP Mesh and Boundary Conventions
==================================

In CSMP, mesh nodes and elements are numbered 0..n-1 following C's array indexing convention. size_t and int64_t are used to represent these numbers.

CSMP supports both, simplex type elements (linear bar, linear triangle and linear tetrahedron), and general element types like the quadrilateral, hexahedron, prism and pyramid. 

For the simplex-type elements the faces sit opposite to their nodes. Thus, the faces of the triangle have the nodes (numbered from the outside looking in), face0 1,2, face1 2,0, face2 0,1. For the tetrahedron, the face nodes (from the outside looking in) are, for face0 ={1,2,3}, face1 ={0,3,2}, face2 ={0,1,3} and for face3 ={0,2,1}. In quadratic elements the mid-side node numbers have higher numbers than the highest corner node. Thus, for the triangle, node 3 is on the 0-1 edge, node 4 on the 1-2 edge and node 5 on the 2-0 edge. For the tetrahedron, mid-side nodes 4,5,6,7,8,9 are on the edges 0-1, 1-2, 2-3, 3-0, 3-1, 3-2, 3-0, respectively.

For the linear hexahedron, face0 with nodes {0,3,2,1} lies in the base plane (labelled bottom), face1 has the nodes {0,1,5,4} in the xy-plane at the front. Face2 with the nodes {1,2,6,5} lies in the yz-plane on the right, Face3 has the nodes {2,3,7,6} and lies in the xy-plane in the back. Face4 with the nodes {0,4,7,3} lies in the yz-plane on the left, and face5 with the nodes {4,5,6,7} lies in the xz-plane at the top.

The linear pyramid, has the same base plane as the hexahedron (nodes 0,1,2,3) but is is face4. Face0 with the nodes {0,1,4} faces the front; face1 with the nodes {1,2,4} faces the right; face2 with the nodes {2,3,4} faces the back. Face3 with nodes {0,4,3}, is on the left.

For the linear prism element, the faces have the nodes: face0 {0,2,1}, face1 {0,1,4,3}, face2 {1,2,5,4}, face3 {0,3,5,2}, and face4  {3,4,5}. Face0 is at the base (bottom), and face4 at the top.

For all elements the neighbor numbering corresponds to the face numbering. The term neighbor refers to the element that sits within the mesh on the opposite side of a Face. The global number of this adjacent element (which can be a volume, surface or line) is stored in the mesh container called Vdata inside of VSet for this face. It is stored in a deque called 'pfverts<int64_t>'. However, elements at the boundaries of the mesh are missing some neighbors. To distinguish this situation, a negative integer value is stored in pfverts, encoded by a boundary flag. 

There are boundary flags for sides, edges, and corners of box-shaped models as well as for internal boundaries, and irregulary shaped boundaries.
In this context, CSMP's right-hand-rule based coordinate system with x pointing to the right (E), y to the top, and z to the front (S) becomes important as it is used for the naming convention of the boundaries. Like the faces of the hexahedron the boundaries are ordered and part of the named enum BOX_BOUNDARY with the values as follows:

BOTTOM: y=min_y (xz-plane), lowest y-coordinate, base of the domain.
FRONT: z=max_z (xy-plane), front face toward viewer.
RIGHT: x=max_x (yz-plane), right side of domain.
BACK: z=min_z (xy-plane), back face away from viewer.
LEFT: x=min_x (yz-plane), left side of domain.
TOP: y=max_y (xz-plane), top of the domain, highest y-coordinate.
IRREGULAR: for arbitrarily shaped boundaries, and INTERNAL for boundaries inside of a model.

Each box-shaped model has 12 edges (parallel to the axes) are:

EDGE1: x-axis, y_min, z_min (e.g., (1,0,0), (2,0,0))
EDGE2: y-axis, x_max, z_min
EDGE3: x-axis, y_max, z_min
EDGE4: y-axis, x_min, z_min
EDGE5: z-axis, x_min, y_min
EDGE6: z-axis, x_max, y_min
EDGE7: z-axis, x_max, y_max
EDGE8: z-axis, x_min, y_max
EDGE9: x-axis, y_min, z_max
EDGE10: y-axis, x_max, z_max
EDGE11: x-axis, y_max, z_max
EDGE12: y-axis, x_min, z_max

These are also contained in BOX_BOUNDARY and get assigned to the nodes of the model. From their localcombination, the position of a face on the side of the model or on an edge can be uniquely identified. The flag 'NOT' is assigned to all nodes that are inside of a model.

The corners of the model mark the endpoints of its external edges:

CNR1: (x_min,y_min,z_min), CNR2: (x_max,y_min,z_min), CNR3: (x_max,y_max,z_min), CNR4: (x_min,y_max,z_min), CNR5: (x_min,y_min,z_max), CNR6: (x_max,y_min,z_max), CNR7: (x_max,y_max,z_max), CNR8: (x_min,y_max,z_max)

Thus, corners 1-4 are in the xy plane at z=0 (back). Corners 5-8 also lie in the xy plane but are in the front of the model at z_max.

Negative integer values distinguish boundary nodes, edges or faces. They are defined (in Box.h) as constexpr int8_t constants with names such as 

BACK_OUTSIDE: (xy-plane)
FRONT_OUTSIDE: (xy-plane)
LEFT_OUTSIDE: (yz-plane)
RIGHT_OUTSIDE: (yz-plane)
BOTTOM_OUTSIDE: (xz-plane)
TOP_OUTSIDE: (xz-plane at y_max)
etc.

By definition, an element has (at most) one equidimensional neighbor per face.
This avoids that line or surface elements in a volumetric mesh form manifolds.
Instead, only one specific face neighbor is chosen. For surface elements, this usually is the one that is closest to the plane of the neighbor on the opposite side of this face. 
As a consequence of these conventions a mesh that consists of volume, surface and line elements will not only have external boundaries, but each series of line elements within it, will have two line elements without a neighbor, and each surface patch within the model volume will have perimeter faces without neighbors. Line elements in the mesh must be oriented such that, in a chain (polyline), the last node matches the first node of the next adjacent element. A special situation arises where line elements surrounding a 2D model form a loop such that they all have neighbors. Similarly, surface elements covering the outside of a box- or arbitrarily shaped model, will not be missing any neighbors.

In addition to this neighbor connectivity, CSMP implements special types of lower-dimensional elements called Face and InterFace which do not only have equidimensional neighbors, but also know their higher-dimensional neighbors.

InterFace elements can only occur inside of models where they are used to connect the boundaries of node-matched but not node-sharing mesh patches. This type of boundary is called SplitBoundary in CSMP. InterFace and Face (node-sharing) elements forming INTERNAL boundaries are oriented such that the normals of them point in a uniform direction marking the inside and outside of the boundary. Face and InterFace thus have access to an InnerParent and OuterParent higher-dimensional element. 

When Face elements form the outer (external) boundaries of the model, they only have their INSIDE higher-dimensional neighbor. This is achieved by orienting them such that the Face normals point to the OUTSIDE.

InterFace elements can also contain lower-dimensional Element objects. In this case, there will be NodeManifold objects that relate three or more nodes to one-another.  

Element, Face and InterFace objects forming a mesh are collected into ModelSubDomain objects. There are three types of these: Element objects form Regions, Face objects form Boundaries and InterFace objects form SplitBoundaries. Each of these consist of cells and Nodes. InterFaces just store NodeManifold objects. Each of these can be unique (space exclusive) or non-unique (overlapping). They can be contiguous (flood-fill works) or discontiguous mesh patches. Each patch has a perimeter and interior. Elements, Faces or InterFaces are considered as perimeter located if they have at least one face/edge o the ModelSubDomain boundary (not necessarily the overall boundary of the mesh). All the elements of a model are grouped into a 'master' Region called "Model".

Polygonal meshes are stored in the VSet : public VData container. The VData base class holds the node coordinates (px,py,pz), cell types ('pelmt' populated with the likes of IsoparametricLinearTriangle etc.), the nodes that make up each element ('plist'), the cell neighbors (pfverts'), node boundary flags ('pbflags'), material IDs of the Elements ('pmtrl'), and some topology flags ('gflags'). In 'plist' and 'pfverts' elements are followed by faces (if any) and interfaces (if any). For Face and InterFace objects, 'pfverts' also contains the indices of the higher-dimensional neighbor elements, following the equidimensional neighbors.

@attention The orientation of the non-Simplex elements (Hexahedron, Prism, Pyramid) relative to the coordinate system is not uniquely defined by these conventions. Thus, many are possible, but 1 of these will give an element orientation where the bottom face is in the xz plane and so forth.

@note The geometry flags 'gflags' that distinguish nodes that define the topology of the model from ones that are mere degrees of freedom in the computational model, are provided only for those VSets, that contain internal boundaries.
 
*/
namespace csmp {

    /// mesh of a 4 linear quadrilateral elements
    VSet<2U>  create_Quadrilateral_VSet();
    
    // 2D single element-type meshes
    // -----------------------------
    
    /// square meshed with 2 triangles, cutting in diagonally downwards in the x direction
    void create_2_Triangle_VSet(VSet<2U>& vset );
    
    /// mesh patch with 10 triangles
    void create_TrianglePatch_VSet( VSet<2U>& );
    
    void create_1Square_VSet(VSet<2U>&, double length_of_sides, bool bSkewed=false );
    
    /// creates a  mesh of rectangles the dimensions of which are defined by the user. The model has a zero width slit in the middle where nodes are duplicated
    void create_SlitRectangle_VSet( VSet<2U>&, size_t x_dimension, size_t y_dimension, double x_length,
                                    double y_length, size_t depth_of_slit, bool bSkewed=false );
                                        
    /// creates either 1 square or 1 split rectangle mesh using the methods above
    void create_Square_VSet( VSet<2U>&, size_t size_sides, double dimension, bool skewed=false );


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
    

    // 3D single-elememt meshes
    // ------------------------
    
    /// single 8-noded hexahedron
    void create_1Hexahedron_VSet( VSet<3U>&, bool bSkewed = false );
    
    /// single 6-noded prism element
    void create_1Prism_VSet( VSet<3U>&, bool bSkewed=false );
    
    
    // 3D poly-elememt meshes
    // ----------------------
    
    /// 6 four-noded tetrahedra filling a cube, no midside nodes
    void create_Tetra_VSet( VSet<3U>& );
    
    /// 27 pyramids packed into cube
    void create_Pyramid_VSet( VSet<3U>&, bool bSkewed=false );
    
    /// 27 hexahedra forming a  cube
    void create_Hexahedra_VSet( VSet<3U>&, bool bSkewed=false );
    
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
    ModelTopology create_FracBox( VSet<3U>& );
}

#endif
