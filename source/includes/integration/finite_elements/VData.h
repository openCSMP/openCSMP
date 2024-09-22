#ifndef CSMP_VDATA_H
#define CSMP_VDATA_H

#include "CSMP_definitions.h"
#include "CSMP_global_enumerations.h"

namespace csmp {

/**

@brief Container for polygonal data implemented on the basis of NCSA 
(Urbana, Champagne, Illinois, USA) VSet
that forms part of their HDF data file storage.

@author Stephan Matthai
@author Stephen G. Roberts
@date 1996

@section motivation Motivation

Mesh connectivity usually termed as VData typically consists of lists
of (1) node and constraint-point coordinates (constraint points can be nodes
between nodes but also integration points), (2) a listing of the nodes
which belong to each element, and (3) a list of the neighbors of each
element. If the element lies on the model boundary and there is no neighbor
then entries in this list will be UNSPECIFIED. (4) A typical VData set 
also specifies which nodes lie at the model boundary and which property
values apply to these nodes. 

The class VData brings all these datasets together into an object 
which holds the information to build this V-type mesh connectivity.

@section applicability Applicability

VData can be used as isolated object but it is also integrated into
the VSet which also holds property data which are assigned to the 
mesh.
 
@section structure Structure
 
The VData object contains coordinate objects, deques of
vectors to store the connectivity and maps to hold the boundary flags
and boundary values.

Detailed content:

(double) px, py, pz[0..n-1]  = coordinates of n nodes (= FE vertices)

(int8_t) bflags[0..n-1] = negative numbers for nodes at boundary, see later

(enum int8_t) pelmt[0..e-1] = CSMP_FEM_TYPE of e finite elements  (1 entry only for single element-type mesh)

(size_t) plist[0..e *[sum npe[e]]] = nodes making up each element

(int64_t) pfverts[0..e *[sum fpe[e]]] = neighbor elements adjacent to the faces of each element

For Face and InterFace objects, 'pfverts' also contains the indices of the higher-dimensional neighbor element.

*/
class VData {
  public:
    VData();

    /// constructor for hybrid element meshes 
    VData( const std::deque<uint32_t>& npes, ///< just the sizes of the different vectors
           const std::deque<uint32_t>& epes,
           size_t nodes );

    /// constructor for meshes that only hold a single element type
    VData( uint32_t nodes_per_element, uint32_t nbors_per_element, size_t nodes, size_t elmts );
    
    /// if the mesh contains only a single element type it can be set with this method
    void SingleElementType( int8_t etype );
    
    /// CSMP_FEM_TYPEs for all the elements stored in the VData
    void ElementTypes( const std::vector<int8_t>& elmt_types );
    
    void AddElementTypes( std::vector<int8_t>::const_iterator first,
                          std::vector<int8_t>::const_iterator last );

    void AddElementTypes( std::deque<int8_t>::const_iterator first,
                          std::deque<int8_t>::const_iterator last );
     
    /// the storage for node manifolds inside of VData
    typedef std::vector< std::pair< std::vector<size_t>, ManifoldType > >  manifoldContainer;

    /// initialises / overwrites the manifold information stored in the VData object
    void AddNodeManifolds( manifoldContainer::const_iterator first,
                           manifoldContainer::const_iterator last );
                          
    /// resizes for a single-element type mesh
    void Resize( size_t nodes_per_element, size_t nbors_per_element, int8_t etype, size_t nodes, size_t elmts );
   
    void Resize( const std::deque<int8_t>& etypes,
                 const std::deque<uint32_t>& npes,  ///< sizes for resizing the member vectors
                 const std::deque<uint32_t>& epes,
                 size_t nodes, size_t faces, size_t interfaces );

    /// for hybrid element meshes, increasing the storage for element types to the new size without invalidating existing types unless the storage is shrunk
    void ResizeElementTypes( size_t elements );

    void ResizeNodes( size_t nodes );
    void ResizePlist( size_t elements );
    void ResizePlist( size_t elements, uint32_t nperelmt );
    void ResizeElementNodes( size_t eid, uint32_t nperelmt );
    void ResizePlist( const std::deque<uint32_t>& mixed_ele_plist );
    void ResizePfverts( size_t elements );
    void ResizePfverts( size_t elements, uint32_t nperelmt );
    void ResizeElementNeighbors( size_t eid, uint32_t nperelmt );
    void ResizePfverts( const std::deque<uint32_t>& mixed_ele_pfverts );
    void ResizeBFlags( /* nodes */ );
    void ResizeBREP_Flags( /* nodes */ );
    
    /// compares 'pelmt', 'plist', 'pvferts' and 'pbflags' among the VData; specific mismatches are reported to std::cerr
    bool operator==( const VData& ) const;
    
    /// using element types, coordinate range, and boundary flags, asesses whether this is a 1D, 2D , or three dimensional model
    uint32_t SpatialDimension() const;

    /// returns number of nodes in the mesh
    size_t Vertices() const;
  
    /// returns number of elements (not including faces and interfaces)
    size_t Elements() const;
  
    /// number of lower-dimensional elements that know their higher-dimensional neighbors
    size_t Faces() const;
  
    /// number of split-boundary faces; faces with multiplicated nodes at boundaries
    size_t Interfaces() const;
  
    /// combined number of any entities: elements + faces + interfaces
    size_t Cells() const;
    
    /// topologically collocated nodes along the node-matched boundaries of mesh patches
    size_t NodeManifolds() const;

    /// number of unique CSMP finite-element/cell types contained in the mesh
    size_t ElementTypes() const;
  
    /// number of element neighbors in a single element type mesh
    size_t ElementNeighbors() const;
  
    /// reports the number of different nde flags stored as boundary flags in the model
    size_t BFlags() const;
  
    /// number of nodes stored for element, face or interface
    uint32_t PlistSize( size_t eidx ) const;

    /// number of neighbors stored for element, face or interface
    uint32_t PfvertsSize( size_t eidx ) const;
    
    /// get type; a mesh is of hybrid-element type if it contains multiple element types
    bool   HybridElementTypeMesh() const;
  
    /// set whether mesh contains multiple element types
    void   HybridElementTypeMesh( bool hybrid_mesh );
    
    /// reports whether the model contains only isoparametric element types
    bool   IsoparametricElementMesh() const;
  
    /// if numbering is not 0..n-1, this method establishes thos
    void   EstablishZeroBasedNumbering();
  
    /// x-coordinate of node i
    void    Px( size_t i, double );
    void    Py( size_t i, double );
    void    Pz( size_t i, double );
    double  Px( size_t i ) const;
    double  Py( size_t i ) const;
    double  Pz( size_t i ) const;
 
    /// range of the X coordinate (increasing to the east)
    std::pair<double,double>  X_Range() const;

    /// range of the Y coordinate (increasing upward and often indicating elevation above sealevel)
    std::pair<double,double>  Y_Range() const;

    /// range of the Z coordinate (from north to south)
    std::pair<double,double>  Z_Range() const;
  
    /// calculates min-max vertex coordinate values of mesh stored in VData
    void CoordinateRange( char coordinate_axis, double& cmin, double& cmax ) const;
  
    /// rescales vertex coordinate values in given spatial direction
    void ScaleCoordinateToRange( char coordinate_axis, double cmin, double cmax ); 
    
    /// to set vertex=node coordinate of node i for user defined coordinate component (x,y, or z)
    void  P( uint32_t coordinate_axis, size_t i, double );

    /// to get vertex=node coordinate of node i for user defined coordinate component (x,y, or z)
    double P( uint32_t coordinate_axis, size_t i ) const;

    /// set CSMP finite element type of element in 'pelmt' container
    void   ElementType( size_t eidx, int8_t type );
  
    /// get CSMP finite element type from 'pelmt' container; stores only 1 elmt in single element-type mesh
    int8_t ElementType( size_t eidx ) const;

    /// assuming that all interpolation functions have same order, returns that order
    uint32_t OrderOfFiniteElementInterpolationFunctions() const;
 
    /// set node index of element in serialised array of node ids; use pelmt to determine how many nodes there shoud be
    void   Plist( size_t eidx, uint32_t node, size_t val );

    /// get node index of element in serialised array of node ids; use pelmt to determine how many nodes there shoud be
    size_t  Plist( size_t eidx, uint32_t node ) const;
  
    /// set neighbor element index (or boundary identifier) for neighbor i of element eidx
    void   Pfvert( size_t eidx, uint32_t i, int64_t val );

    /// get neighbor element index (or boundary identifier) for neighbor i of element eidx
    int64_t  Pfvert( size_t eidx, uint32_t i ) const;

    /// adds id (0..n-1) of boundary node and its BOX_BOUNDARY flag (negative integer)
    void BFlag( size_t node_id, std::int8_t bflag );

    /// returns the boundary flag BOX_BOUNDARY  of the node
    std::int8_t BFlag( size_t node_id ) const;

     /// adds id (0..n-1) of boundary node and its BOX_BOUNDARY flag (negative integer)
    void BREP_Flag( size_t node_id, std::int8_t brep_flag );

    /// returns the boundary flag BOX_BOUNDARY  of the node
    std::int8_t BREP_Flag( size_t node_id ) const;

    /// returns the box boundary identifier of the node if it is located on the model boundary; else returs NOT
    std::int8_t BoundaryFlag( size_t vertex ) const;
  
    /// empties 'pfverts' container if the contained info is flaky so that later code is prompted to recreate it
    void RemovePfverts() { pfverts.clear(); }
  
    /// empties map which stores which nodes lie at boundary and what there BOX_BOUNDARY flag is
    void RemoveBflags();

    // iterators (for any element, face or interface)
    std::deque<std::vector<size_t> >::iterator    PlistBegin();
    std::deque<std::vector<size_t> >::iterator    PlistEnd();
    std::deque<std::vector<int64_t> >::iterator   PfvertsBegin();
    std::deque<std::vector<int64_t> >::iterator   PfvertsEnd();

    std::vector<size_t>::iterator                 PlistBegin( size_t eidx );
    std::vector<size_t>::iterator                 PlistEnd( size_t eidx );
    std::vector<int64_t>::iterator                PfvertsBegin( size_t eidx );
    std::vector<int64_t>::iterator                PfvertsEnd( size_t eidx );

    std::vector<std::int8_t>::iterator            BFlagsBegin();
    std::vector<std::int8_t>::iterator            BFlagsEnd();
    std::vector<std::int8_t>::const_iterator      BFlagsBegin() const;
    std::vector<std::int8_t>::const_iterator      BFlagsEnd() const;
    
    /// labels of the nodes that indicate which point, line or surface of the Boundary Representation (BREP) of the original model it represents
    std::vector<std::int8_t>::iterator            BREP_FlagsBegin();
    std::vector<std::int8_t>::iterator            BREP_FlagsEnd();
    std::vector<std::int8_t>::const_iterator      BREP_FlagsBegin() const;
    std::vector<std::int8_t>::const_iterator      BREP_FlagsEnd() const;

    // const iterators
    std::vector<int8_t>::const_iterator                PelmtBegin() const;
    std::vector<int8_t>::const_iterator                PelmtEnd() const;
    std::deque<std::vector<size_t> >::const_iterator   PlistBegin() const;
    std::deque<std::vector<size_t> >::const_iterator   PlistEnd() const;
    std::deque<std::vector<int64_t> >::const_iterator  PfvertsBegin() const;
    std::deque<std::vector<int64_t> >::const_iterator  PfvertsEnd() const;
  
    /// checks whether pfverts, has right size and contains plausible information (no guarantees!)
    bool WithNeighbourConnectivity() const;

    std::vector<size_t>::const_iterator                PlistBegin( size_t eidx ) const;
    std::vector<size_t>::const_iterator                PlistEnd( size_t eidx ) const;
    std::vector<int64_t>::const_iterator               PfvertsBegin( size_t eidx ) const;
    std::vector<int64_t>::const_iterator               PfvertsEnd( size_t eidx ) const;
    
    // specific element, face and interface iterators
    /// iterator to CSMP finite element type of first face stored in mesh
    std::vector<int8_t>::const_iterator                 PelmtFacesBegin() const;
    /// iterator to CSMP finite element type of first interface stored in mesh
    std::vector<int8_t>::const_iterator                 PelmtInterfacesBegin() const;
  
    // node iterators for subsets of the Plist
    /// Iterator to beginning of elements in the Plist
    std::deque<std::vector<size_t> >::const_iterator    PlistElmtsBegin() const;
    /// Iterator to end of elements in the Plist
    std::deque<std::vector<size_t> >::const_iterator    PlistElmtsEnd() const;
    /// Iterator to beginning of faces in the Plist
    std::deque<std::vector<size_t> >::const_iterator    PlistFacesBegin() const;
    /// Iterator to end of faces in the Plist
    std::deque<std::vector<size_t> >::const_iterator    PlistFacesEnd() const;
    /// Iterator to beginning of interfaces in the Plist
    std::deque<std::vector<size_t> >::const_iterator    PlistInterFacesBegin() const;
    /// Iterator to end of interfaces in the Plist
    std::deque<std::vector<size_t> >::const_iterator    PlistInterFacesEnd() const;

    /// neighbor iterator for first face in plist (use PlistInterFacesBegin() to find last one)
    std::deque<std::vector<int64_t> >::const_iterator    PfvertsFacesBegin() const;
    /// neighbor iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
    std::deque<std::vector<int64_t> >::const_iterator    PfvertsInterfaceBegin() const;
 
    /// read / write access to the stored node manifolds
    manifoldContainer::iterator                          PmanifoldsBegin();
    manifoldContainer::iterator                          PmanifoldsEnd();
 
    manifoldContainer::const_iterator                    PmanifoldsBegin() const;
    manifoldContainer::const_iterator                    PmanifoldsEnd() const;
 
    // EXTRA DATA, MESH MODIFICATION AND REPAIR

    /// flips clockwise-numbered elements, into counter-clockwise right-hand rule compliant orientation; lower dimensional elements are made consistent; returns how many were flipped
    size_t RenumberElementsCounterClockwise2D();
    
    /// computes neighbor element connectivity between line and surface elements, faces and interfaces and replaces existing connectivity with it
    void   EstablishElementConnectivity2D();
    
    /// rebuilds 'pfverts' from scratch, needs correect BOX_BOUNDARY flagging of nodes to get the non-neighbors right
    void   EstablishElementConnectivity3D(); // retested: OK 3/12/21 by SKM
    
    /// eliminates corner tetrahedra with all nodes on the model boundary; extra element degrees of freedom are introduced for boundary condition assignment
    size_t RemeshCornerSpanningTetrahedra();
    
    /// creates an extra array 'pnode' equivalent to a sparsity pattern recording to which nodes each node pnode[i] is connected to
    void   EstablishNodeNeighborConnectivity( std::vector<std::set<size_t>>& pnode ) const;

    /// initialises the TOPOTYPE flags, indicating which elements of the model boundary representation the node forms part of / represents
    void   InitialiseNodeTopologyIdentifiers();
    
    /// eliminate nodes that are not connected to any element, face or interface; report whether there were any
    bool DetectAndEliminateOrphanNodes( bool eliminate_orphan_nodes=true );

    /// eliminate cells with identical nodes, except for interfaces
    bool DetectAndEliminateDuplicateCells( bool verbose=true );
      
    /// clear the container
    void Erase();

    /// wrtie connectivity structure to ASCII text file
    void OutASCII( const char* file ) const;
    
    /// writes initialiser lists for the current VData in C++17 format
    void OutCPP17( std::ofstream& ) const;
  
    /// initialise VData=mesh connectivity structures from binary file
    void InText( std::ifstream& );
  
    /// print connectivity information to screen
    void Out() const;

  protected:
  
    /// remove any gaps in the numbering of nodes and elements
    bool   CheckFix();

    /// aligns potential line elements in a 2D mesh, those at boundary are given the same orientation as the surface-element boundary faces
    void   CreateConsistentLineElementOrientations2D();
    
    /// angle between line elements in degrees
    double AngleBetweenLineElements2D( size_t elmt1, size_t elmt2 );
    double AngleBetweenLineElements3D( size_t elmt1, size_t elmt2 );
    
    /// computes the unit normals to the surfaces and then returns the angle between them
    double AngleBetweenSurfaceElements3D( size_t elmt1, size_t elmt2 );
    
    /// reconnects triangular elements with 3 nodes on the model boundary by switching nodes with their only neighbor; @note needs valid 'pfverts'
    size_t SwitchCornerTriangles2D();
    
    void ReduceTo( const std::map<size_t,size_t>& old_and_new_elmt_ids, std::map<size_t,size_t>& o_n_node_ids );

     /// write mesh to supplied binary file
    void OutBinary( std::fstream& ) const;
  
    /// read mesh from supplied binary file
    void InBinary( std::fstream& );
  
 private:

    bool                               hybrid_mesh_;      ///< mesh that consists of different element types
    std::vector<double>                px, py, pz;        ///< node coordinates
    std::vector<int8_t>                pelmt;             ///< CSMP element type info, needed to read plist & pfverts
    std::deque<std::vector<size_t> >   plist;             ///< nodes of each element, face and interface in that order
    std::deque<std::vector<int64_t> >  pfverts;           ///< element neighbors; same range as eidx, but also negative values possible
    std::vector<std::int8_t>           bflags;            ///< int_8 enumeration flags for those nodes that lie on model boundaries
    std::vector<std::int8_t>           gflags_;           ///< int_8 enumeration flags distinguishing mesh nodes that contribute to the model topology / geometry
    size_t                             first_face_;       ///< faces come after elements; if none this is equal to elements
    size_t                             first_interface_;  ///< interfaces come after faces; if none this is equal to elements
    // collocated nodes connecting mesh patches, and their flags
    manifoldContainer                  pmanifolds_;       ///< node manifolds, are added separately: @todo must be constructed separately

    friend class VData_Test;
};

// RELATED FUNCTIONS

/// converts corner tetrahedron and its neighbor into 3 tetrahedral cells, each with a face on the sides of the box
void splitCornerTetrahedron( VData&, size_t cnr, size_t only_neighbor );

/// returns the 3D bounding box of the element
std::array<double,3>  boundingBox( const VData&, size_t elmt );

/// thus far, only writes tetrahedra to file, appending the element number to the name
void elementToVTK( const VData&, size_t eidx, const char* outfile );

/// prints number of nodes, barycentre, and  node coordinates  placed in the origin and normalised to 1
void printCell( const VData&, size_t cell_id );

/// returns the barycentre of the celll which always has three dimensions but the higher ones will be zero in 2 and 1D models
std::array<double,3> cellBaryCenter( const VData&, size_t cell_id );

} // csmp

#endif 



