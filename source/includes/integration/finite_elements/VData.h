#ifndef CSMP_VDATA_H
#define CSMP_VDATA_H

#include  "CSMP_definitions.h"

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

*/
class VData {

  public:

    VData();

    VData( const std::deque<size_t>& npes, 
           const std::deque<size_t>& epes,
           size_t nodes );
                        
    VData( size_t nodes_per_element, size_t nbors_per_element, size_t nodes, size_t elmts );
    
    void SingleElementType( int8_t etype );
    
    void ElementTypes( const std::vector<int8_t>& elmt_types );
    
    void AddElementTypes( std::vector<int8_t>::const_iterator first,
                          std::vector<int8_t>::const_iterator last );

    void AddElementTypes( std::deque<int8_t>::const_iterator first,
                          std::deque<int8_t>::const_iterator last );

    void Resize( size_t nodes_per_element, size_t nbors_per_element, int8_t etype, size_t nodes, size_t elmts );
   
    void Resize( const std::deque<int8_t>& etypes,
                 const std::deque<size_t>& npes, 
                 const std::deque<size_t>& epes, 
                 size_t nodes, size_t faces, size_t interfaces );
  
    void ResizeNodes( size_t nodes );
    void ResizePlist( size_t elements );
    void ResizePlist( size_t elements, size_t nperelmt );
    void ResizeElementNodes( size_t eid, size_t nperelmt );
    void ResizePlist( const std::deque<size_t>& mixed_ele_plist );
    void ResizePfverts( size_t elements );
    void ResizePfverts( size_t elements, size_t nperelmt );
    void ResizeElementNeighbors( size_t eid, size_t nperelmt );
    void ResizePfverts( const std::deque<size_t>& mixed_ele_pfverts );
    void ResizeBFlags( /* nodes */ );

    virtual ~VData();
    VData( const VData& );
    VData( VData&& ) = default;
    VData& operator=( const VData& );
    VData& operator=( VData&& ) = default;
    bool   operator==( const VData& ) const;

    /// returns number of nodes in the mesh
    size_t Vertices() const;
  
    /// returns number of elements (not including faces and interfaces)
    size_t Elements() const;
  
    /// number of lower-dimensional elements that know their higher-dimensional neighbors
    size_t Faces() const;
  
    /// number of split-boundary faces; faces with multiplicated nodes at boundaries
    size_t InterFaces() const;
  
    /// combined number of any entities: elements + faces + interfaces
    size_t TotalNumberOfCells() const;

    /// CSMP types of the elements, faces and interfaces contained in the mesh
    size_t ElementTypes() const;
  
    /// number of element neighbors in a single element type mesh
    size_t ElementNeighbors() const;
  
    /// number of nodes at model boundaries
    size_t BFlags() const;
  
    /// number of nodes stored for element, face or interface
    size_t PlistSize( size_t eidx ) const;

    /// number of neighbors stored for element, face or interface
    size_t PfvertsSize( size_t eidx ) const;
    
    /// get type; a mesh is of hybrid-element type if it contains multiple element types
    bool   HybridElementTypeMesh() const;
  
    /// set whether mesh contains multiple element types
    void   HybridElementTypeMesh( bool hybrid_mesh );
  
    /// if numbering is not 0..n-1, this method establishes thos
    void   EstablishZeroBasedNumbering();
  
    /// x-coordinate of node i
    void      Px( size_t i, double64 );
    void      Py( size_t i, double64 );
    void      Pz( size_t i, double64 );
    double64  Px( size_t i ) const;
    double64  Py( size_t i ) const;
    double64  Pz( size_t i ) const;
 
    /// range of the X coordinate (increasing to the east)
    std::pair<double64,double64>  X_Range() const;

    /// range of the Y coordinate (increasing upward and often indicating elevation above sealevel)
    std::pair<double64,double64>  Y_Range() const;

    /// range of the Z coordinate (from north to south)
    std::pair<double64,double64>  Z_Range() const;
  
    /// calculates min-max vertex coordinate values of mesh stored in VData
    void CoordinateRange( char coordinate_axis, double64& cmin, double64& cmax ) const;
  
    /// rescales vertex coordinate values in given spatial direction
    void ScaleCoordinateToRange( char coordinate_axis, double64 cmin, double64 cmax ); 
    
    /// to set vertex=node coordinate of node i for user defined coordinate component (x,y, or z)
    void      P( size_t coordinate_axis, size_t i, double64 );

    /// to get vertex=node coordinate of node i for user defined coordinate component (x,y, or z)
    double64  P( size_t coordinate_axis, size_t i ) const;

    /// set CSMP finite element type of element in 'pelmt' container
    void   ElementType( size_t eidx, int8_t type );
  
    /// get CSMP finite element type from 'pelmt' container; stores only 1 elmt in single element-type mesh
    int8_t ElementType( size_t eidx ) const;

    /// assuming that all interpolation functions have same order, returns that order
    size_t OrderOfFiniteElementInterpolationFunctions() const;
 
    /// set node index of element in serialised array of node ids; use pelmt to determine how many nodes there shoud be
    void   Plist( size_t eidx, size_t node, size_t val );

    /// get node index of element in serialised array of node ids; use pelmt to determine how many nodes there shoud be
    size_t Plist( size_t eidx, size_t node ) const;
  
    /// set neighbor element index (or boundary identifier) for neighbor i of element eidx
    void   Pfvert( size_t eidx, size_t i, int32 val );

    /// get neighbor element index (or boundary identifier) for neighbor i of element eidx
    long64 Pfvert( size_t eidx, size_t i ) const;

    /// adds id (0..n-1) of boundary node and its BOX_BOUNDARY flag (negative integer)
    void AddBFlag( size_t node_id, std::int8_t bflag );
  
    /// returns the box boundary identifier of the node if it is located on the model boundary; else returs NOT
    std::int8_t BoundaryFlag( size_t vertex ) const;
  
    /// empties 'pfverts' container if the contained info is flaky so that later code is prompted to recreate it
    void RemovePfverts() { pfverts.clear(); }
  
    /// empties map which stores which nodes lie at boundary and what there BOX_BOUNDARY flag is
    void RemoveBflags();

    // iterators (for any element, face or interface)
    std::deque<std::vector<size_t> >::iterator    PlistBegin();
    std::deque<std::vector<size_t> >::iterator    PlistEnd();
    std::deque<std::vector<long64> >::iterator    PfvertsBegin();
    std::deque<std::vector<long64> >::iterator    PfvertsEnd();

    std::vector<size_t>::iterator                 PlistBegin( size_t eidx );
    std::vector<size_t>::iterator                 PlistEnd( size_t eidx );
    std::vector<long64>::iterator                 PfvertsBegin( size_t eidx );
    std::vector<long64>::iterator                 PfvertsEnd( size_t eidx );

    std::vector<std::int8_t>::iterator            BFlagsBegin();
    std::vector<std::int8_t>::iterator            BFlagsEnd();
    std::vector<std::int8_t>::const_iterator      BFlagsBegin() const;
    std::vector<std::int8_t>::const_iterator      BFlagsEnd() const;

    // const iterators
    std::vector<int8_t>::const_iterator                PelmtBegin() const;
    std::vector<int8_t>::const_iterator                PelmtEnd() const;
    std::deque<std::vector<size_t> >::const_iterator   PlistBegin() const;
    std::deque<std::vector<size_t> >::const_iterator   PlistEnd() const;
    std::deque<std::vector<long64> >::const_iterator   PfvertsBegin() const;
    std::deque<std::vector<long64> >::const_iterator   PfvertsEnd() const;
  
    /// checks whether pfverts, has right size and contains plausible information (no guarantees!)
    bool WithNeighbourConnectivity() const;

    std::vector<size_t>::const_iterator                PlistBegin( size_t eidx ) const;
    std::vector<size_t>::const_iterator                PlistEnd( size_t eidx ) const;
    std::vector<long64>::const_iterator                PfvertsBegin( size_t eidx ) const;
    std::vector<long64>::const_iterator                PfvertsEnd( size_t eidx ) const;
    
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
    std::deque<std::vector<long64> >::const_iterator    PfvertsFacesBegin() const;
    /// neighbor iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
    std::deque<std::vector<long64> >::const_iterator    PfvertsInterfaceBegin() const;
 
 
    // EXTRA DATA, MESH MODIFICATION AND REPAIR 

    /// remove any gaps in the numbering of nodes and elements
    bool   CheckFix();

    /// flips clockwise-numbered elements, into counter-clockwise right-hand rule compliant orientation; lower dimensional elements are made consistent; returns how many were flipped
    size_t RenumberElementsCounterClockwise2D();
    
    /// computes neighbor element connectivity between line and surface elements, faces and interfaces and replaces existing connectivity with it
    void   EstablishElementConnectivity2D();
    
    /// rebuilds 'pfverts' from scratch
    void   EstablishElementConnectivity3D();

    
    // PERSISTANCE (storing mesh in binary file)
    
    /// vertex manifolds: key=-vertex index, value = set of pairs of nodes and their INSIDE,OUTSIDE, MIDDLE classifers
    typedef std::map<size_t,std::set<std::pair<size_t,int8_t> > > vertexManifoldIndices;
    
    /// checks for collocated vertices into transfer data structure
    bool ExtractNodeManifolds( vertexManifoldIndices& ) const;
  
    /// write mesh to supplied binary file
    void OutBinary( std::fstream& ) const;
  
    /// read mesh from supplied binary file
    void InBinary( std::fstream& );
  
    /// wrtie connectivity structure to ASCII text file
    void OutASCII( const char* file ) const;
  
    /// initialise VData=mesh connectivity structures from binary file
    void InText( std::ifstream& );
  
    /// print connectivity information to screen
    void Out() const;
  
    /// eliminate nodes that are not connected to any element, face or interface; report whether there were any
    bool DetectAndEliminateOrphanNodes( bool eliminate_orphan_nodes=true );
  
    /// clear the container
    void Erase();

  protected:
  
    /// aligns potential line elements in a 2D mesh, those at boundary are given the same orientation as the surface-element boundary faces
    void   CreateConsistentLineElementOrientations2D();
    
    /// angle between line elements in degrees
    double64 AngleBetweenLineElements2D( size_t elmt1, size_t elmt2 );
    
    void ReduceTo( const std::map<size_t,size_t>& old_and_new_elmt_ids, std::map<size_t,size_t>& o_n_node_ids );

  private:

    bool                              hybrid_mesh_;      ///< mesh that consists of different element types
    std::vector<double64>             px, py, pz;        ///< node coordinates
    std::vector<int8_t>               pelmt;             ///< CSMP element type info, needed to read plist & pfverts
    // although there's little point to having 64-bit pointers but not 64-bit sizes, after all.
    std::deque<std::vector<size_t> >  plist;             ///< nodes of each element, face and interface in that order
    std::deque<std::vector<long64> >  pfverts;           ///< element neighbors; same range as eidx, but also negative values possible
    // all enums / flags must fit into 8-bit integers
    std::vector<std::int8_t>          bflags;            ///< flags for those nodes that lie on model boundary
    size_t                            first_face_;       ///< faces come after elements; if none this is equal to elements
    size_t                            first_interface_;  ///< interfaces come after faces; if none this is equal to elements 

    friend class VData_Test;
};

} // csmp

#endif 


