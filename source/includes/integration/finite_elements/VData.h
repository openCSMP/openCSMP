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
    
    void SingleElementType( int32 etype );
    
    void ElementTypes( const std::vector<int32>& elmt_types );
    
    void AddElementTypes( std::vector<int32>::const_iterator first,
                          std::vector<int32>::const_iterator last );

    void Resize( size_t nodes_per_element, size_t nbors_per_element, int32 etype, size_t nodes, size_t elmts );
   
    void Resize( const std::deque<int32>& etypes,
                 const std::deque<size_t>& npes, 
                 const std::deque<size_t>& epes, 
                 size_t nodes, size_t faces, size_t interfaces );
  
    void ResizeNodes( size_t nodes );
    // TODO: refactor this dangerous method as it may corrupt vdata if not used wisely
    void ResizeElementTypes( size_t elements );
    void ResizePlist( size_t elements );
    void ResizePlist( size_t elements, size_t nperelmt );
    void ResizeElementNodes( size_t eid, size_t nperelmt );
    void ResizePlist( const std::deque<size_t>& mixed_ele_plist );
    void ResizePfverts( size_t elements );
    void ResizePfverts( size_t elements, size_t nperelmt );
    void ResizeElementNeighbors( size_t eid, size_t nperelmt );
    void ResizePfverts( const std::deque<size_t>& mixed_ele_pfverts );

    virtual ~VData();
    VData( const VData& );
    VData( VData&& );
    VData& operator=( const VData& );
    VData& operator=( VData&& );
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
    size_t Simplices() const;

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
  
    /// remove any gaps in the numbering of nodes and elements
    bool   CheckFix();

    /// x-coordinate of node i
    void      Px( size_t i, double64 );
    void      Py( size_t i, double64 );
    void      Pz( size_t i, double64 );
    double64  Px( size_t i ) const;
    double64  Py( size_t i ) const;
    double64  Pz( size_t i ) const;
  
    /// to set vertex=node coordinate of node i for user defined coordinate component (x,y, or z)
    void      P( size_t coordinate_axis, size_t i, double64 );

    /// to get vertex=node coordinate of node i for user defined coordinate component (x,y, or z)
    double64  P( size_t coordinate_axis, size_t i ) const;

    /// set CSMP finite element type of element in 'pelmt' container
    void   ElementType( size_t eidx, int32 type );
  
    /// get CSMP finite element type from 'pelmt' container; stores only 1 elmt in single element-type mesh
    int32  ElementType( size_t eidx ) const;

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
    void AddBFlag( size_t node_id, long64 bflag );
  
    /// calculates min-max vertex coordinate values of mesh stored in VData
    void CoordinateRange( char coordinate_axis, double64& cmin, double64& cmax ) const;
  
    /// rescales vertex coordinate values in given spatial direction
    void ScaleCoordinateToRange( char coordinate_axis, double64 cmin, double64 cmax ); 
  
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

    std::map<size_t,long64>::iterator             BFlagsBegin();
    std::map<size_t,long64>::iterator             BFlagsEnd();
    std::map<size_t,long64>::const_iterator       BFlagsBegin() const;
    std::map<size_t,long64>::const_iterator       BFlagsEnd() const;

    // const iterators
    std::vector<int32>::const_iterator                 PelmtBegin() const;
    std::vector<int32>::const_iterator                 PelmtEnd() const;
    std::deque<std::vector<size_t> >::const_iterator   PlistBegin() const;
    std::deque<std::vector<size_t> >::const_iterator   PlistEnd() const;
    std::deque<std::vector<long64> >::const_iterator   PfvertsBegin() const;
    std::deque<std::vector<long64> >::const_iterator   PfvertsEnd() const;

    std::vector<size_t>::const_iterator                PlistBegin( size_t eidx ) const;
    std::vector<size_t>::const_iterator                PlistEnd( size_t eidx ) const;
    std::vector<long64>::const_iterator                PfvertsBegin( size_t eidx ) const;
    std::vector<long64>::const_iterator                PfvertsEnd( size_t eidx ) const;
  
    // specific element, face and interface iterators
    /// iterator to CSMP finite element type of first face stored in mesh
    std::vector<int32>::const_iterator                 PelmtFacesBegin() const;
    /// iterator to CSMP finite element type of first interface stored in mesh
    std::vector<int32>::const_iterator                 PelmtInterfacesEnd() const;
  
    /// node iterator for first face in plist (use PlistInterFacesBegin() to find last one)
    std::deque<std::vector<size_t> >::const_iterator    PlistFacesBegin() const;
    /// node iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
    std::deque<std::vector<size_t> >::const_iterator    PlistInterFacesBegin() const;

    /// neighbor iterator for first face in plist (use PlistInterFacesBegin() to find last one)
    std::deque<std::vector<long64> >::const_iterator    PfvertsFacesBegin() const;
    /// neighbor iterator for first interface plist; equivalent to PlistFacesEnd; use PlistEnd() for last one
    std::deque<std::vector<long64> >::const_iterator    PfvertsInterfaceBegin() const;
  

    /// write mesh to supplied binary file
    void OutBinary( std::FILE* fp ) const;
  
    /// read mesh from supplied binary file
    void InBinary( std::FILE* fp );
  
    /// initialise VData=mesh connectivity structures from binary file
    void InText( std::ifstream& ifs );
  
    /// wrtie connectivity structure to ASCII text file
    void OutASCII( const char* file ) const;
  
    /// print connectivity information to screen
    void Out() const;
  
    /// eliminate nodes that are not connected to any element, face or interface; report whether there were any
    bool DetectAndEliminateOrphanNodes( bool eliminate_orphan_nodes=true );
  
    /// clear the container
    void Erase();

  protected:

    void ReduceTo( const std::map<size_t,size_t>& old_and_new_elmt_ids, std::map<size_t,size_t>& o_n_node_ids );
 
    std::vector<int32>                pelmt;            ///< CSMP element type info, needed to read plist & pfverts

  private:

    bool                              hybrid_mesh_;      ///< mesh that consists of different element types
    std::vector<double64>             px, py, pz;        ///< node coordinates
    std::deque<std::vector<size_t> >  plist;             ///< nodes of each element
    std::deque<std::vector<long64> >  pfverts;           ///< element neighbors; same range as eidx, but also negative values possible
    std::map<size_t,long64>           bflags;            ///< flags for those nodes that lie on model boundary
    // NEW
    size_t                            first_interface_;  ///< faces come after elements; if none this is equal to elements
    size_t                            first_face_;       ///< interfaces come after faces; if none this is equal to elements 

    friend class VData_Test;
};

} // csmp

#endif 


