#ifndef CSMP_MODEL_SUB_DOMAIN_H
#define CSMP_MODEL_SUB_DOMAIN_H

#include "CSMP_definitions.h"
#include "LocalVariableStorage.h"

namespace csmp {

/// Parts of ModelSubDomain
enum SUBDOMAIN_PART { COMPLETE,
                      INTERIOR,
                      PERIMETER };

SUBDOMAIN_PART parseSubdomainPart( const char* subdomain );
std::string    parseSubdomainPart( SUBDOMAIN_PART part );

class FiniteElementManager;
template<uint32_t> class PropertyDatabase;
template<uint32_t> class MeshManager;
template<uint32_t> class FiniteVolumeStencilManager;
template<uint32_t> class Point;
template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Interrelation;
template<uint32_t> class Visitor;
class PropertyConstraints;

/**
    Complete index specifications of a ModelSubDomain
    used in the reconstruction of subdomains from binary files.

@author SKM
@date 14/3/2016
*/
struct SubDomainInfo {
   std::string         name;                         ///< unique name
   std::vector<size_t> interior_elmts;               ///< cells that have no face on the perimeter
   std::vector<size_t> perimeter_elmts;              ///< cells that have at least one face on perimeter
   std::vector<std::vector<uint32_t> > perimeter_faces;  ///< local 0..faces-1 identifiers of the faces of the simplices that lie on domain boundary
   std::vector<size_t> interior_nodes;               ///< nodes within the subdomain
   std::vector<size_t> perimeter_nodes;              ///< nodes on the perimeter of the subdomain
};



/** @brief Storage of Node and Cell pointers to subdomains of the overall computational model.
    Distinguishes interior from perimeter nodes (located on subdomain boundary) and interior and perimeter cells
    which have at least a single Face on the subdomain boundary.
    Base class to Region, Boundary and SplitBoundary.
    
    @attention Region, Bounday, SplitBoundary ARE AGNOSTIC ABOUT MESH CONNECTIVITY.
    ModelSubDomain objects do not create Nodes, Elements, Faces or Interfaces. This is the role of the MeshManager
    that is agnostic about subdomains but maintains and updates the mesh.
 
    @note Boundary and SplitBoundare are ModelSubDomain subclasses consisting of lower-dimensional Face and
    InterFace objects.
    - SplitBoundaries are Internal boundaries connecting node-matched but disconnected mesh patches
    - Boundaries tend to be external boundaries with names corresponding to the sides of box-shaped
      models or other unique names.
*/
template<uint32_t dim,template<uint32_t> class CELL>
class ModelSubDomain {
  public:    
    /// constructs incomplete subdomain for later initialisation with suitable methods in subclasses
    ModelSubDomain( const std::string& subdomain_name, const PropertyDatabase<dim>&, bool unique );
    ModelSubDomain( const ModelSubDomain& );
    ModelSubDomain( ModelSubDomain&& );
  
    virtual ~ModelSubDomain() noexcept;
    
    ModelSubDomain<dim,CELL>&  operator=( const ModelSubDomain& ) noexcept;
    ModelSubDomain<dim,CELL>&  operator=( ModelSubDomain&& ) noexcept;

    std::string Name() const noexcept;
    void Name( const std::string& ) noexcept;

    /// local variable storage interface
    virtual PLACEMENT Placement() const = 0;
    virtual bool      ValidVariable( const char* variableName ) const noexcept;

    /// support of the Visitor pattern
    virtual void Accept( Visitor<dim>& ) = 0;
    
    /// modification via Operand-based relations between discretised properties that only modify a single result variable
    void Apply( Interrelation<dim>& );
    
    /// distinguishes PERIMETER simplices that have at least one face on domain boundary from INTERIOR ones; calls PartitionElementVector()
    void IdentifyPerimeter();
    
    /// stores topological  information on the mesh that is associated with PERIMETER, boundary and split boundary information
    void UpdateTopoTypeNodeFlags();
    
    /// creates node pointer vector from cell vector, using a vector to achieve uniqueness via sort, unique, erase algorithms; NOTE: does not distinguish interior and perimeter nodes!
    void CreateNodePointerVector();

    /// creates node pointer vector from the shared face nodes of the supplied range of contacting cells
    void CreateNodePointerVector( std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& contacting_cells );

    /// sorts the node and CELL vectors split into the interior and perimeter ranges (4 sorting operations)
    void SortVectors( size_t interior_cells, size_t interior_nodes );
    
    /// assuming that a partitioned (and sorted) cell vector is in place, constructs the bd_face_vec_ by checking whether neighbor cells belong to the domain or not
    void BuildPerimeterFaceVector( int64_t interior_cells );
    
    /// for rebuilding subdomains when nodes or cells changed:  flag up for a rebuild using RebuildSubDomainAfterChangeOfCellVector
    void ScheduleForRebuild() noexcept;
    bool NeedsRebuild() const noexcept;
    
    /// deletes nullptr cells, rebuilds node vector, sorts everything and re-establishes the perimeter face vectors after modifications of cells
    void RebuildSubDomainAfterChangeOfCellVector();
    
    /// rebuilds node vector, sorts everything and re-establishes interior and perimeter
    void RebuildSubDomainAfterChangeOfNodeVector();
    
    /// rebuilds subdomain on the basis of the cells that will be selected according to the supplied property constraints
    void UpdateCellMembershipApplyingConstraints( typename std::vector<CELL<dim>*>::const_iterator master_domain_start,
                                                  typename std::vector<CELL<dim>*>::const_iterator master_domain_end,
                                                  const PropertyConstraints& );

    /// deletes 'nullptr' cells from cell vector, retaining the sorting into interior and perimeter elements, and rebuilding perimeter faces as necessary
    size_t RebuildCellAndPerimeterFaceVector();

    /// deletes 'nullptr' nodes from Node vector, retaining the sorting into interior and perimeter nodes
    size_t RebuildNodeVector();


    // ----------------------------------------
    // Indexes
    // ----------------------------------------

    /// reference counting-based unique domain identifier  (0..n-1)
    int32_t DomainIndex() const noexcept;

    /// renumbers nodes in domain 0..n-1
    virtual size_t  RenumberNodes() const noexcept;
    
    /// renumbers cells in domain 0..n-1
    size_t  RenumberCells() const noexcept;
    
    /// renumber cells and nodes consecutively from 0..n-1
    void    UpdateMemberIndexes() const noexcept;
    
    /// output current indices to vector
    std::vector<size_t>  MemberCellIndexes() const noexcept;

    // ----------------------------------------
    // access
    // ----------------------------------------

    /// reference to container of finite cell pointers to either Element, Face or InterFace objects; @note used for boolean operations
    const typename std::vector<CELL<dim>*>&  CellVector() const noexcept;
    
    /// do not remove!;  used for boolean operations
    typename std::vector<CELL<dim>*>&        CellVector() noexcept;
  
    /// reference to const Node pointer vector
    const typename std::vector<Node<dim>*>&  NodeVector() const noexcept;

    /// reference to Node pointer vector
    typename std::vector<Node<dim>*>&        NodeVector() noexcept;

    /// const iterators (cell and node pointers cannot be modified but the nodes and cells can!)
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesBegin() const noexcept;
    typename std::vector<csmp::Node<dim>*>::const_iterator  PerimeterNodesBegin() const noexcept;
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesEnd() const noexcept;
    typename std::vector<CELL<dim>*>::const_iterator        CellsBegin() const noexcept;
    typename std::vector<CELL<dim>*>::const_iterator        PerimeterCellsBegin() const noexcept;
    typename std::vector<CELL<dim>*>::const_iterator        CellsEnd() const noexcept;

    /// returns how many of the nodes in the supplied iterator range also form part of the current subdomain's perimeter
    size_t SharedPerimeterNodes( typename std::vector<csmp::Node<dim>*>::const_iterator start,
                                 typename std::vector<csmp::Node<dim>*>::const_iterator end ) const;
                                 
    // see also the non-member functions below

    /// check whether subdomain conatains any cells
    bool              Empty() const noexcept;
    
    /// do not call on SplitBoundary because it only stores NodeManifolds, call nodes=RenumberNodes() on these
    size_t            Nodes() const noexcept;
    
    size_t            InteriorNodes() const noexcept;
    size_t            PerimeterNodes() const noexcept;
    size_t            IntegrationPoints() const noexcept;
    size_t            SectorIntegrationPoints() const noexcept;
    size_t            FacetIntegrationPoints() const noexcept;
    size_t            Cells() const noexcept;
    size_t            InteriorCells() const noexcept;
    size_t            PerimeterCells() const noexcept;
    
    bool              Contains( const CELL<dim>* const ) const noexcept;
    bool              Contains( const Node<dim>* const ) const noexcept;
    
    bool              IsPerimeterNode( const csmp::Node<dim>* const ) const noexcept;
    bool              IsPerimeterCell( const CELL<dim>* const ) const noexcept;
    
    /// for looping over the perimeter faces of perimeter cell with #eid (>=InteriorCells() to Cells()), method is used to travel across subdomain surface / outline
    uint32_t          PerimeterFaces( size_t cell_idx ) const;
    /// returns local face id of face #face that lies on perimeter of model subdomain
    uint32_t          PerimeterFace( size_t cell_idx, uint32_t face ) const;

    /// pointer to node #n in subdomain; @attention node can vary from initialization to initialization
    csmp::Node<dim>*  N( size_t n ) const noexcept;
    /// pointer to cell #n of model subdomain
    CELL<dim>*        E( size_t n ) const noexcept;

    /// is the node located on the surface of the model subdomain?
    bool              IsPerimeterNode( const size_t nidx ) const;
    /// does the cell have at least on face on the surface of the model subdomain?
    bool              IsPerimeterCell( const size_t eidx ) const;

    // ----------------------------------------
    // geometry
    // ----------------------------------------
    
    /// unique versus overlapping domains
    bool IsUnique() const { return is_unique_; }
    void IsUnique( bool unique_domain ) { is_unique_ = unique_domain; }
    
    /// checks whether all cells within the subdomain are interconnected (if the domain has multi-dimensional cells this is never the case)
    bool IsContiguous() const;

    /// if the model subdomain consists of a single cell shape (line, surface, volume) this method returns true revealing the type
    std::pair<CELL_SHAPE,bool>  SingleCellShapeDomain() const;

    /// returns 1) cells of how many different spatial dimensions are contained, and 2) the highest cell spatial dimension in subdomain
    std::pair<int32_t,int32_t>  SpatialDimensions() const;
    
    /// computes the "mid-point of the object
    Point<dim> Centroid() const;
    
    /// computes centre of gravity of a heterogenous subdomain as  volume*density weighted average of the element barycentres.
    Point<dim> CenterOfGravity( const csmp::Index& rho_key ) const;

    /// returns diagonally opposite points of bounding box
    void  MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const;

    /// output coordinates to user-defined vector variable [x1,y1,z1,...xn,yn,zn]
    void AssignNodeCoordinatesTo( const char* vector_prop );
  
    /// assigment of coordinate component 'x','y','z' to scalar variable of choice
    void AssignNodeCoordinatesTo( const char* scalar_prop, char coord );

    /// characteristics like 'length', 'area', 'volume' , 'aspect ratio', 'inner radius' are assigned to user-defined variable
    void AssignCellCharacteristicsTo( const char* characteristic, const char* var );
    
    /// writes a CSV (comma delimited ascii text) file with point coordinates, node-idx, BOX_BOUNDARY flags, and interior vs perimeter information
    void NodeAttributesToCSV();

    // ----------------------------------------
    // manipulation of properties
    // ----------------------------------------

    /// assigns uniform (single) variable value to either the entire subdomain or its interior or perimeter
    template<typename Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );

    /// as InputPropertyValue, but with overwrite protection for variable components that have the flag 'do_not_overwrite'
    template<typename Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value,
                             VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );

    /// changes the flag of the scalar variable 'property' to new value; applied either in the entire subdomain or its interior or perimeter
    void ChangePropertyStatus( const char* property,
                               VARIABLE_FLAG new_status_of_scalar,
                               SUBDOMAIN_PART=COMPLETE );

    /// changes the flag of the scalar variable 'property' to new value if the scalar variable is withing the specified range
    void ChangePropertyStatusWhere( const char* property,
                                    VARIABLE_FLAG new_status_of_scalar,
                                    double min_value_to_change,
                                    double max_value_to_change );

    /// changes the flags of the vector variable 'property' to new value; applied either in the entire subdomain or its interior or perimeter
    void ChangePropertyStatus( const char* property,
                               const std::vector<VARIABLE_FLAG>& new_status,
                               SUBDOMAIN_PART=COMPLETE );

    /// changes the flag of the vector variable 'property' to new value if the scalar variable is withing the specified range
    void ChangePropertyStatusWhere( const char* property,
                                    const std::vector<VARIABLE_FLAG>& new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    /// changes the flag of a particular variable component to new value; applied either in the entire subdomain or its interior or perimeter
    void ChangePropertyStatus( const char* property,
                               uint32_t component,
                               VARIABLE_FLAG new_status,
                               SUBDOMAIN_PART=COMPLETE );

    /// changes the flag of a particular variable component to new value if the scalar variable is withing the specified range
    void ChangePropertyStatusWhere( const char* property,
                                    uint32_t component,
                                    VARIABLE_FLAG new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    /// by default (i=0) returns status of scalar variable or first component of a vector or tensor variable; if i>0 flag of corresponding component is returned
    VARIABLE_FLAG  PropertyStatus( const char* variable, SUBDOMAIN_PART flag=COMPLETE , uint32_t i=0 ) const;

    /// min/max property values (length of vectors and eigenvalues of tensors)
    void MinMaxOf( const char* property,   double& gmin, double& gmax ) const;
    void MinMaxOf( const csmp::Index&,     double& gmin, double& gmax ) const;


    // ----------------------------------------
    // interpolation and extrapolation
    // ----------------------------------------

    void InterpolateNodeToCellProperty( const char* nprop, const char* eprop );
    void InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* eprop );
    void InterpolateIntegrationPointToCellProperty( const char* cprop, const char* eprop );
    void ExtrapolateCellToIntegrationPointProperty( const char* eprop, const char* cprop );
    void ExtrapolateCellToFacetIntegrationPointProperty( const char* eprop, const char* fipprop );
    void ExtrapolateCellToNodeProperty( const char* eprop, const char* nprop, bool by_distance=true );
    void ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop );


    // ----------------------------------------
    // calculations
    // ----------------------------------------

    /// arithmetic (number as opposed to volume weighted) average
    double Average( const char* property ) const;
    bool   CopyGradientOfProperty_A_To_B( const char* node_prop, const char* cell_prop );
    void   CopyReplace( const char* from, const char* to );

    /// returns the averaged unit normal to a lower-dimensional subdomain
    Point<dim> AverageUnitNormal() const;


    // ----------------------------------------
    // output
    // ----------------------------------------

    /// writes complete ModelSubDomain specifications in terms of unique indices as block to binary file
    void WriteDomainIndexesToBinaryFile( std::fstream& ) const;
    // see non-member function readDomainIndexesFromBinaryFile() in this source file for reading the indices back

    void      OutputVariableToScreen( const char* prop ) const;
    void      Out() const;
    
    friend class ModelSubDomain_Test;

  protected:

    /// establishes interior vs. exterior simplices and nodes; returns index of first boundary cell
    size_t  PartitionCellVector();
    size_t  PartitionCellVectorForBoundary();
    size_t  PartitionCellVectorForSplitBoundary();

    const PropertyDatabase<dim>&        pref_; // TODO: should not be member parameter
    std::string                         subdomain_name_="none";  ///< passed down when domain is created so that it can be referred to
    std::vector<CELL<dim>*>             cell_vec_;               ///< doubly sorted, interior cells first
    std::vector<std::vector<uint32_t> > bd_face_vec_;            ///< matching second sorted range of cell_vec_
    std::vector<csmp::Node<dim>*>       node_vec_;               ///< doubly sorted, interior nodes first
    size_t                              first_bd_node_  = std::numeric_limits<size_t>::max(); ///< begin of the perimeter nodes
    inline static int32_t               domain_count_   = 0;     ///< instance-counting to get unique identifiers regions, boundaries and splitboundaries (seperate counts)
    int32_t                             domain_idx_     = 0;     ///< unique subdomain identifier (1..n), created during region, boundary or splitboundary construction
    bool                                rebuilt_needed_ = false; ///< parameter set when mesh gets modified by MeshManager so that update can be prompted
    bool                                is_unique_      = false; ///< distinguishes space-exclusive from potentially overlapping subdomains
    static constexpr bool               verbose_ = false;

  private:
    ModelSubDomain();
};


// NON-MEMBER FUNCTIONS

/// uses subdomain name to distinguish between Region, Boundary and SplitBoundary
PLACEMENT modelSubdomainType( const std::string& subdomain_name ) noexcept;

/// returns number of nodes that are shared by the two subdomains (matches by pointers)
template<uint32_t dim,template<uint32_t> class CELL>
size_t  sharedNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );

/// returns number of nodes on the subdomain perimeters that are shared by the two subdomains (matches by pointers)
template<uint32_t dim,template<uint32_t> class CELL>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );

/// returns the shared perimeter nodes into the argument vector
template<uint32_t dim,template<uint32_t> class CELL>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>&, std::vector<Node<dim>*>& );

/// maps potential contacting faces and cell pointers of cells that are contacting each other in the two regions; @return the number of these cells, pointers to them, and corresponding face number written into the argument map
template<uint32_t dim, template<uint32_t> class CELL>
size_t  sharedPerimeterCells( const ModelSubDomain<dim,CELL>& subdomain1, const ModelSubDomain<dim,CELL>& subdomain2,
                              std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& matching_cells );
                              
/// collects the TOPOTYPEs of the nodes in the subdomain into the (unique) set that is returned
template<uint32_t dim, template<uint32_t> class CELL>
std::set<TOPOTYPE> nodeTopologyFlags( typename std::vector<Node<dim>*>::const_iterator first,
                                      typename std::vector<Node<dim>*>::const_iterator last );
                              
/// collects the TOPOTYPEs of the nodes in the SplitBoundary into the (unique) set that is returned (different method since  SplitBoundary does not store nodes)
template<uint32_t dim>
std::set<TOPOTYPE> nodeTopologyFlags( const SplitBoundary<dim>& );

/// reads ModelSubDomain data block written by writeDomainIndexesToBinaryFile() into the domain info structure
void readDomainIndexesFromBinaryFile( std::fstream&, SubDomainInfo& );












// INLINE FUNCTIONS NOT COVERED BY EXPLICIT TEMPLATE INSTANTIATIONS

template<uint32_t dim, template<uint32_t> class CELL>
inline void  ModelSubDomain<dim,CELL>::ScheduleForRebuild() noexcept
 {
    rebuilt_needed_ = true;
 }

template<uint32_t dim, template<uint32_t> class CELL>
inline bool  ModelSubDomain<dim,CELL>::NeedsRebuild() const noexcept
 {
    return rebuilt_needed_;
 }


template<uint32_t dim, template<uint32_t> class CELL>
inline std::string  ModelSubDomain<dim,CELL>::Name() const noexcept
 {
    return subdomain_name_;
 }

// watch out! - this name is the same as the key for the subdomain in the region,boundary,splitboundary map
template<uint32_t dim, template<uint32_t> class CELL>
inline void  ModelSubDomain<dim,CELL>::Name( const std::string& name ) noexcept
 {
    subdomain_name_ = name;
 }


template<uint32_t dim, template<uint32_t> class CELL>
inline csmp::Node<dim>*  ModelSubDomain<dim,CELL>::N( size_t nd ) const noexcept
 { assert( nd < node_vec_.size() ); return node_vec_[nd]; }


template<uint32_t dim, template<uint32_t> class CELL>
inline CELL<dim>*  ModelSubDomain<dim,CELL>::E( size_t e ) const noexcept
 { assert( e < cell_vec_.size() ); return cell_vec_[e]; }


template<uint32_t dim, template<uint32_t> class CELL>
inline const typename std::vector<CELL<dim>*>&  ModelSubDomain<dim,CELL>::CellVector() const noexcept
 { return cell_vec_; }

template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<CELL<dim>*>&  ModelSubDomain<dim,CELL>::CellVector() noexcept
 { return cell_vec_; }

template<uint32_t dim, template<uint32_t> class CELL>
inline const typename std::vector<Node<dim>*>&  ModelSubDomain<dim,CELL>::NodeVector() const noexcept
  { return node_vec_; }

template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<Node<dim>*>&  ModelSubDomain<dim,CELL>::NodeVector() noexcept
  { return node_vec_; }



template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::NodesBegin() const noexcept
 { return node_vec_.begin(); }

template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::NodesEnd() const noexcept
 { return node_vec_.end(); }

template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::PerimeterNodesBegin() const noexcept
 { return std::next( node_vec_.begin(), static_cast<long>(InteriorNodes()) ); }

template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<CELL<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::CellsBegin() const noexcept
 { return cell_vec_.begin(); }

template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<CELL<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::CellsEnd() const noexcept
 { return cell_vec_.end(); }

template<uint32_t dim, template<uint32_t> class CELL>
inline typename std::vector<CELL<dim>*>::const_iterator  ModelSubDomain<dim,CELL>::PerimeterCellsBegin() const noexcept
  { return std::next( cell_vec_.begin(), static_cast<long>(InteriorCells()) ); }
   



template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::Nodes() const noexcept
  {
//     static_assert( !std::is_same_v<CELL<dim>,InterFace<dim>>, "ModelSubDomain<dim,InterFace>::Nodes: not support for SplitBoundary" );
     return node_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::InteriorNodes() const noexcept
  {
//     static_assert( !std::is_same_v<CELL<dim>,InterFace<dim>>, "ModelSubDomain<dim,InterFace>::InteriorNodes: not support for SplitBoundary" );
     return first_bd_node_;
  }

template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::PerimeterNodes() const noexcept
  {
//     static_assert( !std::is_same_v<CELL<dim>,InterFace<dim>>, "ModelSubDomain<dim,InterFace>::PerimeterNodes: not support for SplitBoundary" );
     return node_vec_.size() - InteriorNodes();
  }

template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::Cells() const noexcept
  {
     return cell_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::InteriorCells() const noexcept
  {
     return cell_vec_.size() - bd_face_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::PerimeterCells() const noexcept
  {
     return bd_face_vec_.size();
  }

template<uint32_t dim, template<uint32_t> class CELL>
inline bool ModelSubDomain<dim,CELL>::Empty() const noexcept
  {
     return cell_vec_.empty();
  }



// PLACEMENT OF PROPERTY ASSIGNMENT

inline SUBDOMAIN_PART parseSubdomainPart( const char* subdomain )
{
    std::string ssubdomain(subdomain);

    std::transform(ssubdomain.begin(),ssubdomain.end(),ssubdomain.begin(),::toupper);

    if (  ssubdomain == "COMPLETE"  )  return COMPLETE;
    if (  ssubdomain == "INTERIOR"  )  return INTERIOR;
    if (  ssubdomain == "PERIMETER" )  return PERIMETER;
    return COMPLETE;
}


inline std::string parseSubdomainPart( SUBDOMAIN_PART ssubdomain )
 {
    if (  ssubdomain == COMPLETE  )  return std::string("COMPLETE");
    if (  ssubdomain == INTERIOR  )  return std::string("INTERIOR");
    if (  ssubdomain == PERIMETER )  return std::string("PERIMETER");
    return std::string("COMPLETE");
}


// INTERRELATIONS INTERFACE




// DOMAIN INDEXES

template<uint32_t dim, template<uint32_t> class CELL>
inline int32_t  ModelSubDomain<dim,CELL>::DomainIndex() const noexcept
{
  return domain_idx_;
}


/**

Returns a vector<double> with the ID numbers of the Elements which belong
to the Region.
*/
template<uint32_t dim, template<uint32_t> class CELL>
inline std::vector<size_t>  ModelSubDomain<dim,CELL>::MemberCellIndexes() const noexcept
 {
    std::vector<size_t> ids;
    ids.reserve( cell_vec_.size() );
    for ( const auto& it : cell_vec_ ) ids.push_back( it->Idx() );
    
    return ids;
 }



/** Renumbers nodes from 0 to n-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::RenumberNodes() const noexcept
 {
    size_t  counter(0U);

    for ( auto& it : node_vec_ ) it->Idx( counter++ );

    return counter;
 }



/** Renumbers cells from 0 to n-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
inline size_t ModelSubDomain<dim,CELL>::RenumberCells() const noexcept
 {
    size_t counter(0U);

    for( auto& it : cell_vec_ ) it->Idx(counter++);
    cell_vec_[0]->FE()->CurrentID( std::numeric_limits<size_t>::max() );

    return counter;
    
 } // end RenumberCells



/** Renumbers cells and nodes from 0 to n-1.
*/
template<uint32_t dim, template<uint32_t> class CELL>
inline void ModelSubDomain<dim,CELL>::UpdateMemberIndexes() const noexcept
 {
    RenumberNodes();
    RenumberCells();
   
 } // end UpdateRegionMemberIndexes


/**
    Uses binary_search on both ranges of the sorted node vector to find the node in question.
    returns true or false.
*/
template<uint32_t dim, template<uint32_t> class CELL>
inline bool ModelSubDomain<dim,CELL>::Contains( const Node<dim>* const nptr ) const noexcept
 {
    assert( nptr != nullptr );

    if ( binary_search( next(node_vec_.begin(), static_cast<long>(InteriorNodes())), node_vec_.end(), nptr ) )
       return true;

    if ( binary_search( node_vec_.begin(), next(node_vec_.begin(), static_cast<long>(InteriorNodes())), nptr ) )
       return true;

    return false;

 } // end


template<uint32_t dim, template<uint32_t> class CELL>
inline bool ModelSubDomain<dim,CELL>::IsPerimeterNode( const csmp::Node<dim>* const nd_ptr ) const noexcept
 {
    assert( nd_ptr != nullptr );
    return std::binary_search( PerimeterNodesBegin(), NodesEnd(), nd_ptr );
 }


template<uint32_t dim, template<uint32_t> class CELL>
inline bool  ModelSubDomain<dim,CELL>::IsPerimeterCell( const CELL<dim>* const e_ptr ) const noexcept
 {
    assert( e_ptr != nullptr );
    return std::binary_search( PerimeterCellsBegin(), CellsEnd(), e_ptr );
 }


// still used by legacy NodeCenteredFiniteVolumeTransport
template<uint32_t dim, template<uint32_t> class CELL>
inline bool  ModelSubDomain<dim,CELL>::IsPerimeterNode( size_t i ) const
 {
    if ( i >= Nodes() ) return false;
    return ( i >=  first_bd_node_ );
 }

/**
    Returns true when the queried cell index (0..n-1)
    is among those of the cells located on the boundary of the region;
    else false.
*/
template<uint32_t dim, template<uint32_t> class CELL>
inline bool  ModelSubDomain<dim,CELL>::IsPerimeterCell( size_t e ) const
 {
    if ( e >= cell_vec_.size() ) return false;
    return (e < InteriorCells()) ? false : true;
 }


} // end namespace

#endif

