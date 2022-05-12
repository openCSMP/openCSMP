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
   std::string         name;                           ///< unique name
   std::vector<uint32_t> interior_elmts;                 ///< cells that have no face on the perimeter
   std::vector<uint32_t> perimeter_elmts;                ///< cells that have at least one face on perimeter
   std::vector<std::vector<int8_t> > perimeter_faces;  ///< local 0..faces-1 identifiers of the faces of the simplices that lie on domain boundary
   std::vector<uint32_t> interior_nodes;                 ///< nodes within the subdomain
   std::vector<uint32_t> perimeter_nodes;                ///< nodes on the perimeter of the subdomain
};



/**
    Blueprint for domains and model boundaries; the latter are
    represented by lower-dimensional cells than the rest of the model.
    - Internal boundaries are named by domains they interface which each other
    - External boundaries have names corresponding to the sides of box-shaped
      models or other unique names.
*/
template<uint32_t dim,template<uint32_t> class CELL>
class ModelSubDomain {
  public:
    typedef CELL<dim> CellType; // used by SteadyStateDiffusor and others
    
    /// constructs incomplete subdomain for later initialisation with suitable methods in subclasses
    ModelSubDomain( const std::string& subdomain_name, const PropertyDatabase<dim>& );
    ModelSubDomain( const ModelSubDomain& );
    ModelSubDomain( ModelSubDomain&& );
  
    virtual ~ModelSubDomain();
    
    ModelSubDomain<dim,CELL>&  operator=( const ModelSubDomain& );
    ModelSubDomain<dim,CELL>&  operator=( ModelSubDomain&& );

    std::string Name() const;
    void Name( const std::string& );

    /// local variable storage interface
    virtual PLACEMENT Placement() const = 0;
    virtual bool      ValidVariable( const char* variableName ) const;

    /// support of the Visitor pattern
    virtual void Accept( Visitor<dim>& );
    
    /// modification via Operand-based relations between discretised properties that only modify a single result variable
    void Apply( Interrelation<dim>& );
    
    /// deletes nullptr cells, rebuilds node vector, sorts everything and re-establishes the perimeter face vectors after modifications of cells
    void RebuildSubDomainAfterChangeOfCellVector();
    
    /// rebuilds subdomain on the basis of the cells that will be selected according to the supplied property constraints
    void UpdateCellMembershipApplyingConstraints( typename std::vector<CELL<dim>*>::const_iterator master_domain_start,
                                                  typename std::vector<CELL<dim>*>::const_iterator master_domain_end,
                                                  const PropertyConstraints& );

    /// distinguishes PERIMETER simplices that have at least one face on domain boundary from INTERIOR ones; calls PartitionElementVector()
    void IdentifyPerimeter();
    
    /// creates node pointer vector from cell vector, using a vector to achieve uniqueness via sort, unique, erase algorithms
    void CreateNodePointerVector();

    /// creates node pointer vector from the shared face nodes of the supplied range of contacting cells
    void CreateNodePointerVector( std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& contacting_cells );

    /// sorts the node and CELL vectors split into the interior and perimeter ranges (4 sorting operations)
    void SortVectors( size_t interior_cells, size_t interior_nodes );
    
    /// assuming that a partitioned (and sorted) cell vector is in place, constructs the bd_face_vec_ by checking whether neighbor cells belong to the domain or not
    void BuildPerimeterFaceVector( size_t interior_cells );
    
    /// removes any cells or node pointers that were set to zero elsewhere; returns number of cells removed
    size_t RemoveNullPointerCells();
    
    /// for rebuilding subdomains when nodes or cells changed:  flag up for a rebuild using RebuildSubDomainAfterChangeOfCellVector
    void ScheduleForRebuilt();
    bool NeedsRebuilt() const;

    // ----------------------------------------
    // Indexes
    // ----------------------------------------

    /// reference counting-based unique domain identifier  (0..n-1)
    size_t DomainIndex() const;

    /// renumbers nodes in domain 0..n-1
    size_t  RenumberNodes() const;
    
    /// renumbers cells in domain 0..n-1
    size_t  RenumberCells() const;
    
    /// renumber cells and nodes consecutively from 0..n-1
    void    UpdateMemberIndexes() const;
    
    /// output current indices to vector
    void    MemberCellIndexes( std::vector<size_t>& ) const;
    
    /// setting all cell indices to a specific value
    void    SetCellIndexes( size_t new_idx );

    // ----------------------------------------
    // access
    // ----------------------------------------

    /// reference to container of finite cell pointers to either Element, Face or InterFace objects; @note used for boolean operations
    const typename std::vector<CELL<dim>*>&  CellVector() const;
    
    /// do not remove!;  used for boolean operations
    typename std::vector<CELL<dim>*>&        CellVector();
  
    /// reference to Node pointer vector
    const typename std::vector<Node<dim>*>&  NodeVector() const;

    /// const iterators (cell and node pointers cannot be modified but the nodes and cells can!)
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator  PerimeterNodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesEnd() const;
    typename std::vector<CELL<dim>*>::const_iterator        CellsBegin() const;
    typename std::vector<CELL<dim>*>::const_iterator        PerimeterCellsBegin() const;
    typename std::vector<CELL<dim>*>::const_iterator        CellsEnd() const;

    /// returns how many of the nodes in the supplied iterator range also form part of the current subdomain's perimeter
    size_t SharedPerimeterNodes( typename std::vector<csmp::Node<dim>*>::const_iterator start,
                                 typename std::vector<csmp::Node<dim>*>::const_iterator end ) const;
                                 
    // see also the non-member functions below

    /// check whether subdomain conatains any cells
    bool              Empty() const;
    
    size_t            Nodes() const;
    size_t            InteriorNodes() const;
    size_t            PerimeterNodes() const;
    size_t            IntegrationPoints() const;
    size_t            SectorIntegrationPoints() const;
    size_t            FacetIntegrationPoints() const;
    size_t            Cells() const;
    size_t            InteriorCells() const;
    size_t            PerimeterCells() const;
    
    bool              Contains( const CELL<dim>* const ) const;
    bool              Contains( const Node<dim>* const ) const;
    
    bool              IsPerimeterNode( const csmp::Node<dim>* const ) const;
    bool              IsPerimeterCell( const CELL<dim>* const ) const;
    
    /// for looping over the perimeter faces of perimeter cell with #eid, method is used to travel across subdomain surface / outline
    uint32_t          PerimeterFaces( size_t eid ) const;
    /// returns local face id of face #face that lies on perimeter of model subdomain
    uint32_t          PerimeterFace( size_t eid, uint32_t face ) const;

    /// pointer to node #n in subdomain; @attention node can vary from initialization to initialization
    csmp::Node<dim>*  N( size_t n ) const;
    /// pointer to cell #n of model subdomain
    CELL<dim>*        E( size_t n ) const;

    /// is the node located on the surface of the model subdomain?
    bool              IsPerimeterNode( const size_t nidx ) const;
    /// does the cell have at least on face on the surface of the model subdomain?
    bool              IsPerimeterCell( const size_t eidx ) const;

    // ----------------------------------------
    // geometry
    // ----------------------------------------
    
    /// checks whether all cells within the subdomain are interconnected (if the domain has multi-dimensional cells this is never the case)
    bool IsContiguous() const;

    /// if the model subdomain consists of a single cell shape (line, surface, volume) this method returns true revealing the type
    std::pair<CELL_SHAPE,bool>  SingleCellShapeDomain() const;

    /// returns 1) cells of how many different spatial dimensions are contained, and 2) the highest cell spatial dimension in subdomain
    std::pair<int32_t,int32_t>  SpatialDimensions() const;

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
    void InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );

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

    // ----------------------------------------
    // output
    // ----------------------------------------

    /// writes complete ModelSubDomain specifications in terms of unique indices as block to binary file
    void WriteDomainIndexesToBinaryFile( std::fstream& ) const;
    // see non-member function readDomainIndexesFromBinaryFile() in this source file for reading the indices back

    void      OutputVariableToScreen( const char* prop ) const;
    void      Out() const;

  protected:

    /// establishes interior vs. exterior simplices and nodes; returns index of first boundary cell
    size_t  PartitionCellVector();

    const PropertyDatabase<dim>&        pref_;
    std::string                         subdomain_name_;         ///< passed down when domain is created so that it can be referred to
    std::vector<CELL<dim>*>             cell_vec_;               ///< doubly sorted, interior cells first
    std::vector<std::vector<uint32_t> > bd_face_vec_;            ///< as in second segment of cell_vec_
    std::vector<csmp::Node<dim>*>       node_vec_;               ///< doubly sorted, interior nodes first
    size_t                              first_bd_node_ = std::numeric_limits<uint32_t>::max(); ///< begin of the perimeter nodes
    inline static int32_t               domain_count_ = 0;       ///<  reference-counting to get unique identifier for subdomains
    int32_t                             domain_idx_;             ///< created during construction from domain_count_
    bool                                rebuilt_needed_ = false; ///< parameter set when mesh gets modified by MeshManager so that update can be prompted
    static constexpr bool               verbose_ = false;

  private:
    ModelSubDomain();
};

/// distinguishes between Region, Boundary and SplitBoundary on the basis of the name string
PLACEMENT modelSubdomainType( const std::string& subdomain_name );

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

/// reads ModelSubDomain data block written by writeDomainIndexesToBinaryFile() into the domain info structure
void readDomainIndexesFromBinaryFile( uint32_t dim, std::fstream&, SubDomainInfo& );


} // end namespace

#endif

