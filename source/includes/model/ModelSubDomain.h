#ifndef CSMP_MODEL_SUB_DOMAIN_H
#define CSMP_MODEL_SUB_DOMAIN_H

#include "CSMP_definitions.h"
#include "LocalVariableStorage.h"
#include "TensorVariable.h"
#include "Box.h"

namespace csmp {

/// Parts of ModelSubDomain
enum SUBDOMAIN_PART { COMPLETE,
                      INTERIOR,
                      PERIMETER };

SUBDOMAIN_PART parseSubdomainPart( const char* subdomain );
std::string    parseSubdomainPart( SUBDOMAIN_PART part );

class FiniteElementManager;
template<size_t> class PropertyDatabase;
template<size_t> class MeshManager;
template<size_t> class FiniteVolumeStencilManager;
template<size_t> class Point;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Interrelation;
template<size_t> class Visitor;

/**
    Complete index specifications of a ModelSubDomain
    used in the reconstruction of subdomains from binary files.

@author SKM
@date 14/3/2016
*/
struct SubDomainInfo {
   std::string         name;                           ///< unique name
   std::vector<size_t> interior_elmts;                 ///< elements that have no face on the perimeter
   std::vector<size_t> perimeter_elmts;                ///< elements that have at least one face on perimeter
   std::vector<std::vector<int8_t> > perimeter_faces;  ///< local 0..faces-1 identifiers of the faces of the simplices that lie on domain boundary
   std::vector<size_t> interior_nodes;                 ///< nodes within the subdomain
   std::vector<size_t> perimeter_nodes;                ///< nodes on the perimeter of the subdomain
};



/**
    Blueprint for regions and model boundaries; the latter are
    represented by lower-dimensional elements than the rest of the model.
    - Internal boundaries are named by regions they interface which each other
    - External boundaries have names corresponding to the sides of box-shaped
      models or other unique names.
*/
template<size_t dim,template<size_t> class CELL>
class ModelSubDomain {
  public:
    // any kind of finite elements; simplex or other types
    typedef CELL<dim>                                  CellType;
    // vertices
    typedef Node<dim>                                  Vertex;
    typedef std::vector<Vertex*>                       VertexContainer;
    typedef std::vector<Node<dim>*>                    NodeContainer;
    typedef typename VertexContainer::iterator         vertexIterator;
    typedef typename VertexContainer::const_iterator   vertexConstIterator;

  public:
    /// constructs incomplete subregion for later initialisation with suitable methods in subclasses
    ModelSubDomain( const std::string& subdomain_name, const PropertyDatabase<dim>& );

    ModelSubDomain( const ModelSubDomain& );
    ModelSubDomain( ModelSubDomain&& );
  
    virtual ~ModelSubDomain();
    ModelSubDomain<dim,CELL>&  operator=( const ModelSubDomain& );

    std::string Name() const;
    void Name( const std::string& );

    /// local variable storage interface
    virtual PLACEMENT Placement() const = 0;
    virtual bool      ValidVariable( const char* variableName ) const;

    virtual void Accept( Visitor<dim>& );
    void Apply( Interrelation<dim>& );

    /// distinguishes PERIMETER simplices that have at least one face on region boundary from INTERIOR ones; calls PartitionElementVector()
    void IdentifyPerimeter();
    
    /// creates node vector from element vector, using a set to achieve uniqueness
    void CreateNodePointerVector1();
    
    /// creates node vector from element vector, using a vector to achieve uniqueness via sort, unique, erase algorithms
    void CreateNodePointerVector2();

    /// sorts the node and CELL vectors split into the interior and perimeter ranges (4 sorting operations)
    void SortVectors( size_t interior_cells, size_t interior_nodes );
    
    /// assuming that a partitioned (and sorted) element vector is in place, constructs the bd_face_vec_ by checking whether neighbor elements belong to the domain or not
    void BuildPerimeterFaceVector( size_t interior_elements );

    // ----------------------------------------
    // Indexes
    // ----------------------------------------

    /// renumbers nodes in domain 0..n-1
    size_t  RenumberNodes() const;
    /// renumbers nodes in domain 0..n-1
    size_t  RenumberElements() const;
    /// renumber elements and nodes
    void    UpdateMemberIndexes() const;
    void    MemberElementIndexes( std::vector<size_t>& ) const;

    // ----------------------------------------
    // access
    // ----------------------------------------

    /// reference to container of finite element pointers to either Element, Face or InterFace objects; @note used for boolean operations
    const typename std::vector<CELL<dim>*>&  CellVector() const;
    /// do not remove!;  used for boolean operations
    typename std::vector<CELL<dim>*>&        CellVector();
  
    /// reference to Node pointer vector
    const typename std::vector<Node<dim>*>&  NodeVector() const;

    // iterators
    typename std::vector<csmp::Node<dim>*>::iterator     NodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator     NodesEnd();
    typename std::vector<csmp::Node<dim>*>::iterator     InteriorNodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator     InteriorNodesEnd();
    typename std::vector<csmp::Node<dim>*>::iterator     PerimeterNodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator     PerimeterNodesEnd();
    typename std::vector<CELL<dim>*>::iterator           ElementsBegin();
    typename std::vector<CELL<dim>*>::iterator           PerimeterElementsBegin();
    typename std::vector<CELL<dim>*>::iterator           ElementsEnd();

    // const iterators (pointer and object that is pointed to cannot be modified)
    typename std::vector<const csmp::Node<dim>* const>::const_iterator  NodesBegin() const;
    typename std::vector<const csmp::Node<dim>* const>::const_iterator  NodesEnd() const;
    typename std::vector<const csmp::Node<dim>* const>::const_iterator  InteriorNodesBegin() const;
    typename std::vector<const csmp::Node<dim>* const>::const_iterator  InteriorNodesEnd() const;
    typename std::vector<const csmp::Node<dim>* const>::const_iterator  PerimeterNodesBegin() const;
    typename std::vector<const csmp::Node<dim>* const>::const_iterator  PerimeterNodesEnd() const;
    typename std::vector<const CELL<dim>* const>::const_iterator  ElementsBegin() const;
    typename std::vector<const CELL<dim>* const>::const_iterator  ElementsEnd() const;
    typename std::vector<const CELL<dim>* const>::const_iterator  PerimeterElementsBegin() const;
    typename std::vector<const CELL<dim>* const>::const_iterator  PerimeterElementsEnd() const;
    typename std::vector<const CELL<dim>* const>::const_iterator  InteriorElementsBegin() const;
    typename std::vector<const CELL<dim>* const>::const_iterator  InteriorElementsEnd() const;

    /// returns the nodes that the region shares with the given range
    size_t SharedPerimeterNodes( typename std::vector<const csmp::Node<dim>* const>::const_iterator start,
                                 typename std::vector<const csmp::Node<dim>* const>::const_iterator end ) const;

    /// check whether subdomain conatains any elements
    bool              Empty() const;

    size_t            Nodes() const;
    size_t            InteriorNodes() const;
    size_t            PerimeterNodes() const;
    size_t            IntegrationPoints() const;
    size_t            SectorIntegrationPoints() const;
    size_t            FacetIntegrationPoints() const;
    size_t            Elements() const;
    size_t            InteriorElements() const;
    size_t            PerimeterElements() const;

    // access via objects and local order in containers
    bool              Contains( const CELL<dim>* const ) const;
    bool              Contains( const Node<dim>* const ) const;
    bool              IsPerimeterNode( const csmp::Node<dim>* const ) const;
    bool              IsPerimeterElement( const CELL<dim>* const ) const;
    /// number of faces of perimeter element #eid, that lie on subdomain surface; @attention member indexes must be are uptodate
    size_t            PerimeterFaces( size_t eid ) const;
    /// returns local face id of face #face that lies on perimeter of model subdomain
    size_t            PerimeterFace( size_t eid, size_t face ) const;
    /// pointer to node #n in subdomain; @attention node can vary from initialization to initialization
    csmp::Node<dim>*  N( size_t n ) const;
    /// pointer to element #n of model subdomain
    CELL<dim>*        E( size_t n ) const;

    // access via object indexes( note: use with caution )
    /// is the node located on the surface of the model subdomain?
    bool              IsPerimeterNode( const size_t nidx ) const;
    /// does the element have at least on face on the surface of the model subdomain?
    bool              IsPerimeterElement( const size_t eidx ) const;

    // ----------------------------------------
    // geometry
    // ----------------------------------------

    /// returns 1) elements of how many different spatial dimensions are contained, and 2) the highest element spatial dimension in subdomain
    std::pair<int32_t,int32_t>  SpatialDimensions() const;

    /// returns diagonally opposite points of bounding box
    void  MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const;

    /// output coordinates to user-defined vector variable [x1,y1,z1,...xn,yn,zn]
    void AssignNodeCoordinatesTo( const char* vector_prop );
  
    /// assigment of coordinate component 'x','y','z' to scalar variable of choice
    void AssignNodeCoordinatesTo( const char* scalar_prop, char coord );

    /// characteristics like 'length', 'area', 'volume' , 'aspect ratio', 'inner radius' are assigned to user-defined variable
    void AssignElementCharacteristicsTo( const char* characteristic, const char* var );

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
                               size_t component,
                               VARIABLE_FLAG new_status,
                               SUBDOMAIN_PART=COMPLETE );

    /// changes the flag of a particular variable component to new value if the scalar variable is withing the specified range
    void ChangePropertyStatusWhere( const char* property,
                                    size_t component,
                                    VARIABLE_FLAG new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    /// by default (i=0) returns status of scalar variable or first component of a vector or tensor variable; if i>0 flag of corresponding component is returned
    VARIABLE_FLAG  PropertyStatus( const char* variable, SUBDOMAIN_PART flag=COMPLETE , size_t i=0 ) const;

    /// min/max property values (length of vectors and eigenvalues of tensors)
    void MinMaxOf( const char* property,   double& gmin, double& gmax ) const;
    void MinMaxOf( const csmp::Index&,     double& gmin, double& gmax ) const;


    // ----------------------------------------
    // interpolation and extrapolation
    // ----------------------------------------

    void InterpolateNodeToElementProperty( const char* nprop, const char* eprop );
    void InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* eprop );
    void InterpolateIntegrationPointToElementProperty( const char* cprop, const char* eprop );
    void ExtrapolateElementToIntegrationPointProperty( const char* eprop, const char* cprop );
    void ExtrapolateElementToFacetIntegrationPointProperty( const char* eprop, const char* fipprop );
    void ExtrapolateElementToNodeProperty( const char* eprop, const char* nprop, bool by_distance=true );
    void ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop );


    // ----------------------------------------
    // calculations
    // ----------------------------------------

    /// arithmetic (number as opposed to volume weighted) average
    double  Average( const char* property ) const;
    bool CopyGradientOfProperty_A_To_B( const char* node_prop, const char* element_prop );
    void CopyReplace( const char* from, const char* to );

    // ----------------------------------------
    // output
    // ----------------------------------------

    /// writes complete ModelSubDomain specifications in terms of unique indices as block to binary file
    void WriteDomainIndexesToBinaryFile( std::fstream& ) const;
    // see non-member function readDomainIndexesFromBinaryFile() in this source file for reading the indices back

    void      OutputVariableToScreen( const char* prop ) const;
    void      Out() const;

    bool      Verbose();
    void      Verbose(bool verbose);

  protected:

    /// establishes interior vs. exterior simplices and nodes; returns index of first boundary element
    size_t  PartitionCellVector();

    const PropertyDatabase<dim>&                pref_;
    std::string                                 subdomain_name_; ///< passed down when region is created so that it can be referred to
    std::vector<CELL<dim>*>                     elmt_vec_;       ///< doubly sorted, interior elements first
    std::vector<std::vector<ONE_BYTE_NUMBER> >  bd_face_vec_;    ///< as in second segment of elmt_vec_
    std::vector<csmp::Node<dim>*>               node_vec_;       ///< doubly sorted, interior nodes first
    size_t                                      first_bd_node_;
    bool                                        verbose_;

  private:
    ModelSubDomain();
};


/// returns number of nodes that are shared by the two subdomains (matches by pointers)
template<size_t dim,template<size_t> class CELL>
size_t  sharedNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );

/// returns number of nodes on the subdomain perimeters that are shared by the two subdomains (matches by pointers)
template<size_t dim,template<size_t> class CELL>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );

/// reads ModelSubDomain data block written by writeDomainIndexesToBinaryFile() into the domain info structure
void readDomainIndexesFromBinaryFile( size_t dim, std::fstream&, SubDomainInfo& );


} // end namespace

#endif

