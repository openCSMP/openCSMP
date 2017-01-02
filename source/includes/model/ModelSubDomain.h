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
   std::string         name;                         ///< unique name
   std::vector<size_t> interior_elmts;               ///< elements that have no face on the perimeter
   std::vector<size_t> perimeter_elmts;              ///< elements that have at least one face on perimeter
   std::vector<std::vector<int8> > perimeter_faces;  ///< local 0..faces-1 identifiers of the faces of the simplices that lie on domain boundary
   std::vector<size_t> interior_nodes;               ///< nodes within the subdomain
   std::vector<size_t> perimeter_nodes;              ///< nodes on the perimeter of the subdomain
};



/**
    Blueprint for regions and model boundaries; the latter are
    represented by lower-dimensional elements than the rest of the model.
    - Internal boundaries are named by regions they interface which each other
    - External boundaries have names corresponding to the sides of box-shaped
      models or other unique names.

      @todo (3) SKM complete switch statements for the new variable types
      @todo (3) Check for redundant inherited and non-inherited methods in subdomains (i.e. OutputVariableTo)
      @todo (1) Test whether Region properties are used correctly (A)
      @todo (2-C) Declare members as virtual if they are

*/
template<size_t dim,template<size_t> class SIMPLEX>
class ModelSubDomain : public LocalVariableStorage<dim,ModelSubDomain<dim,SIMPLEX> > {
  public:
    // simplices
    typedef SIMPLEX<dim>            Simplex;
    typedef std::vector<Simplex*>   SimplexContainer;

    // vertices
    typedef Node<dim>               Vertex;
    typedef std::vector<Vertex*>    VertexContainer;
    typedef std::vector<Node<dim>*> NodeContainer;

    // iterators
    typedef typename SimplexContainer::iterator        simplexIterator;
    typedef typename VertexContainer::iterator         vertexIterator;

    // const iterators
    typedef typename SimplexContainer::const_iterator  simplexConstIterator;
    typedef typename VertexContainer::const_iterator   vertexConstIterator;

  public:

    /// constructs incomplete subregion for later initialisation with suitable methods in subclasses
    ModelSubDomain( const std::string& subdomain_name, const PropertyDatabase<dim>& );
    ModelSubDomain( const ModelSubDomain& );
    ModelSubDomain( ModelSubDomain&& );
    virtual ~ModelSubDomain();
    ModelSubDomain<dim,SIMPLEX>&  operator=( const ModelSubDomain& );

    std::string Name() const;
    void Name( const std::string& );

    /// local variable storage interface
    virtual PLACEMENT Placement() const { return UNDEFINED; }
    virtual bool      ValidVariable( const char* variableName ) const = 0;

    virtual void Accept( Visitor<dim>& );
    void Apply( Interrelation<dim>& );

    /// connects simplices (=cells) with their equidimensional neighbors
    void EstablishNeighborConnectivity();

    /// distinguishes PERIMETER simplices that have at least one face on region boundary from INTERIOR ones; calls PartitionElementVector()
    void IdentifyPerimeter();

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

    /// reference to container of simplices that may be either of pointer to Element, Face or InterFace objects
    typename std::vector<SIMPLEX<dim>*>&  SimplexVector();
  
    /// reference to Node pointer vector
    typename std::vector<Node<dim>*>&     NodeVector();
  
    // TODO: remove this proliferation of names! - if needed put into subclasses
    typename std::vector<SIMPLEX<dim>*>&  ElementVector();
    typename std::vector<SIMPLEX<dim>*>&  FaceVector();
    typename std::vector<SIMPLEX<dim>*>&  InterFaceVector();

    // iterators
    typename std::vector<csmp::Node<dim>*>::iterator        NodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator        PerimeterNodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator        NodesEnd();
    typename std::vector<SIMPLEX<dim>*>::iterator           ElementsBegin();
    typename std::vector<SIMPLEX<dim>*>::iterator           PerimeterElementsBegin();
    typename std::vector<SIMPLEX<dim>*>::iterator           ElementsEnd();

    // const iterators
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator  PerimeterNodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesEnd() const;
    typename std::vector<SIMPLEX<dim>*>::const_iterator     ElementsBegin() const;
    typename std::vector<SIMPLEX<dim>*>::const_iterator     PerimeterElementsBegin() const;
    typename std::vector<SIMPLEX<dim>*>::const_iterator     ElementsEnd() const;

    /// returns the nodes that the region shares with the given range
    size_t SharedPerimeterNodes( typename std::vector<csmp::Node<dim>*>::const_iterator start,
                                 typename std::vector<csmp::Node<dim>*>::const_iterator end ) const;

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
    bool              Contains( const SIMPLEX<dim>* ) const;
    bool              Contains( const Node<dim>* ) const;
    bool              IsPerimeterNode( const csmp::Node<dim>* ) const;
    bool              IsPerimeterElement( const SIMPLEX<dim>* ) const;
    size_t            PerimeterFaces( size_t eid ) const;
    size_t            PerimeterFace( size_t eid, size_t face ) const;
    csmp::Node<dim>*  N( size_t n ) const;  // contained nodes
    SIMPLEX<dim>*     E( size_t n ) const;  // contained elements

    // access via object indexes( note: use with caution )
    bool              IsPerimeterNode( const size_t nidx ) const;
    bool              IsPerimeterElement( const size_t eidx ) const;

    // ----------------------------------------
    // geometry
    // ----------------------------------------

    /// returns 1) elements of how many different spatial dimensions are contained, and 2) the highest element spatial dimension in subdomain
    std::pair<int32,int32>  SpatialDimensions() const;

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
                                    double64 min_value_to_change,
                                    double64 max_value_to_change );

    /// changes the flags of the vector variable 'property' to new value; applied either in the entire subdomain or its interior or perimeter
    void ChangePropertyStatus( const char* property,
                               const std::vector<VARIABLE_FLAG>& new_status,
                               SUBDOMAIN_PART=COMPLETE );

    /// changes the flag of the vector variable 'property' to new value if the scalar variable is withing the specified range
    void ChangePropertyStatusWhere( const char* property,
                                    const std::vector<VARIABLE_FLAG>& new_status,
                                    double64 min_value_to_change,
                                    double64 max_value_to_change );

    /// changes the flag of a particular variable component to new value; applied either in the entire subdomain or its interior or perimeter
    void ChangePropertyStatus( const char* property,
                               size_t component,
                               VARIABLE_FLAG new_status,
                               SUBDOMAIN_PART=COMPLETE );

    /// changes the flag of a particular variable component to new value if the scalar variable is withing the specified range
    void ChangePropertyStatusWhere( const char* property,
                                    size_t component,
                                    VARIABLE_FLAG new_status,
                                    double64 min_value_to_change,
                                    double64 max_value_to_change );

    /// by default (i=0) returns status of scalar variable or first component of a vector or tensor variable; if i>0 flag of corresponding component is returned
    VARIABLE_FLAG  PropertyStatus( const char* variable, SUBDOMAIN_PART flag=COMPLETE , size_t i=0 ) const;

    /// min/max property values (length of vectors and eigenvalues of tensors)
    void MinMaxOf( const char* property,   double64& gmin, double64& gmax ) const;
    void MinMaxOf( const csmp::Index&,     double64& gmin, double64& gmax ) const;


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
    double64  Average( const char* property ) const;
    bool CopyGradientOfProperty_A_To_B( const char* node_prop, const char* element_prop );
    void CopyReplace( const char* from, const char* to );

    // ----------------------------------------
    // output
    // ----------------------------------------

    /// writes complete ModelSubDomain specifications in terms of unique indices as block to binary file
    void WriteDomainIndexesToBinaryFile( FILE* ) const;
    // see non-member function readDomainIndexesFromBinaryFile() to read the indices back

    void      OutputVariableToScreen( const char* prop ) const;
    void      Out() const;

    bool      Verbose();
    void      Verbose(bool verbose);

  protected:

    /// establishes interior vs. exterior simplices and nodes; returns index of first boundary element
    size_t  PartitionCellVector();

    const PropertyDatabase<dim>&                pref_;
    std::string                                 subdomain_name_; ///< passed down when region is created so that it can be referred to
    std::vector<SIMPLEX<dim>*>                  elmt_vec_;       ///< doubly sorted, interior elements first
    std::vector<std::vector<ONE_BYTE_NUMBER> >  bd_face_vec_;    ///< as in second segment of elmt_vec_
    std::vector<csmp::Node<dim>*>               node_vec_;       ///< doubly sorted, interior nodes first
    size_t                                      first_bd_node_;
    bool                                        verbose_;

  private:
    ModelSubDomain();
};


/// returns number of nodes that are shared by the two subdomains (matches by pointers)
template<size_t dim,template<size_t> class SIMPLEX>
size_t  sharedNodes( const ModelSubDomain<dim,SIMPLEX>&, const ModelSubDomain<dim,SIMPLEX>& );

/// returns number of nodes on the subdomain perimeters that are shared by the two subdomains (matches by pointers)
template<size_t dim,template<size_t> class SIMPLEX>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,SIMPLEX>&, const ModelSubDomain<dim,SIMPLEX>& );

/// reads ModelSubDomain data block written by writeDomainIndexesToBinaryFile() into the domain info structure
void readDomainIndexesFromBinaryFile( FILE*, SubDomainInfo& );


/**
    Helper functions that checks complex boundary flags of a variable that shall be assigned.
    They functions only transfer values to it if the variable has not got the specified flag.
    
    Generic version for variables that are placed on the Node, Element, Face etc.
    
    @TODO: make this part of new PDE_Intetgrator class
    
    @author SKM 10/9/2014
*/
template<size_t dim, template<size_t> class SIMPLEX, class Var>
void writeVariableIf( SIMPLEX<dim>*, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
 } // end generic specification


/// write guard for scalar variables
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             const csmp::Index& idx,
                             const ScalarVariable& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(idx) != dont_overwrite )
      ptr->Store( idx, var );
 } // end version for scalars
 

/// write guard for vector variables
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             const csmp::Index& idx,
                             const VectorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( idx, vc );
    for ( size_t i=0U; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             const csmp::Index& idx,
                             const TensorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( idx, ts );
    for ( size_t i=0U; i<dim; i++ )
      if ( ptr->Status(idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( size_t j=0U; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( idx, ts );
 } // end version for tensors



/**
    Helper functions that checks a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
    
    Generic version for variables that are placed on Element/Face/Interface integration points.
*/
template<size_t dim, template<size_t> class SIMPLEX, class Var>
void writeVariableIf( SIMPLEX<dim>*, size_t ip, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
 } // end generic specification

/// write guard for scalar variables
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t ip,
                             const csmp::Index& idx,
                             const ScalarVariable& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for scalars
 

/// write guard for vector variables
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t ip,
                             const csmp::Index& idx,
                             const VectorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( ip, idx, vc );
    for ( size_t i=0U; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( ip, idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t ip,
                             const csmp::Index& idx,
                             const TensorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( ip, idx, ts );
    for ( size_t i=0U; i<dim; i++ )
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( size_t j=0U; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( ip, idx, ts );
 } // end version for tensors



/**
    Helper functions that checks a complex varboundary flags of a variable that shall be assigned
    and only transfers values to it if the variable has not got the specified flag.
    
    Generic version for finite volume-related integration points.
*/
template<size_t dim, template<size_t> class SIMPLEX, class Var>
void writeVariableIf( SIMPLEX<dim>*, size_t sector_or_facet,
                      size_t ip, const csmp::Index&, const Var&, VARIABLE_FLAG )
 {
 } // end generic specification

/// write guard for scalar variables
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t sector_or_facet,
                             size_t ip,
                             const csmp::Index& idx,
                             const ScalarVariable& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    if ( ptr->Status(ip,idx) != dont_overwrite )
      ptr->Store( ip, idx, var );
 } // end version for scalars
 

/// write guard for vector variables
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t sector_or_facet,
                             size_t ip,
                             const csmp::Index& idx,
                             const VectorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    VectorVariable<dim> vc;
    ptr->Read( sector_or_facet, ip, idx, vc );
    for ( size_t i=0U; i<dim; i++ )
      // the component gets overwritten
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           vc.Flag(i) = var.Flag(i);
           vc(i)      = var[i];
        }
    ptr->Store( sector_or_facet, ip, idx, vc );
 } // end version for vector variables


/**
     Write guard for tensor variables
    (where only the diagonal values have flags
     so that only those rows get written where the 
     flag permits this)
*/
template<size_t dim, template<size_t> class SIMPLEX>
inline void writeVariableIf( SIMPLEX<dim>* ptr,
                             size_t sector_or_facet,
                             size_t ip,
                             const csmp::Index& idx,
                             const TensorVariable<dim>& var,
                             VARIABLE_FLAG dont_overwrite )
 {
    TensorVariable<dim> ts;
    ptr->Read( sector_or_facet, ip, idx, ts );
    for ( size_t i=0U; i<dim; i++ )
      if ( ptr->Status(ip,idx,i) != dont_overwrite ) {
           ts.Flag(i) = var.Flag(i);
           for ( size_t j=0U; j<dim; j++ )
             ts(i,j) = var(i,j);
        }
    ptr->Store( sector_or_facet, ip, idx, ts );
 } // end version for tensors



// inlined methods

template<size_t dim, template<size_t> class SIMPLEX>
inline size_t ModelSubDomain<dim,SIMPLEX>::Nodes() const
  {
     return node_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
inline size_t ModelSubDomain<dim,SIMPLEX>::InteriorNodes() const
  {
     return first_bd_node_;
  }

template<size_t dim, template<size_t> class SIMPLEX>
inline size_t ModelSubDomain<dim,SIMPLEX>::PerimeterNodes() const
  {
     return node_vec_.size() - InteriorNodes();
  }

template<size_t dim, template<size_t> class SIMPLEX>
inline size_t ModelSubDomain<dim,SIMPLEX>::Elements() const
  {
     return elmt_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
inline size_t ModelSubDomain<dim,SIMPLEX>::InteriorElements() const
  {
     return elmt_vec_.size() - bd_face_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
inline size_t ModelSubDomain<dim,SIMPLEX>::PerimeterElements() const
  {
     return bd_face_vec_.size();
  }

template<size_t dim, template<size_t> class SIMPLEX>
inline bool ModelSubDomain<dim,SIMPLEX>::Empty() const
  {
     return elmt_vec_.empty();
  }


/// returns how many faces of the target element lie on the subdomain boundary
template<size_t dim, template<size_t> class SIMPLEX>
inline size_t  ModelSubDomain<dim,SIMPLEX>::PerimeterFaces( size_t e ) const
 {
    assert( e >= InteriorElements() );
    assert( e < elmt_vec_.size() );
    return bd_face_vec_[e-InteriorElements()].size();
 }

/// returns the elements local face number of the n'th face that is on the subdomain boundary
template<size_t dim, template<size_t> class SIMPLEX>
inline size_t  ModelSubDomain<dim,SIMPLEX>::PerimeterFace( size_t e, size_t face ) const
 {
    assert( e >= InteriorElements() );
    assert( e < elmt_vec_.size() );
    assert( face < PerimeterFaces(e) );
    return static_cast<size_t>(bd_face_vec_[e-InteriorElements()][face]);
 }


template<size_t dim, template<size_t> class SIMPLEX>
inline csmp::Node<dim>*  ModelSubDomain<dim,SIMPLEX>::N( size_t nd ) const
 { assert( nd < node_vec_.size() ); return node_vec_[nd]; }


template<size_t dim, template<size_t> class SIMPLEX>
inline SIMPLEX<dim>*  ModelSubDomain<dim,SIMPLEX>::E( size_t e ) const
 { assert( e < elmt_vec_.size() ); return elmt_vec_[e]; }


template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::SimplexVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::ElementVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::FaceVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<SIMPLEX<dim>*>&  ModelSubDomain<dim,SIMPLEX>::InterFaceVector()
 { return elmt_vec_; }

template<size_t dim, template<size_t> class SIMPLEX>
typename std::vector<Node<dim>*>&  ModelSubDomain<dim,SIMPLEX>::NodeVector()
  { return node_vec_; }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::NodesBegin() const
 { return node_vec_.begin(); }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterNodesBegin() const
 { return std::next( node_vec_.begin(), InteriorNodes() ); }

template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<csmp::Node<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::NodesEnd() const
 { return node_vec_.end(); }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<SIMPLEX<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::ElementsBegin() const
 { return elmt_vec_.begin(); }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<SIMPLEX<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterElementsBegin() const
 { return std::next( elmt_vec_.begin(), InteriorElements() ); }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<SIMPLEX<dim>*>::const_iterator  ModelSubDomain<dim,SIMPLEX>::ElementsEnd() const
 { return elmt_vec_.end(); }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<csmp::Node<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::NodesBegin()
 { return node_vec_.begin(); }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<csmp::Node<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterNodesBegin()
 { return std::next( node_vec_.begin(), InteriorNodes() ); }


template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<csmp::Node<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::NodesEnd()
 { return node_vec_.end(); }

template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<SIMPLEX<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::ElementsBegin()
 { return elmt_vec_.begin(); }

template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<SIMPLEX<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::PerimeterElementsBegin()
 { return std::next( elmt_vec_.begin(), InteriorElements() ); }

template<size_t dim, template<size_t> class SIMPLEX>
inline typename std::vector<SIMPLEX<dim>*>::iterator  ModelSubDomain<dim,SIMPLEX>::ElementsEnd()
 { return elmt_vec_.end(); }


} // end namespace

#endif

