#ifndef CSMP_SPLIT_BOUNDARY_H
#define CSMP_SPLIT_BOUNDARY_H

#include "ModelSubDomain.h"
#include "InterFace.h"

namespace csmp {

class FiniteElementManager;
template<uint32_t> class PropertyDatabase;
template<uint32_t> class MeshManager;
template<uint32_t> class Boundary;
template<uint32_t> class Element;
template<uint32_t> class Region;
template<uint32_t> class NodeManifold;
//template<typename> class FEM_Data;

/**
    Vector of InterFace  (higher-dim) Element - face idx pairs and Element co-located with InterFace (if present).
    For each InterFace we have a pair or Element pointer - face ID pairs
        1) inner higher-dim element pointer
        2) local face number of element face that is located at interface to outer element
        3) outer higher-dim element pointer
        4) local face number of element face that is located at interface to innner element
        5) pointer to potential  lower-dimensional intervening Element
*/
template<uint32_t dim> ///
struct InterFaceParentElements : public std::vector<std::pair<std::pair<Element<dim>*,uint32_t>, std::pair<Element<dim>*,uint32_t> > > {
    // constructor
    InterFaceParentElements( const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>, std::pair<Element<dim>*,uint32_t> > >& data )
      : std::vector<std::pair<std::pair<Element<dim>*,uint32_t>, std::pair<Element<dim>*,uint32_t> > >(data) {}
      
    typedef typename InterFaceParentElements<dim>::const_iterator ifaceIterator;
    // data members
    Element<dim>* InnerElement( ifaceIterator it ) const { return (*it).first.first; }
    Element<dim>* OuterElement( ifaceIterator it ) const { return (*it).second.first; }
    uint32_t InnerFaceID( ifaceIterator it ) const { return (*it).first.second; }
    uint32_t OuterFaceID( ifaceIterator it ) const { return (*it).second.second; }
};



/** 

@brief Representation of a model boundary across which the elements are node-matching, but have
unique nodes on either side; as many nodes as come together.

@author P. Lang
@author S.K. Matthai
@date 2011

SplitBoundaries are collections of InterFaces (interior material boundaries) singled out for specific
computations.

@section conventions Conventions

The base class node pointer vector contains nodes of the inner side of the split boundary.
If one wants to access neighboring nodes, this would be done via looping over the
boundaries interfaces and use InterFace::N(i,j).

@section implementation Implementation

Split boundaries  split the nodes of the provided Boundary plus those that belong to
another Boundary as well. Hence, use of Boundaries is required where all external
and internal boundaries of the domain are represented as such (exist as a csmp::Boundary).

Opposing parent elements on either side of the csmp::Boundary are removed from each others neighbor list, so they are
from the now duplicated nodes parent list on the opposing side.

@attention Before andy SplitBoundaries may be built all Boundaries have to be set up.
@attention 2D models in 3D space will not allow for proper split boundary creation (unit normal of line element issue)

@section motivation Motivation

To address and carry out computations at material interfaces or model
boundaries.

To avoid conditional processing (if at boundary statements etc.) for all elements of a model, by treating
boundary specific calculations separately.

@section consequences Consequences

SplitBoundary objects permit the implementation of jump discontinuities in continuum models, see Tran et al. (2020, AWR)

@note SplitBoundary has no BOX_BOUNDARY flag because it can only have the value INTERNAL anyway.

*/
template<uint32_t dim>
class SplitBoundary final : public ModelSubDomain<dim,InterFace>,
                            public LocalVariableStorage<dim, SplitBoundary>
 {
  public:
    SplitBoundary() = delete;
    
    /// constructor of split boundary with given name from set of juxtaposed elements; prompts MeshManager to create elements
    SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>&, 
                   MeshManager<dim>&, const InterFaceParentElements<dim>& );
                   
    SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>& );
    
    SplitBoundary( const SplitBoundary& );
    SplitBoundary( SplitBoundary&& );
    
    /// RECONSTRUCTOR  of split boundaries from csmp native file format (call only prior to deleting anythin from colonies)
    SplitBoundary( const PropertyDatabase<dim>&,
                   MeshManager<dim>&,
                   const SubDomainInfo& );  ///< contains correctly partitioned vectors and boundary faces
    
    ~SplitBoundary() = default;
    SplitBoundary<dim>&  operator=( const SplitBoundary& );

    /// applies visitor to split boundary
    virtual void Accept( Visitor<dim>& ) override;
    
    /// methods required for the LocalVariableStorage
    virtual PLACEMENT Placement() const noexcept override final { return SPLIT_BOUNDARY; }
    virtual bool ValidVariable( const char* variableName ) const noexcept override final;


    // --------------------------------------------------
    // building of SplitBoundaries and their modification
    // --------------------------------------------------
    
    /// forms a split boundary from existing interfaces assuming that these are numbered n=elmts+faces .. interfaces-1
    size_t AccumulateByNumber( MeshManager<dim>&, std::vector<size_t>& cell_ids );

    /// creating from supplied vector of interfaces
    bool CreateFrom( const typename std::vector<InterFace<dim>*>::const_iterator ifacesBegin,
                     const typename std::vector<InterFace<dim>*>::const_iterator ifacesEnd );


    // ----------------------------------------
    // user interface
    // ----------------------------------------
    
    /// renumbers nodes in SplitBoundary domain 0..n-1 starting with the inside
    size_t RenumberNodes() const noexcept override final;

    /// generates a vector of NodeManifold pointers that can be iterated over
    std::vector<NodeManifold<dim>*>  NodeManifolds() const;
    
    /// returns sorted node pointer vector partitioned into an interior and perimeter range; perimeter starts at size_t number of nodes in returned pair
    std::pair<std::vector<Node<dim>*>,size_t>  InsideNodes() const;

    /// returns sorted node pointer vector partitioned into an interior and perimeter range; perimeter starts at size_t number of nodes in returned pair
    std::pair<std::vector<Node<dim>*>,size_t>  OutsideNodes() const;

    /// NEW: input constat node variable value on a specific side of the split boundary (options INSIDE or OUTSIDE)
    template<class Var> requires CsmpVariable<dim, Var>
    void InputNodePropertyValue( const char* input_node_prop, const Var&, SUBDOMAIN_PART, INTERFACE_SIDE );

    /// input node variable values on a specific side of the split boundary (options INSIDE or OUTSIDE)
    template<class Var> requires CsmpVariable<dim, Var>
    void InputPropertyValue( const char* input_node_prop, const Var&, SUBDOMAIN_PART, INTERFACE_SIDE );

    /// universal, safe alternative to Store()  to assign SplitBoundary properties (1 value per SplitBoundary()
    template<typename Var> requires CsmpVariable<dim, Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );

    /// same as InputPropertyValue, but with overwrite protection for variable components that have a different flag than 'do_not_overwrite'
    template<typename Var> requires CsmpVariable<dim, Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );

    /// changes the flag of the scalar node variable 'property' to the new value; applied either to the INSIDE or OUTSIDE of the entire splitboundary or its interior or perimeter
    void ChangeNodePropertyStatus( const char* property,
                                   VARIABLE_FLAG new_status_of_scalar,
                                   SUBDOMAIN_PART, INTERFACE_SIDE );

    /// changes the flag of the scalar node variable 'property' to the new value; applied either to the INSIDE or OUTSIDE of the entire splitboundary but only where the variable value is in the target range
    void ChangeNodePropertyStatusWhere( const char*, VARIABLE_FLAG, INTERFACE_SIDE, double, double );

    /// returns facing direction of split boundary relative to adjacent region
    INTERFACE_SIDE  RegionLocation( const Region<dim>& );

    /// output length(2D) or area(3d) of the split boundary=lower dimensional region; middle refers to bisector if nodes are displaced
    double  Area( INTERFACE_SIDE=INSIDE ) const;

    /// output length(of a split boundary that consists of line elements; middle refers to bisector if nodes are displaced
    double  Length( INTERFACE_SIDE=INSIDE ) const;

    /// outputs length of perimeter curve of a 3D split boundary; no meaning in 1 or 2D models
    double  Perimeter( INTERFACE_SIDE=INSIDE ) const;
  
    /// integrates the property over the boundary line or surface
    double  SurfaceIntegral( const PropertyDatabase<dim>&, const char* property, INTERFACE_SIDE=INSIDE ) const;

    /// returns 1) interfaces of how many different spatial dimensions are contained, and 2) the highest interface spatial dimension in subdomain
    std::pair<int32_t,int32_t>  InterFaceSpatialDimensions() const;
    
    /// reports box-boundary flag equivalent which is always INTERNAL because SplitBoundary objects can only exist on the interior of a model
    BOX_BOUNDARY AtBoundary() const { return INTERNAL; }

	/// Method that applies half of displacement to both sides of interface (dotted by unitNormal)
    void PullApartSplitBoundary( double dist );

    /// writes all contained data on the screen
    void Out() const;


    // base class methods that are not available for SplitBoundary because it has no node vector

    /// different version than in the base class that does not attempt to write nodes
    void WriteSplitBoundaryIndexesToBinaryFile( std::fstream& ) const;
 
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesBegin() const = delete;
    typename std::vector<csmp::Node<dim>*>::const_iterator  PerimeterNodesBegin() const = delete;
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesEnd() const = delete;

    // Caution: these methods will not compile only if one attempts to use them on a SplitBoundary
    
    /// all these methods create issues if called on SplitBoundary
    size_t            Nodes() const = delete;
    size_t            InteriorNodes() const = delete;
    size_t            PerimeterNodes() const = delete;
    bool              IsPerimeterNode( const csmp::Node<dim>* const ) const = delete;
    csmp::Node<dim>*  N( size_t n ) const = delete;
    bool              IsPerimeterNode( const size_t nidx ) const = delete;
    void              AssignNodeCoordinatesTo( const char* vector_prop ) = delete;
    void              AssignNodeCoordinatesTo( const char* scalar_prop, char coord ) = delete;
    void              NodeAttributesToCSV() = delete;

    void ChangePropertyStatus( const char*, VARIABLE_FLAG, SUBDOMAIN_PART ) = delete;
    void ChangePropertyStatusWhere( const char*, VARIABLE_FLAG, double, double ) = delete;
    void ChangePropertyStatus( const char*, const std::vector<VARIABLE_FLAG>& ,SUBDOMAIN_PART ) = delete;
    void ChangePropertyStatusWhere( const char*, const std::vector<VARIABLE_FLAG>&, double, double ) = delete;
    void ChangePropertyStatus( const char*, uint32_t, VARIABLE_FLAG, SUBDOMAIN_PART ) = delete;
    void ChangePropertyStatusWhere( const char*, uint32_t, VARIABLE_FLAG, double, double ) = delete;

    void InterpolateNodeToCellProperty( const char*, const char* ) = delete;
    void InterpolateNodeToIntegrationPointProperty( const char*, const char* ) = delete;

  protected:
  
    /// return physical variable count at given integration points
    IntegrationPointVariables  InterFaceIntegrationPointVariables() const;
    LocalVariables             InterFaceVariables() const;
};




// function prototypes

template<uint32_t dim, typename Var>
void inputNodePropertyValue( SplitBoundary<dim>&,
                             const PropertyDatabase<dim>&,
                             const char* input_prop,
                             const Var&,
                             size_t side );
   
  /// matches member function to write domain indices
  void readSplitBoundaryIndexesFromBinaryFile( std::fstream&, SubDomainInfo& info );

} // end csmp

#endif
