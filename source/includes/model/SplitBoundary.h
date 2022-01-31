#ifndef CSMP_SPLIT_BOUNDARY_H
#define CSMP_SPLIT_BOUNDARY_H

#include "ModelSubDomain.h"
#include "InterFace.h"

namespace csmp {

class FiniteElementManager;
template<size_t> class PropertyDatabase;
template<size_t> class MeshManager;
template<size_t> class Boundary;
template<size_t> class Element;
template<size_t> class Region;
template<typename> class FEM_Data;

/**
    Set of InterFace  (higher-dim) Element - face idx pairs and Element co-located with InterFace (if present).
    For each InterFace we have a pair or Element pointer - face ID pairs
        1) inner higher-dim element pointer
        2) local face number of element face that is located at interface to outer element
        3) outer higher-dim element pointer
        4) local face number of element face that is located at interface to innner element
        5) pointer to potential  lower-dimensional intervening Element
*/
template<size_t dim> ///
struct InterFaceSet : public std::set<std::pair<std::pair<Element<dim>*,size_t>, std::pair<Element<dim>*,size_t> > > {
    // constructor
    InterFaceSet( const std::set<std::pair<std::pair<Element<dim>*,size_t>,
                  std::pair<Element<dim>*,size_t> > >& set )
      : std::set<std::pair<std::pair<Element<dim>*,size_t>, std::pair<Element<dim>*,size_t> > >(set) {}
      
    typedef typename InterFaceSet<dim>::const_iterator ifaceIterator;
    // data members
    Element<dim>* InnerElement( ifaceIterator it ) const { return (*it).first.first; }
    Element<dim>* OuterElement( ifaceIterator it ) const { return (*it).second.first; }
    size_t InnerFaceID( ifaceIterator it ) const { return (*it).first.second; }
    size_t OuterFaceID( ifaceIterator it ) const { return (*it).second.second; }
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
template<size_t dim>
class SplitBoundary : public ModelSubDomain<dim,InterFace>,
                      public LocalVariableStorage<dim, SplitBoundary>
 {
  public:
    SplitBoundary() = delete;
    /// constroctor of split boundary with given name from set of juxtaposed elements; prompts MeshManager to create elements
    SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>&, 
                   const FiniteElementManager&, MeshManager<dim>&, const InterFaceSet<dim>& );
                   
    SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>& );
    
    SplitBoundary( const SplitBoundary& );
    SplitBoundary( SplitBoundary&& );
    
    /// RECONSTRUCTOR  of split boundaries from csmp native file format (call only prior to deleting anythin from colonies)
    SplitBoundary( const PropertyDatabase<dim>&,
                   MeshManager<dim>&,
                   const SubDomainInfo& );  ///< contains correctly partitioned vectors and boundary faces
    
    virtual ~SplitBoundary();
    SplitBoundary<dim>&  operator=( const SplitBoundary& );

    /// applies visitor to split boundary
    virtual void Accept( Visitor<dim>& );
    
    /// methods required for the LocalVariableStorage
    virtual PLACEMENT Placement() const { return SPLIT_BOUNDARY; }
    virtual bool ValidVariable( const char* variableName ) const;


    // --------------------------------------------------
    // building of SplitBoundaries and their modification
    // --------------------------------------------------

    /// creates split boundary from boundary assuming that nodes have already been duplicated etc.
    bool CreateFrom( MeshManager<dim>&, Boundary<dim>& );
  
    // ----------------------------------------
    // user interface
    // ----------------------------------------

    /// for the assignment of properties that are unique to the instance of this subclass
    template<typename Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );

    /// as InputPropertyValue, but with overwrite protection for variable components that have the flag 'do_not_overwrite'
    template<typename Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );

    /// input node variable values on specific side of split boundary
    template<class Var>
    void InputNodePropertyValue( const char* input_prop, const Var&, SUBDOMAIN_PART, INTERFACE_SIDE=INSIDE );

    // returns location of split boundary relative to adjacent region
    INTERFACE_SIDE  RegionLocation( const Region<dim>& );

    /// output length(2D) or area(3d) of the split boundary=lower dimensional region; middle refers to bisector if nodes are displaced
    double  Area( INTERFACE_SIDE=MIDDLE ) const;
  
    /// outputs length of perimeter curve of a 3D split boundary; no meaning in 1 or 2D models
    double  Perimeter( INTERFACE_SIDE=MIDDLE ) const;
  
    /// integrates the property over the boundary line or surface
    double  SurfaceIntegral( const PropertyDatabase<dim>&, const char* property, INTERFACE_SIDE=INSIDE  ) const;

    /// returns 1) interfaces of how many different spatial dimensions are contained, and 2) the highest interface spatial dimension in subdomain
    std::pair<int32_t,int32_t>  InterFaceSpatialDimensions() const;
    
    /// reports box-boundary flag equivalent which is always INTERNAL because SplitBoundary objects can only exist on the interior of a model
    BOX_BOUNDARY AtBoundary() const { return INTERNAL; }

    /// writes all contained data on the screen
    void Out() const;

  protected:

    /// reestablishes the pointers to the nodes associated with the stored elements
    void CreateNodePointerVector();

    /// return physical variable count at given integration points
    IntegrationPointVariables  InterFaceIntegrationPointVariables() const;
    LocalVariables             InterFaceVariables() const;
};




// function prototypes

template<size_t dim, typename Var>
void inputNodePropertyValue( SplitBoundary<dim>&,
                             const PropertyDatabase<dim>&,
                             const char* input_prop,
                             const Var&,
                             size_t side );

} // end csmp

#endif
