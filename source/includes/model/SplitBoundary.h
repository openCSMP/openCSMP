#ifndef CSMP_SPLIT_BOUNDARY_H
#define CSMP_SPLIT_BOUNDARY_H

#include "ModelSubDomain.h"
#include "InterFace.h"

namespace csmp {

class FiniteElementManager;
template<size_t> class PropertyDatabase;
template<size_t> class InterFace;
template<size_t> class MeshManager;
template<size_t> class Boundary;
template<size_t> class Element;
template<size_t> class InterFace;
template<size_t> class Region;
template<typename> class FEM_Data;


template<size_t dim> /// set of (higher-dim) Element - (inter)face idx pairs + Element co-located with InterFace (if present)
struct InterFaceSet : public std::set<std::pair<std::pair<Element<dim>*,size_t>, std::pair<Element<dim>*,size_t> > > {
  // constructor
  InterFaceSet( const std::set<std::pair<std::pair<Element<dim>*,size_t>, std::pair<Element<dim>*,size_t> > >& set )
    : std::set<std::pair<std::pair<Element<dim>*,size_t>, std::pair<Element<dim>*,size_t> > >(set) {}
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

Split boundaries will split the interior nodes of the provided Boundary plus those that belong to
another Boundary as well. Hence, use of Boundaries is required where all external
and internal boundaries of the domain are represented as such (exist as a csmp::Boundary).

Opposing parent elements on either side of the csmp::Boundary are removed from each others neighbor list, so they are
from the now duplicated nodes parent list on the opposing side.


@attention Before andy SplitBoundaries may be built all Boundaries have to be set up.
@attention 2D models in 3D space will not allow for proper split boundary creation (unit normal of line element issue)

@section motivation Motivation

To address and carry out computations at material interfaces or model
boundaries.

To avoid conditional processing of all elements of a model, by treating
boundary specific calculations separately.

@section consequences Consequences

Through the use of boundaries essential conditions and coupling terms
arising at material interfaces can be computed with greater ease and
efficiency.

*/
template<size_t dim>
class SplitBoundary : public ModelSubDomain<dim,InterFace> {
  public:
    SplitBoundary() = delete;
    /// constroctor of split boundary with given name from set of juxtaposed elements; prompts MeshManager to create elements
    SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>&, const FiniteElementManager&, MeshManager<dim>&, const InterFaceSet<dim>& );
    SplitBoundary( std::string splitboundaryname, const PropertyDatabase<dim>& );
    SplitBoundary( const SplitBoundary& );
    SplitBoundary( SplitBoundary&& );
    /// re-constructor of split boundaries from csmp native file format
    SplitBoundary( const PropertyDatabase<dim>&,
                   MeshManager<dim>&,       ///< not constant since write access is granted to boundary
                   const SubDomainInfo& );  ///< contains correctly partitioned vectors and boundary faces
    
    virtual ~SplitBoundary();
    SplitBoundary<dim>&  operator=( const SplitBoundary& );

    /// applies visitor to split boundary
    virtual void Accept( Visitor<dim>& );
    
    virtual PLACEMENT Placement() const { return SPLIT_BOUNDARY; }
    virtual bool ValidVariable( const char* variableName ) const;


    // --------------------------------------------------
    // building of SplitBoundaries and their modification
    // --------------------------------------------------

    /// creates split boundary from boundary
    bool CreateFrom( Model<dim>&, Boundary<dim>& );
  
    /// reestablishes the pointers to the nodes associated with the stored elements
    void CreateNodePointerVector();

    // ----------------------------------------
    // user interface
    // ----------------------------------------

    /// input node variable values on specific side of split boundary
    template<class Var>
    void InputNodePropertyValue( const char* input_prop, const Var&, SUBDOMAIN_PART, INTERFACE_SIDE=INSIDE );

    // returns location of split boundary relative to adjacent region
    INTERFACE_SIDE  RegionLocation( const Region<dim>& );

    /// output length(2D) or area(3d) of the split boundary=lower dimensional region; middle refers to bisector if nodes are displaced
    double64  Area( INTERFACE_SIDE=MIDDLE ) const;
  
    /// outputs length of perimeter curve of a 3D split boundary; no meaning in 1 or 2D models
    double64  Perimeter( INTERFACE_SIDE=MIDDLE ) const;
  
    /// integrates the property over the boundary line or surface
    double64  SurfaceIntegral( const PropertyDatabase<dim>&, const char* property, INTERFACE_SIDE=INSIDE  ) const;

    /// returns 1) interfaces of how many different spatial dimensions are contained, and 2) the highest interface spatial dimension in subdomain
    std::pair<int32,int32>  InterFaceSpatialDimensions() const;

    /// writes all contained data on the screen
    void Out() const;

    // -----------------------------------------------
    // Binary input/output
    // -----------------------------------------------

    /// output to binary file
    bool Out( std::fstream& ) const;
  
    /// initializes split boundary from binary file
    bool In( MeshManager<dim>& ,
             const FiniteElementManager& ,
             const Region<dim>&, std::fstream& );
             
             
  protected:
      /// creates split boundary from supplied vectors (used in binary IO - supplied ids used from Model region)
      bool CreateFrom( MeshManager<dim>& ,
                       const FiniteElementManager& ,
                       const std::map<size_t,csmp::Element<dim>*>&        elementIdPtr,
                       const std::vector<CSMP_FEM_TYPE>&                  interfaceTypes,
                       const std::vector<std::vector<size_t> >&           interfaceParents,
                       const std::vector<std::vector<std::pair<size_t,size_t> > >&  interfaceParentNodes );

      /// establishes connectivity and initializes LVS
      void Initialize( bool updateNeighborConnectivity = true, bool updateIndexes = true );

  protected:
    /// return physical variable count at given integration points
    IntegrationPointVariables InterFaceIntegrationPointVariables() const;
    LocalVariables InterFaceVariables() const;

    // TODO: move to the SplitBoundaryInterface
    void Split( Model<dim>& , Boundary<dim>& );
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
