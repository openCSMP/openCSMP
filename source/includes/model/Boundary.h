#ifndef CSMP_BOUNDARY_H
#define CSMP_BOUNDARY_H

#include "ModelSubDomain.h"
#include "Face.h"

namespace csmp {

class FiniteElementManager;
class FaceConstructionData;
template<size_t> class PropertyDatabase;
template<size_t> class MeshManager;
template<typename> class FEM_Data;
template<size_t> class Point;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Region;

/**

@brief Boundary is lower-dimensional subclass and specialisation of ModelSubDomain 
which consist of Face objects.

Boundaries are collections of Faces covering interior and exterior model boundaries,
and are made accessible through BoundaryInterface of the Model class.

@section motivation Motivation

To apply conditions, carry out computations or monitor processes
at material interfaces or model side boundaries.

To distinguish surface elements by turning them into Face objects 
for specific processing, and treating
boundary-specific calculations separately.

@section consequences Consequences

Through the use of boundaries essential conditions and coupling terms
arising at material interfaces can be computed with greater ease and
efficiency.

Boundary adds an extra variable placement BOUNDARY, which permits
targeting of variables like 'basal heat flow' etc.

@author S.K. Matthai
@author P. Lang
@date 2010

*/
template<size_t dim>
class Boundary : public ModelSubDomain<dim, Face>,
                 public LocalVariableStorage<dim, Boundary>
 {
  public:
    // ------------------------------------------------
    // construction of boundaries from scratch
    // ------------------------------------------------
    Boundary() = delete; ///< there is no sensible default contruction

    /// constructor: creating boundary from supplied vector of faces
    Boundary( const std::string& boundary_name,
              const PropertyDatabase<dim>&,
              typename std::vector<Face<dim>*>::iterator facesBegin,
              typename std::vector<Face<dim>*>::iterator facesEnd,
              BOX_BOUNDARY );

    /// create empty boundary with appropriately resized property storage
    Boundary( std::string boundaryname, const PropertyDatabase<dim>&, BOX_BOUNDARY flag );

    Boundary( const Boundary& );
    Boundary( Boundary&& );

    /// gets MeshManager to delete the boundary including Faces
    virtual ~Boundary();

    Boundary<dim>&  operator=( const Boundary& );

    /// methods required for the LocalVariableStorage
    virtual PLACEMENT Placement() const { return BOUNDARY; }
    virtual bool      ValidVariable( const char* variableName ) const;


    // ------------------------------------------------
    // reconstruction of boundaries that existed before
    // ------------------------------------------------

    /// RECONSTRUCTOR of boundary from index data stored in SubDomainInfo
    Boundary( const PropertyDatabase<dim>&,
              const MeshManager<dim>&,       
              const SubDomainInfo&,    ///< contains correctly partitioned vectors and boundary faces
              BOX_BOUNDARY = IRREGULAR );

    /// TODO: still needed? - re-constructor of boundary from index data stored in SubDomainInfo and faces from the MeshManager
    /*
    Boundary( const PropertyDatabase<dim>&,
              const size_t&,
              const std::deque<Node<dim>*>&,
              const std::deque<Face<dim>*>&,
              const SubDomainInfo&,
              BOX_BOUNDARY = IRREGULAR );
    */

    /// Visitors
    virtual void Accept( Visitor<dim>& );


    // -----------------------------------------------
    //  input/output to binary file
    // -----------------------------------------------

    /// for the assignment of properties that are unique to the instance of this subclass
    template<typename Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );

    /// as InputPropertyValue, but with overwrite protection for variable components that have the flag 'do_not_overwrite'
    template<typename Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );

    /// assigning the variable values from the FEM_DATA container to the corresponding property of the Boundary
    template<class Var>
    void InputVariableFrom( const char* property, const FEM_Data<Var>& );

    template<class Var>
    void OutputVariableTo( const char* property, FEM_Data<Var>& ) const;


    // ----------------------------------------
    //  boundary initialisation & modification
    // ----------------------------------------

    /// reestablish nodes based on element container
    void CreateNodePointerVector();
  
    /// if parts of this boundary coincide with subsetRegion, this part will pasted to subsetBoundaryToForm
    bool Divide( const Region<dim>& subsetRegion,
                 Boundary<dim>& subsetBoundaryToForm );

    /// Removes faces from this boundary and assigns to subsetBoundaryToForm
    bool Divide( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                 const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                 Boundary<dim>& subsetBoundaryToForm );

    /// creating from supplied vector of faces
    bool CreateFrom( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                     const typename std::vector<Face<dim>*>::const_iterator facesEnd );

    /// create faces from lower dimensional region
    bool CreateFrom( MeshManager<dim>& meshManager,
                     const Region<dim>& region,
                     const csmp::Index& mtrl_key,
                     BOX_BOUNDARY boxBoundary );
    
    /// creates surface / perimeter line of Faces around the region ( only for volume regions in 3D and surface regions in 2D )
    bool CreateAround( MeshManager<dim>&,
                       const Region<dim>& region,
                       BOX_BOUNDARY boxBoundary = IRREGULAR );

    /// creates surface / perimeter line of Faces between regions (the first is on the inside)
    bool CreateBetween( MeshManager<dim>&,
                        const Region<dim>&,
                        const Region<dim>& );


    // ----------------------------------------
    //  geometric properties and integrals
    // ----------------------------------------

    /// surface area of the boundary
    double64  Area() const;

    /// length of the boundary perimeter in 3D, else NaN
    double64  Perimeter() const;

    /// surface integral over the variable of interest
    double64  SurfaceIntegral( const PropertyDatabase<dim>&, const char* property ) const;
    
    /// reports whether the boundary lies on the outside of the whole model
    bool IsExternal() const;

    /// application of single-value Box boundary flagging
    void AtBoundary( BOX_BOUNDARY boxBoundary );

    /// reports the flagging of the boundary with respect to csmp::Box boundary convention
    BOX_BOUNDARY AtBoundary() const { return boundaryFlag_; }

    /// returns 1) elements of how many different spatial dimensions are contained, and 2) the highest element spatial dimension in boundary
    std::pair<int32, int32>  FaceSpatialDimensions() const;

    void Out() const;

  protected:

    /// establishes connectivity between Faces if not already there, assigns boundary flag and finds perimeter
    void Initialize( BOX_BOUNDARY boxBoundary );

    /// returns local variables stored at face integration points
    IntegrationPointVariables FaceIntegrationPointVariables() const;

    /// returns local variables stored at on the boundary faces
    LocalVariables FaceVariables() const;

  protected:
    // binary IO
    template<class Var>
    bool Out( std::fstream& fp, PLACEMENT place, VARIABLE_TYPE vtype ) const;
    template<class Var>
    bool In( std::fstream& fp, PLACEMENT place, VARIABLE_TYPE vtype );

    BOX_BOUNDARY boundaryFlag_;
};



} // end csmp

#endif
