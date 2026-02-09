#ifndef CSMP_BOUNDARY_H
#define CSMP_BOUNDARY_H

#include "ModelSubDomain.h"
#include "Face.h"

namespace csmp {

class FiniteElementManager;
template<uint32_t> class FaceConstructionData;
template<uint32_t> class PropertyDatabase;
template<uint32_t> class MeshManager;
template<typename> class FEM_Data;
template<uint32_t> class Point;
template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Region;

/**

@brief A Boundary is lower-dimensional ModelSubDomain, consisting of Face objects.
Consisting of Faces endows the boundary with special functionality such
as access to higer-dimensional elements on either side and a distinct
treatment in the accumulation of PDE operators such as surface integrals
for a domain.

Like for the Model or Region, properties can be placed directly on them
saving storage and supporting subdomain specific treatment of processes.

Boundaries have BOX_BOUNDARY flags, telling whether they are internal or
external (box flags including IRREGULAR) boundaries.
Their nodes are flagged correspondingly, which is important because
it helps in the progressive construction of Boundaries to avoid accidential
duplication.

The nodes of a boundary have the same flag as the box boundary.
Nodes within a model are thus flagged INTERNAL unless they already have
a !NOT box-boundary flag.

Boundaries creation and removal are managed via the BoundaryInterface of the Model class.

@note Since Face is a separate variable placement variables placed on the
element are accessible via Read, Store or Status operations.
In practice, however, one always has access to the adjacent higher-dimensional
Element from which element properties can be obtained, for instance, to compute
2D solutions of the governing equations on the Boundary.

@section motivation Motivation

To apply areal boundary conditions, carry out computations or monitor processes
at material interfaces or model side boundaries.

To distinguish surface elements by turning them into Face objects 
for specific processing, and treating
boundary-specific calculations separately.

@author S.K. Matthai
@author P. Lang
@date 2010
@date 2022 major revision

*/
template<uint32_t dim>
class Boundary : public ModelSubDomain<dim, Face>,
                 public LocalVariableStorage<dim, Boundary>
 {
  public:
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
    ~Boundary() override final = default;

    Boundary<dim>&  operator=( const Boundary& );

    /// methods required for the LocalVariableStorage
    PLACEMENT Placement() const noexcept override final { return BOUNDARY; }
    bool      ValidVariable( const char* variableName ) const noexcept override final;

    /// RECONSTRUCTOR of boundary from index data stored in SubDomainInfo (call only prior to deleting anythin from colonies)
    Boundary( const PropertyDatabase<dim>&,
              MeshManager<dim>&,
              const SubDomainInfo&,    ///< contains correctly partitioned vectors and boundary faces
              BOX_BOUNDARY = IRREGULAR );

    /// Visitors
    virtual void Accept( Visitor<dim>& ) override;


    // -----------------------------------------------
    //  input/output to binary file
    // -----------------------------------------------

    /// for the assignment of properties that are unique to the instance of this subclass
    template<typename Var> requires CsmpVariable<dim, Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );

    /// as InputPropertyValue, but with overwrite protection for variable components that have the flag 'do_not_overwrite'
    template<typename Var> requires CsmpVariable<dim, Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );

    /// assigning the variable values from the FEM_DATA container to the corresponding property of the Boundary
    template<class Var> requires CsmpVariable<dim, Var>
    void InputVariableFrom( const char* property, const FEM_Data<Var>& property_values );

    template<class Var> requires CsmpVariable<dim, Var>
    void OutputVariableTo( const char* property, FEM_Data<Var>& property_values ) const;


    // ----------------------------------------
    //  boundary initialisation & modification
    // ----------------------------------------
    
    /// creates pointers to existing Face objects assuming that these have consecutive ID numbers starting with that of the last Element object
    size_t AccumulateByNumber( MeshManager<dim>&, std::vector<size_t>& cell_ids );

    /// creating from supplied vector of faces
    bool CreateFrom( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                     const typename std::vector<Face<dim>*>::const_iterator facesEnd );

    /// create faces from lower dimensional region, deleting these and assigning nullptr values
    bool CreateFrom( Region<dim>& lower_dimensional_region,
                     MeshManager<dim>& meshManager,
                     BOX_BOUNDARY boxBoundary );


    // ----------------------------------------
    //  geometric properties and integrals
    // ----------------------------------------

    /// in case the boundary consists of line segments
    double  Length() const;

    /// surface area of the boundary
    double  Area() const;

    /// length of the boundary perimeter in 3D, else NaN
    double  Perimeter() const;

    /// surface integral over the variable of interest
    double  SurfaceIntegral( const char* property ) const;
    
    /// reports whether the boundary lies on the outside of the whole model
    bool IsExternal() const;

    /// application of single-value Box boundary flagging
    void AtBoundary( BOX_BOUNDARY boxBoundary );

    /// reports the flagging of the boundary with respect to csmp::Box boundary convention
    BOX_BOUNDARY AtBoundary() const { return boundaryFlag_; }

    /// returns 1) elements of how many different spatial dimensions are contained, and 2) the highest element spatial dimension in boundary
    std::pair<int32_t, int32_t>  FaceSpatialDimensions() const;

    void Out() const;

  protected:

    /// assigns box boundary flags to nodes unless they do not have a NOT value
    void InitializeBoundaryFlags( BOX_BOUNDARY boxBoundary );

    /// returns local variables stored at face integration points
    IntegrationPointVariables FaceIntegrationPointVariables() const;

    /// returns local variables stored at on the boundary faces
    LocalVariables FaceVariables() const;

  protected:
    // binary IO
    template<class Var> requires CsmpVariable<dim, Var>
    bool Out( std::fstream& fp, PLACEMENT place, VARIABLE_TYPE vtype ) const;
    template<class Var> requires CsmpVariable<dim, Var>
    bool In( std::fstream& fp, PLACEMENT place, VARIABLE_TYPE vtype );

    BOX_BOUNDARY boundaryFlag_;
};



} // end csmp

#endif
