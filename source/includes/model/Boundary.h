#ifndef CSMP_BOUNDARY_H
#define CSMP_BOUNDARY_H

#include "ModelSubDomain.h"
#include "Face.h"

#include <cstdio>
#include <vector>
#include <utility>

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class MeshManager;
class FiniteElementManager;
template<typename> class FEM_Data;
template<size_t> class Point;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Region;


/** 

@brief Boundaries are lower dimensional regions=subdomains of the Model
which consist of Face objects

@author P. Lang
@author S.K. Matthai
@date 2010

Boundaries are collections of Faces (interior and exterior model boundaries)
made accessible for specific computations.
Like Region objects, they are contructed from the master 
region that contains all elements. 
It must be initialized before boundary operations can proceed.

@section motivation Motivation

To address and carry out computations at material interfaces or model
boundaries.

To avoid conditional processing of all elements of a model, by treating
boundary-specific calculations separately.

@section consequences Consequences

Through the use of boundaries essential conditions and coupling terms
arising at material interfaces can be computed with greater ease and
efficiency.


@todo (3) Shall we flag inner/outer parents (adjacent elements) as at boundary or depricate?
@todo (3) Port any remaining AtBoundary functionality to new framework (EDGE ?) (A)

*/
template<size_t dim>
class Boundary : public ModelSubDomain<dim,Face> {
  public:

    // ------------------------------------------------
    // construction of boundaries from scratch
    // ------------------------------------------------

    /// SKM new constructor: creating boundary from supplied vector of faces
    Boundary( const std::string& boundary_name,
              const PropertyDatabase<dim>&,
              typename std::vector<Face<dim>*>::iterator facesBegin,
              typename std::vector<Face<dim>*>::iterator facesEnd,
              BOX_BOUNDARY );

    /// create empty boundary with appropriately resized property storage
    Boundary( std::string boundaryname, const PropertyDatabase<dim>&, BOX_BOUNDARY flag );

    Boundary( const Boundary& );
    Boundary( Boundary&& );

    virtual ~Boundary();

    Boundary<dim>&  operator=( const Boundary& );

    // ------------------------------------------------
    // reconstruction of boundaries that existed before
    // ------------------------------------------------
  
    /// SKM new constructor: re-constructor of boundary from index data stored in SubDomainInfo
    Boundary( const PropertyDatabase<dim>&,
              MeshManager<dim>&,       ///< not constant since write access is granted to boundary
              const SubDomainInfo&,    ///< contains correctly partitioned vectors and boundary faces
              BOX_BOUNDARY=IRREGULAR );

    /// Local variable storage interface
    virtual PLACEMENT Placement() const { return BOUNDARY; }
    virtual bool      ValidVariable( const char* variableName ) const;

    /// Visitors
    virtual void Accept( Visitor<dim>& );


    // -----------------------------------------------
    // binary input/output
    // -----------------------------------------------

    /// assigning the variable values from the FEM_DATA container to the corresponding property of the Boundary
    template<class Var>
    void InputVariableFrom( const char* property, const FEM_Data<Var>& );

    template<class Var>
    void OutputVariableTo( const char* property, FEM_Data<Var>& ) const;

    // ----------------------------------------
    // building & modification
    // ----------------------------------------

// PHILIP's METHODS
    /// if parts of this boundary coincide with subsetRegion, this part will pasted to subsetBoundaryToForm
    bool Divide( const Region<dim>& subsetRegion,
                 Boundary<dim>& subsetBoundaryToForm );

    /// Removes faces from this boundary and assigns to subsetBoundaryToForm
    bool Divide( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                 const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                 Boundary<dim>& subsetBoundaryToForm );

    /// creating from supplied vector of faces
    bool CreateFrom( const typename std::vector<Face<dim>*>::const_iterator facesBegin,
                     const typename std::vector<Face<dim>*>::const_iterator facesEnd,
                     bool updateFaceConnectivity = true,
                     bool updateIndexes = true);

    /// create faces from lower dimensional region
    bool CreateFrom( MeshManager<dim>& meshManager,
                     const Region<dim>& region,
                     const csmp::Index& mtrl_key,
                     BOX_BOUNDARY boxBoundary );

    /// creates surface / perimeter line of Faces around the region ( only for volume regions in 3D and surface regions in 2D )
    bool CreateAround( MeshManager<dim>&,
                       const FiniteElementManager&,
                       const Region<dim>& region,
                       BOX_BOUNDARY boxBoundary = IRREGULAR );

    /// creates surface / perimeter line of Faces between regions (the first is on the inside)
    bool CreateBetween( MeshManager<dim>&,
                        const FiniteElementManager&,
                        const Region<dim>&,
                        const Region<dim>& );

    /// establishes connectivity, assigns boundary flags and initializes LVS
    void Initialize( BOX_BOUNDARY boxBoundary = IRREGULAR , bool updateNeighborConnectivity = true, bool updateIndexes = true );

    /// reestablish nodes based on element container
    void CreateNodePointerVector();

    // ----------------------------------------
    // various information output
    // ----------------------------------------

    /// surface area of the boundary
    double64  Area() const;
  
    /// length of the perimeter of the boundary
    double64  Perimeter() const;
  
    /// surface integral over the variable of interest
    double64  SurfaceIntegral( const PropertyDatabase<dim>&, const char* property ) const;

    /// application of single-value Box boundary flagging
    // TODO: should this also do the flagging of the edges etc.
    void AtBoundary( BOX_BOUNDARY boxBoundary );
  
    /// reports the flagging of the boundary with respect to csmp::Box boundary convention
    BOX_BOUNDARY AtBoundary() const { return boundaryFlag_; }

    /// returns 1) elements of how many different spatial dimensions are contained, and 2) the highest element spatial dimension in boundary
    std::pair<int32,int32>  FaceSpatialDimensions() const;

    // ----------------------------------------
    // screen output
    // ----------------------------------------
    void Out() const;

  protected:

    /// returns local variables stored at face integration points
    IntegrationPointVariables FaceIntegrationPointVariables() const;
  
    /// returns local variables stored at on the boundary faces
    LocalVariables FaceVariables() const;

  private:

    Boundary();

    // auxilliary binary IO
    template<class Var>
    bool Out( FILE* fp, PLACEMENT place, VARIABLE_TYPE vtype ) const;
    template<class Var>
    bool In( FILE* fp, PLACEMENT place, VARIABLE_TYPE vtype );

    BOX_BOUNDARY boundaryFlag_;
};



} // end csmp

#endif
