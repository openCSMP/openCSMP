#ifndef CSMP_REGION_H
#define CSMP_REGION_H

#include "Element.h"
#include "ModelSubDomain.h"

namespace csmp {

template<size_t> class PropertyDatabase;
class PropertyData;
class PropertyConstraints;
struct LocalVariables;
template<typename> class FEM_Data;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class VSet;
template<size_t> class Region;
template<size_t> class Boundary;

/// returns number of nodes that are shared by the two subdomains
template<size_t dim,template<size_t> class SIMPLEX>
size_t  sharedNodes( const ModelSubDomain<dim,SIMPLEX>&, const ModelSubDomain<dim,SIMPLEX>& );

/// returns number of nodes on the subdomain perimeters that are shared by the two subdomains
template<size_t dim,template<size_t> class SIMPLEX>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim,SIMPLEX>&, const ModelSubDomain<dim,SIMPLEX>& );

// boolean operations (see also region related functions below class declaration)
// -----------------------------------------------------------------------------
template<size_t dim>
size_t  groupUnion( const Region<dim>&, const Region<dim>&, Region<dim>& combined_regions );

template<size_t dim>
size_t  intersection( const Region<dim>&, const Region<dim>&, Region<dim>& intersection_region );

template<size_t dim>
size_t  difference( const Region<dim>&, const Region<dim>&, Region<dim>& difference_region );

template<size_t dim>
size_t  symmetricDifference( const Region<dim>&, const Region<dim>&, Region<dim>& sym_difference_region );

template<size_t dim>
size_t  sharedElements( const Region<dim>&, const Region<dim>& );

template<size_t dim>
bool  isOfLowerDimensionalRepresentation( const Region<dim>& );

template<size_t dim>
bool  containsVolumeElements( const Region<dim>& );

template<size_t dim>
bool  containsSurfaceElements( const Region<dim>& );

template<size_t dim>
bool  containsLineElements( const Region<dim>& );

/**

@brief To associate and access elements which share certain properties.
Once such a distinction has been made these can be drawn
or output to specific visualization software.
Examples: shale horizons, faults etc.
 
@author S.K. Matthai
@author Stephen G. Roberts
@date 1999
 
@section consequences Consequences

When a new group is formed a group internal flag will be assigned to
each node, constraint point and element. 
There are two group-internal object flags, PLAIN and BOUNDARY. When a
new group is formed, the Region method IdentifyBoundaryAs() assigns the 
group-internal object flags, depending on whether the nodes, constraint
points or elements in the group lie at the group boundary or inside of 
the group. Elements are assigned a boundary flag if at least one of 
their faces coincides with the region boundary.

*/
template<size_t dim>
class Region : public ModelSubDomain<dim,Element> {

  public:

    // ---------------------------------------------
    // construction of regions in a new model
    // ---------------------------------------------

    /// construction of named empty region with appropriately sized property storage
    Region( std::string regionname, const PropertyDatabase<dim>& );
  
    Region( const Region& );
    Region( Region&& );
    Region& operator=( const Region& );

    virtual ~Region();

    // ---------------------------------------------
    // reconstruction of regions that existed before
    // ---------------------------------------------
  
    /// re-constructor for regions via the MeshManager
    Region( const PropertyDatabase<dim>&,
                  MeshManager<dim>&,
            const SubDomainInfo& );  ///< contains correctly partitioned vectors and boundary faces

    // --------------------------------------------
    // Property input/output
    // --------------------------------------------

    /// outputs region into VSet polygonal data container; all properties may be output as well
    void OutputTo( VSet<dim>& vset, bool with_properties=true ) const;

    /// outputs region into VSet polygonal data container, including a selected list of properties
    void OutputTo( VSet<dim>& vset, const std::map<std::string,Index>& properties ) const;
  
    /// output all distributed properties into the supplied VSet
    void OutputDataTo( VSet<dim>& vset ) const;
  
    /// outputs finite-volume sector and facet variables to VSet
    void OutputFvDataTo( VSet<dim>& vset ) const;

    /// outputs distributed variable to variable storage container of type PropertyData
    PropertyData  OutputVariableTo( const char* region_property ) const;

    /// outputs distributed variable to FEM_Data container
    template<class Var>
    void OutputVariableTo( const char* property, FEM_Data<Var>& data_local ) const;
  
    /// inputs distributed variable from file into region; @note indices must match
    template<class Var>
    void InputVariableFrom( const char* property, const FEM_Data<Var>& vdata );


    // --------------------------------------------
    // region building & modification
    // --------------------------------------------

    // Accumulate based on the root node of the mesh and its indirect connections to all elements
  
    /// breadth-first traversal of element tree in order to discover all elements in the model
    size_t AccumulateAll( const csmp::Node<dim>* root_node,
                          bool reestablishNeighborConnectivity = true );

    /// accumulate a range of elements into a region defined by iterators over an STL deque container
    void   Accumulate( typename std::deque<csmp::Element<dim> >::iterator start,
                       typename std::deque<csmp::Element<dim> >::iterator end );

    /// accumulate a range of elements into a region defined by constant iterators (accessors only) over an STL vector container
    void   Accumulate( typename std::vector<csmp::Element<dim>*>::const_iterator start,
                       typename std::vector<csmp::Element<dim>*>::const_iterator end );

    /// accumulate a range of elements into a region defined by iterators over the STL set (associative) container
    void   Accumulate( typename std::set<csmp::Element<dim>*>::const_iterator start,
                       typename std::set<csmp::Element<dim>*>::const_iterator end );

    // Accumulate based on pointers to elements

    /// accumulate those elements into a region whose id matches one of the numbers contained in vector 'element_ids'
    void   AccumulateByNumber( typename std::vector<csmp::Element<dim>*>::const_iterator start,
                               typename std::vector<csmp::Element<dim>*>::const_iterator end,
                               std::vector<size_t>& element_ids );

	void AccumulateByNumber( const csmp::MeshManager<dim>& mesh, std::vector<size_t>& element_ids );
    
    // Accumulate based on property values

    /// accumulates region whose elements have properties in the ranges defined inside of the PropertyConstraints object
    void   AccumulateWithinRange( typename std::vector<csmp::Element<dim>*>::const_iterator start, 
                                  typename std::vector<csmp::Element<dim>*>::const_iterator end, 
                                  const PropertyConstraints& );
  
    /// accumulates elements where (at least one node) has 'property' values in the range defined by 'min/max'
    void   AccumulateWithinRange( typename std::vector<csmp::Element<dim>*>::const_iterator start, 
                                  typename std::vector<csmp::Element<dim>*>::const_iterator end,
                                  const char* property,
                                  double64 min, double64 max );

	void AccumulateWithinRange( const csmp::MeshManager<dim>& mesh, const PropertyConstraints& );

	void AccumulateWithinRange( const csmp::MeshManager<dim>& mesh, const char* property, double64 min, double64 max );

    // Accumulate based on the location

    /// accumulates region of elements whose barycenter lies within the defined bounding box
    void   AccumulateRectangularRegion( typename std::vector<csmp::Element<dim>*>::const_iterator start,
                                        typename std::vector<csmp::Element<dim>*>::const_iterator end, 
                                        const Point<dim>& xyz_min, const Point<dim>& xyz_max );

	void AccumulateRectangularRegion( const csmp::MeshManager<dim>& mesh, const Point<dim>& xyz_min, const Point<dim>& xyz_max );
	
    /// Accumulate from the largest component in a mesh
    size_t FromLargestComponent( MeshManager<dim>& mesh, bool reestablishNeighborConnectivity );

    /// merges supplied region with the current one
    void   Add( const Region& );

    /// reestablish nodes based on element container
    void CreateNodePointerVector();

    /// excludes the intersection of elements of the 2 regions from the non-unique region
    bool RemoveRegionFromNonUniqueRegion( const char* non_unique_region, const char* region_to_substract );
    
    /// moves region to from the unique- to the non-unique regions map
    bool MoveRegionRegionToUniqueRegions( const char* unique_region );
  
    // TODO: deprecate since purpose duplicates that of boundaries
    /// creates surface / perimeter line of Elements between regions (the first is on the inside); TODO: @todo check whether this works
    bool CreateBetween( MeshManager<dim>&,
                        const FiniteElementManager&,
                        const Region<dim>&,
                        const Region<dim>& );

    // ----------------------------------------
    // geometry manipulations
    // ----------------------------------------

    /// movement of the nodes in the region
    void MoveNodeCoordinatesBy( const char* vector_variable );

    // ----------------------------------------
    // diagnostics and property calculations
    // ----------------------------------------

    /// Visitors
    virtual void Accept( Visitor<dim>& );

    /// tests whether the current region includes the specified one
    bool   Includes( const Region& ) const;

    /// returns 1) elements of how many different spatial dimensions are contained, and 2) the highest element spatial dimension in subdomain
    std::pair<int32,int32>  ElementSpatialDimensions() const;
  
    /// returns either the lenght of the perimeter (2D) or the surface area of the region (3D)
    double64  SurfaceArea() const;
  
    /// returns either the area of a surface region, the length of a line region or the region's volume
    double64  Volume( bool multiply_with_porosity=false ) const;
  
    /// integration over volume elements
    double64  VolumeIntegral( const char* property, bool multiply_with_porosity, bool verbose=true ) const;
  
    /// integration over volume and surface elements using thickness attribute (=1 for volume elements)
    double64  VolumeIntegral_x_Thickness( const char* property, bool multiply_with_porosity=false ) const;

    // ----------------------------------------
    // screen output
    // ----------------------------------------

    /// no variables in addition to ModelSubDomain
    void Out() const {  std::cout<<"\nRegion:Out:\n"; ModelSubDomain<dim,Element>::Out(); }

  public:
  
    /// Local variable storage interface
    virtual PLACEMENT Placement() const { return REGION; }
	    
    // TODO: trivial and unclear
    virtual bool  ValidVariable( const char* variableName ) const;

  protected:

    /// returns local variables for integration points
    IntegrationPointVariables ElementIntegrationPointVariables() const;
  
    /// returns local variables for elements
    LocalVariables ElementVariables() const;
};

} // csmp   
    
#endif




















