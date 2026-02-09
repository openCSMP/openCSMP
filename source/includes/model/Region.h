#ifndef CSMP_REGION_H
#define CSMP_REGION_H

#include "Element.h"
#include "ModelSubDomain.h"

namespace csmp {

class PropertyData;
class PropertyConstraints;
struct LocalVariables;
template<uint32_t>   class PropertyDatabase;
template<uint32_t>   class Node;
template<uint32_t>   class VSet;
template<uint32_t>   class Region;
template<uint32_t>   class Boundary;
template<uint32_t>   class MeshManager;
template<typename> class FEM_Data;

/// returns number of nodes that are shared by the two subdomains; defined in ModelSubDomain.cpp
template<uint32_t dim, template<uint32_t> class CELL>
size_t  sharedNodes( const ModelSubDomain<dim, CELL>&, const ModelSubDomain<dim, CELL>& );

/// returns number of nodes on the subdomain perimeters that are shared by the two subdomains; defined in ModelSubDomain.cpp
template<uint32_t dim, template<uint32_t> class CELL>
size_t  sharedPerimeterNodes( const ModelSubDomain<dim, CELL>&, const ModelSubDomain<dim, CELL>& );

// boolean operations (see also region related functions below class declaration)
// -----------------------------------------------------------------------------
template<uint32_t dim>
size_t  groupUnion( const Region<dim>&, const Region<dim>&, Region<dim>& combined_regions );

template<uint32_t dim>
size_t  intersection( const Region<dim>&, const Region<dim>&, Region<dim>& intersection_region );

template<uint32_t dim>
size_t  difference( const Region<dim>&, const Region<dim>&, Region<dim>& difference_region );

template<uint32_t dim>
size_t  symmetricDifference( const Region<dim>&, const Region<dim>&, Region<dim>& sym_difference_region );

template<uint32_t dim>
size_t  sharedElements( const Region<dim>&, const Region<dim>& );

template<uint32_t dim>
bool  hasLowerDimensionalRepresentation( const Region<dim>& );

/// if all BOX_BOUDARY flags != NOT region is on an external boundary; method highlights potential inconsistencies with topology
template<uint32_t dim>
bool  formsPartOfExternalBoundary( const Region<dim>&, bool check_topology_as_well=true );

template<uint32_t dim>
bool  containsVolumeElements( const Region<dim>& );

template<uint32_t dim>
bool  containsSurfaceElements( const Region<dim>& );

template<uint32_t dim>
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
template<uint32_t dim>
class Region : public ModelSubDomain<dim, Element>,
               public LocalVariableStorage<dim, Region> {
  public:

    /// construction of named empty region with appropriately sized property storage
    Region( std::string regionname, const PropertyDatabase<dim>&, bool is_unique );

    Region( const Region& );
    Region( Region&& );
    Region& operator=( const Region& );

    virtual ~Region() = default;

    // ---------------------------------------------
    // reconstruction of regions that existed before
    // ---------------------------------------------

    /// RECONSTRUCTOR for regions via the MeshManager (call only prior to deleting anythin from colonies)
    Region( const PropertyDatabase<dim>&,
            MeshManager<dim>&,
            const SubDomainInfo&, bool is_unique );  ///< contains correctly partitioned vectors and boundary faces

    /// RECONSTRUCTOR for regions via the nodes and elements which are explored by the MeshManager
    Region( const PropertyDatabase<dim>&,
            const std::deque<Node<dim>*>&,
            const std::deque<Element<dim>*>&,
            const SubDomainInfo&, bool is_unique );  ///< contains correctly partitioned vectors and boundary faces


   // --------------------------------------------
   // Property input/output
   // --------------------------------------------

    /// Sets the Region ID to ModelSubDomain::domain_idx_ for all elements of unique region; override as is needed
    void SetRegion_ID( int32_t region_id_for_non_unique_regions = UNSPECIFIED );

    /// for the assignment of properties that are unique to the instance of this subclass
    template<typename Var>  requires CsmpVariable<dim, Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );
      
    /// as InputPropertyValue, but with overwrite protection for variable components that have the flag 'do_not_overwrite'
    template<typename Var>  requires CsmpVariable<dim, Var>
    void InputPropertyValue( const char* input_prop, const Var& new_value,
                            VARIABLE_FLAG flag_of_values_to_be_preserved,
                            SUBDOMAIN_PART = COMPLETE );

    /// Local variable storage interface
    virtual PLACEMENT Placement() const noexcept override final { return REGION; }

    /// checks whether variable placement is compatible with placement on a region
    virtual bool  ValidVariable( const char* variableName ) const noexcept override final;

    /// outputs region into VSet polygonal data container; all properties may be output as well
    void OutputTo( VSet<dim>& vset, bool with_properties = true ) const;

    /// outputs region into VSet polygonal data container, including a selected list of properties
    void OutputTo( VSet<dim>& vset, const std::map<std::string,Index>& properties_to_be_output_for_region ) const;

    /// output all distributed properties into the supplied VSet
    void OutputDataTo( VSet<dim>& vset ) const;

    /// outputs finite-volume sector and facet variables to VSet
    void OutputFvDataTo( VSet<dim>& vset ) const;

    /// outputs distributed variable to variable storage container of type PropertyData
    PropertyData  OutputVariableTo( const char* region_property ) const;

    /// outputs distributed variable to FEM_Data container
    template<class Var> requires CsmpVariable<dim, Var>
    void OutputVariableTo( const char* property, FEM_Data<Var>& data_local ) const;

    /// inputs distributed variable from file into region; @note indices must match
    template<class Var> requires CsmpVariable<dim, Var>
    void InputVariableFrom( const char* property, const FEM_Data<Var>& vdata );


    // --------------------------------------------
    // region building & modification
    // --------------------------------------------

    /// creates pointers to non const elements in the mesh directly from MeshManager; returns the number of elements of the new region
    size_t AccumulateAll( MeshManager<dim>& );

    /// creates element pointers assuming that the order in which the elements are stored in the MeshManager matches that in the element_ids vector; no 'idx' searching
    size_t AccumulateByNumber( MeshManager<dim>&, std::vector<size_t>& element_ids );

    // Accumulate based on pointers to elements

    /// accumulate those elements into a region whose id matches one of the numbers contained in vector 'element_ids' using binary_search; @attention use after mesh modification
    size_t AccumulateByNumber( typename std::vector<csmp::Element<dim>*>::const_iterator start,
                               typename std::vector<csmp::Element<dim>*>::const_iterator end,
                               std::vector<size_t>& element_ids );

    /// accumulates range of elements into a region identified by constant pointers created by AccumulateAll; returns # of accumulated elements
    size_t Accumulate( typename std::vector<csmp::Element<dim>*>::const_iterator start,
                       typename std::vector<csmp::Element<dim>*>::const_iterator end );

    /// accumulate a range of elements into a region identified by constant pointers created by AccumulateAll; returns # of accumulated elements
    size_t Accumulate( typename std::set<csmp::Element<dim>*>::const_iterator start,
                       typename std::set<csmp::Element<dim>*>::const_iterator end );

    // Accumulate based on property values

    /// accumulates region whose elements have properties in the ranges defined inside of the PropertyConstraints object; returns number of elements found
    size_t AccumulateWithinRange( MeshManager<dim>&, const PropertyConstraints& );

    /// accumulates elements where (at least one node) has 'property' values in the range defined by 'min/max', returns number of elements found
    size_t AccumulateWithinRange( MeshManager<dim>&, const char* property, double min, double max );

    /// accumulates region of elements whose barycenter lies within the defined bounding box
    size_t AccumulateRectangularRegion( MeshManager<dim>&, const Point<dim>& xyz_min, const Point<dim>& xyz_max );

   /// merges supplied region with the current one
    void  Add( const Region& );

    /// removes those elements in the region whose id matches one of the numbers contained in vector 'element_ids'
    size_t RemoveByNumber( std::vector<size_t>& element_ids, std::vector<csmp::Element<dim>*>& ptrs_to_removed_elements );

    /// removes those elements from the target region whose pointers matches the ones in the range supplied and subsequently rebuilds the region
    size_t RemoveRange( typename std::vector<csmp::Element<dim>*>::iterator begin,
                        typename std::vector<csmp::Element<dim>*>::iterator end );

    // ----------------------------------------
    // geometry manipulations
    // ----------------------------------------

    /// movement of the nodes in the region
    void MoveNodeCoordinatesBy( const char* vector_variable );

    // ----------------------------------------
    // diagnostics and property calculations
    // ----------------------------------------

    /// Visitors
    virtual void Accept( Visitor<dim>& ) override;

    /// tests whether the current region includes the specified one
    bool   Includes( const Region& ) const;

    /// returns 1) elements of how many different spatial dimensions are contained, and 2) the highest element spatial dimension in subdomain
    std::pair<int32_t, int32_t>  ElementSpatialDimensions() const;

    /// returns either the lenght of the perimeter (2D) or the surface area of the region (3D)
    double  SurfaceArea() const;

    /// returns either the area of a surface region, the length of a line region or the region's volume
    double  Volume( bool multiply_with_porosity = false ) const;

    /// integration over volume elements
    double  VolumeIntegral( const char* property, bool multiply_with_porosity, bool verbose = true ) const;

    /// integration over volume and surface elements using thickness attribute (=1 for volume elements)
    double  VolumeIntegral_x_Thickness( const char* property, bool multiply_with_porosity = false ) const;

    // ----------------------------------------
    // screen output
    // ----------------------------------------

    /// no variables in addition to ModelSubDomain
    void Out() const { std::cout << "\nRegion:Out:\n"; ModelSubDomain<dim, Element>::Out(); }

  protected:

    /// returns local variables for integration points
    IntegrationPointVariables ElementIntegrationPointVariables() const;

    /// returns local variables for elements
    LocalVariables ElementVariables() const;
};


/// supporting SplitBoundary creation, this methods finds outside elements that only touch the region parameter with nodes
template<uint32_t dim>
size_t outsideElementsWithNodesTouchingPerimeter( const Region<dim>&,
                                                  const std::vector<Node<dim>*>& perimeter_nodes,
                                                  const std::vector<std::pair<std::pair<Element<dim>*,uint32_t>,
                                                               std::pair<Element<dim>*,uint32_t> > >& perimeter_cells,
                                                  std::map<Node<dim>*,std::map<Element<dim>*,uint32_t>>& touching_elmts );
} // csmp

#endif
