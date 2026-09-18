//
//  ModelSubDomain_Test.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//

#ifndef CSMP_MODEL_SUB_DOMAIN_TEST_H
#define CSMP_MODEL_SUB_DOMAIN_TEST_H

#include "Test.h"
#include "ModelSubDomain.h"
#include "CSMP_mathUtilities.h"
#include "Box.h"

namespace csmp {

/** Tests the functionality of the base class of Region, Boundary and SplitBoundary
    Only tests selected less frequently used, less obvious functionality that requires testing.
    Functions involving variables require testing in 2 and 3D.
    
        /// deletes nullptr cells, rebuilds node vector, sorts everything and re-establishes the perimeter face vectors after modifications of cells
@todo    void RebuildSubDomainAfterChangeOfCellVector();
    
    /// rebuilds subdomain on the basis of the cells that will be selected according to the supplied property constraints
    void UpdateCellMembershipApplyingConstraints( typename std::vector<CELL<dim>*>::const_iterator master_domain_start,
                                                  typename std::vector<CELL<dim>*>::const_iterator master_domain_end,
                                                  const PropertyConstraints& );

    /// removes any cells or node pointers that were set to zero elsewhere; returns number of cells removed
    size_t RemoveNullPointerCells();
    
    /// reference counting-based unique domain identifier  (0..n-1)
    size_t DomainIndex() const;

    bool IsContiguous() const;

    std::pair<CELL_SHAPE,bool>  SingleCellShapeDomain() const;

    /// output coordinates to user-defined vector variable [x1,y1,z1,...xn,yn,zn]
    void AssignNodeCoordinatesTo( const char* vector_prop );
  
    /// assigment of coordinate component 'x','y','z' to scalar variable of choice
    void AssignNodeCoordinatesTo( const char* scalar_prop, char coord );

    /// characteristics like 'length', 'area', 'volume' , 'aspect ratio', 'inner radius' are assigned to user-defined variable
    void AssignCellCharacteristicsTo( const char* characteristic, const char* var );
    
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
                               uint32_t component,
                               VARIABLE_FLAG new_status,
                               SUBDOMAIN_PART=COMPLETE );

    /// changes the flag of a particular variable component to new value if the scalar variable is withing the specified range
    void ChangePropertyStatusWhere( const char* property,
                                    uint32_t component,
                                    VARIABLE_FLAG new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    /// by default (i=0) returns status of scalar variable or first component of a vector or tensor variable; if i>0 flag of corresponding component is returned
    VARIABLE_FLAG  PropertyStatus( const char* variable, SUBDOMAIN_PART flag=COMPLETE , uint32_t i=0 ) const;

    void InterpolateNodeToCellProperty( const char* nprop, const char* eprop );
    void InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* eprop );
    void InterpolateIntegrationPointToCellProperty( const char* cprop, const char* eprop );
    void ExtrapolateCellToIntegrationPointProperty( const char* eprop, const char* cprop );
    void ExtrapolateCellToFacetIntegrationPointProperty( const char* eprop, const char* fipprop );
    void ExtrapolateCellToNodeProperty( const char* eprop, const char* nprop, bool by_distance=true );
    void ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop );

    /// arithmetic (number as opposed to volume weighted) average
    double Average( const char* property ) const;
    bool   CopyGradientOfProperty_A_To_B( const char* node_prop, const char* cell_prop );
    void   CopyReplace( const char* from, const char* to );

*/
class ModelSubDomain_Test : public Test {
  public:
    void run() override final;

    /// compares node locations and connectivity (is used also in other tests)
    template<uint32_t dim,template<uint32_t> class simplicial_complex>
    bool CompareModelSubdomains( const ModelSubDomain<dim,simplicial_complex>&,
                                 const ModelSubDomain<dim,simplicial_complex>&,
                                 bool verbose );
  private:

    // LEGACY TESTS (P.L. 2011)
    /*
        - includes IdentifyPerimeter()
    */
    /// tests method that creates neighbor connectivity inside of CSMP
    bool Test_EstablishNeighborConnectivity();
    
    void Test_RebuildCellAndNodeVectors();



    // TESTING SPECIFIC MODELSUBDOMAIN FUNCTIONALITY

    /**
    std::string Name() const noexcept;
    void Name( const std::string& ) noexcept;
    
    void ScheduleForRebuild() noexcept;
    bool NeedsRebuild() const noexcept;

    size_t RebuildNodeVector();
    void CreateNodePointerVector();
    void CreateNodePointerVector( std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& contacting_cells );
    void SortVectors( size_t interior_cells, size_t interior_nodes );

    int32_t DomainIndex() const noexcept;
    virtual size_t  RenumberNodes() const noexcept;
    size_t  RenumberCells() const noexcept;
    void    UpdateMemberIndexes() const noexcept;
    std::vector<size_t>  MemberCellIndexes() const noexcept;
    typename std::vector<CELL<dim>*>::const_iterator  PerimeterCellsBegin() const noexcept;
    typename std::vector<CELL<dim>*>::const_iterator  CellsEnd() const noexcept;

    bool IsUnique()
    void IsUnique( bool unique_domain )
    bool IsContiguous() const;
    std::pair<CELL_SHAPE,bool>  SingleCellShapeDomain() const;
    */
    void Test_Basics();
  
  
    /**
    void UpdateTopoTypeNodeFlags();
    void BuildPerimeterFaceVector( int64_t interior_cells );
    void RebuildSubDomainAfterChangeOfCellVector();
    void RebuildSubDomainAfterChangeOfNodeVector();

    void UpdateCellMembershipApplyingConstraints( typename std::vector<CELL<dim>*>::const_iterator master_domain_start,
                                                  typename std::vector<CELL<dim>*>::const_iterator master_domain_end,
                                                  const PropertyConstraints& );
    size_t RebuildCellAndPerimeterFaceVector();
     */
    void Test_SubDomainConstructionMethods();


    /**
    bool              Empty() const noexcept;
    size_t            Nodes() const noexcept;
    size_t            InteriorNodes() const noexcept;
    size_t            PerimeterNodes() const noexcept;
    size_t            IntegrationPoints() const noexcept;
    size_t            SectorIntegrationPoints() const noexcept;
    size_t            FacetIntegrationPoints() const noexcept;
    size_t            Cells() const noexcept;
    size_t            InteriorCells() const noexcept;
    size_t            PerimeterCells() const noexcept;
    bool              Contains( const CELL<dim>* const ) const noexcept;
    bool              Contains( const Node<dim>* const ) const noexcept;
    bool              IsPerimeterNode( const csmp::Node<dim>* const ) const noexcept;
    bool              IsPerimeterCell( const CELL<dim>* const ) const noexcept;
    uint32_t          PerimeterFaces( size_t cell_idx ) const;
    uint32_t          PerimeterFace( size_t cell_idx, uint32_t face ) const;
    bool              IsPerimeterNode( const size_t nidx ) const;
    bool              IsPerimeterCell( const size_t eidx ) const;
    size_t            SharedPerimeterNodes( start, end )

    csmp::Node<dim>*  N( size_t n ) const noexcept;
    CELL<dim>*        E( size_t n ) const noexcept;

    */
    void Test_SubDomainDiagnostics();


    /**
    std::pair<int32_t,int32_t>  SpatialDimensions() const;
    Point<dim> Centroid() const;
    Point<dim> CenterOfGravity( const csmp::Index& rho_key ) const;
    void  MinMaxCoordinates( Point<dim>& xyz_min, Point<dim>& xyz_max ) const;
    void AssignNodeCoordinatesTo( const char* vector_prop );
    void AssignNodeCoordinatesTo( const char* scalar_prop, char coord );
    void AssignCellCharacteristicsTo( const char* characteristic, const char* var );
    void NodeAttributesToCSV();
    */
    void Test_GeometricOperations();

    /**
    void InputPropertyValue( const char* input_prop, const Var& new_value, SUBDOMAIN_PART sd=COMPLETE );
    void InputPropertyValue( const char* input_prop, const Var& new_value, VARIABLE_FLAG do_not_overwrite, SUBDOMAIN_PART sd=COMPLETE );
    void ChangePropertyStatus( const char* property, VARIABLE_FLAG new_status_of_scalar, SUBDOMAIN_PART=COMPLETE );
    void ChangePropertyStatusWhere( const char* property,
                                    VARIABLE_FLAG new_status_of_scalar,
                                    double min_value_to_change,
                                    double max_value_to_change );

    void ChangePropertyStatus( const char* property,
                               const std::vector<VARIABLE_FLAG>& new_status,
                               SUBDOMAIN_PART=COMPLETE );

     void ChangePropertyStatusWhere( const char* property,
                                    const std::vector<VARIABLE_FLAG>& new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    void ChangePropertyStatus( const char* property,
                               uint32_t component,
                               VARIABLE_FLAG new_status,
                               SUBDOMAIN_PART=COMPLETE );

    void ChangePropertyStatusWhere( const char* property,
                                    uint32_t component,
                                    VARIABLE_FLAG new_status,
                                    double min_value_to_change,
                                    double max_value_to_change );

    VARIABLE_FLAG  PropertyStatus( const char* variable, SUBDOMAIN_PART flag=COMPLETE , uint32_t i=0 ) const;
    void MinMaxOf( const char* property,   double& gmin, double& gmax ) const;
    void MinMaxOf( const csmp::Index&,     double& gmin, double& gmax ) const;
    void   CopyReplace( const char* from, const char* to );
    */
    void Test_PropertyManipulations();

    /**
    double Average( const char* property ) const;
    Point<dim> AverageUnitNormal() const;
    void InterpolateNodeToCellProperty( const char* nprop, const char* eprop );
    void InterpolateNodeToIntegrationPointProperty( const char* nprop, const char* eprop );
    void InterpolateIntegrationPointToCellProperty( const char* cprop, const char* eprop );
    void ExtrapolateCellToIntegrationPointProperty( const char* eprop, const char* cprop );
    void ExtrapolateCellToFacetIntegrationPointProperty( const char* eprop, const char* fipprop );
    void ExtrapolateCellToNodeProperty( const char* eprop, const char* nprop, bool by_distance=true );
    void ExtrapolateIntegrationPointToNodeProperty( const char* cprop, const char* nprop );
    bool   CopyGradientOfProperty_A_To_B( const char* node_prop, const char* cell_prop );
    */
    void Test_PropertyTransfer();
    
    
    /**
      PLACEMENT modelSubdomainType( const std::string& subdomain_name ) noexcept;
      size_t  sharedNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );
      size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>& );
      size_t  sharedPerimeterNodes( const ModelSubDomain<dim,CELL>&, const ModelSubDomain<dim,CELL>&, std::vector<Node<dim>*>& );
      size_t  sharedPerimeterCells( const ModelSubDomain<dim,CELL>& subdomain1, const ModelSubDomain<dim,CELL>& subdomain2,
                                    std::vector<std::pair<std::pair<CELL<dim>*,uint32_t>,std::pair<CELL<dim>*,uint32_t> > >& matching_cells );
                                    
      std::set<TOPOTYPE> nodeTopologyFlags( typename std::vector<Node<dim>*>::const_iterator first,
                                            typename std::vector<Node<dim>*>::const_iterator last );
                                    
      std::set<TOPOTYPE> nodeTopologyFlags( const SplitBoundary<dim>& );

      void readDomainIndexesFromBinaryFile( std::fstream&, SubDomainInfo& ); // tested elsewhere
    */
    void Test_NonMemberFunctions();
    
    /// helper that prints TOPOTYPE values to VTU
    template<uint32_t dim>
    void BoundaryAndTopoTypeFlagsToVTU( Model<dim>& model );

  private:
    const bool verbose_ = true;
};

} // end csmp

#endif /* CSMP_MODEL_SUB_DOMAIN_TEST_H */
