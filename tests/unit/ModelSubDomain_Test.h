//
//  ModelSubDomain_Test.hpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 25/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#ifndef ModelSubDomain_Test_hpp
#define ModelSubDomain_Test_hpp

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
    virtual void run();

    /// tests method that creates neighbor connectivity inside of CSMP
    bool Test_EstablishNeighborConnectivity();
    
    void Test_RebuildCellAndNodeVectors();

    /// compares node locations and connectivity
    template<uint32_t dim,template<uint32_t> class simplicial_complex>
    bool CompareModelSubdomains( const ModelSubDomain<dim,simplicial_complex>&,
                                 const ModelSubDomain<dim,simplicial_complex>&,
                                 bool verbose );
  private:
    const bool verbose_ = false;
};

} // end csmp

#endif /* ModelSubDomain_Test_hpp */
