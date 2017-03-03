#ifndef CGNS_ELEMENT_SPECIFICATIONS_H
#define CGNS_ELEMENT_SPECIFICATIONS_H

#include "CSMP_definitions.h"
#include "FiniteElement.h"
#include "Box.h"

namespace csmp {

class CGNS_ElementSpecifications {
  public:
    CGNS_ElementSpecifications();
    ~CGNS_ElementSpecifications();

    int  CGNS_Type( const std::string& FEtype ) const;
    std::string CGNS_TypeName( int ) const;

    /// from csmp to cgns
    int  CGNS_TypeFrom_CSMP_Type( int CSMP_finite_element_type ) const;
    int  CGNS_TypeFrom_CSMP_TypeName( const std::string& CSMP_finite_element_type ) const;
    std::string  CGNS_TypeNameFrom_CSMP_Type( int CSMP_finite_element_type ) const;
    std::string  CGNS_TypeNameFrom_CSMP_TypeName( const std::string& CSMP_finite_element_type ) const;

    /// from cgns to csmp
    csmp::CSMP_FEM_TYPE  CSMP_TypeFrom_CGNS_Type( int CGNS_finite_element_type, bool isoparametric, size_t dim ) const;
    csmp::CSMP_FEM_TYPE  CSMP_TypeFrom_CGNS_TypeName( const std::string& CGNS_finite_element_type, bool isoparametric, size_t dim ) const;
    std::string    CSMP_TypeNameFrom_CGNS_Type( int CGNS_finite_element_type, bool isoparametric, size_t dim ) const;
    std::string    CSMP_TypeNameFrom_CGNS_TypeName( const std::string& CGNS_finite_element_type, bool isoparametric, size_t dim ) const;

};


/**
 
@class CGNS_ElementSpecifications  CGNS_ElementSpecifications "interfaces\CGNS_ElementSpecifications.h"
@author R.Manasipov
@date 2015
*/

} // csmp


#endif
