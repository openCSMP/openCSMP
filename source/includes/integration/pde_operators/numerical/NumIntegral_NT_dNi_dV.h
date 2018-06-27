#ifndef NUM_INTEGRAL_NT_DNI_DV_H
#define NUM_INTEGRAL_NT_DNI_DV_H

#include "MathOperatorLHS.h"
#include "CSMP_definitions.h"

namespace csmp {

enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

/**
    Accumulation of scalar coordinate axis aligned test-function * test function derivative products.
    @attention no operand values are required.

    @author Yan Zaretskiy & Sebastian Geiger (2009), coefficient free version by SKM 29/2/16
 
    @remarks refactored by SKM 19/1/2015
*/
template<size_t dim, class CELL=Element<dim> >
class NumIntegral_NT_dNi_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_NT_dNi_dV( const PropertyDatabase<dim>&,
                              /* no operand in this pde operator */
                              const char* basic,
                              const char* test );
    
    virtual ~NumIntegral_NT_dNi_dV();
  
    /// does nothing since there are no operands to read
    virtual void GetOperands( CELL& ) {}
  
    /// calculates the required finite element integral
    virtual void ComputeContribution( CELL& );
  
    /// to chose the spatial derivate direction of interest; default is Y-axis
    void SpatialDerivative( SPATIAL_DERIVATIVE );
  
    /// to transpose the element matrix that will get accumulated; default is false
    void Transposed();
  
    virtual NumIntegral_NT_dNi_dV<dim,CELL>* clone() const { return new NumIntegral_NT_dNi_dV<dim,CELL> (*this); }
  
  private:
    NumIntegral_NT_dNi_dV();
    
    SPATIAL_DERIVATIVE   xyz_;     ///< direction of partial derivative that shall be taken
    bool                 transp_;  ///< transpose element matrix or not(default)
    DenseMatrix<DM_MIN>  TEMP;     ///< temporary matrices
};

} // csmp

#endif
