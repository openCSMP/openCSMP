#ifndef NUM_INTEGRAL_NT_DNI_DV_H
#define NUM_INTEGRAL_NT_DNI_DV_H

#include "MathOperatorLHS.h"
#include "CSMP_definitions.h"

namespace csmp {

enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

/**
    Accumulation of scalar-coordinate-axis-aligned test-function * test function derivative products.
    @attention no operand values are required.

    @author Yan Zaretskiy & Sebastian Geiger (2009), coefficient free version by SKM 29/2/16
 
    @remarks refactored by SKM 19/1/2015
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_NT_dNi_dV : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_NT_dNi_dV( const PropertyDatabase<dim>&,
                           /* no operand 'oper' in this pde operator */
                           const char* basic,
                           const char* test );
    
    virtual ~NumIntegral_NT_dNi_dV();
  
    /// does nothing since there are no operands to read
    virtual void GetOperands( const CELL<dim>& ) { /* no operand to read */ }
  
    /// Computes N_transposed * DN_i product that is weighted by the determinant of Jacobian matrix.
    virtual void ComputeContribution( const CELL<dim>& );
  
    /// chose the desired spatial derivate; default is Y
    void SpatialDerivative( SPATIAL_DERIVATIVE );
  
    /// transposes LHS element matrix before it gets accumulated; default is false
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
