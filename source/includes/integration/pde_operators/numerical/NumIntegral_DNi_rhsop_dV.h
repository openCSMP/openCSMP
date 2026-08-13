#ifndef NUM_INTEGRAL_DNI_RHSOP_DV_H
#define NUM_INTEGRAL_DNI_RHSOP_DV_H

#include "MathOperatorRHS.h"
#include "CSMP_definitions.h"

namespace csmp {

enum SPATIAL_DERIVATIVE { X_DIRECTION=0, Y_DIRECTION=1, Z_DIRECTION=2 };

std::string parse( SPATIAL_DERIVATIVE );

template<uint32_t> class Element;

/**
    Volume integral over the gradient of the Operand in the direction i (i=x,y,z).
    Accumulation into the right-hand side of the linear algebraic system Ax=b
    
    @author SKM 6-12-2016
    
    @note created to compute a grad P right-handside for a 2-step Stokes
    lubrication solver.
*/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_DNi_rhsop_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_DNi_rhsop_dV( const PropertyDatabase<dim>&,
                              const char* oper,
                              const char* test );
    
    ~NumIntegral_DNi_rhsop_dV();
  
    /// reads the operand values from the nodes and stores them in a vector
    void GetOperands( const CELL<dim>& ) override final;
  
    void ComputeContribution( const CELL<dim>& ) override final;
  
    /// to chose the spatial derivate direction of interest; default is Y-axis
    void SpatialDerivative( SPATIAL_DERIVATIVE );
  
    NumIntegral_DNi_rhsop_dV<dim,CELL>* clone() const override final { return new NumIntegral_DNi_rhsop_dV<dim,CELL> (*this); }
  
  private:
    NumIntegral_DNi_rhsop_dV();
    
    SPATIAL_DERIVATIVE           xyz_;     ///< direction of partial derivative of interest
    std::vector<ScalarVariable>  op_vec_;  ///< nodal operand values
    DenseMatrix<DM_MIN>          DERIV_;   ///< shape function derivative matrix
};

} // csmp

#endif
