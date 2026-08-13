#ifndef INTEGRAL_NT_OP_DNI_DV_H
#define INTEGRAL_NT_OP_DNI_DV_H

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// buoyancy for instance
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_NT_op_dNi_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_NT_op_dNi_dV( const PropertyDatabase<dim>&, 
                           const char* oper, const char* mtrl, const char* test );
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    /// the space dimension of the partial derivative which shall be considered
    void SpatialDerivative( uint32_t xyz=2 );
    
    /// in stead of MultiplyWithTimeIncrement() since that would multiply whole contribution
    void MaterialPropertyTimeMultiplier( double time_increment );

    Integral_NT_op_dNi_dV<dim,CELL>* clone() const override final
      { return new Integral_NT_op_dNi_dV<dim,CELL> (*this); }
  
  private:
    std::vector<double>         IPOL;
    DenseMatrix<DM_MIN>         DN;
    csmp::Index                 prop_key;
    ScalarVariable              prop1, prop2;
    VectorVariable<dim>         bcenter;
    const double                gravity;   // acceleration of gravity
    uint32_t                    xyz; // 1=x, 2=y, 3=z
    double                      prop2_time_multiplier;
    std::vector<double>         coord;
    std::vector<ScalarVariable > OP;
};



/**

copyright (c) 1999 by Dr. Stephan K. Matthaei & Stephen G. Roberts */


} // csmp

#endif
















