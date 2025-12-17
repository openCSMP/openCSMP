#ifndef Integral_op_NT_dN_orthogonal_dV_h
#define Integral_op_NT_dN_orthogonal_dV_h

#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/// streamfunction operand for the RHS
template<uint32_t dim, template<uint32_t> class CELL=Element>
class Integral_op_NT_dN_orthogonal_dV : public MathOperatorRHS<dim,CELL> {
  public:
    Integral_op_NT_dN_orthogonal_dV( const PropertyDatabase<dim>&,
                                     const char* basic,            // e.g., fluid pressure
                                     const char* test );           // streamfunction
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;

  private:
    DenseMatrix<DM_MIN>  M, DNORTHO, NT; 
    std::vector<double>         NPROP, IPOL, UNITY, RES;
    const double                zero, one;
};

} // csmp

#endif

























