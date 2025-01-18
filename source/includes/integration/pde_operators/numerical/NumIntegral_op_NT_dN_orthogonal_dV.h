#ifndef NumIntegral_op_NT_dN_orthogonal_dV_h
#define NumIntegral_op_NT_dN_orthogonal_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

template<uint32_t> class Element;

/// streamfunction operand for the RHS
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_op_NT_dN_orthogonal_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_op_NT_dN_orthogonal_dV( const PropertyDatabase<dim>&,
                                        const char* basic,            // e.g., fluid pressure
                                        const char* test );           // streamfunction
    
    virtual void GetOperands( const CELL<dim>& );
    virtual void ComputeContribution( const CELL<dim>& );
    
    virtual NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL>* clone() const { return new NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN>  M, DNORTHO, NT; 
    std::vector<double>  NPROP, IPOL, UNITY, RES;
};

} // csmp

#endif

























