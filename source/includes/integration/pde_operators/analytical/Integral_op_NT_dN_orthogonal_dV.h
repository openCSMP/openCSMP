#ifndef Integral_op_NT_dN_orthogonal_dV_h
#define Integral_op_NT_dN_orthogonal_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/// streamfunction operand for the RHS
template<uint32_t dim,class CELL=Element<dim> >
class Integral_op_NT_dN_orthogonal_dV : public MathOperatorRHS<dim> {
  public:
    Integral_op_NT_dN_orthogonal_dV( const PropertyDatabase<dim>& pref, 
                                     const char* basic,            // e.g., fluid pressure
                                     const char* test );           // streamfunction
    
    virtual void GetOperands( const CELL& );
    virtual void ComputeContribution( const CELL& );

  private:
    DenseMatrix<DM_MIN>  M, DNORTHO, NT; 
    std::vector<double>         NPROP, IPOL, UNITY, RES;
    const double                zero, one;
};

} // csmp

#endif

























