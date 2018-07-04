#ifndef NumIntegral_op_NT_dN_orthogonal_dV_h
#define NumIntegral_op_NT_dN_orthogonal_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"
#include "Operand.h"

namespace csmp {

/// streamfunction operand for the RHS
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_op_NT_dN_orthogonal_dV : public MathOperatorRHS<dim> {
  public:
    NumIntegral_op_NT_dN_orthogonal_dV( const PropertyDatabase<dim>& pref, 
                                        const char* basic,            // e.g., fluid pressure
                                        const char* test );           // streamfunction
    
    virtual void GetOperands( CELL& e );
    virtual void ComputeContribution( CELL& e );
    virtual NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL>* clone() const { return new NumIntegral_op_NT_dN_orthogonal_dV<dim,CELL> (*this); }

  private:
    DenseMatrix<DM_MIN>  M, DNORTHO, NT; 
    std::vector<double64>         NPROP, IPOL, UNITY, RES;
};

} // csmp

#endif

























