#ifndef NUM_INTEGRAL_DNT_OP_GRAVITY_DN_DN_H
#define NUM_INTEGRAL_DNT_OP_GRAVITY_DN_DN_H

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// PDE operator:  by Darcy's law substitution into Laplacian
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_dNT_op_gravity_dN_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_dNT_op_gravity_dN_dV( const PropertyDatabase<dim>& pref, 
                                      const char* oper, ///< e.g. conductivity
                                      const char* body_force_variable,  ///< like fluid density
                                      const char* basic, 
                                      const char* test,
                                      double64 acc_gravity );
  
    virtual void GetOperands( SIMPLEX& element );
  
    virtual void ComputeContribution( SIMPLEX& element );
  
    virtual NumIntegral_dNT_op_gravity_dN_dV<dim,SIMPLEX >* clone() const { return new NumIntegral_dNT_op_gravity_dN_dV<dim,SIMPLEX >(*this); }
    
  protected:
    DenseMatrix<DM_MIN>          B, BT;
    const double64               acc_gravity_;
    csmp::Index                  g_key_;
    std::vector<ScalarVariable>  g_property_;
};

} // csmp

#endif /* NUM_INTEGRAL_DNT_OP_GRAVITY_DN_DN_H */
