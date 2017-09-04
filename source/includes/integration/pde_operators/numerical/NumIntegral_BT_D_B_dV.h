#ifndef NumIntegral_BT_D_B_dV_h
#define NumIntegral_BT_D_B_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/**
    vector solution variable: stiffness matrix with 2(2D) or 3(3D) nodal degrees of freedon
 
    @attention special case: if nu=0.5 a purely viscous (incompressible) fluid is modeled
    this requires a special material property matrix, see Zienkiewicz Vol II
*/
template<size_t dim,class SIMPLEX=Element<dim> >
class NumIntegral_BT_D_B_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_BT_D_B_dV( const PropertyDatabase<dim>& pref, 
                           const char* oper1, const char* oper2, 
                           const char* basic, const char* test,
                           bool plane_strain=true );
    
    virtual void GetOperands( SIMPLEX& e );
    
    virtual void ComputeContribution( SIMPLEX& e );
    
    void PlaneStress( bool yes_no=true ); ///< default is plane strain
    virtual NumIntegral_BT_D_B_dV<dim,SIMPLEX>* clone() const { return new NumIntegral_BT_D_B_dV<dim,SIMPLEX> (*this); }
  private:
    csmp::Index            nu_key_;   ///< Poisson's ratio
    std::vector<double64>  E_, nu_;   ///< variable in which Poisson's ratio will be stored
    DenseMatrix<DM_MIN>    D, B, BT;  ///< material property matrix
    bool                   plane_strain_;
};


} // csmp

#endif
