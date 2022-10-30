#ifndef NumIntegral_BT_D_op_dV_h
#define NumIntegral_BT_D_op_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {

template<uint32_t> class Element;

/**

Vector solution variable: volume strain

@author S.K. Matthaei
@author S. Roberts
@date 1999 */

template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_BT_D_op_dV : public MathOperatorRHS<dim,CELL> {
  public:
    NumIntegral_BT_D_op_dV( const PropertyDatabase<dim>&, const char* oper,     // volume strain
                            const char* youngs,   // Young's modulus
                            const char* poissons, // Poisson's ratio
                            const char* test, bool plane_strain=true ); // displacement

    virtual void GetOperands( const CELL<dim>& );

    virtual void ComputeContribution( const CELL<dim>& );
    
    void PlaneStress();
    
    virtual NumIntegral_BT_D_op_dV<dim,CELL>* clone() const { return new NumIntegral_BT_D_op_dV<dim,CELL> (*this); }
    
  private:  
    csmp::Index                       nu_key_;    ///< Poisson's ratio
    csmp::Index                       Y_key_;     ///< Young's modulus
    std::vector<DenseMatrix<DM_MIN> > E_;         ///< Young's modulus
    std::vector<DenseMatrix<DM_MIN> > nu_;        ///< variable in which Poisson's ratio will be stored
    bool                              plane_strain_;
    DenseMatrix<DM_MIN>               D, B, BT, TEMP, STR;  // material property matrix, stiffness matrix
    ScalarVariable                    sc_;
    VectorVariable<dim>               vc_;
    TensorVariable<dim>               ts_;
};



} // csmp

#endif
