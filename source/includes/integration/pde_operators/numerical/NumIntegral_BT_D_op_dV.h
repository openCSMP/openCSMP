#ifndef NumIntegral_BT_D_op_dV_h
#define NumIntegral_BT_D_op_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorRHS.h"

namespace csmp {
/**
@author S.K. Matthaei
@author S. Roberts
@date 1999 */

/// vector solution variable: volume strain
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_BT_D_op_dV : public MathOperatorRHS<dim> {
  public:
    //                                                                          for instance:                                           
    NumIntegral_BT_D_op_dV( const PropertyDatabase<dim>& pref, const char* oper,     // volume strain
                                                          const char* youngs,   // Young's modulus
                                                          const char* poissons, // Poisson's ratio
                                                          const char* test, bool plane_strain=true );   // displacement

    virtual void GetOperands( CELL& e );

    virtual void ComputeContribution( CELL& e );
    
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
