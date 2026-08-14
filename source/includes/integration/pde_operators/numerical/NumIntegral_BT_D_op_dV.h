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

    void GetOperands( const CELL<dim>& ) override final;

    void ComputeContribution( const CELL<dim>& ) override final;
    
    /// switch from default plane strain to plane stress
    void PlaneStress() noexcept { plane_strain_ = false; }
    
    NumIntegral_BT_D_op_dV<dim,CELL>* clone() const override final { return new NumIntegral_BT_D_op_dV<dim,CELL> (*this); }
    
  private:  
    csmp::Index         nu_key_; ///< Poisson's ratio
    csmp::Index         E_key_;  ///< Young's modulus
    std::vector<double> E_;      ///< Young's modulus
    std::vector<double> nu_;     ///< variable in which Poisson's ratio will be stored
    bool                plane_strain_;
    DenseMatrix<DM_MIN> B_;
    static constexpr auto MATDIM = (dim == 3) ? DM6 : DM3;
    std::vector<std::array<double,MATDIM>> initial_strains_;
    DenseMatrix<MATDIM> D_;
    ScalarVariable      sc_;
    VectorVariable<dim> vc_;
    TensorVariable<dim> ts_;
};



} // csmp

#endif
