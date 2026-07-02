#ifndef NumJacobianIntegral_BT_mN_dV_h
#define NumJacobianIntegral_BT_mN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {
/**
 * @brief Numerically integrated Biot coupling matrix.
 *
 * Computes the element contribution to the coupling block Q in the
 * block-structured poromechanics system:
 *
 *   [ K   Q ] { u }   { f_u }
 *   [ Q^T S ] { p } = { f_p }
 *
 * The element matrix is:
 *
 *   Q^e = alpha * integral( B_u^T * m * N_p ) dV
 *
 * where:
 *   B_u  : strain-displacement matrix        (3 x 2*nodes)
 *   m    : volumetric coupling vector [1,1,0]^T
 *   N_p  : pressure shape function row vector (1 x nodes)
 *   alpha: Biot-Willis coefficient (scalar, element or integration point)
 *
 * Output matrix size: (2*nodes x nodes)
 *
 * The rows correspond to displacement DOFs (interleaved x1,y1,x2,y2,...),
 * the columns correspond to pressure DOFs (one per node).
 
   @attention USE SUBTRACT_ACCUMULATE for pore pressure and stress to get effective stress
 **/
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumJacobianIntegral_BT_m_N_dV : public MathOperatorLHS<dim, CELL>
{
public:
    NumJacobianIntegral_BT_m_N_dV( const PropertyDatabase<dim>& pref,
                                   const char* oper,   // Biot coefficient alpha
                                   const char* basic,  // pore pressure (scalar, NODE)
                                   const char* test ); // displacement (vector, NODE)

    void GetOperands( const CELL<dim>& e ) override;
    void ComputeContribution( const CELL<dim>& e ) override;

    NumJacobianIntegral_BT_m_N_dV<dim,CELL>* clone() const override { return new NumJacobianIntegral_BT_m_N_dV<dim,CELL> (*this); }

private:
    // Biot coefficient storage (element or per integration point)
    std::vector<double> alpha_;   // MTRL

    // Working matrices
    DenseMatrix<DM_MIN>  B_;      // strain-displacement:        3  x 2*nodes
    DenseMatrix<DM_MIN>  BT_;     // transposed:                 2*nodes x 3
    std::vector<double>  Np_;     // pressure shape functions:   1  x nodes
    DenseMatrix<DM_MIN>  BTm_;    // B^T * m:                    2*nodes x 1
    // Compile-time evaluation of the number of strains
    static constexpr uint32_t n_strains_ = (dim == 2) ? 3 : 6;
    
    // Static compile-time constant volumetric coupling vector
    static constexpr std::array<double, n_strains_> m_ = []()
    {
        static_assert( dim == 2 || dim == 3,
            "NumJacobianIntegral_BT_m_N_dV: Only dim=2 and dim=3 are supported." );

        std::array<double, n_strains_> arr{};  // zero-initialise shear components

        if constexpr ( dim == 2 ) {
            arr[0] = 1.0;  // epsilon_xx  — compression positive
            arr[1] = 1.0;  // epsilon_yy  — compression positive
            // arr[2] = 0.0  gamma_xy does not contribute to volumetric strain
        }
        else if constexpr ( dim == 3 ) {
            arr[0] = 1.0;  // epsilon_xx
            arr[1] = 1.0;  // epsilon_yy
            arr[2] = 1.0;  // epsilon_zz
            // arr[3..5] = 0.0  shear components
        }
        return arr;
    }();
 };

} // csmp

#endif
