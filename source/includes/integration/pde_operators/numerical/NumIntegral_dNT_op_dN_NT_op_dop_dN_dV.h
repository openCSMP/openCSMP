#ifndef NumIntegral_dNT_op_dN_NT_op_dop_dN_dV_h
#define NumIntegral_dNT_op_dN_NT_op_dop_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"

namespace csmp {

/// advection-dispersion matrices @note v-term is calculated from 'grad' dop x op multiplier
template<size_t dim,class CELL=Element<dim> >
class NumIntegral_dNT_op_dN_NT_op_dop_dN_dV : public MathOperatorLHS<dim> {
  public:
    NumIntegral_dNT_op_dN_NT_op_dop_dN_dV( const PropertyDatabase<dim>& pref, 
                                          const char* grad_prop,        // e.g., for calc. of v
                                          const char* eprop,            // e.g., property for multiplication with grad
                                          const char* emultiplier,      // e.g., property for multiplication with eprop
                                          const char* oper,             // e.g., conductivity
                                          const char* basic,            // e.g., fluid pressure
                                          const char* test );           // e.g., fluid pressure
                        
    NumIntegral_dNT_op_dN_NT_op_dop_dN_dV( const PropertyDatabase<dim>& pref, 
                                           const char* grad_prop,       // e.g., fluid pressure
                                           const char* eprop,           // e.g., conductivity
                                           const char* emultiplier,     // e.g., property for multiplication with eprop
                                           const char* rrho_prop,       // e.g., for calc. gravity effects
                                           const char* oper,            // e.g., thermal conductivity
                                           const char* basic,           // e.g., temperature
                                           const char* test );          // e.g., temperature
    
    virtual void GetOperands( CELL& e );
    virtual void ComputeContribution( CELL& e );
    
    void SpatialDerivative( size_t num_xyz ); // set gradZ direction to X=1, Y=2, Z=3

  private:
    DenseMatrix<DM_MIN>  DN, DNT, // derivatives of basis functions
                                VIP;     // 'v' (velocity) vector<double> variable
    std::vector<double>             IPOL;    // basis function values (at integration point)
    std::vector<DenseMatrix<DM_MIN> >  NT3;     // IPOL at integration points stored in columns of matrix NTNTNT
    std::vector<double>             NGRAD,   // variable to compute gradient of for calculation of 'v'
                                RDENS;   // relative density if so specified
                    
    std::vector<ScalarVariable >  sc_prop_vec, // vector<double> to read node properties into
                                      rrho_vec; 
    ScalarVariable                emult, econd; 
         
    csmp::Index       grad_key,    // key to variable to compute gradient of for calculation of 'v'
                    rrho_key,    // relative fluid density at nodes or similar operand
                    cond_key,    // e.g., op to be multiplied with dop
                    mult_key;    // for instance heat transport velocity coefficient
    bool            with_gravity;
    const double        gravity;
    size_t       xyz;
};

} // csmp

#endif

























