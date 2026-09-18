// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV_h
#define NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV_h

#include "CSMP_definitions.h"
#include "MathOperatorLHS.h"
#include "Operand.h"
#include "DenseMatrix.h"

namespace csmp {

template<uint32_t> class Element;

/// advection-dispersion matrices @note v-term is calculated from 'grad' test operand x multiplier
template<uint32_t dim, template<uint32_t> class CELL=Element>
class NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV : public MathOperatorLHS<dim,CELL> {
  public:
    NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>&,
                                          const char* emultiplier,      // element v multiplier
                                          const char* nmultiplier,      // nodal v multiplier
                                          const char* grad_prop,        // e.g., for calc. of v
                                          const char* oper,             // e.g., conductivity
                                          const char* basic,            // e.g., fluid pressure
                                          const char* test );           // e.g., fluid pressure
                        
    NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV( const PropertyDatabase<dim>&, 
                                          const char* emultiplier,      // element v multiplier
                                          const char* nmultiplier,      // nodal v multiplier
                                          const char* grad_prop,        // e.g., for calc. of v
                                          const char* rrho_prop,        // e.g., for calc. gravity effects
                                          const char* oper,             // e.g., conductivity
                                          const char* basic,            // e.g., fluid pressure
                                          const char* test );           // e.g., fluid pressure
    
    void GetOperands( const CELL<dim>& ) override final;
    void ComputeContribution( const CELL<dim>& ) override final;
    
    void SpatialDerivative( uint32_t num_xyz ); // set gradZ direction to X=1, Y=2, Z=3

    NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>* clone() const override final { return new NumIntegral_dNT_mixed_op_dN_NT_v_dN_dV<dim,CELL>(*this); }

  private:

    void ReadElementMultiplier( const CELL<dim>&,
                                DenseMatrix<DM_MIN>& MULT );
    
    DenseMatrix<DM_MIN>                DN, DNT, ///< derivatives of basis functions
                                       VIP,     ///< 'v' (velocity) vector<double> variable
                                       EMULT;   ///< element based multiplier for 'v'
    std::vector<double>                IPOL;    ///< basis function values (at integration point)
    std::vector<DenseMatrix<DM_MIN> >  NT3; ///< IPOL at integration points stored in columns of matrix NTNTNT
    std::vector<double>                NMULT,   ///< node based scalar multiplier for 'v'
                                       NGRAD,   ///< variable to compute gradient of for calculation of 'v'
                                       RDENS;   ///< relative density if so specified
                    
    std::vector<ScalarVariable >  sc_prop_vec, ///< vector<double> to read node properties into
                                  rrho_vec;
         
    csmp::Index     grad_key,    ///< key to variable to compute gradient of for calculation of 'v'
                    emulti_key,  ///< key to element based multiplier for 'v'
                    nmulti_key,  ///< key to node based scalar multiplier for 'v'
                    rrho_key;    ///< relative fluid density at nodes or similar operand
    bool            with_gravity;
    const double    gravity;
    uint32_t        xyz;
};

}

#endif
