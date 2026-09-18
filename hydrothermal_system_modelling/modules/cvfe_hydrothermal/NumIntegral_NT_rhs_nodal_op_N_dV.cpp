// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PropertyDatabase.h"
#include "NumIntegral_NT_rhs_nodal_op_N_dV.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "Exception.h"
#include <iostream>
#include <chrono>

using namespace std;

namespace csmp {

/** default constructor */
//template<uint32_t dim, template<uint32_t> class CELL>
//NumIntegral_NT_rhs_nodal_op_N_dV<dim,CELL>::NumIntegral_NT_rhs_nodal_op_N_dV() {}


/** default destructor */
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_rhs_nodal_op_N_dV<dim,CELL>::~NumIntegral_NT_rhs_nodal_op_N_dV() {}

/** custom constructor */
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_rhs_nodal_op_N_dV<dim,CELL>::NumIntegral_NT_rhs_nodal_op_N_dV( const PropertyDatabase<dim>& pref,
                                                                              const char* oper,
                                                                              const char* test,
                                                                              const char* thickness)//Benoit 2025 add
    : MathOperatorRHS<dim,CELL>(pref,oper,test),
    thickness_(pref.StorageKey(thickness))//Benoit 2025 add

{
    MathOperatorRHS<dim,CELL>::Name("NumIntegral_NT_rhs_nodal_op_N_dV", oper, test );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "NumIntegral_NT_rhs_nodal_op_N_dV<dim>::(constructor)",
                               oper, "Operand must be a scalar property." );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != NODE )
        throw csmp::Exception( ERROR, "NumIntegral_NT_rhs_nodal_op_N_dV<dim>::(constructor)",
                               oper, "operand has only been implemented for nodal variables." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
        throw csmp::Exception( ERROR, "NumIntegral_NT_rhs_nodal_op_N_dV<dim>::(constructor)",
                               test, "Dependent variable must be a scalar property placed on the nodes." );
}


/** read in nodal vector of material operand */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_rhs_nodal_op_N_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), VAR );

} // end GetOperands


/** Calculate contributions to the capacitance matrix */
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_rhs_nodal_op_N_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize(e.Nodes());
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(), MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );

    volume = e.Volume();
    volume*=e.Read(thickness_);//Benoit 2025 add

    for ( size_t j{0U}; j<e.Nodes(); j++ )
        MathOperatorRHS<dim,CELL>::RHS[j] = (VAR[j]()*volume) / static_cast<double>(e.Nodes());

} // end ComputeContribution


template class NumIntegral_NT_rhs_nodal_op_N_dV<1U>;
template class NumIntegral_NT_rhs_nodal_op_N_dV<2U>;
template class NumIntegral_NT_rhs_nodal_op_N_dV<3U>;

template class NumIntegral_NT_rhs_nodal_op_N_dV<1U,Face>;
template class NumIntegral_NT_rhs_nodal_op_N_dV<2U,Face>;
template class NumIntegral_NT_rhs_nodal_op_N_dV<3U,Face>;

} // csmp
