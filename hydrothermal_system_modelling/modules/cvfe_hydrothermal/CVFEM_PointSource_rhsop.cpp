// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVFEM_PointSource_rhsop.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_PointSource_rhsop<dim,CELL>::~CVFEM_PointSource_rhsop() {}

template<uint32_t dim, template<uint32_t> class CELL>
CVFEM_PointSource_rhsop<dim,CELL>::CVFEM_PointSource_rhsop( const PropertyDatabase<dim>& pref,
                                                            const char* oper,
                                                            const char* test)
    : CVFEM_MathOperatorRHS<dim>(pref,oper,test),
      SRC_(3)
{
    MathOperatorRHS<dim,CELL>::Name("CVFEM_PointSource_rhsop", oper, test );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR ||
         MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != NODE )
    {
        cout <<"\nMathOperatorRHS->CVFEM_PointSource_rhsop<dim,CELL>::(constructor): ";
        cout <<"Fatal Error: Operand must be a scalar variable placed on the node. Terminating..."<< endl;
        throw invalid_argument("CVFEM_PointSource_rhsop<csp_float,dim>::CVFEM_PointSource_rhsop");
    }
}


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_PointSource_rhsop<dim,CELL>::GetOperands( const CELL<dim>& e )
{
    e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), SRC_ );
    // taking into account that the source contributes to several elements
    for ( uint32_t i{0}; i<e.Nodes(); i++ )
        SRC_[i] /= static_cast<double>(e.N(i)->Parents());

} // end GetOperands


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_PointSource_rhsop<dim,CELL>::GetOperandsCVFEM( const CELL<dim>& e, csmp::Index upwind_var_key  )
{
    e.NodePropertyVector( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), SRC_ );
    e.NodePropertyVector( upwind_var_key, upwind_var_);

    // taking into account that the source contributes to several elements
    for ( uint32_t i{0}; i<e.Nodes(); i++ )
        SRC_[i] /= static_cast<double>(e.N(i)->Parents());

    for ( uint32_t i = 0; i < upwind_var_.size(); i++ )
        SRC_[i]() *= upwind_var_[i]();

} // end GetOperands


template<uint32_t dim, template<uint32_t> class CELL>
void CVFEM_PointSource_rhsop<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorRHS<dim,CELL>::RHS.resize( e.Nodes() );
    for ( auto i{0U}; i<e.Nodes(); i++ )
        MathOperatorRHS<dim,CELL>::RHS[i] = SRC_[i]();

} // end ComputeContribution

template class CVFEM_PointSource_rhsop<1U,Element>;
template class CVFEM_PointSource_rhsop<2U,Element>;
template class CVFEM_PointSource_rhsop<3U,Element>;

} // csmp
