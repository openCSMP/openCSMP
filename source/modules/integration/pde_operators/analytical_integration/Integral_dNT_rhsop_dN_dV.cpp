// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "Integral_dNT_rhsop_dN_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
Integral_dNT_rhsop_dN_dV<dim,CELL>::Integral_dNT_rhsop_dN_dV( const PropertyDatabase<dim>& pref,
                                                              const char*             oper,
                                                              const char*             test )
  : MathOperatorRHS<dim,CELL>(pref,oper,test)
{
    MathOperatorRHS<dim,CELL>::Name("Integral_dNT_rhsop_dN_dV", oper, test );

    // testing the Operands 
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != ELEMENT &&
         MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() != REGION )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_dV::(constructor)",
                    oper, "Operand must be placed on the element or group.");

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE || MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
    throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_dV::(constructor)", 
                    test, "Testfunction (dependent) variable must be a scalar property placed on the nodes.");
}




template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_rhsop_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    assert( e.UsesLocalCoordinates() == false );

    MathOperatorRHS<dim,CELL>::MTRL.resize( 1U );
    MathOperatorRHS<dim,CELL>::MTRL[0].Resize( dim, dim );
    MathOperatorRHS<dim,CELL>::MTRL[0].Zero();

    switch ( MathOperatorRHS<dim,CELL>::MaterialOperandType() )
      {
        case SCALAR:
          MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonalAndZeroOffDiagonal(
              dim, e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey() ) );
          break;
        case VECTOR:
          {
            VectorVariable<dim> vc;
            e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), vc );
            MathOperatorRHS<dim,CELL>::MTRL[0].AssignToDiagonal( vc );
          }
          break;
        case TENSOR:
          {
            TensorVariable<dim> ts;
            e.Read( MathOperatorRHS<dim,CELL>::MaterialOperandKey(), ts );
            MathOperatorRHS<dim,CELL>::MTRL[0] = ts;
          }
          break;
        default:
          throw csmp::Exception( ERROR, "Integral_dNT_rhsop_dN_dV::GetOperands",
                                 "Unsupported material operand type." );
      }

    // read nodal test variable values
    e.NodePropertyVector( MathOperatorRHS<dim,CELL>::TestOperandKey(), u_ );

     e.NodePropertyVector( MathOperatorRHS<dim,CELL>::TestOperandKey(), u_ );

} // end GetOperands




template<uint32_t dim, template<uint32_t> class CELL>
void Integral_dNT_rhsop_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    assert( e.UsesLocalCoordinates() == false );

    MathOperatorRHS<dim,CELL>::RHS.resize( e.Nodes() );

    e.dN( DN_ );
    DN_.Transposed( DNT_ );

    TEMP_ = DNT_;
    TEMP_ *= MathOperatorRHS<dim,CELL>::MTRL[0];
    TEMP_ *= DN_;

    const double volume{ e.Volume() };
    for ( uint32_t j{0U}; j < e.Nodes(); ++j )
      {
        double val{0.0};
        for ( uint32_t k{0U}; k < e.Nodes(); ++k )
          val += TEMP_(j,k) * u_[k]();
        MathOperatorRHS<dim,CELL>::RHS[j] = val * volume;
      }

 } // end ComputeContribution



template class Integral_dNT_rhsop_dN_dV<1U,Element>;
template class Integral_dNT_rhsop_dN_dV<2U,Element>;
template class Integral_dNT_rhsop_dN_dV<3U,Element>;

template class Integral_dNT_rhsop_dN_dV<1U,Face>;
template class Integral_dNT_rhsop_dN_dV<2U,Face>;
template class Integral_dNT_rhsop_dN_dV<3U,Face>;

} // csmp
