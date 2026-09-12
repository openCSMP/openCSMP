#include "NumIntegral_DNT_rhsop_DN_dV.h"
#include "PropertyDatabase.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

// 2-argument constructor — no multiplier
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_rhsop_dN_dV<dim,CELL>::NumIntegral_dNT_rhsop_dN_dV(
    const PropertyDatabase<dim>& pref,
    const char* oper,
    const char* test )
  : MathOperatorRHS<dim,CELL>( pref, oper, test ),
    has_multiplier_( false ),
    multiplier( PLAIN, 1.0 ),
    DN( dim, 3 ),
    DNT( 3, dim ),
    OPMAT( dim, 1 )
 {
    MathOperatorRHS<dim,CELL>::Name( "NumIntegral_dNT_rhsop_dN_dV", oper, test );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE ||
         MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_rhsop_dN_dV<dim>::(constructor)",
                             oper, "Material operand must be a scalar property." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_rhsop_dN_dV<dim>::(constructor)",
                             test, "Test operand must be a scalar property placed on the nodes." );
 }




// 3-argument constructor — with multiplier
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_dNT_rhsop_dN_dV<dim,CELL>::NumIntegral_dNT_rhsop_dN_dV(
    const PropertyDatabase<dim>& pref,
    const char* integral_multiplier,
    const char* oper,
    const char* test )
  : MathOperatorRHS<dim,CELL>( pref, oper, test ),
    has_multiplier_( true ),
    mult_key( pref.StorageKey( integral_multiplier ) ),
    multiplier( PLAIN, 1.0 ),
    DN( dim, 3 ),
    DNT( 3, dim ),
    OPMAT( dim, 1 )
 {
    MathOperatorRHS<dim,CELL>::Name( "NumIntegral_dNT_rhsop_dN_dV", oper, test );

    if ( mult_key.type != SCALAR || mult_key.place != ELEMENT )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_rhsop_dN_dV<dim>::(constructor)",
                             integral_multiplier, "Multiplier must be a scalar element property." );

    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE ||
         MathOperatorRHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_rhsop_dN_dV<dim>::(constructor)",
                             oper, "Material operand must be a scalar property." );

    if ( MathOperatorRHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorRHS<dim,CELL>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_dNT_rhsop_dN_dV<dim>::(constructor)",
                             test, "Test operand must be a scalar property placed on the nodes." );
 }





template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_rhsop_dN_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    assert( e.FE()->Isoparametric() == true );

    // 0. reading the integral multiplier — default to 1 if no multiplier was provided
    if ( has_multiplier_ )
      e.Read( mult_key, multiplier );
    else
      multiplier = ScalarVariable( PLAIN, 1.0 );

    // 1. reading the test variable (e.g. fluid pressure) from nodes
    e.NodePropertyVector( MathOperatorRHS<dim,CELL>::TestOperandKey(), noperand );

    // 2. reading the material operand into MTRL — same logic as MathOperatorLHS::GetOperands
    if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ||
         MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == FACE    ||
         MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == REGION  )
      {
        MathOperatorRHS<dim,CELL>::MTRL.resize( 1U );
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
              throw csmp::Exception( ERROR,
                "NumIntegral_dNT_rhsop_dN_dV::GetOperands",
                "Unsupported MaterialOperandType for ELEMENT/FACE/REGION placement." );
          }
      }
    else
      {
        const auto n_ip{ e.IntegrationPoints() };

        if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT_INTEGRATION_POINT ||
             MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == FACE_INTEGRATION_POINT )
          {
            MathOperatorRHS<dim,CELL>::MTRL.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              {
                switch ( MathOperatorRHS<dim,CELL>::MaterialOperandType() )
                  {
                    case SCALAR:
                      MathOperatorRHS<dim,CELL>::MTRL[i].AssignToDiagonalAndZeroOffDiagonal(
                          dim, e.Read( i, MathOperatorRHS<dim,CELL>::MaterialOperandKey() ) );
                      break;
                    case VECTOR:
                      {
                        VectorVariable<dim> vc;
                        e.Read( i, MathOperatorRHS<dim,CELL>::MaterialOperandKey(), vc );
                        MathOperatorRHS<dim,CELL>::MTRL[i].AssignToDiagonal( vc );
                      }
                      break;
                    case TENSOR:
                      {
                        TensorVariable<dim> ts;
                        e.Read( i, MathOperatorRHS<dim,CELL>::MaterialOperandKey(), ts );
                        MathOperatorRHS<dim,CELL>::MTRL[i] = ts;
                      }
                      break;
                    default:
                      throw csmp::Exception( ERROR,
                        "NumIntegral_dNT_rhsop_dN_dV::GetOperands",
                        "Unsupported MaterialOperandType for ELEMENT_INTEGRATION_POINT/FACE_INTEGRATION_POINT placement." );
                  }
              }
          }
        else if ( MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == NODE )
          {
            MathOperatorRHS<dim,CELL>::MTRL.resize( n_ip );
            for ( uint32_t i{0U}; i < n_ip; ++i )
              MathOperatorRHS<dim,CELL>::PropertyAtIntegrationPoint(
                  e, MathOperatorRHS<dim,CELL>::MaterialOperandKey(), i,
                  MathOperatorRHS<dim,CELL>::MTRL[i] );
          }
        else
          throw csmp::Exception( FATAL_ERROR,
            "NumIntegral_dNT_rhsop_dN_dV::GetOperands",
            "InterFace based operands cannot be accumulated with this method." );
      }

 } // end GetOperands





/**
 
@section arguments Input Arguments 

A reference to the finite-element from which the contribution is 
computed.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_dNT_rhsop_dN_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    MathOperatorRHS<dim,CELL>::RHS.resize( e.Nodes() );
    fill( MathOperatorRHS<dim,CELL>::RHS.begin(),
          MathOperatorRHS<dim,CELL>::RHS.end(), 0.0 );

    const bool piecewise_constant(
        MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT ||
        MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == FACE    ||
        MathOperatorRHS<dim,CELL>::MaterialOperandPlacement() == REGION  );

    for ( uint32_t i{0U}; i < e.FE()->IntegrationPoints(); ++i )
      {
        const double detJ = e.dN_AtIntegrationPoint( DN, i );
        DN.Transposed( DNT );

        // select material matrix at this integration point
        const DenseMatrix<dim>& D = piecewise_constant
                                       ? MathOperatorRHS<dim,CELL>::MTRL[0]
                                       : MathOperatorRHS<dim,CELL>::MTRL[i];

        // DNT(nodes x dim) * D(dim x dim) -> TEMP(nodes x dim)
        // TEMP(nodes x dim) * DN(dim x nodes) -> TEMP(nodes x nodes)
        // then apply nodal test values u_k: RHS[j] += sum_k TEMP(j,k) * u_k
        TEMP = DNT;
        TEMP *= D;
        TEMP *= DN;

        for ( uint32_t j{0U}; j < e.Nodes(); ++j )
          {
            double val{0.0};
            for ( uint32_t k{0U}; k < e.Nodes(); ++k )
              val += TEMP(j,k) * noperand[k]();
            MathOperatorRHS<dim,CELL>::RHS[j] +=
                val * e.WeightAtIntegrationPoint(i) * detJ * multiplier();
          }
      }

 } // end ComputeContribution



// cout <<"\nNumIntegral_dNT_rhsop_dN_dV<dim>::ComputeContribution: Element "<< e.Idx() <<":"<< endl; 
// nicePrint( RHS );

template class NumIntegral_dNT_rhsop_dN_dV<1U,Element>;
template class NumIntegral_dNT_rhsop_dN_dV<2U,Element>;
template class NumIntegral_dNT_rhsop_dN_dV<3U,Element>;

template class NumIntegral_dNT_rhsop_dN_dV<1U,Face>;
template class NumIntegral_dNT_rhsop_dN_dV<2U,Face>;
template class NumIntegral_dNT_rhsop_dN_dV<3U,Face>;

} // csmp
