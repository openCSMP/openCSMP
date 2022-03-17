#include "NumIntegral_NT_dNi_dV_sc.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<uint32_t dim,class CELL>
NumIntegral_NT_dNi_dV_sc<dim,CELL>::~NumIntegral_NT_dNi_dV_sc() {}

//constructor
template<uint32_t dim,class CELL>
NumIntegral_NT_dNi_dV_sc<dim,CELL>::NumIntegral_NT_dNi_dV_sc( const PropertyDatabase<dim>& pref,
                                                                 const char*                  oper,
                                                                 const char*                  basic,
                                                                 const char*                  test )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    DN(3,3), TEMP(3,3), IPOL(3), xyz_(Y_DIRECTION), transp_(false)
{
    MathOperatorLHS<dim>::Name("NumIntegral_NT_dNi_dV_sc",oper, basic, test );
    

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV_sc<dim>::(constructor)",
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim>::TestOperandType() != SCALAR)
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV_sc<dim>::(constructor)",
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}

template<uint32_t dim,class CELL>
void NumIntegral_NT_dNi_dV_sc<dim,CELL>::SpatialDerivative( SPATIAL_DERIVATIVE num_xyz )
 {
    xyz_ = num_xyz;
 }
 
template<uint32_t dim,class CELL>
void NumIntegral_NT_dNi_dV_sc<dim,CELL>::Transposed()
 {
    transp_ = true;
 }
 
//element contribution
template<uint32_t dim,class CELL>
void NumIntegral_NT_dNi_dV_sc<dim,CELL>::ComputeContribution( const CELL& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();
    double det;
    IPOL.resize( e.Nodes() );
    TEMP.Resize( e.Nodes(), e.Nodes() );
    for ( auto i=0; i<e.IntegrationPoints(); i++ )
      {
         e.N_AtIntegrationPoint( i, IPOL );
         det = e.dN_AtIntegrationPoint( DN, i );
         //det = e.det_JINV_AtIntegrationPoint( i );
         for ( auto j=0; j<e.Nodes(); j++ )
            for ( auto k=0; k<e.Nodes(); k++ )
              transp_ ? TEMP(j,k) = IPOL[j] * DN(xyz_,k) : TEMP(j,k) = IPOL[k] * DN(xyz_,j);
         TEMP *= (det * e.WeightAtIntegrationPoint(i));
         MathOperatorLHS<dim>::LHS += TEMP;
      }

} // end ComputeContribution

template class NumIntegral_NT_dNi_dV_sc<1U,Element<1U> >;
template class NumIntegral_NT_dNi_dV_sc<2U,Element<2U> >;
template class NumIntegral_NT_dNi_dV_sc<3U,Element<3U> >;

template class NumIntegral_NT_dNi_dV_sc<1U,Face<1U> >;
template class NumIntegral_NT_dNi_dV_sc<2U,Face<2U> >;
template class NumIntegral_NT_dNi_dV_sc<3U,Face<3U> >;

} //namespace csp
