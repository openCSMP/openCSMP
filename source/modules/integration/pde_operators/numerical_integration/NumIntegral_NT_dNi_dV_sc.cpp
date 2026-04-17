#include "NumIntegral_NT_dNi_dV_sc.h"
#include "Element.h"
#include "Face.h"
#include "Exception.h"

using namespace std;

namespace csmp {

//constructor
template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_NT_dNi_dV_sc<dim,CELL>::NumIntegral_NT_dNi_dV_sc( const PropertyDatabase<dim>& pref,
                                                              const char*                  oper,
                                                              const char*                  basic,
                                                              const char*                  test )
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    DN(3,3), TEMP(3,3), IPOL(3), xyz_(Y_DIRECTION), transp_(false)
{
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_NT_dNi_dV_sc",oper, basic, test );
    

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV_sc<dim>::(constructor)",
                             basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != SCALAR)
      throw csmp::Exception( ERROR, "NumIntegral_NT_dNi_dV_sc<dim>::(constructor)",
                             test, "Operand (test) must be a scalar property placed on the nodes." );
}



template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_dNi_dV_sc<dim,CELL>::SpatialDerivative( SPATIAL_DERIVATIVE num_xyz )
 {
    xyz_ = num_xyz;
 }
 



template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_dNi_dV_sc<dim,CELL>::Transposed()
 {
    transp_ = true;
 }
 




//element contribution
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_NT_dNi_dV_sc<dim,CELL>::ComputeContribution( const CELL<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    // initialize output matrix
    MathOperatorLHS<dim,CELL>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim,CELL>::LHS.Zero();
    double det;
    IPOL.resize( e.Nodes() );
    TEMP.Resize( e.Nodes(), e.Nodes() );
    for ( uint32_t i{0U}; i<e.IntegrationPoints(); i++ )
      {
         e.N_AtIntegrationPoint( i, IPOL );
         det = e.dN_AtIntegrationPoint( DN, i );
         //det = e.det_JINV_AtIntegrationPoint( i );
         for ( uint32_t j{0U}; j<e.Nodes(); j++ )
            for ( uint32_t k{0U}; k<e.Nodes(); k++ )
              transp_ ? TEMP(j,k) = IPOL[j] * DN(xyz_,k) : TEMP(j,k) = IPOL[k] * DN(xyz_,j);
         TEMP *= (det * e.WeightAtIntegrationPoint(i));
         MathOperatorLHS<dim,CELL>::LHS += TEMP;
      }

} // end ComputeContribution

template class NumIntegral_NT_dNi_dV_sc<1U,Element>;
template class NumIntegral_NT_dNi_dV_sc<2U,Element>;
template class NumIntegral_NT_dNi_dV_sc<3U,Element>;

template class NumIntegral_NT_dNi_dV_sc<1U,Face>;
template class NumIntegral_NT_dNi_dV_sc<2U,Face>;
template class NumIntegral_NT_dNi_dV_sc<3U,Face>;

} //namespace csp
