#include "NumIntegral_PT_lhsop_P_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class CELL>
NumIntegral_PT_lhsop_P_dV<dim,CELL>::NumIntegral_PT_lhsop_P_dV( const PropertyDatabase<dim>& pref,
                                                                const char* oper,  // density
                                                                const char* oper2, // porosity
                                                                const char* basic, // interstitial velocity
                                                                const char* test ) // interstitial velocity
  : MathOperatorLHS<dim,CELL>(pref,oper,basic,test),
    phi_key(pref.StorageKey(oper2)),
    nodal_degrees_of_freedom(dim),
    PT(3*2,1), 
    P(1,3*2)
 {
    MathOperatorLHS<dim,CELL>::Name("NumIntegral_PT_lhsop_P_dV", oper, basic, test );
    
    // resize material property matrix 
    if constexpr ( dim == 3U )
      {
         for ( auto it=MathOperatorLHS<dim,CELL>::MTRL.begin();
               it!=MathOperatorLHS<dim,CELL>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }

    if ( phi_key.type != SCALAR && phi_key.place != ELEMENT ) 
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                             oper2, "Must be a scalar variable placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                             oper, "Operand must be of scalar type." );

    if ( MathOperatorLHS<dim,CELL>::BasicOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::BasicOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a vector<double> property placed on the nodes." );

    if ( MathOperatorLHS<dim,CELL>::TestOperandPlacement() != NODE ||
         MathOperatorLHS<dim,CELL>::TestOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a vector<double> property placed on the nodes." );
}




template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_PT_lhsop_P_dV<dim,CELL>::GetOperands( const CELL<dim>& e )
 {
    // this integral is only for numerically integrated isoparametric finite elements
    assert( e.FE()->Isoparametric() == true );

    if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
      {
         MathOperatorLHS<dim,CELL>::MTRL[0].Resize(dim,dim);
         MathOperatorLHS<dim,CELL>::MTRL[0].Zero();
      
         // reading density but using the porosity variable temporarily
         if ( MathOperatorLHS<dim,CELL>::MaterialOperandType() == SCALAR ) {
              e.Read( MathOperatorLHS<dim,CELL>::MaterialOperandKey(), phi );
              for ( auto i{0U}; i<dim; i++ ) 
                MathOperatorLHS<dim,CELL>::MTRL[0](i,i) = phi();
           }
      }
    else for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ )
           MathOperatorLHS<dim,CELL>::PropertyAtIntegrationPoint( e, MathOperatorLHS<dim,CELL>::MaterialOperandKey(),
                                                                  i, MathOperatorLHS<dim,CELL>::MTRL[i] );

   // reading the porosity which must be an element property
   e.Read( phi_key, phi );
   
} // end GetOperands()





template<uint32_t dim, template<uint32_t> class CELL>
void  NumIntegral_PT_lhsop_P_dV<dim,CELL>::N_to_P( const std::vector<double>& N, DenseMatrix<DM_MIN>& mP )
 {
    mP.Resize(1,nodal_degrees_of_freedom * N.size());

    int k{0};
    for ( auto i{0U}; i<N.size(); i++ ) 
      for ( auto j{0U}; j<nodal_degrees_of_freedom; j++ ) mP(0,k++ ) = N[i];
}




/**
 
Computes the volume integral over the element interpolation functions 
times the Operand. If the Operand is 1 over the element, then the volume
integral is naturally 1 as well.  
*/
template<uint32_t dim, template<uint32_t> class CELL>
void NumIntegral_PT_lhsop_P_dV<dim,CELL>::ComputeContribution( const CELL<dim>& e )
{
    MathOperatorLHS<dim,CELL>::LHS.Resize(e.Nodes()*dim,e.Nodes()*dim);
    MathOperatorLHS<dim,CELL>::LHS.Zero();
    
    vector<double>  N( e.Nodes() );
    double          det( 0.0 );
    
    for ( auto i{0U}; i<e.FE()->IntegrationPoints(); i++ )
      {
         e.N_AtIntegrationPoint( i, N );
         det = e.det_JINV_AtIntegrationPoint( i );

         // getting the transposed of P
         N_to_P( N, P );
         P.Transposed( PT );
         
         for ( auto j{0U}; j<P.Cols(); j++ )
           if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == ELEMENT )
             P(0,j) *= MathOperatorLHS<dim,CELL>::MTRL[0](0,0) * phi();
           else if ( MathOperatorLHS<dim,CELL>::MaterialOperandPlacement() == NODE )
             P(0,j) *= MathOperatorLHS<dim,CELL>::MTRL[i](0,0) * phi();
              
         PT   *= P;
         PT   *= (det * e.WeightAtIntegrationPoint(i));
        
         MathOperatorLHS<dim,CELL>::LHS += PT;
      }

    cout <<"\nNumIntegral_PT_lhsop_P_dV<dim>::N_to_P: LHS contribution element: "<< e.Idx() << endl;    
    MathOperatorLHS<dim,CELL>::LHS.Out();
      
} // end ComputeContribution


template class NumIntegral_PT_lhsop_P_dV<1U,Element>;
template class NumIntegral_PT_lhsop_P_dV<2U,Element>;
template class NumIntegral_PT_lhsop_P_dV<3U,Element>;

template class NumIntegral_PT_lhsop_P_dV<1U,Face>;
template class NumIntegral_PT_lhsop_P_dV<2U,Face>;
template class NumIntegral_PT_lhsop_P_dV<3U,Face>;

} // csmp
