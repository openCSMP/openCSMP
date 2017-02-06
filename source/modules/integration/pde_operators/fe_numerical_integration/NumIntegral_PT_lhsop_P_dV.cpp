#include "NumIntegral_PT_lhsop_P_dV.h"
#include "PropertyDatabase.h"
#include "Exception.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {


template<size_t dim,class SIMPLEX>
NumIntegral_PT_lhsop_P_dV<dim,SIMPLEX>::NumIntegral_PT_lhsop_P_dV( const PropertyDatabase<dim>& pref,
                                                      const char* oper,  // density
                                                      const char* oper2, // porosity
                                                      const char* basic, // interstitial velocity
                                                      const char* test ) // interstitial velocity
  : MathOperatorLHS<dim>(pref,oper,basic,test),    
    phi_key(pref.StorageKey(oper2)),
    nodal_degrees_of_freedom(dim),
    PT(3*2,1), 
    P(1,3*2)
 {
    MathOperatorLHS<dim>::Name("NumIntegral_PT_lhsop_P_dV", oper, basic, test );
    
    // resize material property matrix 
    if ( dim == 3 ) 
      {
         typename vector<DenseMatrix<DM_MIN> >::iterator  it;
         for ( it=MathOperatorLHS<dim>::MTRL.begin(); 
               it!=MathOperatorLHS<dim>::MTRL.end(); it++ ) (*it).Resize(3,3);
      }

    if ( phi_key.type != SCALAR && phi_key.place != ELEMENT ) 
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                      oper2, "Must be a scalar variable placed on the nodes." );

    if ( MathOperatorLHS<dim>::MaterialOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                      oper, "Operand must be of scalar type." );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a vector<double64> property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != VECTOR )
      throw csmp::Exception( ERROR, "NumIntegral_PT_lhsop_P_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a vector<double64> property placed on the nodes." );
}




template<size_t dim,class SIMPLEX>
void NumIntegral_PT_lhsop_P_dV<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
 {
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT ) 
      {
         MathOperatorLHS<dim>::MTRL[0].Resize(dim,dim);
         MathOperatorLHS<dim>::MTRL[0].Zero();
      
         // reading density but using the porosity variable temporarily
         if ( MathOperatorLHS<dim>::MaterialOperandType() == SCALAR ) {
              e.Read( MathOperatorLHS<dim>::MaterialOperandKey(), phi );
              for ( size_t i=0; i<dim; i++ ) 
                MathOperatorLHS<dim>::MTRL[0](i,i) = phi();
           }
      }
    else for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
           MathOperatorLHS<dim>::PropertyAtIntegrationPoint( e, MathOperatorLHS<dim>::MaterialOperandKey(), 
                                                                i, MathOperatorLHS<dim>::MTRL[i] );

   // reading the porosity which must be an element property
   e.Read( phi_key, phi );
   
} // end GetOperands()








/**
 
Computes the volume integral over the element interpolation functions 
times the Operand. If the Operand is 1 over the element, then the volume
integral is naturally 1 as well.  
*/
template<size_t dim,class SIMPLEX>
void NumIntegral_PT_lhsop_P_dV<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
    MathOperatorLHS<dim>::LHS.Resize(e.Nodes()*dim,e.Nodes()*dim);
    MathOperatorLHS<dim>::LHS.Zero();
    
    vector<double64>  N( e.Nodes() );
    double64          det( 0.0 );
    
    for ( size_t i=0; i<e.FE()->IntegrationPoints(); i++ )
      {
         e.N_AtIntegrationPoint( i, N );
         det = e.det_JINV_AtIntegrationPoint( i );

         // getting the transposed of P
         N_to_P( N, P );
         P.Transposed( PT );
         
         for ( size_t j=0; j<P.Cols(); j++ )
           if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == ELEMENT )   
             P(0,j) *= MathOperatorLHS<dim>::MTRL[0](0,0) * phi();
           else if ( MathOperatorLHS<dim>::MaterialOperandPlacement() == NODE )
             P(0,j) *= MathOperatorLHS<dim>::MTRL[i](0,0) * phi();
              
         PT   *= P;
         PT   *= (det * e.WeightAtIntegrationPoint(i));
        
         MathOperatorLHS<dim>::LHS += PT;
      }

    cout <<"\nNumIntegral_PT_lhsop_P_dV<dim>::N_to_P: LHS contribution element: "<< e.Idx() << endl;    
    MathOperatorLHS<dim>::LHS.Out(cout);
      
} // end ComputeContribution


template class NumIntegral_PT_lhsop_P_dV<1U,Element<1U> >;
template class NumIntegral_PT_lhsop_P_dV<2U,Element<2U> >;
template class NumIntegral_PT_lhsop_P_dV<3U,Element<3U> >;

template class NumIntegral_PT_lhsop_P_dV<1U,Face<1U> >;
template class NumIntegral_PT_lhsop_P_dV<2U,Face<2U> >;
template class NumIntegral_PT_lhsop_P_dV<3U,Face<3U> >;

} // csp
