#include "CVFE_NumIntegral_dNT_op_dN_dV.h"
#include "TwoPhaseModel.h"
#include "Element.h"
#include "Face.h"

using namespace std;

namespace csmp {

/**
 
The Operand which is used here can be both, an element or a nodal variable
which is then interpolated to the integration points to obtain the 
integral properties.  
*/
template<uint32_t dim,class CELL>
CVFE_NumIntegral_dNT_op_dN_dV<dim,CELL>::CVFE_NumIntegral_dNT_op_dN_dV( const PropertyDatabase<dim>& pref,
                                                                   TwoPhaseModel<dim>&   kri,
                                                                   const char*           oper, 
                                                                   const char*           basic, 
                                                                   const char*           test ) 
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    kri_(kri),
    DN_(dim,3U), FN_(0,dim), LK_(dim,dim)
{
    MathOperatorLHS<dim>::Name("CVFE_NumIntegral_dNT_op_dN_dV", oper, basic, test );
    
    if ( MathOperatorLHS<dim>::MaterialOperandPlacement() != ELEMENT_INTEGRATION_POINT )
      throw csmp::Exception( ERROR, "CVFE_NumIntegral_dNT_op_dN_dV<dim>::(constructor)", 
                      basic, "Operand must be a placed on the integration points." );

    if ( MathOperatorLHS<dim>::BasicOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::BasicOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "CVFE_NumIntegral_dNT_op_dN_dV<dim>::(constructor)", 
                      basic, "Operand (basic) must be a scalar property placed on the nodes." );

    if ( MathOperatorLHS<dim>::TestOperandPlacement() != NODE || 
         MathOperatorLHS<dim>::TestOperandType() != SCALAR )
      throw csmp::Exception( ERROR, "CVFE_NumIntegral_dNT_op_dN_dV<dim>::(constructor)", 
                      test, "Operand (test) must be a scalar property placed on the nodes." );
}





/**
 
Compute "total mobility - permeability products at FV integration points, 
storing the results in the material matrices of the base class.  
*/
template<uint32_t dim,class CELL>
void CVFE_NumIntegral_dNT_op_dN_dV<dim,CELL>::GetOperands( CELL& e )
 {
    // initialising the material constants of the relperm model
    kri_.Initialize(e);

    this->MTRL.resize(e.FV()->Facets());

    // the scalar permeability is stored in the relperm model
    if ( this->MaterialOperandType() == SCALAR ) {
         this->MTRL[0U].AssignToDiagonalAndZeroOffDiagonal( dim, kri_.Permeability() );
      }
    else if ( this->MaterialOperandType() == VECTOR ) {
         // the vector permeability needs to be read again
         VectorVariable<dim> vc;
         e.Read( this->MaterialOperandKey(), vc );
         this->MTRL[0U].AssignToDiagonal( vc );
         for ( size_t i=1U; i<e.FV()->Facets(); i++ )
           this->MTRL[i] =  this->MTRL[0U];
      }
    else if ( this->MaterialOperandType() == TENSOR ) {
         // so does the permeability tensor
         TensorVariable<dim> ts;
         e.Read( this->MaterialOperandKey(), ts );
         this->MTRL[0U] = ts; 
         for ( size_t i=1U; i<e.FV()->Facets(); i++ )
           this->MTRL[i] =  this->MTRL[0U];
      }

    for ( auto i{0}; i<e.FV()->Facets(); i++ )
      {
         kri_.InitializeForFacetIntegrationPoint( i, 0U, e );
         this->MTRL[i] *= kri_.TotalMobility();
      }

 } // end 





/** Laplacian operator of shape function derivatives squared.
*/
template<uint32_t dim,class CELL>
void CVFE_NumIntegral_dNT_op_dN_dV<dim,CELL>::ComputeContribution( CELL& e )
 {
    // initialize output matrix
    MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
    MathOperatorLHS<dim>::LHS.Zero();

    const double fv_integration_weight(1./e.FV()->Facets());

     // for all finite-volume facets
    for ( auto i{0}; i<e.FV()->Facets(); i++ )
      {
         // identifying the finite volumes to which the flux will be distributed
         uint32_t inside_node  = e.FV()->InsideNode(i),
                  outside_node = e.FV()->OutsideNode(i);

         Point<dim> rst  = e.FV()->FacetIntegrationPoint( i, 0U );
         double   detJ = (e).dN_At( rst, DN_ );

         LK_ = MathOperatorLHS<dim>::MTRL[i] * DN_;

         Point<dim>  fn = (e).FacetNormal(i);
         fn *= detJ * fv_integration_weight * (e).FacetArea(i);
         FN_.AssignRow( 0U, fn );
         FN_ *= LK_;

         if ( FN_(0U,0U) < 0. ) {
              // 1. fluxes coming into the sector (fluxes = negative since normals are pointing outward) are added
              //    (outside node = upstream)
              // -------------------------------------------------------------------------------------------------
              MathOperatorLHS<dim>::LHS(inside_node,outside_node) = FN_(0U,0U);  // incoming flux
              
              // 2. outgoing fluxes are added to the matrix diagonal
              // ---------------------------------------------------
              MathOperatorLHS<dim>::LHS(outside_node,inside_node) = -FN_(0U,0U); // outgoing flux
           } 
         else {
              MathOperatorLHS<dim>::LHS(inside_node,outside_node) = -FN_(0U,0U);  // outgoing flux
              MathOperatorLHS<dim>::LHS(outside_node,inside_node) =  FN_(0U,0U);  // incoming flux
           }
      }

} // end ComputeContribution


//cout <<"\nCVFE_NumIntegral_dNT_op_dN_dV: on Element "<< e.Idx() << endl;
//MathOperatorLHS<dim>::LHS.Out();

// for ( auto i{0}; i<this->LHS.Rows(); i++ )
//   if ( this->LHS(i,i) < numeric_limits::epsilon() ) 
//     cout <<"\nCVFE_NumIntegral_dNT_op_dN_dV: zero element in diagonal of element matrix."; 


template class CVFE_NumIntegral_dNT_op_dN_dV<1U,Element<1U> >;
template class CVFE_NumIntegral_dNT_op_dN_dV<2U,Element<2U> >;
template class CVFE_NumIntegral_dNT_op_dN_dV<3U,Element<3U> >;

//template class CVFE_NumIntegral_dNT_op_dN_dV<1U,Face<1U> >;
//template class CVFE_NumIntegral_dNT_op_dN_dV<2U,Face<2U> >;
//template class CVFE_NumIntegral_dNT_op_dN_dV<3U,Face<3U> >;

} // csmp











