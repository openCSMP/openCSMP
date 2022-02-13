#include "CVFEM_PointSource_rhsop.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

using namespace std;

namespace csmp {

template<uint32_t dim, class CELL>
CVFEM_PointSource_rhsop<dim,CELL>::~CVFEM_PointSource_rhsop() {}

template<uint32_t dim, class CELL>
CVFEM_PointSource_rhsop<dim,CELL>::CVFEM_PointSource_rhsop( const PropertyDatabase<dim>& pref, const char* oper, const char* test )
  : CVFEM_MathOperatorRHS<dim>(pref,oper,test),
    SRC_(3)
 {
    MathOperatorRHS<dim>::Name("CVFEM_PointSource_rhsop", oper, test );
 
     if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR || 
          MathOperatorRHS<dim>::MaterialOperandPlacement() != NODE ) 
     {
        cout <<"\nMathOperatorRHS->CVFEM_PointSource_rhsop<dim,CELL>::(constructor): ";
        cout <<"Fatal Error: Operand must be a scalar variable placed on the node. Terminating..."<< endl;
        throw invalid_argument("CVFEM_PointSource_rhsop<csp_float,dim>::CVFEM_PointSource_rhsop");
     }
 }


template<uint32_t dim, class CELL>
void CVFEM_PointSource_rhsop<dim,CELL>::GetOperands( CELL& e )
   { 
      e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), SRC_ );
      // taking into account that the source contributes to several elements
      for ( auto i{0}; i<e.Nodes(); i++ )
        SRC_[i] /= static_cast<double>(e.N(i)->Parents()); 

   } // end GetOperands


template<uint32_t dim, class CELL>
void CVFEM_PointSource_rhsop<dim,CELL>::GetOperandsCVFEM( CELL& e, csmp::Index upwind_var_key  )
   { 

      e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), SRC_ );
      e.NodePropertyVector( upwind_var_key, upwind_var_);
       
      // taking into account that the source contributes to several elements
      for ( auto i{0}; i<e.Nodes(); i++ )
        SRC_[i] /= static_cast<double>(e.N(i)->Parents()); 

      for (auto i = 0; i < upwind_var_.size(); i++)
         SRC_[i]() *= upwind_var_[i]();


   } // end GetOperands


template<uint32_t dim, class CELL>
void CVFEM_PointSource_rhsop<dim,CELL>::ComputeContribution( CELL& e )
{
   MathOperatorRHS<dim>::RHS.resize( e.Nodes() );
   for ( auto i=0; i<e.Nodes(); i++ )
     MathOperatorRHS<dim>::RHS[i] = SRC_[i]();
     
} // end ComputeContribution

template class CVFEM_PointSource_rhsop<1U,Element<1U> >;
template class CVFEM_PointSource_rhsop<2U,Element<2U> >;
template class CVFEM_PointSource_rhsop<3U,Element<3U> >;

} // csmp
