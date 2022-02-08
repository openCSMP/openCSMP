#include "PointSource_lhsop.h"
#include "InterFace.h"
#include "ErrorHandler.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"


using namespace std;

namespace csmp {

template<size_t dim,class SIMPLEX>
PointSource_lhsop<dim,SIMPLEX>::PointSource_lhsop( const PropertyDatabase<dim>& pref, const char* oper, const char* basic, const char* test )
  : MathOperatorLHS<dim>(pref,oper,basic,test),
    SRC_(8)
 {
    MathOperatorLHS<dim>::Name("PointSource_lhsop", oper, test );
 
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
     if ( MathOperatorLHS<dim>::TestOperandType() != SCALAR ||
          MathOperatorLHS<dim>::TestOperandPlacement() != NODE )
       csmp_error.notice( ERROR, "MathOperatorLHS->PointSource_lhsop<dim>::(constructor)",
                          test, "must be a scalar variable placed on the node" ); 
                          
     if ( MathOperatorLHS<dim>::MaterialOperandType() != SCALAR ||
          MathOperatorLHS<dim>::MaterialOperandPlacement() != NODE )
       csmp_error.notice( ERROR, "MathOperatorLHS->PointSource_lhsop<dim>::(constructor)",
                          oper, "must be a scalar variable placed on the node" ); 
 }




/// reads the values of basic operand from element nodes
template<size_t dim,class SIMPLEX>
void PointSource_lhsop<dim,SIMPLEX>::GetOperands( SIMPLEX& e )
   { 
      e.NodePropertyVector( MathOperatorLHS<dim>::MaterialOperandKey(), SRC_ );
       
      // taking into account that the point source contributes to several elements
      for ( size_t i=0U; i<e.Nodes(); i++ )
        SRC_[i]() /= static_cast<double>(e.N(i)->Parents());
        
   } // end GetOperands





/**
 
Applies point (nodal) source terms = Neumann boundary conditions to 
any kind of FEM models. The contribution is divided by the number of 
parent elements of each node because it will be accumulated by each one
of of them. This is already done in the GetOperands() method. 

@param e The current element from which the nodal values are accumulated.
*/
template<size_t dim,class SIMPLEX>
void PointSource_lhsop<dim,SIMPLEX>::ComputeContribution( SIMPLEX& e )
{
   MathOperatorLHS<dim>::LHS.Resize( e.Nodes(), e.Nodes() );
   MathOperatorLHS<dim>::LHS.Zero();
   for ( size_t i=0U; i<e.Nodes(); i++ )
     MathOperatorLHS<dim>::LHS(i,i) = SRC_[i]();
     
} // end ComputeContribution


template class PointSource_lhsop<1U,Element<1U> >;
template class PointSource_lhsop<2U,Element<2U> >;
template class PointSource_lhsop<3U,Element<3U> >;

template class PointSource_lhsop<1U,Face<1U> >;
template class PointSource_lhsop<2U,Face<2U> >;
template class PointSource_lhsop<3U,Face<3U> >;

template class PointSource_lhsop<1U,InterFace<1U> >;
template class PointSource_lhsop<2U,InterFace<2U> >;
template class PointSource_lhsop<3U,InterFace<3U> >;

} // csmp
