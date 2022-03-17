#include "PointSource_rhsop.h"
#include "InterFace.h"
#include "ErrorHandler.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"


using namespace std;

namespace csmp {

template<uint32_t dim,class CELL>
PointSource_rhsop<dim,CELL>::PointSource_rhsop( const PropertyDatabase<dim>& pref, const char* oper, const char* test )
  : MathOperatorRHS<dim>(pref,oper,test),
    SRC_(8)
 {
    MathOperatorRHS<dim>::Name("PointSource_rhsop", oper, test );
 
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
     if ( MathOperatorRHS<dim>::TestOperandType() != SCALAR || 
          MathOperatorRHS<dim>::TestOperandPlacement() != NODE )
       csmp_error.notice( ERROR, "MathOperatorRHS->PointSource_rhsop<dim>::(constructor)", 
                          test, "must be a scalar variable placed on the node" ); 
                          
     if ( MathOperatorRHS<dim>::MaterialOperandType() != SCALAR || 
          MathOperatorRHS<dim>::MaterialOperandPlacement() != NODE )
       csmp_error.notice( ERROR, "MathOperatorRHS->PointSource_rhsop<dim>::(constructor)", 
                          oper, "must be a scalar variable placed on the node" ); 
 }




/// reads the values of basic operand from element nodes
template<uint32_t dim,class CELL>
void PointSource_rhsop<dim,CELL>::GetOperands( const CELL& e )
   { 
      e.NodePropertyVector( MathOperatorRHS<dim>::MaterialOperandKey(), SRC_ );
       
      // taking into account that the point source contributes to several elements
      for ( auto i{0}; i<e.Nodes(); i++ )
        SRC_[i]() /= static_cast<double>(e.N(i)->Parents());
        
   } // end GetOperands





/**
 
Applies point (nodal) source terms = Neumann boundary conditions to 
any kind of FEM models. The contribution is divided by the number of 
parent elements of each node because it will be accumulated by each one
of of them. This is already done in the GetOperands() method. 

@param e The current element from which the value is accumulated.
*/
template<uint32_t dim,class CELL>
void PointSource_rhsop<dim,CELL>::ComputeContribution( const CELL& e )
{
   MathOperatorRHS<dim>::RHS.resize( e.Nodes() );
   for ( auto i{0}; i<e.Nodes(); i++ )
     MathOperatorRHS<dim>::RHS[i] = SRC_[i]();
     
} // end ComputeContribution


template class PointSource_rhsop<1U,Element<1U> >;
template class PointSource_rhsop<2U,Element<2U> >;
template class PointSource_rhsop<3U,Element<3U> >;

template class PointSource_rhsop<1U,Face<1U> >;
template class PointSource_rhsop<2U,Face<2U> >;
template class PointSource_rhsop<3U,Face<3U> >;

template class PointSource_rhsop<1U,InterFace<1U> >;
template class PointSource_rhsop<2U,InterFace<2U> >;
template class PointSource_rhsop<3U,InterFace<3U> >;

} // csmp
