#include "FEM_DataOutputVisitor.h"
#include "Node.h"
#include "Element.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "TensorVariable.h"
#include "FEM_Data.h"

using namespace std;

namespace csmp {

template<typename Var, uint32_t dim>
FEM_DataOutputVisitor<Var,dim>::FEM_DataOutputVisitor( const PropertyDatabase<dim>& p, 
                                                       const char* input_prop,
                                                       FEM_Data<Var>& data  )
    : prop_key_(p.StorageKey(input_prop)),
      input_data_ref_(data),
      counter_(0U)
  { 
     this->ApplicationTarget(prop_key_.place);
  }



template<typename Var, uint32_t dim>
FEM_DataOutputVisitor<Var,dim>::~FEM_DataOutputVisitor() 
 {  
 }

// FEM_DataOutputVisitor Methods ====================

template<typename Var, uint32_t dim>
void FEM_DataOutputVisitor<Var,dim>::Visit( Region<dim>* gptr ) 
 { 
    gptr->Read( prop_key_, input_data_ref_[ counter_++ ] );
 }


template<typename Var, uint32_t dim>
void FEM_DataOutputVisitor<Var,dim>::Visit( Element<dim>* eptr ) 
 { 
    if ( prop_key_.place == ELEMENT ) {
         assert( eptr->Idx() < input_data_ref_.Size() ); 
         eptr->Read( prop_key_, input_data_ref_[ eptr->Idx() ] );
         return;
      }
    // IntegrationPoint properties
    for ( auto i{0}; i<eptr->IntegrationPoints(); i++ ) {
         assert( counter_ < input_data_ref_.Size() ); 
         eptr->Read( i, prop_key_, input_data_ref_[ counter_++ ] );
      }
 }


template<typename Var, uint32_t dim>
void FEM_DataOutputVisitor<Var,dim>::Visit( Node<dim>* nptr ) 
 { 
    assert( nptr->Idx() < input_data_ref_.Size() ); 
    nptr->Read( prop_key_, variable_ );
 }


template<typename Var, uint32_t dim>
void FEM_DataOutputVisitor<Var,dim>::Reset()
 { counter_=0U; }


template class FEM_DataOutputVisitor<ScalarVariable,1U>;
template class FEM_DataOutputVisitor<VectorVariable<1U>,1U>;
template class FEM_DataOutputVisitor<TensorVariable<1U>,1U>;

template class FEM_DataOutputVisitor<ScalarVariable,2U>;
template class FEM_DataOutputVisitor<VectorVariable<2U>,2U>;
template class FEM_DataOutputVisitor<TensorVariable<2U>,2U>;

template class FEM_DataOutputVisitor<ScalarVariable,3U>;
template class FEM_DataOutputVisitor<VectorVariable<3U>,3U>;
template class FEM_DataOutputVisitor<TensorVariable<3U>,3U>;

} // end namespace csmp
