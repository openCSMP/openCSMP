// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "FEM_DataInputVisitor.h"
#include "MeshRenumberVisitor.h"
#include "Node.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "PropertyDatabase.h"
#include "TensorVariable.h"
#include "FEM_Data.h"
#include "Exception.h"

using namespace std;

namespace csmp {

template<typename Var, uint32_t dim>
FEM_DataInputVisitor<Var,dim>::FEM_DataInputVisitor( Model<dim>& sg, 
                                                     const char* input_prop,
                                                     const FEM_Data<Var>& data  )
    : input_data_ref_(data),
      prop_key_(sg.Database().StorageKey(input_prop)),
      counter_(0U)
  { 
     this->ApplicationTarget(prop_key_.place);

     // renumber elements and nodes consecutively
     MeshRenumberVisitor<dim>  renumberer;
     sg.Accept( renumberer );

     if ( prop_key_.place == NODE ) {
         if ( renumberer.VisitedNodes() != data.Size() )
           throw csmp::Exception( ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between input data and number of nodes");
       }
     else if ( prop_key_.place == ELEMENT ) {
         if ( renumberer.VisitedElements() != data.Size() )
           throw csmp::Exception( ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between input data and number of elements");
       }
  }



template<typename Var, uint32_t dim>
FEM_DataInputVisitor<Var,dim>::~FEM_DataInputVisitor() 
 {  
 }

// FEM_DataInputVisitor Methods ====================

template<typename Var, uint32_t dim>
void FEM_DataInputVisitor<Var,dim>::Visit( Region<dim>* gptr ) 
 { 
    gptr->Store( prop_key_, input_data_ref_[ counter_++ ] );
 }


template<typename Var, uint32_t dim>
void FEM_DataInputVisitor<Var,dim>::Visit( Element<dim>* eptr ) 
 { 
    if ( prop_key_.place == ELEMENT ) {
         assert( eptr->Idx() < input_data_ref_.Size() ); 
         eptr->Store( prop_key_, input_data_ref_[ eptr->Idx() ] );
         return;
      }
    // IntegrationPoint properties
    for ( auto i{0U}; i<eptr->IntegrationPoints(); i++ ) {
         assert( counter_ < input_data_ref_.Size() ); 
         eptr->Store( i, prop_key_, input_data_ref_[ counter_++ ] );
      }
 }


template<typename Var, uint32_t dim>
void FEM_DataInputVisitor<Var,dim>::Visit( Node<dim>* nptr ) 
 { 
    assert( nptr->Idx() < input_data_ref_.Size() ); 
    nptr->Store( prop_key_, input_data_ref_[ nptr->Idx() ] );
 }


template<typename Var, uint32_t dim>
void FEM_DataInputVisitor<Var,dim>::Reset()
 { counter_=0U; }


template class FEM_DataInputVisitor<ScalarVariable,1U>;
template class FEM_DataInputVisitor<VectorVariable<1U>,1U>;
template class FEM_DataInputVisitor<TensorVariable<1U>,1U>;

template class FEM_DataInputVisitor<ScalarVariable,2U>;
template class FEM_DataInputVisitor<VectorVariable<2U>,2U>;
template class FEM_DataInputVisitor<TensorVariable<2U>,2U>;

template class FEM_DataInputVisitor<ScalarVariable,3U>;
template class FEM_DataInputVisitor<VectorVariable<3U>,3U>;
template class FEM_DataInputVisitor<TensorVariable<3U>,3U>;

} // end namespace csmp
