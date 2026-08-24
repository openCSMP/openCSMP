#include "DataInputVisitor.h"
#include "MeshRenumberVisitor.h"
#include "Node.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/// visitor for input of variables placed on the node, element and elsewhere
template<typename Var, uint32_t dim>
DataInputVisitor<Var,dim>::DataInputVisitor( Model<dim>& sg, 
                                             const char* input_prop,
                                             const vector<double>& input_data  )
    : input_data_ref_(input_data),
      prop_key_(sg.Database().StorageKey(input_prop)),
      counter_(0U)
  { 
     this->ApplicationLevel(MODEL);
     this->ApplicationTarget(prop_key_.place);

     // renumber elements and nodes consecutively
     MeshRenumberVisitor<dim>  renumberer;
     sg.Accept( renumberer );

     if ( prop_key_.place == NODE ) {
         if ( prop_key_.type == SCALAR and renumberer.VisitedNodes() != input_data.size() )
           throw csmp::Exception( FATAL_ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between scalar input data and number of nodes");
         
         else if ( prop_key_.type == VECTOR and renumberer.VisitedNodes() != input_data.size()*dim )
           throw csmp::Exception( FATAL_ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between vector input data and number of nodes");
         
         else if ( prop_key_.type == TENSOR and renumberer.VisitedNodes() != input_data.size()*dim*dim )
           throw csmp::Exception( FATAL_ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between tensor input data and number of nodes");
       }
     else if ( prop_key_.place == ELEMENT ) {
         if ( prop_key_.type == SCALAR and  renumberer.VisitedElements() != input_data.size() )
           throw csmp::Exception( FATAL_ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between scalar input data and number of elements");
         
         else if ( prop_key_.type == VECTOR and  renumberer.VisitedElements() != input_data.size()*dim )
           throw csmp::Exception( FATAL_ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between vector input data and number of elements");

         else if ( prop_key_.type == TENSOR and  renumberer.VisitedElements() != input_data.size()*dim*dim )
           throw csmp::Exception( FATAL_ERROR, "FEM_DataInputVisitor(constructor)",
                          "size mismatch between tensor input data and number of elements");
       }
  }



// DataInputVisitor Methods ====================

template<typename Var, uint32_t dim>
void DataInputVisitor<Var,dim>::Visit( Model<dim>* gptr )
 { 
    for ( uint32_t i{0U}; i<variable_.Size(); i++ )
      variable_.Component( i, input_data_ref_[ counter_++ ] );
      
    gptr->Store( prop_key_, variable_ );
 }


template<typename Var, uint32_t dim>
void DataInputVisitor<Var,dim>::Visit( Region<dim>* gptr ) 
 { 
    for ( uint32_t i{0U}; i<variable_.Size(); i++ )
      variable_.Component( i, input_data_ref_[ counter_++ ] );
      
    gptr->Store( prop_key_, variable_ );
 }


template<typename Var, uint32_t dim>
void DataInputVisitor<Var,dim>::Visit( Element<dim>* eptr ) 
 { 
    if ( prop_key_.place == ELEMENT ) {
         for ( uint32_t i{0U}; i<variable_.Size(); i++ )
           variable_.Component( i, input_data_ref_[ counter_++ ] );
         eptr->Store( prop_key_, variable_ );
         return;
      }
    // IntegrationPoint properties
    /// TODO: @todo SKM: other placements like FV integration points are not considered here yet
    for ( uint32_t i{0U}; i<eptr->IntegrationPoints(); i++ ) {
         for ( uint32_t j{0U}; j<variable_.Size(); j++ )
           variable_.Component( j, input_data_ref_[ counter_++ ] );
         eptr->Store( i, prop_key_, variable_ );
      }
 }


template<typename Var, uint32_t dim>
void DataInputVisitor<Var,dim>::Visit( Face<dim>* eptr )
  {
    if ( prop_key_.place == FACE ) {
      for ( uint32_t i{0U}; i<variable_.Size(); i++ )
        variable_.Component( i, input_data_ref_[ counter_++ ] );
      eptr->Store( prop_key_, variable_ );
      return;
    }
    // IntegrationPoint properties
    assert( prop_key_.place == FACE_INTEGRATION_POINT );
    /// TODO: @todo SKM: other placements like FV integration points are not considered here yet
    for ( uint32_t i{0U}; i<eptr->IntegrationPoints(); i++ ) {
      for ( uint32_t j{0U}; j<variable_.Size(); j++ )
        variable_.Component( j, input_data_ref_[ counter_++ ] );
      eptr->Store( i, prop_key_, variable_ );
    }
  }

  
template<typename Var, uint32_t dim>
void DataInputVisitor<Var,dim>::Visit( InterFace<dim>* eptr )
  {
    if ( prop_key_.place == INTER_FACE ) {
      for ( uint32_t i{0U}; i<variable_.Size(); i++ )
        variable_.Component( i, input_data_ref_[ counter_++ ] );
      eptr->Store( prop_key_, variable_ );
      return;
    }
    // IntegrationPoint properties
    assert( prop_key_.place == INTER_FACE_INTEGRATION_POINT );
    /// TODO: @todo SKM: other placements like FV integration points are not considered here yet
    for ( uint32_t i{0U}; i<eptr->IntegrationPoints(); i++ ) {
      for ( uint32_t j{0U}; j<variable_.Size(); j++ )
        variable_.Component( j, input_data_ref_[ counter_++ ] );
      eptr->Store( i, prop_key_, variable_ );
    }
  }

  

template<typename Var, uint32_t dim>
void DataInputVisitor<Var,dim>::Visit( Node<dim>* nptr ) 
 {
     const size_t offset(nptr->Idx() * variable_.Size());
 
     for ( uint32_t i{0U}; i<variable_.Size(); i++ )
       variable_.Component( i, input_data_ref_[ offset + i ] );

     nptr->Store( prop_key_, variable_ );
 }


template<typename Var, uint32_t dim>
void DataInputVisitor<Var,dim>::Reset()
 { counter_=0U; }


template class DataInputVisitor<ScalarVariable,1U>;
template class DataInputVisitor<VectorVariable<1U>,1U>;
template class DataInputVisitor<TensorVariable<1U>,1U>;

template class DataInputVisitor<ScalarVariable,2U>;
template class DataInputVisitor<VectorVariable<2U>,2U>;
template class DataInputVisitor<TensorVariable<2U>,2U>;

template class DataInputVisitor<ScalarVariable,3U>;
template class DataInputVisitor<VectorVariable<3U>,3U>;
template class DataInputVisitor<TensorVariable<3U>,3U>;

} // end namespace csmp
