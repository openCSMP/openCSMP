#include "CopyReplaceVisitor_compute.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "Exception.h"
#if defined(_OPENMP )
#include "omp.h"
#endif

using namespace std;

namespace csmp {

  namespace copyReplaceVisitorCompileTimeDispatch {

    template<class V, size_t dim>
    void initializeVariable( V&, const PropertyDatabase<dim>&, const char* )
    {
         // does nothing
    }

    template<size_t dim>
    void initializeVariable( ArrayVariable& var, const PropertyDatabase<dim>& pdb, const char* propertyName )
    {
      var = ArrayVariable( propertyName, pdb );
    }

    template<size_t dim>
    void initializeVariable( FlaggedArrayVariable& var, const PropertyDatabase<dim>& pdb, const char* propertyName )
    {
      var = FlaggedArrayVariable( propertyName, pdb );
    }

  } // copyReplaceVisitorCompileTimeDispatch


/// replaces property b with the value of a
template<typename Var, size_t dim>
CopyReplaceVisitor<Var,dim>::CopyReplaceVisitor( const PropertyDatabase<dim>& p, 
                                                 const char* prop_a, const char* prop_b,
                                                 size_t nodes )
    : prop_key_a_(p.StorageKey(prop_a)), 
      prop_key_b_(p.StorageKey(prop_b)),
      nodes_visited_(nodes)
  { 
     if ( prop_key_a_.place != prop_key_b_.place )
       throw csmp::Exception( FATAL_ERROR, "CopyReplaceVisitor(constructor)",
                             "input and outprop property have different placement.");

     if ( prop_key_a_.type != prop_key_b_.type &&
          (prop_key_a_.type!=ARRAY && prop_key_b_.type!=FLAGGEDARRAY) &&
          (prop_key_a_.type!=FLAGGEDARRAY && prop_key_b_.type!=ARRAY))
       throw csmp::Exception( FATAL_ERROR, "CopyReplaceVisitor(constructor)",
                             "input and outprop property have different types.");

     if ( prop_key_a_.place == NODE and nodes == 0 )
       throw csmp::Exception( FATAL_ERROR, "CopyReplaceVisitor(constructor)",
                             "node count required for node variables" );

     if ( prop_key_a_.dataDepth != prop_key_b_.dataDepth )
       throw csmp::Exception( FATAL_ERROR, "CopyReplaceVisitor(constructor)",
                             "input and outprop property have different data depth (Arrays?).");

     this->ApplicationLevel(MODEL);
     this->ApplicationTarget(prop_key_a_.place);
     
     if ( this->ApplicationTarget() == NODE ) 
       Reset();

     copyReplaceVisitorCompileTimeDispatch::initializeVariable( variable_, p, prop_a );
#if defined(_OPENMP )
     thread_variable_.resize(omp_get_max_threads());
     for (size_t i = 0 ; i < thread_variable_.size() ; i++)
         copyReplaceVisitorCompileTimeDispatch::initializeVariable( thread_variable_[i], p, prop_a );
#endif
  }



template<typename Var, size_t dim>
CopyReplaceVisitor<Var,dim>::~CopyReplaceVisitor() 
 {  
 }

// CopyReplaceVisitor Methods ====================

template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Model<dim>* gptr )
 {
    if ( prop_key_a_.place == MODEL ) {
         assert( prop_key_a_.place == prop_key_b_.place );
         gptr->Read( prop_key_a_, variable_ );
         gptr->Store( prop_key_b_, variable_ );
      }
#if defined(_OPENMP )
    this->Visit(&gptr->Region("Model"));
#endif
 }

template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Region<dim>* gptr ) 
 {
    if ( prop_key_a_.place == REGION ) {
         assert( prop_key_a_.place == prop_key_b_.place );
         gptr->Read( prop_key_a_, variable_ );
         gptr->Store( prop_key_b_, variable_ );
      }
#if defined(_OPENMP )
    assert( prop_key_a_.place == prop_key_b_.place );
    if (prop_key_a_.place == NODE){
#pragma omp parallel
        {
            Node<dim>* np;
#pragma omp for
            for ( long int n= 0 ; n < gptr->Nodes(); n++ ){
                np = gptr->N(n);
                this->ComputeContribution(np);
            }
        }
    }
    if (prop_key_a_.place == ELEMENT){
#pragma omp parallel
        {
            Element<dim>* ep;
#pragma omp for
            for ( long int e= 0 ; e < gptr->Elements(); e++ ){
                ep = gptr->E(e);
                this->ComputeContribution(ep);
            }
        }
    }
#endif
 }

template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Element<dim>* eptr ) 
 { 
#if !defined(_OPENMP)
    if ( prop_key_a_.place == ELEMENT ) {
         eptr->Read( prop_key_a_, variable_ );
         eptr->Store( prop_key_b_, variable_ );
         return;
      }
      
    // IntegrationPoint properties
    if ( prop_key_a_.place == ELEMENT_INTEGRATION_POINT )
    for ( size_t i=0U; i<eptr->IntegrationPoints(); i++ ) {
         eptr->Read( i, prop_key_a_, variable_ );
         eptr->Store( i, prop_key_b_, variable_ );
      }
#endif
}

template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::ComputeContribution( Element<dim>* eptr )
 {
#if defined(_OPENMP )
    size_t tid = omp_get_thread_num();
    if ( prop_key_a_.place == ELEMENT ) {
         eptr->Read( prop_key_a_, thread_variable_[tid] );
         eptr->Store( prop_key_b_, thread_variable_[tid] );
    }

    // IntegrationPoint properties
    if ( prop_key_a_.place == ELEMENT_INTEGRATION_POINT )
        for ( size_t i=0U; i<eptr->IntegrationPoints(); i++ ) {
            eptr->Read( i, prop_key_a_, thread_variable_[tid] );
            eptr->Store( i, prop_key_b_, thread_variable_[tid] );
        }
#endif
 }


template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Node<dim>* nptr ) 
 { 
#if !defined(_OPENMP)
    if ( !nodes_visited_.GetBit( nptr->Idx() ) ) {
        nptr->Read( prop_key_a_, variable_ );
        nptr->Store( prop_key_b_, variable_ );
        nodes_visited_.SetBit( nptr->Idx(), true );
    }
#endif
}
template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::ComputeContribution( Node<dim>* nptr )
 {
#if defined(_OPENMP )
    size_t tid = omp_get_thread_num();
    nptr->Read( prop_key_a_, thread_variable_[tid] );
    nptr->Store( prop_key_b_, thread_variable_[tid] );
#endif
 }


template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Reset()
 { nodes_visited_.SetAll(false); }


template class CopyReplaceVisitor<ScalarVariable,1U>;
template class CopyReplaceVisitor<VectorVariable<1U>,1U>;
template class CopyReplaceVisitor<TensorVariable<1U>,1U>;

template class CopyReplaceVisitor<ScalarVariable,2U>;
template class CopyReplaceVisitor<VectorVariable<2U>,2U>;
template class CopyReplaceVisitor<TensorVariable<2U>,2U>;

template class CopyReplaceVisitor<ScalarVariable,3U>;
template class CopyReplaceVisitor<VectorVariable<3U>,3U>;
template class CopyReplaceVisitor<TensorVariable<3U>,3U>;

template class CopyReplaceVisitor<ArrayVariable,1>;
template class CopyReplaceVisitor<ArrayVariable,2>;
template class CopyReplaceVisitor<ArrayVariable,3>;

template class CopyReplaceVisitor<FlaggedArrayVariable,1>;
template class CopyReplaceVisitor<FlaggedArrayVariable,2>;
template class CopyReplaceVisitor<FlaggedArrayVariable,3>;


} // end namespace csmp
