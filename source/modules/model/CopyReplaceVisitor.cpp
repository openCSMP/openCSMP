#include "CopyReplaceVisitor.h"
#include "Node.h"
#include "Element.h"
#include "Model.h"
#include "Region.h"
#include "PropertyDatabase.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "Exception.h"

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




/// custom constructor: replaces property b with the value of a
template<typename Var, size_t dim>
CopyReplaceVisitor<Var,dim>::CopyReplaceVisitor( const PropertyDatabase<dim>& p, 
                                                 const char* prop_a, const char* prop_b )
    : Visitor<dim>(MODEL,p.Placement(prop_a)),
      prop_key_a_(p.StorageKey(prop_a)),
      prop_key_b_(p.StorageKey(prop_b))
  {
     if ( prop_key_a_.place != prop_key_b_.place )
       throw csmp::Exception( CSMP_FATAL_ERROR, "CopyReplaceVisitor(constructor)",
                             "input and outprop property have different placement.");

     if ( prop_key_a_.type != prop_key_b_.type &&
          (prop_key_a_.type!=ARRAY && prop_key_b_.type!=FLAGGEDARRAY) &&
          (prop_key_a_.type!=FLAGGEDARRAY && prop_key_b_.type!=ARRAY))
       throw csmp::Exception( CSMP_FATAL_ERROR, "CopyReplaceVisitor(constructor)",
                             "input and outprop property have different types.");

     if ( prop_key_a_.dataDepth != prop_key_b_.dataDepth )
       throw csmp::Exception( CSMP_FATAL_ERROR, "CopyReplaceVisitor(constructor)",
                             "input and outprop property have different data depth (Arrays).");

     copyReplaceVisitorCompileTimeDispatch::initializeVariable( variable_, p, prop_a );
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
         gptr->Read( prop_key_a_, variable_ );
         gptr->Store( prop_key_b_, variable_ );
      }
 }


template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Region<dim>* gptr ) 
 {
    if ( prop_key_a_.place == REGION ) {
         gptr->Read( prop_key_a_, variable_ );
         gptr->Store( prop_key_b_, variable_ );
      }
 }


/**
   Considers wide spectrum of placements:
   - element
   - element integration point
   - sector integration point
   - finite volume facet integration point
*/
template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Element<dim>* eptr ) 
 { 
    if ( prop_key_a_.place == ELEMENT ) {
         eptr->Read( prop_key_a_, variable_ );
         eptr->Store( prop_key_b_, variable_ );
         return;
      }
      
    // IntegrationPoint properties
    if ( prop_key_a_.place == ELEMENT_INTEGRATION_POINT ) {
        for ( size_t i=0U; i<eptr->IntegrationPoints(); i++ ) {
             eptr->Read( i, prop_key_a_, variable_ );
             eptr->Store( i, prop_key_b_, variable_ );
          }
        return;
      }

    // IntegrationPoint properties
    if ( prop_key_a_.place == SECTOR_INTEGRATION_POINT ) {
        assert( eptr->FV_Stencil() != nullptr );
        for ( size_t i=0U; i<eptr->Sectors(); i++ )
          for ( size_t j=0U; j<eptr->IntegrationPointsPerSector(); ++j ) {
               eptr->Read( i, j, prop_key_a_, variable_ );
               eptr->Store( i, j, prop_key_b_, variable_ );
            }
        return;
      }

    // IntegrationPoint properties
    if ( prop_key_a_.place == FACET_INTEGRATION_POINT ) {
        assert( eptr->FV_Stencil() != nullptr );
        for ( size_t i=0U; i<eptr->Facets(); i++ )
          for ( size_t j=0U; j<eptr->IntegrationPointsPerFacet(); ++j ) {
               eptr->Read( i, j, prop_key_a_, variable_ );
               eptr->Store( i, j, prop_key_b_, variable_ );
            }
        return;
      }
   
 } // end Visit(Element)




template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Face<dim>* eptr )
 { 
    if ( prop_key_a_.place == FACE ) {
         eptr->Read( prop_key_a_, variable_ );
         eptr->Store( prop_key_b_, variable_ );
         return;
      }
      
    // IntegrationPoint properties
    if ( prop_key_a_.place == FACE_INTEGRATION_POINT ) {
        for ( size_t i=0U; i<eptr->IntegrationPoints(); i++ ) {
             eptr->Read( i, prop_key_a_, variable_ );
             eptr->Store( i, prop_key_b_, variable_ );
          }
        return;
      }

    // IntegrationPoint properties
    if ( prop_key_a_.place == FACE_SECTOR_INTEGRATION_POINT ) {
        assert( eptr->FV_Stencil() != nullptr );
        for ( size_t i=0U; i<eptr->Sectors(); i++ )
          for ( size_t j=0U; j<eptr->IntegrationPointsPerSector(); ++j ) {
               eptr->Read( i, j, prop_key_a_, variable_ );
               eptr->Store( i, j, prop_key_b_, variable_ );
            }
        return;
      }

    // IntegrationPoint properties
    if ( prop_key_a_.place == FACE_FACET_INTEGRATION_POINT ) {
        assert( eptr->FV_Stencil() != nullptr );
        for ( size_t i=0U; i<eptr->Facets(); i++ )
          for ( size_t j=0U; j<eptr->IntegrationPointsPerFacet(); ++j ) {
               eptr->Read( i, j, prop_key_a_, variable_ );
               eptr->Store( i, j, prop_key_b_, variable_ );
            }
        return;
      }
   
 } // end Visit(Face)




template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( InterFace<dim>* eptr )
 { 
    if ( prop_key_a_.place == INTER_FACE ) {
         eptr->Read( prop_key_a_, variable_ );
         eptr->Store( prop_key_b_, variable_ );
         return;
      }
      
    // IntegrationPoint properties
    if ( prop_key_a_.place == INTER_FACE_INTEGRATION_POINT ) {
        for ( size_t i=0U; i<eptr->IntegrationPoints(); i++ ) {
             eptr->Read( i, prop_key_a_, variable_ );
             eptr->Store( i, prop_key_b_, variable_ );
          }
        return;
      }

    // IntegrationPoint properties
    if ( prop_key_a_.place == INTER_FACE_SECTOR_INTEGRATION_POINT ) {
        assert( eptr->FV_Stencil() != nullptr );
        for ( size_t i=0U; i<eptr->Sectors(); i++ )
          for ( size_t j=0U; j<eptr->IntegrationPointsPerSector(); ++j ) {
               eptr->Read( i, j, prop_key_a_, variable_ );
               eptr->Store( i, j, prop_key_b_, variable_ );
            }
        return;
      }

    // IntegrationPoint properties
    if ( prop_key_a_.place == INTER_FACE_FACET_INTEGRATION_POINT ) {
        assert( eptr->FV_Stencil() != nullptr );
        for ( size_t i=0U; i<eptr->Facets(); i++ )
          for ( size_t j=0U; j<eptr->IntegrationPointsPerFacet(); ++j ) {
               eptr->Read( i, j, prop_key_a_, variable_ );
               eptr->Store( i, j, prop_key_b_, variable_ );
            }
        return;
      }
   
 } // end Visit(InterFace)




/**
     @attention watch out that you do not visit the same node
     many times over.
*/
template<typename Var, size_t dim>
void CopyReplaceVisitor<Var,dim>::Visit( Node<dim>* nptr ) 
 { 
    nptr->Read( prop_key_a_, variable_ );
    nptr->Store( prop_key_b_, variable_ );
 }


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
