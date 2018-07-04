#include "EffectiveStressVisitor.h"
#include "Model.h"
#include "Element.h"
#include "Exception.h"


using namespace std;

namespace csmp {

  template<size_t dim>
  EffectiveStressVisitor<dim>::EffectiveStressVisitor( Model<dim>& model, 
                                                       const char* stressTensor,
                                                       const char* effectiveStressTensor,
                                                       const char* fluidPressure )
    : Visitor<dim>( MODEL, ELEMENT ),
      propDB_( model.Database() ),
      sigmaKey_( propDB_.StorageKey(stressTensor) ),
      sigmaEffKey_( propDB_.StorageKey(effectiveStressTensor) ),
      fluidPressureKey_( propDB_.StorageKey(fluidPressure) ),
      meanStressKey_(), outputMeanStress_(false)
  {
    KeyChecks();
  }

  template<size_t dim>
  EffectiveStressVisitor<dim>::EffectiveStressVisitor( Model<dim>& model, 
                                                       const char* stressTensor,
                                                       const char* effectiveStressTensor,
                                                       const char* fluidPressure,
                                                       const char* meanStress )
    : Visitor<dim>( MODEL, ELEMENT ),
    propDB_( model.Database() ),
    sigmaKey_( propDB_.StorageKey(stressTensor) ),
    sigmaEffKey_( propDB_.StorageKey(effectiveStressTensor) ),
    fluidPressureKey_( propDB_.StorageKey(fluidPressure) ),
    meanStressKey_( propDB_.StorageKey(meanStress) ), outputMeanStress_(true)
  {
    KeyChecks();
  }


  template<size_t dim>
  void EffectiveStressVisitor<dim>::KeyChecks() const
  {
    if ( fluidPressureKey_.place != NODE || fluidPressureKey_.type != SCALAR )
      throw csmp::Exception( ERROR, "EffectiveStress_Visitor (constructor)",
      "'fluidPressure' must be a SCALAR variable on the node" );

    if ( sigmaKey_.place != ELEMENT_INTEGRATION_POINT || sigmaKey_.type != TENSOR )
      throw csmp::Exception( ERROR, "EffectiveStress_Visitor (constructor)",
      "'stressTensor' must be a tensor variable on the integration point" );

    if ( sigmaEffKey_.place != ELEMENT_INTEGRATION_POINT || sigmaEffKey_.type != TENSOR )
      throw csmp::Exception( ERROR, "EffectiveStress_Visitor (constructor)",
      "'effectiveStressTensor' must be a tensor variable on the integration point" );

    if(outputMeanStress_)
      if ( meanStressKey_.place != ELEMENT || meanStressKey_.type != SCALAR )
        throw csmp::Exception( ERROR, "EffectiveStress_Visitor (constructor)",
        "'meanStress' must be a SCALAR variable on the element" );

  }


  template<size_t dim>
  EffectiveStressVisitor<dim>::~EffectiveStressVisitor()
  {
  }


/**
    Interpolates fluid pressure to the integration points
    and subtracts it from the diagonal components of 
    the stress tensor.
    
    Optionally, the mean stress can be computed as well.
    
    @attention the input must be a diagonalized stress tensor
    else this will not work.
 
*/
  template<size_t dim>
  void EffectiveStressVisitor<dim>::Visit( Element<dim>* e )
  {
    for ( size_t i(0); i < e->FE()->IntegrationPoints(); ++i )
    {
      e->PropertyValueAtIntegrationPoint( fluidPressureKey_, i, fluidPressure_ );
      e->Read( i, sigmaKey_, sigma_ );

      for( size_t j(0); j < dim; ++j )
		  sigma_(j,j) -= fluidPressure_();

      e->Store( i, sigmaEffKey_, sigma_ );
    }
    
    if(outputMeanStress_)
    {
      meanStress_() = 0.;
      for ( size_t i(0); i < e->FE()->IntegrationPoints(); ++i )
      {
        e->Read( i, sigmaEffKey_, sigma_ );

        for( size_t i(0); i < dim; ++i )
          meanStress_() += sigma_(i,i);
      }
      meanStress_() /= dim*e->FE()->IntegrationPoints();
      e->Store( meanStressKey_, meanStress_ );
    }

  }


  template class EffectiveStressVisitor<1U>;
  template class EffectiveStressVisitor<2U>;
  template class EffectiveStressVisitor<3U>;


} // csmp

