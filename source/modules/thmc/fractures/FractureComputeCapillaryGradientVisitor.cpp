#include "FractureComputeCapillaryGradientVisitor.h"
#include "Model.h"
#include "ScalarVariable.h"

namespace csmp
{


  template<uint32_t dim>
  FractureComputeCapillaryGradientVisitor<dim>::FractureComputeCapillaryGradientVisitor( Model<dim> &model,
                                                                                         TwoPhaseModel<dim>& saturationFunctions,
                                                                                         const char* fracPcGradientTag,
                                                                                         const char* fracPcGradientTermTag,
                                                                                         const char* fracPermTag )

      : Visitor<dim> ( MODEL, ELEMENT ), model_( model ), saturationfunctions_( &saturationFunctions ),
        fracPcGradientKey_( model.Database().StorageKey( fracPcGradientTag ) ),
        fracPcGradientTermKey_( model.Database().StorageKey( fracPcGradientTermTag ) ),
        fracPermKey_( model.Database().StorageKey( fracPermTag ) )
  {
  }

  template<uint32_t dim>
  void FractureComputeCapillaryGradientVisitor<dim>::Visit( Element<dim>* element )
  {
      // setting up the relative permeability model
      saturationfunctions_->Initialize( *element );
      saturationfunctions_->InitializeForBaryCenter( *element );
      saturationfunctions_->EffectiveSaturation();

      element->Read( fracPcGradientKey_, fracPcGradient_ );

      for( uint32_t xyz = 0u; xyz < dim; ++xyz )
        fracPcGradient_(xyz) *= element->Read(fracPermKey_) * saturationfunctions_->MobilityPhase(1);

      element->Store( fracPcGradientTermKey_, fracPcGradient_ );

      //std::cout << "FractureComputeCapillaryGradientVisitor::Visit " << fracPcGradient_ << std::endl;

  }


template class FractureComputeCapillaryGradientVisitor<1>;
template class FractureComputeCapillaryGradientVisitor<2>;
template class FractureComputeCapillaryGradientVisitor<3>;

} //csmp
