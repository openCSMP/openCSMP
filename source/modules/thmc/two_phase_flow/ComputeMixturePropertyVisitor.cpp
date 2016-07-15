#include "ComputeMixturePropertyVisitor.h"
#include "Model.h"


namespace csmp{

  template<size_t dim>
  ComputeMixturePropertyVisitor<dim>::ComputeMixturePropertyVisitor( Model<dim>& model,
                                                                     const char* saturationOilTag ,
                                                                     const char* saturationWaterTag,
                                                                     const char* multiplierOilTag,
                                                                     const char* multiplierWaterTag,
                                                                     const char* divisorOilTag,
                                                                     const char* divisorWaterTag,
                                                                     const char* mixturePropertyTag,
                                                                     const char* externalMultiplierTag )
      : Visitor<dim>( MODEL, ELEMENT ), model_( model ), saturationOilKey_( model.Database().StorageKey( saturationOilTag ) ),
        saturationWaterKey_( model.Database().StorageKey( saturationWaterTag ) ),
        multiplierOilKey_( model.Database().StorageKey( multiplierOilTag ) ),
        multiplierWaterKey_( model.Database().StorageKey( multiplierWaterTag ) ),
        divisorOilKey_( model.Database().StorageKey( divisorOilTag ) ),
        divisorWaterKey_( model.Database().StorageKey( divisorWaterTag ) ),
        mixtureProductKey_( model.Database().StorageKey( mixturePropertyTag ) ),
        externalMultiplierKey_( model.Database().StorageKey( externalMultiplierTag ))
  {
  }

  template<size_t dim>
  void ComputeMixturePropertyVisitor<dim>::Visit( Element<dim>* element )
  {
      element->PropertyValueAtBaryCenter( saturationOilKey_, satOil_ );
      element->PropertyValueAtBaryCenter( saturationWaterKey_, satWater_ );
      element->PropertyValueAtBaryCenter( multiplierOilKey_, multOil_ );
      element->PropertyValueAtBaryCenter( multiplierWaterKey_, multWater_ );
      element->PropertyValueAtBaryCenter( divisorOilKey_, divOil_ );
      element->PropertyValueAtBaryCenter( divisorWaterKey_, divWater_ );
      element->Read( externalMultiplierKey_, externalMultiplier_ );

      mixProperty_ = satOil_.Value() * multOil_.Value() / divOil_.Value() + satWater_.Value() * multWater_.Value() / divWater_.Value();
      mixProperty_ *= externalMultiplier_;

      element->Store(  mixtureProductKey_, mixProperty_ );
  }

  template class ComputeMixturePropertyVisitor<1>;
  template class ComputeMixturePropertyVisitor<2>;
  template class ComputeMixturePropertyVisitor<3>;

} //csmp
