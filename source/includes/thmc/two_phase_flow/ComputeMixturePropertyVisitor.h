#ifndef COMPUTEMIXTUREPROPERTYVISITOR_H
#define COMPUTEMIXTUREPROPERTYVISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"

namespace csmp{

template<uint32_t dim>
class ComputeMixturePropertyVisitor : public Visitor<dim>
  {
    public:
      ComputeMixturePropertyVisitor( Model<dim> &model,
                                     const char* saturationOilTag,
                                     const char* saturationWaterTag,
                                     const char* multiplierOilTag,
                                     const char* multiplierWaterTag,
                                     const char* divisorOilTag,
                                     const char* divisorWaterTag,
                                     const char* mixturePropertyTag,
                                     const char* externalMultiplierTag );

      virtual ~ComputeMixturePropertyVisitor() {}

      virtual void Visit(Element<dim>* element);

    private:
      Model<dim>& model_;
      Index saturationOilKey_, saturationWaterKey_, multiplierOilKey_, multiplierWaterKey_, divisorOilKey_, divisorWaterKey_, mixtureProductKey_, externalMultiplierKey_;
      ScalarVariable satOil_, satWater_, multOil_, multWater_, divOil_, divWater_, externalMultiplier_, mixProperty_;
  };
} //csmp

#endif // COMPUTEMIXTUREPROPERTYVISITOR_H




