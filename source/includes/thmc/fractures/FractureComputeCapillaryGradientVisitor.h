#ifndef FRACTURECOMPUTECAPILLARYGRADIENTVISITOR_H
#define FRACTURECOMPUTECAPILLARYGRADIENTVISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"
#include "TwoPhaseModel.h"

namespace csmp{

  struct Index;
  template<uint32_t> class Model;

  template<uint32_t dim>
  class FractureComputeCapillaryGradientVisitor : public Visitor<dim>
  {
  public:
      FractureComputeCapillaryGradientVisitor( Model<dim>& model, TwoPhaseModel<dim>& saturationFunctions,
                                               const char* fracPcGradientTag,
                                               const char* fracPcGradientTermTag,
                                               const char* fracPermTag );

      void Visit(Element<dim> *element ) override final;

  private:
      Model<dim>& model_;
      TwoPhaseModel<dim>* saturationfunctions_;
      Index fracPcGradientKey_, fracPcGradientTermKey_, fracPermKey_  ;
      VectorVariable<dim> fracPcGradient_;
  };

} //csmp

#endif // FRACTURECOMPUTECAPILLARYGRADIENTVISITOR_H
