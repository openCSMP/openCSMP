#ifndef FRACTUREBROOKSCOREYPARAMETERVISITOR_H
#define FRACTUREBROOKSCOREYPARAMETERVISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"


namespace csmp{

  struct Index;
  template<uint32_t> class Model;


template<uint32_t dim>
class FractureBrooksCoreyParameterVisitor : public Visitor<dim>
  {
  public:
    FractureBrooksCoreyParameterVisitor(Model<dim>& model, const char* frapTag,
                                        const char* brooksCoreyLambdaTag, double meanPoreRadius );
    virtual ~FractureBrooksCoreyParameterVisitor() {}

    virtual void Visit(Element<dim>* element);

  private:
    Model<dim>& model_;
    Index frapKey_, brooksCoreyLambdaKey_;
    const double meanPoreRadius_;
  };

} //csmp


#endif // FRACTUREBROOKSCOREYPARAMETERVISITOR_H
