#ifndef FRACTURECAPILLARYPRESSUREVISITOR_H
#define FRACTURECAPILLARYPRESSUREVISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"

namespace csmp{

  struct Index;
  template<uint32_t> class Model;


template<uint32_t dim>
class FractureCapillaryPressureVisitor : public Visitor<dim>
  {
  public:
    FractureCapillaryPressureVisitor( Model<dim>& model, const char* fracApTag,
                                      const char* fracPcTag, double pd, double ift, double ca);

    virtual ~FractureCapillaryPressureVisitor() {}

    virtual void Visit(Element<dim>* element);

  private:
    Model<dim>& model_;
    Index fracApKey_, fracPcKey_;
    const double Pd_matrix_, ift_, ca_;
    double fracPc_;

  };

} //csmp

#endif // FRACTURECAPILLARYPRESSUREVISITOR_H
