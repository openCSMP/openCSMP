#ifndef CAPILLARYNUMBERVISITOR_H
#define CAPILLARYNUMBERVISITOR_H


#include "Visitor.h"
#include "ScalarVariable.h"

namespace csmp
{

  struct Index;
  template<uint32_t> class Model;


  template<uint32_t dim>
  class CapillaryNumberVisitor : public Visitor<dim>
    {
    public:
      CapillaryNumberVisitor( Model<dim>& model, const char* velTag,
                              const char* viscTag, const char* nCapTag, double ift );

      virtual ~CapillaryNumberVisitor() {}

      virtual void Visit(Element<dim>* element);

      private:
        Model<dim>& model_;
        Index velKey_, nCapKey_, viscKey_;
        const double ift_;
        ScalarVariable capillaryNumber_;
        VectorVariable<dim> velocity_;
      };

} //csmp


#endif // CAPILLARYNUMBERVISITOR_H
