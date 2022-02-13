#ifndef EFFECTIVE_STRESS_VISITOR_H
#define EFFECTIVE_STRESS_VISITOR_H

#include "Visitor.h"
#include "TensorVariable.h"
#include "ScalarVariable.h"

namespace csmp {

  template<uint32_t> class PropertyDatabase;
  class ScalarVariable;
  template<uint32_t> class Model;

  /// effective- and mean stress, by post-processing of a diagonalized stress tensor with fluid pressure as input.
  template<uint32_t dim>
  class EffectiveStressVisitor : public Visitor<dim> {
  public:
    EffectiveStressVisitor( Model<dim>& model, 
                            const char* stressTensor,
                            const char* effectiveStressTensor,
                            const char* fluidPressure );
    EffectiveStressVisitor( Model<dim>& model, 
                            const char* stressTensor,
                            const char* effectiveStressTensor,
                            const char* fluidPressure,
                            const char* meanStress );

    virtual ~EffectiveStressVisitor();

    virtual void Visit(Element<dim>* );

  private:
    EffectiveStressVisitor();
    EffectiveStressVisitor( const EffectiveStressVisitor<dim>& );

    void KeyChecks() const;

  private:
    PropertyDatabase<dim>&   propDB_;
    TensorVariable<dim>      sigma_;     ///< Stress tensor (in diagonal form)
    ScalarVariable      fluidPressure_;  ///< pressure
    ScalarVariable      meanStress_; 
    Index               fluidPressureKey_;
    Index               sigmaKey_, sigmaEffKey_, meanStressKey_;
    bool                outputMeanStress_;
  };



} // csmp

#endif
