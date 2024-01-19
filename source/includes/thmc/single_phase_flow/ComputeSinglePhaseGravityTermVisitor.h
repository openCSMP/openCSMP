#ifndef COMPUTE_SINGLE_PHASE_GRAVITY_TERM_VISITOR_H
#define COMPUTE_SINGLE_PHASE_GRAVITY_TERM_VISITOR_H

#include "Visitor.h"
#include "FiniteElementManager.h"

namespace csmp{

struct Index;
class ScalarVariable;
template<uint32_t> class VectorVariable;
template<uint32_t> class Model;
template<uint32_t> class TwoPhaseModel;

template<uint32_t dim, template<uint32_t> class CELL=Element>
class ComputeSinglePhaseGravityTermVisitor : public Visitor<dim>
{
  public:
    ComputeSinglePhaseGravityTermVisitor(Model<dim>& model,
                                         const char* permeabilityTag,
                                         const char* viscosityTag,
                                         const char* densityTag ,
                                         const char* gravityVectorTag,
                                         const char* model_gravity_vector,
                                         const char *densityTag2 = NULL); // post multiplies gravity term if provided

    virtual ~ComputeSinglePhaseGravityTermVisitor() {}

    virtual void Visit(CELL<dim>* );

  private:
    VectorVariable<dim>  gravity_unit_vector_;
    Index                permeabilityKey_, viscosityKey_, densityKey_,densityKey2_, gravityVectorKey_;
    const double         gravitational_acceleration_;
};

} //csmp

#endif // COMPUTE_SINGLE_PHASE_GRAVITY_TERM_VISITOR_H
