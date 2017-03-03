#ifndef COMPUTESINGLEPHASEGRAVITYTERMVISITOR_H
#define COMPUTESINGLEPHASEGRAVITYTERMVISITOR_H

#include "Visitor.h"
#include "FiniteElementManager.h"

namespace csmp{

struct Index;
class ScalarVariable;
template<size_t> class VectorVariable;
template<size_t> class Model;
template<size_t> class TwoPhaseModel;

template<size_t dim>
class ComputeSinglePhaseGravityTermVisitor : public Visitor<dim>
{
  public:
    ComputeSinglePhaseGravityTermVisitor(Model<dim>& model,
                                         const char* permeabilityTag,
                                         const char* viscosityTag,
                                         const char* densityTag ,
                                         const char* gravityVectorTag,
                                         const char* model_gravity_vector,
                                         const char *densityTag2 = NULL); // this is meant to post multiply the gravity term if provided.

    virtual ~ComputeSinglePhaseGravityTermVisitor() {}

    virtual void Visit(Element<dim>* element);
    void ComputeContribution(Element<dim>* element);
    virtual void Visit(Model<dim>* model);
    virtual void Visit(Region<dim>* region);

  private:
    Model<dim>& model_;
    VectorVariable<dim>           gravity_unit_vector_;
    Index permeabilityKey_, viscosityKey_, densityKey_,densityKey2_, gravityVectorKey_;
    const double gravitational_acceleration_;

#if defined(_OPENMP )
    std::vector<FiniteElementManager> femgrs_; // one manager per thread
#endif

};

} //csmp

#endif // COMPUTESINGLEPHASEGRAVITYTERMVISITOR_H
