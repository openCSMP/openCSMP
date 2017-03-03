#ifndef COMPUTEGRAVITYTERMVISITOR_H
#define COMPUTEGRAVITYTERMVISITOR_H

#include "Visitor.h"
#include "FiniteElementManager.h"

namespace csmp{

struct Index;
class ScalarVariable;
template<size_t> class VectorVariable;
template<size_t> class Model;
template<size_t> class TwoPhaseModel;

template<size_t dim>
class ComputeGravityTermVisitor : public Visitor<dim>
{
  public:
    ComputeGravityTermVisitor( Model<dim>& model,
                               TwoPhaseModel<dim>& saturationFunctions,
                               const char* gravityVectorTag,
                               const char* permeabilityTag );
    
    ComputeGravityTermVisitor( Model<dim>& model,
                               const char* gravityVectorTag,
                               const char* permeabilityTag, 
                               const char* viscosityTag,
                               const char* densityTag );

    ComputeGravityTermVisitor( Model<dim>& model,
                               Index gravityVector,
                               Index permeability,
                               Index viscosity,
                               Index density );

    virtual ~ComputeGravityTermVisitor() {}

    virtual void Visit(Element<dim>* element);
    void ComputeContribution(Element<dim>* element);
    virtual void Visit(Model<dim>* model);
    virtual void Visit(Region<dim>* region);

  private:
    Model<dim>& model_;
    TwoPhaseModel<dim>* saturationFunctions_;
    Index permeabilityKey_, singlePhaseViscosityKey_, singlePhaseDensityKey_, gravityVectorKey_;
    const double gravityAcc_;

#if defined(_OPENMP )
    std::vector<FiniteElementManager> femgrs_; // one manager per thread
#endif

};

} //csmp

#endif // COMPUTEGRAVITYTERMVISITOR_H
