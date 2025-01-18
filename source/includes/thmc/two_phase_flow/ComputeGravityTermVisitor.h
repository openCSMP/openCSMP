#ifndef COMPUTEGRAVITYTERMVISITOR_H
#define COMPUTEGRAVITYTERMVISITOR_H

#include "Visitor.h"
#include "Index.h"
 
namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class TwoPhaseModel;

template<uint32_t dim>
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

    virtual void Visit( Element<dim>* element);

  private:
    Model<dim>& model_;
    TwoPhaseModel<dim>* saturationFunctions_;
    Index permeabilityKey_, singlePhaseViscosityKey_, singlePhaseDensityKey_, gravityVectorKey_;
    const double gravityAcc_;
};

} //csmp

#endif // COMPUTEGRAVITYTERMVISITOR_H
