#ifndef FRACTURECOMPUTEVELOCITYVISITOR_H
#define FRACTURECOMPUTEVELOCITYVISITOR_H

#include "Visitor.h"
#include "ScalarVariable.h"
#include "DenseMatrix.h"
#include "TwoPhaseModel.h"

namespace csmp{

    struct Index;
    template<uint32_t> class Model;


template<uint32_t dim>
class FractureComputeVelocityVisitor : public Visitor<dim>
  {
    public:
        FractureComputeVelocityVisitor( Model<dim>& model, TwoPhaseModel<dim>& saturationFunctions,
                                        const char* vt_Tag,
                                        const char* vn_Tag,
                                        const char* fluidPressureTag,
                                        const char* fractureCapillaryPressureTag,
                                        const char* volFluxTag,
                                        const char* prevVolFluxTag,
                                        const char* gravityVectorTag,
                                        const char* permeabilityTag,
                                        bool withCapillaryGradient,
                                        bool withGravity );

        virtual ~FractureComputeVelocityVisitor() {}

        virtual void Visit(Element<dim>* element);

    private:
        Model<dim>& model_;
        TwoPhaseModel<dim>* saturationFunctions_;
        Index vt_Key_, vn_Key_, fluidPressureKey_, totalMobilityKey_, volumeFluxKey_, previousVolumeFluxKey_,
              gravityVectorKey_, fracCapillaryPressureKey_, permeabilityKey_;
        VectorVariable<dim> velo_, gravityVector_;
        ScalarVariable flux_;
        DenseMatrix<DM_MIN> DERIV_;
        bool withCapillaryGradient_, withGravity_ ;
        const double gravityAcc_;

  };

} //csmp

#endif // FRACTURECOMPUTEVELOCITYVISITOR_H
