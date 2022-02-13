#ifndef RELATIVE_PERMEABILITY_MODEL_TEST_H
#define RELATIVE_PERMEABILITY_MODEL_TEST_H

#include "Test.h"
#include "TwoPhaseModel.h"
#include "Model.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Model;

class RelativePermeabilityModel_Test : public Test {
  public:
    RelativePermeabilityModel_Test( bool verbose );
    virtual ~RelativePermeabilityModel_Test();
    virtual void run();
    void Test( TwoPhaseModel<1U>& relperm, bool extended_property_set );                      
    const PropertyDatabase<1>& Database() const;
      
  private:
    void InitializeFlowProperties();

    Model<1U>*  model_ptr_;
    const bool verbose_;
};

} 

#endif
