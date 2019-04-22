#ifndef PDE_INTEGRATOR_UoM_MOCK_H
#define PDE_INTEGRATOR_UoM_MOCK_H

#include "PDE_Integrator_UoM.h"

namespace csmp {

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN> class PDE_Integrator_UoM_Mock : 
    public PDE_Integrator_UoM <dim, INTEGRATION_DOMAIN> 
  {
  public:
    PDE_Integrator_UoM_Mock();
    virtual ~PDE_Integrator_UoM_Mock();
    SparseMatrix const* GetG() const;
    std::vector<double64> const* GetRH() const;
    std::vector<double64>*  GetX();
    const std::vector<size_t>&    GetDOFIndex() const;
    void EstablishMatrixSetup(const INTEGRATION_DOMAIN<dim>&);
    void EstablishMatrixSetupTest(const INTEGRATION_DOMAIN<dim>&);
    void EnumerateAndFixMatrixSize(const INTEGRATION_DOMAIN<dim>&);
    void Accumulate(const INTEGRATION_DOMAIN<dim>&);
    void AccumulateTest(const INTEGRATION_DOMAIN<dim>&);
    void LateAccumulate(const INTEGRATION_DOMAIN<dim>&);
    void LateAccumulateTest(const INTEGRATION_DOMAIN<dim>&);
    void AssignInitialConditions(const INTEGRATION_DOMAIN<dim>&);
    void AssignInitialConditionsTest(const INTEGRATION_DOMAIN<dim>&);
    void OutputResults(INTEGRATION_DOMAIN<dim>&);
    void OutputResultsTest(INTEGRATION_DOMAIN<dim>&);
	void AssignEssentialConditions(INTEGRATION_DOMAIN<dim>&); //AssignEssentialConditions
	void AssignEssentialConditionsTest(INTEGRATION_DOMAIN<dim>&); //AssignEssentialConditions
	
  };

} // end namespace csmp

#endif
