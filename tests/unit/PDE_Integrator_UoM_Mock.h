#ifndef PDE_INTEGRATOR_UoM_MOCK_H
#define PDE_INTEGRATOR_UoM_MOCK_H

#include "PDE_Integrator_UoM.h"

namespace csmp {

class Solver;

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN> class PDE_Integrator_UoM_Mock : 
    public PDE_Integrator_UoM <dim, INTEGRATION_DOMAIN> 
  {
  public:
    explicit PDE_Integrator_UoM_Mock( Solver& );
    virtual ~PDE_Integrator_UoM_Mock();
    
    virtual void EstablishMatrixSetup(const INTEGRATION_DOMAIN<dim>&);
    virtual void Accumulate(const INTEGRATION_DOMAIN<dim>&);
    virtual void LateAccumulate(const INTEGRATION_DOMAIN<dim>&);
    virtual void AssignInitialConditions(const INTEGRATION_DOMAIN<dim>&);
    virtual void OutputResults(INTEGRATION_DOMAIN<dim>&);
	  virtual void AssignEssentialConditions(INTEGRATION_DOMAIN<dim>&);
   
   // testing
    void EstablishMatrixSetupTest(const INTEGRATION_DOMAIN<dim>&);
    void AccumulateTest(const INTEGRATION_DOMAIN<dim>&);
    void LateAccumulateTest(const INTEGRATION_DOMAIN<dim>&);
    void AssignInitialConditionsTest(const INTEGRATION_DOMAIN<dim>&);
    void OutputResultsTest(INTEGRATION_DOMAIN<dim>&);
    void AssignEssentialConditionsTest(INTEGRATION_DOMAIN<dim>&);

    SparseMatrix const*           GetG() const;
    std::vector<double> const*  GetRH() const;
    std::vector<double>*        GetX();
    const std::vector<size_t>&    GetDOFIndex() const;
    
    friend class Transient_Test;
    friend class PDE_Integrator_UoM_Test;

  private:
    void EnumerateAndFixMatrixSize(const INTEGRATION_DOMAIN<dim>&);
 };

} // end namespace csmp

#endif
