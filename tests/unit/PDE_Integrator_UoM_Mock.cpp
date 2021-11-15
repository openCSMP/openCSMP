#include "PDE_Integrator_UoM_Mock.h"
#include "Solver.h"

using namespace std;

namespace csmp {

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::PDE_Integrator_UoM_Mock( Solver& solver )
  : PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>(solver) {
       
  }
  
  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::~PDE_Integrator_UoM_Mock()  {
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  SparseMatrix const* PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::GetG() const {
    return &this->G_;
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  std::vector<double> const* PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::GetRH() const {
    return &this->rh_;
  }
  
  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  std::vector<double>* PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::GetX() {
	  return &this->x_;
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  const std::vector<size_t>& PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::GetDOFIndex() const {
    return this->DOF_indexes_;
  }			
  
  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::EstablishMatrixSetup(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator<dim, INTEGRATION_DOMAIN>::EstablishMatrixSetup(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::EstablishMatrixSetupTest(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>::EstablishMatrixSetup(domain);
  }  
  
  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::EnumerateAndFixMatrixSize(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>::EnumerateAndFixMatrixSize(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::Accumulate(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator<dim, INTEGRATION_DOMAIN>::Accumulate(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::AccumulateTest(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>::Accumulate(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::LateAccumulate(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator<dim, INTEGRATION_DOMAIN>::LateAccumulate(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::LateAccumulateTest(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>::LateAccumulate(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::AssignInitialConditions(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator<dim, INTEGRATION_DOMAIN>::AssignInitialConditions(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::AssignInitialConditionsTest(const INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>::AssignInitialConditions(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::AssignEssentialConditions(INTEGRATION_DOMAIN<dim>& domain) {
	PDE_Integrator<dim, INTEGRATION_DOMAIN>::AssignEssentialConditions(domain);
  }; //AssignEssentialConditions
  
  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::AssignEssentialConditionsTest(INTEGRATION_DOMAIN<dim>& domain) {
	PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>::AssignEssentialConditions(domain);
  }; //AssignEssentialConditions

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::OutputResults(INTEGRATION_DOMAIN<dim>& domain) {
	  PDE_Integrator<dim, INTEGRATION_DOMAIN>::OutputResults(domain);
  }

  template<size_t dim, template<size_t> class INTEGRATION_DOMAIN>
  void PDE_Integrator_UoM_Mock<dim, INTEGRATION_DOMAIN>::OutputResultsTest(INTEGRATION_DOMAIN<dim>& domain) {
    PDE_Integrator_UoM<dim, INTEGRATION_DOMAIN>::OutputResults(domain);
  }

  
  template class PDE_Integrator_UoM_Mock<1U, Region>;
  template class PDE_Integrator_UoM_Mock<2U, Region>;
  template class PDE_Integrator_UoM_Mock<3U, Region>;

  template class PDE_Integrator_UoM_Mock<1U, Boundary>;
  template class PDE_Integrator_UoM_Mock<2U, Boundary>;
  template class PDE_Integrator_UoM_Mock<3U, Boundary>;

  template class PDE_Integrator_UoM_Mock<1U, SplitBoundary>;
  template class PDE_Integrator_UoM_Mock<2U, SplitBoundary>;
  template class PDE_Integrator_UoM_Mock<3U, SplitBoundary>;


} // end namespace csmp
