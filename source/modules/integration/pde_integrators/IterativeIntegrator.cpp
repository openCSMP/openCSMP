#include "IterativeIntegrator.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
IterativeIntegrator<dim,COMPUTATION_DOMAIN>::IterativeIntegrator()
 : verbose_(false),
   max_iter_(100),
   target_residual_(1.0e-06) 
 {}


template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
double IterativeIntegrator<dim,COMPUTATION_DOMAIN>::Residual()
 {
    cout << "\nIterativeIntegrator< dim>::Residual: this method has not been defined yet."<< endl;
    cout <<" returning 1."<< endl;
    return static_cast<double>(1.);
 }




template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::MaximalIterationNumber(size_t max_iter) {
  max_iter_ = max_iter;
}

template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::TargetResidual(double target_residual) {
  target_residual_ = target_residual;
}

template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::Verbose(bool yesno) {
  verbose_ = yesno;
}

template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::AddPostProcess( Interrelation<dim>* relation ) {
  assert(relation != 0);
  processes_.push_back(std::pair<Interrelation<dim>*, Visitor<dim>*>(relation, 0));
}

template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::AddPostProcess( Visitor<dim>* visitor ) {
  assert(visitor != 0);
  processes_.push_back(std::pair<Interrelation< dim>*, Visitor< dim>*>(0, visitor));
}

template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::SetupEquations( COMPUTATION_DOMAIN<dim>& gref ) {
  this->EstablishMatrixSetup(gref);
  this->Accumulate(gref);
  if (this->Transient()) {
      this->AssignInitialConditions(gref);
      this->LateAccumulate(gref);
   }
  this->AssignEssentialConditions(gref);
}


template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::SolveEquations( COMPUTATION_DOMAIN<dim>& gref ) {
  this->Solve();
  this->OutputResults(gref);
}


template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::ApplyPostProcesses( COMPUTATION_DOMAIN<dim>& gref ) {
  for ( typename list<std::pair<Interrelation<dim>*, Visitor<dim>*> >::iterator 
        it = processes_.begin(); it != processes_.end(); ++it ) {
      if (it->first != 0) gref.Apply(*(it->first));
      else                gref.Accept(*(it->second));
  }
  this->PostProcess(gref);
}


template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
size_t IterativeIntegrator<dim,COMPUTATION_DOMAIN>::MaximalIterationNumber() const {
  return max_iter_;
}


#ifdef CSMP_WITH_SAMG_SOLVER
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
size_t IterativeIntegrator<dim,COMPUTATION_DOMAIN>::Iterations( COMPUTATION_DOMAIN<dim>& sg )
 {
   SAMG_KeepMemory();
   double res = std::numeric_limits<double>::quiet_NaN(); 
   size_t i; 
  
  for ( i{0U}; i < max_iter_; ++i) {
    SetupEquations(sg);
    res = Residual();
    if (res < target_residual_) {
      if (verbose_) {
        cout << "Convergence after " << i << " iterations" << endl; 
        cout << "Final residual: " << res << endl;
      }
      break;
    }
    if (verbose_) {
      cout << "Iteration " << i << ": Residual = " << res << endl;
    }
    SolveEquations(sg);
    ApplyPostProcesses(sg);
    
    SAMG_KeepSettings(); 
  }
  
  ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    
  if ( i == max_iter_ && res > target_residual_ ) {
       csmp_error.notice( ERROR, "IterativeIntegrator<dim>::Iterations", "Solution did not converge! " );
      if ( verbose_ ) {
        cout << "\nIterativeIntegrator< dim>::Iterations: After: " << max_iter_ << " iterations, residual is: " << res << " and target: " << target_residual_ << endl;
        cout << "\n************* SOLUTION DID NOT CONVERGE *************" << endl; 
      }
    }
  return i;
}



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::SAMG_KeepMemory() {
  SAMG_Settings* settings = GetSAMG_Settings();
  if (settings) {
    settings->Set_iswit(4);
  }
}



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void IterativeIntegrator<dim,COMPUTATION_DOMAIN>::SAMG_KeepSettings() {
  SAMG_Settings* settings = GetSAMG_Settings();
  if (settings) {
    settings->Set_iswit(3);
  }
}



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
SAMG_Settings* IterativeIntegrator<dim,COMPUTATION_DOMAIN>::GetSAMG_Settings() {
  SAMG_Settings* result(0);
  if (SAMG_Solver* samg_solver = dynamic_cast<SAMG_Solver*>(this->GetSolver())) {
    result = dynamic_cast<SAMG_Settings*>(samg_solver->GetSolverSettings());
  }
  return result;
}

#endif // with SAMG_Solver

template class IterativeIntegrator<1U,Region>;
template class IterativeIntegrator<2U,Region>;
template class IterativeIntegrator<3U,Region>;

template class IterativeIntegrator<1U,Boundary>;
template class IterativeIntegrator<2U,Boundary>;
template class IterativeIntegrator<3U,Boundary>;

template class IterativeIntegrator<1U,SplitBoundary>;
template class IterativeIntegrator<2U,SplitBoundary>;
template class IterativeIntegrator<3U,SplitBoundary>;

} // end namespace csmp
