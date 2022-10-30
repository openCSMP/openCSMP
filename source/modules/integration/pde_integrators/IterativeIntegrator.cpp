#include "IterativeIntegrator.h"
#include "PDE_Integrator.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "ErrorHandler.h"
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"

using namespace std;

namespace csmp {

template<uint32_t dim,template<uint32_t> class CELLTYPE>
IterativeIntegrator<dim,CELLTYPE>::IterativeIntegrator()
 : // TODO: get the solver from somewhere! PDE_Integrator<dim,CELLTYPE>( solver ),
   verbose_(false),
   max_iter_(100),
   target_residual_(1.0e-06) 
 {}


template<uint32_t dim,template<uint32_t> class CELLTYPE>
double IterativeIntegrator<dim,CELLTYPE>::Residual()
 {
    cout << "\nIterativeIntegrator< dim>::Residual: this method has not been defined yet."<< endl;
    cout <<" returning 1."<< endl;
    return static_cast<double>(1.);
 }




template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::MaximalIterationNumber(size_t max_iter) {
  max_iter_ = max_iter;
}

template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::TargetResidual(double target_residual) {
  target_residual_ = target_residual;
}

template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::Verbose(bool yesno) {
  verbose_ = yesno;
}

template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::AddPostProcess( Interrelation<dim>* relation ) {
  assert(relation != 0);
  processes_.push_back(std::pair<Interrelation<dim>*, Visitor<dim>*>(relation, 0));
}

template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::AddPostProcess( Visitor<dim>* visitor ) {
  assert(visitor != 0);
  processes_.push_back(std::pair<Interrelation< dim>*, Visitor< dim>*>(0, visitor));
}

template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::SetupEquations( ModelSubDomain<dim,CELLTYPE>& gref ) {
  this->EstablishMatrixSetup(gref);
  this->Accumulate(gref);
  if (this->Transient()) {
      this->AssignInitialConditions(gref);
      this->LateAccumulate(gref);
   }
  this->AssignEssentialConditions(gref);
}


template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::SolveEquations( ModelSubDomain<dim,CELLTYPE>& gref ) {
  this->Solve();
  this->OutputResults(gref);
}


template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::ApplyPostProcesses( ModelSubDomain<dim,CELLTYPE>& gref ) {
  for ( typename list<std::pair<Interrelation<dim>*, Visitor<dim>*> >::iterator 
        it = processes_.begin(); it != processes_.end(); ++it ) {
      if (it->first != 0) gref.Apply(*(it->first));
      else                gref.Accept(*(it->second));
  }
  this->PostProcess(gref);
}


template<uint32_t dim,template<uint32_t> class CELLTYPE>
size_t IterativeIntegrator<dim,CELLTYPE>::MaximalIterationNumber() const {
  return max_iter_;
}


#ifdef CSMP_WITH_SAMG_SOLVER
template<uint32_t dim,template<uint32_t> class CELLTYPE>
size_t IterativeIntegrator<dim,CELLTYPE>::Iterations( ModelSubDomain<dim,CELLTYPE>& sg )
 {
   SAMG_KeepMemory();
   double res = std::numeric_limits<double>::quiet_NaN(); 
   uint32_t i{0U}; 
  
  for ( i = 0U; i < max_iter_; ++i) {
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
       csmp_error.Note( ERROR, "IterativeIntegrator<dim>::Iterations", "Solution did not converge! " );
      if ( verbose_ ) {
        cout << "\nIterativeIntegrator< dim>::Iterations: After: " << max_iter_ << " iterations, residual is: " << res << " and target: " << target_residual_ << endl;
        cout << "\n************* SOLUTION DID NOT CONVERGE *************" << endl; 
      }
    }
  return i;
}



template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::SAMG_KeepMemory() {
  SAMG_Settings* settings = GetSAMG_Settings();
  if (settings) {
    settings->Set_iswit(4);
  }
}



template<uint32_t dim,template<uint32_t> class CELLTYPE>
void IterativeIntegrator<dim,CELLTYPE>::SAMG_KeepSettings() {
  SAMG_Settings* settings = GetSAMG_Settings();
  if (settings) {
    settings->Set_iswit(3);
  }
}



template<uint32_t dim,template<uint32_t> class CELLTYPE>
SAMG_Settings* IterativeIntegrator<dim,CELLTYPE>::GetSAMG_Settings() {
  SAMG_Solver& solver = dynamic_cast<SAMG_Solver&>(this->GetSolver());
  return dynamic_cast<SAMG_Settings*>(solver.GetSolverSettings());
}

#endif // with SAMG_Solver

template class IterativeIntegrator<1U>;
template class IterativeIntegrator<2U>;
template class IterativeIntegrator<3U>;

template class IterativeIntegrator<1U,Face>;
template class IterativeIntegrator<2U,Face>;
template class IterativeIntegrator<3U,Face>;

} // end namespace csmp
