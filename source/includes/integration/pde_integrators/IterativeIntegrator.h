#ifndef CSMP_ITERATIVE_INTEGRATOR_H
#define CSMP_ITERATIVE_INTEGRATOR_H

#include "CSMP_definitions.h"
#include "PDE_Integrator.h"

#ifdef CSMP_WITH_SAMG_SOLVER
namespace csmp {
  class SAMG_Settings;
}
#endif

namespace csmp {

template<size_t> class Interrelation;
template<size_t> class Visitor;
template<size_t> class Model;

/**
 
@brief Base class for the iterative itegration of PDEs using the FE method.

@author Adrian Burri
@date 2004
 
Provide a base class for transient algorithms which involve an iterative
solution algorithm on every time-step. 

@section design Design Intent

Make it easier to derive special iterative schemes. Basically, only the
method Residual needs to be overwritten.
 
@section participants Participants

The PDE_Integrator base class.
 
@section collaboration Collaboration

Models, Visitors and Interrelations. Calculations on single Region
objects is not supported.  

*/
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
class IterativeIntegrator : public PDE_Integrator<dim,COMPUTATION_DOMAIN> {
  public:
    IterativeIntegrator();
    
    void MaximalIterationNumber(size_t max_iter);
    size_t MaximalIterationNumber() const;
    void TargetResidual(double64 target_residual);
    void Verbose(bool yesno);
    
    void AddPostProcess( Interrelation<dim>* );
    void AddPostProcess( Visitor<dim>* );
    
    virtual void SetupEquations( COMPUTATION_DOMAIN<dim>& );
    virtual double64  Residual();
    virtual void SolveEquations( COMPUTATION_DOMAIN<dim>& );
    void ApplyPostProcesses( COMPUTATION_DOMAIN<dim>& );

#ifdef CSMP_WITH_SAMG_SOLVER
    virtual size_t Iterations( COMPUTATION_DOMAIN<dim>& );
#else
    /// add extra functionality for alternative solver if needed
#endif
    
  protected:
    /// Solution vector<double64> of last iteration * not yet used
    std::vector<double64> x_old_; 
    bool verbose_;
    
#ifdef CSMP_WITH_SAMG_SOLVER
    void SAMG_KeepMemory();
    void SAMG_KeepSettings();
#else
    /// add extra functionality for alternative solver if needed
#endif

  private:
    std::list<std::pair<Interrelation<dim>*, Visitor<dim>*> > processes_;
    size_t max_iter_;
    double64 target_residual_;
    
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings* GetSAMG_Settings();
#else
    /// add extra functionality for alternative solver if needed
#endif
};

} // end namespace csmp

#endif
