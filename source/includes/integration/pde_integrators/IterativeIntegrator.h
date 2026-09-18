// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CSMP_ITERATIVE_INTEGRATOR_H
#define CSMP_ITERATIVE_INTEGRATOR_H

#include "PDE_Integrator.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Interrelation;
template<uint32_t> class Visitor;
template<uint32_t> class Model;
class SAMG_Settings;

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
template<uint32_t dim,template<uint32_t> class CELLTYPE=Element>
class IterativeIntegrator : public PDE_Integrator<dim,CELLTYPE> {
  public:
    IterativeIntegrator();
    
    void MaximalIterationNumber(size_t max_iter);
    size_t MaximalIterationNumber() const;
    void TargetResidual(double target_residual);
    void Verbose(bool yesno);
    
    void AddPostProcess( Interrelation<dim>* );
    void AddPostProcess( Visitor<dim>* );
    
    virtual void SetupEquations( ModelSubDomain<dim,CELLTYPE>& );
    virtual double  Residual();
    virtual void SolveEquations( ModelSubDomain<dim,CELLTYPE>& );
    void ApplyPostProcesses( ModelSubDomain<dim,CELLTYPE>& );

#ifdef CSMP_WITH_SAMG_SOLVER
    virtual size_t Iterations( ModelSubDomain<dim,CELLTYPE>& );
#else
    /// add extra functionality for alternative solver if needed
#endif
    
  protected:
    /// Solution vector<double> of last iteration * not yet used
    std::vector<double> x_old_; 
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
    double target_residual_;
    
#ifdef CSMP_WITH_SAMG_SOLVER
    SAMG_Settings* GetSAMG_Settings();
#else
    /// add extra functionality for alternative solver if needed
#endif
};

} // end namespace csmp

#endif
