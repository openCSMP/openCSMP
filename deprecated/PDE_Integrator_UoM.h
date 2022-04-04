#ifndef PDE_INTEGRATOR_UOM_H
#define PDE_INTEGRATOR_UOM_H

#include "PDE_Integrator.h"

namespace csmp {

/**
    PDE integrator version that assembles the solution matrix, eliminating Dirichlet constraints
    from the matrix and righthand vector before the accumulation.
 
    @attention works only for scalar solution variables.
 
    @author Khoa Luat Tran (University of Melbourne)
    @date 17/4/2019
*/
template<size_t dim, template<size_t> class COMPUTATION_DOMAIN>
class PDE_Integrator_UoM : public PDE_Integrator<dim, COMPUTATION_DOMAIN> {
  public:
    /// only use this constructor
    explicit PDE_Integrator_UoM( Solver& );
    virtual ~PDE_Integrator_UoM();
    
    virtual void  IntegrateOver( COMPUTATION_DOMAIN<dim>& ); // Luat Khoa Tran
  
  protected:
    virtual void  EstablishMatrixSetup( const COMPUTATION_DOMAIN<dim>& );
    virtual void  AssignInitialConditions( const COMPUTATION_DOMAIN<dim>& );
    virtual void  Accumulate( const COMPUTATION_DOMAIN<dim>& );
    virtual void  LateAccumulate( const COMPUTATION_DOMAIN<dim>& );
	  virtual void  AssignEssentialConditions( const COMPUTATION_DOMAIN<dim>& );
    virtual void  OutputResults( COMPUTATION_DOMAIN<dim>& );
    // For Solve() and PostProcess() the base class methods are used!

  protected:
 	  std::vector<size_t>    DOF_indexes_; ///< for indexing DOFs (only non-Dirichlet BC dofs, enumerated 0 -> maximum DOF
	  std::vector<double64>  pivotVector_; ///< terms recovered from eliminated rows
    void  EnumerateAndFixMatrixSize( const COMPUTATION_DOMAIN<dim>& );  // Luat Khoa Tran

  private:
    double64               scale_factor_; ///< for essential conditions
    bool                   verbose_;
  };

} // csmp

#endif
