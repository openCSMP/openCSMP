#ifndef SAMG_SOLVER_H
#define SAMG_SOLVER_H

#include "Solver.h"
#include "SparseMatrix.h"
#if defined(_OPENMP )
#include "CSRMatrix.h"
#else
#include "CompressedRowMatrix.h"
#endif

namespace csmp {

  class SAMG_Settings;

/** 

@brief Subclass of Solver for linear algebra computations conducted with SAMG, the algebraic multigrid solver for systems of equation from SCAPOS, Fraunhofer, Germany.

  Solves the linear algebraic system A x = b using the 
  algebraic multigrid solver DLL from SCAI, Fraunhofer Gesellschaft (St. Augustin, Bonn),
  which can be obtained from SCAPOS.
  
  Up to 2 istances of the DLL can be communicated with.
  Use another instance of the solver if you have multiple PDEs that are 
  solved sequentially and you want to remember SAMG settings from step to step. 

  @section segFault Segmentation Fault on some linux based os

  We have noticed invalid memory access when using dynamically allocated SAMG_Solver objects
  and their calls to SAMG. To avoid that, use static(stack) instantiation. Right now we have no idea
  why, but we're working on that...
  @code
  // i.e.:
  SAMG_Settings settings;
  settings.Set_napproach(2);
  SAMG_Solver samgSolver(&settings);
  PDE_Integrator<DIM,Region> pressure_diffusion(&samgSolver);
  @endcode

  @todo (1) Check line 254: is resetting res_in allowed if we use relative convergence
  @todo (3) Convergence check 670: 679 check eps or fabs(eps)?
  @todo (1) Look into runtime debug error when solving systems with vector as variable 
  @todo (3) Try mpl::is_volatile  type trait to check whether they should manage the pointers they are supplied with in the constructor (B)

*/
  class SAMG_Solver : public Solver {
    public:
      /// this method will not memory manage the SAMG_Settings object
      explicit SAMG_Solver( SAMG_Settings* settings );
      SAMG_Solver();
      SAMG_Solver( const SAMG_Solver& );
      virtual ~SAMG_Solver();
      SAMG_Solver& operator=( const SAMG_Solver& );
      
      void InputSolverSettings( SAMG_Settings* settings );
      virtual SolverSettings* GetSolverSettings();

      bool Write_SAMG_TextInputFile( const char* filename ) const;
      void Write_SAMG_TextOutput( bool write );
	    	    
      double64 LastSolverResidual() const;

      void CheckSparsityCriterion( int32& levelx ) const;
      void CheckCycleCriterion( int32& iswtch ) const;
      void UpdateCycleCriterion( int32 iswtch );
      bool CheckConvergence( double64 eps ) const;

      void OutputVectors() const;

    protected:
      virtual void SolveMatrixEquation( SparseMatrix& A,
                                        std::vector<double64>& b,
                                        std::vector<double64>& x,
                                        size_t no_unknowns );

      virtual void SolveMatrixEquation( CompressedRowMatrix& A,
                                        std::vector<double64>& b,
                                        std::vector<double64>& x,
                                        size_t no_unknowns );
      
    private:
      SAMG_Settings*        settings_;

      int32                 nsys_;    /// < Number of unknowns
      int32                 npnts_;   /// < Number of points (in mesh)
      int32                 nnu_;	    /// < Number of variables >= 1 (matrix size nnu^2)
      int32                 nna_;	    /// < Number of matrix entries stored in vector a >= nnu

      std::vector<double64> u_;	      /// < [] First guess solution to A.u=f
      std::vector<double64> f_;	      /// < [] Right hand side

      // For coupled systems
      std::vector<int32>    iscale_;  /// < [] which unknowns require scaling (nsys >= 1)
                                      /// < iscale(k)=0	means k-th unknown does NOT require scaling
                                      /// < iscale(k)=1	means k-th unknown DOES require scaling
      std::vector<int32>    iu_;	    /// < [] variable to unknown pointer
      int                   ndiu_;    /// < if nsys>1, m_diu=m_nnu
      std::vector<int32>    ip_;	    /// < []
      int                   ndip_;

#if defined(_OPENMP )
      CSRMatrix crmat_;
#else
      CompressedRowMatrix   crmat_;

#endif

      // Output parameters
      double64  res_in_;      /// < residual of first guess
      double64  res_out_;     /// < residual of final approximation
      int32     ncyc_done_;   /// < total number cycles (iterations) performed
      int32     ncyc_best_;   /// < stores lowest number of iterations achieved

      int32     ierr_;        /// < error inidicator (= 0 no error, > 0 fatal error, < 0 completed with warnings)

      const bool  newed_SAMG_Settings_object; 
      bool output_amg_data_to_text_files_;
};

} // end namespace csmp


void USER_coo( int* i, int* ndim, double* x, double* y, double* z );

#endif
 
