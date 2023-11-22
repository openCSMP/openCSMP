#ifndef SAMGp_Solver_H
#define SAMGp_Solver_H

#ifndef HAVE_MPI_CXX
#define HAVE_MPI_CXX
#endif

#ifdef MPICH_SKIP_MPICXX
#undef MPICH_SKIP_MPICXX
#endif

#include "mpi.h"
#include "Solver.h"
#include "samgp.h"
#include "CSP_SparseMatrix.h"
#include "SAMG_Settings.h"
#include "SAMGp_CommunicationData.h"
#include "CompressedRowMatrix.h"

void USER_coo(int *i, int *ndim, double* x, double* y, double* z);

namespace csp {

template<typename fT, stl_index dim>
  class SAMGp_Solver : public Solver {
    public:
      SAMGp_Solver();
      ~SAMGp_Solver();
      explicit SAMGp_Solver(SAMG_Settings* settings, SAMGp_CommunicationData<fT,dim>* parameters, bool update_halo=false);
      
      bool Write_SAMG_TextInputFile( const char* filename ) const;
	    void Write_SAMG_TextOutput( bool write );
	    
	    void SetSolverSettings(SAMG_Settings* settings);
	    SAMG_Settings* GetSolverSettings() const;
	    
	    csp_float LastSolverResidual() const;

    protected:
      virtual void SolveMatrixEquation(SparseMatrix& A, 
                                       std::vector<csp_float>& b, 
                                       std::vector<csp_float>& x,
                                       stl_index no_unknowns);
    private:
      SAMGp_CommunicationData<fT,dim>* samgp_parameters_;
      SAMG_Settings* settings_;
      
      bool    update_halo_; // if true, after conversion of SAMGp, the data on
                            // outerhalo nodes will be communicated.
      // variables necessary for communication
      std::vector<csp_float>      isnddata;
      std::vector<csp_float>      irecdata;
      int tag;
      int my_rank;
      std::vector<MPI::Request>   mpi_req;

      int     nsys_;    // Number of unknowns
      int     npnts_;   // Number of points (in mesh)
      int     nnu_;	// Number of variables >= 1 (matrix size nnu^2)
      int     nna_;	// Number of matrix entries stored in vector a >= nnu

      std::vector<double> u_;	// [] First guess solution to Au=f
      std::vector<double> f_;	// [] Right hand side
      // For coupled systems
      std::vector<int>    iscale_;  // [] which unknowns require scaling (nsys >= 1) 
                                    // iscale(k)=0	means k-th unknown does NOT require scaling
                                    // iscale(k)=1	means k-th unknown DOES require scaling
      std::vector<int>    iu_;	    // [] variable to unknown pointer
      int                 ndiu_;    // if nsys>1, m_diu=m_nnu
      std::vector<int>    ip_;	    // []
      int                 ndip_;

      int                 nshalo, npsnd, nrhalo, nprec;

      CompressedRowMatrix crmat_;

      // Output parameters 
      double  res_in_;	    // residual of first guess
      double  res_out_;     // residual of final approximation
      int     ncyc_done_;   // total # cycles (iterations) performed
      int     ierr_;	    // error inidicator (=0 no error, >0 fatal error, <0 completed with warnings)
	    
      bool output_amg_data_to_text_files_;
	    
      void CheckConvergence(double64 eps);
      void CommunicateOuterhalo();
      void OutputVectors();

};

} // end namespace csp

#endif
