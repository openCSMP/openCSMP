#ifndef NODE_CENTERED_FINITE_VOLUME_ALGORITHM_H
#define NODE_CENTERED_FINITE_VOLUME_ALGORITHM_H

#include "CSMP_definitions.h"
#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "LinearSolver.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Element;
template<size_t> class TwoPhaseModel;
template<size_t> class Region;
template<size_t> struct StencilProcessor; // base class for element-based FV computations
class Solver;

template<size_t dim>
class NodeCenteredFiniteVolumeAlgorithm {
  public:
    explicit NodeCenteredFiniteVolumeAlgorithm( Region<dim>& );
    ~NodeCenteredFiniteVolumeAlgorithm();

    void ResetLHS( size_t finite_volumes );
    void ResetRHS( size_t finite_volumes );
    void ResetRow( size_t nid);
    void ZeroRHS( size_t nid );

    void AddToLHS( size_t i, size_t j, double addition );
    void AddToRHS( size_t i, double addition );
    void AddToRESULT( size_t i, double addition );
    void AssignRESULT( size_t i, double result );
    double   Result( size_t i ) const;
    std::vector<double>&  ResultVector();
    std::vector<double>&  RHSVector();

    double  LHS_Value( size_t i,  size_t j ) const;


    typename std::vector<double>::iterator  ResultsBegin();
    typename std::vector<double>::iterator  ResultsEnd();

    /// FV volume terms in diagonal & first-order terms=inter-FV fluxes in off-diagonal
    void AccumulateLHS( const StencilProcessor<dim>& es, double time_multiplier );

    /// Accumulation Process for Nonlinear Systems

    void AccumulateMatrix_NonlinearNewtonRaphson( const StencilProcessor<dim>&,
                                                  std::vector<double>& SAT0, double time_multiplier );

    void AccumulateResidual_NonlinearNewtonRaphson( const StencilProcessor<dim>&, std::vector<double>& SAT0,
                                                    double time_multiplier );

    void AccumulateMatrixAtBoundary_NonlinearNewtonRaphson( const StencilProcessor<dim>&,
                                                            std::vector<double>& SAT0, size_t nid,size_t pnid,
                                                            double time_multiplier, const csmp::Index& adv1_key );

    void AccumulateResidualAtBoundary_NonlinearNewtonRaphson( const StencilProcessor<dim>&,
                                                              std::vector<double>& SAT0, size_t nid,size_t pnid,
                                                              double time_multiplier, const csmp::Index& adv1_key );

    void CompensateInflowOutFlowBoundaries( size_t nid,
                                            const double& inflow,
                                            const double& flux_balance);

    /// Backward-Euler righthandside (psi at t0 * pore-volume) / dt for non-iterative computations
    void AccumulateRHS( const StencilProcessor<dim>& es, double time_multiplier );

    /// for iterative computations
    void AccumulateRHS( const StencilProcessor<dim>& es, std::vector<double>& SAT0, double time_multiplier );

    /// higher-order solution is added from RHS (first order in time)
    void AccumulateHigherOrderRHS( const StencilProcessor<dim>& );

    /// higher-order solution is subtracted from RHS (second order in time)
    void AccumulateHigherOrderRHS( const StencilProcessor<dim>&,
                                   const std::vector<std::vector<double> >& FACETFLUXES0,
                                   const std::vector<std::vector<double> >& LTDSATS0 );

      /// first-order Backward Euler solution is added to RHS
    void AddToRHS( const StencilProcessor<dim>& );
    void AddToRHS( const StencilProcessor<dim>&, std::vector<double>& SAT0);

    void SubtractNonConservativeOutFluxesFromRHS( TwoPhaseModel<dim>&,
                                                  const std::vector<double>& FLUX_BALANCE,
                                                  const csmp::Index& ad_key );

    /// diffusion and other divergent fluxes
    void AccumulateIntegral_DNT_op_DN_dV_LHS( const StencilProcessor<dim>& es );

    void AccumulateIntegral_DN_op_dS_LHS( const StencilProcessor<dim>& es );

    void AccumulateSectorSourceTermsInLHS( const StencilProcessor<dim>& );

    /// accumulates FV volume terms and first-order terms=inter-FV fluxes
    void AccumulateFluxUpwindSaturationProducts( const StencilProcessor<dim>& es );

    /// higher-order version for the case when fractional flows are considered
    void AccumulateHigherOrderFluxSaturationProducts( const StencilProcessor<dim>& es );

    void SolveMatrixEquation();

    /// results are stored back in the Region where condition flag is PLAIN
    double OutputResults(const PropertyDatabase<dim>&,
                           const csmp::Index& adv_key,
                           bool show_range ,
                           const size_t var_comp_nr=0) const;

    /// multiphase version for range 0..1 (adv1=wetting phase=1)
    double OutputResults( const PropertyDatabase<dim>&,
                            size_t result_phase, // enter either 1(w) or 2(nw) here
                            const csmp::Index& adv1_key,
                            const csmp::Index& adv2_key,
                            bool show_range ) const;


    /// testing only: results are stored back in the Region where condition flag is PLAIN
    double OutputResultsWithL2NormRes( const PropertyDatabase<dim>&, const csmp::Index& adv_key, bool show_range ) const;

    /// multiphase version for range 0..1 (adv1=wetting phase=1)
    double OutputResults_NonlinearNewtonRaphson( const PropertyDatabase<dim>&,
                              size_t result_phase, // enter either 1(w) or 2(nw) here
                              const csmp::Index& adv1_key,
                              const csmp::Index& adv2_key,
                              bool show_range );

    /// multiphase version for range 0..1 (adv1=wetting phase=1)
    int32_t OutputResults_NonlinearNewtonRaphson( const PropertyDatabase<dim>&,
                              size_t result_phase, // enter either 1(w) or 2(nw) here
                              const csmp::Index& adv1_key,
                              const csmp::Index& adv2_key,
                              bool show_range,
                              std::vector<double>&,
                              std::vector<double>&);

    /// interface to solver settings and instance
    CSMP_DEFAULT_LINEAR_SOLVER_SETTINGS&  GetSolverSettings();
    CSMP_DEFAULT_LINEAR_SOLVER* GetSolver();

    // Added Julian 23-09-2013
    void Verbose(bool verbose){this->verbose_=verbose;}
    bool Verbose(){return this->verbose_;}
    bool Verbose() const {return this->verbose_;}

  private:
    Region<dim>&           gref_;
    SparseMatrix           LHS;
    std::vector<double>  RHS, RESULT;
    DenseMatrix<DM_MIN>    DN, DNT;
    Solver*                solver_;
    const uint32_t           MAX_NODES_GAUSS_SOLVER;
    bool                   firstCall_;
    bool                   verbose_;
};


template<size_t dim>
inline void NodeCenteredFiniteVolumeAlgorithm<dim>::AddToLHS( size_t i, size_t j, double addition )
 { LHS.Add( i, j, addition ); }

template<size_t dim>
inline void NodeCenteredFiniteVolumeAlgorithm<dim>::AddToRHS( size_t i, double addition )
 { RHS[i] += addition; }

template<size_t dim>
inline void NodeCenteredFiniteVolumeAlgorithm<dim>::AddToRESULT( size_t i, double addition )
 { RESULT[i] += addition; }

template<size_t dim>
inline void NodeCenteredFiniteVolumeAlgorithm<dim>::AssignRESULT( size_t i, double fvresult )
 { RESULT[i] = fvresult; }

template<size_t dim>
inline double NodeCenteredFiniteVolumeAlgorithm<dim>::Result( size_t i ) const
 { return RESULT[i]; }

template<size_t dim>
inline std::vector<double>&  NodeCenteredFiniteVolumeAlgorithm<dim>::ResultVector()
 { return RESULT; }

template<size_t dim>
inline std::vector<double>&  NodeCenteredFiniteVolumeAlgorithm<dim>::RHSVector()
 { return RHS; }

template<size_t dim>
inline typename std::vector<double>::iterator  NodeCenteredFiniteVolumeAlgorithm<dim>::ResultsBegin()
 { return RESULT.begin(); }

template<size_t dim>
inline typename std::vector<double>::iterator  NodeCenteredFiniteVolumeAlgorithm<dim>::ResultsEnd()
 { return RESULT.end(); }

template<size_t dim>
inline double NodeCenteredFiniteVolumeAlgorithm<dim>::LHS_Value( size_t i,  size_t j ) const
 { return LHS(i,j); }



} // end namespace csmp

#endif























