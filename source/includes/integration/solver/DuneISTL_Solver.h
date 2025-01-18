#ifndef DUNE_ISTL_SOLVER_H
#define DUNE_ISTL_SOLVER_H

#ifdef CSMP_WITH_DUNE_ISTL

#include <dune/common/exceptions.hh>
#include <dune/common/deprecated.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/typetraits.hh>

#include <dune/istl/basearray.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/vbvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/operators.hh>

#include <dune/istl/io.hh>
#include <dune/common/timer.hh>

#include <dune/istl/owneroverlapcopy.hh>
#include <dune/istl/solvercategory.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/paamg/amg.hh>
#include <dune/istl/paamg/pinfo.hh>
#include <dune/istl/superlu.hh>

// csmp related
#include "SparseMatrix.h"
#include "Solver.h"

namespace csmp {

/** Interface to ISTL ( Iterative Solver Template Library ) of Dune ( Distributed and Unified Numerics Environment )
 *  https://www.dune-project.org
 *  https://www.dune-project.org/coremodules.html
 *  Part of interface to istl solvers is taken from dune-pdelab
 *  https://www.dune-project.org/pdelab/index.html
 *  in particular from corresponding ISTLBackend's:
 *  dune/dune-pdelab/dune/pdelab/backend/seqistlsolverbackend.hh
 */

class DuneISTL_Settings
{
public:

    typedef Dune::Amg::Parameters AmgParameters;

    DuneISTL_Settings()
        :verbose(1),
         maxiter(5000),
         reduction(1E-12),
         n_(1),
         w_(1.0),
         params(15,2000),
         reuse(false),
         usesuperlu(true)
    {}

    DuneISTL_Settings( const DuneISTL_Settings& settings )
        :verbose   ( settings.verbose   ),
         maxiter   ( settings.maxiter   ),
         reduction ( settings.reduction ),
         n_        ( settings.n_ ),
         w_        ( settings.w_ ),
         params    ( settings.params ),
         reuse     ( settings.reuse ),
         usesuperlu( settings.usesuperlu )
    {
    }

    ~DuneISTL_Settings()
    {}

    DuneISTL_Settings& operator=( const DuneISTL_Settings& settings )
    {
        if( this != &settings )
        {
            verbose    = settings.verbose;
            maxiter    = settings.maxiter;
            reduction  = settings.reduction;
            n_         = settings.n_;
            w_         = settings.w_;
            params     = settings.params;
            reuse      = settings.reuse;
            usesuperlu = settings.usesuperlu;
        }
        return *this;
    }

    int      Get_verbose() const { return verbose; }
    unsigned Get_maxiter() const { return maxiter; }
    double   Get_reduction() const { return reduction; }
    double   Get_w() const { return w_; }
    int      Get_n() const { return n_; }
    const    AmgParameters& Get_amg_params() const { return params; }
    bool     Get_reuse() const { return reuse; }
    bool     Get_usesuperlu() const { return usesuperlu; }

    void     Set_verbose( int _verbose ) { verbose = _verbose; }
    void     Set_maxiter( unsigned _maxiter ) { maxiter = _maxiter; }
    void     Set_reduction( double _reduction ) { reduction = _reduction; }
    void     Set_w( double _w ){ w_ = _w; }
    void     Set_n( int _n ){ n_ = _n; }
    void     Set_amg_params( const AmgParameters& _params ){ params = _params; }
    void     Set_reuse( bool _reuse ) { reuse = _reuse; }
    void     Set_usesuperlu( bool _usesuperlu ) { usesuperlu = _usesuperlu; }

private:

    int      verbose;       // verbosity level
    unsigned maxiter;       // max number of iterations
    double   reduction;     // reduction
    double   w_;            // relaxation factor (ILUn)
    int      n_;            // level (ILUn)
    AmgParameters params;   // specific parameters for AMG
    bool     reuse;         // resue matrix
    bool     usesuperlu;    // use superlu solver
};

/**
@class DuneISTL_Settings DuneISTL_Settings "solver/DuneISTL_Solver.h"
@author Roman M.N.
@date 2015
*/

/*! \brief csmp solver interface to dune-istl solvers

  \param M the given matrix
  \param V the solution vector to be computed
  \param W right hand side
  \param DuneISTL_SolverInterface specific solver iterface
*/

template<class M, class V, class W,
         class DuneISTL_SolverInterface>
class DuneISTL_Solver : public Solver
{
  public:

    explicit DuneISTL_Solver( DuneISTL_Settings* settings ):
        settings_(settings),
        newed_DUNE_Settings_object_(false),
        solver_(settings_)
    {
        assert(settings!=NULL);
    }

    DuneISTL_Solver()
        :settings_(new DuneISTL_Settings()),
         newed_DUNE_Settings_object_(false),
         solver_(settings_)
    {}

    DuneISTL_Solver( const DuneISTL_Solver& solver )
        :settings_(new DuneISTL_Settings(*solver.settings_)),
         newed_DUNE_Settings_object_(true),
         solver_(settings_)
    {
    }

    virtual ~DuneISTL_Solver()
    {
        if ( newed_DUNE_Settings_object_ )
            delete settings_;
    }

    DuneISTL_Solver& operator=(const DuneISTL_Solver& solver)
    {
        if( this != &solver )
        {
            if ( newed_DUNE_Settings_object_ ) {
                delete settings_;
                settings_ = new DuneISTL_Settings(*solver.settings_);
            }
            else settings_ = solver.settings_;
        }
        return *this;
    }

    void SolverSettings( DuneISTL_Settings* settings )
    {
        assert(settings!=NULL);
        settings_ = settings;
    }

    DuneISTL_Settings* SolverSettings() const
    {
        return settings_;
    }

    void ConvertSparseMatrixToBCRSMatrix( const csmp::SparseMatrix& A, Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> >& B)
    {
        typedef Dune::FieldMatrix<double,1,1> MB;
        typedef Dune::BCRSMatrix<MB> Matrix;
        typedef MB::size_type size_type;

        // build matrix
        const size_type cols( A.Cols() );
        const size_type rows( A.Rows() );
        const size_type entries( A.Entries() );
        B.setSize( rows, cols, entries );
        B.setBuildMode( Dune::BCRSMatrix<MB>::row_wise );

        // create sparsity pattern
        typedef Dune::BCRSMatrix<MB>::CreateIterator Iter;
        for(Iter row=B.createbegin(); row!=B.createend(); ++row ){
            const auto i = row.index();
            // add nonzeros elements
            for(csmp::SparseMatrix::colsConstIterator cit = A.RowBegin( i ); cit != A.RowEnd( i ); ++cit ){
                const size_t j = (*cit).first;
                row.insert( j );
            }
        }
        // assign data
        for (Matrix::RowIterator i=B.begin(); i!=B.end(); ++i)
          for (Matrix::ColIterator j=(*i).begin(); j!=(*i).end(); ++j )
              (*j) = A.At( i.index(), j.index() );
    }
    void ConvertVectorToBlockVector( const std::vector<double>& v, Dune::BlockVector<Dune::FieldVector<double,1> >& bv)
    {
        typedef Dune::FieldVector<double,1> VB;
        typedef Dune::BlockVector<VB> Vector;
        typedef VB::size_type size_type;

        // build vector
        const size_type entries( v.size() );
        bv.resize(entries);
        // assign data
        for (size_t i{0U}; i<entries; ++i)
            bv[i] = v[i];
    }
    void ConvertBlockVectorToVector( const Dune::BlockVector<Dune::FieldVector<double,1> >& bv, std::vector<double>& v )
    {
        typedef Dune::FieldVector<double,1> VB;
        typedef Dune::BlockVector<VB> Vector;
        typedef VB::size_type size_type;

        // build vector
        const size_type entries( bv.size() );
        v.resize(entries);
        // assign data
        for (size_t i{0U}; i<entries; ++i)
            v[i] = bv[i];
    }

protected:

    virtual void SolveMatrixEquation( SparseMatrix& A,
                                    std::vector<double>& b,
                                    std::vector<double>& x,
                                    size_t no_unknowns )
    {
        // Todo: optimise conversion
        M lhs;
        ConvertSparseMatrixToBCRSMatrix(A,lhs);
        V sol;
        ConvertVectorToBlockVector(x,sol);
        W rhs;
        ConvertVectorToBlockVector(b,rhs);
        solver_.apply( lhs, sol, rhs, settings_->Get_reduction() );
        ConvertBlockVectorToVector(sol,x);
    }

  private:

    DuneISTL_Settings*        settings_;
    DuneISTL_SolverInterface  solver_;
    bool                      newed_DUNE_Settings_object_;

};

/**
@class DuneISTL_Solver DuneISTL_Solver "solver/DuneISTL_Solver.h"
@author Roman M.N.
@date 2015
*/



namespace istl {

// ********************************************************************************
// Helper functions for uniform access to ISTL containers
//
// The following suite of raw() functions should be used in places where an
// algorithm might work on either the bare ISTL container or the PDELab
// wrapper and has to access the bare container.
// ********************************************************************************

//! Returns the raw ISTL object associated with v, or v itself it is already an ISTL object.
template<typename V>
V& raw(V& v)
{
  return v;
}

//! Returns the raw ISTL object associated with v, or v itself it is already an ISTL object.
template<typename V>
const V& raw(const V& v)
{
  return v;
}

//! Returns the raw ISTL type associated with C, or C itself it is already an ISTL type.
template<typename C>
struct raw_type
{
  typedef C type;
};

}

struct SequentialNorm
{/*! \brief compute global norm of a vector

    \param[in] v the given vector
  */
  template<class V>
  typename Dune::template FieldTraits<typename V::ElementType >::real_type norm(const V& v) const
  {
    return istl::raw(v).two_norm();
  }
};

// Status information of a linear solver
template<class RFType>
struct LinearSolverResult
{
  bool converged;            // Solver converged
  unsigned int iterations;   // number of iterations
  double elapsed;            // total user time in seconds
  RFType reduction;          // defect reduction
  RFType conv_rate;          // convergence rate (average reduction per step)

  LinearSolverResult() :
    converged(false),
    iterations(0),
    elapsed(0.0),
    reduction(0.0),
    conv_rate(0.0)
  {}
};

class LinearResultStorage
{
public:
  /*! \brief Return access to result data */
  const LinearSolverResult<double>& result() const
  {
    return res;
  }

protected:

  LinearSolverResult<double> res;
};

//==============================================================================
// Here we add some standard linear solvers conforming to the linear solver
// interface required to solve linear and nonlinear problems.
//==============================================================================

//! Solver to be used for explicit time-steppers with (block-)diagonal mass matrix
class DuneISTL_SEQ_ExplicitDiagonal
  : public SequentialNorm, public LinearResultStorage
{
public:
  /*! \brief make a linear solver object
  */
  DuneISTL_SEQ_ExplicitDiagonal ()
  {}

  DuneISTL_SEQ_ExplicitDiagonal (DuneISTL_Settings* sets )
  {}

  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename W::value_type reduction)
  {
    Dune::SeqJac<M,V,W> jac(istl::raw(A),1,1.0);
    jac.pre(z,r);
    jac.apply(z,r);
    jac.post(z);
    res.converged  = true;
    res.iterations = 1;
    res.elapsed    = 0.0;
    res.reduction  = reduction;
    res.conv_rate  = reduction; // pow(reduction,1.0/1)
  }
};



/////////////////////////////////////////////////////////
///! Base class for Iterative Solver with Preconditioner
/////////////////////////////////////////////////////////

template<template<class,class,class,int> class Preconditioner,
         template<class> class Solver>
class DuneISTL_SEQ_Base
  : public SequentialNorm, public LinearResultStorage
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_Base(unsigned maxiter_=5000, int verbose_=1)
    : maxiter(maxiter_),
      verbose(verbose_)
  {}

  explicit DuneISTL_SEQ_Base(DuneISTL_Settings* sets )
    : maxiter(sets->Get_maxiter()),
      verbose(sets->Get_verbose())
  {}

  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename W::value_type reduction)
  {
    Dune::MatrixAdapter<M,V,W> opa(istl::raw(A));
    Preconditioner<M,V,W,1> prec(istl::raw(A), 3, 1.0);
    Solver<V> solver(opa, prec, reduction, maxiter, verbose);
    Dune::InverseOperatorResult stat;
    solver.apply(istl::raw(z), istl::raw(r), stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
    res.conv_rate  = stat.conv_rate;
  }

private:
  unsigned maxiter;
  int verbose;
};






/////////////////////////////////////////////////////////
///! Jacobi Solver
/////////////////////////////////////////////////////////

/**
 * @brief Interface for sequential loop solver with Jacobi preconditioner.
 */
class DuneISTL_SEQ_LOOP_Jac
  : public DuneISTL_SEQ_Base<Dune::SeqJac, Dune::LoopSolver>
{
public:
  /*! \brief make a linear solver object
    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_LOOP_Jac (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_Base<Dune::SeqJac, Dune::LoopSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_LOOP_Jac (DuneISTL_Settings* sets )
    : DuneISTL_SEQ_Base<Dune::SeqJac, Dune::LoopSolver>(sets)
  {}
};





/////////////////////////////////////////////////////////
///! ILU Solvers
/////////////////////////////////////////////////////////

template<template<typename> class Solver>
class DuneISTL_SEQ_ILU0
  :  public SequentialNorm, public LinearResultStorage
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_ILU0 (unsigned maxiter_=5000, int verbose_=1)
    : maxiter(maxiter_),
      verbose(verbose_)
   {}
  explicit DuneISTL_SEQ_ILU0 (DuneISTL_Settings* sets )
    : maxiter(sets->Get_maxiter()),
      verbose(sets->Get_verbose())
   {}
  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename W::value_type reduction)
  {
    Dune::MatrixAdapter<M,V,W> opa(istl::raw(A));
    Dune::SeqILU0<M,V,W> ilu0(istl::raw(A), 1.0);
    Solver<V> solver(opa, ilu0, reduction, maxiter, verbose);
    Dune::InverseOperatorResult stat;
    solver.apply(istl::raw(z), istl::raw(r), stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
    res.conv_rate  = stat.conv_rate;
   }
private:
  unsigned maxiter;
  int verbose;
};

template<template<typename> class Solver>
class DuneISTL_SEQ_ILUn
  :  public SequentialNorm, public LinearResultStorage
{
public:
  /*! \brief make a linear solver object
    \param[in] n The number of levels to be used.
    \param[in] w The relaxation factor.
    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  DuneISTL_SEQ_ILUn (int n=1, double w=1.0, unsigned maxiter_=5000, int verbose_=1)
    : n_(n),
      w_(w),
      maxiter(maxiter_),
      verbose(verbose_)
   {}
  DuneISTL_SEQ_ILUn (DuneISTL_Settings* sets)
    : n_(sets->Get_n()),
      w_(sets->Get_w()),
      maxiter(sets->Get_maxiter()),
      verbose(sets->Get_verbose())
   {}
  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename W::value_type reduction)
  {
    Dune::MatrixAdapter<M,V,W> opa(istl::raw(A));
    Dune::SeqILUn<M,V,W> ilun(istl::raw(A), n_, w_);
    Solver<V> solver(opa, ilun, reduction, maxiter, verbose);
    Dune::InverseOperatorResult stat;
    solver.apply(istl::raw(z), istl::raw(r), stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
    res.conv_rate  = stat.conv_rate;
   }
private:
  int n_;
  double w_;

  unsigned maxiter;
  int verbose;
};






/////////////////////////////////////
///! MINRES solvers
/////////////////////////////////////

/**
 * @brief Interface using a MINRes solver preconditioned by SSOR.
 */
class DuneISTL_SEQ_MINRES_SSOR
  : public DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::MINRESSolver>
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_MINRES_SSOR (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::MINRESSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_MINRES_SSOR (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::MINRESSolver>(sets)
  {}
};





/////////////////////////////////////
///! CG - conjugate gradient solvers
/////////////////////////////////////

/**
 * @brief Interface for conjugate gradient solver with Jacobi preconditioner.
 */
class DuneISTL_SEQ_CG_Jac
  : public DuneISTL_SEQ_Base<Dune::SeqJac, Dune::CGSolver>
{
public:
  /*! \brief make a linear solver object
    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_CG_Jac (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_Base<Dune::SeqJac, Dune::CGSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_CG_Jac (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_Base<Dune::SeqJac, Dune::CGSolver>(sets)
  {}
};

/**
 * @brief Interface for sequential conjugate gradient solver with SSOR preconditioner.
 */
class DuneISTL_SEQ_CG_SSOR
  : public DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::CGSolver>
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_CG_SSOR (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::CGSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_CG_SSOR (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::CGSolver>(sets)
  {}
};

/**
 * @brief Interface for sequential conjugate gradient solver with ILU0 preconditioner.
 */
class DuneISTL_SEQ_CG_ILU0
  : public DuneISTL_SEQ_ILU0<Dune::CGSolver>
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_CG_ILU0 (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_ILU0<Dune::CGSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_CG_ILU0 (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_ILU0<Dune::CGSolver>(sets)
  {}
};

//! \brief Sequential congute gradient solver with ILU0 preconditioner
class DuneISTL_SEQ_CG_ILUn
  : public DuneISTL_SEQ_ILUn<Dune::CGSolver>
{
public:
  /*! \brief make a linear solver object


    \param[in] n_ The number of levels to be used.
    \param[in] w_ The relaxation factor.
    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_CG_ILUn (int n_=1, double w_=1.0, unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_ILUn<Dune::CGSolver>(n_, w_, maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_CG_ILUn (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_ILUn<Dune::CGSolver>(sets)
  {}
};









/////////////////////////////////////
///! BiCGSTAB solvers
/////////////////////////////////////

/**
 * @brief Interface for sequential BiCGSTAB solver with Jacobi preconditioner.
 */
class DuneISTL_SEQ_BCGS_Jac
  : public DuneISTL_SEQ_Base<Dune::SeqJac, Dune::BiCGSTABSolver>
{
public:
  /*! \brief make a linear solver object
    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_BCGS_Jac (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_Base<Dune::SeqJac, Dune::BiCGSTABSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_BCGS_Jac (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_Base<Dune::SeqJac, Dune::BiCGSTABSolver>(sets)
  {}
};

/**
 * @brief Interface for sequential BiCGSTAB solver with SSOR preconditioner.
 */
class DuneISTL_SEQ_BCGS_SSOR
  : public DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::BiCGSTABSolver>
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_BCGS_SSOR (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::BiCGSTABSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_BCGS_SSOR (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_Base<Dune::SeqSSOR, Dune::BiCGSTABSolver>(sets)
  {}
};

 /**
 * @brief Interface for sequential BiCGSTAB solver with ILU0 preconditioner.
 */
class DuneISTL_SEQ_BCGS_ILU0
  : public DuneISTL_SEQ_ILU0<Dune::BiCGSTABSolver>
{
public:
  /*! \brief make a linear solver object

    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_BCGS_ILU0 (unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_ILU0<Dune::BiCGSTABSolver>(maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_BCGS_ILU0 (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_ILU0<Dune::BiCGSTABSolver>(sets)
  {}
};

//! \brief Sequential BiCGStab solver with ILU0 preconditioner
class DuneISTL_SEQ_BCGS_ILUn
  : public DuneISTL_SEQ_ILUn<Dune::BiCGSTABSolver>
{
public:
  /*! \brief make a linear solver object


    \param[in] n_ The number of levels to be used.
    \param[in] w_ The relaxation factor.
    \param[in] maxiter_ maximum number of iterations to do
    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_BCGS_ILUn (int n_=1, double w_=1.0, unsigned maxiter_=5000, int verbose_=1)
    : DuneISTL_SEQ_ILUn<Dune::BiCGSTABSolver>(n_, w_, maxiter_, verbose_)
  {}
  explicit DuneISTL_SEQ_BCGS_ILUn (DuneISTL_Settings* sets)
    : DuneISTL_SEQ_ILUn<Dune::BiCGSTABSolver>(sets)
  {}
};







/////////////////////////////////////
///! SUPERLU Solver
/////////////////////////////////////

#if HAVE_SUPERLU
/**
 * @brief Interface for SuperLU as a direct solver.
 */
class DuneISTL_SEQ_SuperLU
  : public SequentialNorm, public LinearResultStorage
{
public:
  /*! \brief make a linear solver object

    \param[in] verbose_ print messages if true
  */
  explicit DuneISTL_SEQ_SuperLU (int verbose_=1)
    : verbose(verbose_)
  {}

  /*! \brief make a linear solver object

    \param[in] maxiter Maximum number of allowed steps (ignored)
    \param[in] verbose_ print messages if true
  */
  DuneISTL_SEQ_SuperLU (int maxiter, int verbose_)
    : verbose(verbose_)
  {}

  explicit DuneISTL_SEQ_SuperLU (DuneISTL_Settings* sets )
    : verbose(sets->Get_verbose())
  {}

  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  template<class M, class V, class W>
  void apply(M& A, V& z, W& r, typename W::value_type reduction)
  {
    typedef typename M::Container ISTLM;
    Dune::SuperLU<ISTLM> solver(istl::raw(A), verbose);
    Dune::InverseOperatorResult stat;
    solver.apply(istl::raw(z), istl::raw(r), stat);
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
    res.conv_rate  = stat.conv_rate;
  }

private:
  int verbose;
};
#endif // HAVE_SUPERLU


/////////////////////////////////////
///! AMG Solver basic
/////////////////////////////////////

/**
 * @brief Class providing some statistics of the AMG solver.
 *
 */
struct ISTLAMGStatistics
{
  /**
   * @brief The needed for computing the parallel information and
   * for adapting the linear system.
   */
  double tprepare;
  /** @brief the number of levels in the AMG hierarchy. */
  int levels;
  /** @brief The time spent in solving the system (without building the hierarchy. */
  double tsolve;
  /** @brief The time needed for building the AMG hierarchy (coarsening). */
  double tsetup;
  /** @brief The number of iterations performed until convergence was reached. */
  int iterations;
  /** @brief True if a direct solver was used on the coarset level. */
  bool directCoarseLevelSolver;
};

template<uint32_t dim, class M, class V,
         template<class,class,class,int> class Preconditioner,
         template<class> class Solver,
         bool skipBlocksizeCheck = false>
class DuneISTL_SEQ_AMG : public LinearResultStorage
{
  typedef M MatrixType;
  typedef V VectorType;
  typedef Preconditioner<MatrixType,VectorType,VectorType,1> Smoother;
  typedef Dune::MatrixAdapter<MatrixType,VectorType,VectorType> Operator;
  typedef typename Dune::Amg::SmootherTraits<Smoother>::Arguments SmootherArgs;
  typedef Dune::Amg::AMG<Operator,VectorType,Smoother> AMG;
  typedef Dune::Amg::Parameters Parameters;

public:
  DuneISTL_SEQ_AMG(unsigned maxiter_=5000, int verbose_=1,
                   bool reuse_=false, bool usesuperlu_=true)
    : maxiter(maxiter_),
      params(15,2000),
      verbose(verbose_),
      reuse(reuse_),
      firstapply(true),
      usesuperlu(usesuperlu_)
  {
    params.setDefaultValuesIsotropic(dim);
    params.setDebugLevel(verbose_);
#if !HAVE_SUPERLU
    if (usesuperlu == true)
      {
        std::cout << "WARNING: You are using AMG without SuperLU!"
                  << " Please consider installing SuperLU,"
                  << " or set the usesuperlu flag to false"
                  << " to suppress this warning." << std::endl;
      }
#endif
  }

  DuneISTL_SEQ_AMG(DuneISTL_Settings* sets)
    : maxiter(sets->Get_maxiter()),
      params(sets->Get_amg_params()),
      verbose(sets->Get_verbose()),
      reuse(sets->Get_reuse()),
      usesuperlu(sets->Get_usesuperlu()),
      firstapply(true)
  {
#if !HAVE_SUPERLU
    if (usesuperlu == true)
      {
        std::cout << "WARNING: You are using AMG without SuperLU!"
                  << " Please consider installing SuperLU,"
                  << " or set the usesuperlu flag to false"
                  << " to suppress this warning." << std::endl;
      }
#endif
  }

   /*! \brief set AMG parameters

    \param[in] params_ a parameter object of Type Dune::Amg::Parameters
  */
  void setparams(Parameters params_)
  {
    params = params_;
  }

  /*! \brief compute global norm of a vector

    \param[in] v the given vector
  */
  typename V::value_type norm (const V& v) const
  {
    return istl::raw(v).two_norm();
  }

  /*! \brief solve the given linear system

    \param[in] A the given matrix
    \param[out] z the solution vector to be computed
    \param[in] r right hand side
    \param[in] reduction to be achieved
  */
  void apply(M& A, V& z, V& r, typename V::value_type reduction)
  {
    Dune::Timer watch;
    MatrixType& mat=istl::raw(A);
    typedef Dune::Amg::CoarsenCriterion<Dune::Amg::SymmetricCriterion<MatrixType,
      Dune::Amg::FirstDiagonal> > Criterion;
    SmootherArgs smootherArgs;
    smootherArgs.iterations = 1;
    smootherArgs.relaxationFactor = 1;

    Criterion criterion(params);
    Operator oop(mat);
    //only construct a new AMG if the matrix changes
    if (reuse==false || firstapply==true){
      amg.reset(new AMG(oop, criterion, smootherArgs));
      firstapply = false;
      stats.tsetup = watch.elapsed();
      stats.levels = amg->maxlevels();
      stats.directCoarseLevelSolver=amg->usesDirectCoarseLevelSolver();
    }
    watch.reset();
    Dune::InverseOperatorResult stat;

    Solver<VectorType> solver(oop,*amg,reduction,maxiter,verbose);
    solver.apply(istl::raw(z),istl::raw(r),stat);
    stats.tsolve= watch.elapsed();
    res.converged  = stat.converged;
    res.iterations = stat.iterations;
    res.elapsed    = stat.elapsed;
    res.reduction  = stat.reduction;
    res.conv_rate  = stat.conv_rate;
  }


  /**
   * @brief Get statistics of the AMG solver (no of levels, timings).
   * @return statistis of the AMG solver.
   */
  const ISTLAMGStatistics& statistics() const
  {
    return stats;
  }

private:
  unsigned maxiter;
  Parameters params;
  int verbose;
  bool reuse;
  bool firstapply;
  bool usesuperlu;
  Dune::shared_ptr<AMG> amg;
  ISTLAMGStatistics stats;
};






/////////////////////////////////////////////////////////////////////
///! Iterative Solver with AMG Preconditioner smoothed by SSOR or SOR
/////////////////////////////////////////////////////////////////////

/**
 * @brief Sequential Loop solver preconditioned with AMG smoothed by SSOR
 * @tparam GO The type of the grid operator
 * (or the fakeGOTraits class for the old grid operator space).
 */
template<uint32_t dim, class M, class V>
class DuneISTL_SEQ_LS_AMG_SSOR
  : public DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::LoopSolver>
{

public:
  /**
   * @brief Constructor
   * @param maxiter_ The maximum number of iterations allowed.
   * @param verbose_ The verbosity level to use.
   * @param reuse_ Set true, if the Matrix to be used is always identical
   * (AMG aggregation is then only performed once).
   * @param usesuperlu_ Set false, to suppress the no SuperLU warning
   */
  DuneISTL_SEQ_LS_AMG_SSOR(unsigned maxiter_=5000, int verbose_=1,
                           bool reuse_=false, bool usesuperlu_=true)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::LoopSolver>
      (maxiter_, verbose_, reuse_, usesuperlu_)
  {}
  DuneISTL_SEQ_LS_AMG_SSOR(DuneISTL_Settings* sets)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::LoopSolver>(sets)
  {}
};

/**
 * @brief Sequential Loop solver preconditioned with AMG smoothed by SOR
 * @tparam GO The type of the grid operator
 * (or the fakeGOTraits class for the old grid operator space).
 */
template<uint32_t dim, class M, class V>
class DuneISTL_SEQ_LS_AMG_SOR
  : public DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSOR, Dune::LoopSolver>
{

public:
  /**
   * @brief Constructor
   * @param maxiter_ The maximum number of iterations allowed.
   * @param verbose_ The verbosity level to use.
   * @param reuse_ Set true, if the Matrix to be used is always identical
   * (AMG aggregation is then only performed once).
   * @param usesuperlu_ Set false, to suppress the no SuperLU warning
   */
  DuneISTL_SEQ_LS_AMG_SOR(unsigned maxiter_=5000, int verbose_=1,
                          bool reuse_=false, bool usesuperlu_=true)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSOR, Dune::LoopSolver>
      (maxiter_, verbose_, reuse_, usesuperlu_)
  {}
  DuneISTL_SEQ_LS_AMG_SOR(DuneISTL_Settings* sets)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSOR, Dune::LoopSolver>(sets)
  {}
};




///////////////////////////////////////////////////////////////
///! CG Solver with AMG Preconditioner smoothed by SSOR
///////////////////////////////////////////////////////////////


/**
 * @brief Sequential conjugate gradient solver preconditioned with AMG smoothed by SSOR
 * @tparam GO The type of the grid operator
 * (or the fakeGOTraits class for the old grid operator space).
 */
template<uint32_t dim, class M, class V>
class DuneISTL_SEQ_CG_AMG_SSOR
  : public DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::CGSolver>
{

public:
  /**
   * @brief Constructor
   * @param maxiter_ The maximum number of iterations allowed.
   * @param verbose_ The verbosity level to use.
   * @param reuse_ Set true, if the Matrix to be used is always identical
   * (AMG aggregation is then only performed once).
   * @param usesuperlu_ Set false, to suppress the no SuperLU warning
   */
  DuneISTL_SEQ_CG_AMG_SSOR(unsigned maxiter_=5000, int verbose_=1,
                           bool reuse_=false, bool usesuperlu_=true)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::CGSolver>
      (maxiter_, verbose_, reuse_, usesuperlu_)
  {}
  DuneISTL_SEQ_CG_AMG_SSOR(DuneISTL_Settings* sets)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::CGSolver>(sets)
  {}
};


///////////////////////////////////////////////////////////////
///! BiCGSTAB Solver with AMG Preconditioner smoothed by SSOR
///////////////////////////////////////////////////////////////


/**
 * @brief Sequential BiCGStab solver preconditioned with AMG smoothed by SSOR
 * @tparam GO The type of the grid operator
 * (or the fakeGOTraits class for the old grid operator space).
 */
template<uint32_t dim, class M, class V>
class DuneISTL_SEQ_BCGS_AMG_SSOR
  : public DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::BiCGSTABSolver>
{

public:
  /**
   * @brief Constructor
   * @param maxiter_ The maximum number of iterations allowed.
   * @param verbose_ The verbosity level to use.
   * @param reuse_ Set true, if the Matrix to be used is always identical
   * (AMG aggregation is then only performed once).
   * @param usesuperlu_ Set false, to suppress the no SuperLU warning
   */
  DuneISTL_SEQ_BCGS_AMG_SSOR(unsigned maxiter_=5000, int verbose_=1,
                             bool reuse_=false, bool usesuperlu_=true)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::BiCGSTABSolver>
      (maxiter_, verbose_, reuse_, usesuperlu_)
  {}
  DuneISTL_SEQ_BCGS_AMG_SSOR(DuneISTL_Settings* sets)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSSOR, Dune::BiCGSTABSolver>(sets)
  {}
};

/**
 * @brief Sequential BiCGSTAB solver preconditioned with AMG smoothed by SOR
 * @tparam GO The type of the grid operator
 * (or the fakeGOTraits class for the old grid operator space).
 */
template<uint32_t dim, class M, class V>
class DuneISTL_SEQ_BCGS_AMG_SOR
  : public DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSOR, Dune::BiCGSTABSolver>
{

public:
  /**
   * @brief Constructor
   * @param maxiter_ The maximum number of iterations allowed.
   * @param verbose_ The verbosity level to use.
   * @param reuse_ Set true, if the Matrix to be used is always identical
   * (AMG aggregation is then only performed once).
   * @param usesuperlu_ Set false, to suppress the no SuperLU warning
   */
  DuneISTL_SEQ_BCGS_AMG_SOR(unsigned maxiter_=5000, int verbose_=1,
                            bool reuse_=false, bool usesuperlu_=true)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSOR, Dune::BiCGSTABSolver>
      (maxiter_, verbose_, reuse_, usesuperlu_)
  {}
  DuneISTL_SEQ_BCGS_AMG_SOR(DuneISTL_Settings* sets)
    : DuneISTL_SEQ_AMG<dim,M,V, Dune::SeqSOR, Dune::BiCGSTABSolver>(sets)
  {}
};


} // namespace csmp

#endif // end CSMP_WITH_DUNE_ISTL

#endif
