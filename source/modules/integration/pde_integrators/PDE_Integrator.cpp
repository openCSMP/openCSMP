#include "PDE_Integrator.h"
#if defined(_OPENMP )
#include "omp.h"
#include "FiniteElementManager.h"
#endif

using namespace std;

namespace csmp {


/** Default Constructor

The default constructor sets the computational domain equivalent to the
Model.
The default Solver is SAMG, without license, you only have access
to a single level solver.
The default computation is a steady-state computation. Thus, if you
want to carry out a transient calculation you have to use the method
Transient().

@attention SIMPLICIAL_COMPLEX is used here because DOMAIN caused a clash
with DOMAIN defined in <cmath>
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::PDE_Integrator()
 :
#ifdef CSMP_WITH_SAMG_SOLVER
   solver_(new SAMG_Solver()),
#else
   /// add extra functionality for alternative solver if needed
   solver_(new CSMP_DEFAULT_LINEAR_SOLVER()),
#endif
   dim2_(dim*dim),
   dof_per_node_(0),
   setup_established_(false),
   retain_matrix_(false),
   newed_Solver_object(true),
   time_increment_(0.),
   scale_factor_(1.),
   verbose_(true)
{
   target_.nodes = target_.elements = 0U;
#if defined(_OPENMP )
   this->femgrs_.resize(omp_get_max_threads());
   this->thread_G_.resize(omp_get_max_threads());
   this->thread_rh_.resize(omp_get_max_threads());
   this->thread_lhs_operators_.resize(omp_get_max_threads());
   this->thread_rhs_operators_.resize(omp_get_max_threads());
#pragma omp parallel
    {
        size_t tid = omp_get_thread_num();
        this->femgrs_[tid].InitializeElements(dim,1,true);
    }
#endif
}



/**
    Use this constructor wherever possible.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::PDE_Integrator( Solver& solver )
 : solver_(&solver),
   dim2_(dim*dim),
   dof_per_node_(0),
   setup_established_(false),
   retain_matrix_(false),
   newed_Solver_object(false),
   time_increment_(0.),
   scale_factor_(1.),
   verbose_(true)
{
  target_.nodes = target_.elements = 0U;
#if defined(_OPENMP )
    femgrs_.resize(omp_get_max_threads());
    this->thread_G_.resize(omp_get_max_threads());
    this->thread_rh_.resize(omp_get_max_threads());
    this->thread_lhs_operators_.resize(omp_get_max_threads());
    this->thread_rhs_operators_.resize(omp_get_max_threads());
#pragma omp parallel
    {
        size_t tid = omp_get_thread_num();
        this->femgrs_[tid].InitializeElements(dim,1,true);
    }
#endif
}




/** custom constructor.

  @attention Potential memory leak! - if you write     integrator(new Solver())
 
  The constructor will not manage the supplied pointer, i.e. delete the supplied object instance.
 
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::PDE_Integrator( Solver* solver )
 : solver_(solver),
   dim2_(dim*dim),
   dof_per_node_(0),
   setup_established_(false),
   retain_matrix_(false),
   newed_Solver_object(false),
   time_increment_(0.),
   scale_factor_(1.),
   verbose_(true)
{
  assert( solver != NULL );
  target_.nodes = target_.elements = 0U;
#if defined(_OPENMP )
    femgrs_.resize(omp_get_max_threads());
    this->thread_G_.resize(omp_get_max_threads());
    this->thread_rh_.resize(omp_get_max_threads());
    this->thread_lhs_operators_.resize(omp_get_max_threads());
    this->thread_rhs_operators_.resize(omp_get_max_threads());
#pragma omp parallel
    {
        size_t tid = omp_get_thread_num();
        this->femgrs_[tid].InitializeElements(dim,1,true);
    }
#endif
}








template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::~PDE_Integrator()
 {
    if ( newed_Solver_object ) delete solver_;
#if defined(_OPENMP )
    for ( typename map<string,MathOperatorLHS<dim>*>::iterator
          it_lhs=lhs_operators_.begin(); it_lhs!=lhs_operators_.end(); it_lhs++ )
        for (size_t tid = 0 ; tid < omp_get_num_threads() ;tid++)
            if (this->thread_lhs_operators_[tid].find(it_lhs->second->Name()) != this->thread_lhs_operators_[tid].end())
                delete(this->thread_lhs_operators_[tid].at(it_lhs->second->Name()));

    for ( typename map<string,MathOperatorRHS<dim>*>::iterator
          it_rhs=rhs_operators_.begin(); it_rhs!=rhs_operators_.end(); it_rhs++ )
        for (size_t tid = 0 ; tid < omp_get_num_threads() ;tid++)
            if (this->thread_rhs_operators_[tid].find(it_rhs->second->Name()) != this->thread_rhs_operators_[tid].end())
                delete(this->thread_rhs_operators_[tid].at(it_rhs->second->Name()));
#endif
 }



/** 
@attention A new Solver object is created here and it is of type LUdcmp_Solver, hence not honoring provided instance
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::PDE_Integrator( const PDE_Integrator<dim,SIMPLICIAL_COMPLEX>& a )
 :
   #ifdef CSMP_WITH_SAMG_SOLVER
   solver_(new SAMG_Solver()),
   #else
   /// add extra functionality for alternative solver if needed
   solver_(new CSMP_DEFAULT_LINEAR_SOLVER()),
   #endif
   rh_(a.rh_),
   x_(a.x_),
   G_(a.G_),
   lhs_operators_(a.lhs_operators_),
   rhs_operators_(a.rhs_operators_),
   basic_operands_(a.basic_operands_),
   test_operands_(a.test_operands_),
   postpro_operators_(a.postpro_operators_),
   dim2_(dim*dim),
   dof_per_node_(a.dof_per_node_),
   setup_established_(a.setup_established_),
   retain_matrix_(a.retain_matrix_),
   newed_Solver_object(true),
   scale_factor_(1.),
   time_increment_(a.time_increment_),
   target_(a.target_),
   verbose_(a.verbose_)
{
   assert( solver_ != NULL );
   cout <<"\nPDE_Integrator: copy constructor: ";
   cout <<"New algorithm uses new instance of solver."<< endl;
#if defined(_OPENMP )
    femgrs_.resize(omp_get_max_threads());
    this->thread_G_.resize(omp_get_max_threads());
    this->thread_rh_.resize(omp_get_max_threads());
    this->thread_lhs_operators_.resize(omp_get_max_threads());
    this->thread_rhs_operators_.resize(omp_get_max_threads());
#pragma omp parallel
    {
        size_t tid = omp_get_thread_num();
        this->femgrs_[tid].InitializeElements(dim,1,true);
    }
#endif
}



/** 
@attention A new Solver object is created here and it is of type LUdcmp_Solver, hence not honoring provided instance

@todo SKM: fix so that the copy construction also involves the Solver
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
PDE_Integrator<dim,SIMPLICIAL_COMPLEX>& PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::operator=( const PDE_Integrator<dim,SIMPLICIAL_COMPLEX>& a )
{
  if ( &a != this ) {
      G_                 = a.G_;
      rh_                = a.rh_;
      x_                 = a.x_;

      if ( a.solver_ != NULL ) {
           assert( solver_ != NULL or solver_ == 0 );
           if ( newed_Solver_object ) {
               delete solver_;
               /// SKM fix:  this is not a clean solution. One should bring over solver from other integrator instance
               #ifdef CSMP_WITH_SAMG_SOLVER
               solver_ = new SAMG_Solver();
               #else
               /// add extra functionality for alternative solver if needed
               solver_ = new CSMP_DEFAULT_LINEAR_SOLVER();
               #endif
            }
           else solver_ = a.solver_;
        }

      lhs_operators_     = a.lhs_operators_;
      rhs_operators_     = a.rhs_operators_;
      postpro_operators_ = a.postpro_operators_;
      basic_operands_    = a.basic_operands_;
      test_operands_     = a.test_operands_;
      time_increment_    = a.time_increment_;
      retain_matrix_     = a.retain_matrix_;
      setup_established_ = a.setup_established_;
      target_            = a.target_;
      dof_per_node_      = a.dof_per_node_;
      scale_factor_      = a.scale_factor_;
      verbose_           = a.verbose_;

#if defined(_OPENMP )
      this->femgrs_=a.femgrs_;
      this->thread_G_=a.thread_G_;
      this->thread_rh_=a.thread_rh_;
      this->thread_lhs_operators_=a.thread_lhs_operators_;
      this->thread_rhs_operators_=a.thread_rhs_operators_;
#endif
    }
  return *this;
}



template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::RetainGlobalSolutionMatrix( bool retain ) 
 { retain_matrix_=retain; }




/**
     Deletes the current solver replacing it with the supplied one.
     Does not assume responsibility for the supplield solver 
     which must be deleted elsewhere if created on the heap.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::SetSolver( Solver* new_solver ) {
  assert( new_solver != NULL );
  if ( solver_ != new_solver ) {
       delete solver_;
       solver_ = new_solver;
    }
  newed_Solver_object = false;
}



template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
Solver* PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::GetSolver() const {
  return solver_;
}



template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AdjustSolverSettings()
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    csmp_error.notice( WARNING, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AdjustSolverSettings",
                      "this call to the method did nothing, the method is defined only in the subclasses.");

 } // end AdjustSolverSettings




/** Transient PDE_Integrator

The accumulation process for a transient PDE_Integrator is different from a
steady-state PDE_Integrator in that initial conditions must be assigned and
matrices must be multiplied with the time_increment. This method
sets the PDE_Integrator state to transient or retrieves its state.

@section application Application

Every PDE_Integrator which you define in your program uses storage for G, x, and
rhs. You should therefore always limit the number of A. to those essential.
Use this method if you want to re-use a previously defined steady-state
PDE_Integrator.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
bool  PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::Transient() const
 { return !(time_increment_ < numeric_limits<double64>::epsilon()); }




/** Sets the time-icrement in a transient calculation.

If one forget to set
the PDE_Integrator to transient, this is done as well.

@section application Application

Set the time-increment of an PDE_Integrator before you Apply() it to the
Model or target Region objects.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void   PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::TimeIncrement( double64 dt )
 {
    time_increment_ = dt;
 }




/** Prints out the sparse solution matrix and the righthand vector

...which were accumulated by the PDE_Integrator. If you execute this method after passing
the PDE_Integrator to the Model, the matrix will have been modified by the
Solver object.

@section application  Application

To test the accumulation process by visual examination of the matrices,
you must call it directly after executing Accumulate(), see below.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void  PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::OutputGlobals( int32 precision )
 {
   cout <<"\nGlobal solution matrix: "<< G_.Rows() <<" x "<< G_.Cols() << endl;

   cout.setf(ios::scientific);

   for ( size_t i=0; i<G_.Rows(); i++ )
      {
         for ( size_t j=0; j<G_.Cols(); j++ )
           {
              cout.precision(precision);
              if ( G_.At(i,j) >= 0 ) cout <<" ";
              cout <<  G_.At(i,j) <<" ";
           }
         cout << endl;
      }

   cout <<"\n\nGlobal righthand vector of length: "<< rh_.size() << endl;
   for ( size_t i=0U; i<rh_.size(); i++ )
     {
        cout.precision(precision);
        if ( rh_[i] >= 0. ) cout <<" ";
        cout << rh_[i] <<" ";
     }
   cout << endl;

   cout <<"\n\nGlobal solution vector of length: "<< x_.size() << endl;
   for ( size_t i=0U; i<x_.size(); i++ )
     {
        cout.precision(precision);
        if ( x_[i] >= 0 ) cout <<" ";
        cout << x_[i] <<" ";
     }
   cout << endl;

   cout.unsetf(ios::scientific);

 } // end OutputGlobals







/**

The PDE_Integrator object maintains a list of pointers to MathOperator instances
which represent parts of the PDE which is being modeled by the PDE_Integrator.
Add() allows you to associate such operators_ (like dove^2) with the
PDE_Integrator.

@param op PDE operator usually representing the integral form of an equation
term that shall be added to the integrator.

Pointer to PDE operator objects (=Operands) which you defined
earlier in your program. These are always subclasses derived
from the MathOperatorLHSptr or MathOperatorRHSptr base classes.

@section implementation Implementation

The MathOperator pointers are stored in the lhs_operators_ and rhs_operators_
maps. Their functionality is accessed in the accumulation process of the
solution matrix G and the righthand vector rh. The maps prevent, that
multiple instances of the same MathOperator object with the same
name are added.

@section application Application

Each PDE_Integrator needs at least one left and one righthand math operator.
Define these before you pass the PDE_Integrator to the Model.

*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::Add( MathOperatorLHS<dim>* op )
 {
    // The name for the algorithm is combined out of its operands
    lhs_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
#if defined(_OPENMP )
    for (size_t tid = 0 ; tid < omp_get_max_threads() ;tid++)
        thread_lhs_operators_[tid][op->Name()]=op->clone();
#endif
 }

template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::Add( MathOperatorRHS<dim>* op )
 {
    rhs_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
#if defined(_OPENMP )
    for (size_t tid = 0 ; tid < omp_get_max_threads() ;tid++)
        thread_rhs_operators_[tid][op->Name()]=op->clone();
#endif
 }


template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AddPostProcess( MathOperatorLHS<dim>* op )
 {
    postpro_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;

#if defined(_OPENMP )

#endif
 }




/**

Lists all MathOperator objects which have been defined for the PDE_Integrator.
These are output by name in alphabetical order. For the output, the
MathOperator interface method Out() is used. Thus, you influence the
information that is output when you define your own PDE operator
subclasses.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::ListMathOperatorsLHS() const
 {
    typename map<string,MathOperatorLHS<dim>*>::const_iterator it;

    for ( it=lhs_operators_.begin(); it!=lhs_operators_.end(); it++ )
      {
         cout << (*it).first <<":  ";
         (*it).second->Out();
      }
 }


template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::ListMathOperatorsRHS() const
 {
    typename map<string,MathOperatorRHS<dim>*>::const_iterator it;

    for ( it=rhs_operators_.begin(); it!=rhs_operators_.end(); it++ )
      {
         cout << (*it).first <<":  ";
         (*it).second->Out();
      }
 }






/**

@section application Application

Use this method to extract the solution vector from the algorithm, for
instance to use it as initial guess in another time step. (Use method
FirstGuess)
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void  PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::SolutionVector( vector<double64>& sol ) const {
  sol.resize(x_.size());
  vector<double64>( sol ).swap( sol );
  copy(x_.begin(), x_.end(), sol.begin());
}




/**

Registers vector guess as first guess in the solution procedure for the
algebraic multigrid solver.

@section application  Application

Use this method to supply a first guess solution different from the previous
solution obtained with this algorithm.

NOTE: Make sure that you actually use an algebraic multigrid solver object and that the
parameter ifirst in its SAMG_Settings object is set to 0.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void  PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::FirstGuess( const vector<double64>& guess ) {
  assert(guess.size() == x_.size());
  copy(guess.begin(), guess.end(), x_.begin());
}





/**

After an PDE_Integrator has been applied, the solution matrix G and the righthand
vector rh contain processed accumulation data which need to be removed before
a new accumulation can take place. This is done automatically by the
PDE_Integrator, unless it has been set to retain the solution matrix.

In addition, you may want to re-use the PDE_Integrator and define
new PDE operators, a new computational domain and/or computational
method. In this case, Reset() allows to restore the PDE_Integrator's default state
and zero's the elements of its matrix and vector storage.

@section application Application

Typically, you will reset an PDE_Integrator if you can just re-use it for
another task. This is always recommended, since this permits to make use
of the already initialized solution matrix, because the memory of it
(owing to the Meschach implementation) is not
de-allocated before the programme terminates.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void  PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::Reset( bool delete_math_operators )
 {
    if ( delete_math_operators ) {
         lhs_operators_.erase(  lhs_operators_.begin(),  lhs_operators_.end() );
         rhs_operators_.erase(  rhs_operators_.begin(),  rhs_operators_.end() );
         basic_operands_.erase( basic_operands_.begin(), basic_operands_.end() );
         test_operands_.erase(  test_operands_.begin(),  test_operands_.end() );
         postpro_operators_.erase( postpro_operators_.begin(), postpro_operators_.end() );
         setup_established_ = false;
      }

    // restoring defaults
    time_increment_ = 0.;
    retain_matrix_  = false;
    target_.nodes   = target_.elements = 0U;

} // end Reset






/**

Calls the Solver interface SolveMatrixEquation() and hands G, rh, x to the
Solver which returns the solution into the vector x.

@section application Application

The method Solve() is called internally, when the PDE_Integrator object is
passed to the Model.

The solver reports a range of parameters in the solution process. These
are described in detail in the CSMP User's guide. Important is that
the final residue of the solution (which could be compared to a
signal-to-noise ratio) is several orders of magnitude lower than the
initial residue. Otherwise the obtained solution is useless.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::Solve()
 {
    solver_->Solve( G_, rh_, x_, dof_per_node_ );
 } // end Solve








/**

EstablishMatrixSetup() involves the following steps:

1. The basic sizes of the solution matrix, solution vector and righthand
   vector are determined.

2. Each MathOperator has a basic, and a test variable operand. Basic operands
   are checked for their type and compared among operations to identify what
   size the final solution matrix G should have.

3. Offset are assigned to individual MathOperators to indicate to what
   positions in the solution matrix they must accumulate to. This applies
   when multiple dependent variables are being solved for in the solution
   process.

4. The solution matrix G, the solution vector x and the righthand vector
   rh are resized according to the operator specifications and the
   application target (Region or Model).

EstablishMatrixSetup() takes information from the mesh and property
managers and may query the property database for the type of variables.
In the case of a group-restricted computation, the group interface is
accessed via a pointer (gptr).

@section application Application

EstablishMatrixSetup() is applied internally when the PDE_Integrator is
passed to the Model.

@section messages Messages
EstablishMatrixSetup() can terminate emitting messages when:

@code
    computational domain not specified"
@endcode

If the target group could not be found.

@code
PDE_Integrator cannot add new variable to database
@endcode

When you mispelled a variable name or if you tried to use a non-existant
variable in your computation. .

@code
righthand basic Operand's offset could not be resolved...
@endcode

Again, when a physical variable is unknown to the PropertyDatabase object
such that it cannot resolve its csmp::Index.

The following Error won't terminate PDE_Integrator application:

@code
There is no test operand for basic operand: x csmp::Index: idx
@endcode

Here, the finite-element interpolation functions do not model the
dependent variable.

Warning messages that can arise during the execution of
EstablishMatrixSetup() are the following:

@code
PDE_Integrator was used before, using same specifications again...
@endcode

You forgot to call Reset() before re-using an existing PDE_Integrator. You will
now accumulate into the inverted solution matrix.

@code
lefthand basic operand not found
@endcode

Each PDE_Integrator needs a lefthand basic operand in order to define
a PDE problem. Define one before applying the PDE_Integrator.

@code
There is no test operand for basic operand:
@endcode

Your element interpolation functions do not model your
dependent variable.

@code
Basic operand is an Element variable. Is this intended ?
@endcode

Typically, the dependent variable in a finite-element computation, e.g.,
the basic operand is placed on the nodes. This lies in the very nature
of the finite-element method. If you are doing a finite-volume or
other computation, just ignore this warning.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup( const SIMPLICIAL_COMPLEX<dim>& gref )
{

    // -------------------------------------------------------------------
    // 0. If the algorithm is just re-used, (and has not been reset by
    //    the user, the righthand vector and
    //    the solution vector are zeroed and nothing else is done.
    // -------------------------------------------------------------------
    // * Change: x is not reset anymore
    if ( setup_established_ &&
         ( !basic_operands_.empty() && !test_operands_.empty() ) )
    {
        if ( rh_.size() > 0 ) fill( rh_.begin(), rh_.end(), 0. );
#if defined(_OPENMP )
        for (size_t tid = 0 ; tid < omp_get_max_threads() ; tid++){
           if ( thread_rh_[tid].size() > 0 ) fill( thread_rh_[tid].begin(), thread_rh_[tid].end(), 0. );
        }
#endif
        if ( G_.Rows()  > 0 && !retain_matrix_ ) {
            G_.Erase();
            G_.Resize( rh_.size() );
#if defined(_OPENMP )
            for (size_t tid = 0 ; tid < omp_get_max_threads() ; tid++){
                this->thread_G_[tid].Erase();
                this->thread_G_[tid].Resize(rh_.size());
            }
#endif
        }
    }


   // -------------------------------------------
   // 1. determine basic sizes for G, x, rh
   // -------------------------------------------
    dof_per_node_    = 0U;
    target_.nodes    = gref.Nodes();
    target_.elements = gref.Elements();


   // ------------------------------------------------------
   // 2. if A.-setup for first time or if rebuild is necessary:
   //    checking basic operands and comparing MathOperators
   // ------------------------------------------------------
   const Index  unspecified;

   // lefthand MathOperators first
   // ----------------------------
   for ( typename map<string,MathOperatorLHS<dim>*>::iterator
         lhs_it=lhs_operators_.begin(); lhs_it!=lhs_operators_.end(); lhs_it++ )
     {
        // making list of unique basic operands
        Index pkey = (*lhs_it).second->BasicOperandKey();
        if (verbose_)
            cout <<"\nFor: '"<< (*lhs_it).first <<"' PDE operator is added to lefthand term list."<< endl;

        // checking whether the intended variables exist in the database
        if ( unspecified != pkey )
            basic_operands_[(*lhs_it).second->BasicOperand()] = 0U;
        else
            throw csmp::Exception( WARNING,
                                   "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
                                   "lefthand basic operand not found");
        // test function operands are picked up when the righthandside is accumulated
        // since they must also be present in there
     }

   // righthand MathOperators
   // -----------------------
   for ( typename map<string,MathOperatorRHS<dim>*>::iterator
         rhs_it=rhs_operators_.begin(); rhs_it!=rhs_operators_.end(); rhs_it++ )
     {
        // making a list of unique test operands
        Index pkey = (*rhs_it).second->TestOperandKey();
        if (verbose_)
            cout <<"\nFor: '"<< (*rhs_it).first <<"' PDE operator is added to righthand term list."<< endl;

        if ( unspecified != pkey )
            test_operands_[(*rhs_it).second->TestOperand()] = 0U;
        else
            throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
                                         "righthand test operand not found");
     }

   // ---------------------------------------------------------------------------------------
   // 3. Testing: for each righthand operand there must be a basic or test operand on the LHS
   // ---------------------------------------------------------------------------------------
   for ( typename map<string,MathOperatorLHS<dim>*>::iterator
         lhs_it=lhs_operators_.begin(); lhs_it!=lhs_operators_.end(); lhs_it++ )
     if ( basic_operands_.find( ((*lhs_it).second->BasicOperand()) ) == basic_operands_.end() &&
          test_operands_.find( ((*lhs_it).second->TestOperand()) ) == test_operands_.end() )
       {
          cout <<"\nPDE_Integrator<"<<  dim <<">::EstablishMatrixSetup: ";
          cout <<"There is no lefthand operand corresponding to righthand operand. ";
          cout <<"\nThe system of equations is undefined. ";
          cout <<"\nCreate corresponding LHS basic or test Operand for: ";
          cout << (*lhs_it).first << endl;
          throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
                                         "lefthand basic or test operand missing");
       }

   // ----------------------------------------------------------------
   //   4. Offsets are assigned to test Operands
   //      indicating positions, i,j in solution matrix.
   // ----------------------------------------------------------------
   size_t  offset(0U);
   for ( operandsIterator
         iter=test_operands_.begin(); iter!=test_operands_.end(); iter++ )
     {
        // calculating the degrees of freedom per node
        if      ( (*iter).first.key.type == SCALAR )
          dof_per_node_ += 1U;
        else if ( (*iter).first.key.type == VECTOR )
          dof_per_node_ += dim;
        else if ( (*iter).first.key.type == TENSOR )
          dof_per_node_ += dim2_;
        else if ( (*iter).first.key.type == ARRAY )
          dof_per_node_ += (*iter).first.key.dataDepth;
        else if ( (*iter).first.key.type == FLAGGEDARRAY )
          dof_per_node_ += (*iter).first.key.dataDepth;

        // the matrix/righthand offsets are stored with the operands
        // offset starts out as zero.
        switch ( (*iter).first.key.place )
         {
            case NODE:
              (*iter).second = offset;
              switch ( (*iter).first.key.type ) {
                   case SCALAR: offset += target_.nodes;
                     break;
                   case VECTOR: offset += target_.nodes * dim;
                     break;
                   case TENSOR: offset += target_.nodes * dim2_;
                     break;
                   case ARRAY:  offset += target_.nodes * (*iter).first.key.dataDepth;
                     break;
                   case FLAGGEDARRAY:  offset += target_.nodes * (*iter).first.key.dataDepth;
                }
              break;
            case ELEMENT:
              cout <<"\nPDE_Integrator<"<< dim <<">::EstablishMatrixSetup: ";
              cout <<" Test operand is an Element variable. Is this intended ?"<< endl;
              (*iter).second = offset;
              switch ( (*iter).first.key.type ) {
                   case SCALAR: offset += target_.elements;
                     break;
                   case VECTOR: offset += target_.elements * dim;
                     break;
                   case TENSOR: offset += target_.elements * dim2_;
                     break;
                   case ARRAY:  offset += target_.elements * (*iter).first.key.dataDepth;
                      break;
                   case FLAGGEDARRAY:  offset += target_.elements * (*iter).first.key.dataDepth;
                }
               break;
            default:
              throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
                                           "Test operand placement unresolved");
         }
      }

   // --------------------------------------------------------------------------
   // 5. communicating offsets to MathOperators
   // --------------------------------------------------------------------------
   operandsIterator  iter;

   for ( typename map<string,MathOperatorLHS<dim>*>::iterator
         lhs_it=lhs_operators_.begin(); lhs_it!=lhs_operators_.end(); lhs_it++ )
     {
       // since basic and test operand offsets must be the same, but only the test operands
       // have been assigned an offset basic operand offsets are derived from test operand offsets
       if ( (iter=test_operands_.find((*lhs_it).second->BasicOperand())) != test_operands_.end() ) {
           (*lhs_it).second->BasicOperandOffset( (*iter).second );
           if ( (iter=basic_operands_.find((*lhs_it).second->BasicOperand())) != basic_operands_.end() )
               (*iter).second = (*lhs_it).second->BasicOperandOffset();
       }
       else
           throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
                                  "LHS basic operand matrix placement i unresolved");

        if ( (iter=test_operands_.find((*lhs_it).second->TestOperand())) != test_operands_.end() )
          (*lhs_it).second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
                                       "LHS test function operand matrix placement j unresolved");
     }
   for ( typename map<string,MathOperatorRHS<dim>*>::iterator
         rhs_it=rhs_operators_.begin(); rhs_it!=rhs_operators_.end(); rhs_it++ )
     {
        if ( (iter=test_operands_.find((*rhs_it).second->TestOperand())) != test_operands_.end() )
          (*rhs_it).second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::EstablishMatrixSetup",
                          (*rhs_it).first.c_str(), "RHS operand vector^T placement i unresolved...");
     }

   // ----------------------------------------------------------------------
   // 6. Resizing 'G' and righthand vector 'rh' which is initialised to zero
   // ----------------------------------------------------------------------
   // * Change: resize x as well
   G_.Resize( offset );
   rh_.resize( offset );
   vector<double64>( rh_ ).swap( rh_ );
   fill( rh_.begin(), rh_.end(), 0. );
   x_.resize( offset );
   vector<double64>( x_ ).swap( x_ );
#if defined(_OPENMP )
   for (size_t tid = 0 ; tid < omp_get_max_threads() ; tid++){
       this->thread_G_[tid].Resize(offset);
       this->thread_rh_[tid].resize(offset);
       vector<double64>( thread_rh_[tid] ).swap( thread_rh_[tid] );
       fill( thread_rh_[tid].begin(), thread_rh_[tid].end(), 0. );
   }
#endif
   setup_established_ = true;

 } // end EstablishMatrixSetup()






/**

AssignInitialConditions() writes the present value of the dependent
variable(s) of the computation into the righthand vector. It is applied
only in transient calculations. AssignInitialConditions() makes use
of the mesh and property storage managers.

@section implementation Implementation

AssignInitialConditions() searches the dependent variable lists
(basic operands) to retrieve those to accumulate them into the
righthand vector. If there are multiple ones, the 'offset' variable
is used to distinguish their placement. Dependent on the variable
type (scalar, vector, tensor) the accumulation process differs
since the degrees of freedom differ. Here the model dimension is used
to determine the DOFs.

@section application Application

AssignInitialConditions() is applied internally, when the PDE_Integrator
is passed to the Model object. Initial conditions are only assigned
if the model is transient.

@section messages Messages

The following errors and warnings will be reported:

@code
please call EstablishMatrixSetup() prior to this method.
@endcode

If EstablishMatrixSetup() has not been executed so far.

@code
No basic operands have been specified...
@endcode

If there are no dependent variables defined !!!

@code
Warning: So far no conditions are assigned to elements, faces, segments
@endcode

The dependent variable must be placed on the nodes:
 */
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignInitialConditions( const SIMPLICIAL_COMPLEX<dim>& gref )
 {
    size_t                  i, j;
    size_t                  position, offset;
    ScalarVariable          sc;
    VectorVariable<dim>     vc;
    TensorVariable<dim>     ts;
    Index                   prop_key;

    if ( !setup_established_ )
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignInitialConditions",
                             "please call EstablishMatrixSetup() prior to this method.");

    if ( basic_operands_.empty() ) {
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignInitialConditions",
                             "No basic operands have been specified...");
         return;
      }

     for ( operandsConstIterator
           it=test_operands_.begin(); it!=test_operands_.end(); it++ )
         {
            prop_key = (*it).first.key;
            offset   = (*it).second;
            typename vector<csmp::Node<dim>*>::const_iterator  niter(gref.NodesBegin());

            if ( prop_key.place != NODE ) {
                 throw csmp::Exception( WARNING, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignInitialConditions",
                 "So far no conditions are assigned to elements, faces, segments");
                 return;
              }

            switch( prop_key.type )
              {
                 case SCALAR:
                    while ( niter != gref.NodesEnd() ) {
                         // here the number in the new list is required to facilitate input
                         // into the size-restricted computation matrix
                         (*niter)->Read( prop_key, sc );
                         position = (*niter)->Idx() + offset;
                         rh_[ position ] *= sc();
                         niter++;
                      }
                   break;
                 case VECTOR:
                    while ( niter != gref.NodesEnd() ) {
                         (*niter)->Read( prop_key, vc );
                         for ( i=0; i<dim; i++ ) {
                             position = (*niter)->Idx() * dim + i + offset;
                             rh_[ position ] *= vc(i);
                           }
                         niter++;
                      }
                   break;
                 case TENSOR:
                    while ( niter != gref.NodesEnd() ) {
                         (*niter)->Read( prop_key, ts );
                         for ( i=0; i<dim; i++ )
                           for ( j=0; j<dim; j++ ) {
                               position = (*niter)->Idx() * dim2_ + i * dim + j + offset;
                               rh_[ position ] *= ts(i,j);
                             }
                         niter++;
                      }
                    break;
                 case ARRAY:{
                   ArrayVariable  ar(prop_key.dataDepth);
                   while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, ar );
                        for ( i=0; i<prop_key.dataDepth; i++ ) {
                            position = (*niter)->Idx() * prop_key.dataDepth  + i + offset;
                            rh_[ position ] *= ar(i);
                          }
                        niter++;
                     }
                   }
                   break;
                 case FLAGGEDARRAY:{
                  FlaggedArrayVariable far(prop_key.dataDepth);
                  while ( niter != gref.NodesEnd() ) {
                       (*niter)->Read( prop_key, far );
                       for ( i=0; i<prop_key.dataDepth; i++ ) {
                           position = (*niter)->Idx() * prop_key.dataDepth  + i + offset;
                           rh_[ position ] *= far(i);
                         }
                       niter++;
                    }
                  }
             }
         }

 } // end AssignInitialConditions





/** 
    Applies scale factor to Dirichlet matrix diagonal entries as applied by AssignEssentialConditions() and the rhs entries.
    
    Use, for instance, in mechanics problems where the unscaled matrix diagonal entry of (1) 
    is >10 orders of magnitude greater than the other matrix entries and often, this is aggrevated
    by extremely small entries in the RHS vector.
    
    @author SKM 1/10/2014
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::ScaleEssentialConditions( double64 scale_factor )
 {
    scale_factor_ = scale_factor;
   
 } // end ScaleEssentialConditions




/**

Two versions, for Modeland for Region computations are defined.
This method accumulates the basic 'Operands' of the finite-element equations
into the global solution matrix and the righthand vector. This is
done only if their condition flag (VARIABLE_FLAG) is equal to the essential
condition flag defined for this basic operand. (by default,
essential-condition flags of Operands are DIRICH(let)).

The method retrieves information from the mesh and the property managers.
If a group computation is carried out, the target group is accessed by
a pointer.

@section implementation Implementation

AssignEssentialConditions() does not decrease the size of the global
matrix 'G' or the righthand vector 'rh'. In stead, it zeros rows
in 'G' and then sets the G'ith element in question to one while the
Operand value is then placed into 'rh' (this technique is described for
fixed displacements (Dirichlet) in Segerlind, 1984, p. 417ff).

@section application Application

Because it simultaneously affects the solution matrix and the righthand
vector, AssignEssentialConditions() must be applied AFTER the accumulation
process has been completed.

When you assemble vector or tensor variables, you also have the option
of only assembling one of their components. To do this just flag the
components that you want to assemble as DIRICH.

@section messages Messages

The following errors and warnings will be reported:

@code
please call EstablishMatrixSetup() prior to this method.
@endcode

If EstablishMatrixSetup() has not been executed so far.

@code
No basic operands have been specified...
@endcode

If there are no dependent variables defined !!!

@code
Warning: So far no conditions are assigned to elements, faces, segments
@endcode

The dependent variable must be placed on the nodes:
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions( const SIMPLICIAL_COMPLEX<dim>& gref )
 {
    if ( !setup_established_ )
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions",
                      "please call EstablishMatrixSetup() prior to this method.");

    if ( basic_operands_.empty() )
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions",
                             "No (basic) operands have been specified...");

    // ----------------------------------------------------
    // 2. if the PDE_Integrator applies to a Region
    // ----------------------------------------------------
    //  This method accumulates basic 'Operands' in the finite-element equations
    //  into the global solution matrix and the righthand vector, if their condition
    //  flag is equal to the essential condition flag of the basic operand.
    //  (by default, the essential-condition flag of Operand's is DIRICH(let)).

    //  AssignEssentialConditions() must be applied AFTER accumulation process.
    //
    //  When you assemble vector or tensor variables, you also have the option
    //  of only assembling one of their components. To do this just set the
    //  components that you do not want to assemble to DBL_MAX.
    // -------------------------------------------------------------
     for ( operandsConstIterator
           it=test_operands_.begin(); it!=test_operands_.end(); it++ )
       {
         typename vector<csmp::Node<dim>*>::const_iterator  niter(gref.NodesBegin());
         csmp::Index prop_key = (*it).first.key;
         size_t      offset   = (*it).second;

         if ( prop_key.place != NODE )
           throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions",
                                 "So far no conditions are assigned to elements, faces, segments");

          switch( prop_key.type )
           {
              case SCALAR:
                while ( niter != gref.NodesEnd() ) {
                      if ( (*niter)->Status( prop_key ) == DIRICH )
                        {
                           size_t  position = (*niter)->Idx() + offset;
                           G_.ZeroRow( position );
                           G_.Add( position, position, scale_factor_ );
                           rh_[ position ] = scale_factor_ * (*niter)->Read( prop_key );
                        }
                       niter++;
                    }
                break;
              case VECTOR: {
                  VectorVariable<dim>  vc;
                  while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, vc );
                        for ( size_t i=0U; i<dim; i++ )
                          if ( vc.Flag(i) == DIRICH ) {
                               size_t  position = (*niter)->Idx() * dim + i + offset;
                               G_.ZeroRow( position );
                               G_.Add( position, position, scale_factor_ );
                               rh_[ position ] = scale_factor_ * vc(i);
                            }
                        niter++;
                     }
                  }
                break;
              case TENSOR: {
                 TensorVariable<dim>  ts;
                 while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, ts );
                        for ( size_t i=0U; i<dim; i++ )
                          if ( ts.Flag(i) == DIRICH ) 
                            for ( size_t j=0U; j<dim; j++ )
                              {
                                size_t  position = (*niter)->Idx() * dim2_ + i * dim + j + offset;
                                G_.ZeroRow( position );
                                G_.Add( position, position, scale_factor_ );
                                rh_[ position ] = scale_factor_ * ts(i,j);
                              }
                         ++niter;
                      }
                   }
                 break;
              case ARRAY: {
                  ArrayVariable  ar(prop_key.dataDepth);
                  while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, ar );
                        if ( ar.Flag() == DIRICH )
                        {
                            for ( size_t i=0U; i<prop_key.dataDepth; i++ )
                            {
                               size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                               G_.ZeroRow( position );
                               G_.Add( position, position, scale_factor_ );
                               rh_[ position ] = scale_factor_ * ar(i);
                            }
                        }
                        niter++;
                     }
                  }
                 break;
              case FLAGGEDARRAY: {
                  FlaggedArrayVariable  ar(prop_key.dataDepth);
                  while ( niter != gref.NodesEnd() ) {
                        (*niter)->Read( prop_key, ar );
                        for ( size_t i=0U; i<prop_key.dataDepth; i++ )
                          if ( ar.Flag(i) == DIRICH ) {
                               size_t  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                               G_.ZeroRow( position );
                               G_.Add( position, position, scale_factor_ );
                               rh_[ position ] = scale_factor_ * ar(i);
                            }
                        niter++;
                     }
                  }
                 break;
               default:
                 throw csmp::Exception( FATAL_ERROR,
                                       "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::AssignEssentialConditions",
                                       "Variable type not recognised by this method" );
        }
    } // end for

} // end AssignEssentialConditions






/**

Accumulate() loops over the finite-elements in the Model or target Region
and accumulates their contributions to the solution matrix and the righthand
vector. If the computation is part of a finite-difference time-stepping
scheme this is accounted for by multiplying time-dependent contributions
like source or sink rates with the time-icrement.
The method retrieves information from the mesh and the property managers.
If a group computation is carried out, the target group is accessed by
a pointer.

@section implementation Implementation

The generation and addition of the element contributions to the solution
matrix and righthand vector requires the execution of the math operator
methods:

@code
GetOperands();
ComputeContribution();
AssignToGlobal();
@endcode

If a transient calculation is chosen the "grad-test function" terms
and the Neumann boundary condition terms are multiplied
with the time-increment.

@section application  Application

Accumulate() is executed internally when the PDE_Integrator is passed to the
Model.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::Accumulate( const SIMPLICIAL_COMPLEX<dim>& gref )
 {
    // setting up the index mapping from global to local node ID numbers
    gref.RenumberNodes();
#if !defined(_OPENMP)
    // accumulating into the sparse matrix 'G'
    // ---------------------------------------
     for ( typename map<string,MathOperatorLHS<dim>*>::iterator
           it_lhs=lhs_operators_.begin(); it_lhs!=lhs_operators_.end(); it_lhs++ )
       if ( !(*it_lhs).second->AddLater() && !(*it_lhs).second->SubtractLater() )
         for ( typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
               git=gref.ElementsBegin(); git!=gref.ElementsEnd(); git++ )
           {
             (*it_lhs).second->GetOperands( *(*git) );
             (*it_lhs).second->ComputeContribution( *(*git) );
             if ( (*it_lhs).second->MultiplyWithTimeIncrement() )
               (*it_lhs).second->MultiplyWithTimeFactor( time_increment_ );
             if ( (*it_lhs).second->DivideByTimeIncrement() )
                 (*it_lhs).second->MultiplyWithTimeFactor( 1./time_increment_ );
             (*it_lhs).second->AssignToGlobal( *(*git), G_ );
           }

    // accumulating via multiplication into the righthand vector 'rhs'
    // --------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::iterator
           it_rhs=rhs_operators_.begin(); it_rhs!=rhs_operators_.end(); it_rhs++ )
       if ( !(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater() )
         for ( typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
               git=gref.ElementsBegin(); git!=gref.ElementsEnd(); git++ )
           {
             (*it_rhs).second->GetOperands( *(*git) );
             (*it_rhs).second->ComputeContribution( *(*git) );
             if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
               (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
             if ( (*it_rhs).second->DivideByTimeIncrement() )
                 (*it_rhs).second->MultiplyWithTimeFactor( 1./time_increment_ );
             (*it_rhs).second->AssignToGlobal( *(*git), rh_ );
           }
#else
     // accumulating into the sparse matrix 'G' for each thread
     // ---------------------------------------
#pragma omp parallel
    {
        FiniteElement* fe_tmp;
        size_t tid = omp_get_thread_num();
        for ( typename map<string,MathOperatorLHS<dim>*>::iterator
              it_lhs=this->thread_lhs_operators_[tid].begin(); it_lhs!=this->thread_lhs_operators_[tid].end(); it_lhs++ ){
            if ( !(*it_lhs).second->AddLater() && !(*it_lhs).second->SubtractLater() ){
#pragma omp for
                for ( int32 e = 0 ; e < gref.Elements(); e++ )
                {
//                    cout<<"element: "<<e<<endl;
                    typename SIMPLICIAL_COMPLEX<dim>::Simplex* eit = gref.E(e);
                    fe_tmp=eit->FE(); //save old pointer.
                    // change pointer here
                    eit->Assign(femgrs_[tid].E(eit->FE_Type()));

                    // do needed work.
                    (*it_lhs).second->GetOperands( (*eit) );
                    (*it_lhs).second->ComputeContribution( (*eit) );
                    if ( (*it_lhs).second->MultiplyWithTimeIncrement() )
                        (*it_lhs).second->MultiplyWithTimeFactor( time_increment_ );
                    if ( (*it_lhs).second->DivideByTimeIncrement() )
                        (*it_lhs).second->MultiplyWithTimeFactor( 1./time_increment_ );
                    (*it_lhs).second->AssignToGlobal( (*eit), thread_G_[tid] );

                    // put it back here
                    eit->Assign(fe_tmp);
                }
            }
        }
    }


#pragma omp parallel
    {
        FiniteElement* fe_tmp;
        size_t tid = omp_get_thread_num();
        // accumulating via multiplication into the righhand vector 'rhs' for each thread
        // --------------------------------------------------------------
        for ( typename map<string,MathOperatorRHS<dim>*>::iterator
              it_rhs=this->thread_rhs_operators_[tid].begin(); it_rhs!=this->thread_rhs_operators_[tid].end(); it_rhs++ ){
            if ( !(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater() ){
#pragma omp for
                for ( int32 e = 0 ; e < gref.Elements(); e++ )
                {
                    typename SIMPLICIAL_COMPLEX<dim>::Simplex* eit = gref.E(e);
                    fe_tmp=eit->FE(); //save old pointer.
                    // change pointer here
                    eit->Assign(femgrs_[tid].E(eit->FE_Type()));

                    // do needed work.
                    (*it_rhs).second->GetOperands( (*eit) );
                    (*it_rhs).second->ComputeContribution( (*eit) );
                    if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
                        (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
                    if ( (*it_rhs).second->DivideByTimeIncrement() )
                        (*it_rhs).second->MultiplyWithTimeFactor( 1./time_increment_ );

                    (*it_rhs).second->AssignToGlobal( (*eit), thread_rh_[tid] );

                    // put it back here
                    eit->Assign(fe_tmp);
                }
            }
        }
    }

    for (size_t tid = 0 ; tid < omp_get_max_threads(); tid++)
        this->G_+=this->thread_G_[tid];

    for (size_t tid = 0 ; tid < omp_get_max_threads(); tid++)
        for (size_t i = 0 ; i < this->rh_.size(); i++)
            this->rh_[i]=this->rh_[i]+this->thread_rh_[tid][i];

    //    G_.Out();
    //    for (size_t i = 0; i<rh_.size();i++)
    //        cout<<"rh:["<<i<<"]: "<<rh_[i]<<endl;
    //    exit(1);

#endif
} // end Accumulate







/**

LateAccumulate() adds source and sink terms into the righthand
vector after all other accumulations have occurred. This timing is
required if an implicit time-stepping scheme is used.
The method retrieves information from the mesh and the property managers.
If a group computation is carried out, the target group is accessed by
a pointer.

@section implementation Implementation

Those math operators which have been ear-marked for late accumulation
are added by this method. Like in Accumulate() this addition requires
the execution of the math operator methods:

@code
GetOperands();
ComputeContribution();
AssignToGlobal();
@endcode

@section application Application

LateAccumulate() is executed only in time-dependent calculations. The
execution is invoked internally, when the PDE_Integrator is passed to the
Model.
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void  PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::LateAccumulate( const SIMPLICIAL_COMPLEX<dim>& gref )
 {
    // accumulating as late addition into the righhand vector 'rhs'
    // ------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::const_iterator
           it_rhs=rhs_operators_.begin(); it_rhs!=rhs_operators_.end(); it_rhs++ )
       if ( (*it_rhs).second->AddLater() || (*it_rhs).second->SubtractLater() )
         for ( typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
               git=gref.ElementsBegin(); git!=gref.ElementsEnd(); git++ )
           {
             (*it_rhs).second->GetOperands( *(*git) );
             (*it_rhs).second->ComputeContribution( *(*git) );
             if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
               (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
             if ( (*it_rhs).second->DivideByTimeIncrement() )
               (*it_rhs).second->MultiplyWithTimeFactor( 1./time_increment_ );
             (*it_rhs).second->AssignToGlobal( *(*git), rh_ );
           }

 } // end Late Accumulate







/**

PostProcess() is a virtual PDE_Integrator interface to be implemented in a
subclass. You can use it for any kind of post-processing operations especially
those which use data from the solution matrix. In sublclasses of the older
PDE_Integrator class, post-processing is used to calculate flow velocities
(Velocity, InterstitialVelocity subclasses) or stresses and strains from
displacements (LinearElasticity subclass).
The method retrieves information from the mesh and the property managers.
If a group computation is carried out, the target group is accessed by
a pointer.

@section implementation Implementation

Define your own using the interfaces of the MeshManger, the
MemoryManager, the PropertyDatabase and the Region.

@section application Application

PostProcess() is executed internally when the PDE_Integrator is passed to
the Model.

@section message Messages

When you use the PDE_Integrator baseclass in your computations or if you
do not define PostProcess(), you will get the info message:

"no post-processing operations for REGION were defined in derived PDE_Integrator"
*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void  PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::PostProcess( const SIMPLICIAL_COMPLEX<dim>& gref )
 {
    if ( postpro_operators_.empty() ) return;

    for ( typename map<string,MathOperatorLHS<dim>*>::iterator
          it=postpro_operators_.begin(); it!=postpro_operators_.end(); it++ )
      {
         for ( size_t i=1; i<=(*it).second->ApplicationCycles(); i++ )
           {
              // setting application cylce such that it can be used by PDE operator
              (*it).second->ApplicationCycle(i);
              if (verbose_) cout <<"\nPDE_Integrator<"<<  dim;
//              if (verbose_) cout <<">::PostProcess: Computing: "<< (*it).first <<" in region'"<< gref.Name() <<"'\n";
              if (verbose_) cout <<">::PostProcess: Computing: "<< (*it).first <<"\n";
              for ( typename vector<typename SIMPLICIAL_COMPLEX<dim>::Simplex*>::const_iterator
                    git=gref.ElementsBegin(); git!=gref.ElementsEnd(); git++ )
                {
                   (*it).second->GetOperands( *(*git) );
                   (*it).second->ComputeContribution( *(*git) );
                   (*it).second->WriteOperands( *(*git) );
                }
          }
      }

 } // end PostProcess





/**

Once a solution has been computed by the current Solver object, the results
stored in the solution vector are written back to the Model or
subregions thereof. The original variable flags are preserved. 

@section implementation Implementation

The solution vector indices are re-translated into the temporary node ids.

@section application Application

OutputResults() is executed internally, when a PDE_Integrator object is passed
to the Model.

 */
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::OutputResults( SIMPLICIAL_COMPLEX<dim>& gref ) 
 {
   Index   prop_key;
   size_t  offset;

    for ( operandsIterator
          it=basic_operands_.begin(); it!=basic_operands_.end(); it++ )
     {
        typename vector<Node<dim>*>::iterator  gfirst(gref.NodesBegin());
        prop_key = (*it).first.key;
        offset   = (*it).second;

        if ( prop_key.place != NODE )
            throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::OutputResults(Model)",
                                           "only nodal properties can be output by this method.");
        switch ( prop_key.type  )
         {
            case SCALAR:
                 while ( gfirst != gref.NodesEnd() ) {
                      const double64 sc = x_[ (*gfirst)->Idx() + offset ];
                      (*gfirst)->Store( prop_key, makeScalar((*gfirst)->Status(prop_key),sc) );
                      gfirst++;
                   }
               break;
            case VECTOR: {
                 VectorVariable<dim>  vc;
                 while ( gfirst != gref.NodesEnd() ) {
                      (*gfirst)->Read( prop_key, vc );
                      for ( size_t i=0U; i<dim; i++ ) vc(i) = x_[ (*gfirst)->Idx() * dim + i + offset ];
                      (*gfirst)->Store( prop_key, vc );
                      gfirst++;
                   }
                 }
               break;
            case TENSOR: {
                 TensorVariable<dim>  ts;
                 while ( gfirst != gref.NodesEnd() ) {
                      (*gfirst)->Read( prop_key, ts );
                      for ( size_t i=0U; i<dim; i++ )
                        for ( size_t k=0U; k<dim; k++ )
                          ts(i,k) = x_[ (*gfirst)->Idx() * dim2_ + i * dim + k + offset ];
                      (*gfirst)->Store( prop_key, ts );
                      gfirst++;
                   }
                }
              break;
            case ARRAY: {
                 ArrayVariable  ar(prop_key.dataDepth);
                 while ( gfirst != gref.NodesEnd() ) {
                      (*gfirst)->Read( prop_key, ar );
                      for ( size_t i=0U; i<prop_key.dataDepth; i++ ) ar(i) = x_[ (*gfirst)->Idx() * prop_key.dataDepth + i + offset ];
                      (*gfirst)->Store( prop_key, ar );
                      gfirst++;
                   }
                 }
               break;
            case FLAGGEDARRAY: {
                 FlaggedArrayVariable  ar(prop_key.dataDepth);
                 while ( gfirst != gref.NodesEnd() ) {
                      (*gfirst)->Read( prop_key, ar );
                      for ( size_t i=0U; i<prop_key.dataDepth; i++ ) ar(i) = x_[ (*gfirst)->Idx() * prop_key.dataDepth + i + offset ];
                      (*gfirst)->Store( prop_key, ar );
                      gfirst++;
                   }
                 }
               break;
        default:
              throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::OutputResults(Model)",
                                            "Output to ARRAY type variables is not supported by this method yet.");
            
         } // end switch(type)

 } // end for

 } // end OutputResults







template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::IntegrateOver( SIMPLICIAL_COMPLEX<dim>& domain, bool debug )
 {
    // 1. configure algorithm
	clock_t loc_t = clock();
	EstablishMatrixSetup( domain );
	
    // 2. Accumulation: Note that the conditions that pertain to the group must be input ! 
    Accumulate( domain );

    // 3. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) AssignInitialConditions( domain );

    // 4. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) LateAccumulate( domain );    

    // 5. assign conditions like Dirichlet or Neumann boundary conditions etc.
    AssignEssentialConditions( domain );
	loc_t = clock() - loc_t;
	cout << "\n\n Accumulation CLOCK \t" << loc_t << "\n\n";
    // 6. diagnostics
    if ( debug ) {
         OutputGlobals();
         // OutputInput();
      }
	loc_t = clock();
    // 7. invert global matrix
    Solve();
	loc_t = clock() - loc_t;
	cout << "\n\n Solving CLOCK \t" << loc_t << "\n\n";
    // 8. write results back into Model
    OutputResults( domain );
                           
    // 9. Calculation of result-dependent properties                                 
    PostProcess( domain );

 } // end IntegrateOver




template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
void PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::WriteGlobalMatrixBitMapToText( const char* file )
 {
    char  outfile[INFO_STRING];
    strcpy( outfile, file );
    strcat( outfile, ".txt" );

     // 1. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       throw csmp::Exception( ERROR, "PDE_Integrator<dim,SIMPLICIAL_COMPLEX>::WriteGlobalMatrixBitMapToText",
                                      "Output file could not be opened" );

     // 2. writing G matrix to file
     for ( size_t i=0U; i<G_.Rows(); i++ )
       {
          for ( size_t j=0U; j<G_.Cols(); j++ )
            if ( G_.At(i,j) != 0. ) ofs << 1 <<" ";
            else                   ofs << 0 <<" ";
          ofs << endl;
       }
     ofs << endl;
     ofs.close();
     cout <<"\nPDE_Integrator<"<< dim;
     cout <<">::WriteGlobalMatrixBitMapToText: '";
     cout << outfile <<"' has been written successfully." << endl;

 } // end WriteGlobalMatrixMap




template class PDE_Integrator<1U,Region>;
template class PDE_Integrator<2U,Region>;
template class PDE_Integrator<3U,Region>;

template class PDE_Integrator<1U,Boundary>;
template class PDE_Integrator<2U,Boundary>;
template class PDE_Integrator<3U,Boundary>;

template class PDE_Integrator<1U,SplitBoundary>;
template class PDE_Integrator<2U,SplitBoundary>;
template class PDE_Integrator<3U,SplitBoundary>;

} // end namespace csmp
