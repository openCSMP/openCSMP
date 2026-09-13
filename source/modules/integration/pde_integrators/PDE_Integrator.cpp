#include "PDE_Integrator.h"
#include "Model.h"
#include "ModelSubDomain.h"
#include "Region.h"
#include "Boundary.h"
#include "Face.h"
#include "InterFace.h"
#include "SplitBoundary.h"
#include "NimbleRegion.h"
#include "ErrorHandler.h"
#include "Exception.h"

#ifdef CSMP_WITH_SAMG_SOLVER
#include "SAMG_Solver.h"
#include "SAMG_Settings.h"
#else
#include "LinearSolver.h"
#endif

#include "CompressedRowMatrix.h"
#include "SparseMatrix.h"

// #define DEBUG_CSMP_PDE_INTEGRATOR

using namespace std;

namespace csmp {

// ============================================================================
// PDE_INTEGRATOR: CONSTRUCTORS AND DESTRUCTORS
// ============================================================================

/** Default Constructor

The default constructor sets the computational domain equivalent to the
Model.
The default Solver is SAMG, without license, you only have access
to a single level solver.
The default computation is a steady-state computation. Thus, if you
want to carry out a transient calculation you have to use the method
Transient().

@attention CELLTYPE is used here because DOMAIN caused a clash
with DOMAIN defined in <cmath>
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::PDE_Integrator()
 :
#ifdef CSMP_WITH_SAMG_SOLVER
   solver_(new SAMG_Solver()),
#else
   /// add extra functionality for alternative solver if needed
   solver_(new CSMP_DEFAULT_LINEAR_SOLVER()),
#endif
   newed_Solver_object_(true),
   dof_per_node_(0),
   setup_established_(false),
   retain_matrix_(false),
   trim_vectors_(false),
   time_increment_(0.),
   verbose_(true)
{
   static_assert( !is_same<CELLTYPE<dim>,Region<dim>>::value,
                  "PDE_Integrator: template template parameter must be Element or Face");

   static_assert( !is_same<CELLTYPE<dim>,Boundary<dim>>::value,
                  "PDE_Integrator: template template parameter must be Element or Face");

   static_assert( !is_same<CELLTYPE<dim>,SplitBoundary<dim>>::value,
                  "PDE_Integrator: template template parameter must be Element or Face");

   target_.nodes = target_.elements = 0U;
}


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::PDE_Integrator( Solver& solver )
 : solver_(&solver),
   newed_Solver_object_(false),
   dof_per_node_(0),
   setup_established_(false),
   retain_matrix_(false),
   trim_vectors_(false),
   time_increment_(0.),
   verbose_(true)
{
  target_.nodes = target_.elements = 0U;
}


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::~PDE_Integrator()
 {
    if ( newed_Solver_object_ ) delete solver_;
 }

// ============================================================================
// SIMPLE ACCESSOR METHODS
// ============================================================================

template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Verbose( bool verbose )
{ verbose_=verbose; }


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
bool PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Verbose() const
{ return verbose_; }


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::RetainGlobalSolutionMatrix( bool retain )
 { retain_matrix_=retain; }


/**
    Default = false, switch on if size matters more than speed.
    TODO: is this really needed? - it should be done automatically
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::TrimExcessCapacityOfVectors( bool trim )
 {
     trim_vectors_ = trim;
 }


/**
    If a solver was allocated earlier it is deleted before the new_solver is 
    connected to the Integrator.
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::SetSolver( Solver& new_solver ) {
   if ( solver_ != &new_solver and newed_Solver_object_ ) delete solver_;
   newed_Solver_object_ = false;
   solver_ = &new_solver;
}


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
Solver& PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::GetSolver() const {
  return *solver_;
}


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
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
bool  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Transient() const
 { return !(time_increment_ < numeric_limits<double>::epsilon()); }


/** Sets the time-icrement in a transient calculation.

If one forget to set
the PDE_Integrator to transient, this is done as well.

@section application Application

Set the time-increment of an PDE_Integrator before you Apply() it to the
Model or target Region objects.
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void   PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::TimeIncrement( double dt )
 {
    time_increment_ = dt;
 }

// ============================================================================
// OUTPUT METHODS
// ============================================================================

/** Prints out the sparse solution matrix and the righthand vector

...which were accumulated by the PDE_Integrator. If you execute this method after passing
the PDE_Integrator to the Model, the matrix will have been modified by the
Solver object.

@param precision defines number of decimal places the matrix entries shall be printed with.

@attention the output can be forced to print evenly spaced integers by setting the precision to -1.
In this case, rounding is performed accordingly.

@section application  Application

To test the accumulation process by visual examination of the matrices,
you must call it directly after executing Accumulate(), see below.

*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::OutputGlobals( int32_t precision )
{
    // -----------------------------------------------------------------------
    // Path 1: Scientific notation with user-specified decimal places
    // -----------------------------------------------------------------------
    if ( precision > 0 )
    {
        // Use a local format guard so cout state is always restored
        const std::ios_base::fmtflags old_flags = cout.flags();
        const std::streamsize         old_prec  = cout.precision();

        cout << std::scientific << std::setprecision(precision);

        cout << "\nGlobal solution matrix: "
             << G_.Rows() << " x " << G_.Cols() << "\n";
        G_.Out( precision );

        cout << "\n\nGlobal righthand vector of length: "
             << rh_.size() << "\n";
        for ( size_t i{0U}; i < rh_.size(); ++i )
        {
            if ( rh_[i] >= 0.0 ) cout << " ";
            cout << rh_[i] << " ";
        }
        cout << "\n";

        cout << "\n\nGlobal solution vector of length: "
             << x_.size() << "\n";
        for ( size_t i{0U}; i < x_.size(); ++i )
        {
            if ( x_[i] >= 0.0 ) cout << " ";
            cout << x_[i] << " ";
        }
        cout << "\n";

        // Restore cout state
        cout.flags( old_flags );
        cout.precision( old_prec );
        return;
    }

    // -----------------------------------------------------------------------
    // Path 2: Integer (rounded) format  — precision <= 0
    // -----------------------------------------------------------------------
    constexpr int32_t col_stride{ 3 };
    constexpr int32_t label_len { static_cast<int32_t>(
                                      std::string_view("column").size()) };

    // Helper: safe integer print with NaN/Inf guard
    auto print_int = [&]( double v )
    {
        if ( !std::isfinite(v) ) {
            cout << " ???";
            return;
        }
        const int64_t rounded = std::lround(v);
        const int32_t digits  = static_cast<int32_t>(
                                    std::to_string(std::abs(rounded)).length());
        for ( int32_t k{0}; k < col_stride - digits; ++k ) cout << " ";
        if ( v >= 0.0 ) cout << " ";
        cout << rounded << " ";
    };

    // --- Matrix ---
    cout << "\nGlobal solution matrix (integer format): "
         << G_.Rows() << " x " << G_.Cols() << "\n";

    // Column labels
    cout << "column:     ";
    for ( size_t i{0U}; i < G_.Cols(); ++i )
    {
        const int32_t digits = static_cast<int32_t>(
                                   std::to_string(i).length());
        for ( int32_t k{0}; k < col_stride - digits; ++k ) cout << " ";
        cout << " " << i << " ";
    }
    cout << "\n";

    // Rows
    for ( size_t i{0U}; i < G_.Rows(); ++i )
    {
        // Row label with alignment
        const int32_t idx_digits = static_cast<int32_t>(
                                       std::to_string(i).length());
        const int32_t padding    = col_stride + label_len
                                 - idx_digits - 2;
        cout << "row ";
        for ( int32_t k{0}; k < padding; ++k ) cout << " ";
        cout << i << ":";

        // Matrix entries
        for ( size_t j{0U}; j < G_.Cols(); ++j )
            print_int( G_(i,j) );

        cout << "\n";
    }

    // --- RHS vector ---
    cout << "\n\nGlobal righthand vector (integer format), length "
         << rh_.size() << ":\n";
    for ( size_t i{0U}; i < rh_.size(); ++i )
    {
        if ( rh_[i] >= 0.0 ) cout << " ";
        cout << std::lround( rh_[i] ) << " ";
    }
    cout << "\n";

    // --- Solution vector ---
    cout << "\n\nGlobal solution vector (integer format), length "
         << x_.size() << ":\n";
    for ( size_t i{0U}; i < x_.size(); ++i )
    {
        if ( x_[i] >= 0.0 ) cout << " ";
        cout << std::lround( x_[i] ) << " ";
    }
    cout << "\n";

} // end OutputGlobals


// ============================================================================
// OPERATOR MANAGEMENT
// ============================================================================

/**

The PDE_Integrator object maintains a list of pointers to MathOperator instances
which represent parts of the PDE which is being modeled by the PDE_Integrator.
Add() allows you to associate such operators_ (like dove^2) with the
PDE_Integrator.

@param op pointer to PDE operator objects (=Operands) which you defined
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
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Add( MathOperatorLHS<dim,CELLTYPE>* op )
 {
    // The name for the algorithm is combined out of its operands
    lhs_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }

template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Add( MathOperatorRHS<dim,CELLTYPE>* op )
 {
    rhs_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AddBoundaryIntegral( MathOperatorRHS<dim,Face>* op )
 {
    rhs_boundary_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AddSplitBoundaryIntegral( MathOperatorRHS<dim,InterFace>* op )
 {
    rhs_split_boundary_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AddSplitBoundaryIntegral( MathOperatorLHS<dim,InterFace>* op )
 {
    lhs_split_boundary_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AddPostProcess( MathOperatorLHS<dim,CELLTYPE>* op )
 {
    postpro_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }


/**

Lists all MathOperator objects which have been defined for the PDE_Integrator.
These are output by name in alphabetical order. For the output, the
MathOperator interface method Out() is used. Thus, you influence the
information that is output when you define your own PDE operator
subclasses.
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::ListMathOperatorsLHS() const
 {
    for ( const auto& it : lhs_operators_ )
      {
         cout << it.first <<":  ";
         it.second->Out();
      }
 }


template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::ListMathOperatorsRHS() const
 {
    for ( const auto& it : rhs_operators_ )
      {
         cout << it.first <<":  ";
         it.second->Out();
      }
 }

// ============================================================================
// SOLUTION VECTOR MANAGEMENT
// ============================================================================

/**

@section application Application

Use this method to extract the solution vector from the algorithm, for
instance to use it as initial guess in another time step. (Use method
FirstGuess)
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::SolutionVector( vector<double>& sol ) const {
  sol.resize(x_.size());
  vector<double>( sol ).swap( sol );
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
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::FirstGuess( const vector<double>& guess ) {
  assert(guess.size() == x_.size());
  copy(guess.begin(), guess.end(), x_.begin());
}

// ============================================================================
// RESET METHOD
// ============================================================================

/**

After an PDE_Integrator has been applied, the solution matrix G and the righthand
vector rh contain data which need to be removed before
a new accumulation can take place. This is done automatically by the
PDE_Integrator, unless it has been asked to retain the solution matrix.

However, you may want to use new PDE operators, a new computational domain and/or computational
method. In this case, Reset() allows to restore the PDE_Integrator's default state
and zero's the elements of its matrix and vector storage.

@section application Application

Typically, you will reset an PDE_Integrator if you can just re-use it for
another task. This is always recommended, since this permits to make use
of the already initialized solution matrix, because the memory of it
(owing to the Meschach implementation) is not
de-allocated before the programme terminates.
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Reset( bool delete_math_operators )
 {
    if ( delete_math_operators ) {
         lhs_operators_.erase(  lhs_operators_.begin(),  lhs_operators_.end() );
         rhs_operators_.erase(  rhs_operators_.begin(),  rhs_operators_.end() );
         basic_operands_.erase( basic_operands_.begin(), basic_operands_.end() );
         test_operands_.erase(  test_operands_.begin(),  test_operands_.end() );
         postpro_operators_.erase( postpro_operators_.begin(), postpro_operators_.end() );
         setup_established_ = false;
      }

    // cleaning up matrix and vectors
    if ( rh_.size()          > 0 ) rh_.clear();
    if ( x_.size()           > 0 ) x_.clear();
    if ( G_.Rows()           > 0 ) G_.Erase();
    if ( DOF_indexes_.size() > 0 ) DOF_indexes_.clear();
    if ( pivotVector_.size() > 0 ) pivotVector_.clear();
    
    // boundary integrals
    boundary_faces_.clear();
    splitboundary_interfaces_.clear();

    // restoring defaults
    time_increment_ = 0.;
    retain_matrix_  = false;
    target_.nodes   = target_.elements = 0U;

} // end Reset

// ============================================================================
// SOLVER INTERFACE
// ============================================================================

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
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Solve()
 {
    solver_->Solve( G_, rh_, x_, dof_per_node_ );
 } // end Solve

// ============================================================================
// MATRIX SETUP - CRITICAL FOR MATRIX TYPE HANDLING
// ============================================================================

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

@return true if the setup was newly established and (false) if it was merely reused.

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

TODO: consider case where one might want to retain the right-hand vector, but not the matrix
TODO: rather than throwing the entire matrix away, one might just remove off-diagonal elements
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
bool PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup( const ModelSubDomain<dim,CELLTYPE>& gref, const Index& periodic_key )
 {
   // -------------------------------------------------------------------
   // 0. PDE_Integrator re-use: (and has not been reset by
   //    the user), the righthand vector and
   //    the solution vector are zeroed and nothing else is done.
   // -------------------------------------------------------------------
   if ( setup_established_ &&
        target_.nodes == gref.Nodes() &&
       !basic_operands_.empty() && !test_operands_.empty() )
     {
        if ( rh_.size() > 0 ) {
             fill( rh_.begin(), rh_.end(), 0. );
             // if any of the essential condition values changed, this must be updated
             fill( pivotVector_.begin(), pivotVector_.end(), 0. );
          }
        // provisions for the matrix class
        if ( G_.Rows() > 0 && retain_matrix_ == false ) {
             if constexpr ( is_same<MATRIXTYPE, CompressedRowMatrix>::value ) G_.Zero();
             else G_.Erase();
          }
        return false;
     }
   else gref.RenumberNodes();

   // -------------------------------------------
   // 1. determine basic sizes for G, x, rh
   // -------------------------------------------
    dof_per_node_    = 0U;
    target_.nodes    = gref.Nodes();

   // ---------------------------------------------------------------------------------
   // 2. if the PDE_Integrator is setup for first time or if its rebuild is necessary,
   //    new maps of basic and test function operands are established
   // ----------------------------------------------------------------------------------
   if ( basic_operands_.empty() || test_operands_.empty() )
     {
       const Index  unspecified;

       // lefthand MathOperators first
       // ----------------------------
       for ( auto& lhs_it : lhs_operators_ )
         {
            // making list of unique basic operands
            Index pkey = lhs_it.second->BasicOperandKey();
            if (verbose_) cout <<"\nFor: '"<< lhs_it.first <<"' PDE operator is added to lefthand operator list."<< endl;
            // checking whether the intended variables exist in the database
            //                                   Index,   calculation offset
            if ( pkey != unspecified ) basic_operands_[ lhs_it.second->BasicOperand() ] = 0U;
            else
                throw csmp::Exception( WARNING,
                                       "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup:",
                                       "lefthand basic operand not found.");
            // test function operands are picked up when the righthandside is accumulated
            // since they must also be present in there
         }

       // splitboundaries
       for ( auto& lhs_it : lhs_split_boundary_operators_ )
         {
            // making list of unique basic operands
            Index pkey = lhs_it.second->BasicOperandKey();
            if (verbose_) cout <<"\nFor: '"<< lhs_it.first <<"' PDE splitboundary operator is added to lefthand operator list."<< endl;
            // checking whether the intended variables exist in the database
            //                                   Index,   calculation offset
            if ( pkey != unspecified ) basic_operands_[ lhs_it.second->BasicOperand() ] = 0U;
            else
                throw csmp::Exception( WARNING,
                                       "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup:",
                                       "lefthand basic operand for splitboundary integral not found.");
            // test function operands are picked up when the righthandside is accumulated
            // since they must also be present in there
         }

       // righthand MathOperators
       // -----------------------
       for ( auto& rhs_it : rhs_operators_ )
         {
            // making a list of unique test operands
            Index pkey = rhs_it.second->TestOperandKey();
            if (verbose_) cout <<"\nFor: '"<< rhs_it.first <<"' PDE operator is added to righthand operator list."<< endl;
            if ( pkey != unspecified ) test_operands_[ rhs_it.second->TestOperand() ] = 0;
            else
              throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup:",
                                             "righthand test operand not found.");
         }
       // including pde operators on the model boundary
       if ( !rhs_boundary_operators_.empty() ) {
           for ( auto& rhs_it : rhs_boundary_operators_ )
             {
                // making a list of unique test operands
                Index pkey = rhs_it.second->TestOperandKey();
                if (verbose_) cout <<"\nFor: '"<< rhs_it.first <<"' PDE boundary operator is added to righthand operator list."<< endl;
                if ( pkey != unspecified ) test_operands_[ rhs_it.second->TestOperand() ] = 0;
                else
                  throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup:",
                                                 "righthand test operand for boundary integral not found.");
             }
         }
       // split boundaries
       if ( !rhs_split_boundary_operators_.empty() ) {
           for ( auto& rhs_it : rhs_split_boundary_operators_ )
             {
                // making a list of unique test operands
                Index pkey = rhs_it.second->TestOperandKey();
                if (verbose_) cout <<"\nFor: '"<< rhs_it.first <<"' PDE splitboundary operator is added to righthand operator list."<< endl;
                if ( pkey != unspecified ) test_operands_[ rhs_it.second->TestOperand() ] = 0;
                else
                  throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup:",
                                                 "righthand test operand for splitboundary integral not found.");
             }
         }

       // Testing: for each lefthand operand there must be a basic or test operand
       // ------------------------------------------------------------------------
       for ( auto& lhs_it : lhs_operators_ )
         if ( basic_operands_.find( lhs_it.second->BasicOperand() ) == basic_operands_.end() &&
              test_operands_.find( lhs_it.second->TestOperand() ) == test_operands_.end() )
           {
              cout <<"\nPDE_Integrator<"<<  dim <<">::EstablishMatrixSetup: ";
              cout <<"There is no lefthand operand corresponding to righthand operand. ";
              cout <<"\nThe system of equations is undefined. ";
              cout <<"\nCreate corresponding LHS basic or test Operand for: ";
              cout << lhs_it.first << endl;
              throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup",
                                             "lefthand basic or test operand missing");
           }

      } // end: establishing basic and test function operators
   

   // ----------------------------------------------------------------
   //   3. Offsets are assigned to test Operands
   //      indicating positions, i,j in solution matrix.
   // ----------------------------------------------------------------
   size_t         offset(0U);
   const uint32_t dim2(dim * dim);
   for ( auto& iter : test_operands_ )
     {
        // calculating the degrees of freedom per node
        if      ( iter.first.key.type == SCALAR )
          dof_per_node_ += 1U;
        else if ( iter.first.key.type == VECTOR )
          dof_per_node_ += dim;
        else if ( iter.first.key.type == TENSOR )
          dof_per_node_ += dim2;
        else if ( iter.first.key.type == ARRAY )
          dof_per_node_ += iter.first.key.dataDepth;
        else if ( iter.first.key.type == FLAGGEDARRAY )
          dof_per_node_ += iter.first.key.dataDepth;

        // the matrix/righthand offsets are stored with the operands
        // offset starts out as zero.
        switch ( iter.first.key.place )
         {
            case NODE:
              iter.second = offset;
              switch ( iter.first.key.type ) {
                   case SCALAR: offset += target_.nodes;
                     break;
                   case VECTOR: offset += target_.nodes * dim;
                     break;
                   case TENSOR: offset += target_.nodes * dim2;
                     break;
                   case ARRAY:  offset += target_.nodes * iter.first.key.dataDepth;
                     break;
                   case FLAGGEDARRAY:  offset += target_.nodes * iter.first.key.dataDepth;
                }
              break;
            case ELEMENT:
              cout <<"\nPDE_Integrator<"<< dim <<">::EstablishMatrixSetup: ";
              cout <<" Test operand is an Element variable. Is this intended ?"<< endl;
              iter.second = offset;
              switch ( iter.first.key.type ) {
                   case SCALAR: offset += target_.elements;
                     break;
                   case VECTOR: offset += target_.elements * dim;
                     break;
                   case TENSOR: offset += target_.elements * dim2;
                     break;
                   case ARRAY:  offset += target_.elements * iter.first.key.dataDepth;
                      break;
                   case FLAGGEDARRAY:  offset += target_.elements * iter.first.key.dataDepth;
                }
               break;
            default:
              throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup",
                                           "Test operand placement unresolved");
         }
      }


   // --------------------------------------------------------------------------
   // 4. communicating offsets to MathOperators
   // --------------------------------------------------------------------------
   operandsIterator  iter;

   // LEFT-HANDSIDE
   for ( auto& lhs_it : lhs_operators_ )
     {
        // since basic and test operand offsets must be the same, but only the test operands
        // have been assigned an offset basic operand offsets are derived from test operand offsets
        if ( (iter=test_operands_.find( lhs_it.second->BasicOperand())) != test_operands_.end() ) {
             lhs_it.second->BasicOperandOffset( (*iter).second );
             if ( (iter=basic_operands_.find( lhs_it.second->BasicOperand())) != basic_operands_.end() )
               (*iter).second = lhs_it.second->BasicOperandOffset();
          }
        else
          throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup",
                                       "LHS basic operand matrix placement i unresolved");

        if ( (iter=test_operands_.find(lhs_it.second->TestOperand())) != test_operands_.end() )
          lhs_it.second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup",
                                       "LHS test function operand matrix placement j unresolved");
     }
   
     //LEFT HANDSIDE - SPLIT BOUNDARY INTEGRALS
     for (auto& lhs_it : lhs_split_boundary_operators_){
       //Same passing of basic and test offsets are done for pde operators defined on a split boundary
       //Passing the basic operands to pde operators in split boundary
       if( (iter = test_operands_.find( lhs_it.second->BasicOperand() )) != test_operands_.end()  ){
         lhs_it.second->BasicOperandOffset( iter->second); //if basic operand in operator is found in the stored test operands, then we pass the offset to the operators
         if ( (iter=basic_operands_.find( lhs_it.second->BasicOperand() )) != basic_operands_.end() )
           iter->second = lhs_it.second->BasicOperandOffset();     //also assing offset calibrated in test operands to the basic operands
       } else
         throw csmp::Exception(ERROR, "PDE_Integrator::EstablishMatrixSetup",
                               "Basic operand of pde operator not found in test operands list!");

       //passing test operators to pde operators in split boundary
       if ( (iter=test_operands_.find(lhs_it.second->TestOperand() ) ) != test_operands_.end() ){
         lhs_it.second->TestOperandOffset( iter->second );
       } else
         throw csmp::Exception(ERROR, "PDE_Integrator::EstablishMatrixSetup",
                               "No matching test operand found from pde operator");

     }

   // RIGHT-HAND SIDE
   for ( auto& rhs_it : rhs_operators_ )
     {
        if ( (iter=test_operands_.find(rhs_it.second->TestOperand())) != test_operands_.end() )
          rhs_it.second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup",
                          rhs_it.first.c_str(), "RHS operand vector^T placement i unresolved...");
     }

   // RIGHT-HAND SIDE: BOUNDARY INTEGRALS
   for ( auto& rhs_it : rhs_boundary_operators_ )
     {
        if ( (iter=test_operands_.find(rhs_it.second->TestOperand())) != test_operands_.end() )
          rhs_it.second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EstablishMatrixSetup",
                                 rhs_it.first.c_str(), "RHS boundary-integral operand vector^T placement i unresolved...");
     }

   // RIGHT-HAND SIDE: SPLIT BOUNDARY INTEGRALS
     for ( auto& rhs_it : rhs_split_boundary_operators_ ){
       if ( (iter=test_operands_.find(rhs_it.second->TestOperand() ) ) != test_operands_.end()){
         rhs_it.second->TestOperandOffset( iter->second );
       } else
         throw csmp::Exception(ERROR, "PDE_Integrator::EstablishMatrixSetup",
                               "RHS split boundary pde operator contains test operand not initialised");
     }

   setup_established_ = true;

   // ------------------------------------------------------------------------------
   // 5. Sizing 'G' and righthand vector 'rh' eliminating Dirichlet rows and columns
   // ------------------------------------------------------------------------------
   EliminateEssentialConditions( gref, offset, periodic_key );

   return true;
   
 } // end EstablishMatrixSetup



// ============================================================================
// INITIAL CONDITIONS AND ESSENTIAL CONDITIONS
// ============================================================================

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

TODO: write alternative method that deals with the case when a boundary is only partially overlapping with the model domain.

 */
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AssignInitialConditions( const ModelSubDomain<dim,CELLTYPE>& gref )
 {
    if ( !setup_established_ )
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AssignInitialConditions",
                             "please call EstablishMatrixSetup() prior to this method.");

    if ( basic_operands_.empty() )
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AssignInitialConditions",
                             "No basic operands have been specified...");

    size_t position(NULL_IDX);

     for ( const auto& it : test_operands_ )
       {
          Index prop_key = it.first.key;
          if ( prop_key.place != NODE )
            throw csmp::Exception( WARNING, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AssignInitialConditions",
                                            "So far no conditions are assigned to elements, faces, segments");
          size_t offset = it.second;
          auto   niter(gref.NodesBegin());

          switch (prop_key.type)
            {
              case SCALAR:
                while (niter != gref.NodesEnd()) {
                  // here the number in the new list is required to facilitate input
                  // into the size-restricted computation matrix
                  position = (*niter)->Idx() + offset;
                  position = DOF_indexes_[position];
                  if ( position != NULL_IDX )
                    this->rh_[position] *= (*niter)->Read(prop_key);
                  niter++;
                }
                break;
              case VECTOR: {
                  VectorVariable<dim> vc;
                  while (niter != gref.NodesEnd()) {
                      (*niter)->Read(prop_key, vc);
                      for ( uint32_t i{0U}; i < dim; i++) {
                          position = (*niter)->Idx() * dim + i + offset;
                          position = DOF_indexes_[position];
                          if (position != NULL_IDX) this->rh_[position] *= vc(i);
                        }
                      niter++;
                    }
                  }
                break;
              case TENSOR: {
                  TensorVariable<dim> ts;
                  const size_t        dim2(dim * dim);
                  while (niter != gref.NodesEnd()) {
                    (*niter)->Read(prop_key, ts);
                    for ( uint32_t i{0U}; i < dim; i++)
                      for ( uint32_t j{0U}; j < dim; j++) {
                          position = (*niter)->Idx() * dim2 + i * dim + j + offset;
                          position = DOF_indexes_[position];
                          if (position != NULL_IDX) this->rh_[position] *= ts(i, j);
                        }
                    niter++;
                    }
                  }
                break;
              case ARRAY: {
                  ArrayVariable  ar(prop_key.dataDepth);
                  while (niter != gref.NodesEnd()) {
                      (*niter)->Read(prop_key, ar);
                      for ( uint32_t i{0U}; i < prop_key.dataDepth; i++) {
                          position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                          position = DOF_indexes_[position];
                          if (position != NULL_IDX) this->rh_[position] *= ar(i);
                        }
                      niter++;
                    }
                  }
                break;
              case FLAGGEDARRAY: {
                  FlaggedArrayVariable far(prop_key.dataDepth);
                  while (niter != gref.NodesEnd()) {
                      (*niter)->Read(prop_key, far);
                      for ( uint32_t i{0U}; i < prop_key.dataDepth; i++) {
                          position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                          position = DOF_indexes_[position];
                          if (position != NULL_IDX) this->rh_[position] *= far(i);
                        }
                      niter++;
                    }
                 }
               break;
             default:
                 throw csmp::Exception( WARNING, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AssignInitialConditions",
                                                 "Variable placement not recognised; nothing was done.");
        } // end switch
        
    } // end for

 } // end AssignInitialConditions



/**
    The Dirichlet constraints have already been eliminated,
    but their contributions to the non-Dirichlet rows must be added to RHS.
 
    Here the unreduced pivotVector_ indices are mapped into the reduced rh_ layout using your DOF_indexes mapping array.

    The elimination and matrix size reduction was already performed in   ReduceSystemSizeEliminatingEssentialConditions().
    
    @attention the use of DOF_indexes_ requires knowledge of variable type and offset
 
    @author Luat Khoa Tran
    @author SKM - revised for systems of equations
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim, CELLTYPE, MATRIXTYPE>::AssignEssentialConditions( const ModelSubDomain<dim,CELLTYPE>&  )
 {
    const size_t rh_size(this->rh_.size());
    for ( size_t i{0U}; i < rh_size; ++i ) {
         this->rh_[i] += pivotVector_[i];
      }
 }



// ============================================================================
// HELPER FUNCTIONS FOR BOUNDARY AND INTERFACE CHECKING
// ============================================================================

/**
    determining whether the nodes of the supplied element are contained in the computational domain
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE>
bool isContainedIn( const ModelSubDomain<dim,CELLTYPE>& comp_domain, const Face<dim>& face )
 {
    const auto nodes(face.Nodes());
    for ( uint32_t i{0U}; i<nodes; ++i )
      // since we are dealing with a Face it will always be on the domain perimeter
      if ( !comp_domain.IsPerimeterNode( face.N(i) ) ) return false;
    return true;
 }


// same but for SplitBoundary objects
template<uint32_t dim, template<uint32_t> class CELLTYPE>
bool isContainedIn( const ModelSubDomain<dim,CELLTYPE>& comp_domain, const InterFace<dim>& interface )
 {
    if constexpr ( is_same<CELLTYPE<dim>,Element<dim>>::value ) {
        if ( !comp_domain.Contains( interface.InnerParent() ) ||
             !comp_domain.Contains( interface.OuterParent() ) ) return false;
      }
    else return false;
    return true;
 }


/**
       Assumptions:
       - region "Model" contains all split boundaries
       - a split boundary is considered contained if the higher-dimensional neighbor elements of all interfaces are contained in it
*/
template<uint32_t dim>
bool isSplitBoundaryContainedInRegion( const SplitBoundary<dim>& split_boundary, const ModelSubDomain<dim,Element>& region )
 {
    // all splitboundaries are contained region "Model"
    if ( region.Name() == "Model" ) return true;
    
    // checking whether all interfaces of split_boundary are contained in comp_domain
    for ( const auto& it : split_boundary.CellVector() ) {
         if ( !region.Contains( it->InnerParent() ) ) return false;
         if ( !region.Contains( it->OuterParent() ) ) return false;
      }
    return true;
 }

// ============================================================================
// ACCUMULATION METHODS
// ============================================================================

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
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Accumulate( const ModelSubDomain<dim,CELLTYPE>& gref )
 {
    // accumulating into the sparse matrix 'G'
    // ---------------------------------------
     for ( const auto& it_lhs : lhs_operators_ )
       if ( !it_lhs.second->AddLater() && !it_lhs.second->SubtractLater() )
         for ( auto git = gref.CellsBegin(); git!=gref.CellsEnd(); ++git  )
           {
             it_lhs.second->GetOperands( *(*git) );
             it_lhs.second->ComputeContribution( *(*git) );
             if ( it_lhs.second->MultiplyWithTimeIncrement() )
               it_lhs.second->MultiplyWithTimeFactor( time_increment_ );
             it_lhs.second->AssignToGlobal(*(*git), this->G_, pivotVector_, DOF_indexes_);
           }

    // accumulating into the righhand vector 'rhs'
    // --------------------------------------------------------------
     for ( const auto& it_rhs : rhs_operators_ )
       if ( !it_rhs.second->AddLater() && !it_rhs.second->SubtractLater() )
         for ( auto git=gref.CellsBegin(); git!=gref.CellsEnd(); git++ )
           {
             it_rhs.second->GetOperands( *(*git) );
             it_rhs.second->ComputeContribution( *(*git) );
             if ( it_rhs.second->MultiplyWithTimeIncrement() )
               it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
             it_rhs.second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_);
           }

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

    @attention Method relies on a 0..n contiguous numbering of the finite element nodes.

*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::LateAccumulate( const ModelSubDomain<dim,CELLTYPE>& gref )
 {
    // accumulating as late addition into the righhand vector 'rhs'
    // ------------------------------------------------------------
     for ( const auto& it_rhs : rhs_operators_ )
       if ( it_rhs.second->AddLater() || it_rhs.second->SubtractLater() )
         for ( auto git=gref.CellsBegin(); git!=gref.CellsEnd(); git++ )
           {
             it_rhs.second->GetOperands( *(*git) );
             it_rhs.second->ComputeContribution( *(*git) );
             if ( it_rhs.second->MultiplyWithTimeIncrement() )
               it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
             it_rhs.second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_ );
           }

 } // end LateAccumulate

// ============================================================================
// BOUNDARY INTEGRAL ACCUMULATION
// ============================================================================

/**
    For RHS vector accumulation of Neumann-flagged element integrals evaluated on Face objects.
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AccumulateBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>& comp_domain,
                                                                            const Boundary<dim>& boundary )
 {
    bool did_accumulation{ false };

    // accumulating into the righhand vector 'rhs'
    // --------------------------------------------------------------
    for ( const auto& it_rhs : rhs_boundary_operators_ ) {
        if ( it_rhs.second->AddLater() || it_rhs.second->SubtractLater() )
            continue;

        const auto& key       = it_rhs.second->MaterialOperandKey();
        const bool  is_vector = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool  at_ip     = it_rhs.second->MaterialOperandPlacement() == FACE_INTEGRATION_POINT;

        for ( auto& git : boundary.CellVector() ) {

            // check NEUMANN status — most selective and cheap test first.
            // for vector variables, any component flagged NEUMANN is sufficient.
            // isContainedIn() is the most expensive check so it runs last.
            const bool is_neumann = [&]() -> bool {
                if ( is_vector ) {
                    const uint32_t depth = it_rhs.second->MaterialOperandDataDepth();
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? git->Status( 0U, key, c ) == NEUMANN
                                   : git->Status(     key, c ) == NEUMANN )
                            return true;
                    return false;
                }
                // scalar
                return at_ip ? git->Status( 0U, key ) == NEUMANN
                             : git->Status(     key ) == NEUMANN;
            }();

            if ( !is_neumann || !isContainedIn( comp_domain, *git ) )
                continue;

            it_rhs.second->GetOperands( *git );
            it_rhs.second->ComputeContribution( *git );
            if ( it_rhs.second->MultiplyWithTimeIncrement() )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *git, rh_, DOF_indexes_ );
            did_accumulation = true;
        }
    }
   
    if ( did_accumulation ) {
        cout <<"\nPDE_Integrator<"<< dim <<",CELLTYPE,MATRIXTYPE>::AccumulateBoundaryIntegrals: ";
        cout <<"accumulated integrals into RHS from '"<< boundary.Name() <<"'"<< endl;
      }

 } // AccumulateBoundaryIntegrals




/**
    Version for model subdomains other than "Model"
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AccumulateBoundaryIntegrals()
{
    // accumulate into the right-hand vector from boundary faces
    for ( const auto& it_rhs : rhs_boundary_operators_ ) {
        if ( it_rhs.second->AddLater() || it_rhs.second->SubtractLater() )
            continue;

        const auto&  key       = it_rhs.second->MaterialOperandKey();
        const bool   is_vector = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool   at_ip     = it_rhs.second->MaterialOperandPlacement() == FACE_INTEGRATION_POINT;
        const uint32_t depth   = is_vector
                                 ? it_rhs.second->MaterialOperandDataDepth()
                                 : 0U;

        for ( const auto& bit : boundary_faces_ ) {

            // lambda checks NEUMANN status — for vector variables any component suffices
            const bool is_neumann = [&]() -> bool {
                if ( is_vector ) {
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? bit->Status( 0U, key, c ) == NEUMANN
                                   : bit->Status(     key, c ) == NEUMANN )
                            return true;
                    return false;
                }
                return at_ip ? bit->Status( 0U, key ) == NEUMANN
                             : bit->Status(     key ) == NEUMANN;
            }();

            if ( !is_neumann )
                continue;

            it_rhs.second->GetOperands( *bit );
            it_rhs.second->ComputeContribution( *bit );
            if ( it_rhs.second->MultiplyWithTimeIncrement() )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *bit, rh_, DOF_indexes_ );
        }
    }

} // AccumulateBoundaryIntegrals




/**
    For source terms on the surface that need to be added to the righthand side after time or other
    conditions were multiplied in the righthand vector.

    @attention Method relies on a 0..n contiguous numbering of the finite element nodes.
 
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::LateAccumulateBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>& comp_domain,
                                                                               const Boundary<dim>& boundary )
{
    bool did_accumulation{ false };

    // late accumulation into the right-hand vector from boundary faces
    for ( const auto& it_rhs : rhs_boundary_operators_ ) {
        if ( !it_rhs.second->AddLater() && !it_rhs.second->SubtractLater() )
            continue;

        const auto&    key       = it_rhs.second->MaterialOperandKey();
        const bool     is_vector = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool     at_ip     = it_rhs.second->MaterialOperandPlacement() == FACE_INTEGRATION_POINT;
        const uint32_t depth     = is_vector
                                   ? it_rhs.second->MaterialOperandDataDepth()
                                   : 0U;

        for ( auto& git : boundary.CellVector() ) {

            // NEUMANN status check — cheap, runs before isContainedIn
            // for vector variables any component flagged NEUMANN is sufficient
            const bool is_neumann = [&]() -> bool {
                if ( is_vector ) {
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? git->Status( 0U, key, c ) == NEUMANN
                                   : git->Status(      key, c ) == NEUMANN )
                            return true;
                    return false;
                }
                return at_ip ? git->Status( 0U, key ) == NEUMANN
                             : git->Status(      key ) == NEUMANN;
            }();

            // isContainedIn is the most expensive check — runs last
            if ( !is_neumann || !isContainedIn( comp_domain, *git ) ) continue;

            it_rhs.second->GetOperands( *git );
            it_rhs.second->ComputeContribution( *git );
            if ( it_rhs.second->MultiplyWithTimeIncrement() )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *git, this->rh_, DOF_indexes_ );
            did_accumulation = true;
        }
    }

    if ( did_accumulation ) {
        cout << "\nPDE_Integrator<" << dim
             << ",CELLTYPE,MATRIXTYPE>::LateAccumulateBoundaryIntegrals: "
             << "late accumulating integrals into RHS from '"
             << boundary.Name() << "'\n";
    }

} // end LateAccumulateBoundaryIntegrals



/**
    Version for model subdomains other than "Model"
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::LateAccumulateBoundaryIntegrals()
{
    // late accumulation into the right-hand vector from boundary faces
    for ( const auto& it_rhs : rhs_boundary_operators_ ) {
        if ( !it_rhs.second->AddLater() && !it_rhs.second->SubtractLater() )
            continue;

        const auto&    key       = it_rhs.second->MaterialOperandKey();
        const bool     is_vector = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool     at_ip     = it_rhs.second->MaterialOperandPlacement() == FACE_INTEGRATION_POINT;
        const uint32_t depth     = is_vector
                                   ? it_rhs.second->MaterialOperandDataDepth()
                                   : 0U;

        for ( const auto& bit : boundary_faces_ ) {

            // NEUMANN status check — for vector variables any component suffices
            const bool is_neumann = [&]() -> bool {
                if ( is_vector ) {
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? bit->Status( 0U, key, c ) == NEUMANN
                                   : bit->Status(     key, c ) == NEUMANN )
                            return true;
                    return false;
                }
                return at_ip ? bit->Status( 0U, key ) == NEUMANN
                             : bit->Status(     key ) == NEUMANN;
            }();

            if ( !is_neumann ) continue;

            it_rhs.second->GetOperands( *bit );
            it_rhs.second->ComputeContribution( *bit );
            if ( it_rhs.second->MultiplyWithTimeIncrement() )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *bit, rh_, DOF_indexes_ );
        }
    }

} // LateAccumulateBoundaryIntegrals




// ============================================================================
// SPLIT BOUNDARY INTEGRAL ACCUMULATION
// ============================================================================

template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AccumulateSplitBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>& /*comp_domain*/,
                                                                                const SplitBoundary<dim>& splitboundary )
{
    bool did_lhs_accumulation{ false };
    bool did_rhs_accumulation{ false };

    // --- LHS: accumulate into sparse matrix G ---
    for ( const auto& it_lhs : lhs_split_boundary_operators_ ) {
        if ( it_lhs.second->AddLater() || it_lhs.second->SubtractLater() )
            continue;

        for ( InterFace<dim>* git : splitboundary.CellVector() ) {
            it_lhs.second->GetOperands( *git );
            it_lhs.second->ComputeContribution( *git );
            if ( it_lhs.second->MultiplyWithTimeIncrement() )
                it_lhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_lhs.second->AssignToGlobal( *git, this->G_, pivotVector_, DOF_indexes_ );
            did_lhs_accumulation = true;
        }
    }

    // --- RHS: accumulate into right-hand vector ---
    for ( auto& it_rhs : rhs_split_boundary_operators_ ) {
        if ( it_rhs.second->AddLater() || it_rhs.second->SubtractLater() )
            continue;

        const bool     time_increment = it_rhs.second->MultiplyWithTimeIncrement();
        const bool     is_vector      = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool     at_ip          = it_rhs.second->MaterialOperandPlacement()
                                        == INTER_FACE_INTEGRATION_POINT;
        const auto&    key            = it_rhs.second->MaterialOperandKey();
        const uint32_t depth          = is_vector
                                        ? it_rhs.second->MaterialOperandDataDepth()
                                        : 0U;

        for ( InterFace<dim>* git : splitboundary.CellVector() ) {

            // check ROBIN status — for vector variables any component suffices
            const bool is_robin = [&]() -> bool {
                if ( is_vector ) {
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? git->Status( 0U, key, c ) == ROBIN
                                   : git->Status(      key, c ) == ROBIN )
                            return true;
                    return false;
                }
                return at_ip ? git->Status( 0U, key ) == ROBIN
                             : git->Status(      key ) == ROBIN;
            }();

            if ( !is_robin )
                continue;

            it_rhs.second->GetOperands( *git );
            it_rhs.second->ComputeContribution( *git );
            if ( time_increment )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *git, pivotVector_, DOF_indexes_ );
            did_rhs_accumulation = true;
        }
    }

    if ( did_lhs_accumulation )
        cout << "\nPDE_Integrator<" << dim
             << ",CELLTYPE,MATRIXTYPE>::AccumulateSplitBoundaryIntegrals: "
             << "accumulated integrals from '"
             << splitboundary.Name() << "' into LHS.\n";

    if ( did_rhs_accumulation )
        cout << "\nPDE_Integrator<" << dim
             << ",CELLTYPE,MATRIXTYPE>::AccumulateSplitBoundaryIntegrals: "
             << "accumulated integrals from '"
             << splitboundary.Name() << "' into RHS.\n";

} // AccumulateSplitBoundaryIntegrals




template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::AccumulateSplitBoundaryIntegrals()
{
    // --- LHS: accumulate into sparse matrix G ---
    for ( const auto& it_lhs : lhs_split_boundary_operators_ ) {
        if ( it_lhs.second->AddLater() || it_lhs.second->SubtractLater() )
            continue;

        for ( const auto& git : splitboundary_interfaces_ ) {
            it_lhs.second->GetOperands( *git );
            it_lhs.second->ComputeContribution( *git );
            if ( it_lhs.second->MultiplyWithTimeIncrement() )
                it_lhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_lhs.second->AssignToGlobal( *git, this->G_, pivotVector_, DOF_indexes_ );
        }
    }

    // --- RHS: accumulate into right-hand vector ---
    for ( const auto& it_rhs : rhs_split_boundary_operators_ ) {
        if ( it_rhs.second->AddLater() || it_rhs.second->SubtractLater() )
            continue;

        const auto&    key       = it_rhs.second->MaterialOperandKey();
        const bool     is_vector = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool     at_ip     = it_rhs.second->MaterialOperandPlacement()
                                   == INTER_FACE_INTEGRATION_POINT;
        const uint32_t depth     = is_vector
                                   ? it_rhs.second->MaterialOperandDataDepth()
                                   : 0U;

        for ( const auto& git : splitboundary_interfaces_ ) {

            // check ROBIN status — for vector variables any component suffices
            const bool is_robin = [&]() -> bool {
                if ( is_vector ) {
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? git->Status( 0U, key, c ) == ROBIN
                                   : git->Status(     key, c ) == ROBIN )
                            return true;
                    return false;
                }
                return at_ip ? git->Status( 0U, key ) == ROBIN
                             : git->Status(     key ) == ROBIN;
            }();

            if ( !is_robin ) continue;

            it_rhs.second->GetOperands( *git );
            it_rhs.second->ComputeContribution( *git );
            if ( it_rhs.second->MultiplyWithTimeIncrement() )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *git, rh_, DOF_indexes_ );
        }
    }

} // AccumulateSplitBoundaryIntegrals


// ----------------------------------------------------------------------------
//  LateAccumulateSplitBoundaryIntegrals( comp_domain, splitboundary )
// ----------------------------------------------------------------------------

template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::LateAccumulateSplitBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>& /*comp_domain*/,
                                                                                    const SplitBoundary<dim>&  splitboundary )
{
    bool did_lhs_accumulation{ false };
    bool did_rhs_accumulation{ false };

    // --- LHS: late accumulate into sparse matrix G ---
    for ( const auto& it_lhs : lhs_split_boundary_operators_ ) {
        if ( !it_lhs.second->AddLater() && !it_lhs.second->SubtractLater() )
            continue;

        for ( InterFace<dim>* git : splitboundary.CellVector() ) {
            it_lhs.second->GetOperands( *git );
            it_lhs.second->ComputeContribution( *git );
            if ( it_lhs.second->MultiplyWithTimeIncrement() )
                it_lhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_lhs.second->AssignToGlobal( *git, this->G_, pivotVector_, DOF_indexes_ );
            did_lhs_accumulation = true;
        }
    }

    // --- RHS: late accumulate into right-hand vector ---
    for ( const auto& it_rhs : rhs_split_boundary_operators_ ) {
        if ( !it_rhs.second->AddLater() && !it_rhs.second->SubtractLater() ) continue;

        const bool     is_vector = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool     at_ip     = it_rhs.second->MaterialOperandPlacement()
                                   == INTER_FACE_INTEGRATION_POINT;
        const auto&    key       = it_rhs.second->MaterialOperandKey();
        const uint32_t depth     = is_vector
                                   ? it_rhs.second->MaterialOperandDataDepth()
                                   : 0U;

        for ( InterFace<dim>* git : splitboundary.CellVector() ) {

            // check ROBIN status — for vector variables any component suffices
            const bool is_robin = [&]() -> bool {
                if ( is_vector ) {
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? git->Status( 0U, key, c ) == ROBIN
                                   : git->Status(     key, c ) == ROBIN )
                            return true;
                    return false;
                }
                return at_ip ? git->Status( 0U, key ) == ROBIN
                             : git->Status(     key ) == ROBIN;
            }();

            if ( !is_robin ) continue;

            it_rhs.second->GetOperands( *git );
            it_rhs.second->ComputeContribution( *git );
            if ( it_rhs.second->MultiplyWithTimeIncrement() )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *git, pivotVector_, DOF_indexes_ );
            did_rhs_accumulation = true;
        }
    }

    if ( did_lhs_accumulation )
        cout << "\nPDE_Integrator<" << dim
             << ",CELLTYPE,MATRIXTYPE>::LateAccumulateSplitBoundaryIntegrals: "
             << "late accumulated integrals from '"
             << splitboundary.Name() << "' into LHS.\n";

    if ( did_rhs_accumulation )
        cout << "\nPDE_Integrator<" << dim
             << ",CELLTYPE,MATRIXTYPE>::LateAccumulateSplitBoundaryIntegrals: "
             << "late accumulated integrals from '"
             << splitboundary.Name() << "' into RHS.\n";

} // LateAccumulateSplitBoundaryIntegrals





template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::LateAccumulateSplitBoundaryIntegrals()
{
    // --- LHS: late accumulate into sparse matrix G ---
    for ( const auto& it_lhs : lhs_split_boundary_operators_ ) {
        if ( !it_lhs.second->AddLater() && !it_lhs.second->SubtractLater() )
            continue;

        for ( const auto& git : splitboundary_interfaces_ ) {
            it_lhs.second->GetOperands( *git );
            it_lhs.second->ComputeContribution( *git );
            if ( it_lhs.second->MultiplyWithTimeIncrement() )
                it_lhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_lhs.second->AssignToGlobal( *git, this->G_, pivotVector_, DOF_indexes_ );
        }
    }

    // --- RHS: late accumulate into right-hand vector ---
    for ( const auto& it_rhs : rhs_split_boundary_operators_ ) {
        if ( !it_rhs.second->AddLater() && !it_rhs.second->SubtractLater() ) continue;

        const auto&    key       = it_rhs.second->MaterialOperandKey();
        const bool     is_vector = it_rhs.second->MaterialOperandType() == VECTOR;
        const bool     at_ip     = it_rhs.second->MaterialOperandPlacement()
                                   == INTER_FACE_INTEGRATION_POINT;
        const uint32_t depth     = is_vector
                                   ? it_rhs.second->MaterialOperandDataDepth()
                                   : 0U;

        for ( const auto& git : splitboundary_interfaces_ ) {

            // check ROBIN status — for vector variables any component suffices
            const bool is_robin = [&]() -> bool {
                if ( is_vector ) {
                    for ( uint32_t c{ 0U }; c < depth; ++c )
                        if ( at_ip ? git->Status( 0U, key, c ) == ROBIN
                                   : git->Status(     key, c ) == ROBIN )
                            return true;
                    return false;
                }
                return at_ip ? git->Status( 0U, key ) == ROBIN
                             : git->Status(     key ) == ROBIN;
            }();

            if ( !is_robin ) continue;

            it_rhs.second->GetOperands( *git );
            it_rhs.second->ComputeContribution( *git );
            if ( it_rhs.second->MultiplyWithTimeIncrement() )
                it_rhs.second->MultiplyWithTimeFactor( time_increment_ );
            it_rhs.second->AssignToGlobal( *git, rh_, DOF_indexes_ );
        }
    }

} // LateAccumulateSplitBoundaryIntegrals





// ============================================================================
// DOMAIN COUPLING
// ============================================================================

/**

Loops over the nodes of the Region Boundary.
If these are manifolds, the coupling is applied.

@note whether the test function variable should be coupled across the interface is determined from the value of the test-function operand at the SplitBoudary.
If it is flagged ANY or PLAIN, lke at any no-flow boundary, no coupling is created, if it is ROBIN, the nodes in the manifold are coupled.

@todo perhaps introduce new flag called VARIABLE_FLAG : COUPLED to express the state of the variable at the internal SplitBoundary.

*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::CoupleDomainsAcrossSplitBoundary( ModelSubDomain<dim,CELLTYPE>& subdomain )
{
  static_assert( !is_same<CELLTYPE<dim>,SplitBoundary<dim> >::value,
                 "PDE_Integrator::CoupleDomainsAcrossSplitBoundary: only works for Region or Boundary objects" );
   
  if ( test_operands_.size() > 1U )
    throw csmp::Exception( ERROR, "PDE_Integrator::CoupleDomainsAcrossSplitBoundary",
                           "method implemented for only one test-function operand so far");

  // TODO: get this info from the solution variable
  const INDEX<SCALAR,NODE> var_key = INDEX<SCALAR,NODE>{ (*test_operands_.begin()).first.key };
  const VARIABLE_FLAG couple_if{ ROBIN };

  // manifolds can only be present at the subdomain perimeter
  for ( auto nit = subdomain.PerimeterNodesBegin(); nit != subdomain.NodesEnd(); ++nit )
    if ( (*nit)->IsManifold() )
    {
      const auto md = (*nit)->Manifold();
      const auto n_branches{ md->Branches() };

      // 1. RHS: compute total load at master dof then copy that value to slave node
      size_t masterIDX = md->N(0)->Idx();
      size_t masterPOS = DOF_indexes_[masterIDX];
      if(masterPOS != NULL_IDX) {
        for ( uint32_t n{1U}; n < n_branches; n++) {
          auto slave_node = md->N(n);
          // here the assumption is made that if the control variable value = 1, the interface should be coupled
          if ( slave_node->Status(var_key) == couple_if ) {
            size_t slaveIDX = slave_node->Idx();
            size_t slavePOS = DOF_indexes_[slaveIDX];
            if (slavePOS != NULL_IDX) rh_[masterPOS] += rh_[slavePOS];
          }
        }

        // 2. RHS: apply reciprocal coupling
        for ( uint32_t n{1U}; n < n_branches; n++) {
          auto slave_node = md->N(n);
          if ( slave_node->Status(var_key) == couple_if ) {
            size_t slaveIDX = slave_node->Idx();
            size_t slavePOS = DOF_indexes_[slaveIDX];
            if (slavePOS != NULL_IDX) rh_[slavePOS] = rh_[masterPOS];
          }
        }

        // 3. LHS: adding all the element on slave row to master dof - except slave dof
        for ( uint32_t n{1U}; n < n_branches; n++) {
          auto slave_node = md->N(n);
          if ( slave_node->Status(var_key) == couple_if ) {
            size_t slaveIDX = slave_node->Idx();
            size_t slavePOS = DOF_indexes_[slaveIDX];
            if (slavePOS != NULL_IDX) {
              for ( uint32_t j(0U); j < G_.Cols(); ++j)
                if (j != slavePOS && j != masterPOS)
                  G_.Add( static_cast<uint32_t>(masterPOS), j, G_(static_cast<uint32_t>(slavePOS), j));

              // 4. adding diagonal value to master dof
              G_.Add( static_cast<uint32_t>(masterPOS), static_cast<uint32_t>(masterPOS),
                      G_(static_cast<uint32_t>(slavePOS), static_cast<uint32_t>(slavePOS)));
            }
          }
        }

        // 5. copy that value from master dof to slave dof - except slave and master dofs position
        for ( uint32_t n{1U}; n < n_branches; n++) {
          auto slave_node = md->N(n);
          if ( slave_node->Status(var_key) == couple_if ) {
            size_t slaveIDX = slave_node->Idx();
            size_t slavePOS = DOF_indexes_[slaveIDX];
            if (slavePOS != NULL_IDX) {
              for ( uint32_t j{0}; j < G_.Cols(); ++j)
                if (j != slavePOS && j != masterPOS)
                  G_.Assign( static_cast<uint32_t>(slavePOS), j, G_( static_cast<uint32_t>(masterPOS), j));

              // 6.  copy diagonal value from master dof to slave dof
              G_.Assign( static_cast<uint32_t>(slavePOS), static_cast<uint32_t>(slavePOS),
                         G_( static_cast<uint32_t>(masterPOS), static_cast<uint32_t>(masterPOS)));
            }
          }
        }
      }
    }
} // end CoupleDomainsAcrossSplitBoundary



// ============================================================================
// POST-PROCESSING AND OUTPUT
// ============================================================================

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
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::PostProcess( ModelSubDomain<dim,CELLTYPE>& gref )
 {
    if ( postpro_operators_.empty() ) return;

    for ( const auto& it : postpro_operators_ )
      {
         for ( uint32_t i{1u}; i<=it.second->ApplicationCycles(); i++ )
           {
              // setting application cylce such that it can be used by PDE operator
              it.second->ApplicationCycle(i);
              if (verbose_) cout <<"\nPDE_Integrator<"<<  dim;
              if (verbose_) cout <<">::PostProcess: Computing: "<< it.first <<" in region'"<< gref.Name() <<"'\n";
              for ( auto git=gref.CellsBegin(); git!=gref.CellsEnd(); git++ )
                {
                   it.second->GetOperands( *(*git) );
                   it.second->ComputeContribution( *(*git) );
                   it.second->WriteOperands( *(*git) );
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

@attention Method relies on an 0..n contiguous node numbering in the region (computational domain).

 */
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::OutputResults( ModelSubDomain<dim,CELLTYPE>& gref )
 {
   size_t         position;
   const uint32_t dim2(dim * dim);

    for ( const auto& it : basic_operands_ )
     {
        auto         gfirst(gref.NodesBegin());
        const Index  prop_key = it.first.key;
        if ( prop_key.place != NODE )
            throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::OutputResults(Model)",
                                           "only nodal properties can be output by this method.");
        size_t  offset = it.second;
       
        switch ( prop_key.type  )
         {
            case SCALAR:
              while (gfirst != gref.NodesEnd()) {
                  position = (*gfirst)->Idx() + offset;
                  position = DOF_indexes_[position];
                  if (position != NULL_IDX) {
                    const double sc = this->x_[position];
                    (*gfirst)->Store(prop_key, makeScalar((*gfirst)->Status(prop_key), sc));
                  }
                  gfirst++;
                }
              break;
            case VECTOR: {
              VectorVariable<dim>  vc;
              while (gfirst != gref.NodesEnd()) {
                  (*gfirst)->Read(prop_key, vc);
                  for ( uint32_t i{0U}; i < dim; i++) {
                      position = (*gfirst)->Idx() * dim + i + offset;
                      position = DOF_indexes_[position];
                      if (position != NULL_IDX) vc(i) = this->x_[position];
                    }
                  (*gfirst)->Store(prop_key, vc);
                  gfirst++;
                }
              }
            break;
          case TENSOR: {
            TensorVariable<dim>  ts;
            while (gfirst != gref.NodesEnd()) {
                  (*gfirst)->Read(prop_key, ts);
                  for ( uint32_t i{0U}; i < dim; i++)
                    for ( uint32_t k{0U}; k < dim; k++) {
                         position = (*gfirst)->Idx() * dim2 + i * dim + k + offset;
                         position = DOF_indexes_[position];
                         if (position != NULL_IDX) ts(i, k) = this->x_[position];
                      }
                  (*gfirst)->Store(prop_key, ts);
                  gfirst++;
                }
              }
            break;
          case ARRAY: {
            ArrayVariable  ar(prop_key.dataDepth);
            while (gfirst != gref.NodesEnd()) {
                 (*gfirst)->Read(prop_key, ar);
                 for ( uint32_t i{0U}; i < prop_key.dataDepth; i++) {
                     position = (*gfirst)->Idx() * prop_key.dataDepth + i + offset;
                     position = DOF_indexes_[position];
                     if (position != NULL_IDX) ar(i) = this->x_[position];
                   }
                 (*gfirst)->Store(prop_key, ar);
                 gfirst++;
              }
            }
          break;
          case FLAGGEDARRAY: {
              FlaggedArrayVariable  ar(prop_key.dataDepth);
              while (gfirst != gref.NodesEnd()) {
                  (*gfirst)->Read(prop_key, ar);
                  for ( uint32_t i{0U}; i < prop_key.dataDepth; i++) {
                      position = (*gfirst)->Idx() * prop_key.dataDepth + i + offset;
                      position = DOF_indexes_[position];
                      if (position != NULL_IDX) ar(i) = this->x_[position];
                    }
                  (*gfirst)->Store(prop_key, ar);
                  gfirst++;
                }
            }
          break;
          default:
              throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::OutputResults(Model)",
                                            "Output to ARRAY type variables is not supported by this method yet.");
             
         } // end switch(type)

 } // end for

 } // end OutputResults

// ============================================================================
// SYSTEM REDUCTION AND MATRIX SETUP
// ============================================================================

/**
    Assembly and solution of the system of linear algebraic equations, but without consideration of boundary integrals:
    
    1. Dimensionsing of the solution matrix (n-variables, scalar or vector etc.)
    
    2. Accumulation of the FE or FV integrals
    
    3. Multiplication of the initial conditions into the RHS in the case of the transient problem.
    
    4. Late accumulation to RHS
    
    5. Essential conditions (either with or without elimination of the Dirichlet constraints from the matrix).
    
    6. Output matrix and vectors to get diagnostics
    
    7. Solution
    
    8. Output results back to moel
    
    9. Postprocessing (if respective pde operators were added to the PDE_Integrator) and writing related results to model.
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::IntegrateOver( ModelSubDomain<dim,CELLTYPE>& domain, bool debug, const Index& periodic_key )
 {
    // 1. configure algorithm
    EstablishMatrixSetup( domain, periodic_key );
    
    if ( domain.Name() != "Model" ) {
         if ( !rhs_boundary_operators_.empty() ||
              !rhs_split_boundary_operators_.empty() ||
              !lhs_split_boundary_operators_.empty() )
            throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::IntegrateOver( ModelSubDomain)",
                                  "this method only works if there are no split boundary or boundary integrals");
 
         if ( !rhs_boundary_operators_.empty() )
           throw csmp::Exception( ERROR, "PDE_Integrator<>::IntegrateOver(domain):",
                                 "integrator contains Boundary object integrals; call IntegrateOver(model,domain), to assure that these get considered." );
      }
    // 2. Accumulation of finite element integrals
    domain.RenumberNodes();
    Accumulate( domain );

    // 3. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) AssignInitialConditions( domain );

    // 4. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) LateAccumulate( domain );    

    // 5. Assign conditions like Dirichlet or Neumann boundary conditions etc.
    AssignEssentialConditions( domain );
                                                                          
    // 6. Diagnostics
    if ( debug ) {
         Out();
         OutputGlobals();
      }
 
    // 7. Solve linear algebraic system of equations
    Solve();

    // 8. Write results from the solution vector back to Model
    OutputResults( domain );
                           
    // 9. Calculation of result-dependent properties
    PostProcess( domain );

 } // end IntegrateOver

// ============================================================================
// MAIN INTEGRATION METHOD WITH BOUNDARY AND SPLIT BOUNDARY SUPPORT
// ============================================================================

/** 
    Accumulates element integrals and boundary integrals that might arise from potential Boundary and SplitBoundary objects associated.
    The domain is the computational domain to which the PDE_Integrator is applied, but this method also collects data from the entire model.
    
    @attention costly element search; use only if there are boundary integrals on other domains present
    
    @author SKM 7/7/2015
    
    TODO: detections of intersections between box boundaries and computational domains is missing
    TODO: avoid repeated check of whether faces or interfaces are contained in comp-domain, by remembering subset of these after first check
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::IntegrateOver( Model<dim>& model,
                                                  ModelSubDomain<dim,CELLTYPE>& domain,
                                                  bool debug )
 {
    // 1. configure operators and linear algebraic system
    if ( EstablishMatrixSetup( domain, model.Database().StorageKey("master_node_id") ) == true ) {
         // collecting Faces and InterFaces if boundary integrals need to be accumulated
         if ( !rhs_boundary_operators_.empty() ) {
              if constexpr ( is_same<CELLTYPE<dim>,Element<dim>>::value )
               collectBorderFacesOfComputationRegion( model, domain, boundary_faces_ );
           }
         if ( !lhs_split_boundary_operators_.empty() || !rhs_split_boundary_operators_.empty() ) {
              if constexpr ( is_same<CELLTYPE<dim>,Element<dim>>::value )
               collectInterfacesInComputationRegion( model, domain, splitboundary_interfaces_ );
           }
       }

    // 2. Accumulation: Note that the conditions that pertain to the group must be input !
    domain.RenumberNodes();
    Accumulate( domain );
    
     // 2.1 Accumulation of potential boundary integrals from Boundary objects that share nodes with the model subdomain of interest
    if ( !rhs_boundary_operators_.empty() ) {
        if ( domain.Name() == "Model" ) {
            for ( auto bit=model.BoundariesBegin(); bit!=model.BoundariesEnd(); ++bit ) {
                 const Boundary<dim>& domain_boundary = model.Boundary( (*bit).first.c_str() );
                 AccumulateBoundaryIntegrals( domain, domain_boundary );
              }
            cout << endl;
          }
        // the prestored Face objects are used
        else AccumulateBoundaryIntegrals();
      }

    // 2.2 Accumulation of potential split-boundary integrals from SplitBoundary objects inside of model subdomain of interest
    if ( !lhs_split_boundary_operators_.empty() || !rhs_split_boundary_operators_.empty() ) {
        if ( domain.Name() == "Model" ) {
            for ( auto sb=model.SplitBoundariesBegin(); sb!=model.SplitBoundariesEnd(); ++sb ) {
                 const SplitBoundary<dim>& split_boundary = model.SplitBoundary( (*sb).first.c_str() );
                 AccumulateSplitBoundaryIntegrals( domain, split_boundary );
              }
           }
         // the prestored InterFace objects are used
         else AccumulateSplitBoundaryIntegrals();
      }
   

    // 3. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) AssignInitialConditions( domain );

    // 4. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) {
         LateAccumulate( domain );
        if ( !rhs_boundary_operators_.empty() ) {
            if ( domain.Name() == "Model" ) {
                for ( auto bit=model.BoundariesBegin(); bit!=model.BoundariesEnd(); ++bit ) {
                     const Boundary<dim>& domain_boundary = model.Boundary( (*bit).first.c_str() );
                     LateAccumulateBoundaryIntegrals( domain, domain_boundary );
                  }
                cout << endl;
              }
            // the prestored Face objects are used
            else LateAccumulateBoundaryIntegrals();
          }

        // 2.2 Accumulation of potential split-boundary integrals from SplitBoundary objects inside of model subdomain of interest
        if ( !lhs_split_boundary_operators_.empty() || !rhs_split_boundary_operators_.empty() ) {
            if ( domain.Name() == "Model" ) {
                for ( auto sb=model.SplitBoundariesBegin(); sb!=model.SplitBoundariesEnd(); ++sb ) {
                     const SplitBoundary<dim>& split_boundary = model.SplitBoundary( (*sb).first.c_str() );
                     LateAccumulateSplitBoundaryIntegrals( domain, split_boundary );
                  }
               }
             // the prestored InterFace objects are used
             else LateAccumulateSplitBoundaryIntegrals();
          }
      }

    // 5. assign conditions like Dirichlet or Neumann boundary conditions etc.
    AssignEssentialConditions( domain );
    
    // Couple domains across split boundaries if continuity of the solution variable(s) is desired
    // TODO: check whether this diagnostic is the correct one?
    // ALT - EP's check is !contained_splitboundaries.empty(), contained_splitboundaries
    // having been returned by IdentifySharedSplitBoundaries( model, domain)
    if ( !lhs_split_boundary_operators_.empty() || !rhs_split_boundary_operators_.empty() )
      CoupleDomainsAcrossSplitBoundary( domain );

    // 6. diagnostics
    if ( debug ) {
         Out();
         OutputGlobals();
         // OutputInput();
      }
   
    // 7. invert global matrix
    Solve();

    // 8. write results back into Model
    OutputResults( domain );
                           
    // 9. Calculation of result-dependent properties
    PostProcess( domain );

 } // end IntegrateOver

// ============================================================================
// BOUNDARY AND SPLIT BOUNDARY IDENTIFICATION
// ============================================================================

/**
    Identifies the names of any model boundaries (standard and split ones),
    which are shared with the computational domain 'subdomain'
    by searching for the corresponding string in the boundaries.
    The identified boundaries  are returned in the argument list.
    
    @attention method assumes that the names of the boundaries are composed of the strings
    BOUNDARY, the name of the volumetric boundary that is adjacent to them and their name,
    for instance,  BOUNDARY_MATRIX_LANDSURFACE. The only other names that are handled are those
    of the boundaries in a box-shaped model.
    
    @note method is costly so that it should not be repeatedly applied, but results should be stored.
    
    @attention method also handles case, where the computational domain is the region 'Model'.
    
    @attention internal model boundaries or edge boundaries are not handled yet.
    
    @author SKM
*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
list<string> PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::IdentifySharedBoundaries( const Model<dim>& model,
                                                                     const ModelSubDomain<dim,CELLTYPE>& subdomain ) const
 {
    list<string> shared_boundaries;

    // 1. enlist any internal or free-form boundaries
    for ( auto it=model.BoundariesBegin(); it!=model.BoundariesEnd(); it++ )
      if ( (*it).first.find(subdomain.Name()) != string::npos )
        shared_boundaries.push_back( (*it).first );
   
    // 2. box-shaped model and computational domain = "Model"
    if ( subdomain.Name() == "Model" ) {
         shared_boundaries.push_back( "BOTTOM" );
         if constexpr ( dim != 1U ) {
              shared_boundaries.push_back( "LEFT" );
              shared_boundaries.push_back( "RIGHT" );
              shared_boundaries.push_back( "TOP" );
           }
         if constexpr ( dim == 3U ) {
              shared_boundaries.push_back( "FRONT" );
              shared_boundaries.push_back( "BACK" );
           }
         return shared_boundaries;
      }
      
    // 3. if the computational domain is any other subregion of the model
    else {
        for ( auto it=model.BoundariesBegin(); it!=model.BoundariesEnd(); it++ )
          if ( subdomain.SharedPerimeterNodes( (*it).second.NodesBegin(), (*it).second.NodesEnd() ) > 0 )
            shared_boundaries.push_back( (*it).first );    
      }
   
    // 3. for an entire irregularly shaped model find boundaries which touch each other
    model.Region("Model").UpdateMemberIndexes();
   
    for ( auto it=model.BoundariesBegin(); it!=model.BoundariesEnd(); it++ )
      for ( auto fit=(*it).second.CellsBegin(); fit!=(*it).second.CellsEnd(); fit++ )
        if ( subdomain.IsPerimeterCell( (*fit)->InnerParent()->Idx() ) ) {
             shared_boundaries.push_back( (*it).first );
             break;
          }
   
    // restoring a node numbering that is unique to the computational domain
    subdomain.UpdateMemberIndexes();
   
    return shared_boundaries;

 } // end IdentifySharedBoundaries


/// Encapsulated method for finding all relevant splitboundaries touching with subdomain
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
list<string>  PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::IdentifySharedSplitBoundaries( const Model<dim>& model,
                                                                                      const ModelSubDomain<dim,CELLTYPE>& subdomain ) const
 {
    list<string> contained_splitboundaries;
    if constexpr ( is_same<CELLTYPE<dim>,Element<dim>>::value )
        if ( model.SplitBoundaries() > 0U )
          for ( auto sb=model.SplitBoundariesBegin(); sb!=model.SplitBoundariesEnd(); ++sb )
            if ( isSplitBoundaryContainedInRegion( (*sb).second, subdomain ) )
              contained_splitboundaries.push_back( (*sb).first );

    return contained_splitboundaries;
}

void SetDOF(bool iscorner, csmp::VARIABLE_FLAG status, size_t position, size_t master_position, size_t offset, vector<size_t>& DOF_indexes_, size_t& DOF, vector<size_t>& DOF_masters_, vector<std::tuple<size_t,size_t>>& slave_and_master, csmp::Index periodic_key) {
  assert( position < DOF_indexes_.size() );
  if(iscorner && status == PERIODIC){
    DOF_indexes_[position] = NULL_IDX;
  }
  else{
    if (status == DIRICH){
      DOF_indexes_[position] = NULL_IDX;
    } 
    else if ( status == PERIODIC ){
      if(position == master_position){
        DOF_indexes_[position] = DOF++;
      }
      else{
        assert( master_position < DOF_indexes_.size() );
        slave_and_master.emplace_back(position, master_position);
      }
    }
    else {
      DOF_indexes_[position] = DOF++;
    }
  }
}


// ============================================================================
// SYSTEM REDUCTION - CRITICAL FOR MATRIX TYPE HANDLING
// ============================================================================

/**
    Method must be applied AFTER establish matrix setup process.
    and BEFORE the accumulation process.
 
    @attention G matrix is erased by this method
    @attention Method relies on a 0..n contiguous numbering of the finite element nodes.
 
    @author Luat Khoa Tran
    
    @attention SKM - some refactoring 31/8/22
    
    @note CRITICAL: For CompressedRowMatrix, this method calls generateSparsityPatternEliminatingEssentialConditions
          to set up the sparsity pattern. For SparseMatrix, this is not needed.

*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim, CELLTYPE, MATRIXTYPE>::EliminateEssentialConditions( const ModelSubDomain<dim,CELLTYPE>& gref,
                                                                              size_t total_degrees_of_freedom, const csmp::Index& periodic_key )
 {
    // Uses offset = total_degrees_of_freedom of system of equations, called after EstablishMatrixSetup
    DOF_indexes_.resize(total_degrees_of_freedom);
    DOF_masters_.resize(total_degrees_of_freedom);
    fill( DOF_indexes_.begin(), DOF_indexes_.end(), NULL_IDX );
    fill( DOF_masters_.begin(), DOF_masters_.end(), NULL_IDX );
    if ( trim_vectors_ ) {
      DOF_indexes_.shrink_to_fit();
      DOF_masters_.shrink_to_fit();
    }

    if (!this->setup_established_)
      throw csmp::Exception(ERROR, 
        "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EliminateEssentialConditions",
        "please call EstablishMatrixSetup() prior to this method.");

    if (this->basic_operands_.empty())
      throw csmp::Exception(ERROR, 
        "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EliminateEssentialConditions",
        "No basic operands have been specified...");

    // ATTENTION: unique node-numbers from 0..n-1 are required for this condensation
    size_t DOF{0U};

    // ========================================================================
    // BUILD DOF MAPPING
    // ========================================================================
    
    for ( auto& it : test_operands_ )
      {
        auto        niter(gref.NodesBegin());
        csmp::Index prop_key = it.first.key;
        size_t      offset = it.second;
        if (prop_key.place != NODE)
          throw csmp::Exception(ERROR, 
            "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EliminateEssentialConditions",
            "So far no conditions are assigned to elements, faces, segments");

        size_t  position{0U};
        std::vector<std::tuple<size_t,size_t>> slave_and_master;

        switch (prop_key.type) {
            case SCALAR:
              while (niter != gref.NodesEnd()) {
                position = (*niter)->Idx() + offset;
                const auto master_node_id = periodic_key.place==UNDEFINED? std::numeric_limits<size_t>::quiet_NaN():(*niter)->Read(periodic_key);
                size_t master_position = master_node_id + offset;
                DOF_masters_[position] = master_position;
                SetDOF(isCorner((*niter)->AtBoundary()), (*niter)->Status(prop_key), position, master_position, offset, DOF_indexes_, DOF, DOF_masters_, slave_and_master, periodic_key);
                niter++;
              }
              break;
            case VECTOR:
              while ( niter != gref.NodesEnd()) {
                const auto is_corner = isCorner((*niter)->AtBoundary());
                const auto master_node_id = periodic_key.place==UNDEFINED? std::numeric_limits<size_t>::quiet_NaN():(*niter)->Read(periodic_key);
                for ( uint32_t i{0U}; i < dim; ++i ) {
                  position = (*niter)->Idx() * dim + i + offset;
                  size_t master_position = master_node_id * dim + i + offset;
                  DOF_masters_[position] = master_position;
                  SetDOF(is_corner, (*niter)->Status(prop_key, i), position, master_position, offset, DOF_indexes_, DOF, DOF_masters_, slave_and_master, periodic_key);
                }
                niter++;
              }
              break;
            case TENSOR: {
                constexpr uint32_t dim2(dim * dim);
                // tensors have flags only for their diagonal elements
                while ( niter != gref.NodesEnd() ){
                  const auto is_corner = isCorner((*niter)->AtBoundary());
                  const auto master_node_id = periodic_key.place==UNDEFINED? std::numeric_limits<size_t>::quiet_NaN():(*niter)->Read(periodic_key);
                  for ( uint32_t i{0U}; i < dim; i++) {
                    for ( uint32_t j{0U}; j < dim; j++ ) {
                      position = (*niter)->Idx() * dim2 + i * dim + j + offset;
                      size_t master_position = master_node_id * dim2 + i * dim + j + offset;
                      DOF_masters_[position] = master_position;
                      SetDOF(is_corner, (*niter)->Status(prop_key, i), position, master_position, offset, DOF_indexes_, DOF, DOF_masters_, slave_and_master, periodic_key);
                    }
                  }
                  niter++;
                }
              }
              break;
            case ARRAY:
              // array variables only have a single flag
              while (niter != gref.NodesEnd()) {
                const auto is_corner = isCorner((*niter)->AtBoundary());
                const auto master_node_id = periodic_key.place==UNDEFINED? std::numeric_limits<size_t>::quiet_NaN():(*niter)->Read(periodic_key);
                for ( uint32_t i{0U}; i < prop_key.dataDepth; i++) {
                  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                  size_t master_position = master_node_id * prop_key.dataDepth + i + offset;
                  DOF_masters_[position] = master_position;
                  SetDOF(is_corner, (*niter)->Status(prop_key), position, master_position, offset, DOF_indexes_, DOF, DOF_masters_, slave_and_master, periodic_key);
                }
                niter++;
              } // end while
              break;
            case FLAGGEDARRAY:
              while ( niter != gref.NodesEnd() ) {
                const auto is_corner = isCorner((*niter)->AtBoundary());
                const auto master_node_id = periodic_key.place==UNDEFINED? std::numeric_limits<size_t>::quiet_NaN():(*niter)->Read(periodic_key);
                for ( uint32_t i{0U}; i < prop_key.dataDepth; i++ ) {
                  position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                  size_t master_position = master_node_id * prop_key.dataDepth + i + offset;
                  DOF_masters_[position] = master_position;
                  SetDOF(is_corner, (*niter)->Status(prop_key, i), position, master_position, offset, DOF_indexes_, DOF, DOF_masters_, slave_and_master, periodic_key);
                }
                niter++;
                }
           break;
            default:
              throw csmp::Exception(FATAL_ERROR,
                "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::EliminateEssentialConditions",
                "Variable type not recognised by this method");
            }
          // now slave indexes are known - we can copy the DOF indexes from the respective masters
          // this is needed for periodic boundary conditions
          for(const auto& [slave_position, master_position] : slave_and_master){
            switch (prop_key.type) {
              case SCALAR:{
                DOF_indexes_[slave_position] = DOF_indexes_[master_position];
                break;
              }
              case VECTOR:{
                for ( auto i{0U}; i < dim; ++i ) {
                  DOF_indexes_[slave_position] = DOF_indexes_[master_position];
                }
                break;
              }
              case TENSOR:{
                constexpr uint32_t dim2(dim * dim);
                for (auto i{0U}; i < dim; i++) {
                  for ( auto j{0U}; j < dim; j++ ) {
                    DOF_indexes_[slave_position] = DOF_indexes_[master_position];
                  }
                }
                break;
              }
              case ARRAY:case FLAGGEDARRAY:{
                for (auto i{0U}; i < prop_key.dataDepth; i++) {
                  DOF_indexes_[slave_position] = DOF_indexes_[master_position];
                }
                break;
              }
              default:{
                throw csmp::Exception(FATAL_ERROR,
                  "PDE_Integrator<dim,CELLTYPE>::ReduceSystemSizeEliminatingEssentialConditions",
                  "Variable type not recognised by this method");
            }
          }
        }
      } // end for (all Dirichlet flagged variables)

    // ========================================================================
    // ERASE GLOBAL SOLUTION MATRIX (SparseMatrix or CompressedRowMatrx)
    // ========================================================================
    
    this->G_.Erase();

    // ========================================================================
    // DIFFERENTIATED TREATMENT FOR MATRIX TYPES
    // ========================================================================
    
    if constexpr ( is_same<MATRIXTYPE, CompressedRowMatrix>::value )
    {
#ifndef NDEBUG
#ifdef DEBUG_CSMP_PDE_INTEGRATOR
        if ( verbose_ ) {
            cout << "\n=== PDE_Integrator::EliminateEssentialConditions Debug ===" << endl;
            cout << "test_operands size: " << test_operands_.size() << endl;
            cout << "DOF_indexes size: " << DOF_indexes_.size() << endl;
            cout << "nodes in domain: " << gref.Nodes() << endl;

            for (const auto& test_operand : test_operands_) {
                cout << "\nProcessing operand: " << test_operand.first.name << endl;
                cout << "  Offset: " << test_operand.second << endl;
                cout << "  Type: " << parseType(test_operand.first.key.type) << endl;
            }
         }
#endif
#endif
        // CompressedRowMatrix: Generate sparsity pattern
        // NOTE: This function expects DOF_indexes_ to be indexed by (node_idx + operand_offset)
        this->G_ = CompressedRowMatrix( generateSparsityPattern<dim,CELLTYPE>( test_operands_, DOF_indexes_, gref ) );

#ifndef NDEBUG
#ifdef DEBUG_CSMP_PDE_INTEGRATOR
        if ( verbose_ ) {
            cout << "\n=== PDE_Integrator::EliminateEssentialConditions (CompressedRowMatrix) ===" << endl;
            cout << "Original DOFs: " << total_degrees_of_freedom << endl;
            cout << "Reduced DOFs (non-Dirichlet): " << DOF << endl;
            cout << "Matrix dimensions: " << G_.Rows() << " x " << G_.Cols() << endl;
            cout << "Non-zero entries: " << G_.Entries() << endl;
        }
#endif
#endif
    }
    else
    {
        // SparseMatrix and other types: Just resize
        this->G_.Resize(DOF);
        
#ifndef NDEBUG
#ifdef DEBUG_CSMP_PDE_INTEGRATOR
        if ( verbose_ ) {
            cout << "\n=== EliminateEssentialConditions (SparseMatrix) ===" << endl;
            cout << "Total original DOFs: " << total_degrees_of_freedom << endl;
            cout << "Reduced DOFs (non-Dirichlet): " << DOF << endl;
            cout << "Matrix will build sparsity pattern dynamically during accumulation" << endl;
        }
#endif
#endif
    }

    // ========================================================================
    // INITIALIZE VECTORS (same for all matrix types)
    // ========================================================================
    
    this->rh_.resize(DOF);
    fill(this->rh_.begin(), this->rh_.end(), 0.);
    this->x_.resize(DOF);
    pivotVector_.resize(DOF);
    fill(pivotVector_.begin(), pivotVector_.end(), 0.);

    if ( trim_vectors_ ) {
         rh_.shrink_to_fit();
         x_.shrink_to_fit();
         pivotVector_.shrink_to_fit();
      }

 } // end EliminateEssentialConditions




// ============================================================================
// MATRIX OUTPUT AND DIAGNOSTICS
// ============================================================================

template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::WriteGlobalMatrixBitMapToText( const char* file )
 {
    char  outfile[INFO_STRING];
    strcpy( outfile, file );
    strcat( outfile, ".txt" );

     // 1. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       throw csmp::Exception( ERROR, "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::WriteGlobalMatrixBitMapToText",
                                      "Output file could not be opened" );

     // 2. writing G matrix to file
     for ( uint32_t i{0U}; i<G_.Rows(); i++ )
       {
          for ( uint32_t j{0U}; j<G_.Cols(); j++ )
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


/**
    Prints the state of the integrator to the screen.
 
    std::map<std::string,MathOperatorLHS<dim>*>  lhs_operators_;
    std::map<std::string,MathOperatorRHS<dim>*>  rhs_operators_;
    std::map<std::string,MathOperatorRHS<dim>*>  rhs_boundary_operators_;
    std::map<std::string,MathOperatorLHS<dim>*>  postpro_operators_;
    std::map<Parameter,size_t>                   basic_operands_;
    std::map<Parameter,size_t>                   test_operands_;

    MATRIXTYPE              G_;
    std::vector<double>   rh_;
    std::vector<double>   x_;
    Solver*                 solver_;

    const size_t            dim2_;
    size_t                  dof_per_node_;
    bool                    setup_established_, retain_matrix_;
    const bool              newed_Solver_object;
    double                time_increment_;

    struct SIZES {
        size_t nodes;
        size_t elements;
    } target_;

  private:

    double                scale_factor_; ///< for essential conditions
    bool                    verbose_;

*/
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::Out() const
 {
     cout <<"\n\n"<<"PDE_Integrator<"<< dim <<",CELLTYPE,MATRIXTYPE>::Out:\n";
   
     if ( basic_operands_.empty() || test_operands_.empty() ) {
          cout <<"\n\t"<<"integrator has not been initialized yet.\n";
          return;
       }
   
     cout <<"\n\t"<<"solution variables, their type and their offsets in the righthand vector:";
     //  std::map<Parameter,size_t>  test_operands_;
     for ( const auto& it : test_operands_ )
       cout <<"\n\t\t"<< it.first.name <<" ("<< parseType(it.first.key.type) <<"), offset: "<< it.second;

     cout <<"\n\n\t"<<"lefthand element integrals that will be accumulated (material, basic and test operands):";
     for ( const auto& it : lhs_operators_ ) {
          cout <<"\n\t\t"<< it.first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }
   
     cout <<"\n\n\t"<<"righthand element integrals that will be accumulated:";
     for ( const auto& it : rhs_operators_ ) {
          cout <<"\n\t\t"<< it.first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }

     cout <<"\n\n\t"<<"surface element integrals that will be accumulated:";
     for ( const auto& it : rhs_boundary_operators_ ) {
          cout <<"\n\t\t"<< it.first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }

     cout <<"\n\n\t"<<"post-processing operators:";
     for ( const auto& it : postpro_operators_ ) {
          cout <<"\n\t\t"<< it.first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }
   
      // -----------------------------------------------------------------------
      // System diagnostics
      // -----------------------------------------------------------------------
      cout << "\n\nSystem dimensions: "
           << G_.Rows() << " x " << G_.Cols() << "\n";

      cout << std::boolalpha
           << "\n\tdegrees of freedom per node:        " << dof_per_node_
           << "\n\tsetup established:                  " << setup_established_
           << "\n\tretain matrix between steps:        " << retain_matrix_
           << "\n\ttime increment:                     " << time_increment_
           << std::noboolalpha
           << "\n";

      if ( solver_ != nullptr )
          cout << "\nSolver: " << solver_->Name() << "\n";
      else
          cout << "\nSolver: (none assigned)\n";
    
 } // end Out

// ============================================================================
// HELPER FUNCTIONS FOR COLLECTING BOUNDARY FACES AND INTERFACES
// ============================================================================

/**
    Collects pointers to Faces that may be needed to accumulate surface integrals on the computational domain.
    Reports whether any faces were identified.
    
    Two cases are considered:
        1. If the computational domain is Model, all boundaries will be considered automatically and the 'boundary_faces' vector stays uninitialsed
        2. When computations are performed on a model subdomain the Faces that share higher-dimensional elements with it are recorded.
        
    @return reports on whether a Face*  vector was initialised
    
*/
template<uint32_t dim>
bool collectBorderFacesOfComputationRegion( const Model<dim>& model,
                                            const ModelSubDomain<dim,Element>& comp_domain,
                                            vector<const Face<dim>*>& boundary_faces )
  {
     // forgetting older settings
     boundary_faces.clear();
     
     // in computations on the "Model" all boundaries are considered automatically
     if ( comp_domain.Name() == "Model" ) return false;
     
     // searching existing faces to see whether their parents are situated on the boundary of
     // the computational domain
     for ( auto ft=model.Mesh().FacesBegin(); ft!=model.Mesh().FacesEnd(); ++ft ) {
          // all Faces have higher-dim element parents on the inside
          if ( comp_domain.IsPerimeterCell( (*ft).InnerParent() ) ) {
               // cast needed because of issue with constness in plf::colony
               boundary_faces.push_back( &(*ft) );
            }
          // the outside is considered, but only if there are elements there
          else if ( (*ft).OuterParent() && comp_domain.IsPerimeterCell( (*ft).OuterParent() ) ) {
                boundary_faces.push_back( &(*ft) );
            }
        }
     // vector will be unique already, but may be too big
     boundary_faces.shrink_to_fit();
       
     // reporting results
     if ( boundary_faces.empty() ) return false;
     return true;
     
  } // end collectBorderFacesOfComputationRegion

template bool collectBorderFacesOfComputationRegion( const Model<3U>&, const ModelSubDomain<3U,Element>&, vector<const Face<3U>*>& );
template bool collectBorderFacesOfComputationRegion( const Model<2U>&, const ModelSubDomain<2U,Element>&, vector<const Face<2U>*>& );
template bool collectBorderFacesOfComputationRegion( const Model<1U>&, const ModelSubDomain<1U,Element>&, vector<const Face<1U>*>& );


/**
      Same as above but for computations on boundaries that are rimmed by line faces (Edge objects).
      But method has to identify the higher-dimensional face neighbors of the line faces without being able to use inner or outer parent elements.
      
      @todo we must make sure that these line faces are not included into the computational domain! - Create new type ?
      @todo an edge could be a lower-dimensional face that knows its higher dimensional faces
       
*/
/*
template<uint32_t dim>
bool collectBorderFacesOfComputationBoundary( const Model<dim>& model, const Boundary<dim>& comp_domain,
                                              vector<const csmp::Edge<dim>*>& boundary_faces )
  {
     // forgetting older settings
     boundary_faces.clear();
     
     // in computations on the "Model" all boundaries are considered automatically
     if ( comp_domain.Name() == "Model" ) return false;
     
     // searching existing faces to see whether their parents are situated on the boundary of
     // the computational domain
     for ( auto ft=model.Mesh().EdgesBegin(); ft!=model.Mesh().EdgesEnd(); ++ft ) {
          // all Faces have higher-dim element parents on the inside
          if ( comp_domain.IsPerimeterCell( (*ft).InnerParent() ) ) {
               boundary_faces.push_back( &(*ft) );
            }
          // the outside is considered, but only if there are elements there
          else if ( (*ft).OuterParent() && comp_domain.IsPerimeterCell( (*ft).OuterParent() ) ) {
                boundary_faces.push_back( &(*ft) );
            }
        }
       
     // reporting results
     if ( boundary_faces.empty() ) return false;
     return true;
     
  } // end collectBorderFacesOfComputationDomain

template bool collectBorderFacesOfComputationBoundary( const Model<3U>&, const Boundary<3U>&, vector<const Edge<3U>*>& );
template bool collectBorderFacesOfComputationBoundary( const Model<2U>&, const Boundary<2U>&, vector<const Edge<2U>*>& );
template bool collectBorderFacesOfComputationBoundary( const Model<1U>&, const Boundary<1U>&, vector<const Edge<1U>*>& );
*/


/**
    Collects pointers to  InterFaces that may be needed to accumulate splitboundary integrals on a computational domain consisting of elements.
    Only those interfaces are considered are considered that have both higher-dim parents inside of the computational domain.
      
*/
template<uint32_t dim>
bool collectInterfacesInComputationRegion( const Model<dim>& model,
                                           const ModelSubDomain<dim,Element>& comp_domain,
                                           vector<const InterFace<dim>*>& splitboundary_interfaces )
 {
     // forgetting older settings
     splitboundary_interfaces.clear();
     
     // in computations on the "Model" all splitboundaries will be considered automatically
     if ( comp_domain.Name() == "Model" ) return false;
     
     // searching existing faces to see whether their parents are situated on the boundary of
     // the computational domain
     for ( auto ift=model.Mesh().InterfacesBegin(); ift!=model.Mesh().InterfacesEnd(); ++ift ) {
          // all Faces have higher-dim element parents on the inside
          if ( comp_domain.IsPerimeterCell( (*ift).InnerParent() ) &&
               comp_domain.IsPerimeterCell( (*ift).OuterParent() ) ) {
               splitboundary_interfaces.push_back( const_cast<InterFace<dim>*>(&(*ift)) );
            }
        }
     // vector will be unique already, but may be too big
     splitboundary_interfaces.shrink_to_fit();
       
     // reporting results
     if ( splitboundary_interfaces.empty() ) return false;
     return true;

 } // end collectInterfacesInComputationRegion

template bool collectInterfacesInComputationRegion( const Model<3U>&, const ModelSubDomain<3U,Element>&, vector<const InterFace<3U>*>& );
template bool collectInterfacesInComputationRegion( const Model<2U>&, const ModelSubDomain<2U,Element>&, vector<const InterFace<2U>*>& );
template bool collectInterfacesInComputationRegion( const Model<1U>&, const ModelSubDomain<1U,Element>&, vector<const InterFace<1U>*>& );

// ============================================================================
// EXPLICIT TEMPLATE INSTANTIATIONS
// ============================================================================

// Instantiate for all supported dimensions with default CompressedRowMatrix
template class PDE_Integrator<1U>;
template class PDE_Integrator<2U>;
template class PDE_Integrator<3U>;

// Instantiate for Face cell type with default CompressedRowMatrix
template class PDE_Integrator<1U,Face>;
template class PDE_Integrator<2U,Face>;
template class PDE_Integrator<3U,Face>;

// Instantiate for all supported dimensions with SparseMatrix
template class PDE_Integrator<1U,Element,SparseMatrix>;
template class PDE_Integrator<2U,Element,SparseMatrix>;
template class PDE_Integrator<3U,Element,SparseMatrix>;

// Instantiate for Face cell type with SparseMatrix
template class PDE_Integrator<1U,Face,SparseMatrix>;
template class PDE_Integrator<2U,Face,SparseMatrix>;
template class PDE_Integrator<3U,Face,SparseMatrix>;


// COVERED BY DEFAULT INITIALISATION

// Instantiate for all supported dimensions with CompressedRowMatrix
//template class PDE_Integrator<1U,Element,CompressedRowMatrix>;
//template class PDE_Integrator<2U,Element,CompressedRowMatrix>;
//template class PDE_Integrator<3U,Element,CompressedRowMatrix>;

// Instantiate for Face cell type with CompressedRowMatrix
//template class PDE_Integrator<1U,Face,CompressedRowMatrix>;
//template class PDE_Integrator<2U,Face,CompressedRowMatrix>;
//template class PDE_Integrator<3U,Face,CompressedRowMatrix>;

} // end namespace csmp

