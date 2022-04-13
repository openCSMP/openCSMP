#include "PDE_Integrator.h"
#include "Model.h"
#include "Region.h"
#include "Boundary.h"
#include "Node.h"
#include "NimbleRegion.h"

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

@attention COMPUTATION_DOMAIN is used here because DOMAIN caused a clash
with DOMAIN defined in <cmath>
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
PDE_Integrator<dim,COMPUTATION_DOMAIN>::PDE_Integrator()
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
   target_.nodes = target_.elements = 0U;
}





template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
PDE_Integrator<dim,COMPUTATION_DOMAIN>::PDE_Integrator( Solver& solver )
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



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
PDE_Integrator<dim,COMPUTATION_DOMAIN>::~PDE_Integrator()
 {
    if ( newed_Solver_object_ ) delete solver_;
 }



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::Verbose( bool verbose )
{ verbose_=verbose; }



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
bool PDE_Integrator<dim,COMPUTATION_DOMAIN>::Verbose() const
{ return verbose_; }



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::RetainGlobalSolutionMatrix( bool retain )
 { retain_matrix_=retain; }


/**
    Default = false, switch on if size matters more than speed.
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::TrimExcessCapacityOfVectors( bool trim )
 {
     trim_vectors_ = trim;
 }



/**
    If a solver was allocated earlier it is deleted before the new_solver is 
    connected to the Integrator.
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::SetSolver( Solver& new_solver ) {
   if ( solver_ != &new_solver and newed_Solver_object_ ) delete solver_;
   newed_Solver_object_ = false;
   solver_ = &new_solver;
}



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
Solver& PDE_Integrator<dim,COMPUTATION_DOMAIN>::GetSolver() const {
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
bool  PDE_Integrator<dim,COMPUTATION_DOMAIN>::Transient() const
 { return !(time_increment_ < numeric_limits<double>::epsilon()); }




/** Sets the time-icrement in a transient calculation.

If one forget to set
the PDE_Integrator to transient, this is done as well.

@section application Application

Set the time-increment of an PDE_Integrator before you Apply() it to the
Model or target Region objects.
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void   PDE_Integrator<dim,COMPUTATION_DOMAIN>::TimeIncrement( double dt )
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::OutputGlobals( int32_t precision )
 {
   cout <<"\nGlobal solution matrix: "<< G_.Rows() <<" x "<< G_.Cols() << endl;
   G_.Out( precision );

   cout.setf(ios::scientific);
   cout <<"\n\nGlobal righthand vector of length: "<< rh_.size() << endl;
   for ( size_t i{0U}; i<rh_.size(); i++ )
     {
        cout.precision(precision);
        if ( rh_[i] >= 0. ) cout <<" ";
        cout << rh_[i] <<" ";
     }
   cout << endl;

   cout <<"\n\nGlobal solution vector of length: "<< x_.size() << endl;
   for ( size_t i{0U}; i<x_.size(); i++ )
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::Add( MathOperatorLHS<dim>* op )
 {
    // The name for the algorithm is combined out of its operands
    lhs_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }

template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::Add( MathOperatorRHS<dim>* op )
 {
    rhs_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }


template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::AddBoundaryIntegral( MathOperatorRHS<dim>* op )
 {
    rhs_boundary_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::AddSplitBoundaryIntegral( MathOperatorRHS<dim>* op )
 {
    rhs_split_boundary_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }



template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::AddSplitBoundaryIntegral( MathOperatorLHS<dim>* op )
 {
    lhs_split_boundary_operators_[ op->Name() ] = op;
    // force update during next application
    setup_established_ = false;
 }


template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::AddPostProcess( MathOperatorLHS<dim>* op )
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::ListMathOperatorsLHS() const
 {
    typename map<string,MathOperatorLHS<dim>*>::const_iterator it;

    for ( it=lhs_operators_.begin(); it!=lhs_operators_.end(); it++ )
      {
         cout << (*it).first <<":  ";
         (*it).second->Out();
      }
 }


template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::ListMathOperatorsRHS() const
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::SolutionVector( vector<double>& sol ) const {
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::FirstGuess( const vector<double>& guess ) {
  assert(guess.size() == x_.size());
  copy(guess.begin(), guess.end(), x_.begin());
}





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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::Reset( bool delete_math_operators )
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::Solve()
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup( const COMPUTATION_DOMAIN<dim>& gref )
 {
   // -------------------------------------------------------------------
   // 0. If the algorithm is just re-used, (and has not been reset by
   //    the user), the righthand vector and
   //    the solution vector are zeroed and nothing else is done.
   // -------------------------------------------------------------------
   if ( setup_established_ &&
        target_.nodes == gref.Nodes() &&
       !basic_operands_.empty() && !test_operands_.empty() )
     {
        // TODO: one might want to retain the right-hand vector, but not the matrix
        if ( rh_.size() > 0 ) fill( rh_.begin(), rh_.end(), 0. );
        if ( G_.Rows()  > 0 && !retain_matrix_ ) {
            // TODO: rather than throwing the entire matrix away, one might just remove off-diagonal elements
            G_.Erase();
            G_.Resize( rh_.size() );
          }
        return;
     }

   // -------------------------------------------
   // 1. determine basic sizes for G, x, rh
   // -------------------------------------------
    dof_per_node_    = 0U;
    target_.nodes    = gref.Nodes();
    target_.elements = gref.Cells();


   // ---------------------------------------------------------------------------------
   // 2. if the PDE_Integrator is setup for first time or if its rebuild is necessary,
   //    new maps of basic and test function operands are established
   // ----------------------------------------------------------------------------------
   if ( basic_operands_.empty() || test_operands_.empty() )
     {
       const Index  unspecified;

       // lefthand MathOperators first
       // ----------------------------
       for ( typename map<string,MathOperatorLHS<dim>*>::iterator
             lhs_it=lhs_operators_.begin(); lhs_it!=lhs_operators_.end(); lhs_it++ )
         {
            // making list of unique basic operands
            Index pkey = (*lhs_it).second->BasicOperandKey();
            if (verbose_) cout <<"\nFor: '"<< (*lhs_it).first <<"' PDE operator is added to lefthand term list."<< endl;
            // checking whether the intended variables exist in the database
            //                                   Index,   calculation offset
            if ( pkey != unspecified ) basic_operands_[(*lhs_it).second->BasicOperand()] = 0U;
            else
                throw csmp::Exception( WARNING,
                                       "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup:",
                                       "lefthand basic operand not found.");
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
            if (verbose_) cout <<"\nFor: '"<< (*rhs_it).first <<"' PDE operator is added to righthand term list."<< endl;
            if ( pkey != unspecified ) test_operands_[(*rhs_it).second->TestOperand()] = 0;
            else
              throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup:",
                                             "righthand test operand not found.");
         }
       // including pde operators on the model boundary
       if ( !rhs_boundary_operators_.empty() ) {
           for ( typename map<string,MathOperatorRHS<dim>*>::iterator
                 rhs_it=rhs_boundary_operators_.begin(); rhs_it!=rhs_boundary_operators_.end(); rhs_it++ )
             {
                // making a list of unique test operands
                Index pkey = (*rhs_it).second->TestOperandKey();
                if (verbose_) cout <<"\nFor: '"<< (*rhs_it).first <<"' PDE boundary operator is added to righthand term list."<< endl;
                if ( pkey != unspecified ) test_operands_[(*rhs_it).second->TestOperand()] = 0;
                else
                  throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup:",
                                                 "righthand test operand not found.");
             }
         }

       // Testing: for each righthand operand there must be a basic or test operand on the LHS
       // ------------------------------------------------------------------------------------
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
              throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                                             "lefthand basic or test operand missing");
           }

      } // end: establishing basic and test function operators
   

   // ----------------------------------------------------------------
   //   3. Offsets are assigned to test Operands
   //      indicating positions, i,j in solution matrix.
   // ----------------------------------------------------------------
   size_t       offset(0U);
   const uint32_t dim2(dim * dim);
   for ( operandsIterator
         iter=test_operands_.begin(); iter!=test_operands_.end(); iter++ )
     {
        // calculating the degrees of freedom per node
        if      ( (*iter).first.key.type == SCALAR )
          dof_per_node_ += 1U;
        else if ( (*iter).first.key.type == VECTOR )
          dof_per_node_ += dim;
        else if ( (*iter).first.key.type == TENSOR )
          dof_per_node_ += dim2;
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
                   case TENSOR: offset += target_.nodes * dim2;
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
                   case TENSOR: offset += target_.elements * dim2;
                     break;
                   case ARRAY:  offset += target_.elements * (*iter).first.key.dataDepth;
                      break;
                   case FLAGGEDARRAY:  offset += target_.elements * (*iter).first.key.dataDepth;
                }
               break;
            default:
              throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                                           "Test operand placement unresolved");
         }
      }


   // --------------------------------------------------------------------------
   // 4. communicating offsets to MathOperators
   // --------------------------------------------------------------------------
   operandsIterator  iter;

   // LEFT-HANDSIDE
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
          throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                                       "LHS basic operand matrix placement i unresolved");

        if ( (iter=test_operands_.find((*lhs_it).second->TestOperand())) != test_operands_.end() )
          (*lhs_it).second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( FATAL_ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                                       "LHS test function operand matrix placement j unresolved");
     }
   
   // RIGHT-HAND SIDE
   for ( typename map<string,MathOperatorRHS<dim>*>::iterator
         rhs_it=rhs_operators_.begin(); rhs_it!=rhs_operators_.end(); rhs_it++ )
     {
        if ( (iter=test_operands_.find((*rhs_it).second->TestOperand())) != test_operands_.end() )
          (*rhs_it).second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                          (*rhs_it).first.c_str(), "RHS operand vector^T placement i unresolved...");
     }

   // RIGHT-HAND SIDE: BOUNDARY INTEGRALS
   for ( typename map<string,MathOperatorRHS<dim>*>::iterator
         rhs_it=rhs_boundary_operators_.begin(); rhs_it!=rhs_boundary_operators_.end(); rhs_it++ )
     {
        if ( (iter=test_operands_.find((*rhs_it).second->TestOperand())) != test_operands_.end() )
          (*rhs_it).second->TestOperandOffset( (*iter).second );
        else
          throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::EstablishMatrixSetup",
                                 (*rhs_it).first.c_str(), "RHS boundary-integral operand vector^T placement i unresolved...");
     }


   // ----------------------------------------------------------------------
   // 5. Resizing 'G' and righthand vector 'rh' which is initialised to zero
   // ----------------------------------------------------------------------
   G_.Resize( offset );
   rh_.resize( offset );
   if ( trim_vectors_ ) vector<double>( rh_ ).swap( rh_ );
   fill( rh_.begin(), rh_.end(), 0. );
   x_.resize( offset );
   if ( trim_vectors_ ) vector<double>( x_ ).swap( x_ );
   setup_established_ = true;

 } // end EstablishMatrixSetup






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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::AssignInitialConditions( const COMPUTATION_DOMAIN<dim>& gref )
 {
    if ( !setup_established_ )
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::AssignInitialConditions",
                             "please call EstablishMatrixSetup() prior to this method.");

    if ( basic_operands_.empty() )
      throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::AssignInitialConditions",
                             "No basic operands have been specified...");

    size_t position(NULL_IDX);

     for ( operandsConstIterator it=test_operands_.begin(); it!=test_operands_.end(); it++ )
       {
          Index prop_key = (*it).first.key;
          if ( prop_key.place != NODE )
            throw csmp::Exception( WARNING, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::AssignInitialConditions",
                                            "So far no conditions are assigned to elements, faces, segments");
          size_t offset = (*it).second;
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
                      for ( auto i{0U}; i < dim; i++) {
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
                    for ( auto i{0U}; i < dim; i++)
                      for ( auto j{0U}; j < dim; j++) {
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
                      for ( auto i{0U}; i < prop_key.dataDepth; i++) {
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
                      for ( auto i{0U}; i < prop_key.dataDepth; i++) {
                          position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                          position = DOF_indexes_[position];
                          if (position != NULL_IDX) this->rh_[position] *= far(i);
                        }
                      niter++;
                    }
                 }
               break;
             default:
                 throw csmp::Exception( WARNING, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::AssignInitialConditions",
                                                 "Variable placement not recognised; nothing was done.");
        } // end switch
        
    } // end for

 } // end AssignInitialConditions






/**
    The Dirichlet constraints have already been eliminated,
    but their contributions to the non-Dirichlet rows have to be added to the RHS.
 
    The guts of the elimation now live in EnumerateAndFixMatrixSize().

    @author Luat Khoa Tran
*/
template<uint32_t dim, template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim, COMPUTATION_DOMAIN>::AssignEssentialConditions(const COMPUTATION_DOMAIN<dim>& domain)
 {
    const size_t rh_size(this->rh_.size());
    for ( size_t i{0U}; i < rh_size; ++i ) {
         this->rh_[i] += pivotVector_[i];
      }
 }





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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::Accumulate( const COMPUTATION_DOMAIN<dim>& gref )
 {
    // accumulating into the sparse matrix 'G'
    // ---------------------------------------
     for ( typename map<string,MathOperatorLHS<dim>*>::iterator
           it_lhs=lhs_operators_.begin(); it_lhs!=lhs_operators_.end(); it_lhs++ )
       if ( !(*it_lhs).second->AddLater() && !(*it_lhs).second->SubtractLater() )
         for ( auto git=gref.CellsBegin(); git!=gref.CellsEnd(); git++ )
           {
             (*it_lhs).second->GetOperands( *(*git) );
             (*it_lhs).second->ComputeContribution( *(*git) );
             if ( (*it_lhs).second->MultiplyWithTimeIncrement() )
               (*it_lhs).second->MultiplyWithTimeFactor( time_increment_ );
             // LUAT KHOA TRAN - (*it_lhs).second->AssignToGlobal( *(*git), G_ );
             (*it_lhs).second->AssignToGlobal(*(*git), this->G_, pivotVector_, DOF_indexes_);
           }

    // accumulating into the righhand vector 'rhs'
    // --------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::iterator
           it_rhs=rhs_operators_.begin(); it_rhs!=rhs_operators_.end(); it_rhs++ )
       if ( !(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater() )
         for ( auto git=gref.CellsBegin(); git!=gref.CellsEnd(); git++ )
           {
             (*it_rhs).second->GetOperands( *(*git) );
             (*it_rhs).second->ComputeContribution( *(*git) );
             if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
               (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
             // LUAT KHOA TRAN - (*it_rhs).second->AssignToGlobal( *(*git), rh_ );
             (*it_rhs).second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_);
           }

 } // end Accumulate





/**
    determining whether the nodes of the supplied element are contained in the computational domain
*/
template<uint32_t dim, template<uint32_t> class COMPUTATION_DOMAIN>
bool isContainedIn( const COMPUTATION_DOMAIN<dim>& comp_domain, const Face<dim>& face )
 {
    const size_t nodes(face.Nodes());
    for ( auto i{0U}; i<nodes; ++i )
      if ( !comp_domain.IsPerimeterNode( face.N(i) ) ) return false;
    return true;
 }

// same but for SplitBoundary objects
template<uint32_t dim, template<uint32_t> class COMPUTATION_DOMAIN>
bool isContainedIn( const COMPUTATION_DOMAIN<dim>& comp_domain, const InterFace<dim>& interface )
 {
    const size_t nodes(interface.Nodes());
    for ( auto i{0U}; i<nodes; ++i )
      if ( !comp_domain.IsPerimeterNode( interface.N(i) ) ) return false;
    return true;
 }



/**
    for the accumulation of Neumann-flagged element integrals evaluated on Face objects
    
    TODO: include boundary conditions applied to LHS
    TODO: adopt method to handle InterFace objects (in split boundaries) as well
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::AccumulateBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>& comp_domain,
                                                                                       const Boundary<dim>& boundary )
 {
    // accumulating into the righhand vector 'rhs'
    // --------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::iterator
           it_rhs=rhs_boundary_operators_.begin(); it_rhs!=rhs_boundary_operators_.end(); it_rhs++ )
       if ( !(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater() )
         for ( auto git=boundary.CellsBegin(); git!=boundary.CellsEnd(); git++ )
           // if the material operand is flagged Neumann, an accumulation will be performed
           // @attention it is assumed that the material operand has the same status at all integration points
           if ( isContainedIn( comp_domain, *(*git) ) == true && (
                ( (*it_rhs).second->MaterialOperandPlacement() == FACE && (*git)->Status( (*it_rhs).second->MaterialOperandKey() ) == NEUMANN ) ||
                ( (*it_rhs).second->MaterialOperandPlacement() == FACE_INTEGRATION_POINT && (*git)->Status( 0U, (*it_rhs).second->MaterialOperandKey() ) == NEUMANN ) ) )
             {
               (*it_rhs).second->GetOperands( *(*git) );
               (*it_rhs).second->ComputeContribution( *(*git) );
               if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
                 (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
               (*it_rhs).second->AssignToGlobal( *(*git), rh_, DOF_indexes_ );
             }
   
 } // AccumulateBoundaryIntegrals



/**
    Coupling conditions / surface integrals applied to SplitBoundary objects
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::AccumulateSplitBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>& comp_domain,
                                                                                            const SplitBoundary<dim>& boundary )
 {
    // TODO: include accumulation for LHS operators
    
    // accumulating into the righhand vector 'rhs'
    // --------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::iterator
           it_rhs=rhs_split_boundary_operators_.begin(); it_rhs!=rhs_split_boundary_operators_.end(); it_rhs++ )
       if ( !(*it_rhs).second->AddLater() && !(*it_rhs).second->SubtractLater() )
         for ( auto git=boundary.CellsBegin(); git!=boundary.CellsEnd(); git++ )
           // if the material operand is flagged Robin, the accumulation will be performed
           // @attention it is assumed that the material operand has the same status at all integration points
           if ( isContainedIn( comp_domain, *(*git) ) == true && (
                ( (*it_rhs).second->MaterialOperandPlacement() == INTER_FACE && 
                  (*git)->Status( (*it_rhs).second->MaterialOperandKey() ) == ROBIN ) ||
                ( (*it_rhs).second->MaterialOperandPlacement() == INTER_FACE_INTEGRATION_POINT && 
                  (*git)->Status( 0U, (*it_rhs).second->MaterialOperandKey() ) == ROBIN ) ) )
             {
               (*it_rhs).second->GetOperands( *(*git) );
               (*it_rhs).second->ComputeContribution( *(*git) );
               if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
                 (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
               (*it_rhs).second->AssignToGlobal( *(*git), rh_, DOF_indexes_ );
             }
   
 } // AccumulateSplitBoundaryIntegrals




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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::LateAccumulate( const COMPUTATION_DOMAIN<dim>& gref )
 {
    // accumulating as late addition into the righhand vector 'rhs'
    // ------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::const_iterator
           it_rhs=rhs_operators_.begin(); it_rhs!=rhs_operators_.end(); it_rhs++ )
       if ( (*it_rhs).second->AddLater() || (*it_rhs).second->SubtractLater() )
         for ( auto git=gref.CellsBegin(); git!=gref.CellsEnd(); git++ )
           {
             (*it_rhs).second->GetOperands( *(*git) );
             (*it_rhs).second->ComputeContribution( *(*git) );
             if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
               (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
             // LUAT KHOA TRAN - (*it_rhs).second->AssignToGlobal( *(*git), rh_ );
             (*it_rhs).second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_ );
           }

 } // end Late Accumulate





/**
    For source terms on the surface that need to be added to the righthand side after time or other
    conditions were multiplied in the righthand vector.

    @attention Method relies on a 0..n contiguous numbering of the finite element nodes.
 
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::LateAccumulateBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>& comp_domain,
                                                                                           const Boundary<dim>& boundary )
 {
    // accumulating as late addition into the righhand vector 'rhs'
    // ------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::const_iterator
           it_rhs=rhs_boundary_operators_.begin(); it_rhs!=rhs_boundary_operators_.end(); it_rhs++ )
       if ( (*it_rhs).second->AddLater() || (*it_rhs).second->SubtractLater() )
         for ( auto git=boundary.CellsBegin(); git!=boundary.CellsEnd(); git++ )
           // accumulations need to be performed only where material operands are flagged Neumann
           if ( isContainedIn( comp_domain, *(*git) ) == true && (
                ( (*it_rhs).second->MaterialOperandPlacement() == FACE && 
                  (*git)->Status( (*it_rhs).second->MaterialOperandKey() ) == NEUMANN ) ||
                ( (*it_rhs).second->MaterialOperandPlacement() == FACE_INTEGRATION_POINT && (*git)->Status( 0U, (*it_rhs).second->MaterialOperandKey() ) == NEUMANN ) ) )
             {
               (*it_rhs).second->GetOperands( *(*git) );
               (*it_rhs).second->ComputeContribution( *(*git) );
               if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
                 (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
               // LUAT KHOA TRAN - (*it_rhs).second->AssignToGlobal( *(*git), rh_ );
               (*it_rhs).second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_ );
             }

 } // end LateAccumulateBoundaryIntegrals


// same but for SplitBoundary objects
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::LateAccumulateSplitBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>& comp_domain,
                                                                                                const SplitBoundary<dim>& boundary )
 {
    // accumulating as late addition into the righhand vector 'rhs'
    // ------------------------------------------------------------
     for ( typename map<string,MathOperatorRHS<dim>*>::const_iterator
           it_rhs=rhs_boundary_operators_.begin(); it_rhs!=rhs_boundary_operators_.end(); it_rhs++ )
       if ( (*it_rhs).second->AddLater() || (*it_rhs).second->SubtractLater() )
         for ( auto git=boundary.CellsBegin(); git!=boundary.CellsEnd(); git++ )
           // accumulations need to be performed only where material operands are flagged Neumann
           if ( isContainedIn( comp_domain, *(*git) ) == true && (
                ( (*it_rhs).second->MaterialOperandPlacement() == INTER_FACE && 
                  (*git)->Status( (*it_rhs).second->MaterialOperandKey() ) == ROBIN ) ||
                ( (*it_rhs).second->MaterialOperandPlacement() == INTER_FACE_INTEGRATION_POINT &&  
                  (*git)->Status( 0U, (*it_rhs).second->MaterialOperandKey() ) == NEUMANN ) ) )
             {
               (*it_rhs).second->GetOperands( *(*git) );
               (*it_rhs).second->ComputeContribution( *(*git) );
               if ( (*it_rhs).second->MultiplyWithTimeIncrement() )
                 (*it_rhs).second->MultiplyWithTimeFactor( time_increment_ );
               (*it_rhs).second->AssignToGlobal(*(*git), this->rh_, DOF_indexes_ );
             }

 } // end LateAccumulateBoundaryIntegrals



/**

Loops over the nodes of the Region Boundary.
If these are manifolds, the coupling is applied.

@note whether the test function variable should be coupled across the interface is determined from the value of the test-function operand at the SplitBoudary.
If it is flagged ANY or PLAIN, lke at any no-flow boundary, no coupling is created, if it is ROBIN, the nodes in the manifold are coupled.

@todo perhaps introduce new flag called VARIABLE_FLAG : COUPLED to express the state of the variable at the internal SplitBoundary.

*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::CoupleDomainsAcrossSplitBoudary( COMPUTATION_DOMAIN<dim>& subdomain )
{
   if ( test_operands_.size() > 1U )
     throw csmp::Exception( ERROR, "PDE_Integrator::CoupleDomainsAcrossSplitBoudary",
                           "method implemented for only one test-function operand so far");
     
    // TODO: get this info from the solution variable
    // const INDEX<SCALAR,NODE> key_continuous_p = INDEX<SCALAR,NODE>( model.Database().StorageKey("pressure continuity status") );
    const INDEX<SCALAR,NODE> var_key = INDEX<SCALAR,NODE>{ (*test_operands_.begin()).first.key };
    const VARIABLE_FLAG couple_if{ ROBIN };

    // manifolds can only be present at the subdomain perimeter
    for ( auto mit = subdomain.PerimeterNodesBegin(); mit != subdomain.NodesEnd(); ++mit )
      if ( (*mit)->IsManifold() )
        {
          const auto n_branches{ (*mit)->Manifold()->Branches() };
          
          // 1. RHS: compute total load at master dof then copy that value to slave node
          size_t masterIDX = (*mit)->Manifold()->N(0)->Idx();
          for ( auto n{1U}; n < n_branches; n++ ) {
                auto slave_node = (*mit)->Manifold()->N(n);
                // here the assumption is made that if the control variable value = 1, the interface should be coupled
                // int continuous_p = static_cast<int>(slave_node->Read(key_continuous_p)); - use ROBIN status instead
                if ( slave_node->Status(var_key) == couple_if )
                  rh_[masterIDX] += rh_[slave_node->Idx()];
            }
            
          // 2. RHS: apply reciprocal coupling
          for ( auto n{1U}; n < n_branches; n++ ) {
                auto slave_node = (*mit)->Manifold()->N(n);
                if ( slave_node->Status(var_key) == couple_if )
                  rh_[slave_node->Idx()] = rh_[masterIDX];
             }

          // 3. LHS: adding all the element on slave row to master dof - except slave dof
          for( auto n{1U}; n_branches; n++ ) {
                auto slave_node = (*mit)->Manifold()->N(n);
                if ( slave_node->Status(var_key) == couple_if ) {
                    size_t slaveIDX = slave_node->Idx();
                    for (size_t j(0U); j < G_.Cols(); ++j )
                      if (j != slaveIDX && j != masterIDX)
                        G_.Add(masterIDX, j, G_.At(slaveIDX, j));
                    
                    // 4. adding diagonal value to master dof
                    G_.Add(masterIDX, masterIDX, G_.At(slaveIDX, slaveIDX));
                 }
        }

      // 5. copy that value from master dof to slave dof - except slave and master dofs position
      for ( auto n{1U}; n_branches; n++ ) {
            auto slave_node = (*mit)->Manifold()->N(n);
            if ( slave_node->Status(var_key) == couple_if ) {
                size_t slaveIDX = slave_node->Idx();
                for ( size_t j(0); j < G_.Cols(); ++j )
                   if ( j != slaveIDX && j != masterIDX)
                     G_.Assign( slaveIDX, j, G_.At(masterIDX, j) );
 
                // 6.  copy diagonal value from master dof to slave dof
                G_.Assign(slaveIDX, slaveIDX, G_.At(masterIDX, masterIDX) );
             }
        }
    }
    
} // end CoupleContacts




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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void  PDE_Integrator<dim,COMPUTATION_DOMAIN>::PostProcess( COMPUTATION_DOMAIN<dim>& gref )
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
              for ( typename vector<typename COMPUTATION_DOMAIN<dim>::CellType*>::const_iterator
                    git=gref.CellsBegin(); git!=gref.CellsEnd(); git++ )
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

@attention Method relies on an 0..n contiguous node numbering in the region (computational domain).

 */
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::OutputResults( COMPUTATION_DOMAIN<dim>& gref )
 {
   size_t       position;
   const uint32_t dim2(dim * dim);

    for ( operandsIterator it=basic_operands_.begin(); it!=basic_operands_.end(); it++ )
     {
        auto         gfirst(gref.NodesBegin());
        const Index  prop_key = (*it).first.key;
        if ( prop_key.place != NODE )
            throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::OutputResults(Model)",
                                           "only nodal properties can be output by this method.");
        size_t  offset = (*it).second;
       
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
                  for (auto i{0U}; i < dim; i++) {
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
                  for (auto i{0U}; i < dim; i++)
                    for ( auto k{0U}; k < dim; k++) {
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
                 for (auto i{0U}; i < prop_key.dataDepth; i++) {
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
                  for (auto i{0U}; i < prop_key.dataDepth; i++) {
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
              throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::OutputResults(Model)",
                                            "Output to ARRAY type variables is not supported by this method yet.");
             
         } // end switch(type)

 } // end for

 } // end OutputResults













/**
    Solution of the system of linear algebraic equations:
    
    1. Dimensionsing of the solution matrix (n-variables, scalar or vector etc.)
    
    2. Accumulation of the FE or FV integrals
    
    3. Multiplication of the initial conditions into the RHS in the case of the transient problem.
    
    4. Late accumulation to RHS
    
    5. Essential conditions (either with or without elimination of the Dirichlet constraints from the matrix).
    
    6. Solution
    
    7. Postprocessing (if respective pde operators were added to the PDE_Integrator).
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::IntegrateOver( COMPUTATION_DOMAIN<dim>& domain, bool debug )
 {
    // 1. configure algorithm
    EstablishMatrixSetup( domain );
    if ( !rhs_boundary_operators_.empty() )
      throw csmp::Exception( ERROR, "PDE_Integrator<>::IntegrateOver(domain):",
                            "integrator contains Boundary object integrals; call IntegrateOver(model,domain), such that boundary objects can be considered." );
   
    // setting up 0..n contiguous node numbering for index mapping
    domain.RenumberNodes();

    ReduceSystemSizeEliminatingEssentialConditions( domain );
 
    // 2. Accumulation: Note that the conditions that pertain to the group must be input !                                 
    Accumulate( domain );

    // 3. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) AssignInitialConditions( domain );

    // 4. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) LateAccumulate( domain );    

    // 5. assign conditions like Dirichlet or Neumann boundary conditions etc.
    AssignEssentialConditions( domain );
    
    // 6. diagnostics
    if ( debug ) {
         Out();
         OutputGlobals();
      }
 
    // 7. invert global matrix
    Solve();

    // 8. write results back into Model
    OutputResults( domain );
                           
    // 9. Calculation of result-dependent properties                                 
    PostProcess( domain );

 } // end IntegrateOver









/** 
    Accumulates element integrals, but
    simultaneously considering potential Boundary objects associated with the simplicial complex.
    The domain is the computational domain to which the PDE_Integrator is applied.
    
    @author SKM 7/7/2015
*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::IntegrateOver( Model<dim>& model,
                                                                        COMPUTATION_DOMAIN<dim>& domain,
                                                                        bool debug )
 {
    // 1. configure algorithm
    EstablishMatrixSetup( domain );

    // setting up 0..n contiguous node numbering for index mapping
    domain.RenumberNodes();

    ReduceSystemSizeEliminatingEssentialConditions( domain );

    // 2. Accumulation: Note that the conditions that pertain to the group must be input !                                 
    Accumulate( domain );
   
    // 3. Accumulation: potential boundary integrals from Boundary objects that share nodes with the simplicial
    //    complex of interest
    list<string> shared_boundaries;
    // enlist all boundaries: split- and regular ones
    IdentifySharedBoundaries( model, domain, shared_boundaries );
   
    // TODO: distinguish between split boundaries and boundaries  
    if ( !shared_boundaries.empty() ) {
         for ( list<string>::const_iterator
               it=shared_boundaries.begin(); it!=shared_boundaries.end(); it++ ) {
              // handling the split boundaries
              if ( (*it).find("SPLIT_BOUNDARY") != std::string::npos ) {
                   const SplitBoundary<dim>& domain_boundary = model.SplitBoundary( (*it).c_str() );
                   cout <<"\nPDE_Integrator<dim,COMPUTATION_DOMAIN>::IntegrateOver: ";
                   cout <<" accumulating split boundary: "<< (*it) <<"\n";
                   AccumulateSplitBoundaryIntegrals( domain, domain_boundary );
                }
              // handing the normal boundaries
              else {
                   const Boundary<dim>& domain_boundary = model.Boundary( (*it).c_str() );
                   cout <<"\nPDE_Integrator<dim,COMPUTATION_DOMAIN>::IntegrateOver: ";
                   cout <<" accumulating boundary: "<< (*it) <<"\n";
                   AccumulateBoundaryIntegrals( domain, domain_boundary );
                }
           }
      }

    // 4. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) AssignInitialConditions( domain );

    // 5. If the computation is transient initial conditions must be input into the righthand vector
    if ( Transient() == true ) {
         LateAccumulate( domain );
         if ( !shared_boundaries.empty() ) {
              for ( list<string>::const_iterator
                    it=shared_boundaries.begin(); it!=shared_boundaries.end(); it++ ) {
                  if ( (*it).find("SPLIT_BOUNDARY") != std::string::npos ) {
                        const SplitBoundary<dim>& domain_boundary = model.SplitBoundary( (*it).c_str() );
                        LateAccumulateSplitBoundaryIntegrals( domain, domain_boundary );
                     }
                   else {
                        const Boundary<dim>& domain_boundary = model.Boundary( (*it).c_str() );
                        LateAccumulateBoundaryIntegrals( domain, domain_boundary );
                     }
                }
           }
      }
   
    // 6. assign conditions like Dirichlet or Neumann boundary conditions etc.
    AssignEssentialConditions( domain );
    
    // 7. diagnostics
    if ( debug ) {
         Out();
         OutputGlobals();
         // OutputInput();
      }
   
// SKM_TEST checks whether the contributions are written to the correct nodes
//csmp::Index test_key = model.Database().StorageKey("test variable");
//for ( size_t i{0U}; i<domain.Nodes(); i++ )
// domain.N(i)->Store( test_key, makeScalar(PLAIN, rh_[i] ) );
 
    // 8. invert global matrix
    Solve();

    // 9. write results back into Model
    OutputResults( domain );
                           
    // 10. Calculation of result-dependent properties
    PostProcess( domain );

 } // end IntegrateOver
 
 
 
/**
    Identifies the names of any model boundaries (standard and split ones),
    which are shared with the computational domain 'subdomain'
    by searching for the corresponding string in the boundaries. The identified boundaries 
    are returned in the argument list.
    
    @attention method assumes that the names of the boundaries are composed of the strings
    BOUNDARY, the name of the volumetric boundary that is adjacent to them and their name,
    for instance,  BOUNDARY_MATRIX_LANDSURFACE. The only other names that are handled are those
    of the boundaries in a box-shaped model.
    
    @attention method also handles case, where the computational domain is the region 'Model'.
    
    @attention internal model boundaries are not handled yet.
    
    @author SKM
       
    TODO: think about meaningful PDE operators for internal model boundaries

*/
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
bool PDE_Integrator<dim,COMPUTATION_DOMAIN>::IdentifySharedBoundaries( const Model<dim>& model,
                                                                                   const COMPUTATION_DOMAIN<dim>& subdomain,
                                                                                   list<string>& shared_boundaries )
 {
    shared_boundaries.clear();

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
         return true;
      }
      
    // 3. if the computational domain is any other subregion of the model
    else {
        for ( auto it=model.BoundariesBegin(); it!=model.BoundariesEnd(); it++ )
          if ( subdomain.SharedPerimeterNodes( (*it).second.NodesBegin(), (*it).second.NodesEnd() ) > 0 )
            shared_boundaries.push_back( (*it).first );    
      }
   
    // 3. for an entire irregularly shaped model find boundaries which touch each other
    model.Region("Model").UpdateMemberIndexes();
   
    for ( typename BoundaryInterface<dim,Boundary>::boundaryConstIterator
          it=model.BoundariesBegin(); it!=model.BoundariesEnd(); it++ )
      for ( auto fit=(*it).second.CellsBegin(); fit!=(*it).second.CellsEnd(); fit++ )
        if ( subdomain.IsPerimeterCell( (*fit)->InnerParent()->Idx() ) ) {
             shared_boundaries.push_back( (*it).first );
             break;
          }
   
    // restoring a node numbering that is unique to the computational domain
    subdomain.UpdateMemberIndexes();
   
    if ( !shared_boundaries.empty() ) return true;
    return false;
   
 } // end IdentifySharedBoundaries





/**
    Enumerate method must be applied AFTER establish matrix setup process.
    and BEFORE the accumulated proecess.
    When you assemble vector or tensor variables, you also have the option
    of only assembling one of their components. To do this just set the
    components that you do not want to assemble to DBL_MAX.
 
    @attention Method relies on a 0..n contiguous numbering of the finite element nodes.
 
    @author Luat Khoa Tran

*/
template<uint32_t dim, template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim, COMPUTATION_DOMAIN>::ReduceSystemSizeEliminatingEssentialConditions( const COMPUTATION_DOMAIN<dim>& gref )
 {
    const uint32_t dim2(dim * dim);
    DOF_indexes_.resize(this->rh_.size());
    fill(DOF_indexes_.begin(), DOF_indexes_.end(), 0);

    if (!this->setup_established_)
      throw csmp::Exception(ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::ReduceSystemSizeEliminatingEssentialConditions",
        "please call EstablishMatrixSetup() prior to this method.");

    if (this->basic_operands_.empty())
      throw csmp::Exception(ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::ReduceSystemSizeEliminatingEssentialConditions",
        "No (basic) operands have been specified...");


    size_t DOF(0);

    for ( typename PDE_Integrator<dim,COMPUTATION_DOMAIN>::operandsConstIterator
          it = this->test_operands_.begin(); it != this->test_operands_.end(); it++)
        {
          auto        niter(gref.NodesBegin());
          csmp::Index prop_key = (*it).first.key;
          size_t      offset = (*it).second;

          if (prop_key.place != NODE)
            throw csmp::Exception(ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::ReduceSystemSizeEliminatingEssentialConditions",
              "So far no conditions are assigned to elements, faces, segments");

          size_t  position(0);
          switch (prop_key.type) {
              case SCALAR:
                while (niter != gref.NodesEnd()) {
                    position = (*niter)->Idx() + offset;
                    if ((*niter)->Status(prop_key) == DIRICH) DOF_indexes_[position] = NULL_IDX;
                    else {
                        DOF_indexes_[position] = DOF;
                        DOF = DOF + 1U;
                      }
                    niter++;
                  }
                break;
              case VECTOR:
                while ( niter != gref.NodesEnd()) {
                       for ( auto i{0U}; i < dim; ++i ) {
                            position = (*niter)->Idx() * dim + i + offset;
                            if ( (*niter)->Status(prop_key,i) == DIRICH ) DOF_indexes_[position] = NULL_IDX;
                            else {
                                 DOF_indexes_[position] = DOF;
                                 DOF = DOF + 1U;
                              }
                         }
                      niter++;
                   }
                break;
              case TENSOR:
                // tensors have flags only for their diagonal elements
                while ( niter != gref.NodesEnd() )
                  {
                     for (auto i{0U}; i < dim; i++) {
                        if ( (*niter)->Status(prop_key,i) == DIRICH )
                          for (size_t j{0U}; j < dim; j++) {
                               position = (*niter)->Idx() * dim2 + i * dim + j + offset;
                               DOF_indexes_[position] = NULL_IDX;
                            }
                        else for (size_t j{0U}; j < dim; j++) {
                                  position = (*niter)->Idx() * dim2 + i * dim + j + offset;
                                  DOF_indexes_[position] = DOF;
                                  DOF = DOF + 1U;
                               }

                       }
                    niter++;
                  }
                break;
              case ARRAY:
                // array variables only have a single flag
                while (niter != gref.NodesEnd()) {
                    if ( (*niter)->Status(prop_key) == DIRICH )
                      {
                        for (auto i{0U}; i < prop_key.dataDepth; i++) {
                          position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                          DOF_indexes_[position] = NULL_IDX;
                        }
                      }
                    else {
                        for (auto i{0U}; i < prop_key.dataDepth; i++) {
                          position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                          DOF_indexes_[position] = DOF;
                          DOF = DOF + 1U;
                        }
                      }// end if
                    niter++;
                  } // end while
                break;
              case FLAGGEDARRAY:
                while ( niter != gref.NodesEnd() )
                  {
                     for (auto i{0U}; i < prop_key.dataDepth; i++ )
                        if ( (*niter)->Status(prop_key,i) == DIRICH ) {
                             position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                             DOF_indexes_[position] = NULL_IDX;
                          }
                        else {
                             position = (*niter)->Idx() * prop_key.dataDepth + i + offset;
                             DOF_indexes_[position] = DOF;
                             DOF = DOF + 1U;
                          }
                     niter++;
                  }
             break;
              default:
                throw csmp::Exception(FATAL_ERROR,
                  "PDE_Integrator<dim,COMPUTATION_DOMAIN>::ReduceSystemSizeEliminatingEssentialConditions",
                  "Variable type not recognised by this method");
              }
          
      } // end for (all Dirichlet flagged variables)

    this->G_.Resize(DOF);
    this->rh_.resize(DOF);
    if ( trim_vectors_ ) vector<double>(this->rh_).swap(this->rh_);
    fill(this->rh_.begin(), this->rh_.end(), 0.);
    this->x_.resize(DOF);
    if ( trim_vectors_ ) vector<double>(this->x_).swap(this->x_);

    pivotVector_.resize(DOF);
    fill(pivotVector_.begin(), pivotVector_.end(), 0.);
    if ( trim_vectors_ ) vector<double>(this->x_).swap(this->x_);

 } // end ReduceSystemSizeEliminatingEssentialConditions








template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::WriteGlobalMatrixBitMapToText( const char* file )
 {
    char  outfile[INFO_STRING];
    strcpy( outfile, file );
    strcat( outfile, ".txt" );

     // 1. opening data output file in ascii format
     ofstream ofs;
     ofs.open( outfile, ios::out|ios::trunc );
     if ( !ofs )
       throw csmp::Exception( ERROR, "PDE_Integrator<dim,COMPUTATION_DOMAIN>::WriteGlobalMatrixBitMapToText",
                                      "Output file could not be opened" );

     // 2. writing G matrix to file
     for ( size_t i{0U}; i<G_.Rows(); i++ )
       {
          for ( size_t j{0U}; j<G_.Cols(); j++ )
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

    SparseMatrix            G_;
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
template<uint32_t dim,template<uint32_t> class COMPUTATION_DOMAIN>
void PDE_Integrator<dim,COMPUTATION_DOMAIN>::Out() const
 {
     cout <<"\nPDE_Integrator<dim,COMPUTATION_DOMAIN>::Out:\n";
   
     if ( basic_operands_.empty() || test_operands_.empty() ) {
          cout <<"\nintegrator has not been initialized yet.\n";
          return;
       }
   
     cout <<"\nsolution variables, their type and their offsets in the righthand vector:\n";
     //  std::map<Parameter,size_t>  test_operands_;
     for ( auto it=test_operands_.begin(); it!=test_operands_.end(); it++ )
       cout <<"\n"<< (*it).first.name <<" ("<< parseType((*it).first.key.type) <<"), offset: "<< (*it).second;

     cout <<"\n\nlefthand element integrals that will be accumulated (material, basic and test operands):\n";
     for ( auto it=lhs_operators_.begin(); it!=lhs_operators_.end(); it++ ) {
          cout <<"\n"<< (*it).first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }
   
     cout <<"\n\nrighthand element integrals that will be accumulated:\n";
     for ( auto it=rhs_operators_.begin(); it!=rhs_operators_.end(); it++ ) {
          cout <<"\n"<< (*it).first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }

     cout <<"\n\nsurface element integrals that will be accumulated:\n";
     for ( auto it=rhs_boundary_operators_.begin(); it!=rhs_boundary_operators_.end(); it++ ) {
          cout <<"\n"<< (*it).first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }

     cout <<"\n\npost-processing operators:\n";
     for ( auto it=postpro_operators_.begin(); it!=postpro_operators_.end(); it++ ) {
          cout <<"\n"<< (*it).first <<": ";
          //cout << (*it).second->MaterialOperand().name <<", "<< (*it).second->BasicOperand().name <<", "<< (*it).second->TestOperand().name;
       }
   
     cout <<"\n\nsystem dimensions: "<< G_.Rows() <<" x "<< G_.Cols() <<"\n";
     cout <<"\n\tdegrees of freedom (dof) per node:  "<< dof_per_node_;
     cout <<"\n\tset up etablished:                  "<< setup_established_;
     cout <<"\n\tkeep solution matrix between steps: "<< retain_matrix_;
     cout <<"\n\ttime increment:                     "<< time_increment_;
   
     if ( solver_ != NULL ) cout <<"\n\nSolver: "<< typeid(solver_).name() << endl;
   
 } // end Out





template class PDE_Integrator<1U,Region>;
template class PDE_Integrator<2U,Region>;
template class PDE_Integrator<3U,Region>;

template class PDE_Integrator<1U,Boundary>;
template class PDE_Integrator<2U,Boundary>;
template class PDE_Integrator<3U,Boundary>;

template class PDE_Integrator<1U,SplitBoundary>;
template class PDE_Integrator<2U,SplitBoundary>;
template class PDE_Integrator<3U,SplitBoundary>;

template class PDE_Integrator<1U,NimbleRegion>;
template class PDE_Integrator<2U,NimbleRegion>;
template class PDE_Integrator<3U,NimbleRegion>;

} // end namespace csmp
