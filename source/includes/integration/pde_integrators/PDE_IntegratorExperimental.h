#ifndef PDE_IntegratorExperimental_EXPERIMENTAL_H
#define PDE_IntegratorExperimental_EXPERIMENTAL_H

#include "ErrorHandler.h"

#include "CSMP_definitions.h"
#include "LinearSolver.h"

#include "Index.h"
#include "Parameter.h"

#include "DenseMatrix.h"
#include "SparseMatrix.h"

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"

#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

#include "ModelSubDomain.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

#include "MathOperatorRHS.h"
#include "MathOperatorLHS.h"

#include <string>
#include <iostream>
#include <cassert>

namespace csmp {

template<size_t> class Model;
template<size_t> class Boundary;

/**
@brief Accumulator and integrator for FEM or FVM or mixed global solution
of linear algebraic sets of equations (system Ax=b) on a Model or its subregions.

@author S.K. Matthai
@author Stephen G. Roberts
@date 1997

@attention By default(if not explicitly specified otherwise), this uses the LUdcmp_Solver as default to circumvent 3rd party dependence

@section motivation Motivation

To encapsulate the integration and numerical solution of
partial-differential equations (PDE) in their finite-element form. This
should entail (1) the assembly of a sparse matrix with the contributions
from the finite elements and (2) the assembly of a right-hand vector with
essential and initial conditions. (3) The matrix equation Ax=b is inverted
to find the solution vector x. The dependent physical variable stored in
the solution vector must be stored back into property storage.


@section design Design Intent

The intent of the algorithm
design was to separate the PDE operator implementation from the
implementation of the type of finite element which is used. Thus, the same
PDE operator can now be used for a range of finite element types with
different node numbers and so forth.

The PDE_IntegratorExperimental class is a base class from which specific PDE_IntegratorExperimentals can
be inherited to test improvements of its functionality.


@section applicability Applicability

PDE_IntegratorExperimental objects can be used to carry out any 2D and 3D computation which
are possible with CSP. These computations may either appy to entire
Model objects or to Regions which form subsets thereof.


@section participant Participants

Each PDE_IntegratorExperimental contains maps of pointers to left and righthand operators
that are inherited from the MathOperatorLHS and MathOperatorRHS base
classes. Other variables set the state of the algorithm, the time_increment
in transient calculations and so forth.

PDE_IntegratorExperimentals contain a Solver object which inverts the sparse solution matrix.



@section collaborations Collaborations

PDE_IntegratorExperimentals interact with SuperGrous and Regions. They query these objects
for the data which they need to setup the solution matrices and vectors.


They also use references to the MeshManager, MemoryManager, and the
PropertyDatabase to gain access to the variables and the finite elements
(Element class instances).

PDE_IntegratorExperimentals obtain the finite-element contributions from the Element class
that in turn accesses the specific finite-element type which appears hidden
behind the specific element that is assembled.


@section consequences Consequences

If you want to create a new PDE operator you now only have to write a
single addition inherited class. None of the existing code has to be
modified. You do this by inheriting from MathOperatorLHS or
MathOperatorRHS.


@section implementation Implementation

PDE_IntegratorExperimentals are stand-alone objects which are passed to the Model. The
interaction with the S.G. follows then a visitor pattern. The Model
passes the PDE_IntegratorExperimental on to a group if its application has been restricted
to Regions. Inside the Model, the following steps take place when an
PDE_IntegratorExperimental is applied by calling Model::Apply():


First, the Model checks whether the PDE_IntegratorExperimental must be handed down to
a specific group. If not, the PDE_IntegratorExperimental is allowed to build or update its
vectors and matrices:

@code
PDE_IntegratorExperimental::EstablishMatrixSetup( mesh, phys_vars );
@endcode

Now the PDE_IntegratorExperimental accumulates the contributions to the global solution
matrices seqentially:

@code
PDE_IntegratorExperimental::Accumulate ( mesh, property_collection );
@endcode

If the computation is a transient one, initial conditions are added into
the righthand vector:

@code
PDE_IntegratorExperimental::AssignInitialConditions( mesh, property_collection );
@endcode

Then, the boundary conditions are applied to the lefthand and righthand
side:

@code
PDE_IntegratorExperimental::AssignEssentialConditions( mesh, property_collection );
@endcode

If the computation is a transient one, and the solution method is Backward
Euler implicit timestepping, righthand-side source terms (typically fluid
source or other rates) are now multiplied with the time increment and added
to the righthand vector.

@code
PDE_IntegratorExperimental::LateAccumulate( mesh, property_collection );
@endcode

Now the solver is invoked to invert the solution matrix:

@code
PDE_IntegratorExperimental::Solve();
@endcode

The results from the computation are mapped back into property storage:

@code
PDE_IntegratorExperimental::OutputResults( mesh, property_collection );
@endcode

and post-processing operations are applied:

@code
PDE_IntegratorExperimental::PostProcess( mesh, property_collection );
@endcode

Post-processing may be for instance, the computation of flow velocities
from computed fluid-pressure gradients and permeability values.

Important is also the mapping from global node ID numbers to entry
positions in the global solution matrix: In the simplemost case, of
a single degree-of-freedom per node and a global computation, the
mapping is:

G(n,n), rh(n) -> n = node-ID-1


For vector variables:

G(n,n), rh(n) -> n = (node-ID-1) * dimensions + j

where j is the index variable for the spatial dimensions x, y, z.

For tensor variables:

G(n,n), rh(n) -> n = (node-ID-1) * dimensions^2 i*dimensions + j

where i, j are the index variables corresponding to the rows and
column indices of the tensor variable with the i_max = j_max =
spatial dimensions.

In the case of Region computations, an instance of the class
IndexMapper transforms the global node ID numbers into the scaled
ID numbers 1...n-nodes in group. Subsequently, the same transformations
are applied for vector and tensor variables.


@section examples Application Examples


The simplemost computation with an PDE_IntegratorExperimental is a 3D steady-state
computation of fluid pressure. This involves the definition and
instantiation of some PDE operators as first step:

@code
  // compute steady-state fluid pressure from [K]{p} = {q}
  PDE_IntegratorExperimental<2U>  fluid_pressure;

  NumIntegral_dNT_op_dN_dV<2U>  conductance( example_model.Database(), "conductivity",   "fluid pressure", "fluid pressure" );
  NumIntegral_NT_op_N_dV<2U>    source( example_model.Database(), "fluid volume source", "fluid pressure" );
  VelocityAndVolumeFlux<2U>     velo( example_model, "conductivity", "porosity", "fluid pressure", true );

  fluid_pressure.Add( &conductance );
  fluid_pressure.Add( &source );
  fluid_pressure.AddPostProcess( &velo );
  example_model.Apply( fluid_pressure );
@endcode


A more advanced calculation of displacements, strains and stresses is
done in the following example:

@code
  PDE_IntegratorExperimental<2U>  deformation;
  if ( restricted_to_rock ) deformation.RestrictApplicationTo("rock");

  PT_op<2U>                    bforces( p_ref, "force", "displacement" );
  NumIntegral_BT_D_B_dV<2U>    stiffness( p_ref, "Young's modulus", "Poisson's ratio", "displacement", "displacement");
  if ( with_plane_stress ) stiffness.PlaneStress();
  NumIntegral_PT_op_dS<2U>     bstresses( p_ref, "Neumann stress", "displacement" );
  NumIntegral_PT_op_dV<2U>     bodyforce( p_ref, "gravity force", "displacement");
  NumIntegral_BT_D_op_dV<2U>   volstrain( p_ref, "dilatation", "Young's modulus", "Poisson's ratio", "displacement");
  NumIntegral_BT_op_dV<2U>     porepressure( p_ref, "fluid pressure", "displacement");

  deformation.Add( &stiffness );
  deformation.Add( &bforces );  // must always be there so that Dirichlet conditions are accumulated
  if ( with_boundary_stresses ) deformation.Add( &bstresses );
  if ( with_body_forces )       deformation.Add( &bodyforce );
  if ( with_volume_strains )    deformation.Add( &volstrain );
  if ( with_pore_pressure )     deformation.Add( &porepressure );

  StressesAndStrains<2U>  postpro( example_model, "Young's modulus",
                                                           "Poisson's ratio", "displacement", true, true );

  if ( with_plane_stress ) postpro.PlaneStress();
   deformation.AddPostProcess( &postpro );

  example_model.Apply( deformation );
@endcode

@section outlook Future Implementations

@todo !!! SKM: Implement the automatic integration over boundaries in the case where surface integrals are present
      (design approved: Garmisch and Colleoli)

*/
template<size_t dim,template<size_t> class SIMPLICIAL_COMPLEX>
class PDE_IntegratorExperimental {

  public:

    typedef typename std::map<Parameter,size_t>::const_iterator operandsConstIterator; ///< constant iterator over solution variables
    typedef typename std::map<Parameter,size_t>::iterator       operandsIterator;      ///< iterator over solution variables

    PDE_IntegratorExperimental();
  
    /// create an integrator with a solver managed elsewhere
    explicit PDE_IntegratorExperimental( Solver& );
    virtual ~PDE_IntegratorExperimental();

  protected:

    /// make a copy of this PDE integrator with its own unique solver
    PDE_IntegratorExperimental( const PDE_IntegratorExperimental& );
  
    /// assign setup to another PDE integrator
    PDE_IntegratorExperimental& operator=( const PDE_IntegratorExperimental& );

  public:

    /// add desired integral terms (pde operators) to the left-hand side matrix of the algebraic system of equations
    void          Add( MathOperatorLHS<dim>* );

    /// add desired integral terms (pde operators) to the right-hand side vector of the algebraic system of equations
    void          Add( MathOperatorRHS<dim>* );
  
    /// adds integral terms on the boundary of the computational domain if any
    void          AddBoundaryIntegrals( MathOperatorRHS<dim>* );
  
    /// adds math operators that will be applied in a second loop after the matrix has been inverted
    void          AddPostProcess( MathOperatorLHS<dim>* );

    /// sets the time-increment for- and triggers a transient calculation (use 1/t if problem has been set up that way)
    void          TimeIncrement( double64 dt );
  
    /// returns whether a finite-difference time increment has been set
    bool          Transient() const;
  
    /// accumulates, assembles and solves PDEs in domain of interest; @param debug prompts output of solution matrices to file
    void          IntegrateOver( SIMPLICIAL_COMPLEX<dim>&, bool debug=false );
  
    /// performs an elimination of the Dirichlet degrees-of-freedom prior to solving the problem
    void          IntegrateOver1( SIMPLICIAL_COMPLEX<dim>&, bool debug=false );

    /// simultaneously considers potential Boundary objects sharing nodes with the simplicial complex on which the solution is obtained
    void          IntegrateOver( Model<dim>&, SIMPLICIAL_COMPLEX<dim>&, bool debug=false );

    /// switch to another solver deleting any dynamically allocated solver that was associated with integrator
    void          SetSolver( Solver& );
  
    /// returns a reference to the current solver object
    Solver&       GetSolver() const;
  
    /// applies scale factor to Dirichlet matrix-diagonal entries as applied by AssignEssentialConditions() and the rhs entries
    void          ScaleEssentialConditions( double64 scale_factor);
  
    /// returns PDE integrator into the state created by default constructor ( @todo is Reset is this needed?)
    virtual void  Reset( bool delete_math_operators=true );

    void          ListMathOperatorsLHS() const;
    void          ListMathOperatorsRHS() const;
  
    /// outputs solution vector; @attention this is useful only after an application of the PDE integrator
    void          SolutionVector( std::vector<double64>& ) const;
  
    /// optional input of a solution vector with an initial guess of the result; for SAMG set ifirst=0 so that this vector is used
    void          FirstGuess( const std::vector<double64>& );
  
    /// if an evolutionary problem where only the right-hand side changes, the matrix needs to be assembled only once
    void          RetainGlobalSolutionMatrix( bool yes_or_no );
  
    /// writes out the sparsity pattern of the solution matrix
    void          WriteGlobalMatrixBitMapToText( const char* file_name );
  
    /// outputs the system A x = b to the screen
    void          OutputGlobals( int precision=1 );
  
    /// prints current parameters settings and pde operators to screen
    void          Out() const;
  
    /// extra detailed screen output about the assembly  solution progress
    void          Verbose( bool );
    bool          Verbose() const;

  protected:

    /// checks whether (returns true) any Boundary object in the model is a surface of the computational domain
    bool IdentifySharedBoundaries( const Model<dim>&, const SIMPLICIAL_COMPLEX<dim>&, std::list<std::string>& shared_boundaries );
  
    /// resizes sparse solution matrix and establishes variable offsets if a system of equations will be solved
    virtual void  EstablishMatrixSetup( const SIMPLICIAL_COMPLEX<dim>& );

    /// if this is a time-dependent calculation, this method assigns initial conditions to the RHS
    virtual void  AssignInitialConditions( const SIMPLICIAL_COMPLEX<dim>& );

    /// zeroes out Dirichlet matrix rows, puts 1's into its diagonal, and overwrites RHS with condition value
    virtual void  AssignEssentialConditions( const SIMPLICIAL_COMPLEX<dim>& );

    /// eliminates essential (Dirichlet) conditions, condensing the the solution matrix, rhs etc. to that of the remaining DOF
    virtual void  EliminateEssentialConditions( const SIMPLICIAL_COMPLEX<dim>& );

    /// accumulates the finite element integrals into the solution matrix and right-hand side
    virtual void  Accumulate( const SIMPLICIAL_COMPLEX<dim>& );

    /// accumulates surface integrals from those parts of the supplied boundary that are shared with the computational domain
    virtual void  AccumulateBoundaryIntegrals( const SIMPLICIAL_COMPLEX<dim>&, const Boundary<dim>& );

    /// permits to add finite-element integrals to matrix and vector after terms were multiplied into them
    virtual void  LateAccumulate( const SIMPLICIAL_COMPLEX<dim>& );

    /// late accumulates surface integrals from those parts of the supplied boundary that are shared with the computational domain
    virtual void  LateAccumulateBoundaryIntegrals( const SIMPLICIAL_COMPLEX<dim>&, const Boundary<dim>& );

    /// calls connected solver object to find x in A x = b problem
    virtual void  Solve();

    /// allows to apply pde operators to post-process the solution
    virtual void  PostProcess( const SIMPLICIAL_COMPLEX<dim>& );

    /// transfers the results stored in solution vector onto the nodes of the computational domain
    virtual void  OutputResults( SIMPLICIAL_COMPLEX<dim>& );
  
    // SKM FIX output method that takes into account the Dirichlet_index_mapping_ from the reduced solution matrix to the model
    virtual void  OutputResults1( SIMPLICIAL_COMPLEX<dim>& );

  protected:

    std::map<std::string,MathOperatorLHS<dim>*>  lhs_operators_;           ///< pde operators for solution matrix
    std::map<std::string,MathOperatorRHS<dim>*>  rhs_operators_;           ///< pde operators for righthand vector
    std::map<std::string,MathOperatorRHS<dim>*>  rhs_boundary_operators_;  ///< potential surface integrals for accumulation over boundary
    std::map<std::string,MathOperatorLHS<dim>*>  postpro_operators_;       ///< post-processing: Darcy velocities, stresses from strains etc.
    std::map<Parameter,size_t>                   basic_operands_;
    std::map<Parameter,size_t>                   test_operands_;           ///< dependent variables in the solved system of equations

    SparseMatrix            G_;                        ///< solution matrix
    std::vector<double64>   rh_;                       ///< righthand vector
    std::vector<double64>   x_;                        ///< solution vector
    std::vector<long64>     Dirichlet_index_mapping_;  ///< mapping from the unknowns in the condensed solution matrix back to the indexed mode
    Solver*                 solver_;

    size_t                  dof_per_node_;
    bool                    setup_established_, retain_matrix_;
    const bool              newed_Solver_object_;
    double64                time_increment_;

    struct SIZES {
        size_t nodes;
        size_t elements;
    } target_;

  private:

    double64                scale_factor_; ///< for essential conditions
    bool                    verbose_;
};

} // csmp

#endif





