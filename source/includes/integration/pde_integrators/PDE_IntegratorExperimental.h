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

@brief Accumulator and integrator for FEM or FVM-derived integral forms of PDEs
and their solution in the form of linear algebraic sets of equations (system Ax=b).
Applicable to CSMP++ Model objects or their subregions.

@author S.K. Matthai
@author Stephen G. Roberts
@author various extra constributions (Adrian Burri, Julian Mindel, Kho Luat Tran).
@date 1997

@section motivation Motivation

To encapsulate the integration and numerical solution of
partial-differential equations (PDE) in their integral form. This
process entails (1) the assembly of a sparse matrix with the contributions
from finite elements or finite volumes, (2) the assembly of a right-hand vector with
essential and initial conditions; elimination of Dirichlet constraints.
(3) Inversion of the matrix equation Ax=b to find the solution vector x.
The dependent physical variable(s) stored in
the solution vector are updated on the mesh.


@section design Design Intent

The intent of using a 'strategy' pattern for the
design was to separate the PDE operator / integrator implementation from the
the type of finite element(s) or FV's used and the solution method chosen.
As a consequence, the same PDE operator objects and solution methods can
be used for a range of computational cell types with
different node numbers and so forth.

The PDE_IntegratorExperimental class is a base class from which specific PDE_IntegratorExperimental subclasses can
be inherited to test improvements of its functionality.


@section applicability Applicability

PDE_IntegratorExperimental objects can be used to carry out any 2D and 3D computation which
are possible with CSMP++. These computations may either appy to entire
Model objects or to Regions which form subsets thereof.


@section participant Participants

Each PDE_IntegratorExperimental contains maps of pointers to left and righthand operators
that are inherited from the MathOperatorLHS and MathOperatorRHS base
classes. Other variables set the state of the algorithm, the time_increment
in transient calculations and so forth.

PDE_IntegratorExperimental also contains or connects to a Solver object which inverts the sparse solution matrix.



@section collaborations Collaborations

PDE_IntegratorExperimental interact with Model, Region, Boundary and SplitBoundary objects.
It queries these objects for the data needed to setup the solution matrices and vectors.

PDE_IntegratorExperimental also collaborates with the MeshManager and the
PropertyDatabase to gain access to variables and finite elements
(Element class instances).

PDE_IntegratorExperimental obtains the finite-element contributions from the Element class
that in turn accesses the specific finite-element type which appears hidden
behind the specific element that is assembled.
This is implemented as a 'bridge' pattern.


@section consequences Consequences

The design as a full-function base class implies that to create a PDE operator subclass one only has to
overwrite those methods that one wants to modify.
Other code does not need to be
modified. To create additional PDE operators, inherit these from MathOperatorLHS or
MathOperatorRHS.


@section implementation Implementation

The PDE_IntegratorExperimental is a stand-alone object, which is passed to the Model.
The interaction with the Model follows a visitor like pattern. The Model
passes the PDE_IntegratorExperimental on to a region to which its application has been restricted.
Inside the Model, the following steps take place when a
PDE_IntegratorExperimental is applied by calling Model::Apply():

To start with, the Model checks whether the PDE_IntegratorExperimental must be handed down to
a specific region, boundary or split boundary. Then the PDE_IntegratorExperimental builds or updates its
vectors and matrices so that they can store a sufficient number of entries:

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
Here an elimination is performed, retaining, if so, the symmetric shape of the solution matrix.
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

The simple most computation with an PDE_IntegratorExperimental is a 3D steady-state
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

*/
template<size_t dim,template<size_t> class COMPUTATION_DOMAIN>
class PDE_IntegratorExperimental {
  public:

    typedef typename std::map<Parameter,size_t>::const_iterator operandsConstIterator; ///< constant iterator over solution variables
    typedef typename std::map<Parameter,size_t>::iterator       operandsIterator;      ///< iterator over solution variables

    PDE_IntegratorExperimental();
  
    /// create an integrator with a solver managed elsewhere
    explicit PDE_IntegratorExperimental( Solver& );
    virtual ~PDE_IntegratorExperimental();
    /// do not copy because of all the pointer ramificatons
    PDE_IntegratorExperimental( const PDE_IntegratorExperimental& ) = delete;
    /// do not assign
    PDE_IntegratorExperimental& operator=( const PDE_IntegratorExperimental& ) = delete;

    /// add desired integral terms (pde operators) to the left-hand side matrix of the algebraic system of equations
    void          Add( MathOperatorLHS<dim>* );

    /// add desired integral terms (pde operators) to the right-hand side vector of the algebraic system of equations
    void          Add( MathOperatorRHS<dim>* );
  
    /// adds integral terms on the boundary of the computational domain if any
    void          AddBoundaryIntegral( MathOperatorRHS<dim>* );

    /// for the computation of fluxes across internal split boundaries
    void          AddSplitBoundaryIntegral( MathOperatorRHS<dim>* );
    void          AddSplitBoundaryIntegral( MathOperatorLHS<dim>* );

    /// adds math operators that will be applied in a second loop after the matrix has been inverted
    void          AddPostProcess( MathOperatorLHS<dim>* );

    /// sets the time-increment for- and triggers a transient calculation (use 1/t if problem has been set up that way)
    void          TimeIncrement( double dt );
  
    /// returns whether a finite-difference time increment has been set
    bool          Transient() const;
  
    /// accumulates, assembles and solves PDEs in domain of interest; @param debug prompts output of solution matrices to file; uses node numbering
    void          IntegrateOver( COMPUTATION_DOMAIN<dim>&, bool debug=false );
  
    /// simultaneously considers potential Boundary objects sharing nodes with the model subdomain on which the solution is obtained; uses node numbering
    void          IntegrateOver( Model<dim>&, COMPUTATION_DOMAIN<dim>&, bool debug=false );

    /// switch to another solver deleting any dynamically allocated solver that was associated with integrator
    void          SetSolver( Solver& );
  
    /// returns a reference to the current solver object
    Solver&       GetSolver() const;
  
    /// returns PDE integrator into the state created by default constructor ( @todo is Reset is this needed?)
    virtual void  Reset( bool delete_math_operators=true );

    /// outputs solution vector; @attention this produces meaningful results only after application of the PDE integrator
    void          SolutionVector( std::vector<double>& ) const;
  
    /// optional input of a solution vector with an initial guess of the result; for SAMG set ifirst=0 so that this vector is used
    void          FirstGuess( const std::vector<double>& );
  
    /// if an evolutionary problem where only the right-hand side changes, the matrix needs to be assembled only once
    void          RetainGlobalSolutionMatrix( bool yes_or_no );
  
    /// 'swap trick' keeps size and capacity of solution and righthand vectors equal
    void          TrimExcessCapacityOfVectors( bool trim );
  
    /// writes out the sparsity pattern of the solution matrix
    void          WriteGlobalMatrixBitMapToText( const char* file_name );
  
    void          ListMathOperatorsLHS() const;
    void          ListMathOperatorsRHS() const;
  
    /// outputs the system A x = b to the screen
    void          OutputGlobals( int precision=1 );
  
    /// prints current parameters settings and pde operators to screen
    void          Out() const;
  
    /// extra detailed screen output about the assembly  solution progress
    void          Verbose( bool );
    bool          Verbose() const;

  protected:

    /// checks whether (returns true) any Boundary object in the model is a surface of the computational domain
    bool IdentifySharedBoundaries( const Model<dim>&, const COMPUTATION_DOMAIN<dim>&, std::list<std::string>& shared_boundaries );
  
    /// for the elimination of Dirichlet constraints (Luat Khoa Tran)
    void  ReduceSystemSizeEliminatingEssentialConditions( const COMPUTATION_DOMAIN<dim>& );
  
    /// resizes sparse solution matrix and establishes variable offsets if a system of equations will be solved
    virtual void  EstablishMatrixSetup( const COMPUTATION_DOMAIN<dim>& );

    /// in time-dependent calculations this method assigns initial conditions to the RHS; uses node numbering
    virtual void  AssignInitialConditions( const COMPUTATION_DOMAIN<dim>& );

    /// eliminates Dirichlet conditions from the solution matrix and right-hand vector; uses node numbering
    virtual void  AssignEssentialConditions( const COMPUTATION_DOMAIN<dim>& );

    /// accumulates finite element integrals into solution matrix and right-hand side; uses node numbering
    virtual void  Accumulate( const COMPUTATION_DOMAIN<dim>& );

    /// accumulates surface integrals from Neumann-flagged Face object variables representing those parts of all boundaries that delimit the computational domain
    virtual void  AccumulateBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>&, const Boundary<dim>& );
    
    /// accumulation of Robin-type boundary conditions to SplitBoundary interfaces
    virtual void  AccumulateSplitBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>&, const SplitBoundary<dim>& );

    /// accumulates finite-element integrals to matrix and vector after the corresponding entries were already multiplied with the initial conditions; uses node numbering
    virtual void  LateAccumulate( const COMPUTATION_DOMAIN<dim>& );

    /// late accumulates surface integrals from Neumann-flagged Face object variables representing those parts of all boundaries that delimit the computational domain
    virtual void  LateAccumulateBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>&, const Boundary<dim>& );
    virtual void  LateAccumulateSplitBoundaryIntegrals( const COMPUTATION_DOMAIN<dim>&, const SplitBoundary<dim>& );

    /// calls connected solver object to find x in G x = rh problem
    virtual void  Solve();

    /// allows to apply pde operators to post-process the newly computed solution
    virtual void  PostProcess( const COMPUTATION_DOMAIN<dim>& );

    /// transfers the results stored in solution vector onto the nodes of the computational domain; uses node numbering
    virtual void  OutputResults( COMPUTATION_DOMAIN<dim>& );

    std::map<std::string,MathOperatorLHS<dim>*>  lhs_operators_;           ///< stencils for lefthand solution matrix
    std::map<std::string,MathOperatorRHS<dim>*>  rhs_operators_;           ///< stencils for righthand vector
    std::map<std::string,MathOperatorRHS<dim>*>  rhs_boundary_operators_;  ///< surface integrals for accumulation over boundary
    std::map<std::string,MathOperatorLHS<dim>*>  lhs_split_boundary_operators_;  ///< implicit integral coupling terms for SplitBoundary
    std::map<std::string,MathOperatorRHS<dim>*>  rhs_split_boundary_operators_;  ///< explicit integral coupling terms for SplitBoundary
    std::map<std::string,MathOperatorLHS<dim>*>  postpro_operators_;       ///< post-processing: Darcy velocities, stresses from strains etc.
    std::map<Parameter,size_t>                   basic_operands_;
    std::map<Parameter,size_t>                   test_operands_;           ///< dependent variables in the solved system of equations

    SparseMatrix            G_;                        ///< solution matrix
    std::vector<double>   rh_;                       ///< righthand vector
    std::vector<double>   x_;                        ///< solution vector
    std::vector<size_t>     DOF_indexes_;              ///< for indexing DOFs (only non-Dirichlet dofs, enumerated 0 -> maximum DOF
    std::vector<double>   pivotVector_;              ///< terms recovered from eliminated rows

    Solver*                 solver_;

    size_t                  dof_per_node_;
    bool                    setup_established_, retain_matrix_, trim_vectors_; ///< false, false, false to start with
    bool                    newed_Solver_object_;
    double                time_increment_;

    struct SIZES {
        size_t nodes;
        size_t elements;
      } target_;

  private:

    bool verbose_ = true;
};

} // csmp

#endif





