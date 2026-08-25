#ifndef CSMP_PDE_INTEGRATOR_H
#define CSMP_PDE_INTEGRATOR_H

#include "MathOperatorRHS.h"
#include "MathOperatorLHS.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Face;
template<uint32_t> class InterFace;
template<uint32_t> class NimbleRegion;
template<uint32_t> class Model;
template<uint32_t,template<uint32_t> class> class ModelSubDomain;
template<uint32_t> class Region;
template<uint32_t> class Boundary;
template<uint32_t> class SplitBoundary;
class Solver;

class PDE_Integrator_Test;
class PDE_Integrator_Transient_Test;
class PDE_Integrator_Computation_Test;

/**

@brief Accumulator and integrator for FEM or FVM-derived integral forms of PDEs
and their solution by solving linear sets of linear algebraic equations (system Ax=b).
Applicable to CSMP++ Model objects or their subregions.

@author S.K. Matthai
@author Stephen G. Roberts
@author various extra contributions (Adrian Burri, Julian Mindel, Kho Luat Tran).
@date 1997
@date replaced by merge of PDE_Integrator 1/4/_22
@date working toward a rigorous test of the functionality and better user feedback

@section motivation Motivation

To encapsulate the integration / numerical solution of
partial-differential equations (PDE). This
process entails

(0) the creation of a sparsity pattern
(1) the assembly of a sparse solution matrix, G, using integral contributions from finite elements or finite volumes,
(2) the assembly of a right-hand vector with essential and initial conditions; before, Dirichlet conditions are eliminated.
(3) Inversion of the matrix equation Ax=b to find the solution vector x.
(4) Mapping the results from the result vector 'x' back onto the mesh, performing range checks.


@section design Design

The 'strategy' design pattern is used to separate the PDE operator / integrator implementation from the
the type of finite element(s) or FV's used and the solution method chosen.
As a consequence, the same PDE operator objects and solution methods can
be used for a range of computational cell types with
different node numbers and so forth.
The Bridge pattern is used to decouple specific cell types and interpolation orders from the
way the contributions from the finite elements are integrated. These implementations
are found in collections of subclasses of the MathOperatorLHS and MathOperatorRHS,
which are grouped together in specific PDE operator libraries.
These are further subdivided into ones for analytic and others for numeric integration.

The PDE_Integrator class is a base class from which specific PDE_Integrator subclasses can
be inherited. These subclasses  are found in the pde_integrators directory in the source.
Subclasses are implemented to support nonlinear iteration loops
and more user-friendly high-level forms of integrators like the SteadyStateDiffusor.


@section applicability Applicability

PDE_Integrator supports the solution of elliptic, parabolic and even hyperbolic PDEs in 2D and 3D .
These computations can be targeted on the entire Model, Regions thereof and Boundary objects..


@section participant Participants

Each PDE_Integrator contains maps of pointers to left and righthand operators
that are inherited from the MathOperatorLHS and MathOperatorRHS base
classes. Other variables set the state of the algorithm, the time_increment
in transient calculations and so forth.

To solve the system of equations, the PDE_Integrator uses a Solver object.
Different solvers can be configured, applied and managed via the Strategy pattern
as described in the users guide.


@section collaborations Collaborations

PDE_Integrator interact with Model, Region, Boundary and SplitBoundary objects.
It queries these objects for the data needed to setup the solution matrices and vectors.
It also uses them to managed evolving boundary conditions, including ones
needed to dynamically couple domains together.

PDE_Integrator also collaborates with the MeshManager and the
PropertyDatabase to gain access to variables and finite elements
(Element class instances).

PDE_Integrator obtains the finite-element contributions from the Element class
that in turn has access to specific finite-element functionality via its FE policy.
The interior workings of this remain hidden behind the interface of the specific element that is assembled.
This is implemented as based on a 'bridge' pattern.


@section consequences Consequences

The design as a full-function base class implies that to create a PDE operator subclass one only has to
overwrite those methods that one wants to modify.
Other code does not need to be
modified. To create additional PDE operators, inherit these from MathOperatorLHS or
MathOperatorRHS.


@section implementation Implementation

The PDE_Integrator is a stand-alone object, which is passed to the Model or applied
to one of its regions or boundaries with the IntegrateOver() method.
The interaction with the Model follows a visitor-like pattern. The Model
passes the PDE_Integrator on to a region to which its application has been restricted.
Inside the Model, the following steps take place when a
PDE_Integrator is applied by calling Model::Apply():

To start with, the Model checks whether the PDE_Integrator must be handed down to
a specific region, boundary or split boundary. Then the PDE_Integrator builds or updates its
vectors and matrices so that they can store a sufficient number of entries:

@code
PDE_Integrator::EstablishMatrixSetup( mesh, phys_vars );
@endcode

Now the PDE_Integrator accumulates the contributions to the global solution
matrix seqentially:

@code
PDE_Integrator::Accumulate ( mesh, property_collection );
@endcode

If the computation is a transient one, initial conditions are added into
the righthand vector:

@code
PDE_Integrator::AssignInitialConditions( mesh, property_collection );
@endcode

Then, the boundary conditions are applied to the lefthand and righthand
side:

@code
PDE_Integrator::AssignEssentialConditions( mesh, property_collection );
Here an elimination is performed, retaining, if so, the symmetric shape of the solution matrix.
@endcode

If the computation is a transient one, and the time-stepping solution method is Backward
Euler implicit, then righthand-side source terms (typically fluid
source or other rates) are added to to the righthand vector after initial conditions were multiplied into it.
For such integrals LateAccumulate needs to be specified.

@code
PDE_Integrator::LateAccumulate( mesh, property_collection );
@endcode

Now the solver is invoked to solve the ensuing linear algebraic system  LHS * x = RHS.

@code
PDE_Integrator::Solve();
@endcode

The results from the computation are mapped back into property storage:

@code
PDE_Integrator::OutputResults( mesh, property_collection );
@endcode

and post-processing operations are applied:

@code
PDE_Integrator::PostProcess( mesh, property_collection );
@endcode

Post-processing may be, for instance, the computation of flow velocities
from computed fluid-pressure gradients and permeability values.

Important also is the mapping from global node ID numbers to entry
positions in the global solution matrix troughout the solution process:
In the simplest case, of a single degree-of-freedom per node and a global computation, the
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

The simple most computation with an PDE_Integrator is a 3D steady-state
computation of fluid pressure. This involves the definition and
instantiation of some PDE operators as first step:

@code
  // compute steady-state fluid pressure from [K]{p} = {q}
  PDE_Integrator<2U>  fluid_pressure;

  NumIntegral_dNT_lhsop_dN_dV<2U>  conductance( example_model.Database(), "conductivity",   "fluid pressure", "fluid pressure" );
  NumIntegral_NT_rhsop_N_dV<2U> source( example_model.Database(), "fluid volume source", "fluid pressure" );
  VelocityAndVolumeFlux<2U>     velo( example_model, "conductivity", "porosity", "fluid pressure", true );

  fluid_pressure.Add( &conductance );
  fluid_pressure.Add( &source );
  fluid_pressure.AddPostProcess( &velo );
  example_model.Apply( fluid_pressure );
@endcode


A more advanced calculation of displacements, strains and stresses is
done in the following example:

@code
  PDE_Integrator<2U>  deformation;
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
template<uint32_t dim, template<uint32_t> class CELLTYPE=Element, class MATRIXTYPE=CompressedRowMatrix>
class PDE_Integrator {
  public:
    typedef typename std::map<Parameter,size_t>::const_iterator operandsConstIterator; ///< constant iterator over solution variables
    typedef typename std::map<Parameter,size_t>::iterator       operandsIterator;      ///< iterator over solution variables

    PDE_Integrator();
  
    /// create an integrator with a solver managed elsewhere
    explicit PDE_Integrator( Solver& );
    virtual ~PDE_Integrator();
    /// do not copy because of all the pointer ramificatons
    PDE_Integrator( const PDE_Integrator& ) = delete;
    /// do not assign
    PDE_Integrator& operator=( const PDE_Integrator& ) = delete;

    /// add desired integral terms (pde operators) to the left-hand side matrix of the algebraic system of equations
    void          Add( MathOperatorLHS<dim,CELLTYPE>* );

    /// add desired integral terms (pde operators) to the right-hand side vector of the algebraic system of equations
    void          Add( MathOperatorRHS<dim,CELLTYPE>* );
  
    /// adds math operators that will be applied in a second loop after the matrix has been inverted
    void          AddPostProcess( MathOperatorLHS<dim,CELLTYPE>* );

    /// adds integral terms on the boundary of the computational domain if any
    void          AddBoundaryIntegral( MathOperatorRHS<dim,Face>* );
    
    /// to perform checks during the application of the PDE_Integrator
    bool          HasBoundaryIntegrals() const { return rhs_boundary_operators_.empty(); }

    /// for the computation of fluxes across internal split boundaries
    void          AddSplitBoundaryIntegral( MathOperatorRHS<dim,InterFace>* );
    void          AddSplitBoundaryIntegral( MathOperatorLHS<dim,InterFace>* );

    /// sets the time-increment for- and triggers a transient calculation (use 1/t if problem has been set up that way)
    void          TimeIncrement( double dt );
  
    /// returns whether a finite-difference time increment has been set
    bool          Transient() const;
  
    /// accumulates, assembles, and solves PDEs in domain of interest; @param debug prompts output of solution matrices to file; uses node numbering
    void          IntegrateOver( ModelSubDomain<dim,CELLTYPE>&, bool debug=false );

    /// also considers  "dS" pde operators from Boundary or SplitBoundary objects if these share nodes with domain on which the solution is obtained
    void          IntegrateOver( Model<dim>&, ModelSubDomain<dim,CELLTYPE>&, bool debug=false );

    /// switch to another solver deleting any dynamically allocated solver that was associated with integrator
    void          SetSolver( Solver& );
  
    /// returns a reference to the current solver object
    Solver&       GetSolver() const;
  
    /// returns PDE integrator into the state created by default constructor
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

    /// returns the equation number in the reduced (Dirichlet DOF-eliminated) global system or size_t::max if 'global_dof' is none because it is a Dirichlet constraint
    size_t EquationNumber(size_t global_dof ) const noexcept { return DOF_indexes_[global_dof]; }
    
    void ExpandNodeDOFs( const Node<dim>* node, const csmp::Index& top_key, size_t test_operand_offset, std::vector<size_t>& out );

    /// checks whether any Boundary object in the model is a surface of the computational domain
    std::list<std::string> IdentifySharedBoundaries( const Model<dim>&, const ModelSubDomain<dim,CELLTYPE>& ) const;

    /// checks whether any SplitBoundary object in the model is a surface of the computational domain
    std::list<std::string> IdentifySharedSplitBoundaries( const Model<dim>&, const ModelSubDomain<dim,CELLTYPE>& ) const;
  
    /// resizes sparse solution matrix and establishes variable offsets if a system of equations will be solved, but without boundary or splitboundary integrals
    virtual bool EstablishMatrixSetup( const ModelSubDomain<dim,CELLTYPE>& );

    /// for the elimination of Dirichlet constraints from the solution matrix; called after EstablishMatrixSetup but before accumulation
    void  EliminateEssentialConditions( const ModelSubDomain<dim,CELLTYPE>&, size_t max_offset );
  
    /// in time-dependent calculations this method assigns initial conditions to the RHS; uses node numbering
    virtual void  AssignInitialConditions( const ModelSubDomain<dim,CELLTYPE>& );

    /// modifies right-hand vector; adding Dirichlet condition terms that were elimitated before
    virtual void  AssignEssentialConditions( const ModelSubDomain<dim,CELLTYPE>& );

    /// accumulates finite element integrals into solution matrix and right-hand side; uses node numbering
    virtual void  Accumulate( const ModelSubDomain<dim,CELLTYPE>& );

    /// accumulates finite-element integrals to matrix and vector after the corresponding entries were already multiplied with the initial conditions; uses node numbering
    virtual void  LateAccumulate( const ModelSubDomain<dim,CELLTYPE>& );

    /// accumulates surface integrals from Neumann-flagged Face object variables representing those parts of all boundaries that delimit the computational domain
    virtual void  AccumulateBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>&, const Boundary<dim>& );
    /// version that uses 'boundary_faces_' to accumulate over
    virtual void  AccumulateBoundaryIntegrals();
    
    /// accumulation of Robin-type boundary conditions to SplitBoundary interfaces
    virtual void  AccumulateSplitBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>&, const SplitBoundary<dim>& );
    /// version that uses 'splitboundary_interfaces_' to accumulate over
    virtual void  AccumulateSplitBoundaryIntegrals();

    /// late accumulates surface integrals from Neumann-flagged Face object variables representing those parts of all boundaries that delimit the computational domain
    virtual void  LateAccumulateBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>&, const Boundary<dim>& );
    virtual void  LateAccumulateBoundaryIntegrals();
    virtual void  LateAccumulateSplitBoundaryIntegrals( const ModelSubDomain<dim,CELLTYPE>&, const SplitBoundary<dim>& );
    virtual void  LateAccumulateSplitBoundaryIntegrals();
    
    /// couples domains separated by SplitBoundaries using the information from NodeManifolds @todo:  refactor this method from Luat to increase efficiency
    void CoupleDomainsAcrossSplitBoundary( ModelSubDomain<dim,CELLTYPE>& );

    /// calls connected solver object to find x in G x = rh problem
    virtual void  Solve();

    /// allows to apply pde operators to post-process the newly computed solution
    virtual void  PostProcess( ModelSubDomain<dim,CELLTYPE>& );

    /// transfers the results stored in solution vector onto the nodes of the computational domain; uses node numbering
    virtual void  OutputResults( ModelSubDomain<dim,CELLTYPE>& );

    std::map<std::string,MathOperatorLHS<dim,CELLTYPE>*>  lhs_operators_;           ///< stencils for lefthand solution matrix
    std::map<std::string,MathOperatorRHS<dim,CELLTYPE>*>  rhs_operators_;           ///< stencils for righthand vector
    std::map<std::string,MathOperatorLHS<dim,CELLTYPE>*>  postpro_operators_;       ///< post-processing: Darcy velocities, stresses from strains etc.
    std::map<std::string,MathOperatorRHS<dim,Face>*>      rhs_boundary_operators_;  ///< surface integrals for accumulation over boundary
    std::map<std::string,MathOperatorLHS<dim,InterFace>*> lhs_split_boundary_operators_;  ///< implicit integral coupling terms for SplitBoundary
    std::map<std::string,MathOperatorRHS<dim,InterFace>*> rhs_split_boundary_operators_;  ///< explicit integral coupling terms for SplitBoundary
    std::map<csmp::Parameter,size_t>                      basic_operands_;          ///< column operators ordered alphabetical by corresponding primary variable name
    std::map<csmp::Parameter,size_t>                      test_operands_;           ///< row operators ordered alphabetical by corresponding primary variable name
  
    // TODO: recheck the circumstances under which these extra pointers must be stored
    std::vector<const csmp::Face<dim>*>       boundary_faces_;           ///< empty if computation applies to entire model; else Faces needed for dS integrals
    std::vector<const csmp::InterFace<dim>*>  splitboundary_interfaces_; ///< empty if computation applies to entire model; else InterFaces needed for domain coupling

    MATRIXTYPE            G_;           ///< solution matrix
    std::vector<double>   rh_;          ///< righthand vector
    std::vector<double>   x_;           ///< solution vector
    std::vector<size_t>   DOF_indexes_; ///< indices of DOFs, but only of the non-Dirichlet dofs, size enumerated 0 - DOF-1 (including Dirich DOF)
    std::vector<double>   pivotVector_; ///< full-system DOF (including Dirich); accumulates products of eliminated Dirichlet rows and RHS DIrich entries

    Solver*               solver_ = nullptr;

    size_t                dof_per_node_;
    bool                  setup_established_, retain_matrix_, trim_vectors_; ///< false, false, false to start with
    bool                  newed_Solver_object_;
    double                time_increment_;

    struct SIZES {
        size_t nodes;
        size_t elements;
      } target_;

  private:

    bool verbose_ = true;
    
    friend PDE_Integrator_Test;
    friend PDE_Integrator_Transient_Test;
    friend PDE_Integrator_Computation_Test;
};


/// collects pointers to Faces that may be needed to accumulate surface integrals on a computational domain consisting of elements
template<uint32_t dim>
bool collectBorderFacesOfComputationRegion( const Model<dim>&, const ModelSubDomain<dim,Element>&,
                                            std::vector<const csmp::Face<dim>*>& boundary_faces );
    
/// collects pointers to  InterFaces that may be needed to accumulate splitboundary integrals on a computational domain consisting of elements
template<uint32_t dim>
bool collectInterfacesInComputationRegion( const Model<dim>&, const ModelSubDomain<dim,Element>&,
                                           std::vector<const csmp::InterFace<dim>*>& splitboundary_interfaces );

// INLINE FUNCTIONS

/// expands node indices into degrees of freedom of the solution variable in the global system
template<uint32_t dim, template<uint32_t> class CELLTYPE, class MATRIXTYPE>
inline void PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::ExpandNodeDOFs( const Node<dim>* node,
                                                                     const csmp::Index& top_key, ///< test operand key
                                                                     size_t offset,
                                                                     std::vector<size_t>& expandedDOF_vec )
{
    const size_t base = node->Idx();

    auto push_eq = [&](size_t global_dof) {
        const size_t eq = DOF_indexes_[global_dof];
        if (eq != NULL_IDX)
            expandedDOF_vec.push_back(eq);
    };

    switch (top_key.type) {
        case SCALAR:
             push_eq(base + offset);
          break;
        case VECTOR:
             for (uint32_t i = 0; i < dim; ++i)
                push_eq(base * dim + i + offset);
          break;
        case TENSOR: {
             constexpr uint32_t dim2 = dim * dim;
             for (uint32_t i = 0; i < dim; ++i)
                for (uint32_t j = 0; j < dim; ++j)
                    push_eq(base * dim2 + i * dim + j + offset);
          break;
        }
        case ARRAY:
        case FLAGGEDARRAY:
             for ( uint32_t i = 0; i < top_key.dataDepth; ++i )
               push_eq(base * top_key.dataDepth + i + offset);
          break;
        default:
            throw std::runtime_error("PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::ExpandNodeDOFs: Unsupported variable type");
    }
}



} // csmp

#endif /* CSMP_PDE_INTEGRATOR_H */

