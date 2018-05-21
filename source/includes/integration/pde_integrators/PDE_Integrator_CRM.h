#ifndef PDE_INTEGRATOR_CRM_H
#define PDE_INTEGRATOR_CRM_H

#include "ErrorHandler.h"

#include "CSMP_definitions.h"
#include "LinearSolver.h"

#include "Index.h"
#include "Parameter.h"

#include "DenseMatrix.h"
#include "SparseMatrix.h"
#include "CompressedRowMatrix.h"

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

  template<size_t dim, template<size_t> class SIMPLICIAL_COMPLEX>
  class PDE_Integrator_CRM {
  public:

    typedef typename std::map<Parameter, size_t>::const_iterator operandsConstIterator;
    typedef typename std::map<Parameter, size_t>::iterator       operandsIterator;

    PDE_Integrator_CRM();
    explicit PDE_Integrator_CRM(Solver&);

    /// accepts a pointer to a Solver object that is managed somewhere else
    explicit PDE_Integrator_CRM(Solver*);
    virtual ~PDE_Integrator_CRM();

  protected:

    PDE_Integrator_CRM(const PDE_Integrator_CRM&);
    PDE_Integrator_CRM& operator=(const PDE_Integrator_CRM&);

  public:

    void          Add(MathOperatorLHS<dim>*);
    void          Add(MathOperatorRHS<dim>*);
    void          AddPostProcess(MathOperatorLHS<dim>*);
    void      AddBoundaryIntegrals(MathOperatorRHS<dim>*);

    void          TimeIncrement(double64 dt);
    bool          Transient() const;

    void      IntegrateOver(Model<dim>& model, SIMPLICIAL_COMPLEX<dim>& domain, bool debug = false);

    void          IntegrateOver(SIMPLICIAL_COMPLEX<dim>&, bool debug = false);

    bool        IdentifySharedBoundaries(const Model<dim>& model, const SIMPLICIAL_COMPLEX<dim>& subdomain, std::list<std::string>& shared_boundaries);

    void      AccumulateBoundaryIntegrals(const SIMPLICIAL_COMPLEX<dim>& comp_domain, const Boundary<dim>& boundary);

    void      LateAccumulateBoundaryIntegrals(const SIMPLICIAL_COMPLEX<dim>& comp_domain, const Boundary<dim>& boundary);

    void          SetSolver(Solver* new_solver);
    Solver*       GetSolver() const;
    virtual void  AdjustSolverSettings();
    void      outputCRMfile(std::string& );
    /// applies scale factor to Dirichlet matrix-diagonal entries as applied by AssignEssentialConditions() and the rhs entries
    void          ScaleEssentialConditions(double64 scale_factor);

    virtual void  Reset(bool delete_math_operators = true);

    void          ListMathOperatorsLHS() const;
    void          ListMathOperatorsRHS() const;
    void          SolutionVector(std::vector<double64>&) const;
    void          FirstGuess(const std::vector<double64>&);
    void          RetainGlobalSolutionMatrix(bool yes_or_no);
    void          WriteGlobalMatrixBitMapToText(const char* file_name);
    void          OutputGlobals(int precision = 1);
    void          Out() const { std::cout << "PDE_Integrator_CRM 2 \n"; std::getchar(); };
    void          Verbose(bool verbose) { verbose_ = verbose; }
    bool          GetVerbose() { return verbose_; }

  protected:

    virtual void  EstablishMatrixSetup(const SIMPLICIAL_COMPLEX<dim>&);

    virtual void  AssignInitialConditions(const SIMPLICIAL_COMPLEX<dim>&);

    /// zeroes out Dirichlet matrix rows, puts 1's into its diagonal, and overwrites RHS with condition value
    virtual void  AssignEssentialConditions(const SIMPLICIAL_COMPLEX<dim>&);

    virtual void  Accumulate(const SIMPLICIAL_COMPLEX<dim>&);

    virtual void  LateAccumulate(const SIMPLICIAL_COMPLEX<dim>&);

    virtual void  Solve();

    virtual void  PostProcess(const SIMPLICIAL_COMPLEX<dim>&);

    virtual void  OutputResults(SIMPLICIAL_COMPLEX<dim>&);

  protected:

    std::map<std::string, MathOperatorLHS<dim>*>  lhs_operators_;
    std::map<std::string, MathOperatorRHS<dim>*>  rhs_operators_;
    std::map<std::string, MathOperatorRHS<dim>*>  rhs_boundary_operators_;  ///< potential surface integrals for accumulation over boundary
    std::map<std::string, MathOperatorLHS<dim>*>  postpro_operators_;
    std::map<Parameter, size_t>                   basic_operands_;
    std::map<Parameter, size_t>                   test_operands_;

    CompressedRowMatrix     G_;
    std::vector<double64>   rh_;
    std::vector<double64>   RH_;
    std::vector<double64>   x_;
    Solver*                 solver_;
    std::vector<Entry>      dirich_;

    const size_t            dim2_;
    size_t                  dof_per_node_;
    bool                    setup_established_, retain_matrix_;
    bool                    newed_Solver_object;
    double64                time_increment_;

    struct SIZES {
      size_t nodes;
      size_t elements;
    } target_;

  private:

    double64                scale_factor_; ///< for essential conditions
    bool                    verbose_;
  };


  /**
  @class PDE_Integrator PDE_Integrator "main_library/PDE_Integrator.h"

  @author S.K. Matthaei
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

  The PDE_Integrator class is a base class from which specific PDE_Integrators can
  be inherited to test improvements of its functionality.


  @section applicability Applicability

  PDE_Integrator objects can be used to carry out any 2D and 3D computation which
  are possible with CSP. These computations may either appy to entire
  Model objects or to Regions which form subsets thereof.


  @section participant Participants

  Each PDE_Integrator contains maps of pointers to left and righthand operators
  that are inherited from the MathOperatorLHS and MathOperatorRHS base
  classes. Other variables set the state of the algorithm, the time_increment
  in transient calculations and so forth.

  PDE_Integrators contain a Solver object which inverts the sparse solution matrix.



  @section collaborations Collaborations

  PDE_Integrators interact with SuperGrous and Regions. They query these objects
  for the data which they need to setup the solution matrices and vectors.


  They also use references to the MeshManager, MemoryManager, and the
  PropertyDatabase to gain access to the variables and the finite elements
  (Element class instances).

  PDE_Integrators obtain the finite-element contributions from the Element class
  that in turn accesses the specific finite-element type which appears hidden
  behind the specific element that is assembled.


  @section consequences Consequences

  If you want to create a new PDE operator you now only have to write a
  single addition inherited class. None of the existing code has to be
  modified. You do this by inheriting from MathOperatorLHS or
  MathOperatorRHS.



  @section implementation Implementation

  PDE_Integrators are stand-alone objects which are passed to the Model. The
  interaction with the S.G. follows then a visitor pattern. The Model
  passes the PDE_Integrator on to a group if its application has been restricted
  to Regions. Inside the Model, the following steps take place when an
  PDE_Integrator is applied by calling Model::Apply():


  First, the Model checks whether the PDE_Integrator must be handed down to
  a specific group. If not, the PDE_Integrator is allowed to build or update its
  vectors and matrices:

  @code
  PDE_Integrator::EstablishMatrixSetup( mesh, phys_vars );
  @endcode

  Now the PDE_Integrator accumulates the contributions to the global solution
  matrices seqentially:

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
  @endcode

  If the computation is a transient one, and the solution method is Backward
  Euler implicit timestepping, righthand-side source terms (typically fluid
  source or other rates) are now multiplied with the time increment and added
  to the righthand vector.

  @code
  PDE_Integrator::LateAccumulate( mesh, property_collection );
  @endcode

  Now the solver is invoked to invert the solution matrix:

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


  The simplemost computation with an PDE_Integrator is a 3D steady-state
  computation of fluid pressure. This involves the definition and
  instantiation of some PDE operators as first step:

  @code
  // compute steady-state fluid pressure from [K]{p} = {q}
  PDE_Integrator<2U>  fluid_pressure;

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

  @section outlook Future Implementations

  @todo (1) Implement AssignNonEssentialConditions( const SIMPLICIAL_COMPLEX<SIMPLEX<dim> >& );

  @todo !!! SKM: Implement the automatic integration over boundaries in the case where surface integrals are present
  (design approved: Garmisch and Colleoli)

  */

} // csmp

#endif






