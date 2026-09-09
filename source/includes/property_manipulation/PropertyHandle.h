// VERSION THAT HANDLES aother variable placements

#ifndef CSMP_PROPERTY_HANDLE4_H
#define CSMP_PROPERTY_HANDLE4_H

#include "ScalarVariable.h"
#include "VectorVariable.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "Exception.h"
#include "CSMP_mathUtilities.h"
#include "Region.h"
#include "Boundary.h"
#include "SplitBoundary.h"

namespace csmp {

/**
@brief PropertyHandle - accessor and mutator of distributed data on
       Region, Boundary, and SplitBoundary subdomains.

@author S.K. Matthaei
@author Stephen G. Roberts
@author (refactored, extended) 2024

@section motivation Motivation

Extends PropertyHandle2 to support all CSMP variable placements,
including those on Boundary (Face-based) and SplitBoundary
(InterFace-based) subdomains, as well as finite-volume sector and
facet integration point placements.

@section subdomains Supported subdomains

The constructor accepts any subdomain name and detects the type
using the primary detection mechanism:

  - model_.ContainsSplitBoundary(name) → SplitBoundary<dim>
  - model_.ContainsBoundary(name)      → Boundary<dim>
  - otherwise                          → Region<dim>

SplitBoundary is checked first because a SplitBoundary name may
also satisfy the Boundary name pattern. Each subdomain is unique
within a Model.

@section placements Supported placements

  Region:        NODE, ELEMENT, ELEMENT_INTEGRATION_POINT,
                 SECTOR_INTEGRATION_POINT, FACET_INTEGRATION_POINT,
                 REGION

  Boundary:      NODE, FACE, FACE_INTEGRATION_POINT,
                 FACE_SECTOR_INTEGRATION_POINT,
                 FACE_FACET_INTEGRATION_POINT

  SplitBoundary: NODE (see @ref interFaceNodes),
                 INTER_FACE, INTER_FACE_INTEGRATION_POINT,
                 INTER_FACE_SECTOR_INTEGRATION_POINT,
                 INTER_FACE_FACET_INTEGRATION_POINT

The legacy BOUNDARY placement is rejected at construction time with
a descriptive error. Use FACE or INTER_FACE instead.

@section fv Finite-volume integration points

Sector and facet integration points always have exactly one IP per
sector or facet (ip = 0U always). The sector count per cell equals
the number of nodes in the underlying finite element; the facet
count equals Segments().

@section interFaceNodes InterFace node access

A SplitBoundary does not maintain a NodeVector() because its
InterFace objects have two distinct sets of collocated but
independent nodes: one set belonging to the higher-dimensional
element on the INSIDE (accessible via InterFace::Parent(INSIDE))
and one set belonging to the higher-dimensional element on the
OUTSIDE (accessible via InterFace::Parent(OUTSIDE)). When an
InterveningElement is present, a third set of MIDDLE nodes connects
to the lower-dimensional mesh sandwiched between the two sides.

For NODE-placed operations on a SplitBoundary, PropertyHandle
iterates the nodes of the underlying finite element of each
InterFace via:

@code
for ( uint32_t i = 0; i < interface.FE()->Nodes(); ++i )
{
    node_inside  = interface.N( i, INSIDE );
    node_outside = interface.N( i, OUTSIDE );
    if ( interface.HasInterveningElement() )
        node_middle = interface.N( i, MIDDLE );
}
@endcode

Since nodes are shared between adjacent InterFace objects, an
unordered_set is used to deduplicate and ensure each node is
visited exactly once. INSIDE, OUTSIDE, and MIDDLE nodes are all
visited and operated on independently.

@section interpolation Automatic placement alignment

When two handles with different placements are combined in an
arithmetic expression, the right-hand side is automatically
interpolated onto the placement of the left-hand side using the
ModelSubDomain interpolation methods. The following placement
pairs are supported:

  Region / Boundary:
    NODE                    <-> ELEMENT / FACE
    NODE                    <-> ELEMENT_INTEGRATION_POINT /
                                FACE_INTEGRATION_POINT
    ELEMENT / FACE          <-> ELEMENT_INTEGRATION_POINT /
                                FACE_INTEGRATION_POINT
    ELEMENT / FACE           -> FACET_INTEGRATION_POINT /
                                FACE_FACET_INTEGRATION_POINT

  SplitBoundary:
    INTER_FACE              <-> INTER_FACE_INTEGRATION_POINT
    INTER_FACE               -> INTER_FACE_FACET_INTEGRATION_POINT

All interpolation is scoped to the handle's subdomain, not the
whole Model.

@section splitBoundaryInterpolation SplitBoundary interpolation limitation

For SplitBoundary handles, interpolation paths involving NODE
placement are not supported. Although each node belongs
unambiguously to one side of the InterFace — INSIDE nodes belong
to Parent(INSIDE) and OUTSIDE nodes belong to Parent(OUTSIDE) —
PropertyHandle cannot determine which side the user intends to
interpolate from without explicit guidance. Furthermore, the
interpolation from a chosen set of nodes to the InterFace
barycentre or integration points requires evaluating finite element
shape functions at the parametric coordinates of the InterFace
within the chosen parent element face, which is beyond the scope
of PropertyHandle.

To perform this interpolation:
  (1) Access the nodes explicitly via InterFace::N(i, INSIDE) or
      InterFace::N(i, OUTSIDE), choose the appropriate side, and
      use InterFace::Parent(INSIDE) or InterFace::Parent(OUTSIDE)
      with the finite element interpolation machinery directly.
  (2) Perform the interpolation on the higher-dimensional Region
      first using a whole-Model PropertyHandle, then copy the
      result to the SplitBoundary handle.
      
@section crosssubdomain Cross-subdomain assignment

The only supported cross-subdomain operation is assigning from a
whole-Model NODE-placed handle into a Boundary or SplitBoundary
NODE-placed handle. This is valid because boundary and split-boundary
nodes are a subset of the model nodes.

All other cross-subdomain combinations are rejected:
  - FACE placement exists exclusively on Boundary subdomains.
  - INTER_FACE placement exists exclusively on SplitBoundary subdomains.
  - ELEMENT and integration point placements exist exclusively on
    Region subdomains.
  - Assigning from a sub-domain into the whole model is always
    rejected because it would only partially overwrite the model.

@section angles Angles

All trigonometric methods (Sin, Cos, Tan, Acos, Asin, Atan) expect
and produce values in degrees.

@section arrays Array variables

ArrayVariable and FlaggedArrayVariable are not supported.
Attempting to construct a PropertyHandle for such a variable
throws an ERROR exception at construction time with a descriptive
message. For array-typed variables, iterate directly over the
subdomain's cells or nodes using Read / Store.

@section ownership Ownership and lifetime

A handle that creates a new variable owns it and deletes it on
destruction, provided the variable still exists in the model
database (guarded by IsDefined to handle model teardown order).
A copy-constructed handle never owns the variable; the original
retains ownership. Move construction and move assignment are
deleted because the class holds references that cannot be rebound.

@section thread_safety Thread safety

The iteration helpers are single-threaded. Parallelisation, if
required, should be added at the call site.
*/
template<uint32_t dim>
class PropertyHandle
{
public:

    // -----------------------------------------------------------------------
    // Construction / destruction
    // -----------------------------------------------------------------------

    /**
    Binds to (or creates) @p var_name scoped to the whole Model.
    */
    PropertyHandle( Model<dim>&   model,
                     const char*   var_name,
                     VARIABLE_TYPE type  = SCALAR,
                     PLACEMENT     place = NODE,
                     uint32_t      vsize = 1U );

    /**
    Binds to (or creates) @p var_name restricted to the subdomain
    identified by @p subdomain_name.

    The subdomain type is detected automatically:
      - ContainsBoundary(subdomain_name)      → Boundary
      - ContainsSplitBoundary(subdomain_name) → SplitBoundary
      - otherwise                             → Region
    */
    PropertyHandle( Model<dim>&   model,
                     const char*   subdomain_name,
                     const char*   var_name,
                     VARIABLE_TYPE type  = SCALAR,
                     PLACEMENT     place = NODE,
                     uint32_t      vsize = 1U );

    /**
    Copy construction: shares the same underlying variable.
    The copy does NOT take ownership.
    */
    PropertyHandle( const PropertyHandle& other );

    /** Move construction is deleted: references cannot be rebound. */
    PropertyHandle( PropertyHandle&& ) = delete;

    /**
    Destroys any variable created by this handle, provided it still
    exists in the model database.
    */
    ~PropertyHandle();

    // -----------------------------------------------------------------------
    // Assignment
    // -----------------------------------------------------------------------

    /**
    Copies distributed values from @p other's variable into this
    handle's variable, interpolating between placements as needed.
    Copying from a whole-Model handle into a sub-domain handle is
    permitted. The reverse is not.
    */
    PropertyHandle& operator=( const PropertyHandle& other );

    /** Move assignment is deleted. */
    PropertyHandle& operator=( PropertyHandle&& ) = delete;

    /** Assigns a uniform scalar value to every qualifying point. */
    PropertyHandle& operator=( double val );

    /** Assigns the value of @p s to every qualifying point. */
    PropertyHandle& operator=( const ScalarVariable& s );

    /**
    Assigns the vector @p v to every qualifying point.
    If this handle holds a SCALAR variable, the vector length is
    assigned with a diagnostic message.
    */
    PropertyHandle& operator=( const VectorVariable<dim>& v );

    /**
    Assigns the tensor @p t to every qualifying point.
    If this handle holds a SCALAR variable, the tensor determinant
    is assigned with a diagnostic message.
    */
    PropertyHandle& operator=( const TensorVariable<dim>& t );

    /**
    Assigns a pre-built vector of VectorVariables directly into
    storage via FEM_Data.
    */
    PropertyHandle& operator=( const std::vector<VectorVariable<dim>>& vc );

    /**
    Assigns a pre-built vector of TensorVariables directly into
    storage via FEM_Data.
    */
    PropertyHandle& operator=( const std::vector<TensorVariable<dim>>& ts );
    
        /**
    Assigns the values of @p av to every qualifying point.
    The array is broadcast uniformly — every point receives the same
    values. The flag_output_ guard is applied via ArrayVariable::Flag().
    Throws ERROR if this handle does not hold an ARRAY variable.
    */
    PropertyHandle& operator=( const ArrayVariable& av );

    /**
    Assigns the values of @p fav to every qualifying point.
    The array is broadcast uniformly — every point receives the same
    values and flags. Throws ERROR if this handle does not hold a
    FLAGGEDARRAY variable.
    */
    PropertyHandle& operator=( const FlaggedArrayVariable& fav );


    // -----------------------------------------------------------------------
    // Compound assignment — scalar right-hand side
    // -----------------------------------------------------------------------

    /** Adds @p val to every component of every qualifying point. */
    PropertyHandle& operator+=( double val );

    /** Subtracts @p val from every component of every qualifying point. */
    PropertyHandle& operator-=( double val );

    /** Multiplies every component of every qualifying point by @p val. */
    PropertyHandle& operator*=( double val );

    /**
    Divides every component of every qualifying point by @p val.
    Throws ERROR if @p val is zero.
    */
    PropertyHandle& operator/=( double val );

    // -----------------------------------------------------------------------
    // Compound assignment — PropertyHandle right-hand side
    // -----------------------------------------------------------------------

    /**
    Adds the values of @p other's variable to this handle's variable.
    Placement alignment is performed automatically; any temporary
    property created is deleted immediately after the operation.
    */
    PropertyHandle& operator+=( const PropertyHandle& other );

    /**
    Subtracts the values of @p other's variable from this handle's
    variable. Placement alignment is performed automatically.
    */
    PropertyHandle& operator-=( const PropertyHandle& other );

    /**
    Multiplies this handle's variable by @p other's variable.
    Placement alignment is performed automatically.
    Mixed-type operations (VECTOR *= SCALAR, TENSOR *= SCALAR) are
    supported component-wise.
    */
    PropertyHandle& operator*=( const PropertyHandle& other );

    /**
    Divides this handle's variable by @p other's variable.
    Placement alignment is performed automatically.
    */
    PropertyHandle& operator/=( const PropertyHandle& other );

    // -----------------------------------------------------------------------
    // Component-wise mathematical operations
    // -----------------------------------------------------------------------

    /**
    Applies any callable @p f with signature void(double&) to every
    component of every qualifying point.

    Example:
    @code
    phi.Apply( []( double& x ) { x = 1.0 - std::exp( -x / 0.3 ); } );
    @endcode
    */
    template<typename F>
    void Apply( F&& f );

    /**
    Applies callable @p f with signature void(ArrayVariable&) to the
    ArrayVariable at every qualifying point.

    The whole array is read, passed to @p f, and written back in a
    single read-modify-write cycle per point. The flag_output_ guard
    is applied once per point via ArrayVariable::Flag() — if the
    array's flag does not match flag_output_, the point is skipped
    entirely.

    Only valid for ARRAY-typed handles. Throws ERROR for other types.

    Example:
    @code
    ph.ApplyArray( []( ArrayVariable& av )
    {
        for ( uint32_t i = 0; i < av.Size(); ++i )
            av(i) = std::log( av(i) );
    });
    @endcode
    */
    template<typename F>
    void ApplyArray( F&& f );

    /**
    Applies callable @p f with signature void(FlaggedArrayVariable&)
    to the FlaggedArrayVariable at every point.

    The whole array is read, passed to @p f, and written back in a
    single read-modify-write cycle per point. No whole-variable flag
    check is applied — @p f receives the complete FlaggedArrayVariable
    including its per-element flags and is responsible for per-element
    flag checking. Use ApplyToFlaggedElements() as a helper.

    Only valid for FLAGGEDARRAY-typed handles. Throws ERROR for other
    types.

    Example:
    @code
    ph.ApplyFlaggedArray( [&]( FlaggedArrayVariable& fav )
    {
        PropertyHandle3::ApplyToFlaggedElements( fav, ph.OutputCondition(),
            []( double& x ) { x = std::exp(x); } );
    });
    @endcode
    */
    template<typename F>
    void ApplyFlaggedArray( F&& f );

    /**
    Helper for use inside ApplyFlaggedArray lambdas.
    Applies @p op(double&) to each element of @p fav where the
    element flag equals @p flag_output.

    This provides the per-element flag checking that ApplyFlaggedArray
    delegates to the caller, consistent with the design intent of
    FlaggedArrayVariable (some elements fixed, others free).

    @param fav         The FlaggedArrayVariable to operate on.
    @param flag_output Only elements whose flag matches this value
                       are passed to @p op.
    @param op          Callable with signature void(double&).
    */
    template<typename Op>
    static void ApplyToFlaggedElements( FlaggedArrayVariable& fav,
                                         VARIABLE_FLAG         flag_output,
                                         Op&&                  op );

    /**
    Component-wise squaring: each component T_ij becomes T_ij^2.
    Equivalent to the Hadamard product T ∘ T.
    Valid for SCALAR, VECTOR, TENSOR, ARRAY, FLAGGEDARRAY.
    */
    void Squared();

    /**
    Matrix product of a tensor with itself: result_ij = sum_k T_ik * T_kj.
    Only valid for TENSOR handles. Throws ERROR for other types.
    */
    void MatrixSquared();

    /**
    Double contraction T:T = sum_ij T_ij^2, producing a scalar result.
    Only valid for TENSOR handles. Throws ERROR for other types.
    Note: this reduces the tensor to a scalar — the result is stored
    as a SCALAR variable. The handle must therefore be a SCALAR handle
    receiving the result, not the TENSOR handle itself.
    This operation cannot be performed in-place on a TENSOR handle.

    Double contraction T:T = sum_ij T_ij^2 stored into a SCALAR handle.
    Only valid when this handle holds SCALAR and @p tensor holds TENSOR
    on the same placement. Throws ERROR otherwise.
    */
    void DoubleContraction( const PropertyHandle& tensor );

    /** Component-wise square root. */
    void Sqrt();

    /**
    Natural logarithm (base e), component-wise.
    Non-positive components are left unchanged.
    */
    void Ln();

    /**
    Decadic logarithm (base 10), component-wise.
    Non-positive components are left unchanged.
    */
    void Log10();

    /** Exponential function e^x, component-wise. */
    void Exp();

    /** Raises each component to the power @p exponent. */
    void Pow( double exponent );

    /** Replaces every NaN component with @p replacement. */
    void ZapNAN( double replacement );

    /** Sine of each component. Input values must be in degrees. */
    void Sin();

    /** Cosine of each component. Input values must be in degrees. */
    void Cos();

    /** Tangent of each component. Input values must be in degrees. */
    void Tan();

    /**
    Arc cosine of each component. Output is in degrees.
    */
    void Acos();

    /**
    Arc sine of each component. Output is in degrees.
    */
    void Asin();

    /**
    Arc tangent of each component. Output is in degrees.
    */
    void Atan();

    /** Absolute value of each component. */
    void Abs();

    /**
    Clamps every component to [@p lo, @p hi].
    Throws ERROR if @p lo > @p hi.
    */
    void Clamp( double lo, double hi );

    // -----------------------------------------------------------------------
    // Queries
    // -----------------------------------------------------------------------

    /** Returns the name of the associated variable. */
    const char* VariableName() const;

    /** Returns the storage key of the associated variable. */
    const csmp::Index& Key() const;

    /**
    Returns the global minimum and maximum of the associated variable
    into @p omin and @p omax.
    */
    void Range( double& omin, double& omax ) const;

    /**
    Returns true if all values lie within the range registered in
    the PropertyDatabase. For TENSOR variables, eigenvalues are used
    where applicable (see IsWithinRange implementation).
    */
    bool IsWithinRange() const;

    /**
    Returns the VARIABLE_FLAG that a point's flag must match for the
    handle to modify it. Default is ANY.
    */
    VARIABLE_FLAG OutputCondition() const;

    /** Sets the VARIABLE_FLAG guard to @p c. */
    void OutputCondition( VARIABLE_FLAG c );

    // -----------------------------------------------------------------------
    // Output
    // -----------------------------------------------------------------------

    /** Prints variable metadata and distributed values to stdout. */
    void Out() const;

    /**
    Writes variable metadata and distributed values to @p filename.
    Throws ERROR if the file cannot be opened.
    */
    void Out( const char* filename ) const;

private:

    // -----------------------------------------------------------------------
    // Private helpers
    // -----------------------------------------------------------------------

    /**
    Shared constructor body. Detects subdomain type, validates
    placement and variable type, binds to an existing variable or
    creates a new one, and sets key_.
    */
    void Initialise( const char*   subdomain_name,
                     const char*   var_name,
                     VARIABLE_TYPE type,
                     PLACEMENT     place,
                     uint32_t      vsize );

    /**
    Dispatches to the appropriate ModelSubDomain interpolation method.
    Throws ERROR if no path exists between @p src and @p dst.

    @note For SplitBoundary handles, interpolation paths involving NODE
    placement are not supported. Interpolating between NODE and
    INTER_FACE requires finite element shape function evaluation at
    the parametric coordinates of the InterFace within its parent
    element face, which is beyond the scope of PropertyHandle.
    Use InterFace::Parent(INSIDE) or InterFace::Parent(OUTSIDE) with
    the finite element interpolation machinery directly, or perform
    the interpolation on the higher-dimensional Region first.
    */
    void Interpolate( const char* src_name,
                      const char* dst_name,
                      PLACEMENT   src,
                      PLACEMENT   dst );

    /**
    Assigns the length of the VECTOR variable identified by @p rkey
    to this handle's SCALAR variable at every qualifying point.
    */
    void AssignVectorLengthToScalar( const csmp::Index& rkey );

    /**
    Assigns the determinant of the TENSOR variable identified by
    @p rkey to this handle's SCALAR variable at every qualifying point.
    */
    void AssignTensorDetToScalar( const csmp::Index& rkey );

    /**
    Copies values from the FLAGGEDARRAY variable identified by @p rkey
    into this handle's ARRAY variable at every qualifying point,
    discarding per-element flags. The flag_output_ guard is applied
    once per point via the whole-array flag.
    Used by operator=(PropertyHandle3) for the ARRAY = FLAGGEDARRAY
    cross-type path.
    */
    void AssignFlaggedArrayToArray( const csmp::Index& rkey );

    /**
    Returns an iterable range over the cells of the active subdomain.
    Dispatches to Region, Boundary, or SplitBoundary CellVector().
    */
    auto CellRange() const;

    /**
    Returns an iterable range over the nodes of the active subdomain.
    Only valid for Region and Boundary; throws for SplitBoundary
    (use InterFace FE node iteration instead).
    */
    auto NodeRange() const;
    
    /**
    Calls f( active_subdomain ) where active_subdomain is whichever of
    region_, boundary_, or split_boundary_ is non-null.
    F must be a generic lambda accepting any of the three types.
    */
    template<typename F>
    void WithSubdomain( F&& f ) const;

    template<typename F>
    void WithSubdomain( F&& f );

    /**
    Applies a binary operation between *this and @p rhs, which must
    already be on the same placement as *this.

      SF: void( double&, double )
      VF: void( VectorVariable<dim>&, const VectorVariable<dim>& )
      TF: void( TensorVariable<dim>&, const TensorVariable<dim>& )

    Mixed-type combinations (VECTOR/SCALAR, TENSOR/SCALAR) are
    handled by routing through SF component-wise.
    */
    template<typename SF, typename VF, typename TF>
    void ApplyBinaryOpSamePlacement( const PropertyHandle& rhs,
                                     SF&& sop,
                                     VF&& vop,
                                     TF&& top );

    /**
    Applies a binary operation between *this and @p other, creating
    and immediately deleting a temporary property if the placements
    differ. Delegates to ApplyBinaryOpSamePlacement once aligned.
    */
    template<typename SF, typename VF, typename TF>
    PropertyHandle& ApplyBinaryOpAligned( const PropertyHandle& other,
                                           SF&& sop,
                                           VF&& vop,
                                           TF&& top );

    /**
    Iterates over every qualifying point and calls f(double&) in a
    read-modify-write cycle on the scalar value stored there.
    The flag_output_ guard is applied at every point.
    Handles all placements including FV sector/facet IPs and
    InterFace node access.
    */
    template<typename F>
    void ApplyScalar( F&& f );

    /**
    Iterates over every qualifying point and calls
    f(VectorVariable<dim>&) in a read-modify-write cycle.
    */
    template<typename F>
    void ApplyVector( F&& f );

    /**
    Iterates over every qualifying point and calls
    f(TensorVariable<dim>&) in a read-modify-write cycle.
    */
    template<typename F>
    void ApplyTensor( F&& f );

    /**
    Dispatches to ApplyScalar, ApplyVector, or ApplyTensor based on
    key_.type. Throws ERROR for ARRAY and FLAGGED_ARRAY types.
    */
    template<typename SF, typename VF, typename TF>
    void ApplyByType( SF&& sf, VF&& vf, TF&& tf );
    
    /**
    Private helper used by all named math methods (Pow, Sqrt, Ln, etc.)
    to dispatch the correct operation across all five variable types.

      ScalarOp:       void(double&)              — used for SCALAR/VECTOR/TENSOR
      ArrayOp:        void(ArrayVariable&)       — used for ARRAY
      FlaggedArrayOp: void(FlaggedArrayVariable&) — used for FLAGGEDARRAY

    For SCALAR/VECTOR/TENSOR, delegates to Apply(scalar_op) which
    iterates components via applyScalar/applyVector/applyTensor.
    For ARRAY, delegates to applyArray(array_op).
    For FLAGGEDARRAY, delegates to applyFlaggedArray(farray_op).
    Calls IsWithinRange() after the operation.
    */
    template<typename ScalarOp, typename ArrayOp, typename FlaggedArrayOp>
    void ApplyMathOp( ScalarOp&&       scalar_op,
                      ArrayOp&&        array_op,
                      FlaggedArrayOp&& farray_op );
                  
    // -----------------------------------------------------------------------
    // Data members
    // -----------------------------------------------------------------------

    /** Reference to the Model that owns the property database. */
    Model<dim>&   model_;

    /**
    Exactly one of these three pointers is non-null, depending on
    which subdomain type was detected at construction time.
    All are non-owning; always valid for the lifetime of the model.
    */
    Region<dim>*        region_         = nullptr;
    Boundary<dim>*      boundary_       = nullptr;
    SplitBoundary<dim>* split_boundary_ = nullptr;;

    /** Name of the subdomain; used for consistency checks. */
    std::string          subdomain_name_;

    /**
    True if the subdomain is a SplitBoundary.
    Required to select the correct node-access strategy for
    NODE-placed operations (InterFace FE nodes, deduplicated).
    */
    bool                 is_split_boundary_;

    /** Name of the associated variable in the PropertyDatabase. */
    std::string          var_name_;

    /**
    Storage key for the associated variable. Kept current via the
    IndexTracker / Observer pattern in the PropertyDatabase.
    */
    csmp::Index          key_;

    /**
    A point is only modified if its VARIABLE_FLAG matches this value.
    Default is ANY.
    */
    VARIABLE_FLAG        flag_output_;

    /**
    True if and only if this handle created the variable. Controls
    whether the destructor calls model_.DeleteProperty.
    */
    bool                 owns_variable_;
};


// ============================================================================
//  Inline template member definitions
//  All inner template members must be defined here so that the lambda
//  types deduced at each call site are visible to the compiler.
// ============================================================================

template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::WithSubdomain( F&& f ) const
{
    if      ( region_        ) f( *region_ );
    else if ( boundary_      ) f( *boundary_ );
    else if ( split_boundary_) f( *split_boundary_ );
    else
        throw csmp::Exception( ERROR,
            "PropertyHandle::WithSubdomain",
            var_name_.c_str(),
            "No subdomain is bound to this handle." );
}

template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::WithSubdomain( F&& f )
{
    if      ( region_        ) f( *region_ );
    else if ( boundary_      ) f( *boundary_ );
    else if ( split_boundary_) f( *split_boundary_ );
    else
        throw csmp::Exception( ERROR,
            "PropertyHandle::WithSubdomain",
            var_name_.c_str(),
            "No subdomain is bound to this handle." );
}

// ============================================================================
//  Apply
// ============================================================================
template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::Apply( F&& f )
{
    ApplyByType(
        [&]( double& sc )              { f( sc ); },
        [&]( VectorVariable<dim>& vc ) { for ( uint32_t j=0; j<dim; ++j ) f( vc(j) ); },
        [&]( TensorVariable<dim>& ts ) { for ( uint32_t i=0; i<dim; ++i )
                                             for ( uint32_t j=0; j<dim; ++j ) f( ts(i,j) ); }
    );
}





// ============================================================================
//  ApplyScalar
// ============================================================================

template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::ApplyScalar( F&& f )
{
    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                // Region or Boundary — NodeVector() is available.
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                        if ( nit->Status( key_ ) == flag_output_ )
                        {
                            double sc = nit->Read( key_ );
                            f( sc );
                            nit->Store( key_, makeScalar( flag_output_, sc ) );
                        }
                });
            }
            else
            {
                // SplitBoundary — iterate InterFace FE nodes, deduplicated.
                // Call split_boundary_ directly to avoid instantiating
                // InterFace-specific methods (N, HasInterveningElement)
                // for Element and Face types inside WithSubdomain.
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                    {
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                            {
                                double sc = nd->Read( key_ );
                                f( sc );
                                nd->Store( key_, makeScalar( flag_output_, sc ) );
                            }
                        }
                        if ( ifit->HasInterveningElement() )
                        {
                            auto* nd = ifit->N( i, MIDDLE );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                            {
                                double sc = nd->Read( key_ );
                                f( sc );
                                nd->Store( key_, makeScalar( flag_output_, sc ) );
                            }
                        }
                    }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    if ( cit->Status( key_ ) == flag_output_ )
                    {
                        double sc = cit->Read( key_ );
                        f( sc );
                        cit->Store( key_, makeScalar( flag_output_, sc ) );
                    }
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        if ( cit->Status( ip, key_ ) == flag_output_ )
                        {
                            double sc = cit->Read( ip, key_ );
                            f( sc );
                            cit->Store( ip, key_, makeScalar( flag_output_, sc ) );
                        }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        if ( cit->Status( s, 0U, key_ ) == flag_output_ )
                        {
                            double sc = cit->Read( s, 0U, key_ );
                            f( sc );
                            cit->Store( s, 0U, key_, makeScalar( flag_output_, sc ) );
                        }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        if ( cit->Status( fac, 0U, key_ ) == flag_output_ )
                        {
                            double sc = cit->Read( fac, 0U, key_ );
                            f( sc );
                            cit->Store( fac, 0U, key_, makeScalar( flag_output_, sc ) );
                        }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                if ( sd.Status( key_ ) == flag_output_ )
                {
                    double sc = sd.Read( key_ );
                    f( sc );
                    sd.Store( key_, makeScalar( flag_output_, sc ) );
                }
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle::ApplyScalar",
                var_name_.c_str(),
                "Unsupported placement." );
    }
}



template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::ApplyVector( F&& f )
{
    VectorVariable<dim> vc;

    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                    {
                        nit->Read( key_, vc );
                        f( vc );
                        nit->Store( key_, vc );
                    }
                });
            }
            else
            {
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                    {
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second )
                            {
                                nd->Read( key_, vc );
                                f( vc );
                                nd->Store( key_, vc );
                            }
                        }
                        if ( ifit->HasInterveningElement() )
                        {
                            auto* nd = ifit->N( i, MIDDLE );
                            if ( visited.insert( nd ).second )
                            {
                                nd->Read( key_, vc );
                                f( vc );
                                nd->Store( key_, vc );
                            }
                        }
                    }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                {
                    cit->Read( key_, vc );
                    f( vc );
                    cit->Store( key_, vc );
                }
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                    {
                        cit->Read( ip, key_, vc );
                        f( vc );
                        cit->Store( ip, key_, vc );
                    }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                    {
                        cit->Read( s, 0U, key_, vc );
                        f( vc );
                        cit->Store( s, 0U, key_, vc );
                    }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                    {
                        cit->Read( fac, 0U, key_, vc );
                        f( vc );
                        cit->Store( fac, 0U, key_, vc );
                    }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                sd.Read( key_, vc );
                f( vc );
                sd.Store( key_, vc );
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle::ApplyVector",
                var_name_.c_str(),
                "Unsupported placement." );
    }
}



template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::ApplyTensor( F&& f )
{
    TensorVariable<dim> ts;

    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                    {
                        nit->Read( key_, ts );
                        f( ts );
                        nit->Store( key_, ts );
                    }
                });
            }
            else
            {
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                    {
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second )
                            {
                                nd->Read( key_, ts );
                                f( ts );
                                nd->Store( key_, ts );
                            }
                        }
                        if ( ifit->HasInterveningElement() )
                        {
                            auto* nd = ifit->N( i, MIDDLE );
                            if ( visited.insert( nd ).second )
                            {
                                nd->Read( key_, ts );
                                f( ts );
                                nd->Store( key_, ts );
                            }
                        }
                    }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                {
                    cit->Read( key_, ts );
                    f( ts );
                    cit->Store( key_, ts );
                }
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                    {
                        cit->Read( ip, key_, ts );
                        f( ts );
                        cit->Store( ip, key_, ts );
                    }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                    {
                        cit->Read( s, 0U, key_, ts );
                        f( ts );
                        cit->Store( s, 0U, key_, ts );
                    }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                    {
                        cit->Read( fac, 0U, key_, ts );
                        f( ts );
                        cit->Store( fac, 0U, key_, ts );
                    }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                sd.Read( key_, ts );
                f( ts );
                sd.Store( key_, ts );
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle::ApplyTensor",
                var_name_.c_str(),
                "Unsupported placement." );
    }
}



// ============================================================================
//  ApplyToFlaggedElements  (static helper)
// ============================================================================
template<uint32_t dim>
template<typename Op>
inline void PropertyHandle<dim>::ApplyToFlaggedElements(
    FlaggedArrayVariable& fav,
    VARIABLE_FLAG         flag_output,
    Op&&                  op )
{
    for ( uint32_t i = 0; i < fav.Size(); ++i )
        if ( fav.Flag(i) == flag_output )
            op( fav(i) );
}

// ============================================================================
//  applyArray
// ============================================================================
template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::ApplyArray( F&& f )
{
    ArrayVariable av;

    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                        if ( nit->Status( key_ ) == flag_output_ )
                        {
                            nit->Read( key_, av );
                            f( av );
                            nit->Store( key_, av );
                        }
                });
            }
            else
            {
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                    {
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                            {
                                nd->Read( key_, av );
                                f( av );
                                nd->Store( key_, av );
                            }
                        }
                        if ( ifit->HasInterveningElement() )
                        {
                            auto* nd = ifit->N( i, MIDDLE );
                            if ( visited.insert( nd ).second &&
                                 nd->Status( key_ ) == flag_output_ )
                            {
                                nd->Read( key_, av );
                                f( av );
                                nd->Store( key_, av );
                            }
                        }
                    }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    if ( cit->Status( key_ ) == flag_output_ )
                    {
                        cit->Read( key_, av );
                        f( av );
                        cit->Store( key_, av );
                    }
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        if ( cit->Status( ip, key_ ) == flag_output_ )
                        {
                            cit->Read( ip, key_, av );
                            f( av );
                            cit->Store( ip, key_, av );
                        }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        if ( cit->Status( s, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( s, 0U, key_, av );
                            f( av );
                            cit->Store( s, 0U, key_, av );
                        }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        if ( cit->Status( fac, 0U, key_ ) == flag_output_ )
                        {
                            cit->Read( fac, 0U, key_, av );
                            f( av );
                            cit->Store( fac, 0U, key_, av );
                        }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                if ( sd.Status( key_ ) == flag_output_ )
                {
                    sd.Read( key_, av );
                    f( av );
                    sd.Store( key_, av );
                }
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle3::applyArray",
                var_name_.c_str(),
                "Unsupported placement." );
    }
}

// ============================================================================
//  applyFlaggedArray
// ============================================================================
template<uint32_t dim>
template<typename F>
inline void PropertyHandle<dim>::ApplyFlaggedArray( F&& f )
{
    FlaggedArrayVariable fav;

    switch ( key_.place )
    {
        case NODE:
            if ( !is_split_boundary_ )
            {
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& nit : sd.NodeVector() )
                    {
                        nit->Read( key_, fav );
                        f( fav );
                        nit->Store( key_, fav );
                    }
                });
            }
            else
            {
                std::unordered_set<Node<dim>*> visited;
                for ( auto& ifit : split_boundary_->CellVector() )
                    for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                    {
                        for ( auto side : { INSIDE, OUTSIDE } )
                        {
                            auto* nd = ifit->N( i, side );
                            if ( visited.insert( nd ).second )
                            {
                                nd->Read( key_, fav );
                                f( fav );
                                nd->Store( key_, fav );
                            }
                        }
                        if ( ifit->HasInterveningElement() )
                        {
                            auto* nd = ifit->N( i, MIDDLE );
                            if ( visited.insert( nd ).second )
                            {
                                nd->Read( key_, fav );
                                f( fav );
                                nd->Store( key_, fav );
                            }
                        }
                    }
            }
            break;

        case ELEMENT:
        case FACE:
        case INTER_FACE:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                {
                    cit->Read( key_, fav );
                    f( fav );
                    cit->Store( key_, fav );
                }
            });
            break;

        case ELEMENT_INTEGRATION_POINT:
        case FACE_INTEGRATION_POINT:
        case INTER_FACE_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                    {
                        cit->Read( ip, key_, fav );
                        f( fav );
                        cit->Store( ip, key_, fav );
                    }
            });
            break;

        case SECTOR_INTEGRATION_POINT:
        case FACE_SECTOR_INTEGRATION_POINT:
        case INTER_FACE_SECTOR_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                    {
                        cit->Read( s, 0U, key_, fav );
                        f( fav );
                        cit->Store( s, 0U, key_, fav );
                    }
            });
            break;

        case FACET_INTEGRATION_POINT:
        case FACE_FACET_INTEGRATION_POINT:
        case INTER_FACE_FACET_INTEGRATION_POINT:
            WithSubdomain( [&]( auto& sd )
            {
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                    {
                        cit->Read( fac, 0U, key_, fav );
                        f( fav );
                        cit->Store( fac, 0U, key_, fav );
                    }
            });
            break;

        case REGION:
            WithSubdomain( [&]( auto& sd )
            {
                sd.Read( key_, fav );
                f( fav );
                sd.Store( key_, fav );
            });
            break;

        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle3::applyFlaggedArray",
                var_name_.c_str(),
                "Unsupported placement." );
    }
}








// ============================================================================
//  ApplyByType
// ============================================================================

template<uint32_t dim>
template<typename SF, typename VF, typename TF>
inline void PropertyHandle<dim>::ApplyByType( SF&& sf, VF&& vf, TF&& tf )
{
    switch ( key_.type )
    {
        case SCALAR:      ApplyScalar( std::forward<SF>( sf ) ); break;
        case VECTOR:      ApplyVector( std::forward<VF>( vf ) ); break;
        case TENSOR:      ApplyTensor( std::forward<TF>( tf ) ); break;
        case ARRAY:
            throw csmp::Exception( ERROR,
                "PropertyHandle3::applyByType",
                var_name_.c_str(),
                "Use ApplyArray() for ARRAY-typed handles." );
        case FLAGGEDARRAY:
            throw csmp::Exception( ERROR,
                "PropertyHandle3::applyByType",
                var_name_.c_str(),
                "Use ApplyFlaggedArray() for FLAGGEDARRAY-typed handles." );
        default:
            throw csmp::Exception( ERROR,
                "PropertyHandle3::applyByType",
                var_name_.c_str(),
                "Unrecognised variable type." );
    }
}



// ============================================================================
//  ApplyBinaryOpSamePlacement
// ============================================================================
template<uint32_t dim>
template<typename SF, typename VF, typename TF>
inline void PropertyHandle<dim>::ApplyBinaryOpSamePlacement(
    const PropertyHandle<dim>& rhs,
    SF&& sop,
    VF&& vop,
    TF&& top )
{
    const csmp::Index& lkey = key_;
    const csmp::Index& rkey = rhs.key_;

    if ( lkey.type == rkey.type )
    {
        switch ( lkey.type )
        {
            case SCALAR:
                switch ( lkey.place )
                {
                    case NODE:
                        if ( !is_split_boundary_ )
                        {
                            WithSubdomain( [&]( auto& sd )
                            {
                                for ( auto& nit : sd.NodeVector() )
                                    if ( nit->Status( lkey ) == flag_output_ )
                                    {
                                        double a = nit->Read( lkey );
                                        double b = nit->Read( rkey );
                                        sop( a, b );
                                        nit->Store( lkey, makeScalar( flag_output_, a ) );
                                    }
                            });
                        }
                        else
                        {
                            std::unordered_set<Node<dim>*> visited;
                            for ( auto& ifit : split_boundary_->CellVector() )
                                for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                                {
                                    for ( auto side : { INSIDE, OUTSIDE } )
                                    {
                                        auto* nd = ifit->N( i, side );
                                        if ( visited.insert( nd ).second &&
                                             nd->Status( lkey ) == flag_output_ )
                                        {
                                            double a = nd->Read( lkey );
                                            double b = nd->Read( rkey );
                                            sop( a, b );
                                            nd->Store( lkey, makeScalar( flag_output_, a ) );
                                        }
                                    }
                                    if ( ifit->HasInterveningElement() )
                                    {
                                        auto* nd = ifit->N( i, MIDDLE );
                                        if ( visited.insert( nd ).second &&
                                             nd->Status( lkey ) == flag_output_ )
                                        {
                                            double a = nd->Read( lkey );
                                            double b = nd->Read( rkey );
                                            sop( a, b );
                                            nd->Store( lkey, makeScalar( flag_output_, a ) );
                                        }
                                    }
                                }
                        }
                        break;

                    case ELEMENT:
                    case FACE:
                    case INTER_FACE:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                if ( cit->Status( lkey ) == flag_output_ )
                                {
                                    double a = cit->Read( lkey );
                                    double b = cit->Read( rkey );
                                    sop( a, b );
                                    cit->Store( lkey, makeScalar( flag_output_, a ) );
                                }
                        });
                        break;

                    case ELEMENT_INTEGRATION_POINT:
                    case FACE_INTEGRATION_POINT:
                    case INTER_FACE_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                                    if ( cit->Status( ip, lkey ) == flag_output_ )
                                    {
                                        double a = cit->Read( ip, lkey );
                                        double b = cit->Read( ip, rkey );
                                        sop( a, b );
                                        cit->Store( ip, lkey, makeScalar( flag_output_, a ) );
                                    }
                        });
                        break;

                    case SECTOR_INTEGRATION_POINT:
                    case FACE_SECTOR_INTEGRATION_POINT:
                    case INTER_FACE_SECTOR_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                                    if ( cit->Status( s, 0U, lkey ) == flag_output_ )
                                    {
                                        double a = cit->Read( s, 0U, lkey );
                                        double b = cit->Read( s, 0U, rkey );
                                        sop( a, b );
                                        cit->Store( s, 0U, lkey, makeScalar( flag_output_, a ) );
                                    }
                        });
                        break;

                    case FACET_INTEGRATION_POINT:
                    case FACE_FACET_INTEGRATION_POINT:
                    case INTER_FACE_FACET_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                                    if ( cit->Status( fac, 0U, lkey ) == flag_output_ )
                                    {
                                        double a = cit->Read( fac, 0U, lkey );
                                        double b = cit->Read( fac, 0U, rkey );
                                        sop( a, b );
                                        cit->Store( fac, 0U, lkey, makeScalar( flag_output_, a ) );
                                    }
                        });
                        break;

                    case REGION:
                        WithSubdomain( [&]( auto& sd )
                        {
                            if ( sd.Status( lkey ) == flag_output_ )
                            {
                                double a = sd.Read( lkey );
                                double b = sd.Read( rkey );
                                sop( a, b );
                                sd.Store( lkey, makeScalar( flag_output_, a ) );
                            }
                        });
                        break;

                    default:
                        throw csmp::Exception( ERROR,
                            "PropertyHandle::ApplyBinaryOpSamePlacement",
                            var_name_.c_str(), "Unsupported placement." );
                }
                break;

            case VECTOR:
            {
                VectorVariable<dim> va, vb;
                switch ( lkey.place )
                {
                    case NODE:
                        if ( !is_split_boundary_ )
                        {
                            WithSubdomain( [&]( auto& sd )
                            {
                                for ( auto& nit : sd.NodeVector() )
                                {
                                    nit->Read( lkey, va );
                                    nit->Read( rkey, vb );
                                    vop( va, vb );
                                    nit->Store( lkey, va );
                                }
                            });
                        }
                        else
                        {
                            std::unordered_set<Node<dim>*> visited;
                            for ( auto& ifit : split_boundary_->CellVector() )
                                for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                                {
                                    for ( auto side : { INSIDE, OUTSIDE } )
                                    {
                                        auto* nd = ifit->N( i, side );
                                        if ( visited.insert( nd ).second )
                                        {
                                            nd->Read( lkey, va );
                                            nd->Read( rkey, vb );
                                            vop( va, vb );
                                            nd->Store( lkey, va );
                                        }
                                    }
                                    if ( ifit->HasInterveningElement() )
                                    {
                                        auto* nd = ifit->N( i, MIDDLE );
                                        if ( visited.insert( nd ).second )
                                        {
                                            nd->Read( lkey, va );
                                            nd->Read( rkey, vb );
                                            vop( va, vb );
                                            nd->Store( lkey, va );
                                        }
                                    }
                                }
                        }
                        break;

                    case ELEMENT:
                    case FACE:
                    case INTER_FACE:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                            {
                                cit->Read( lkey, va );
                                cit->Read( rkey, vb );
                                vop( va, vb );
                                cit->Store( lkey, va );
                            }
                        });
                        break;

                    case ELEMENT_INTEGRATION_POINT:
                    case FACE_INTEGRATION_POINT:
                    case INTER_FACE_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                                {
                                    cit->Read( ip, lkey, va );
                                    cit->Read( ip, rkey, vb );
                                    vop( va, vb );
                                    cit->Store( ip, lkey, va );
                                }
                        });
                        break;

                    case SECTOR_INTEGRATION_POINT:
                    case FACE_SECTOR_INTEGRATION_POINT:
                    case INTER_FACE_SECTOR_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                                {
                                    cit->Read( s, 0U, lkey, va );
                                    cit->Read( s, 0U, rkey, vb );
                                    vop( va, vb );
                                    cit->Store( s, 0U, lkey, va );
                                }
                        });
                        break;

                    case FACET_INTEGRATION_POINT:
                    case FACE_FACET_INTEGRATION_POINT:
                    case INTER_FACE_FACET_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                                {
                                    cit->Read( fac, 0U, lkey, va );
                                    cit->Read( fac, 0U, rkey, vb );
                                    vop( va, vb );
                                    cit->Store( fac, 0U, lkey, va );
                                }
                        });
                        break;

                    case REGION:
                        WithSubdomain( [&]( auto& sd )
                        {
                            sd.Read( lkey, va );
                            sd.Read( rkey, vb );
                            vop( va, vb );
                            sd.Store( lkey, va );
                        });
                        break;

                    default:
                        throw csmp::Exception( ERROR,
                            "PropertyHandle::ApplyBinaryOpSamePlacement",
                            var_name_.c_str(), "Unsupported placement." );
                }
                break;
            }

            case TENSOR:
            {
                TensorVariable<dim> ta, tb;
                switch ( lkey.place )
                {
                    case NODE:
                        if ( !is_split_boundary_ )
                        {
                            WithSubdomain( [&]( auto& sd )
                            {
                                for ( auto& nit : sd.NodeVector() )
                                {
                                    nit->Read( lkey, ta );
                                    nit->Read( rkey, tb );
                                    top( ta, tb );
                                    nit->Store( lkey, ta );
                                }
                            });
                        }
                        else
                        {
                            std::unordered_set<Node<dim>*> visited;
                            for ( auto& ifit : split_boundary_->CellVector() )
                                for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                                {
                                    for ( auto side : { INSIDE, OUTSIDE } )
                                    {
                                        auto* nd = ifit->N( i, side );
                                        if ( visited.insert( nd ).second )
                                        {
                                            nd->Read( lkey, ta );
                                            nd->Read( rkey, tb );
                                            top( ta, tb );
                                            nd->Store( lkey, ta );
                                        }
                                    }
                                    if ( ifit->HasInterveningElement() )
                                    {
                                        auto* nd = ifit->N( i, MIDDLE );
                                        if ( visited.insert( nd ).second )
                                        {
                                            nd->Read( lkey, ta );
                                            nd->Read( rkey, tb );
                                            top( ta, tb );
                                            nd->Store( lkey, ta );
                                        }
                                    }
                                }
                        }
                        break;

                    case ELEMENT:
                    case FACE:
                    case INTER_FACE:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                            {
                                cit->Read( lkey, ta );
                                cit->Read( rkey, tb );
                                top( ta, tb );
                                cit->Store( lkey, ta );
                            }
                        });
                        break;

                    case ELEMENT_INTEGRATION_POINT:
                    case FACE_INTEGRATION_POINT:
                    case INTER_FACE_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                                {
                                    cit->Read( ip, lkey, ta );
                                    cit->Read( ip, rkey, tb );
                                    top( ta, tb );
                                    cit->Store( ip, lkey, ta );
                                }
                        });
                        break;

                    case SECTOR_INTEGRATION_POINT:
                    case FACE_SECTOR_INTEGRATION_POINT:
                    case INTER_FACE_SECTOR_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                                {
                                    cit->Read( s, 0U, lkey, ta );
                                    cit->Read( s, 0U, rkey, tb );
                                    top( ta, tb );
                                    cit->Store( s, 0U, lkey, ta );
                                }
                        });
                        break;

                    case FACET_INTEGRATION_POINT:
                    case FACE_FACET_INTEGRATION_POINT:
                    case INTER_FACE_FACET_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                                {
                                    cit->Read( fac, 0U, lkey, ta );
                                    cit->Read( fac, 0U, rkey, tb );
                                    top( ta, tb );
                                    cit->Store( fac, 0U, lkey, ta );
                                }
                        });
                        break;

                    case REGION:
                        WithSubdomain( [&]( auto& sd )
                        {
                            sd.Read( lkey, ta );
                            sd.Read( rkey, tb );
                            top( ta, tb );
                            sd.Store( lkey, ta );
                        });
                        break;

                    default:
                        throw csmp::Exception( ERROR,
                            "PropertyHandle::ApplyBinaryOpSamePlacement",
                            var_name_.c_str(), "Unsupported placement." );
                }
                break;
            }

            case ARRAY:
            {
                ArrayVariable av1, av2;

                // Verify sizes match — programming error if they differ.
                WithSubdomain( [&]( auto& sd )
                {
                    if ( !sd.CellVector().empty() )
                    {
                        sd.CellVector().front()->Read( lkey, av1 );
                        sd.CellVector().front()->Read( rkey, av2 );
                        assert( av1.Size() == av2.Size()
                            && "PropertyHandle3::ApplyBinaryOpSamePlacement: "
                               "ArrayVariable size mismatch between operands." );
                    }
                });

                switch ( lkey.place )
                {
                    case NODE:
                        if ( !is_split_boundary_ )
                        {
                            WithSubdomain( [&]( auto& sd )
                            {
                                for ( auto& nit : sd.NodeVector() )
                                    if ( nit->Status( lkey ) == flag_output_ )
                                    {
                                        nit->Read( lkey, av1 );
                                        nit->Read( rkey, av2 );
                                        for ( uint32_t i = 0; i < av1.Size(); ++i )
                                            sop( av1(i), av2[i] );
                                        nit->Store( lkey, av1 );
                                    }
                            });
                        }
                        else
                        {
                            std::unordered_set<Node<dim>*> visited;
                            for ( auto& ifit : split_boundary_->CellVector() )
                                for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                                    for ( auto side : { INSIDE, OUTSIDE } )
                                    {
                                        auto* nd = ifit->N( i, side );
                                        if ( visited.insert( nd ).second &&
                                             nd->Status( lkey ) == flag_output_ )
                                        {
                                            nd->Read( lkey, av1 );
                                            nd->Read( rkey, av2 );
                                            for ( uint32_t j = 0; j < av1.Size(); ++j )
                                                sop( av1(j), av2[j] );
                                            nd->Store( lkey, av1 );
                                        }
                                    }
                        }
                        break;

                    case ELEMENT:
                    case FACE:
                    case INTER_FACE:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                if ( cit->Status( lkey ) == flag_output_ )
                                {
                                    cit->Read( lkey, av1 );
                                    cit->Read( rkey, av2 );
                                    for ( uint32_t i = 0; i < av1.Size(); ++i )
                                        sop( av1(i), av2[i] );
                                    cit->Store( lkey, av1 );
                                }
                        });
                        break;

                    case ELEMENT_INTEGRATION_POINT:
                    case FACE_INTEGRATION_POINT:
                    case INTER_FACE_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t ip = 0;
                                      ip < cit->IntegrationPoints(); ++ip )
                                    if ( cit->Status( ip, lkey ) == flag_output_ )
                                    {
                                        cit->Read( ip, lkey, av1 );
                                        cit->Read( ip, rkey, av2 );
                                        for ( uint32_t i = 0; i < av1.Size(); ++i )
                                            sop( av1(i), av2[i] );
                                        cit->Store( ip, lkey, av1 );
                                    }
                        });
                        break;

                    case SECTOR_INTEGRATION_POINT:
                    case FACE_SECTOR_INTEGRATION_POINT:
                    case INTER_FACE_SECTOR_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                                    if ( cit->Status( s, 0U, lkey ) == flag_output_ )
                                    {
                                        cit->Read( s, 0U, lkey, av1 );
                                        cit->Read( s, 0U, rkey, av2 );
                                        for ( uint32_t i = 0; i < av1.Size(); ++i )
                                            sop( av1(i), av2[i] );
                                        cit->Store( s, 0U, lkey, av1 );
                                    }
                        });
                        break;

                    case FACET_INTEGRATION_POINT:
                    case FACE_FACET_INTEGRATION_POINT:
                    case INTER_FACE_FACET_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                                    if ( cit->Status( fac, 0U, lkey ) == flag_output_ )
                                    {
                                        cit->Read( fac, 0U, lkey, av1 );
                                        cit->Read( fac, 0U, rkey, av2 );
                                        for ( uint32_t i = 0; i < av1.Size(); ++i )
                                            sop( av1(i), av2[i] );
                                        cit->Store( fac, 0U, lkey, av1 );
                                    }
                        });
                        break;

                    case REGION:
                        WithSubdomain( [&]( auto& sd )
                        {
                            if ( sd.Status( lkey ) == flag_output_ )
                            {
                                sd.Read( lkey, av1 );
                                sd.Read( rkey, av2 );
                                for ( uint32_t i = 0; i < av1.Size(); ++i )
                                    sop( av1(i), av2[i] );
                                sd.Store( lkey, av1 );
                            }
                        });
                        break;

                    default:
                        throw csmp::Exception( ERROR,
                            "PropertyHandle3::ApplyBinaryOpSamePlacement",
                            var_name_.c_str(), "Unsupported placement." );
                }
                break;
            }

            case FLAGGEDARRAY:
            {
                FlaggedArrayVariable fav1, fav2;

                // Verify sizes match.
                WithSubdomain( [&]( auto& sd )
                {
                    if ( !sd.CellVector().empty() )
                    {
                        sd.CellVector().front()->Read( lkey, fav1 );
                        sd.CellVector().front()->Read( rkey, fav2 );
                        assert( fav1.Size() == fav2.Size()
                            && "PropertyHandle3::ApplyBinaryOpSamePlacement: "
                               "FlaggedArrayVariable size mismatch." );
                    }
                });

                switch ( lkey.place )
                {
                    case NODE:
                        if ( !is_split_boundary_ )
                        {
                            WithSubdomain( [&]( auto& sd )
                            {
                                for ( auto& nit : sd.NodeVector() )
                                {
                                    nit->Read( lkey, fav1 );
                                    nit->Read( rkey, fav2 );
                                    for ( uint32_t i = 0; i < fav1.Size(); ++i )
                                        if ( fav1.Flag(i) == flag_output_ )
                                            sop( fav1(i), fav2[i] );
                                    nit->Store( lkey, fav1 );
                                }
                            });
                        }
                        else
                        {
                            std::unordered_set<Node<dim>*> visited;
                            for ( auto& ifit : split_boundary_->CellVector() )
                                for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                                    for ( auto side : { INSIDE, OUTSIDE } )
                                    {
                                        auto* nd = ifit->N( i, side );
                                        if ( visited.insert( nd ).second )
                                        {
                                            nd->Read( lkey, fav1 );
                                            nd->Read( rkey, fav2 );
                                            for ( uint32_t j = 0; j < fav1.Size(); ++j )
                                                if ( fav1.Flag(j) == flag_output_ )
                                                    sop( fav1(j), fav2[j] );
                                            nd->Store( lkey, fav1 );
                                        }
                                    }
                        }
                        break;

                    case ELEMENT:
                    case FACE:
                    case INTER_FACE:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                            {
                                cit->Read( lkey, fav1 );
                                cit->Read( rkey, fav2 );
                                for ( uint32_t i = 0; i < fav1.Size(); ++i )
                                    if ( fav1.Flag(i) == flag_output_ )
                                        sop( fav1(i), fav2[i] );
                                cit->Store( lkey, fav1 );
                            }
                        });
                        break;

                    case ELEMENT_INTEGRATION_POINT:
                    case FACE_INTEGRATION_POINT:
                    case INTER_FACE_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t ip = 0;
                                      ip < cit->IntegrationPoints(); ++ip )
                                {
                                    cit->Read( ip, lkey, fav1 );
                                    cit->Read( ip, rkey, fav2 );
                                    for ( uint32_t i = 0; i < fav1.Size(); ++i )
                                        if ( fav1.Flag(i) == flag_output_ )
                                            sop( fav1(i), fav2[i] );
                                    cit->Store( ip, lkey, fav1 );
                                }
                        });
                        break;

                    case SECTOR_INTEGRATION_POINT:
                    case FACE_SECTOR_INTEGRATION_POINT:
                    case INTER_FACE_SECTOR_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                                {
                                    cit->Read( s, 0U, lkey, fav1 );
                                    cit->Read( s, 0U, rkey, fav2 );
                                    for ( uint32_t i = 0; i < fav1.Size(); ++i )
                                        if ( fav1.Flag(i) == flag_output_ )
                                            sop( fav1(i), fav2[i] );
                                    cit->Store( s, 0U, lkey, fav1 );
                                }
                        });
                        break;

                    case FACET_INTEGRATION_POINT:
                    case FACE_FACET_INTEGRATION_POINT:
                    case INTER_FACE_FACET_INTEGRATION_POINT:
                        WithSubdomain( [&]( auto& sd )
                        {
                            for ( auto& cit : sd.CellVector() )
                                for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                                {
                                    cit->Read( fac, 0U, lkey, fav1 );
                                    cit->Read( fac, 0U, rkey, fav2 );
                                    for ( uint32_t i = 0; i < fav1.Size(); ++i )
                                        if ( fav1.Flag(i) == flag_output_ )
                                            sop( fav1(i), fav2[i] );
                                    cit->Store( fac, 0U, lkey, fav1 );
                                }
                        });
                        break;

                    case REGION:
                        WithSubdomain( [&]( auto& sd )
                        {
                            sd.Read( lkey, fav1 );
                            sd.Read( rkey, fav2 );
                            for ( uint32_t i = 0; i < fav1.Size(); ++i )
                                if ( fav1.Flag(i) == flag_output_ )
                                    sop( fav1(i), fav2[i] );
                            sd.Store( lkey, fav1 );
                        });
                        break;

                    default:
                        throw csmp::Exception( ERROR,
                            "PropertyHandle3::ApplyBinaryOpSamePlacement",
                            var_name_.c_str(), "Unsupported placement." );
                }
                break;
            }

            default:
                throw csmp::Exception( ERROR,
                    "PropertyHandle::ApplyBinaryOpSamePlacement",
                    var_name_.c_str(), "Unknown variable type." );
        }
    }
    // ------------------------------------------------------------------
    //  Mixed type: VECTOR op SCALAR
    // ------------------------------------------------------------------
    else if ( lkey.type == VECTOR && rkey.type == SCALAR )
    {
        VectorVariable<dim> va;
        switch ( lkey.place )
        {
            case NODE:
                if ( !is_split_boundary_ )
                {
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& nit : sd.NodeVector() )
                        {
                            nit->Read( lkey, va );
                            double b = nit->Read( rkey );
                            for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                            nit->Store( lkey, va );
                        }
                    });
                }
                else
                {
                    std::unordered_set<Node<dim>*> visited;
                    for ( auto& ifit : split_boundary_->CellVector() )
                        for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                        {
                            for ( auto side : { INSIDE, OUTSIDE } )
                            {
                                auto* nd = ifit->N( i, side );
                                if ( visited.insert( nd ).second )
                                {
                                    nd->Read( lkey, va );
                                    double b = nd->Read( rkey );
                                    for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                                    nd->Store( lkey, va );
                                }
                            }
                            if ( ifit->HasInterveningElement() )
                            {
                                auto* nd = ifit->N( i, MIDDLE );
                                if ( visited.insert( nd ).second )
                                {
                                    nd->Read( lkey, va );
                                    double b = nd->Read( rkey );
                                    for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                                    nd->Store( lkey, va );
                                }
                            }
                        }
                }
                break;

            case ELEMENT:
            case FACE:
            case INTER_FACE:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                    {
                        cit->Read( lkey, va );
                        double b = cit->Read( rkey );
                        for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                        cit->Store( lkey, va );
                    }
                });
                break;

            case ELEMENT_INTEGRATION_POINT:
            case FACE_INTEGRATION_POINT:
            case INTER_FACE_INTEGRATION_POINT:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                        for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        {
                            cit->Read( ip, lkey, va );
                            double b = cit->Read( ip, rkey );
                            for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                            cit->Store( ip, lkey, va );
                        }
                });
                break;

            case SECTOR_INTEGRATION_POINT:
            case FACE_SECTOR_INTEGRATION_POINT:
            case INTER_FACE_SECTOR_INTEGRATION_POINT:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                        for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        {
                            cit->Read( s, 0U, lkey, va );
                            double b = cit->Read( s, 0U, rkey );
                            for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                            cit->Store( s, 0U, lkey, va );
                        }
                });
                break;

            case FACET_INTEGRATION_POINT:
            case FACE_FACET_INTEGRATION_POINT:
            case INTER_FACE_FACET_INTEGRATION_POINT:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                        for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        {
                            cit->Read( fac, 0U, lkey, va );
                            double b = cit->Read( fac, 0U, rkey );
                            for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                            cit->Store( fac, 0U, lkey, va );
                        }
                });
                break;

            case REGION:
                WithSubdomain( [&]( auto& sd )
                {
                    sd.Read( lkey, va );
                    double b = sd.Read( rkey );
                    for ( uint32_t j = 0; j < dim; ++j ) sop( va(j), b );
                    sd.Store( lkey, va );
                });
                break;

            default:
                throw csmp::Exception( ERROR,
                    "PropertyHandle::ApplyBinaryOpSamePlacement",
                    var_name_.c_str(), "Unsupported placement." );
        }
    }
    // ------------------------------------------------------------------
    //  Mixed type: TENSOR op SCALAR
    // ------------------------------------------------------------------
    else if ( lkey.type == TENSOR && rkey.type == SCALAR )
    {
        TensorVariable<dim> ta;
        switch ( lkey.place )
        {
            case NODE:
                if ( !is_split_boundary_ )
                {
                    WithSubdomain( [&]( auto& sd )
                    {
                        for ( auto& nit : sd.NodeVector() )
                        {
                            nit->Read( lkey, ta );
                            double b = nit->Read( rkey );
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j ) sop( ta(i,j), b );
                            nit->Store( lkey, ta );
                        }
                    });
                }
                else
                {
                    std::unordered_set<Node<dim>*> visited;
                    for ( auto& ifit : split_boundary_->CellVector() )
                        for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                        {
                            for ( auto side : { INSIDE, OUTSIDE } )
                            {
                                auto* nd = ifit->N( i, side );
                                if ( visited.insert( nd ).second )
                                {
                                    nd->Read( lkey, ta );
                                    double b = nd->Read( rkey );
                                    for ( uint32_t ii = 0; ii < dim; ++ii )
                                        for ( uint32_t j = 0; j < dim; ++j ) sop( ta(ii,j), b );
                                    nd->Store( lkey, ta );
                                }
                            }
                            if ( ifit->HasInterveningElement() )
                            {
                                auto* nd = ifit->N( i, MIDDLE );
                                if ( visited.insert( nd ).second )
                                {
                                    nd->Read( lkey, ta );
                                    double b = nd->Read( rkey );
                                    for ( uint32_t ii = 0; ii < dim; ++ii )
                                        for ( uint32_t j = 0; j < dim; ++j ) sop( ta(ii,j), b );
                                    nd->Store( lkey, ta );
                                }
                            }
                        }
                }
                break;

            case ELEMENT:
            case FACE:
            case INTER_FACE:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                    {
                        cit->Read( lkey, ta );
                        double b = cit->Read( rkey );
                        for ( uint32_t i = 0; i < dim; ++i )
                            for ( uint32_t j = 0; j < dim; ++j ) sop( ta(i,j), b );
                        cit->Store( lkey, ta );
                    }
                });
                break;

            case ELEMENT_INTEGRATION_POINT:
            case FACE_INTEGRATION_POINT:
            case INTER_FACE_INTEGRATION_POINT:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                        for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        {
                            cit->Read( ip, lkey, ta );
                            double b = cit->Read( ip, rkey );
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j ) sop( ta(i,j), b );
                            cit->Store( ip, lkey, ta );
                        }
                });
                break;

            case SECTOR_INTEGRATION_POINT:
            case FACE_SECTOR_INTEGRATION_POINT:
            case INTER_FACE_SECTOR_INTEGRATION_POINT:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                        for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        {
                            cit->Read( s, 0U, lkey, ta );
                            double b = cit->Read( s, 0U, rkey );
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j ) sop( ta(i,j), b );
                            cit->Store( s, 0U, lkey, ta );
                        }
                });
                break;

            case FACET_INTEGRATION_POINT:
            case FACE_FACET_INTEGRATION_POINT:
            case INTER_FACE_FACET_INTEGRATION_POINT:
                WithSubdomain( [&]( auto& sd )
                {
                    for ( auto& cit : sd.CellVector() )
                        for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        {
                            cit->Read( fac, 0U, lkey, ta );
                            double b = cit->Read( fac, 0U, rkey );
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j ) sop( ta(i,j), b );
                            cit->Store( fac, 0U, lkey, ta );
                        }
                });
                break;

            case REGION:
                WithSubdomain( [&]( auto& sd )
                {
                    sd.Read( lkey, ta );
                    double b = sd.Read( rkey );
                    for ( uint32_t i = 0; i < dim; ++i )
                        for ( uint32_t j = 0; j < dim; ++j ) sop( ta(i,j), b );
                    sd.Store( lkey, ta );
                });
                break;

            default:
                throw csmp::Exception( ERROR,
                    "PropertyHandle::ApplyBinaryOpSamePlacement",
                    var_name_.c_str(), "Unsupported placement." );
        }
    }
    else
    {
        throw csmp::Exception( ERROR,
            "PropertyHandle::ApplyBinaryOpSamePlacement",
            var_name_.c_str(),
            "Unsupported type combination for binary operation." );
    }
}




// ============================================================================
//  ApplyBinaryOpAligned
// ============================================================================

template<uint32_t dim>
template<typename SF, typename VF, typename TF>
inline PropertyHandle<dim>& PropertyHandle<dim>::ApplyBinaryOpAligned(
    const PropertyHandle<dim>& other,
    SF&& sop,
    VF&& vop,
    TF&& top )
{
    const bool same_subdomain = ( subdomain_name_ == other.subdomain_name_ );
    const bool other_is_model = ( other.subdomain_name_ == "Model" );

    // The only valid cross-subdomain assignment is:
    //   whole-model NODE → Boundary NODE
    //   whole-model NODE → SplitBoundary NODE
    // All other cross-subdomain combinations are invalid because
    // FACE, INTER_FACE and their IP variants only exist on
    // Boundary and SplitBoundary respectively, and ELEMENT,
    // ELEMENT_INTEGRATION_POINT etc. only exist on Region.

    if ( !same_subdomain && other_is_model )
    {
        if ( key_.place != NODE || other.key_.place != NODE )
            throw Exception( ERROR,
                "PropertyHandle::operator=(PropertyHandle)",
                var_name_.c_str(),
                ( std::string(
                    "Cross-subdomain assignment is only supported between "
                    "NODE-placed variables. "
                    "FACE and INTER_FACE placements exist exclusively on "
                    "Boundary and SplitBoundary subdomains respectively. "
                    "ELEMENT and integration point placements exist "
                    "exclusively on Region subdomains. "
                    "Cannot assign '"
                    + std::string( other.VariableName() )
                    + "' (placement "
                    + std::to_string( other.key_.place )
                    + ") into '"
                    + var_name_
                    + "' (placement "
                    + std::to_string( key_.place )
                    + ")." )
                ).c_str() );
    }

    const PLACEMENT dst = key_.place;
    const PLACEMENT src = other.key_.place;

    if ( dst == src )
    {
        ApplyBinaryOpSamePlacement( other,
            std::forward<SF>( sop ),
            std::forward<VF>( vop ),
            std::forward<TF>( top ) );
        IsWithinRange();
        return *this;
    }

    // ------------------------------------------------------------------
    //  Placements differ — create a temporary property on the
    //  destination placement, interpolate into it, operate, then
    //  delete it immediately.
    // ------------------------------------------------------------------
    static std::atomic<uint32_t> tmp_counter{0};
    const std::string tmp_name = var_name_
                               + "__ph3_tmp_"
                               + std::to_string( tmp_counter.fetch_add(1) );

    model_.CreateProperty( tmp_name.c_str(), tmp_name.c_str(), "SI",
                           other.key_.type, dst, 1U );
    try
    {
        Interpolate( other.VariableName(), tmp_name.c_str(), src, dst );

        // Construct a non-owning view onto the temporary.
        // Because the property already exists when this constructor
        // runs, Initialise sets owns_variable_ = false, so its
        // destructor will NOT call DeleteProperty.
        // We delete explicitly below.
        PropertyHandle<dim> aligned_view( model_,
                                           subdomain_name_.c_str(),
                                           tmp_name.c_str(),
                                           other.key_.type,
                                           dst );

        ApplyBinaryOpSamePlacement( aligned_view,
            std::forward<SF>( sop ),
            std::forward<VF>( vop ),
            std::forward<TF>( top ) );
    }
    catch ( ... )
    {
        model_.DeleteProperty( tmp_name.c_str() );
        throw;
    }

    model_.DeleteProperty( tmp_name.c_str() );
    IsWithinRange();
    return *this;
}


template<uint32_t dim>
template<typename ScalarOp, typename ArrayOp, typename FlaggedArrayOp>
inline void PropertyHandle<dim>::ApplyMathOp( ScalarOp&& scalar_op, ArrayOp&& array_op, FlaggedArrayOp&& farray_op )
{
    switch ( key_.type )
    {
        case ARRAY:
            ApplyArray( std::forward<ArrayOp>( array_op ) );
            break;
        case FLAGGEDARRAY:
            ApplyFlaggedArray( std::forward<FlaggedArrayOp>( farray_op ) );
            break;
        default:
            // SCALAR, VECTOR, TENSOR — Apply dispatches via applyByType.
            Apply( std::forward<ScalarOp>( scalar_op ) );
            break;
    }
    IsWithinRange();
}


// ============================================================================
//  Named math methods (all inline, all delegate to Apply)
// ============================================================================

template<uint32_t dim> inline void PropertyHandle<dim>::Squared()
{
    // Component-wise squaring for all types.
    // For TENSOR this is the Hadamard product T_ij^2, NOT the matrix
    // product T*T. Use MatrixSquared() for the matrix product.
    ApplyMathOp(
        []( double& x )                     { x *= x; },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i ) av(i) *= av(i); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x *= x; } ); }
    );
    // Note: for VECTOR and TENSOR, Apply(f(double&)) applies f
    // component-wise via applyVector/applyTensor, so this is
    // always Hadamard squaring regardless of type.
}

template<uint32_t dim> inline void PropertyHandle<dim>::MatrixSquared()
{
    if ( key_.type != TENSOR )
        throw csmp::Exception( ERROR,
            "PropertyHandle3::MatrixSquared",
            var_name_.c_str(),
            "MatrixSquared is only defined for TENSOR handles. "
            "Use Squared() for component-wise squaring." );
    ApplyTensor( []( TensorVariable<dim>& ts ) { ts *= ts; } );
    IsWithinRange();
}


template<uint32_t dim>
inline void PropertyHandle<dim>::DoubleContraction( const PropertyHandle<dim>& tensor )
{
    if ( key_.type != SCALAR )
        throw csmp::Exception( ERROR,
            "PropertyHandle3::DoubleContraction",
            var_name_.c_str(),
            "DoubleContraction requires this handle to be SCALAR "
            "(to receive the result)." );
    if ( tensor.key_.type != TENSOR )
        throw csmp::Exception( ERROR,
            "PropertyHandle3::DoubleContraction",
            var_name_.c_str(),
            "DoubleContraction requires the argument handle to be TENSOR." );
    if ( key_.place != tensor.key_.place )
        throw csmp::Exception( ERROR,
            "PropertyHandle3::DoubleContraction",
            var_name_.c_str(),
            "DoubleContraction requires both handles to have the same "
            "placement." );

    // Use applyBinaryOpAligned with a scalar functor that reads the
    // tensor and computes T:T = sum_ij T_ij^2.
    // We cannot use the standard binary op path since the types differ,
    // so we iterate manually.
    const csmp::Index& skey = key_;
    const csmp::Index& tkey = tensor.key_;

    WithSubdomain( [&]( auto& sd )
    {
        TensorVariable<dim> ts;
        switch ( skey.place )
        {
            case ELEMENT:
            case FACE:
            case INTER_FACE:
                for ( auto& cit : sd.CellVector() )
                    if ( cit->Status( skey ) == flag_output_ )
                    {
                        cit->Read( tkey, ts );
                        double sum = 0.0;
                        for ( uint32_t i = 0; i < dim; ++i )
                            for ( uint32_t j = 0; j < dim; ++j )
                                sum += ts(i,j) * ts(i,j);
                        cit->Store( skey, makeScalar( flag_output_, sum ) );
                    }
                break;

            case NODE:
                if ( !is_split_boundary_ )
                {
                    for ( auto& nit : sd.NodeVector() )
                        if ( nit->Status( skey ) == flag_output_ )
                        {
                            nit->Read( tkey, ts );
                            double sum = 0.0;
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j )
                                    sum += ts(i,j) * ts(i,j);
                            nit->Store( skey, makeScalar( flag_output_, sum ) );
                        }
                }
                else
                {
                    std::unordered_set<Node<dim>*> visited;
                    for ( auto& ifit : split_boundary_->CellVector() )
                        for ( uint32_t i = 0; i < ifit->FE()->Nodes(); ++i )
                            for ( auto side : { INSIDE, OUTSIDE } )
                            {
                                auto* nd = ifit->N( i, side );
                                if ( visited.insert( nd ).second &&
                                     nd->Status( skey ) == flag_output_ )
                                {
                                    nd->Read( tkey, ts );
                                    double sum = 0.0;
                                    for ( uint32_t ii = 0; ii < dim; ++ii )
                                        for ( uint32_t jj = 0; jj < dim; ++jj )
                                            sum += ts(ii,jj) * ts(ii,jj);
                                    nd->Store( skey, makeScalar( flag_output_, sum ) );
                                }
                            }
                }
                break;

            case ELEMENT_INTEGRATION_POINT:
            case FACE_INTEGRATION_POINT:
            case INTER_FACE_INTEGRATION_POINT:
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t ip = 0; ip < cit->IntegrationPoints(); ++ip )
                        if ( cit->Status( ip, skey ) == flag_output_ )
                        {
                            cit->Read( ip, tkey, ts );
                            double sum = 0.0;
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j )
                                    sum += ts(i,j) * ts(i,j);
                            cit->Store( ip, skey, makeScalar( flag_output_, sum ) );
                        }
                break;

            case SECTOR_INTEGRATION_POINT:
            case FACE_SECTOR_INTEGRATION_POINT:
            case INTER_FACE_SECTOR_INTEGRATION_POINT:
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t s = 0; s < cit->Sectors(); ++s )
                        if ( cit->Status( s, 0U, skey ) == flag_output_ )
                        {
                            cit->Read( s, 0U, tkey, ts );
                            double sum = 0.0;
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j )
                                    sum += ts(i,j) * ts(i,j);
                            cit->Store( s, 0U, skey, makeScalar( flag_output_, sum ) );
                        }
                break;

            case FACET_INTEGRATION_POINT:
            case FACE_FACET_INTEGRATION_POINT:
            case INTER_FACE_FACET_INTEGRATION_POINT:
                for ( auto& cit : sd.CellVector() )
                    for ( uint32_t fac = 0; fac < cit->Facets(); ++fac )
                        if ( cit->Status( fac, 0U, skey ) == flag_output_ )
                        {
                            cit->Read( fac, 0U, tkey, ts );
                            double sum = 0.0;
                            for ( uint32_t i = 0; i < dim; ++i )
                                for ( uint32_t j = 0; j < dim; ++j )
                                    sum += ts(i,j) * ts(i,j);
                            cit->Store( fac, 0U, skey, makeScalar( flag_output_, sum ) );
                        }
                break;

            case REGION:
                if ( sd.Status( skey ) == flag_output_ )
                {
                    sd.Read( tkey, ts );
                    double sum = 0.0;
                    for ( uint32_t i = 0; i < dim; ++i )
                        for ( uint32_t j = 0; j < dim; ++j )
                            sum += ts(i,j) * ts(i,j);
                    sd.Store( skey, makeScalar( flag_output_, sum ) );
                }
                break;

            default:
                throw csmp::Exception( ERROR,
                    "PropertyHandle3::DoubleContraction",
                    var_name_.c_str(), "Unsupported placement." );
        }
    });
    IsWithinRange();
}


template<uint32_t dim> inline void PropertyHandle<dim>::Sqrt()
{
    ApplyMathOp(
        []( double& x )                     { x = std::sqrt(x); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i ) av(i) = std::sqrt(av(i)); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = std::sqrt(x); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Ln()
{
    ApplyMathOp(
        []( double& x )                     { if (x>0.0) x = std::log(x); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  if (av(i)>0.0) av(i) = std::log(av(i)); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { if (x>0.0) x = std::log(x); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Log10()
{
    ApplyMathOp(
        []( double& x )                     { if (x>0.0) x = std::log10(x); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  if (av(i)>0.0) av(i) = std::log10(av(i)); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { if (x>0.0) x = std::log10(x); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Exp()
{
    ApplyMathOp(
        []( double& x )                     { x = std::exp(x); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i ) av(i) = std::exp(av(i)); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = std::exp(x); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Pow( double exponent )
{
    ApplyMathOp(
        [exponent]( double& x )                     { x = std::pow(x, exponent); },
        [exponent]( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                          av(i) = std::pow(av(i), exponent); },
        [exponent, this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                          [exponent]( double& x ) { x = std::pow(x, exponent); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::ZapNAN( double replacement )
{
    ApplyMathOp(
        [replacement]( double& x )                     { if (std::isnan(x)) x = replacement; },
        [replacement]( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                             if (std::isnan(av(i))) av(i) = replacement; },
        [replacement, this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                             [replacement]( double& x ) { if (std::isnan(x)) x = replacement; } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Abs()
{
    ApplyMathOp(
        []( double& x )                     { x = std::abs(x); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i ) av(i) = std::abs(av(i)); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = std::abs(x); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Clamp( double lo, double hi )
{
    if ( lo > hi )
        throw csmp::Exception( ERROR, "PropertyHandle3::Clamp",
            var_name_.c_str(), "lo must be <= hi." );
    ApplyMathOp(
        [lo,hi]( double& x )                     { x = std::max(lo, std::min(x, hi)); },
        [lo,hi]( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                       av(i) = std::max(lo, std::min(av(i), hi)); },
        [lo,hi,this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                      [lo,hi]( double& x ) { x = std::max(lo, std::min(x, hi)); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Sin()
{
    ApplyMathOp(
        []( double& x )                     { x = std::sin( degreesToRadians(x) ); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  av(i) = std::sin( degreesToRadians(av(i)) ); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = std::sin( degreesToRadians(x) ); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Cos()
{
    ApplyMathOp(
        []( double& x )                     { x = std::cos( degreesToRadians(x) ); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  av(i) = std::cos( degreesToRadians(av(i)) ); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = std::cos( degreesToRadians(x) ); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Tan()
{
    ApplyMathOp(
        []( double& x )                     { x = std::tan( degreesToRadians(x) ); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  av(i) = std::tan( degreesToRadians(av(i)) ); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = std::tan( degreesToRadians(x) ); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Acos()
{
    // Range check only meaningful for SCALAR — for VECTOR, TENSOR,
    // ARRAY and FLAGGEDARRAY the range is over lengths or components
    // which may legitimately exceed [-1,1] even when individual
    // component values are valid inputs to acos.
    if ( key_.type == SCALAR )
    {
        double omin, omax;
        Range( omin, omax );
        if ( omin < -1.0 || omax > 1.0 )
            throw csmp::Exception( ERROR, "PropertyHandle3::Acos",
                var_name_.c_str(), "Values must lie in [-1,1]." );
    }
    ApplyMathOp(
        []( double& x )                     { x = radiansToDegrees( std::acos(x) ); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  av(i) = radiansToDegrees( std::acos(av(i)) ); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = radiansToDegrees( std::acos(x) ); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Asin()
{
    if ( key_.type == SCALAR )
    {
        double omin, omax;
        Range( omin, omax );
        if ( omin < -1.0 || omax > 1.0 )
            throw csmp::Exception( ERROR, "PropertyHandle3::Asin",
                var_name_.c_str(), "Values must lie in [-1,1]." );
    }
    ApplyMathOp(
        []( double& x )                     { x = radiansToDegrees( std::asin(x) ); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  av(i) = radiansToDegrees( std::asin(av(i)) ); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = radiansToDegrees( std::asin(x) ); } ); }
    );
}

template<uint32_t dim> inline void PropertyHandle<dim>::Atan()
{
    if ( key_.type == SCALAR )
    {
        double omin, omax;
        Range( omin, omax );
        if ( omin < -1.0 || omax > 1.0 )
            throw csmp::Exception( ERROR, "PropertyHandle3::Atan",
                var_name_.c_str(), "Values must lie in [-1,1]." );
    }
    ApplyMathOp(
        []( double& x )                     { x = radiansToDegrees( std::atan(x) ); },
        []( ArrayVariable& av )             { for ( uint32_t i=0; i<av.Size(); ++i )
                                                  av(i) = radiansToDegrees( std::atan(av(i)) ); },
        [this]( FlaggedArrayVariable& fav ) { ApplyToFlaggedElements( fav, flag_output_,
                                                []( double& x ) { x = radiansToDegrees( std::atan(x) ); } ); }
    );
}


} // namespace csmp

#endif // CSMP_PROPERTY_HANDLE3_H

