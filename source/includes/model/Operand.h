#ifndef CSMP_OPERAND_H
#define CSMP_OPERAND_H

#include "Index.h"
#include "TensorVariable.h"
#include "ArrayVariable.h"
#include "FlaggedArrayVariable.h"
#include "ErrorHandler.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;

/**
 
@brief Interface for all of CSMP's variable types facilitating that they
can be used together in computations; currently used only by Interrelation; outdated design, refactoring required.

@author S.K. Matthai
@author Stephen G. Roberts
@date 1999
 
@section motivation Motivation

Forget for a moment all about the rest of the CSMP code except for the basic 
variables ScalarVariable, VectorVariable and TensorVariable and think of 
the possibility of creating objects that are overloaded such that you can
write expressions for the whole mesh object or sub-regions thereof, which 
read like mathematical equations in which the variables refer to local
or distributed data.  

This quest was the design intent of the Operand class and it brings CSP
users close to this goal when they are prepared to write their own
Interrelation subclasses. The goal is fully realized through the new
PropertyHandle classes for use in main() programs. These are however only 
complementary to the Operand and Interrelation mechanism, because they
do not encapsulate the calculations which Interrelations do. 

 
@section design Design Intent

Apart from being able to write mathematical expressions which are applied to
the whole Model or sub-Regions thereof (as is facilitated through
the PropertyHandle class), there is a need to use distributed instances of 
properties in calculations inside objects which are still efficient and 
allow the user to partially ignore variable IDs and the variable placement 
in the mesh (Node, IntegrationPoint or Element).  This requires (1) the use 
of a notation to represent the variable in mathematical expressions,
and (2) an association of the csmp::Index of the variable with the notation such
that it does not need to be requested from the PropertyDatabase every time
the expression is evaluated. Finally (3), temporary variable storage is 
needed to allow that a value can be range-checked before
it is written back into the Model storage. This buffer is essential 
when it comes to checking where in a chain of expressions, a physically 
meaningless value was produced.
 
 
@section applicability Applicability

Operands are used in Interrelations, Algorithms, Algorithms, and
Visitor objects.
 
 
@section structure Structure

The Operand is a simple class object with a String containing its
name, a csmp::Index describing the represented variable, and 
each one of the basic CSMP variables as well as an output condition
flag which specifies the VARIABLE_FLAG of the property locations which
the variable is permitted to overwrite if the calculation yielded a 
result which is within the permitted range.
 
@section participants Participants
 
A csmp::Index, a String, an output VARIABLE_FLAG and CSP basic variables
(ScalarVariable, VectorVariable, and TensorVariable).
 
 
@section collaborations Collaborations

Operands collaborate with CSP basic variables and the
enumerations which are defined in CSMP_definitions.h (e.g., VARIABLE_FLAG...).
 
 
@section consequences Consequences
 

The operators \f$+, -, *, /, +=, -=, *=, /= \f$ are overloaded for Operands such
that they can interact with eachother, double data, and the CSP basic
variables. One should not write expressions like A = C. Following the logic
of an assignment, A = B means: The Operand A now 
refers to exactly the same variable as B and it will also contain the same 
values. If Operand A represents a property stored in the Model, the 
variable A will no longer exist following the assignment.  

Importantly, scalar-, vector-, and tensor variables are more than 
just mathematical variables in that they also contain a STATUS flag which 
specifies how they are treated in the assembly of the finite-element 
equations. This requires some rules about the flag assignment: If you say
 
@code
A = C + B;
@endcode

the temporary variable that is created will be flagged like the lefthand 
operator before the assigment. If \f$C = PLAIN\f$ and \f$B=DIRICH\f$, A will be
\f$PLAIN\f$. If the combined operator += is used, A will retain its original flag
because the assignment here uses the assignment as defined for scalar
variables.  

If you compare an Operand with a VARIABLE_FLAG enumeration, the comparison
will be carried out using all flag(s) of the underlying variable. Thus, 
if a vector variable in a 2D calculation has two different flags (for instance, 
a DIRICH flag for the X displacement and a PLAIN flag for the Y displacement) 
a comparison will evaluate as false. In this case one should use 
VectorVariable::Flag( int i) to carry out a more specific comparison. 
Any comparisons with not initialized Operands or Operands that 
store NAN values in the variables will evaluate as false. 
 
 
@section implementation Implementation

By default, Operands have storage for a single variable value of type 
ScalarVariable, VectorVariable, or TensorVariable, but only if this variable 
exists in the PropertyDatabase (in contrast, PropertyHandles spawn the
generation of a new variable in the PropertyDatabase if the variable does 
not exist).   

The need for large computations imposes stringend constraints on the properties
of Operands. They may be big in size because there will only be a few, but their 
overloaded operators have to be efficient since single calculations may involve 
millions of calls to those. Already a function call overhead is significant for 
the overloaded operators. To avoid this overhead auto-inlining is used in the
compilation. 

The construction of temporary 'Operands' in expressions like:

@code
C = A + 32.5;
@endcode

should also be avoided. In those cases the construction by the compiler
of temporary objects (object x which will be assigned to C) can be avoided 
only in some cases through construction of the Operands in the return 
statements of the overloaded operators.  

To avoid the creation of temporaries you should rather write:  
 
@code
C = 32.5;\
C += A;
@endcode

Operands also store a 'calculation offset' which is used by Algorithms
to modify the placement of dependent variables in computations 
involving the Operand variable, when coupled systems of equations 
are solved. 

 
@section application Application Examples
 

Operands are typically used in Algorithms and Interrelations. When you
derive a subclass from an Interrelation, you first create a Operand 
reference to a physical variable in your model:

@code
Operand&  K;
@endcode

You then associate the Operand with a property which has been defined
before as a member of the PropertyDatabase. This happens in your 
interrelation constructor.
 
@code
K( GlobalProperty("conductivity") )
@endcode

This defines K as a notation for 'conductivity'. You can 
now use K in expressions using CSMP basic variables and other Operands:

@code
K += PHI / 322.1 - 45.0;
@endcode

*/
template<uint32_t dim>
class Operand {
  public:
    Operand();
    explicit Operand( const csmp::Index& );
    Operand( const Operand&  ); 
    Operand( const char* var, const PropertyDatabase<dim>& );
    Operand( const csmp::Index&, const ScalarVariable& ); 
    Operand( const csmp::Index&, const VectorVariable<dim>& );
    Operand( const csmp::Index&, const TensorVariable<dim>& );
    Operand( const csmp::Index&, const ArrayVariable&  );
    Operand( const csmp::Index&, const FlaggedArrayVariable&  );
    ~Operand();
    // these calculate using 'val' on right and object specs on the left
    Operand   operator+( double ) const;
    Operand   operator-( double ) const;
    Operand   operator*( double ) const;
    Operand   operator/( double ) const;
    // these will create copies of the specifications of object on the right 
    Operand   operator+( const Operand& ) const;
    Operand   operator-( const Operand& ) const;
    Operand   operator*( const Operand& ) const;
    Operand   operator/( const Operand& ) const;
    // these will "self-assign" calculation result to object on the left
    Operand&  operator+=( double );
    Operand&  operator-=( double );
    Operand&  operator*=( double );
    Operand&  operator/=( double );
    // with ScalarVariables
    Operand&  operator+=( const ScalarVariable& );
    Operand&  operator-=( const ScalarVariable& );
    Operand&  operator*=( const ScalarVariable& );
    Operand&  operator/=( const ScalarVariable& );
    // with VectorVariables
    Operand&  operator+=( const VectorVariable<dim>& );
    Operand&  operator-=( const VectorVariable<dim>& );
    Operand&  operator*=( const VectorVariable<dim>& );
    Operand&  operator/=( const VectorVariable<dim>& );
    // with TensorVariables
    Operand&  operator+=( const TensorVariable<dim>& );
    Operand&  operator-=( const TensorVariable<dim>& );
    Operand&  operator*=( const TensorVariable<dim>& );
    Operand&  operator/=( const TensorVariable<dim>& );
    // with ArrayVariables
    Operand&  operator+=( const ArrayVariable& );
    Operand&  operator-=( const ArrayVariable& );
    Operand&  operator*=( const ArrayVariable& );
    Operand&  operator/=( const ArrayVariable& );
    // with FlaggedArrayVariables
    Operand&  operator+=( const FlaggedArrayVariable& );
    Operand&  operator-=( const FlaggedArrayVariable& );
    Operand&  operator*=( const FlaggedArrayVariable& );
    Operand&  operator/=( const FlaggedArrayVariable& );
    // with Operands
    Operand&  operator+=( const Operand& );
    Operand&  operator-=( const Operand& );
    Operand&  operator*=( const Operand& );
    Operand&  operator/=( const Operand& );
    // logical
    bool       operator>(  double val ) const;
    bool       operator<(  double val ) const;
    bool       operator>=( double val ) const;
    bool       operator<=( double val ) const;
    bool       operator==( double val ) const;
    bool       operator!=( double val ) const;

    bool       operator==( const ScalarVariable& )      const;
    bool       operator!=( const ScalarVariable& )      const;
    bool       operator==( const VectorVariable<dim>& ) const;
    bool       operator!=( const VectorVariable<dim>& ) const;
    bool       operator==( const TensorVariable<dim>& ) const;
    bool       operator!=( const TensorVariable<dim>& ) const;
    bool       operator==( const ArrayVariable& )       const;
    bool       operator!=( const ArrayVariable& )       const;
    bool       operator==( const FlaggedArrayVariable& )const;
    bool       operator!=( const FlaggedArrayVariable& )const;

    bool       operator>(  Operand& ) const;
    bool       operator<(  Operand& ) const;
    bool       operator>=( Operand& ) const;
    bool       operator<=( Operand& ) const;
    bool       operator==( Operand& ) const;
    bool       operator!=( Operand& ) const;
    // assignment
    Operand&   operator=( double );
    Operand&   operator=( const ScalarVariable& );
    Operand&   operator=( const VectorVariable<dim>& );
    Operand&   operator=( const TensorVariable<dim>& );
    Operand&   operator=( const ArrayVariable& );
    Operand&   operator=( const FlaggedArrayVariable& );
    Operand&   operator=( const Operand& );
    
    // methods
    void                Name( const char* );
    std::string         Name()              const;
    size_t              Size( );
    bool                IsWithinRange()     const;
    template<class cspT> bool IsWithinRange( const cspT& val ) const;
    double            MaxValue()          const;
    double            MinValue()          const;
    
    const csmp::Index&  Key()               const;
    size_t              Index()             const;
    PLACEMENT           Placement()         const;
    VARIABLE_TYPE       Type()              const;
    size_t              CalculationOffset() const;
    VARIABLE_FLAG       OutputCondition()   const;
    VARIABLE_FLAG       EssentialCondition()const;
    void                CalculationOffset( size_t );
    void                OutputCondition( VARIABLE_FLAG );
    void                EssentialCondition( VARIABLE_FLAG );
    
    void                AssignTo( ScalarVariable& )        const;
    void                AssignTo( VectorVariable<dim>& )   const;
    void                AssignTo( TensorVariable<dim>& )   const;
    void                AssignTo( ArrayVariable& )         const;
    void                AssignTo( FlaggedArrayVariable& )  const;

    void                PrintValue()        const;
    void                Out()               const;

  private:
    std::string          name_;
    csmp::Index          prop_key_;
    double               omin, omax; 
    VARIABLE_FLAG        flag_essential_, flag_output_;
    size_t               calc_offset_;
    ScalarVariable       scalar_storage_;
    VectorVariable<dim>  vector_storage_;
    TensorVariable<dim>  tensor_storage_;
    ArrayVariable        array_storage_;
    FlaggedArrayVariable flagged_array_storage_;
};

} // csmp

#endif
