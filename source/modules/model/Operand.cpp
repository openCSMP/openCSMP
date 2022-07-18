#include "Operand.h"
#include "PropertyDatabase.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

/**
 
The default constructor of the Operand sets the csmp::Index type to SCALAR,
the placement to MODEL and the index to ULONG_MAX. The Operand will
have no name, its output flag will be PLAIN, the essential flag DIRICH,
and the value rangle will be  -1.0e+50 to 1.0e+50. The calculation offset
will be zero. 
*/
template<uint32_t dim>
Operand<dim>::Operand()
 : name_("undefined"),
   omin(-1.0e+50),
   omax(1.0e+50),
   flag_output_(PLAIN),
   flag_essential_(DIRICH),
   calc_offset_(0)
 {
 }




/**
 
Operands can be constructed efficiently using only a csmp::Index as argument.
In this case they will have no name, a variable range from -1.0e+50 to
1.0e+50, and a calculation offset of zero. The output flag will be PLAIN
and the essential flag will be DIRICH. 

@param idx name of the csmp::Index of the variable with which
the Operand shall be associated. This index must correspond to a physical
variable which is already known to the PropertyDatabase. 
*/
template<uint32_t dim>
Operand<dim>::Operand( const csmp::Index& idx )
 : name_("undefined"),
   prop_key_(idx),
   omin(-1.0e+50),
   omax(1.0e+50),
   flag_output_(PLAIN),
   flag_essential_(DIRICH),
   calc_offset_(0)
 {
 }


template<uint32_t dim>
Operand<dim>::Operand( const Operand<dim>& op ) 
 {
    *this = op;
 } 


/**
 
This constructor will fully initialize an Operand using the current
variable data from the PropertyDatabase. The Operand name will be that of
the supplied variable, the calculation offset will be zero. The output flag 
will be PLAIN and the essential flag will be DIRICH. 

@param var the name of the variable with which the Operand
shall be associated. This name must be already known to the PropertyDatabase
a reference of which is passed as a second argument to the constructor. 

@section messages Messages

If the property database cannot identify the variable by its name, it
will report an error and the Operand will not be initialized properly. 
*/
template<uint32_t dim>
Operand<dim>::Operand( const char* var, const PropertyDatabase<dim>& p )
 : name_(var),
   prop_key_(p.StorageKey(var)),
   flag_output_(PLAIN),
   flag_essential_(DIRICH),
   calc_offset_(0)
 {
    p.RangeOf ( var, omin, omax );
 }



/**
 
Efficient constructor for temporary operands needed in mathematical 
expressions like:

@code
Op_A = Op_B + Op_C;
@endcode

@param idx the name of the csmp::Index of the variable with which
the Operand shall be associated. This index must correspond to a physical
variable which is already known to the PropertyDatabase. 

@param sc an initialized CSP variable, the type of which must
be equivalent to the type specified in the csmp::Index> 

@section implementation Implementation

These seemingly 'weird' constructors have been created only because
they allow a smart compiler to create temporaries in the place of their
application (return value optimization).

The resulting temporary Operands have the name 'tmp' (also remember that, in 
large computations, one should write x = 0; x += Operand + Operand; 
if one can, to avoid the creation of temporaries all together. 

@section application Application

Temporary variables in mathematical expressions involving Operands. 
*/
template<uint32_t dim>
Operand<dim>::Operand( const csmp::Index& idx, const ScalarVariable& sc )
    : name_("tmp"), 
      prop_key_(idx),
      omin(-1.0e+50), omax(1.0e+50),
      flag_essential_(DIRICH),
      flag_output_(PLAIN),
      calc_offset_(0),
      scalar_storage_(sc)
 {
 } // end 


template<uint32_t dim>
Operand<dim>::Operand( const csmp::Index& idx, const VectorVariable<dim>& vc )
    : name_("tmp"), 
      prop_key_(idx),
      omin(-1.0e+50), omax(1.0e+50),
      flag_essential_(DIRICH),
      flag_output_(PLAIN),
      calc_offset_(0),
      vector_storage_(vc)
{
} // end


template<uint32_t dim>
Operand<dim>::Operand( const csmp::Index& idx, const TensorVariable<dim>& ts )
    : name_("tmp"), 
      prop_key_(idx),
      omin(-1.0e+50), omax(1.0e+50),
      flag_essential_(DIRICH),
      flag_output_(PLAIN),
      calc_offset_(0),
      tensor_storage_(ts)
 {
 } // end 
 

template<uint32_t dim>
Operand<dim>::Operand( const csmp::Index& idx, const ArrayVariable& ar )
    : name_("tmp"),
      prop_key_(idx),
      omin(-1.0e+50), omax(1.0e+50),
      flag_essential_(DIRICH),
      flag_output_(PLAIN),
      calc_offset_(0),
      array_storage_(ar)
 {
 } // end

template<uint32_t dim>
Operand<dim>::Operand( const csmp::Index& idx, const FlaggedArrayVariable& ar )
    : name_("tmp"),
      prop_key_(idx),
      omin(-1.0e+50), omax(1.0e+50),
      flag_essential_(DIRICH),
      flag_output_(PLAIN),
      calc_offset_(0),
      flagged_array_storage_(ar)
 {
 } // end

template<uint32_t dim>
 Operand<dim>::~Operand()
 {
 }


 /**
 Inputs the CSMP basic variable associated with the Operand into the
 supplied variable.

 @param s A CSMP variable, ie., a ScalarVariable, VectorVariable or TensorVariable.

 @section application Application .

 When you cannot use any of the overloaded operators of the Operand class
 for your specific calculation, just output the Operand variable into a
 basic CSMP variable and manipulate it directly.

 @section messages Messages

 AssignTo() will cause the program to terminate if you try to assign an
 Operand variable to a variable of the wrong type.
 
 */
template<uint32_t dim>
void Operand<dim>::AssignTo( ScalarVariable& s ) const
  {
     if ( prop_key_.type == SCALAR ) {
          s = scalar_storage_;
          return;
       }

     if ( prop_key_.type == VECTOR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of VECTOR type.\n";
     else if ( prop_key_.type == TENSOR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of TENSOR type.\n";
     else if ( prop_key_.type == ARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of ARRAY type.\n";
     else if ( prop_key_.type == FLAGGEDARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of FLAGGEDARRAY type.\n";
  }


template<uint32_t dim>
void Operand<dim>::AssignTo( VectorVariable<dim>& v ) const
  {
     if ( prop_key_.type == VECTOR ) {
          v = vector_storage_;
          return;
       }

     if ( prop_key_.type == SCALAR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of SCALAR type.\n";
     else if ( prop_key_.type == TENSOR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of TENSOR type.\n";
     else if ( prop_key_.type == ARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of ARRAY type.\n";
     else if ( prop_key_.type == FLAGGEDARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of FLAGGEDARRAY type.\n";
  }


template<uint32_t dim>
void Operand<dim>::AssignTo( TensorVariable<dim>& t ) const
  {
     if ( prop_key_.type == TENSOR ) {
          t = tensor_storage_;
          return;
       }

     if ( prop_key_.type == SCALAR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of SCALAR type.\n";
     else if ( prop_key_.type == VECTOR )
      std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of VECTOR type.\n";
     else if ( prop_key_.type == ARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of ARRAY type.\n";
     else if ( prop_key_.type == FLAGGEDARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of FLAGGEDARRAY type.\n";
  }


template<uint32_t dim>
void Operand<dim>::AssignTo( ArrayVariable& a ) const
  {
     if ( prop_key_.type == ARRAY ) {
          a = array_storage_;
          return;
       }

     if ( prop_key_.type == SCALAR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of SCALAR type.\n";
     else if ( prop_key_.type == VECTOR )
      std::cerr <<"\nOperand<dim>::AssignTo:  Expects argument of VECTOR type.\n";
     else if ( prop_key_.type == TENSOR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of TENSOR type.\n";
     else if ( prop_key_.type == FLAGGEDARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of FLAGGEDARRAY type.\n";
  }


template<uint32_t dim>
void Operand<dim>::AssignTo( FlaggedArrayVariable& a ) const
  {
     if ( prop_key_.type == FLAGGEDARRAY ) {
          a = flagged_array_storage_;
          return;
       }

     if ( prop_key_.type == SCALAR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of SCALAR type.\n";
     else if ( prop_key_.type == VECTOR )
      std::cerr <<"\nOperand<dim>::AssignTo:  Expects argument of VECTOR type.\n";
     else if ( prop_key_.type == TENSOR )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of TENSOR type.\n";
     else if ( prop_key_.type == ARRAY )
       std::cerr <<"\nOperand<dim>::AssignTo: Expects argument of ARRAY type.\n";
  }


/**
 
This assignment operator creates an Operand that is exactly identical to
the Operand on the left side of the expression.  Thus, it will overwrite 
whatever definition the Operand on the left side of the expression had and
using this operator for 
 
K = Pf
 
While this is exactly the functionality which one expects from a C++
assignment operator, it means that K will refer to Pf from now on and 
the variable with which K was associated, will not be accessible anymore 
by this Operand. 

@param op A reference to the Operand which will be assigned.

@section implementation Implementation

The method checks for self-assigment before setting the internal data. 

@section application Application

The copy constructor and for duplication of Operands. 
*/
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator=( const Operand<dim>& op ) 
 {
    // check for self assigment
    if ( &op != this ) 
      {
         name_              = op.name_;
         prop_key_          = op.prop_key_;
         omin               = op.omin;
         omax               = op.omax;
         flag_output_       = op.flag_output_;
         flag_essential_    = op.flag_essential_;
         calc_offset_       = op.calc_offset_;
         if      ( op.prop_key_.type == SCALAR )        scalar_storage_         = op.scalar_storage_;
         else if ( op.prop_key_.type == VECTOR )        vector_storage_         = op.vector_storage_;
         else if ( op.prop_key_.type == TENSOR )        tensor_storage_         = op.tensor_storage_;
         else if ( op.prop_key_.type == ARRAY  )        array_storage_          = op.array_storage_;
         else if ( op.prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  = op.flagged_array_storage_;
      }
    return *this;
 }



template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator=( const ScalarVariable& s )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ = s;
    else if ( prop_key_.type == VECTOR ) vector_storage_ = s;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ = s;
    else if ( prop_key_.type == ARRAY  ) array_storage_  = s;
    else if ( prop_key_.type == FLAGGEDARRAY  ) flagged_array_storage_  = s;


    return *this;
 }


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator=( const VectorVariable<dim>& vc )
 {
    if ( prop_key_.type != VECTOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator=(vector)", name_.c_str(),
                            "can only assign vectors to VectorVariable Operands");
    vector_storage_ = vc;

    return *this;
 }


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator=( const TensorVariable<dim>& ts )
 {
    if ( prop_key_.type != TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator=(tensor)", name_.c_str(),
                                   "can only assign tensors to TensorVariable Operands");
    tensor_storage_ = ts;

    return *this;
 }

template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator=( const ArrayVariable& ar )
 {
    if ( prop_key_.type != ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator=(array)", name_.c_str(),
                                   "can only assign arrays to ArrayVariable Operands");
    array_storage_ = ar;

    return *this;
 }

template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator=( const FlaggedArrayVariable& ar )
 {
    if ( prop_key_.type != FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator=(flagged array)", name_.c_str(),
                                   "can only assign flagged arrays to FlaggedArrayVariable Operands");
    flagged_array_storage_ = ar;

    return *this;
 }





template<uint32_t dim>
std::string   Operand<dim>::Name() const { return name_; }

template<uint32_t dim>
void          Operand<dim>::Name( const char* s ) { name_ = s; }

template<uint32_t dim>
const csmp::Index& Operand<dim>::Key() const { return prop_key_; }

template<uint32_t dim>
size_t     Operand<dim>::Index() const { return prop_key_.index; }

template<uint32_t dim>
PLACEMENT     Operand<dim>::Placement() const { return prop_key_.place; }

template<uint32_t dim>
VARIABLE_TYPE Operand<dim>::Type() const { return prop_key_.type; }

template<uint32_t dim>
size_t     Operand<dim>::CalculationOffset() const { return calc_offset_; }

template<uint32_t dim>
void          Operand<dim>::CalculationOffset( size_t o ) { calc_offset_=o; }

template<uint32_t dim>
void          Operand<dim>::OutputCondition( VARIABLE_FLAG c ) { flag_output_ = c; }

template<uint32_t dim>
VARIABLE_FLAG    Operand<dim>::OutputCondition() const { return flag_output_; }

template<uint32_t dim>
VARIABLE_FLAG    Operand<dim>::EssentialCondition() const { return flag_essential_; }

template<uint32_t dim>
void          Operand<dim>::EssentialCondition( VARIABLE_FLAG c ) { flag_essential_ = c; }




/**

Assignments of CSP basic variables to Operands are implemented as
assignment operators in which the type of the variable which is associated
with the Operand determines the outcome of the assignment. This restricts
the number of meaningful assignments for CSP basic variables: Scalars can
be assigned to vectors and tensors but not the other way around. Vectors
cannot be assigned to tensors and vice versa.

@param  val Constant references to double, ScalarVariable, VectorVariable, and
TensorVariable.

@section implementation Implementation

The assignments are implemented using the overloaded operators of the
CSP basic variables.

@section messages Messages

Erratic assignments may halt the program and lead to an error report
including the name of the Operand to which the assignment was
attempted.
*/
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator=( double val )
{
    switch ( prop_key_.type )
    {
         case SCALAR:
              scalar_storage_ = val;
           break;
         case VECTOR:
              vector_storage_ = val;
           break;
         case TENSOR:
              tensor_storage_ = val;
              break;
         case ARRAY:
              array_storage_ = val;
              break;
         case FLAGGEDARRAY:
              flagged_array_storage_ = val;
              break;
        default:
           throw csmp::Exception( ERROR, "Operand<dim>::operator=(double)",
                                  "Array type variables are not handled yet.");
    }
    return *this;
}



/**
 
Create temporary Operands with the same csmp::Index as the Operand on the
left of the operator in the mathematical expression. The basic CSP
variable of the temporary Operands will have been transformed using the
operator an the right-side argument of it. 

@param val A double value typically typedef'ed as a fT. The value will appear
on the right of the operator in the mathematical expression involving
the Operand. 

@section application Application

NOTE: These overloaded Operands will return Operands which have exactly
the same csmp::Index as the Operand on the lefthandside of the operator
sign in the expression in which they are used. Thus beware of writing
something like:

@code
Op1 = Op2 * 322.4;
@endcode

This will overwrite the Op1 with the csmp::Index of Op2 because the 
assignment operator is used. The overloaded operators are still useful 
for expressions like:

@code
Op1 += Op2 * 322.4;
@endcode
 */
template<uint32_t dim>
Operand<dim>  Operand<dim>::operator+( double val ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ + val) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ + val) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ + val) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  + val) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  + val) );

    return Operand<dim>();
 }
 

template<uint32_t dim>
Operand<dim>  Operand<dim>::operator-( double val ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ - val) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ - val) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ - val) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  - val) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  - val) );

    return Operand<dim>();
 }
 
 
template<uint32_t dim>
Operand<dim>  Operand<dim>::operator*( double val ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ * val) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ * val) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ * val) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  * val) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  * val) );

    return Operand<dim>();
 }
 
 
template<uint32_t dim>
Operand<dim>  Operand<dim>::operator/( double val ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ / val) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ / val) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ / val) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  / val) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  / val) );

    return Operand<dim>();
 }




/**
 
Create temporary Operands with the same csmp::Index as the Operand on the
left of the operator in the mathematical expression. The basic CSP
variable of the temporary Operands will have been transformed using the
operator an the right-side argument of it. 

@param op An Operand which will appear on the right of the operator in the mathematical
expression involving the lefthand Operand. 

@section application Application

NOTE: These overloaded Operands will return Operands which have exactly
the same csmp::Index as the Operand on the lefthandside of the operator
sign in the expression in which they are used. Thus beware of writing
something like:
  
@code
Op1 = Op2 * 322.4;
@endcode
  
This will overwrite the Op1 with the csmp::Index of Op2 because the 
assignment operator is used. The overloaded operators are still useful 
for expressions like:
  
@code
Op1 += Op2 * 322.4;
@endcode

@section implementation Implementation

These overloaded operators use the overloaded operators of CSP basic
variables. 

@section messages Messages

Depending on the basic variable types of the left and rigth Operands 
errors may be reported, for instance, when you try to add a tensor to
a scalar variable. Mathematically meaningful operations, like multiplying
a right-side vector with a left-side tensor, are enabled. 
 */
template<uint32_t dim>
Operand<dim>  Operand<dim>::operator+( const Operand<dim>& op ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ + op.scalar_storage_) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ + op.vector_storage_) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ + op.tensor_storage_) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  + op.array_storage_) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  + op.flagged_array_storage_) );

    return Operand<dim>();
 }
 
 
template<uint32_t dim>
Operand<dim>  Operand<dim>::operator-( const Operand<dim>& op ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ - op.scalar_storage_) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ - op.vector_storage_) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ - op.tensor_storage_) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  - op.array_storage_) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  - op.flagged_array_storage_) );

    return Operand<dim>();
 }
 
 
template<uint32_t dim>
Operand<dim>  Operand<dim>::operator*( const Operand<dim>& op ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ * op.scalar_storage_) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ * op.vector_storage_) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ * op.tensor_storage_) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  * op.array_storage_) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  * op.flagged_array_storage_) );

    return Operand<dim>();
 }
 
 
template<uint32_t dim>
Operand<dim>  Operand<dim>::operator/( const Operand<dim>& op ) const
 {
    if ( prop_key_.type == SCALAR )
      return Operand<dim>( prop_key_, (scalar_storage_ / op.scalar_storage_) );
    if ( prop_key_.type == VECTOR )
      return Operand<dim>( prop_key_, (vector_storage_ / op.vector_storage_) );
    if ( prop_key_.type == TENSOR )
      return Operand<dim>( prop_key_, (tensor_storage_ / op.tensor_storage_) );
    if ( prop_key_.type == ARRAY )
      return Operand<dim>( prop_key_, (array_storage_  / op.array_storage_) );
    if ( prop_key_.type == FLAGGEDARRAY )
      return Operand<dim>( prop_key_, (flagged_array_storage_  / op.flagged_array_storage_) );

    return Operand<dim>();
 }




/**
 
Transforms the basic variable associated with the lefthand Operand by the 
right-side double using the arithmetic operators +, -, *, or / . The C++
short form for writing a = a + b; as a += b; is used in such expressions. 

@param val A double value typically typedef'ed as a fT. The value will appear
on the right of the operator in the mathematical expression involving
the Operand. 

@section application Application

Avoids the creation of temporaries and should therefore be prefered in
mathematical expressions involving Operands.  

@section implementation Implementation

The overloaded operators use the overloaded operators of CSP basic
variables. 

@section messages Messages

Depending on the basic variable types of the left and rigth Operands 
errors may be reported, for instance, when you try to add a tensor to
a scalar variable. Mathematically meaningful operations, like multiplying
a right-side vector with a left-side tensor, are enabled. 
 */
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator+=( double val )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ += val;
    else if ( prop_key_.type == VECTOR ) vector_storage_ += val;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ += val;
    else if ( prop_key_.type == ARRAY )  array_storage_  += val;
    else if ( prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  += val;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator-=( double val )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ -= val;
    else if ( prop_key_.type == VECTOR ) vector_storage_ -= val;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ -= val;
    else if ( prop_key_.type == ARRAY  ) array_storage_  -= val;
    else if ( prop_key_.type == FLAGGEDARRAY  ) flagged_array_storage_  -= val;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator*=( double val )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ *= val;
    else if ( prop_key_.type == VECTOR ) vector_storage_ *= val;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ *= val;
    else if ( prop_key_.type == ARRAY )  array_storage_  *= val;
    else if ( prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  *= val;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator/=( double val )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ /= val;
    else if ( prop_key_.type == VECTOR ) vector_storage_ /= val;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ /= val;
    else if ( prop_key_.type == ARRAY )  array_storage_  /= val;
    else if ( prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  /= val;

    return *this;
 } // end



/**
 
Transforms the basic variable associated with the lefthand Operand by the 
right-side  ScalarVariable using the arithmetic operators +, -, *, or / . 
The C++ short form for writing a = a + b; as a += b; is used in such 
expressions. 

@param sc A ScalarVariable value typically typedef'ed as a fT. The value will
appear on the right of the operator in the mathematical expression involving
the Operand. 

@section application Application

Avoids the creation of temporaries and should therefore be prefered in
mathematical expressions involving Operands.  

@section implementation Implementation

The overloaded operators use the overloaded operators of CSP basic
variables. 

@section messages Messages

Depending on the basic variable types of the left and rigth Operands 
errors may be reported, for instance, when you try to add a tensor to
a scalar variable. Mathematically meaningful operations, like multiplying
a right-side vector with a left-side tensor, are enabled. 

*/
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator+=( const ScalarVariable& sc )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ += sc;
    else if ( prop_key_.type == VECTOR ) vector_storage_ += sc;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ += sc;
    else if ( prop_key_.type == ARRAY )  array_storage_  += sc;
    else if ( prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  += sc;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator-=( const ScalarVariable& sc )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ -= sc;
    else if ( prop_key_.type == VECTOR ) vector_storage_ -= sc;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ -= sc;
    else if ( prop_key_.type == ARRAY )  array_storage_  -= sc;
    else if ( prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  -= sc;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator*=( const ScalarVariable& sc )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ *= sc;
    else if ( prop_key_.type == VECTOR ) vector_storage_ *= sc;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ *= sc;
    else if ( prop_key_.type == ARRAY )  array_storage_  *= sc;
    else if ( prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  *= sc;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator/=( const ScalarVariable& sc )
 {
    if      ( prop_key_.type == SCALAR ) scalar_storage_ /= sc;
    else if ( prop_key_.type == VECTOR ) vector_storage_ /= sc;
    else if ( prop_key_.type == TENSOR ) tensor_storage_ /= sc;
    else if ( prop_key_.type == ARRAY )  array_storage_  /= sc;
    else if ( prop_key_.type == FLAGGEDARRAY )  flagged_array_storage_  /= sc;

    return *this;
 } // end



/**
 
Transforms the basic variable associated with the lefthand Operand by the 
right-side VectorVariable using the arithmetic operators +, -, *, or / . 
The C++ short form for writing a = a + b; as a += b; is used in such 
expressions. 

@param vc A VectorVariable which will appear on the right of the operator in
the mathematical expression involving the Operand. 

@section application Application

Avoids the creation of temporaries and should therefore be prefered in
mathematical expressions involving Operands.  

@section implementation Implementation

The overloaded operators use the overloaded operators of CSP basic
variables. 

@section messages Messages

Depending on the basic variable types of the left and rigth Operands 
errors may be reported, for instance, when you try to add a tensor to
a scalar variable. Mathematically meaningful operations, like multiplying
a right-side vector with a left-side tensor, are enabled. 

For the addition (+) and subtraction (-) operators the application is 
restricted to the case where the variable associated with the Operand
is also a vector. If this is not the case the program will be 
terminated and a fatal error will be reported. 

 */
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator+=( const VectorVariable<dim>& vc )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding vector to scalar");
    else if ( prop_key_.type == VECTOR ) vector_storage_ += vc;
    else if ( prop_key_.type == TENSOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(), 
                                   "no rule for adding vector to tensor");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding vector to array");
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding vector to flagged array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator-=( const VectorVariable<dim>& vc )
 {
    if      ( prop_key_.type == SCALAR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting vector from scalar");
    else if ( prop_key_.type == VECTOR ) vector_storage_ -= vc;
    else if ( prop_key_.type == TENSOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting vector from tensor");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting vector from array");
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting vector from flagged array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator*=( const VectorVariable<dim>& vc )
 {
    if      ( prop_key_.type == SCALAR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying scalar with vector");
    else if ( prop_key_.type == VECTOR ) vector_storage_ *= vc;
    else if ( prop_key_.type == TENSOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "result of multiplying tensor with vector would be a vector");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying vector with array");
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying vector with flagged array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator/=( const VectorVariable<dim>& vc )
 {
    if      ( prop_key_.type == SCALAR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing scalar by vector");
    else if ( prop_key_.type == VECTOR ) vector_storage_ /= vc;
    else if ( prop_key_.type == TENSOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for for dividing tensor by vector");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing vector by array");
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing vector by flagged array");
    return *this;
 } // end

/// Array
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator+=( const ArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding array to scalar");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                     "no rule for adding array to vector");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding array to tensor");
    else if ( prop_key_.type == ARRAY ) array_storage_ += ar;
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding array to flagged array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator-=( const ArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting array from scalar");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                     "no rule for subtractiong array from vector");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting array from tensor");
    else if ( prop_key_.type == ARRAY ) array_storage_ -= ar;
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtractiong array from flagged array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator*=( const ArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying scalar with array");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                     "no rule for multiplying vector with array");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying scalar with array");
    else if ( prop_key_.type == ARRAY ) array_storage_ *= ar;
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying flagged array with array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator/=( const ArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing scalar by array");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                     "no rule for dividing vector by array");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing tensor by array");
    else if ( prop_key_.type == ARRAY ) array_storage_ /= ar;
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing flagged array by array");

    return *this;
 } // end

/// Flagged Array
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator+=( const FlaggedArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding flagged array to scalar");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                     "no rule for adding flagged array to vector");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding flagged array to tensor");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding flagged array to array");
    else if ( prop_key_.type == FLAGGEDARRAY ) flagged_array_storage_ += ar;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator-=( const FlaggedArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting flagged array from scalar");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                     "no rule for subtractiong flagged array from vector");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting flagged array from tensor");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtractiong flagged array from array");
    else if ( prop_key_.type == FLAGGEDARRAY ) flagged_array_storage_ -= ar;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator*=( const FlaggedArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying scalar with flagged array");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                     "no rule for multiplying vector with flagged array");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying scalar with flagged array");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying array with flagged array");
    else if ( prop_key_.type == FLAGGEDARRAY ) flagged_array_storage_ *= ar;

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator/=( const FlaggedArrayVariable& ar )
 {
    if      ( prop_key_.type == SCALAR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing scalar by flagged array");
    else if ( prop_key_.type == VECTOR )
        throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                     "no rule for dividing vector by flagged array");
    else if ( prop_key_.type == TENSOR )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing tensor by flagged array");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing array by flagged array");
    else if ( prop_key_.type == FLAGGEDARRAY ) flagged_array_storage_ /= ar;

    return *this;
 } // end

/**
 
Transforms the basic variable associated with the lefthand Operand by the
right-side TensorVariable using the arithmetic operators +, -, *, or / . 
The C++ short form for writing a = a + b; as a += b; is used in such 
expressions. 

@param ts A TensorVariable which will appear on the right of the operator in the
mathematical expression involving the Operand. 

@section application Application

Avoids the creation of temporaries and should therefore be prefered in
mathematical expressions involving Operands.  

@section implementation Implementation

The overloaded operators use the overloaded operators of CSP basic
variables. 

@section messages Messages

Depending on the basic variable types of the left and rigth Operands 
errors may be reported, for instance, when you try to add a tensor to
a scalar variable. Mathematically meaningful operations, like multiplying
a right-side vector with a left-side tensor, are enabled. 

Right-side TensorVariables can only be added to, subtracted from, multiplied 
with or divide other TensorVariables. The left-side Operand must therefore
also be a TensorVariable. If this is not the case the program will be 
terminated and a fatal error will be reported. 

*/
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator+=( const TensorVariable<dim>& ts )
 {
    if      ( prop_key_.type == TENSOR ) tensor_storage_ += ts;
    else if ( prop_key_.type == SCALAR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding tensor to scalar");
    else if ( prop_key_.type == VECTOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(), 
                                   "no rule for adding tensor to vector");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding tensor to array");
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for adding tensor to flagged array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator-=( const TensorVariable<dim>& ts )
 {
    if      ( prop_key_.type == TENSOR ) tensor_storage_ -= ts;
    else if ( prop_key_.type == SCALAR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting tensor from scalar");
    else if ( prop_key_.type == VECTOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting tensor from vector");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting tensor from array");
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator-=", name_.c_str(),
                                   "no rule for subtracting tensor from flagged array");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator*=( const TensorVariable<dim>& ts )
 {
    if      ( prop_key_.type == TENSOR ) tensor_storage_ *= ts;
    else if ( prop_key_.type == SCALAR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying scalar with tensor");
    else if ( prop_key_.type == VECTOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multiplying vector with tensor");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multying array with tensor");
    else if ( prop_key_.type == FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator*=", name_.c_str(),
                                   "no rule for multying flagged array with tensor");

    return *this;
 } // end


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator/=( const TensorVariable<dim>& ts )
 {
    if      ( prop_key_.type == TENSOR ) tensor_storage_ /= ts;
    else if ( prop_key_.type == SCALAR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing scalar by tensor");
    else if ( prop_key_.type == VECTOR ) 
      throw csmp::Exception( FATAL_ERROR, "Operand::operator/=", name_.c_str(),
                                   "no rule for dividing vector by tensor");
    else if ( prop_key_.type == ARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for dividing array by tensor");
    else if ( prop_key_.type ==FLAGGEDARRAY )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator+=", name_.c_str(),
                                   "no rule for dividing flagged array by tensor");

    return *this;
 } // end



/**
 
Transforms the basic variable associated with the lefthand Operand by the 
right-side Operand variable using the arithmetic operators +, -, *, or / . 
The C++ short form for writing a = a + b; as a += b; is used in such 
expressions. Since only the mathematical transformation is defined,
the other internal variables of the lefthand Operand will not be
affected by a call to one of these operators. 

@param op A constant reference to an Operand which will appear on the right of the
operator in the mathematical expression involving the Operand. 

@section application Application

Avoids the creation of temporaries and should therefore be prefered in
mathematical expressions involving Operands. This is also the safe way
of writing expressions involving several Operands, because the 
possibility of accidentially overwriting an Operand with another 
is ruled out.  

@section implementation Implementation

The overloaded operators use the overloaded operators of CSP basic
variables. 

@section messages Messages

Depending on the basic variable types of the left and rigth Operands 
errors may be reported, for instance, when you try to add a tensor to
a scalar variable. Mathematically meaningful operations, like multiplying
a right-side vector with a left-side tensor, are enabled. 
*/
template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator+=( const Operand<dim>& op )
 {
    if      ( prop_key_.type == SCALAR )
      {
         if ( op.prop_key_.type == SCALAR ) return *this += op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this += op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this += op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this += op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this += op.flagged_array_storage_;
      }
    else if ( prop_key_.type == VECTOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this += op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this += op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this += op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this += op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this += op.flagged_array_storage_;
      }
    else if ( prop_key_.type == TENSOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this += op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this += op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this += op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this += op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this += op.flagged_array_storage_;
      }
    else if ( prop_key_.type == ARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this += op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this += op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this += op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this += op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this += op.flagged_array_storage_;
      }
    else if ( prop_key_.type == FLAGGEDARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this += op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this += op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this += op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this += op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this += op.flagged_array_storage_;
      }


    return *this;
 }


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator-=( const Operand<dim>& op )
 {
    if      ( prop_key_.type == SCALAR )
      {
         if ( op.prop_key_.type == SCALAR ) return *this -= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this -= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this -= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this -= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this -= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == VECTOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this -= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this -= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this -= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this -= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this -= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == TENSOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this -= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this -= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this -= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this -= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this -= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == ARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this -= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this -= op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this -= op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this -= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this -= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == FLAGGEDARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this -= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this -= op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this -= op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this -= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this -= op.flagged_array_storage_;
      }

    return *this;
 }


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator*=( const Operand<dim>& op )
 {
    if      ( prop_key_.type == SCALAR )
      {
         if ( op.prop_key_.type == SCALAR ) return *this *= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this *= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this *= op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this *= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this *= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == VECTOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this *= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this *= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this *= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this *= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this *= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == TENSOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this *= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this *= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this *= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this *= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this *= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == ARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this *= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this *= op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this *= op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this *= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this *= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == FLAGGEDARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this *= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this *= op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this *= op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this *= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this *= op.flagged_array_storage_;
      }
    return *this;
 }


template<uint32_t dim>
Operand<dim>&  Operand<dim>::operator/=( const Operand<dim>& op )
 {
    if      ( prop_key_.type == SCALAR )
      {
         if ( op.prop_key_.type == SCALAR ) return *this /= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this /= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this /= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this /= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this /= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == VECTOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this /= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this /= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this /= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this /= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this /= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == TENSOR ) 
      {
         if ( op.prop_key_.type == SCALAR ) return *this /= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this /= op.vector_storage_; 
         if ( op.prop_key_.type == TENSOR ) return *this /= op.tensor_storage_; 
         if ( op.prop_key_.type == ARRAY )  return *this /= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this /= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == ARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this /= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this /= op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this /= op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this /= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this /= op.flagged_array_storage_;
      }
    else if ( prop_key_.type == FLAGGEDARRAY )
      {
         if ( op.prop_key_.type == SCALAR ) return *this /= op.scalar_storage_;
         if ( op.prop_key_.type == VECTOR ) return *this /= op.vector_storage_;
         if ( op.prop_key_.type == TENSOR ) return *this /= op.tensor_storage_;
         if ( op.prop_key_.type == ARRAY )  return *this /= op.array_storage_;
         if ( op.prop_key_.type == FLAGGEDARRAY )  return *this /= op.flagged_array_storage_;
      }

    return *this;
 }
 
 

/**
 
Operands can be compared with other Operands, double values, and CSP
basic variables. Whether this comparison is meaningful, depends on the 
variable types on the left and right-side of the comparitors,
and errors will be reported for obviously nonsensical comparisons.
Comparisons of variables of the same type are straightforward and the 
comparison of a ScalarVariable with a double is as well. 

If a vector variable is compared with a scalar, this implementation of 
comparitors will compare the length of the vector with the scalar,
since in most cases this is a physically sensible thing to do.  

If you compare a tensor variable with a scalar, this comparison will
use the largest or the smallest element in the tensor, but the 
comparitors '==' and '!=' are not defined. Vectors can  
not be compared with tensors.  

Only Operands with the same type of associated variables can be compared. 

@param val A constant reference to any of the types, double, Operand, ScalarVariable,
VectorVariable, or TensorVariable.  

@return A boolean variable which evaluates to true or false.

@section implementation Implementation

In as far as possible, the overloaded operators of the basic CSP 
variables are used to ascertain consistency of CSP variable behaviour. 
No provisions are made for the comparison of enumeration
VARIABLE_FLAG flags. This would open a confusing array of possibilities.
All comparitors are constant methods which means that they do not 
change the state of Operands. 

@section application Application

Typically, you want to compare Operands with floating-point numbers in 'if' 
statements inside of Calculate() methods of Interrelation subclasses. 
This may involve statements like if ( T > 500.0 )... 
This type of comparisons is fully supported by the 
overloaded comparitors of the Operand class. 

@section messages Messages

The most common runtime message that you may get using a comparitor is
 
'cannot compare * variable with * variable'
 
To make it easier to track down where this error occurred, the name of 
the Operand involved in the comparison will be reported as well. 
Remember, however, that temporary Operands have no distinct names. 
*/
template<uint32_t dim>
bool  Operand<dim>::operator>( double val ) const
 {
    if ( prop_key_.type == SCALAR ) return (scalar_storage_() > val);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator>", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with float variable");
    return false;
 }
 
 
template<uint32_t dim>
bool  Operand<dim>::operator<( double val ) const
 {
    if ( prop_key_.type == SCALAR ) return (scalar_storage_() < val);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator<", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with float variable");
    return false;
 }
 
 
template<uint32_t dim>
bool  Operand<dim>::operator>=( double val ) const
 {
    if ( prop_key_.type == SCALAR ) return (scalar_storage_() >= val);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator>=", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with float variable");
    return false;
 }
 
 
template<uint32_t dim>
bool  Operand<dim>::operator<=( double val ) const
 {
    if ( prop_key_.type == SCALAR ) return (scalar_storage_() <= val);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator<=", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with float variable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator==( double val ) const
 {
    if ( prop_key_.type == SCALAR ) return approximatelyEqual( scalar_storage_(), val );

    throw csmp::Exception( FATAL_ERROR, "Operand::operator==", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with float variable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator!=( double val ) const
 {
    if ( prop_key_.type == SCALAR ) return !approximatelyEqual( scalar_storage_(), val );


    throw csmp::Exception( FATAL_ERROR, "Operand::operator!=", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with float variable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator==( const ScalarVariable& sc ) const 
 {
    if ( prop_key_.type == SCALAR ) return (scalar_storage_ == sc);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator==", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with ScalarVariable");
    return false;
 }

template<uint32_t dim>
bool  Operand<dim>::operator!=( const ScalarVariable& sc ) const 
 {
    if ( prop_key_.type == SCALAR ) return (scalar_storage_ != sc);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator!=", name_.c_str(),
                           "cannot compare vector, tensor, array or flagged array variable with ScalarVariable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator==( const VectorVariable<dim>& vc ) const
 {
    if ( prop_key_.type == VECTOR ) return (vector_storage_ == vc);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator==", name_.c_str(),
                           "cannot compare scalar, tensor, array or flagged array variable with VectorVariable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator!=( const VectorVariable<dim>& vc ) const
 {
    if ( prop_key_.type == VECTOR ) return (vector_storage_ != vc);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator!=", name_.c_str(),
                           "cannot compare scalar, tensor, array or flagged array variable with VectorVariable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator==( const TensorVariable<dim>& ts ) const
 {
    if ( prop_key_.type == TENSOR ) return (tensor_storage_ == ts);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator==", name_.c_str(),
                           "cannot compare scalar, vector, array or flagged array variable with TensorVariable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator!=( const TensorVariable<dim>& ts ) const 
 {
    if ( prop_key_.type == TENSOR ) return (tensor_storage_ != ts);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator!=", name_.c_str(),
                           "cannot compare scalar, vector, array or flagged array variable with TensorVariable");
    return false;
 }

template<uint32_t dim>
bool  Operand<dim>::operator==( const ArrayVariable& ar ) const
 {
    if ( prop_key_.type == ARRAY ) return (array_storage_ == ar);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator==", name_.c_str(),
                           "cannot compare scalar, vector, tensor or flagged array variable with ArrayVariable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator!=( const ArrayVariable& ar ) const
 {
    if ( prop_key_.type == ARRAY ) return (array_storage_ != ar);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator!=", name_.c_str(),
                           "cannot compare scalar, vector, tensor or flagged array variable with ArrayVariable");
    return false;
 }

template<uint32_t dim>
bool  Operand<dim>::operator==( const FlaggedArrayVariable& ar ) const
 {
    if ( prop_key_.type == FLAGGEDARRAY ) return (flagged_array_storage_ == ar);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator==", name_.c_str(),
                           "cannot compare scalar, vector, tensor or array variable with FlaggedArrayVariable");
    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator!=( const FlaggedArrayVariable& ar ) const
 {
    if ( prop_key_.type == FLAGGEDARRAY ) return (flagged_array_storage_ != ar);

    throw csmp::Exception( FATAL_ERROR, "Operand::operator!=", name_.c_str(),
                           "cannot compare scalar, vector, tensor or array variable with FlaggedArrayVariable");
    return false;
 }

template<uint32_t dim>
bool  Operand<dim>::operator>( Operand<dim>& op ) const
 {
    if ( prop_key_.type != op.prop_key_.type )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator>", name_.c_str(),
                                   "cannot compare Operands with different variable types");
     
    if ( prop_key_.type == SCALAR ) 
      return (scalar_storage_ > op.scalar_storage_);
      
    if ( prop_key_.type == VECTOR )
      return (vector_storage_.Length() > op.vector_storage_.Length());
      
    if ( prop_key_.type == TENSOR )
      return (tensor_storage_.MaxElement() > op.tensor_storage_.MaxElement());   

    if ( prop_key_.type == ARRAY )
        return (array_storage_.Size() > op.array_storage_.Size());

    if ( prop_key_.type == FLAGGEDARRAY )
        return (flagged_array_storage_.Size() > op.flagged_array_storage_.Size());

    return false;
 }
 
 
template<uint32_t dim>
bool  Operand<dim>::operator<( Operand<dim>& op ) const
 {
    if ( prop_key_.type != op.prop_key_.type )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator<", name_.c_str(),
                                   "cannot compare Operands with different variable types");
     
    if ( prop_key_.type == SCALAR ) 
      return (scalar_storage_ < op.scalar_storage_);
      
    if ( prop_key_.type == VECTOR )
      return (vector_storage_.Length() < op.vector_storage_.Length());
      
    if ( prop_key_.type == TENSOR )
      return (tensor_storage_.MinElement() < op.tensor_storage_.MinElement());   

    if ( prop_key_.type == ARRAY )
        return (array_storage_.Size() < op.array_storage_.Size());

    if ( prop_key_.type == FLAGGEDARRAY )
        return (flagged_array_storage_.Size() < op.flagged_array_storage_.Size());

    return false;
 }
 
 
template<uint32_t dim>
bool  Operand<dim>::operator>=( Operand<dim>& op ) const
 {
    if ( prop_key_.type != op.prop_key_.type )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator>=", name_.c_str(),
                                   "cannot compare Operands with different variable types");
     
    if ( prop_key_.type == SCALAR ) 
      return (scalar_storage_ >= op.scalar_storage_);
      
    if ( prop_key_.type == VECTOR )
      return (vector_storage_.Length() >= op.vector_storage_.Length());
      
    if ( prop_key_.type == TENSOR )
      return (tensor_storage_.MaxElement() >= op.tensor_storage_.MaxElement());   

    if ( prop_key_.type == ARRAY )
        return (array_storage_.Size() >= op.array_storage_.Size());

    if ( prop_key_.type == FLAGGEDARRAY )
        return (flagged_array_storage_.Size() >= op.flagged_array_storage_.Size());

    return false;
 }

 
template<uint32_t dim>
bool  Operand<dim>::operator<=( Operand<dim>& op ) const
 {
    if ( prop_key_.type != op.prop_key_.type )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator<=", name_.c_str(),
                                   "cannot compare Operands with different variable types");
     
    if ( prop_key_.type == SCALAR ) 
      return (scalar_storage_ <= op.scalar_storage_);
      
    if ( prop_key_.type == VECTOR )
      return (vector_storage_.Length() <= op.vector_storage_.Length());
      
    if ( prop_key_.type == TENSOR )
      return (tensor_storage_.MinElement() <= op.tensor_storage_.MinElement());

    if ( prop_key_.type == ARRAY )
        return (array_storage_.Size() <= op.array_storage_.Size());

    if ( prop_key_.type == FLAGGEDARRAY )
        return (flagged_array_storage_.Size() <= op.flagged_array_storage_.Size());

    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator==( Operand<dim>& op ) const 
 {
    if ( prop_key_.type != op.prop_key_.type )
      throw csmp::Exception( FATAL_ERROR, "Operand::operator==", name_.c_str(),
                                   "cannot compare Operands with different variable types");
    if ( prop_key_.type == SCALAR ) 
      return (scalar_storage_ == op.scalar_storage_);
      
    if ( prop_key_.type == VECTOR )
      return (vector_storage_ == op.vector_storage_);
      
    if ( prop_key_.type == TENSOR )
      return (tensor_storage_ == op.tensor_storage_);   

    if ( prop_key_.type == ARRAY )
      return (array_storage_ == op.array_storage_);

    if ( prop_key_.type == FLAGGEDARRAY )
      return (flagged_array_storage_ == op.flagged_array_storage_);

    return false;
 }


template<uint32_t dim>
bool  Operand<dim>::operator!=( Operand<dim>& op ) const 
 {
    return !(*this == op);
 }


/**

Tests whether the value of the physical variable is within the bounds that
were specified for that variable inside the PropertyDatabase record. If
IsWithinRange() is called for a temporary Operand, the positions of the
variable value relative to the default variable range of -1e+50 to 1e+50
will be tested for. If the variable type is a vector, its length will be
compared with the range, and if the variable type is a tensor its largest
and its smallest element will be range-checked.

@section arguments Input Arguments

Optionally, rather than testing the variable which is associated with
the Operand, external CSP basic variables can be tested for the range
specified for the Operand. Just supply a constant reference to the
test variable to IsWithinRange().

@section application Application .

A range check is performed automatically on output variables in subclasses
of Interrelations, before the calculated values are stored in the
MemoryManager. This is a fundamental CSP feature, prohibiting error
propagation in complex simulations.
*/
template<uint32_t dim>
bool Operand<dim>::IsWithinRange() const
 {
    if ( prop_key_.type == SCALAR )
      return (scalar_storage_() >= omin && scalar_storage_() <= omax);

    if ( prop_key_.type == VECTOR )
      return (vector_storage_.Length() >= omin && vector_storage_.Length() <= omax);

    if ( prop_key_.type == TENSOR )
      return (tensor_storage_.MinElement() >= omin && tensor_storage_.MaxElement() <= omax);

    if ( prop_key_.type == ARRAY )
        return array_storage_.IsWithinRange(omin,omax);

    if ( prop_key_.type == FLAGGEDARRAY )
      return flagged_array_storage_.IsWithinRange(omin,omax);

    return false;

 } // end IsWithinRange





/**

Tests the supplied value or components thereof for whether it lies within
the physical meaninful range which was assigned to the Operand when it
was constructed.

@section arguments Input Arguments

The value which shall be tested is supplied as first argument to the
overloaded method IsWithinRange(). The supported argument types are
double, and the basic CSP variables ScalarVariable, VectorVariable, and
TensorVariable.

@return The method returns the boolean variable 'true' if the argument value
is in the physical meaningful range; else it returns 'false'.

@section implementation Implementation

IsWithinRange() tests each component of its argument for whether it lies
in the specified range. Thus, for a tensor argument, 9 checks are
performed if the model is three dimensional.

@section application Application .

Before a property value is returned to the Model storage, it is
tested for whether it lies in the permitted range.

@section messages Messages

If the value is out of range the method prints a message to 'stderr'
and outputs the violating variable value(s).
*/
template<uint32_t dim>
template<typename cspT>
bool Operand<dim>::IsWithinRange( const cspT& val ) const
 {
    ErrorHandler& csmp_error( ErrorHandler::Instance() );

    if ( !val.IsWithinRange( omin, omax ) ) {
          std::cerr <<"\n'"<< name_ <<"' has erratic value:";
          val.Out();
          std::cerr << std::endl << std::endl;
          csmp_error.Note( ERROR, "Operand::IsWithinRange", "Property value out of range.");
          return false;
       }

     return true;

 } // end IsWithinRange


template<uint32_t dim>
double Operand<dim>::MinValue() const
{
    return omin;

} // end MinValue

template<uint32_t dim>
double Operand<dim>::MaxValue() const
{
    return omax;

} // end MaxValue












/**
 
Calls the Out() method of the current variable of the Operand. This will
print the CSP basic variable values and the VARIABLE_FLAG flags which is 
useful for testing a derived Interrelation before using it in
computations.  

*/
template<uint32_t dim>
void Operand<dim>::PrintValue() const
 {
    if      ( prop_key_.type == SCALAR )        scalar_storage_.Out();
    else if ( prop_key_.type == VECTOR )        vector_storage_.Out();
    else if ( prop_key_.type == TENSOR )        tensor_storage_.Out();
    else if ( prop_key_.type == ARRAY  )        array_storage_.Out();
    else if ( prop_key_.type == FLAGGEDARRAY  ) flagged_array_storage_.Out();
    
 } // end PrintValue

 

 
/** Standard CSP method to check the state of an Operand by printing it to stdout.

*/
template<uint32_t dim>
void  Operand<dim>::Out() const 
 {
    string stype    = parseType(prop_key_.type);
    string place    = parsePlacement(prop_key_.place);
    string esstatus = parseStatus(flag_essential_);
    string oustatus = parseStatus(flag_output_);
           
    cout <<"\n\nOperand<dim>: '" << name_ <<"'"<< endl;
    cout <<"\nproperty type:       "<< stype << endl;
    cout <<"placement:           "<< place << endl;
    cout <<"index:               "<< prop_key_.index << endl;
    cout <<"property range, min: "<< omin <<",  max: "<< omax << endl;
    cout <<"flag_essential:      "<< esstatus << endl;
    cout <<"flag_output:         "<< oustatus << endl;
    cout <<"calculation offset:  "<< calc_offset_ << endl;

    PrintValue();
     
 } // end Out



template class Operand<1U>;
template class Operand<2U>;
template class Operand<3U>;


} // end namespace csp


