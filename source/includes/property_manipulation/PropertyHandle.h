#ifndef CSMP_PROPERTY_HANDLE_H
#define CSMP_PROPERTY_HANDLE_H

#include "ScalarVariable.h"
#include "Element.h"

namespace csmp {

template<size_t> class VectorVariable;
template<size_t> class TensorVariable;
template<size_t> class Region;
template<size_t> class Model;

<<<<<<< HEAD
=======
template<size_t dim>
class PropertyHandle {
  public:
    PropertyHandle( Model<dim>& sg, const char* var_name,
                    VARIABLE_TYPE type=SCALAR, PLACEMENT place=NODE, size_t vsize=1U );

    PropertyHandle( Model<dim>& sg, const char* target_group, 
                    const char* var_name,
                    VARIABLE_TYPE type=SCALAR, PLACEMENT place=NODE, size_t vsize=1U );

    PropertyHandle( const PropertyHandle& op ); 

    ~PropertyHandle();

    // assignment
    PropertyHandle&  operator=( const PropertyHandle& op );
    PropertyHandle&  operator=( double64 val );
    PropertyHandle&  operator=( const ScalarVariable& s );
    PropertyHandle&  operator=( const VectorVariable<dim>& v );
    PropertyHandle&  operator=( const TensorVariable<dim>& t );
    PropertyHandle&  operator=( const std::vector<VectorVariable<dim> >& vc ); 
    PropertyHandle&  operator=( const std::vector<TensorVariable<dim> >& ts ); 

    // "self-assign" calculation result to PropertyHandle on the left
    PropertyHandle&  operator+=( double64 val );
    PropertyHandle&  operator-=( double64 val );
    PropertyHandle&  operator*=( double64 val );
    PropertyHandle&  operator/=( double64 val );

    // "self-assign" results from spatially variable calculations
    PropertyHandle&  operator+=( const PropertyHandle& op );
    PropertyHandle&  operator-=( const PropertyHandle& op );
    PropertyHandle&  operator*=( const PropertyHandle& op );
    PropertyHandle&  operator/=( const PropertyHandle& op );
    
    // standard math
    void          Squared();
    void          Sqrt();
    void          Ln();      // natural logarithm
    void          Log10();   // decadic logarithm
    void          Exp();
    void          Pow( double64 raised_to );
    void          ZapNAN( double64 with );
    void          Sin();
    void          Cos();
    void          Tan();
    void          Acos();
    void          Asin();
    void          Atan();
    
    // various operations
    const char*   VariableName() const;
    const csmp::Index&  Key() const;
    void           Range( double64& omin, double64& omax ) const;
    bool           IsWithinRange() const;
    VARIABLE_FLAG  OutputCondition() const;
    void           OutputCondition( VARIABLE_FLAG c );
    
    void Out() const { Out(std::cout); }
    void Out(std::ostream& os) const;
    void Out( const char* text_file_name ) const;
  
  private:
    Model<dim>&           super_group;
    std::string           group_name; 
    Region<dim>&          group;
    std::string           var_name;
    csmp::Index           key_;
    VARIABLE_FLAG         flag_output;
    bool                  new_variable_created;
};

>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3
/**
 
@brief PropertyHandle - a global accessor and mutator of distributed data

@author S.K. Matthaei
@author Stephen G. Roberts
@date 1999

 
@section motivation Motivation

Especially for test computations, one would like write arithmetic
expressions which apply to distributed physical variables. To make
these expressions concise it is useful to replace lengthy variable 
names by a notation of choice, for instance, to define a variation
of porosity with depth:
 
phi = phi - 1 / z^2
 
Finally, there is a need to create variables for such mathematical 
expressions at runtime.
 
 
@section design Design Intent

The PropertyHandle was designed following the motivation from above, but 
also by using the existing Operand classs as a model. Thus, an PropertyHandle shares
many of the Operand interfaces, but as a difference, it will currently
apply to the associated variable on all the elements, nodes, or constraint 
points of a mesh, depending on where the PropertyHandle variable is placed. 

The desired capability of writing arithmetic expressions is implemented
via C++ operator overloading. Notably, there is no assignment operator for
operands to operands, since the expected result of a C++ assignment 
would imply a complete overwriting of the lefthand physical variable with 
the one on the right of the assignment operator. 

Self-assignment is provided to write expressions which are efficient,
since they do not require the creation of temporary variables. Thus, the
expression 
 
phi = phi - 1 / z^2
 
can be replaced by
 
phi -= 1 / z^2
 
which avoids the creation of a temporary variable for the righthand side
which would subsequently be assigned to 'phi' on the left.
 
 
@section applicability Applicability

PropertyHandles can be applied anywhere in a main program or a high-level 
class object or function which has access to the CSMP library. To create
a PropertyHandle its constructor takes arguments which allow it to associate 
itself with an existing variable in the Model or to request a new 
Model variable which will be globally known until the PropertyHandle is
destructed again. If a new variable is created, it will be destroyed when 
the destructor of the PropertyHandle is called.
 
 
@section participants Participants

The PropertyHandle has references to the Mesh and Memory managers, and
the PropertyDatabase and it stores the csmp::Index of the associated 
physical variable for rapid variable access. 
 
 
@section collaborations Collaborations

The PropertyHandle collaborates with the Model and its contained 
managers in order to carry out its arithmetic operations.
 
 
@section consequences Consequences

PropertyHandles allow the user to add, subtract, multiply, and divide
distributed physical variables among eachother, to assign values to
them, or to add, subtract, multiply, and divide CSP basic variables to
the distributed physical variables.
 
 
@section implementation Implementation

For clarity, the placement and type of the physical variables which
are associated with PropertyHandles are used to impose restrictions on 
whether PropertyHandles can be used together in arithmetic expressions.
Only if the placement and variable type is the same this can be
done. 

The user is requested to make recommendations of how PropertyHandles can
be made more flexibility while avoiding ambiguities of assigning 
variables of different type to one another. 
 
 
@section examples Application Examples

Below the full source code for the porosity example from above is listed:
 
@code
// existing property using default arguments 
PropertyHandle  phi( current_supergroup, "porosity", SCALAR, ELEMENT );
phi = 0.4;

// new property 
PropertyHandle  z( current_supergroup, "depth", SCALAR, ELEMENT );

// initialize the new property using an Interrelation subclass
AssignDepthToProperty  make_depth( reference_to_property_database );
current_supergroup.Apply( make_depth );

// desired arithmetic expression
phi -= 1 / z^2;
phi.OutMatlab();
@endcode
 
In this example, the existing scalar element variable 'porosity' is assigned
to the PropertyHandle 'phi' and intialized to 40%. The construction of the 
PropertyHandle 'z' prompts the generation of a new variable "depth" as
a scalar placed on the element and initialized to NAN (not a number). 
The interrelation AssignDepthToProperty assigns the z-coordinate of the 
element centers to the new variable 'depth'. Now the following arithmetic 
expression has the necessary input values and can be applied. The result
is output to a text file in Matlab format. 

PropertyHandles also allow to create temporary variables, then apply
standard mathematical operations on these, and subsequently output 
the results to a file for visualization. In the following example, the 
natural logarithm of the variable porosity is output to a textfile: 
 
@code
PropertyHandle  phi( current_supergroup, "porosity", SCALAR, ELEMENT );
PropertyHandle  log_phi( current_supergroup, "log porosity", SCALAR, ELEMENT );

// assign the porosity to the new variable without direct assignment
log_phi  = 0.0;
log_phi += phi;
// logarithmitize the values

log_phi.Log();
log_phi.Out( logarithmitized_porosity );

// delete the temporary variable
log_phi.~PropertyHandle();
@endcode

@todo !!! SKM: fix math function calls and catch IEEE exceptions 
@todo SKM: deal with new variable placements
@todo (1) Crashes when it wants to create a REGION property (A)
@todo (3) A range of switch statements must be extended for type REGION, FACE, INTER_FACE (A)
 
*/
template<size_t dim>
class PropertyHandle {
  public:
    PropertyHandle( Model<dim>& sg, const char* var_name,
                    VARIABLE_TYPE type=SCALAR, PLACEMENT place=NODE, size_t vsize=1U );

    PropertyHandle( Model<dim>& sg, const char* target_group, 
                    const char* var_name,
                    VARIABLE_TYPE type=SCALAR, PLACEMENT place=NODE, size_t vsize=1U );

    PropertyHandle( const PropertyHandle& op ); 

    ~PropertyHandle();

    // assignment
    PropertyHandle&  operator=( const PropertyHandle& op );
    PropertyHandle&  operator=( double64 val );
    PropertyHandle&  operator=( const ScalarVariable& s );
    PropertyHandle&  operator=( const VectorVariable<dim>& v );
    PropertyHandle&  operator=( const TensorVariable<dim>& t );
    PropertyHandle&  operator=( const std::vector<VectorVariable<dim> >& vc ); 
    PropertyHandle&  operator=( const std::vector<TensorVariable<dim> >& ts ); 

    // "self-assign" calculation result to PropertyHandle on the left
    PropertyHandle&  operator+=( double64 val );
    PropertyHandle&  operator-=( double64 val );
    PropertyHandle&  operator*=( double64 val );
    PropertyHandle&  operator/=( double64 val );

    // "self-assign" results from spatially variable calculations
    PropertyHandle&  operator+=( const PropertyHandle& op );
    PropertyHandle&  operator-=( const PropertyHandle& op );
    PropertyHandle&  operator*=( const PropertyHandle& op );
    PropertyHandle&  operator/=( const PropertyHandle& op );
    
    // standard math
    void          Squared();
    void          Sqrt();
    void          Ln();      // natural logarithm
    void          Log10();   // decadic logarithm
    void          Exp();
    void          Pow( double64 raised_to );
    void          ZapNAN( double64 with );
    void          Sin();
    void          Cos();
    void          Tan();
    void          Acos();
    void          Asin();
    void          Atan();
    
    // various operations
    const char*   VariableName() const;
    const csmp::Index&  Key() const;
    void           Range( double64& omin, double64& omax ) const;
    bool           IsWithinRange() const;
    VARIABLE_FLAG  OutputCondition() const;
    void           OutputCondition( VARIABLE_FLAG c );
    
    void Out() const;
    void Out( const char* text_file_name ) const;
  
  private:
    Model<dim>&           super_group;
    std::string           group_name; 
    Region<dim>&          group;
    std::string           var_name;
    csmp::Index           key_;
    VARIABLE_FLAG         flag_output;
    bool                  new_variable_created;
};

} // csmp

#endif




