#ifndef CSMP_INTERRELATION_H
#define CSMP_INTERRELATION_H

#include "CSMP_definitions.h"
#include "Operand.h"
#include "MeshManager.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Model;
template<uint32_t,template<uint32_t> class> class ModelSubDomain;

template<uint32_t dim>
class Interrelation {
  public:
    explicit Interrelation( const PropertyDatabase<dim>& p );
    virtual ~Interrelation(); 
    
    void          Name( const char* s );

    PLACEMENT     ApplicationLevel() const;

    void          Apply( Model<dim>& );                     
                         
    /// a simplex is any geometrical association of points
    template<template<uint32_t> class CELL>
    void          Apply( ModelSubDomain<dim,CELL>& );                     

    virtual void  Calculate();

    bool          ResultWithinRange();

    typename std::map<std::string,Operand<dim> >::const_iterator  OperandsBegin() const;

    typename std::map<std::string,Operand<dim> >::const_iterator  OperandsEnd() const;
      
    typename std::map<std::string,Operand<dim> >::iterator  Result();
      
      
  protected:
    Interrelation();
    Operand<dim>&   GlobalProperty( const char* var );
    void            ResultProperty( const char* var );
    void            OutputCondition( Operand<dim>& var, VARIABLE_FLAG output );

    const PropertyDatabase<dim>&  p_ref;
    std::string                   name_, result_property_;
    PLACEMENT                     application_level_;
    typename std::map<std::string,Operand<dim> >  operand_list_;
    typename std::map<std::string,Operand<dim> >::iterator result_;
};



/**
@class Interrelation Interrelation "main_library/Interrelation.h"

@author S.K. Matthaei
@author Stephen G. Roberts
@date 1999


@section motivation Motivation

Traditionally, locally resolvable interdependencies of dependent model 
variables are computed in so-called 'post-processing' operations which 
involve looping through arrays of nodes and elements, applying mathematical 
relationships to compute the variable(s) of interest. This procedure is
error prone, firstly, because of the danger of writing beyond array 
bounds which the top level user is not necessarily familiar with; 
secondly, the results of such calculations are not tested before they are
written into the global arrays; and, thirdly, if one has many of such
calculation procedures inside of a time-stepping loop, it becomes hard
to track which post-processing procedure depends on a previous one,
to determine an execution order which must be adhered to. 

The Interrelation class represents an attempt to overcome these dangers
by encapsulating post-processing operations into modules. This design 
achieves the following objectives: 

Firstly, the loop operations and variable access is embedded into
the Interrelation class. The user must only associate an Operand
object with the variable which he can subsequently use in
the mathematical expressions. 

Secondly, The result of such a calculation is tested before it is 
stored back into the MemoryManager. This test determines whether
the result lies in a physically meaningful range for which other 
Interrelations that use the variable as an input, have been tested
(This is essential to ascertain the correctness of a complex
simulation).  

Thirdly, Interrelations always modify only a single variable. This
restriction allows them to be managed by an InterrelationManager
which automatically determines an Interrelation computation
sequence in time-stepping loops. 

There is also a need to used node variables in the computation of 
element properties and so forth. For these type of calculations, the
Interrelation class has the capability of taking property averages
or apply other useful strategies through the use of Operand
objects. 

 
 
@section design Design Intent
 

In the CSMP modeling code, all calculations that do not involve the 
inversion of a global solution matrix are done by passing 
Interrelation subclass instances to Model objects.  An example of 
an Interrelation would be the calculation of the hydraulic conductivity 
from the permeability by dividing permeability by viscosity and assigning 
the result to the variable conductivity.  The code that does this is 
provided as the class 'ConductivityFromPerm' as defined in a header file 
with the same name.  
  
 
@section applicability Applicability
 

All calculations in a CSP model which can be done in an element by element
fashion or just by looking at node or constraint point variables only,
should be encapsulated into Interrelations. Exceptions from this rule 
are global expressions which one may want to write using PropertyHandle
objects. Also, for the special case where a well-understood process 
modifies several variables simultaneously, which would require multiple
interrelations, a visitor object may be better suited to achieve the 
computation goal.  

Importantly, Interrelations can be restricted such that they apply
only to subregions of a model identified as Region objects. This 
capability allows, for instance, to compute a magmatic fluid source 
within a cooling granite without the overhead of traversing the 
entire mesh.
 
 
@section structure Structure
 

The Interrelation class is a base class from which the user inherits 
a subclass object for instance 'ConductivityFromPerm' which 
contains the mathematical expression which the user wants to apply.
Inheriting such a subclass involves the following steps: 

Firstly, the user must define a notation for his mathematical expression
by associating a Operand objects with the global physical variables of 
interest (for the above example): 
 
@code
class ConductivityFromPerm : public Interrelation {
Operand& K; // hydraulic conductivity
Operand& k; // permeability
Operand& V; // dynamic viscosity
...
@endcode

Secondly, the user must register the variables with the Model,
define a name of the Interrelation, the output property, and the write 
permission flag. This is done in the constructor of the Interrelation 
subclass: 

@code
ConductivityFromPerm::ConductivityFromPerm( const Property& p ) 
      : Interrelation(p),
        k( GlobalProperty("permeability") ),
        K( GlobalProperty("conductivity") ),
        V( GlobalProperty("viscosity") )
 {
    Name(" ConductivityFromPerm");
    ResultProperty("conductivity");
    OutputConditions( K, PLAIN );
 }
 @endcode

k, K, and V are associated with the global physical variables
identified by their spelled out name, 'conductivity' is defined as 
the property which will be modified by the Interrelation, but only
where it has the flag PLAIN.  

Thirdly, the user must overwrite the pure virtual method Calculate()
of the Interrelation base class, defining the mathematical expression
he wants to apply n order to compute the target variable: 
 
@code
inline void ConductivityFromPerm::Calculate()
 {
    K = 1.0 / viscosity;
    K *= k;
 }
@endcode
  
Here Calculate() is defined as  a C++ inline function (to reduce 
computation overhead), which first assigns 1/viscosity to the
global physical variable conductivity. Then the latter is 
multiplied by the permeability according to Darcy's law:  
 
K = k/v
 
Please refer to the documentation of the Operand class to see 
how to write mathematical expressions involving Operand objects.
 
 
@section participant Participants

An Interrelation object always contains Operands and a reference to the 
PropertyDatabase in order to register and manage the physical variables.
 
 
@section collaborations Collaborations
 
The Interrelation collaborates with Model and Region objects in 
order to compute the desired properties. To apply an Interrelation to
a model just submitt it to the Model: 
 
ConductivityFromPerm  my_interrelation( var_database_reference );
model.Apply( my_interrelation );
 
 
 
@section consequences Consequences

To make use of the Interrelation capability of CSP one must define and 
compile Interrelation subclasses. A small library of these is provided
but the intend is to build a large moderated internet repository of 
Interrelation objects from which all CSP users can profit and to which
they can submit their own Interrelation subclasses.   
 
 
@section implementation Implementation
 

The Interrelation uses a conditional loop within its Apply() interface
to establish the input values for each calculation. It tests the 
calculation results with ResultWithinRange(), before these are 
communicated back to the Model object. 

If an interrelation subclass involves scalar, vector, and/or tensor
variables simultaneously, you can access these variables directly
by using the AssignTo() method of the Operand class. This allows you
to make use of the rich interfaces of the Scalar-, Vector-, and 
TensorVariable classes.
 
 
@section examples Application Examples

As already discussed above, one must define and compile Interrelation
subclasses to make use of the Interrelation capability of CSP. The steps
involved in this process are described in more detail below: 

   1. Test the interrelation among variables that you would like to express. 
   This may involve graphing and designing the Interrelation in matlab or 
   Mathematica or doing a curve fit through some experimental data within
   Excel. You will also have to identify the value range for each of the 
   involved variables. Within this value range your mathematical 
   formulation must produce meaningful results. This value range should be 
   at least as wide as the value range that you specified in the variable 
   database text file from which the programm initializes its variable 
   storage. This is extremely important since it is the only way of 
   ascertaining that the CSP code produces correct results once your model 
   becomes complex and contains several interrelations.  

   
   2. Define your notation and make each symbol a reference to an 
   Operand object. These references should 
   be private data members of your subclass which you define through public 
   derivation from the Interrelation class. Importantly, if you don't make 
   Operands references, you will create new variables and your calculation will 
   either not read any of the variable values stored in the Model 
   object, or your calculation result will not be transmitted back into 
   the Model. 
   
   3.  Register your variables in the constructor of your derived class.  In 
   C++ all references must be initialized in the constructor.  First, you 
   construct the base class, then you use the base class method 
   GlobalProperty() to register the involved variables and finally you 
   define which of these variables will store the result of your 
   calculation.  This is done with the base class method 
   DefineResultProperty().  You also have to provide the Interrelation 
   with output conditions (use method OutputConditions().  
   The use of these conditions is to prevent the Interrelation from 
   overwriting variable values that you have fixed (e.g.  DIRICH conditions) 
   or to read variables that have not yet a meaningful value.  
   
   4. Define your own version of the pure virtual member function 
   Calculate() of the Interrelation base class. This method will contain 
   your mathematical formulation of the interrelation written using a 
   normal notation involving * / + -  < > <= >= and other mathematical 
   operations as are
   defined as overloaded operators of the Operand class. Avoid the 
   construction of too many temporaries within this method. If you need 
   many, make them static such that they do not get destructed when the 
   method goes out of scope. Also, if you write something like: 
 
       K += b + c - d * 24.5;
 
   you will automatically create a 4 temporary Operands (why ? - make this 
   a little exercise).  
   
   5. Decide where to apply your Interrelation and pass it to the 
   Model which will call the Interrelation method Apply(), that 
   performs a range check on the values that your calculation produces.  

   6. You have to define the application domain, if you want to have your
   interrelation applied only to a group of elements rather than the 
   Model. Use the method RestrictApplicationTo() to achieve this. 
   
The calculations that you can perform with Interrelation 
subclasses are limited.  You can only calculate result variables that 
are on the same or a higher level in the Model storage hierarchy.  For 
instance you can calculate an element property from node, or constraint point 
property, but not the other way around. Use a visitor and the Parent references
of the Node or IntegrationPoint classes to do this. Also if you calculate an element 
property from node property values, the interrelation will use an 
average over the nodes rather than individual values in the calculation.  
If this is too limiting for you than you can overrule the baseclass methods 
with your own method.  Before you do that, however, consider to do 
your calculation as postprocessng operations of an Algorithm or via 
a Visitor object.  

@todo (1) Remove interrelations after porting functionality to PropertyHandle and visitors (A)
 */



/**

Give constant STL iterators to the first and last elements of a map of
Operands stored inside of the Interrelation. In this map, the name of the
physical variable which is associated with the Operand is used as a key
to the later. This implies that there may only be a single Operand per
physical variable. 

@return A constant iterator to the Operand map.

@section application Application

Gives the Model access to the list of Operands inside an 
Interrelation. The constant iterators are also used inside of the 
Interrelation to avoid accidential modification of any property other
than the result property. 
*/
template<uint32_t dim>
inline typename std::map<std::string,Operand<dim> >::const_iterator  
                          Interrelation<dim>::OperandsBegin() const
 { 
     return operand_list_.begin(); 
 }

template<uint32_t dim>
inline typename std::map<std::string,Operand<dim> >::const_iterator 
                            Interrelation<dim>::OperandsEnd() const
 {   
     return operand_list_.end();
 }


/**

Returns an iterator to the result Operand / property of a calculation. 

@return An iterator to an STL map where the name of each physical variable is
used as a key to store the Operands. 

@section application Application

Is used inside Apply() to modify the result property. 
*/
template<uint32_t dim>
inline typename std::map<std::string,Operand<dim> >::iterator  Interrelation<dim>::Result()
 { 
     return result_; 
 }
      


/**
 
Tells you on what level in the Model hierarchy your 
interrelation is going to be applied. The default is ELEMENT. If you do 
however carry out a calculation that only involves node variables, the 
application level is set to NODE automatically, to avoid that averages of 
variable values are taken.  

@section arguments Input Arguments

To modify the application level manually, supply the desired target object as 
first method argument. 

@return Without first argument, ApplicationLevel() returns the current setting
of the private variable that defines the application level. 

@section implementation Implementation

The application level defines the target set of objects in the Model 
hierarchy through which the Interrelation will loop in order to carry out 
the calculation. The default level is ELEMENT, but if you define 
Interrelations that only use node variables, the application 
level is automatically set to NODE. 
*/
template<uint32_t dim>
inline PLACEMENT   Interrelation<dim>::ApplicationLevel()  const 
{ return application_level_; }





/**

Uses the value range that Operands get for the associated variable from the 
PropertyDatabase to test whether the calculated result is within range. 
To modify the range inside the Operand, call on the public interface of 
the Operand class inside the Interrelation subclass constructor to avoid 
redefining the range during each call to Calculate(). 

@return The method returns the boolean variable 'true' if the value ofthe tested
physical variable is within the specified range. 

@section implementation Implementation

Calls the Operand method IsWithinRange(). 

@section application Application

Is called inside the Interrelation method Apply() on each computed 
variable value. 

@section messages Messages

If the computed value is out of range, an error is reported. 
*/
template<uint32_t dim>
inline bool Interrelation<dim>::ResultWithinRange()
 {
    return (*result_).second.IsWithinRange();
 }




/**

Pure virtual function that the user overrides to express the interrelation
between variables mathematically. A more detailed description of this 
inheritance process is given in the introduction to the class above. 

@section implementation Implementation

Within Calculate() the user has access to Operands defined
in the private data section of the Interrelation subclass. These Operands 
can be used in mathematical and logical expressions, or associated
variables can be assigned to local variables of double, scalar, vector
or tensor type. Please refer to the documentation of the Operand 
interface to learn about the overloaded operators.  

@section application Application

To calculate local (element-restricted) interdependencies of physical
variables inside of the CSP model. 

@section messages Messages

Use the method PrintValue() of the Operand class to output computed
property values to screen. 
*/
template<uint32_t dim>
inline void Interrelation<dim>::Calculate()
 {
    std::cout<<"\nInterrelation::Calculate: Please overwrite this virtual function by own implementation..."<< std::endl;
 }

} // csp

#endif
