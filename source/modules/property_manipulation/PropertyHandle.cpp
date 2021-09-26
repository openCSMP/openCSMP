#include "PropertyHandle.h"
#include "FEM_Data.h"
#include "TensorVariable.h"
#include "Node.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"
#include "CSMP_highLevelUtilities.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/**
 
Constructor associates a distributed physical variable with the
constructed PropertyHandle. If the variable already exists in the 
PropertyDatabase, the new PropertyHandle will be a handle to it. If the 
variable name is not known to the PropertyDatabase, it will be registered
as a new variable and property storage will be allocated accordingly
by the MemoryManager. 

@section arguments Input Arguments

There are four constructor arguments. The first is a reference to the 
current Model, the second is the variable name, and the remaining
two arguments specify the variable type and placement, respectively.
These last two arguments have the default values SCALAR, and NODE. 
If the variable type or placement differs from these values, they must be
specified. If a new variable shall be created the 
user must specify its type and placement. 

@section implementation Implementation 

The constructor initializes the output variable flag to PLAIN and requests 
the csmp::Index for the associated variable from the PropertyDatabase. It 
prompts the creation of a new variable if the property does not already 
exist. 
*/
template<size_t dim>
PropertyHandle<dim>::PropertyHandle( Model<dim>& sg, 
                                     const char*   var_name, 
                                     VARIABLE_TYPE ptype, 
                                     PLACEMENT     place,
                                     size_t        psize )
 : super_group(sg),
   group(sg.Region("Model")),
   group_name("Model"),
   flag_output(ANY),
   var_name(var_name)
 {
    if ( place == BOUNDARY )
      throw Exception( ERROR, "PropertyHandle<dim>(constructor)",
                       var_name, "Property handles cannot handle BOUNDARY properties" );

    // if the variable exists the Operand is associated with it
    if ( sg.Database().IsDefined( var_name ) ) 
      {
         csmp::Index  key = sg.Database().StorageKey( var_name );
         if ( ptype != SCALAR && key.type != ptype ) {
              cout <<"\nPropertyHandle<"<< dim <<">(constructor): ";
              cout <<"Type of existing variable is different from what you have anticipated: ";
              sg.Database().FlushToScreen( var_name ); 
           }
         if ( place != NODE && key.place != place ) {
              cout <<"\nPropertyHandle<"<< dim <<">(constructor): ";
              cout <<"Placement of existing variable is different from what you have anticipated: ";
              sg.Database().FlushToScreen( var_name ); 
           }
         new_variable_created = false;
      }
    else { // a new variable is created 
         sg.CreateProperty( var_name, "SI", ptype, place, psize );
         new_variable_created = true;
      }
      
    key_ = sg.Database().StorageKey(var_name);
      
 } // end







template<size_t dim>
PropertyHandle<dim>::PropertyHandle( Model<dim>& sg, 
                                     const char*   group,
                                     const char*   var_name, 
                                     VARIABLE_TYPE ptype, 
                                     PLACEMENT     place,
                                     size_t        psize )
 : super_group(sg),
   group(sg.Region(group)),
   group_name(group),
   flag_output(ANY),
   var_name(var_name)
 {
    if ( place == BOUNDARY )
      throw Exception( ERROR, "PropertyHandle<dim>(constructor)",
                       var_name, "Property handles cannot handle BOUNDARY properties" );

    // if the variable exists the Operand is associated with it
    if ( sg.Database().IsDefined( var_name ) ) 
      {
         csmp::Index  key = sg.Database().StorageKey( var_name );
         if ( ptype != SCALAR && key.type != ptype ) {
              cout <<"\nPropertyHandle<"<< dim <<">(constructor): ";
              cout <<"Type of existing variable is different from what you have anticipated: ";
              sg.Database().FlushToScreen( var_name ); 
           }
         if ( place != NODE && key.place != place ) {
              cout <<"\nPropertyHandle<"<< dim <<">(constructor): ";
              cout <<"Placement of existing variable is different from what you have anticipated: ";
              sg.Database().FlushToScreen( var_name ); 
           }
         new_variable_created = false;
      }
    else { // a new variable is created 
         sg.CreateProperty( var_name, "SI", ptype, place, psize );
         new_variable_created = true;
      }

    key_ = sg.Database().StorageKey(var_name);

 } // end ct




/** A copy-constructed PropertyHandle is associated exactly with the same variable as its input argument.
*/
template<size_t dim>
PropertyHandle<dim>::PropertyHandle( const PropertyHandle<dim>& op )
 : super_group(op.super_group),
   group_name(op.group_name),
   group(op.group),
   var_name(op.var_name),
   key_(op.key_),
   flag_output(op.flag_output),
   new_variable_created(op.new_variable_created)
{
}





/**

The destructor of a PropertyHandle linked to dynamically created variable,
will prompt the PropertyDatabase to delete the variable entry and call the
MemoryManager to delete the storage which was allocated for that 
variable. 
*/
template<size_t dim>
PropertyHandle<dim>::~PropertyHandle()
 {
    if ( new_variable_created ) 
      super_group.DeleteProperty( var_name.c_str() );
    
 } // end destructor





/** Returns the storage specification of the variable which is associated with the csmp::Operand into a corresponding structure.
*/    
template<size_t dim>
const csmp::Index&  PropertyHandle<dim>::Key() const { 
    return key_;
 }

template<size_t dim>
const char*  PropertyHandle<dim>::VariableName() const
 { return var_name.c_str(); }


/**

OutputCondition() returns  the VARIABLE_FLAG flag which the
distributed physical variable which is associated with the PropertyHandle
must have in order to allow modification by the PropertyHandle. 

@return The flag that the physical variables flag must correspond to if the
variable shall be modified.

@section application Application 

Since each CSMP scalar, vector, or tensor variable has a single or a set
of flags which indicate to the solver whether it may modify or use this 
variable instance as a constraint, the same rule applies to PropertyHandles.
Thus, only if the output flag matches the local flag of the distributed
variable, it will modify the latter. Accordingly, the variable flag may
be used to protect certain variable values from modification by 
PropertyHandles.*/    
template<size_t dim>
VARIABLE_FLAG  PropertyHandle<dim>::OutputCondition() const { return flag_output; }
    
    
/**

OutputCondition() assigns the VARIABLE_FLAG flag which the
distributed physical variable which is associated with the PropertyHandle
must have in order to allow modification by the PropertyHandle. 

@param c The flag which the output flag shall be changed into.

@attention The flag of the physical variables flag must correspond,
else the variable is not modified.

@section application Application 

Since each CSMP scalar, vector, or tensor variable has a single or a set
of flags which indicate to the solver whether it may modify or use this 
variable instance as a constraint, the same rule applies to PropertyHandles.
Thus, only if the output flag matches the local flag of the distributed
variable, it will modify the latter. Accordingly, the variable flag may
be used to protect certain variable values from modification by 
PropertyHandles. 
*/
template<size_t dim>
void  PropertyHandle<dim>::OutputCondition( VARIABLE_FLAG c ) { flag_output=c; }




   


template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator=( const PropertyHandle& op )
 {
    if ( group_name != op.group_name )
      throw Exception( ERROR, "PropertyHandle<dim>::operator=(PropertyHandle)",
                      VariableName(), "Property handles are associated with different model subdomains" );
 
    if ( &op != this ) {
         // key is not touched, thus placement and variable type must remain
         // the same as well and all other variables are retained
         if ( group_name != op.group_name )
           throw csmp::Exception( WARNING, "PropertyHandle<dim>::operator=", 
                          "operands have different targets using that of the assigned operand"  );
         group_name = op.group_name; 
 
          if ( var_name != op.var_name ) {
               throw csmp::Exception( WARNING, "PropertyHandle<dim>::operator=", 
                              "target operands have different names, no assignment was made"  );
               return *this;
            }
    
       VectorVariable<dim>  vc(flag_output,std::numeric_limits<double64>::quiet_NaN());
       TensorVariable<dim>  ts(flag_output,std::numeric_limits<double64>::quiet_NaN());
       csmp::Index  key   = super_group.Database().StorageKey( var_name.c_str() );
       csmp::Index  opkey = super_group.Database().StorageKey( op.VariableName() );
    
       // these combinations of properties are not possible   
       if ( key.place != opkey.place )
         throw Exception( ERROR, "PropertyHandle<dim>::operator=(PropertyHandle)",
                          VariableName(), "the property that shall be assigned must have the same placement as this one" );
       
        switch( key.type ) {
	         case SCALAR: 
	              switch( key.place ) {
	                   case NODE:
	                        for ( typename vector<Node<dim>*>::iterator
	                              nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
	                          if ( (*nit)->Status( key ) == flag_output ) 
	                            (*nit)->Store( key, makeScalar( flag_output, (*nit)->Read( opkey )) );
	                     break;
	                   case ELEMENT_INTEGRATION_POINT:
	                        for ( typename vector<Element<dim>*>::iterator
	                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
	                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
	                            if ( (*eit)->Status( i, key ) == flag_output ) 
	                              (*eit)->Store( i, key, makeScalar( flag_output, (*eit)->Read( i, opkey )) );
	                     break;
	                   case ELEMENT:        
	                        for ( typename vector<Element<dim>*>::iterator
	                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
	                          if ( (*eit)->Status( key ) == flag_output ) 
	                            (*eit)->Store( key, makeScalar( flag_output, (*eit)->Read( opkey )) );
	                      break;
	                   case REGION:        
                          if ( group.Status( key ) == flag_output ) 
                            group.Store( key, makeScalar( flag_output, group.Read( opkey )) );
	                      break;
	                   default:
                         throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                              "BOUNDARY and MODEL properties are not handled by this class yet."  );
	                }
	           break;
	         case VECTOR: 
	              switch( key.place )
	                {
	                   case NODE:
	                        for ( typename vector<Node<dim>*>::iterator 
	                              nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
	                            {
	                               (*nit)->Read( opkey, vc );	                          
	                               (*nit)->Store( key, vc );
	                            }
	                     break;
	                   case ELEMENT_INTEGRATION_POINT:
	                        for ( typename vector<Element<dim>*>::iterator
	                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
	                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
	                          if ( (*eit)->Status( i, key ) == flag_output ) {
	                               (*eit)->Read( i, opkey, vc );	                          
	                               (*eit)->Store( i, key, vc );
	                            }
	                     break;
	                   case ELEMENT:        
	                        for ( typename vector<Element<dim>*>::iterator
	                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
	                            {
	                               (*eit)->Read( opkey, vc );	                          
	                               (*eit)->Store( key, vc );
	                            }
	                      break;
	                   case REGION:        
                          if ( group.Status( key ) == flag_output ) { 
                               group.Read( opkey, vc );
	                             group.Store( key, vc );
                            }     
	                      break;
	                   default:
                         throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                              "BOUNDARY and MODEL properties are not handled by this class yet."  );
	                }
	           break;
	         case TENSOR: 
	              switch( key.place )
	                {
	                   case NODE:
	                        for ( typename vector<Node<dim>*>::iterator 
	                              nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
	                            {
	                               (*nit)->Read( opkey, ts );	                          
	                               (*nit)->Store( key, ts );
	                            }
	                     break;
	                   case ELEMENT_INTEGRATION_POINT:
	                        for ( typename vector<Element<dim>*>::iterator
	                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
	                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
  	                            {
  	                               (*eit)->Read( i, opkey, ts );	                          
  	                               (*eit)->Store( i, key, ts );
  	                            }
	                     break;
	                   case ELEMENT:        
	                        for ( typename vector<Element<dim>*>::iterator
	                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
	                            {
	                               (*eit)->Read( opkey, ts );	                          
	                               (*eit)->Store( key, ts );
	                            }
	                      break;
	                   case REGION:        
                            { 
                               group.Read( opkey, ts );
	                             group.Store( key, ts );
                            }     
	                      break;
	                   default:
                         throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                              "BOUNDARY and MODEL properties are not handled by this class yet."  );
	                }
	            break;
	          default:
              throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                             "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
          }
      }
    IsWithinRange();

    return *this;
 }






/**
 
Assigns the value of its argument to the physical variable associated
with the PropertyHandle at where the variable has a flag which is equal
to the output flag. 

@section arguments Input Arguments

Variables of the types double64, ScalarVariable, VectorVariable, or
TensorVariable may be assigned. 

@return The assignment returns the value into the PropertyHandle to the left of the
equal sign. The original values are overwritten irrespective of their
original flag.  

@section implementation Implementation 

The assigment operator relies on the overloaded operators of the basic
CSMP variables (scalar, vector, and tensor) to perform the assignments. 

In the special case, where an assignment of a vector to a scalar is 
attempted, the length of the vector will be assigned. 

The attempt to assign a tensor to a vector is undefined and 
will therefore provoke an error. 

If a tensor shall be assigned to a scalar, the determinant of the tensor
will be assigned. 

@section application Application 

To write PropertyHandle expressions like: 

@code
csp_specific_op = 12.5;
@endcode

@section messages Messages 

Undefined assigments, like that of a tensor to a vector will invoke
error messages. Because these are special operations, if the length of
a vector or the determinant of a tensor is assigned to a scalar,
a warning will be issued which specifies the type of assignment made. 
*/
template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator=( double64 val )
 {
    switch( key_.type )
      {
         case SCALAR: {
              ScalarVariable  sc(flag_output,val);
              switch( key_.place )
                {
                   case NODE:
                        for ( typename vector<Node<dim>*>::iterator
                              nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          (*nit)->Store( key_, sc );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( typename vector<Element<dim>*>::iterator
                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            (*eit)->Store( i, key_, sc );
                     break;
                   case ELEMENT:        
                        for ( typename vector<Element<dim>*>::iterator
                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          (*eit)->Store( key_, sc );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
             }
           break;
         case VECTOR: {
              VectorVariable<dim>  vc(flag_output,std::numeric_limits<double64>::quiet_NaN());
              vc = val;
              switch( key_.place )
                {
                   case NODE:
                        for ( typename vector<Node<dim>*>::iterator
                              nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          (*nit)->Store( key_, vc );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( typename vector<Element<dim>*>::iterator
                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            (*eit)->Store( i, key_, vc );
                     break;
                   case ELEMENT:        
                        for ( typename vector<Element<dim>*>::iterator
                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          (*eit)->Store( key_, vc );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
             }
           break;
         case TENSOR: {
              TensorVariable<dim>  ts(flag_output,std::numeric_limits<double64>::quiet_NaN());
              ts = val;  
              switch( key_.place )
                {
                   case NODE:
                        for ( typename vector<Node<dim>*>::iterator
                              nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          (*nit)->Store( key_, ts );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( typename vector<Element<dim>*>::iterator
                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            (*eit)->Store( i, key_, ts );
                     break;
                   case ELEMENT:        
                        for ( typename vector<Element<dim>*>::iterator
                              eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          (*eit)->Store( key_, ts );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
              }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    return *this;
 }



template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator=( const ScalarVariable& s )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                         sc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    VectorVariable<dim>                     vc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    TensorVariable<dim>                     ts(flag_output,std::numeric_limits<double64>::quiet_NaN());
    sc = s;
    vc = s;
    ts = s;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output ) (*nit)->Store( key, sc );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) (*eit)->Store( i, key, sc );
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) (*eit)->Store( key, sc );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case VECTOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          (*nit)->Store( key, vc );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) (*eit)->Store( i, key, vc );
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          (*eit)->Store( key, vc );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          (*nit)->Store( key, ts );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            (*eit)->Store( i, key, ts );
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          (*eit)->Store( key, ts );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    return *this;
 }


template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator=( const VectorVariable<dim>& vc )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                         sc(flag_output,vc.Length());

    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR: 
              cout <<"\nPropertyHandle::operator=: Assigning vector length to scalar !"<< endl;
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output ) (*nit)->Store( key, sc );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) (*eit)->Store( i, key, sc );
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) (*eit)->Store( key, sc );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case VECTOR: {
              VectorVariable<dim>  temp;
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ ) {
                             (*nit)->Read( key, temp );
                             for ( size_t j=0U; j<dim; j++ )
                               if ( (*nit)->Status( key, j ) == flag_output ) temp(j) = vc[j];
                             (*nit)->Store( key, temp );
                          }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ ) {
                               (*eit)->Read( i, key, temp );
                               for ( size_t j=0U; j<dim; j++ )
                                 if ( (*eit)->Status(i, key, j ) == flag_output ) temp(j) = vc[j];
                               (*eit)->Store( i, key, vc );
                            }
                            
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ ) {
                             (*eit)->Read( key, temp );
                             for ( size_t j=0U; j<dim; j++ )
                               if ( (*eit)->Status( key, j ) == flag_output ) temp(j) = vc[j];
                             (*eit)->Store( key, vc );
                          }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
             }
           break;
         case TENSOR: 
              cout <<"\nPropertyHandle::operator=: No rule exists to assign 'vector' to 'tensor'. ";
              cout <<"No assignments were made."<< endl;
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    return *this;
 }


/**
    Assigns the values of the tensor variable to the distributed property.
*/
template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator=( const TensorVariable<dim>& ts )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable  sc(flag_output,strtod("NAN",NULL));
    sc() = ts.Determinant();
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              cout <<"\nPropertyHandle::operator=: Assigning determinant of tensor variable to scalar."<< endl; 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output ) (*nit)->Store( key, sc );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) (*eit)->Store( i, key, sc );
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) (*eit)->Store( key, sc );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case VECTOR: 
              cout <<"\nPropertyHandle::operator=: No rule was supplied to assign tensor to vector. "; 
              cout <<"Nothing was done."<< endl; 
           break;
         case TENSOR: 
              // NB: using only the status of the first diagonal element of the tensor
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key, 0 ) == flag_output ) (*nit)->Store( key, ts );
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key, 0 ) == flag_output ) (*eit)->Store( i, key, ts );
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key, 0 ) == flag_output ) (*eit)->Store( key, ts );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    return *this;
 }






template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator=( const std::vector<VectorVariable<dim> >& vc ) 
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    if ( key.type != VECTOR ) {
         throw csmp::Exception( ERROR, "PropertyHandle<dim>::operator=( vector of vectors )", 
                                    "Type mismatch; PropertyHandle does not contain vectors" );
      }
    switch ( key.place ) {
         case NODE: 
              if ( group.Nodes() != vc.size() )
                throw csmp::Exception( FATAL_ERROR, "PropertyHandle<dim>::operator=( vector of vectors )",
                                             "vector size does not match number of Nodes" ); 
           break;
         case ELEMENT_INTEGRATION_POINT:
              if ( group.IntegrationPoints() != vc.size() )
                throw csmp::Exception( FATAL_ERROR, "PropertyHandle<dim>::operator=( vector of vectors )",
                                             "vector size does not match number of IntegrationPoints" ); 
           break;
         case ELEMENT:
              if ( group.Elements() != vc.size() )
                throw csmp::Exception( FATAL_ERROR, "PropertyHandle<dim>::operator=( vector of vectors )",
                                             "vector size does not match number of Elements" ); 
            break;
         default:
             throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                  "REGION and MODEL properties are not handled by this class yet."  );
      }
  
    FEM_Data<VectorVariable<dim> >  var_data( key.place, vc );

    super_group.InputVariableFrom( super_group.Database().Name(key), var_data );
      
    return *this; 

 } // end operator



template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator=( const std::vector<TensorVariable<dim> >& ts ) 
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    if ( key.type != TENSOR ) {
         throw csmp::Exception( ERROR, "PropertyHandle<dim>::operator=( vector of tensors )", 
                                    "Type mismatch; PropertyHandle does not contain tensors" );
      }
    switch ( key.place ) {
         case NODE: 
              if ( group.Nodes() != ts.size() )
                throw csmp::Exception( FATAL_ERROR, "PropertyHandle<dim>::operator=( vector of tensors )",
                                             "vector size does not match number of Nodes" ); 
           break;
         case ELEMENT_INTEGRATION_POINT:
              if ( group.IntegrationPoints() != ts.size() )
                throw csmp::Exception( FATAL_ERROR, "PropertyHandle<dim>::operator=( vector of tensors )",
                                             "vector size does not match number of IntegrationPoints" ); 
           break;
         case ELEMENT:
              if ( group.Elements() != ts.size() )
                throw csmp::Exception( FATAL_ERROR, "PropertyHandle<dim>::operator=( vector of tensors )",
                                             "vector size does not match number of Elements" ); 
            break;
         default:
             throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                  "REGION and MODEL properties are not handled by this class yet."  );
      }
  
    FEM_Data<TensorVariable<dim> >  var_data( key.place, ts );

    super_group.InputVariableFrom( super_group.Database().Name(key), var_data );
      
    return *this; 

 } // end operator







/**
 
The operators '+=', '-=', '*=', and '/=' perform combined arithmetic and
assignment operations but avoid the generation of temporary variables as
would be generated, if the operatios were performed separately. 

@param val The argument to the right of the combined operators must be a floating point
value. 

@return The floating point value on the right modifies the PropertyHandle value on
the left according to the arithmetic operator used.

@section implementation Implementation 

The operators are based on the overloaded operators of the CSP basic 
variables, the ScalarVariable, VectorVariable, and TensorVariable. 

@section application Application 

To add a floating point value to, subtract from, multiply, or divide by the 
physical variable which is associated with the PropertyHandle, write
expressions like: 

@code
my_operand += 5.2;
@endcode
 */
template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator+=( double64 val )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                         sc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    VectorVariable<dim>                     vc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    TensorVariable<dim>                     ts(flag_output,std::numeric_limits<double64>::quiet_NaN());
    size_t i, j; 
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() += val;
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key, sc );
                                 sc() += val;
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() += val;
                               (*eit)->Store( key, sc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case VECTOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) += val;
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i) += val;
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                  case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) += val;
                               (*eit)->Store( key, vc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) += val;
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) += val;
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) += val;
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    
    return *this;

 } // end +=



template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator-=( double64 val )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable      sc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    VectorVariable<dim>  vc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    TensorVariable<dim>  ts(flag_output,std::numeric_limits<double64>::quiet_NaN());
    size_t i, j; 
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() -= val;
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key, sc );
                                 sc() -= val;
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() -= val;
                               (*eit)->Store( key, sc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case VECTOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) -= val;
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i) -= val;
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) -= val;
                               (*eit)->Store( key, vc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) -= val;
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) -= val;
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) -= val;
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    return *this;

 } // end -=



/**
    Where the variable is flagged that the output flag, the default value of which is ANY,
    its value is multiplied by the argument value.
*/
template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator*=( double64 val )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable      sc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    VectorVariable<dim>  vc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    TensorVariable<dim>  ts(flag_output,std::numeric_limits<double64>::quiet_NaN());
    size_t i, j; 
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() *= val;
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key, sc );
                                 sc() *= val;
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() *= val;
                               (*eit)->Store( key, sc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case VECTOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) *= val;
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i) *= val;
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) *= val;
                               (*eit)->Store( key, vc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) *= val;
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) *= val;
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) *= val;
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    return *this;

 } // end *=



template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator/=( double64 val )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable      sc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    VectorVariable<dim>  vc(flag_output,std::numeric_limits<double64>::quiet_NaN());
    TensorVariable<dim>  ts(flag_output,std::numeric_limits<double64>::quiet_NaN());   
    size_t i, j; 
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() /= val;
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key, sc );
                                 sc() /= val;
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() /= val;
                               (*eit)->Store( key, sc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case VECTOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) /= val;
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i) /= val;
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i) /= val;
                               (*eit)->Store( key, vc );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR: 
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) /= val;
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) /= val;
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ ) for ( j=0; j<dim; j++ ) ts(i,j) /= val;
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
      }
    IsWithinRange();
    return *this;

 } // end operator/=




/**
 
The operators '+=', '-=', '*=', and '/=' perform combined arithmetic and
assignment operations involving two PropertyHandles but avoid the generation 
of temporary variables as would be generated, if the operatios were 
performed separately. 

@param op The argument to the right of the combined operators must be a
PropertyHandle.

@return The value of the physical variable of the PropertyHandle on the right modifies
the PropertyHandle value on the left according to the arithmetic operator
used. 

@section implementation Implementation 

The operators are based on the overloaded operators of the CSP basic 
variables, the ScalarVariable, VectorVariable, and TensorVariable. 

@section application Application 

To add a floating point value to, subtract from, multiply, or divide by the 
physical variable which is associated with the PropertyHandle, write
expressions like: 

@code
my_operand += 5.2;
@endcode

@section messages Messages 

Operations involving variables of different type or different placement 
cannot be performed. In this case an error will be reported. 
*/
template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator+=( const PropertyHandle<dim>& op )
 {
    if ( group_name != op.group_name )
      throw Exception( ERROR, "PropertyHandle<dim>::operator+=(PropertyHandle)",
                      VariableName(), "Property handles are associated with different model subdomains" );

    csmp::Index  key   = super_group.Database().StorageKey( var_name.c_str() );
    csmp::Index  opkey = super_group.Database().StorageKey( op.VariableName() );
 
    if ( key.type  != opkey.type ) {
         throw csmp::Exception( ERROR, "PropertyHandle::operator+=", "Operands are not of the same type",
                                "Nothing was done.");
         return *this;
      }

    //if ( key.place != opkey.place ) {
    //     throw csmp::Exception( ERROR, "PropertyHandle::operator+=", "Operands have not the same placement",
    //                            "Nothing was done.");
    //     return *this;
    //  }

    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc1, sc2;
    VectorVariable<dim>                     vc1, vc2;
    TensorVariable<dim>                     ts1, ts2;

    if (key.place == opkey.place)
      switch( opkey.type )
      {
         case SCALAR:    
              switch( opkey.place )
              {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                            {
                               (*nit)->Read( key,    sc1 );
                               (*nit)->Read( opkey, sc2 );
                               sc1 += sc2;
                               (*nit)->Store( key, sc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key,    sc1 );
                                 (*eit)->Read( i, opkey, sc2 );
                                 sc1 += sc2;
                                 (*eit)->Store( i, key, sc1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key,    sc1 );
                               (*eit)->Read( opkey, sc2 );
                               sc1 += sc2;
                               (*eit)->Store( key, sc1 );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
               }
            break;
         case VECTOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc1 );
                               (*nit)->Read( opkey, vc2 );
                               vc1 += vc2;
                               (*nit)->Store( key, vc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc1 );
                                 (*eit)->Read( i, opkey, vc2 );
                                 vc1 += vc2;
                                 (*eit)->Store( i, key, vc1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc1 );
                               (*eit)->Read( opkey, vc2 );
                               vc1 += vc2;
                               (*eit)->Store( key, vc1 );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
              break;
         case TENSOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts1 );
                               (*nit)->Read( opkey, ts2 );
                               ts1 += ts2;
                               (*nit)->Store( key, ts1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts1 );
                                 (*eit)->Read( i, opkey, ts2 );
                                 ts1 += ts2;
                                 (*eit)->Store( i, key, ts1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts1 );
                               (*eit)->Read( opkey, ts2 );
                               ts1 += ts2;
                               (*eit)->Store( key, ts1 );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
              break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       }

    else if (key.place == ELEMENT && opkey.place == NODE)
    {
        PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateNodeToElementProperty(op.VariableName(), (var_name + "_temp").c_str() );
        (*this) += tempElementVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempNodeVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToNodeProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) += tempNodeVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT_INTEGRATION_POINT)
    {
        PropertyHandle<dim> tempNodeVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateIntegrationPointToNodeProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) += tempNodeVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == NODE)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateNodeToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) += tempIPVariable;
    }
    else if (key.place == ELEMENT && opkey.place == ELEMENT_INTEGRATION_POINT)
    {
        PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateIntegrationPointToElementProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) += tempElementVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) += tempIPVariable;
    }
    else if (key.place == FACET_INTEGRATION_POINT && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempFIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToFacetIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) += tempFIPVariable;
    }
    else
    {
        throw csmp::Exception( ERROR, "PropertyHandle::operator+=", "This placement combination has not been implemented yet",
                                      "Nothing was done.");
        return *this;
    }

    IsWithinRange();
    return *this;
 
 } // end operator+=
 


template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator-=( const PropertyHandle<dim>& op )
 {
     if ( group_name != op.group_name )
      throw Exception( ERROR, "PropertyHandle<dim>::operator-=(PropertyHandle)",
                      VariableName(), "Property handles are associated with different model subdomains" );

    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    csmp::Index  opkey = super_group.Database().StorageKey( op.VariableName() );

    if ( key.type  != opkey.type ) {
         throw csmp::Exception( ERROR, "PropertyHandle::operator+=", "Operands are not of the same type",
                                "Nothing was done.");
         return *this;
      }
    //if ( key.place != opkey.place ) {
    //     throw csmp::Exception( ERROR, "PropertyHandle::operator+=", "Operands have not the same placement",
    //                            "Nothing was done.");
    //     return *this;
    //  }
    
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc1, sc2;
    VectorVariable<dim>                     vc1, vc2;
    TensorVariable<dim>                     ts1, ts2;

    if (key.place == opkey.place)
      switch( opkey.type )
      {
         case SCALAR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key,    sc1 );
                               (*nit)->Read( opkey, sc2 );
                               sc1 -= sc2;
                               (*nit)->Store( key, sc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,    sc1 );
                                 (*eit)->Read( i, opkey, sc2 );
                                 sc1 -= sc2;
                                 (*eit)->Store( i, key, sc1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key,    sc1 );
                               (*eit)->Read( opkey, sc2 );
                               sc1 -= sc2;
                               (*eit)->Store( key, sc1 );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc1 );
                               (*nit)->Read( opkey, vc2 );
                               vc1 -= vc2;
                               (*nit)->Store( key, vc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc1 );
                                 (*eit)->Read( i, opkey, vc2 );
                                 vc1 -= vc2;
                                 (*eit)->Store( i, key, vc1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc1 );
                               (*eit)->Read( opkey, vc2 );
                               vc1 -= vc2;
                               (*eit)->Store( key, vc1 );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts1 );
                               (*nit)->Read( opkey, ts2 );
                               ts1 -= ts2;
                               (*nit)->Store( key, ts1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts1 );
                                 (*eit)->Read( i, opkey, ts2 );
                                 ts1 -= ts2;
                                 (*eit)->Store( i, key, ts1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts1 );
                               (*eit)->Read( opkey, ts2 );
                               ts1 -= ts2;
                               (*eit)->Store( key, ts1 );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       }
    else if (key.place == ELEMENT && opkey.place == NODE)
    {
        PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateNodeToElementProperty(op.VariableName(), (var_name + "_temp").c_str() );
        (*this) -= tempElementVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempNodeVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToNodeProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) -= tempNodeVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT_INTEGRATION_POINT)
    {
        PropertyHandle<dim> tempNodeVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateIntegrationPointToNodeProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) -= tempNodeVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == NODE)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateNodeToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) -= tempIPVariable;
    }
    else if (key.place == ELEMENT && opkey.place == ELEMENT_INTEGRATION_POINT)
    {
        PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateIntegrationPointToElementProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) -= tempElementVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) -= tempIPVariable;
    }
    else
    {
        throw csmp::Exception( ERROR, "PropertyHandle::operator-=", "This placement combination has not been implemented yet",
                                      "Nothing was done.");
        return *this;
    }

    IsWithinRange();
    return *this;
 
 } // end operator-=
 


template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator*=( const PropertyHandle<dim>& op )
 {
    if ( group_name != op.group_name )
      throw Exception( ERROR, "PropertyHandle<dim>::operator*=(PropertyHandle)",
                      VariableName(), "Property handles are associated with different model subdomains" );

    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    csmp::Index  opkey = super_group.Database().StorageKey( op.VariableName() );


    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                     sc1, sc2;
    VectorVariable<dim>                     vc1, vc2;
    TensorVariable<dim>                     ts1, ts2;

    if (key.place == opkey.place)
      switch( opkey.type )
      {
         case SCALAR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key,    sc1 );
                               (*nit)->Read( opkey, sc2 );
                               sc1 *= sc2;
                               (*nit)->Store( key, sc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,    sc1 );
                                 (*eit)->Read( i, opkey, sc2 );
                                 sc1 *= sc2;
                                 (*eit)->Store( i, key, sc1 );
                              }
                     break;
                   case ELEMENT:  
                     if( key.type == opkey.type )
                     {
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key,    sc1 );
                               (*eit)->Read( opkey, sc2 );
                               sc1 *= sc2;
                               (*eit)->Store( key, sc1 );
                           }
                     }
                     else if( key.type == VECTOR && opkey.type == SCALAR )
                     {
                       for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ ) // no flag manipulations for props of diff type
                         {
                           (*eit)->Read( key,    vc1 );
                           (*eit)->Read( opkey,  sc1 );
                           vc1 *= sc1;
                           (*eit)->Store( key, vc1 );
                         }
                     }
                     else
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(*=)", "Type combination not implemented"  );
                     break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc1 );
                               (*nit)->Read( opkey, vc2 );
                               vc1 *= vc2;
                               (*nit)->Store( key, vc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc1 );
                                 (*eit)->Read( i, opkey, vc2 );
                                 vc1 *= vc2;
                                 (*eit)->Store( i, key, vc1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc1 );
                               (*eit)->Read( opkey, vc2 );
                               vc1 *= vc2;
                               (*eit)->Store( key, vc1 );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts1 );
                               (*nit)->Read( opkey, ts2 );
                               ts1 *= ts2;
                               (*nit)->Store( key, ts1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts1 );
                                 (*eit)->Read( i, opkey, ts2 );
                                 ts1 *= ts2;
                                 (*eit)->Store( i, key, ts1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts1 );
                               (*eit)->Read( opkey, ts2 );
                               ts1 *= ts2;
                               (*eit)->Store( key, ts1 );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       }
    else if (key.place == ELEMENT && opkey.place == NODE)
    {
       PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), opkey.type, key.place );
       group.InterpolateNodeToElementProperty(op.VariableName(), (var_name + "_temp").c_str() );
       (*this) *= tempElementVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT) // P. Lang EDIT
    {
        string opVarName(  op.VariableName() );
        PropertyHandle<dim> tempNodeVariable( super_group, (opVarName + "_temp").c_str(), opkey.type, key.place );
        group.ExtrapolateElementToNodeProperty( op.VariableName(), (opVarName + "_temp").c_str() );
        (*this) *= tempNodeVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT_INTEGRATION_POINT) /// @todo (1-F) key.type, key.place should be opkey.type, key.place
    {
        PropertyHandle<dim> tempNodeVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateIntegrationPointToNodeProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) *= tempNodeVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == NODE)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateNodeToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) *= tempIPVariable;
    }
    else if (key.place == ELEMENT && opkey.place == ELEMENT_INTEGRATION_POINT)
    {
        PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateIntegrationPointToElementProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) *= tempElementVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) *= tempIPVariable;
    }
    else
    {
        throw csmp::Exception( ERROR, "PropertyHandle::operator*=", "This placement combination has not been implemented yet",
                                      "Nothing was done.");
        return *this;
    }

    IsWithinRange();
    return *this;
 
 } // end operator*=
 


template<size_t dim>
PropertyHandle<dim>&  PropertyHandle<dim>::operator/=( const PropertyHandle<dim>& op )
 {
    if ( group_name != op.group_name )
      throw Exception( ERROR, "PropertyHandle<dim>::operator/=(PropertyHandle)",
                      VariableName(), "Property handles are associated with different model subdomains" );

    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    csmp::Index  opkey = super_group.Database().StorageKey( op.VariableName() );

    if ( key.type  != opkey.type ) {
         throw csmp::Exception( ERROR, "PropertyHandle::operator+=", "Operands are not of the same type",
                                "Nothing was done.");
         return *this;
      }
    //if ( key.place != opkey.place ) {
    //     throw csmp::Exception( ERROR, "PropertyHandle::operator+=", "Operands have not the same placement",
    //                            "Nothing was done.");
    //     return *this;
    //  }
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                     sc1, sc2;
    VectorVariable<dim>                     vc1, vc2;
    TensorVariable<dim>                     ts1, ts2;

    if (key.place == opkey.place)
      switch( opkey.type )
      {
         case SCALAR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key,    sc1 );
                               (*nit)->Read( opkey, sc2 );
                               sc1 /= sc2;
                               (*nit)->Store( key, sc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,    sc1 );
                                 (*eit)->Read( i, opkey, sc2 );
                                 sc1 /= sc2;
                                 (*eit)->Store( i, key, sc1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key,    sc1 );
                               (*eit)->Read( opkey, sc2 );
                               sc1 /= sc2;
                               (*eit)->Store( key, sc1 );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc1 );
                               (*nit)->Read( opkey, vc2 );
                               vc1 /= vc2;
                               (*nit)->Store( key, vc1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc1 );
                                 (*eit)->Read( i, opkey, vc2 );
                                 vc1 /= vc2;
                                 (*eit)->Store( i, key, vc1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc1 );
                               (*eit)->Read( opkey, vc2 );
                               vc1 /= vc2;
                               (*eit)->Store( key, vc1 );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( opkey.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts1 );
                               (*nit)->Read( opkey, ts2 );
                               ts1 /= ts2;
                               (*nit)->Store( key, ts1 );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts1 );
                                 (*eit)->Read( i, opkey, ts2 );
                                 ts1 /= ts2;
                                 (*eit)->Store( i, key, ts1 );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts1 );
                               (*eit)->Read( opkey, ts2 );
                               ts1 /= ts2;
                               (*eit)->Store( key, ts1 );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       }
    else if (key.place == ELEMENT && opkey.place == NODE)
    {
        PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateNodeToElementProperty(op.VariableName(), (var_name + "_temp").c_str() );
        (*this) /= tempElementVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempNodeVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToNodeProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) /= tempNodeVariable;
    }
    else if (key.place == NODE && opkey.place == ELEMENT_INTEGRATION_POINT)
    {
        PropertyHandle<dim> tempNodeVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateIntegrationPointToNodeProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) /= tempNodeVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == NODE)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateNodeToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) /= tempIPVariable;
    }
    else if (key.place == ELEMENT && opkey.place == ELEMENT_INTEGRATION_POINT)
    {
        PropertyHandle<dim> tempElementVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.InterpolateIntegrationPointToElementProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) /= tempElementVariable;
    }
    else if (key.place == ELEMENT_INTEGRATION_POINT && opkey.place == ELEMENT)
    {
        PropertyHandle<dim> tempIPVariable( super_group, ( var_name + "_temp").c_str(), key.type, key.place );
        group.ExtrapolateElementToIntegrationPointProperty( op.VariableName(), (var_name + "_temp").c_str() );
        (*this) /= tempIPVariable;
    }
    else
    {
        throw csmp::Exception( ERROR, "PropertyHandle::operator/=", "This placement combination has not been implemented yet",
                                      "Nothing was done.");
        return *this;
    }

    IsWithinRange();
    return *this;
 
 } // end operator/=
 


/**
 
Multiplies the physical variable which is associated with the PropertyHandle
by itself. If this variable is vector, or a tensor, the corresponding 
vector or matrix products will be assigned, respectively. 

@section implementation Implementation 

The method relies on the interface of the CSP basic variables 
ScalarVariable, VectorVariable, and TensorVariable. 

@section application Application 

Apart from offering an efficient way to raise variable values to 
the power of 2, Squared() implements vector and matrix multiplication
for vector and tensor variables, respectively. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Squared()
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc *= sc;
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc *= sc;
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc *= sc;
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               vc *= vc;
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 vc *= vc;
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               vc *= vc;
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               ts *= ts;
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 ts *= ts;
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               ts *= ts;
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end Squared




/**
 
Takes the square root of the values of the physical variable which is
associated with the PropertyHandle. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Sqrt()
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = sqrt( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = sqrt( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = sqrt( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=sqrt( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=sqrt( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=sqrt( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = sqrt( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = sqrt( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = sqrt( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
    
 } // end Sqrt



/**
 

Calculates the natural logarithm (base 2) of each component of the
CSP basic variable which is associated with the PropertyHandle. 

@section messages Messages 

Since the method cannot take the natural logarithm of 0 or a negative 
number, variables with such a value are not modified. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Ln()
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               if ( sc() > 0.0 ) sc() = log( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 if ( sc() > 0.0 ) sc() = log( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               if ( sc() > 0.0 ) sc() = log( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) 
                                 if ( vc(i) > 0.0 ) vc(i)=log( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) 
                                   if ( vc(i) > 0. ) vc(i)=log( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ )
                                 if ( vc(i) > 0.0 ) vc(i)=log( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   if ( ts(i,j) > 0.0 ) ts(i,j) = log( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     if ( ts(i,j) > 0.0 ) ts(i,j) = log( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   if ( ts(i,j) > 0.0 ) ts(i,j) = log( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 

/**
 
Calculates the decadic logarithm (base 10) of each component of the
CSP basic variable which is associated with the PropertyHandle. 

@section messages Messages 

Since the method cannot take the logarithm of 0 or a negative 
number, variables with such a value are not modified. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Log10()
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               if ( sc() > 0.0 ) sc() = log10( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 if ( sc() > 0.0 ) sc() = log10( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               if ( sc() > 0.0 ) sc() = log10( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) 
                                 if ( vc(i) > 0.0 ) vc(i)=log10( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ )
                                   if ( vc(i) > 0.0 ) vc(i)=log10( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ )
                                 if ( vc(i) > 0.0 ) vc(i)=log10( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   if ( ts(i,j) > 0.0 ) ts(i,j) = log10( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     if ( ts(i,j) > 0.0 ) ts(i,j) = log10( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   if ( ts(i,j) > 0.0 ) ts(i,j) = log10( ts(i,j) );
                               (*eit)->Store( key, ts );
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                            }
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**

Calculates E raised to the power of each component of the physical
variable associated with the PropertyHandle. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Exp()
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = exp( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = exp( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = exp( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=exp( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=exp( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=exp( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = exp( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = exp( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = exp( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**
 
Raises each component of the CSP basic variable associated with the
PropertyHandle to the power of its first argument. 

@param raised_to The power to which the variable shall be raised.
*/
template<size_t dim>
void  PropertyHandle<dim>::Pow( double64 raised_to )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = pow( sc(), raised_to );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = pow( sc(), raised_to );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = pow( sc(), raised_to );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=pow( vc(i), raised_to );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=pow( vc(i), raised_to );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=pow( vc(i), raised_to );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = pow( ts(i,j), raised_to );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = pow( ts(i,j), raised_to );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = pow( ts(i,j), raised_to );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**
 
Replaces each NAN (not a number) value of the physical variable which is
associated with the PropertyHandle with the supplied floating point 
argument. 

@param with The floating point number with which NAN values shall be replaced.
*/
template<size_t dim>
void  PropertyHandle<dim>::ZapNAN( double64 with )
 {
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               if ( isnan(sc()) ) sc() = with;
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 if ( isnan(sc()) ) sc() = with;
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               if ( isnan(sc()) ) sc() = with;
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) 
                                 if ( isnan(sc()) ) vc(i) = with;
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) 
                                   if ( isnan(vc(i)) ) vc(i) = with;
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) 
                                 if ( isnan(vc(i)) ) vc(i) = with;
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   if ( isnan(ts(i,j)) ) ts(i,j) = with;
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     if ( isnan(ts(i,j)) ) ts(i,j) = with;
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   if ( isnan(ts(i,j)) ) ts(i,j) = with;
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**
 
Returns the sin of the physical variable component value which must be in
radians. A radian (rad) is defined as the ratio of the corresponding
arc length to the radius of the circle. Thus, 1 radian = 57.2958o and 
the angle in rad is (grad * pi) / 180. 

@section messages Messages 

It is tested whether the input variable is in the legitimate range of zero 
to 2 Pi. If not, an error message is returned and no calculation is
performed. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Sin()
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    double64        omin, omax;
    Range( omin, omax );
    if ( omin < 0.0 || omax > (3.1415927*2.0) )
      {
         throw csmp::Exception( ERROR, "PropertyHandle::Sin", "Operand value not in range of 0 to 2 Pi",
                                "Nothing was done.");
         return;
      }

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = sin( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = sin( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = sin( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=sin( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=sin( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=sin( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = sin( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = sin( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = sin( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**
 
Returns the cos of the physical variable component value which must be in
radians. A radian (rad) is defined as the ratio of the corresponding
arc length to the radius of the circle. Thus, 1 radian = 57.2958o and 
the angle in rad is (grad * pi) / 180. 

@section messages Messages 

It is tested whether the input variable is in the legitimate range of zero 
to 2 Pi. If not, an error message is returned and no calculation is
performed. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Cos()
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    double64                             omin, omax;
    Range( omin, omax );
    if ( omin < 0.0 || omax > (3.1415927*2.0) )
      {
         throw csmp::Exception( ERROR, "PropertyHandle::Cos", "Operand value not in range of 0 to 2 Pi",
                                "Nothing was done.");
         return;
      }

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = cos( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = cos( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = cos( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=cos( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=cos( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=cos( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = cos( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = cos( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = cos( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**
 
Returns the tan of the physical variable component value which must be in
radians. A radian (rad) is defined as the ratio of the corresponding
arc length to the radius of the circle. Thus, 1 radian = 57.2958o and 
the angle in rad is (grad * pi) / 180. 

@section messages Messages 

It is tested whether the input variable is in the legitimate range of zero 
to 2 Pi. If not, an error message is returned and no calculation is
performed. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Tan()
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    double64                             omin, omax;
    Range( omin, omax );
    if ( omin < 0.0 || omax > (3.1415927*2.0) )
      {
         throw csmp::Exception( ERROR, "PropertyHandle::Tan", "Operand value not in range of 0 to 2 Pi",
                                "Nothing was done.");
         return;
      }

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = tan( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = tan( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = tan( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=tan( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=tan( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=tan( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = tan( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = tan( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = tan( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**
 
Returns the arc cos of the physical variable component value. This value
must be between -1 and 1. 

@section messages Messages 

If the physical variable component value is not within a range between
-1 and 1, an error will be reported.
*/
template<size_t dim>
void  PropertyHandle<dim>::Acos()
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    double64                             omin, omax;
    Range( omin, omax );
    if ( omin < -1.0 || omax > 1.0 )
      {
         throw csmp::Exception( ERROR, "PropertyHandle::Acos", "Operand value not in range of -1 to 1",
                                "Nothing was done.");
         return;
      }

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = acos( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = acos( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = acos( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=acos( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=acos( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=acos( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = acos( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = acos( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = acos( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 



/**
 
Returns the arc sin of the physical variable component value. This value
must be between -1 and 1. 

@section messages Messages 

If the physical variable component value is not within a range between
-1 and 1, an error will be reported. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Asin()
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    typename vector<Node<dim>*>::iterator             nit;
    typename vector<Element<dim>*>::iterator          eit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    double64                             omin, omax;
    Range( omin, omax );
    if ( omin < -1.0 || omax > 1.0 )
      {
         throw csmp::Exception( ERROR, "PropertyHandle::Asin", "Operand value not in range of -1 to 1",
                                "Nothing was done.");
         return;
      }

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = asin( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = asin( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = asin( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=asin( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=asin( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=asin( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = asin( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = asin( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = asin( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end ln
 
 




/**
 
Returns the arc tan of the physical variable component value. This value
must be between -1 and 1. 

@section messages Messages 

If the physical variable component value is not within a range between
-1 and 1, an error will be reported. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Atan()
 {
    csmp::Index  key = super_group.Database().StorageKey( var_name.c_str() );
    typename vector<Element<dim>*>::iterator eit;
    typename vector<Node<dim>*>::iterator    nit;
    ScalarVariable                          sc;
    VectorVariable<dim>                      vc;
    TensorVariable<dim>                      ts;
    size_t i, j;
    double64                             omin, omax;
    Range( omin, omax );
    if ( omin < -1.0 || omax > 1.0 )
      {
         throw csmp::Exception( ERROR, "PropertyHandle::Atan", "Operand value not in range of -1 to 1",
                                "Nothing was done.");
         return;
      }

    switch( key.type )
      {
         case SCALAR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                          if ( (*nit)->Status( key ) == flag_output )
                            {
                               (*nit)->Read( key, sc );
                               sc() = atan( sc() );
                               (*nit)->Store( key, sc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                            if ( (*eit)->Status( i, key ) == flag_output ) 
                              {
                                 (*eit)->Read( i, key,  sc );
                                 sc() = atan( sc() );
                                 (*eit)->Store( i, key, sc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          if ( (*eit)->Status( key ) == flag_output ) 
                            {
                               (*eit)->Read( key, sc );
                               sc() = atan( sc() );
                               (*eit)->Store( key, sc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
         case VECTOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=atan( vc(i) );
                               (*nit)->Store( key, vc );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, vc );
                                 for ( i=0; i<dim; i++ ) vc(i)=atan( vc(i) );
                                 (*eit)->Store( i, key, vc );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, vc );
                               for ( i=0; i<dim; i++ ) vc(i)=atan( vc(i) );
                               (*eit)->Store( key, vc );
                           }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
           break;
         case TENSOR:
              switch( key.place )
                {
                   case NODE:
                        for ( nit=group.NodesBegin(); nit!=group.NodesEnd(); nit++ )
                           {
                               (*nit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = atan( ts(i,j) );
                               (*nit)->Store( key, ts );
                            }
                     break;
                   case ELEMENT_INTEGRATION_POINT:
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                          for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
                              {
                                 (*eit)->Read( i, key, ts );
                                 for ( i=0; i<dim; i++ )
                                   for ( j=0; j<dim; j++ )
                                     ts(i,j) = atan( ts(i,j) );
                                 (*eit)->Store( i, key, ts );
                              }
                     break;
                   case ELEMENT:        
                        for ( eit=group.ElementsBegin(); eit!=group.ElementsEnd(); eit++ )
                            {
                               (*eit)->Read( key, ts );
                               for ( i=0; i<dim; i++ )
                                 for ( j=0; j<dim; j++ )
                                   ts(i,j) = atan( ts(i,j) );
                               (*eit)->Store( key, ts );
                            }
                      break;
                   default:
                       throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                            "REGION and MODEL properties are not handled by this class yet."  );
                }
            break;
          default:
            throw csmp::Exception( ERROR, "PropertyHandle<dim>::(method)", 
                           "property type could not be resolved (valid options are: SCALAR, VECTOR or TENSOR)"  );
       } 
    IsWithinRange();
 } // end Atan
 
 



/**
 
Performs a range check for the physical variable which is associated
with the PropertyHandle. 

Range() will return the minimum and the maximum value of the physical 
variable which is associated with the PropertyHandle into its first and its
second argument. 

@section implementation Implementation 

Range() calls the MinMaxOf() interface of the Model. 
*/
template<size_t dim>
void  PropertyHandle<dim>::Range( double64& omin, double64& omax ) const
 {
    super_group.MinMaxOf( var_name.c_str(), omin, omax );
 }





/**
 
Tests whether the physical variable which is associated with the
PropertyHandle has a value which complies with the legitimate range defined
in the PropertyDatabase object. 

If the variable was defined at runtime, it is tested against a default 
range from -1.0e+30 to 1.0e+30. 

@return IsWithinRange() returns the boolean variable 'true' if the variable
values lie within the expected range.

@section implementation Implementation 

sWithinRange() compares the actual value range of the variable with 
that specified in the PropertyDatabase. 

@section application Application 

If so desired by the user, IsWithinRange() will be applied after each 
operation on a PropertyHandle. 
*/
template<size_t dim>
bool  PropertyHandle<dim>::IsWithinRange() const
 {
    if ( key_.type == TENSOR ) {
         ErrorHandler&  csmp_error(ErrorHandler::Instance());
         csmp_error.notice( WARNING, "PropertyHandle<dim>::IsWithinRange:",
                           "cannot return range of arbitrary tensors for which Eigenvalues cannot be found.");
         return true;
      }
    double64  omin, omax, pmin, pmax;
    group.MinMaxOf( var_name.c_str(), omin, omax );
    super_group.Database().RangeOf(  var_name.c_str(), pmin, pmax );
    if ( omin >= pmin && omax <= pmax ) return true;

    return false;
 }



/**
 
The standard CSP object interface Out() will print the state of the
PropertyHandle to stdout. This output comprises the name of the physical
variable which is associated with the PropertyHandle, the output flag,
and the values of the distributed physical variable. 

@section implementation Implementation 

The method Out() uses the Model method OutputVariableToScreen() to
print the variable values. 
*/ 
template<size_t dim>
void  PropertyHandle<dim>::Out() const 
 {
    cout <<"\nPropertyHandle('"<< group_name <<"')::Out: "<< endl;
    cout <<"Physical variable: "<< var_name << endl;
    cout <<"flag_output:         "<< parseStatus(flag_output) << endl;
    group.OutputVariableToScreen( var_name.c_str() );

 } // end Out




/**
 
The CSP object interface Out() will print the values of the PropertyHandle
associated variable to a textfile. 

@param text_file_name The name of the textfile is supplied as first argument to Out(). 

@section implementation Implementation 

Out() relies on the Model method OutputDataAsTextColumns() in the 
creation of the output textfile.  

@section application Application 

To save a PropertyHandle temporary variable to file. 

@section messages Messages 

Refer to the documentation of the Model to learn about the 
footprint of OutputDataAsTextColumns(). 
*/
template<size_t dim>
void  PropertyHandle<dim>::Out( const char* text_file_name ) const
 {
    ofstream  ofs(text_file_name);
 
    ofs <<"\nPropertyHandle::Out: "<< endl;
    ofs <<"Physical variable: "<< var_name << endl;
    ofs <<"flag_output:         "<< parseStatus(flag_output) << endl;

 } // end Out



template class PropertyHandle<1U>;
template class PropertyHandle<2U>;
template class PropertyHandle<3U>;


} // end namespace csp








