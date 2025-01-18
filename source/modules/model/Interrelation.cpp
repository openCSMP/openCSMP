#include "Interrelation.h"
#include "PropertyDatabase.h"
#include "Node.h"
#include "Face.h"
#include "InterFace.h"
#include "Element.h"
#include "Boundary.h"
#include "SplitBoundary.h"
#include "Region.h"
#include "Model.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**

Builds Interrelation instance and stores a reference to the variable 
database. The private variables of the Interrelation are set as
follows. The application target is the Model, the application level is
the element, and the name of the interrelation and the result property 
are 'unspecified'.  

@param p A constant reference to the PropertyDatabase inside of the Model
object. Use the Model method Database() to obtain
this reference. 

@section application Application

Is called inside each interrelation subclass. 
*/
template<uint32_t dim>
Interrelation<dim>::Interrelation( const PropertyDatabase<dim>& p ) 
  : p_ref(p),
    name_("unspecified"),
    result_property_("unspecified"),
    application_level_(ELEMENT)
 {
 }
  
  
  
  
template<uint32_t dim>
Interrelation<dim>::~Interrelation()
 {
 } 
 



 








/**

Sets the name of an Interrelation subclass such that it can be output
in an error message of the base class. 

@param s The name of the Interrelation subclass is input as first method argument.
While one can input arbitrary names it will be easier to detect errors 
which may arise, if the name is exactly identical to the name of the 
Interrelation subclass. 

@section application Application

Assign a name to an interrelation subclass inside of its constructor. 
*/
template<uint32_t dim>
void  Interrelation<dim>::Name( const char* s )             
{ name_ = s; }



 


/**

Associates a reference to an Operand object with a globally known 
physical variable registered in the PropertyDatabase. 

@param var The name of the property which shall be associated with the Operand.

@return A reference to an Operand which is stored inside an Operand map inside
of the Interrelation object. 

@section implementation Implementation

Looks for existing property in the database and assigns it to an
Operand by storing the two as a pair in an Operand map.  

@section application Application

Call inside of the constructor of an Interrelation subclass to initialize
the corresponding Operand reference.  

@section messages Messages

Errors will be reported if the target property is unknown to the 
PropertyDatabase or if the Operand already exists in the Operand map. 
*/
template<uint32_t dim>
Operand<dim>&  Interrelation<dim>::GlobalProperty( const char* var )
 {
    if ( p_ref.IsDefined(var) == false )
      throw csmp::Exception( ERROR, "Interrelation<dim>::GlobalProperty", 
                             "Operand does not exist in variable database");
    
    Operand<dim> op( var, p_ref );

    operand_list_[ op.Name() ] = op;
    
    typename map<string,Operand<dim> >::iterator it(operand_list_.find(op.Name())); 
     
    if ( it == operand_list_.end() )
      throw csmp::Exception( ERROR, "Interrelation<dim>::GlobalProperty", 
                             "Operand assignment to operand list failed");
    
    return (*it).second;
       
 } // end GlobalProperty
 


 
 
                                                                       
/**

Specifies which property will be written back to the Model once the 
calculation has been performed.  

@param var the name of the physical variable which shall be modified by the
calculation. 

@section implementation Implementation

Method defines which Operand is associated with the result property of 
the calculation and sets the application level accordingly, ie. if the 
result property is an element variable the target level will be the 
element and so forth. 

@section application Application

Call inside of the constructor of the Interrelation subclass to define
which property shall be modified by the Interrelation. 

@section messages Messages

Errors will be reported if the target property is unknown to the 
PropertyDatabase or if the Operand already exists in the Operand map. 
*/
template<uint32_t dim>
void  Interrelation<dim>::ResultProperty( const char* var )
 {
    result_ = operand_list_.find( string(var) );
    
    if ( result_ == operand_list_.end() ) 
         throw csmp::Exception( ERROR, "Interrelation<dim>::DefineResultProperty", 
                                "Result property was not found in Operand list");
    result_property_   = var;                            
    application_level_ = (*result_).second.Placement();                            
 }






/**

OutputCondition() determines the flag of the locations in the mesh where 
the global physical variable which is modified by the Interrelation is 
allowed to be stored back into the Model object. 

@param var a reference to the Operand which
is associated with the output variable
@param output the output flag, controls at which locations the calculated variable value is output
to the Model and where the flag differs from the local variable 
flag and output is therefore prohibited. 

@section implementation Implementation

Calls the Operand class interface OutputCondition(). 

@section application Application

OutputCondition() must be called inside the constructor of the 
Interrelation subclass. 
*/
template<uint32_t dim>
void Interrelation<dim>::OutputCondition( Operand<dim>& var, VARIABLE_FLAG output )
 {
    var.OutputCondition( output );  
 }




/**
     Method for the application of Interrelations to the entire model 
     including subregions and boundaries. 
     
     If the result variable of an interrelation is placed on the REGION or 
     the BOUNDARY, other variables that are involved in the calculation 
     are volume integrated over each region, then the result is normalized
     by the region volume/surface/length.
*/
template<uint32_t dim>
void  Interrelation<dim>::Apply( Model<dim>& sg )
 {
     if ( application_level_ == REGION ) 
       {
         // unique regions first
         for ( typename map<string,csmp::Region<dim> >::iterator
               git=sg.UniqueRegionsBegin(); git!=sg.UniqueRegionsEnd(); git++ )
           {
              // 1. Get the Operands
              // -------------------
              for ( typename map<string,Operand<dim> >::iterator
                    oiter=operand_list_.begin(); oiter!=operand_list_.end(); oiter++ )
                {
                   csmp::Index prop_key = (*oiter).second.Key();
                   
                   if ( prop_key.place == BOUNDARY )
                    throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                    "BOUNDARY variables cannot be involved in REGION calculations ");
                   
                   // REGION properties
                   if ( prop_key.place == REGION ) {
                       switch( prop_key.type ) {
                            case SCALAR:
                                 (*oiter).second = (*git).second.Read( prop_key ); 
                              break;
                            case VECTOR: {
                                   VectorVariable<dim>  vc;
                                   (*git).second.Read( prop_key, vc ); 
                                   (*oiter).second = vc;
                                 }
                              break;
                            case TENSOR: {
                                   TensorVariable<dim>  ts;
                                   (*git).second.Read( prop_key, ts ); 
                                   (*oiter).second = ts;
                                 }
                              break;
                            default: throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                                                   "variable type could not be resolved");
                         }
                     }
                   // OTHER properties are volume averaged
                   else {
                        if ( prop_key.type != SCALAR )
                          throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                          "only SCALAR non-REGION variables can be involved in REGION calculations ");
                       
                        (*oiter).second = (*git).second.VolumeIntegral( sg.Database().Name(prop_key) ,false) / (*git).second.Volume();
                     }
                   
                 } // end looping through operands

              // 2. Doing the calculation of the result
              // --------------------------------------
              Calculate();
                        
              // 3. Mapping the result back to the current element or the new variable container
              // -------------------------------------------------------------------------------
              if ( ResultWithinRange() == false ) {
                  (*result_).second.Out();
                   string message("calculated result out of range; nothing is done; this affects the result property '");
                   message += result_property_.c_str();
                   message +="'";
                   throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(), message.c_str()  );
                }
              else
              switch( (*result_).second.Type() ) {
                   case SCALAR: {
                          ScalarVariable  sc;
                          (*result_).second.AssignTo( sc );
                          (*git).second.Store( (*result_).second.Key(), sc );
                        }
                     break;
                   case VECTOR: {
                          VectorVariable<dim>  vc;
                          (*result_).second.AssignTo( vc );
                          (*git).second.Store( (*result_).second.Key(), vc );
                        }
                     break;
                   case TENSOR: {
                          TensorVariable<dim>  ts;
                          (*result_).second.AssignTo( ts );
                          (*git).second.Store( (*result_).second.Key(), ts );
                        }
                      break;
                    default: 
                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                           "variable type could not be resolved");
                 }
            } // end regions loop
            
         // non-unique regions second
         for ( typename map<string,csmp::Region<dim> >::iterator
               git=sg.RegionsBegin(); git!=sg.RegionsEnd(); git++ )
           {
              // 1. Get the Operands
              // -------------------
              for ( typename map<string,Operand<dim> >::iterator
                    oiter=operand_list_.begin(); oiter!=operand_list_.end(); oiter++ )
                {
                   csmp::Index prop_key = (*oiter).second.Key();
                   
                   if ( prop_key.place == BOUNDARY )
                    throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                    "BOUNDARY variables cannot be involved in REGION calculations ");
                   
                   // REGION properties
                   if ( prop_key.place == REGION ) {
                       switch( prop_key.type ) {
                            case SCALAR:
                                 (*oiter).second = (*git).second.Read( prop_key ); 
                              break;
                            case VECTOR: {
                                   VectorVariable<dim>  vc;
                                   (*git).second.Read( prop_key, vc ); 
                                   (*oiter).second = vc;
                                 }
                              break;
                            case TENSOR: {
                                   TensorVariable<dim>  ts;
                                   (*git).second.Read( prop_key, ts ); 
                                   (*oiter).second = ts;
                                 }
                              break;
                            default: throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                                                   "variable type could not be resolved");
                         }
                     }
                   // OTHER properties are volume averaged
                   else {
                        if ( prop_key.type != SCALAR )
                          throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                          "only SCALAR non-REGION variables can be involved in REGION calculations ");
                       
                        (*oiter).second = (*git).second.VolumeIntegral( sg.Database().Name(prop_key) ,false) / (*git).second.Volume();
                     }
                   
                 } // end looping through operands

              // 2. Doing the calculation of the result
              // --------------------------------------
              Calculate();
                        
              // 3. Mapping the result back to the current element or the new variable container
              // -------------------------------------------------------------------------------
              if ( ResultWithinRange() == false ) {
                  (*result_).second.Out();
                   string message("calculated result out of range; nothing is done; this affects the result property '");
                   message += result_property_.c_str();
                   message +="'";
                   throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(), message.c_str()  );
                }
              else
              switch( (*result_).second.Type() ) {
                   case SCALAR: {
                          ScalarVariable  sc;
                          (*result_).second.AssignTo( sc );
                          (*git).second.Store( (*result_).second.Key(), sc );
                        }
                     break;
                   case VECTOR: {
                          VectorVariable<dim>  vc;
                          (*result_).second.AssignTo( vc );
                          (*git).second.Store( (*result_).second.Key(), vc );
                        }
                     break;
                   case TENSOR: {
                          TensorVariable<dim>  ts;
                          (*result_).second.AssignTo( ts );
                          (*git).second.Store( (*result_).second.Key(), ts );
                        }
                      break;
                    default: 
                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                           "variable type could not be resolved");
                 }
            } // end non-unique regions loop
            
          return;
          
       } // end for region
     
     if (application_level_ == BOUNDARY ) 
       {
         for ( typename map<std::string,csmp::Boundary<dim> >::iterator
               git=sg.BoundariesBegin(); git!=sg.BoundariesEnd(); git++ )
           {
              // 1. Get the Operands
              // -------------------
              for ( typename map<string,Operand<dim> >::iterator
                    oiter=operand_list_.begin(); oiter!=operand_list_.end(); oiter++ )
                {
                   csmp::Index prop_key = (*oiter).second.Key();
                   
                   if ( prop_key.place == REGION )
                    throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                    "REGION variables cannot be involved in BOUNDARY calculations ");
                   
                   // REGION properties
                   if ( prop_key.place == REGION ) {
                       switch( prop_key.type ) {
                            case SCALAR:
                                 (*oiter).second = (*git).second.Read( prop_key ); 
                              break;
                            case VECTOR: {
                                   VectorVariable<dim>  vc;
                                   (*git).second.Read( prop_key, vc ); 
                                   (*oiter).second = vc;
                                 }
                              break;
                            case TENSOR: {
                                   TensorVariable<dim>  ts;
                                   (*git).second.Read( prop_key, ts ); 
                                   (*oiter).second = ts;
                                 }
                              break;
                            default: throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                                                   "variable type could not be resolved");
                         }
                     }
                   // OTHER properties are volume averaged
                   else {
                        if ( prop_key.type != SCALAR )
                          throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                          "only SCALAR non-BOUNDARY variables can be involved in BOUNDARY calculations ");
                       
                        (*oiter).second = (*git).second.SurfaceIntegral( sg.Database(), sg.Database().Name(prop_key) ) / (*git).second.Area(); 
                     }
                   
                 } // end looping through operands

              // 2. Doing the calculation of the result
              // --------------------------------------
              Calculate();
                        
              // 3. Mapping the result back to the current element or the new variable container
              // -------------------------------------------------------------------------------
              if ( ResultWithinRange() == false ) {
                  (*result_).second.Out();
                   string message("calculated result out of range; nothing is done; this affects the result property '");
                   message += result_property_.c_str();
                   message +="'";
                   throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(), message.c_str()  );
                }
              else
              switch( (*result_).second.Type() ) {
                   case SCALAR: {
                          ScalarVariable  sc;
                          (*result_).second.AssignTo( sc );
                          (*git).second.Store( (*result_).second.Key(), sc );
                        }
                     break;
                   case VECTOR: {
                          VectorVariable<dim>  vc;
                          (*result_).second.AssignTo( vc );
                          (*git).second.Store( (*result_).second.Key(), vc );
                        }
                     break;
                   case TENSOR: {
                          TensorVariable<dim>  ts;
                          (*result_).second.AssignTo( ts );
                          (*git).second.Store( (*result_).second.Key(), ts );
                        }
                      break;
                    default: 
                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                           "variable type could not be resolved");
                 }
            } // end boundaries loop

         for ( typename map<std::string,csmp::SplitBoundary<dim> >::iterator
               git=sg.SplitBoundariesBegin(); git!=sg.SplitBoundariesEnd(); git++ )
           {
              // 1. Get the Operands
              // -------------------
              for ( typename map<string,Operand<dim> >::iterator
                    oiter=operand_list_.begin(); oiter!=operand_list_.end(); oiter++ )
                {
                   csmp::Index prop_key = (*oiter).second.Key();
                   
                   if ( prop_key.place == REGION )
                    throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                    "REGION variables cannot be involved in BOUNDARY calculations ");
                   
                   // REGION properties
                   if ( prop_key.place == REGION ) {
                       switch( prop_key.type ) {
                            case SCALAR:
                                 (*oiter).second = (*git).second.Read( prop_key ); 
                              break;
                            case VECTOR: {
                                   VectorVariable<dim>  vc;
                                   (*git).second.Read( prop_key, vc ); 
                                   (*oiter).second = vc;
                                 }
                              break;
                            case TENSOR: {
                                   TensorVariable<dim>  ts;
                                   (*git).second.Read( prop_key, ts ); 
                                   (*oiter).second = ts;
                                 }
                              break;
                            default: throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                                                   "variable type could not be resolved");
                         }
                     }
                   // OTHER properties are volume averaged
                   else {
                        if ( prop_key.type != SCALAR )
                          throw Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                          "only SCALAR non-BOUNDARY variables can be involved in BOUNDARY calculations ");
                       
                        (*oiter).second = (*git).second.SurfaceIntegral( sg.Database(), sg.Database().Name(prop_key) ) / (*git).second.Area(); 
                     }
                   
                 } // end looping through operands

              // 2. Doing the calculation of the result
              // --------------------------------------
              Calculate();
                        
              // 3. Mapping the result back to the current element or the new variable container
              // -------------------------------------------------------------------------------
              if ( ResultWithinRange() == false ) {
                  (*result_).second.Out();
                   string message("calculated result out of range; nothing is done; this affects the result property '");
                   message += result_property_.c_str();
                   message +="'";
                   throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(), message.c_str()  );
                }
              else
              switch( (*result_).second.Type() ) {
                   case SCALAR: {
                          ScalarVariable  sc;
                          (*result_).second.AssignTo( sc );
                          (*git).second.Store( (*result_).second.Key(), sc );
                        }
                     break;
                   case VECTOR: {
                          VectorVariable<dim>  vc;
                          (*result_).second.AssignTo( vc );
                          (*git).second.Store( (*result_).second.Key(), vc );
                        }
                     break;
                   case TENSOR: {
                          TensorVariable<dim>  ts;
                          (*result_).second.AssignTo( ts );
                          (*git).second.Store( (*result_).second.Key(), ts );
                        }
                      break;
                    default: 
                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply(Model)", name_.c_str(),
                                          "variable type could not be resolved");
                 }
            } // end split-boundaries loop
            
         return;
         
       } // end BOUNDARY
 
    // if the result variable is Node, IntegrationPoint, or Element
    Region<dim>& gref(sg.Region("Model"));
    Apply( gref );
    
 } // end Apply (model)




 
/**
 
Apply() loops over the elements, nodes, or constraint points in a 
subregion of the Model mesh which has been identified as a group. It
collects the values of the input variables of the interrelation, calls the 
method Calculate(), checks the computed results, and maps these back
into the Model, ie. the MemoryManager. The output variable is
stored back into the Model, if and only if it has a flag equivalent 
to the output flag of the interrelation and under the proviso that its
value is inside the legitimate range. 

@param gref a pointer to the target group to which the interrelation shall be
applied to. 

@section implementation Implementation

A few restrictions apply with regard to the placement of the output 
variable and the placement of the variables which can be employed in an 
interrelation: 

If the output variable is a node or constraint point variable, the 
interrelation must only involve node or constraint point variables,
respectively. To explain why, let's assume that one of the input variables
was an element variable but the output variable would be a node variable.
Since the node simultaneously belongs to several elements whose properties
may vary greatly, it is undefined from which element the interrelation
should attain the element property input value. For this reason this type 
of interrelation is not supported. 

If the output variable is an element variable but the calculation employs
node and/or constraint point variables, Apply() uses the Average() interface
of the Element to obtain a mean value of the node property. In many cases,
this is legitimate but there are cases where this is not so. For instance,
if you have a depth-dependent nodal property and two adjacent triangular
elements that break up a square, one of the elements will have two nodes
at the top of the square and the other one has two at the bottom. The
computed mean value will be different in these two elements which may
be reflected in the result property. 

Importantly, interrelations are but one way of computing local 
interdependencies of physical variables. Please refer to the documentation
of the Visitor base class to see how one can achieve more 
flexibility. Also, in a CSMP mesh each node or constraint point knows 
its parent elements which are stored in counterclockwise (righthandrule)
fashion in a corresponding vector (see Node and IntegrationPoint 
documetation). This knowledge may be employed to compute node properties 
from element properties. 

Finally, since each element knows its neighbors, mesh sweeps may be 
instrumentalized via visitor objects if so desired. 

@section application Application

Is called by the Model when the interrelation is passed to it using 
the Apply() interface. 

@section messages Messages

Apply() will report an error if the target group cannot be accessed (it may
not have been defined in this case) or if the output variable of the 
interrelation is a Node or ConstrainPoint property, but the calculation involves 
Element properties. If the result property value of an interrelation
is outside the legitimate range another error is reported together
with the name of the property which is concerned.  
 */
template<uint32_t dim>
template<template<uint32_t> class CELL>
void  Interrelation<dim>::Apply( ModelSubDomain<dim,CELL>& gref )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
   
    if ( application_level_ == REGION or application_level_ == BOUNDARY )
      throw Exception( ERROR, "Interrelation<dim>::Apply(ModelSubDomain)", name_.c_str(),
                      "attempt to apply this on the application level REGION or BOUNDARY" );
 
     auto nst1 = gref.NodesBegin(), nst2 = gref.NodesEnd();
     auto est1 = gref.CellsBegin(), est2 = gref.CellsEnd();

     ScalarVariable       sc;
     VectorVariable<dim>  vc;
     TensorVariable<dim>  ts;
 
     // ----------------------
     // RESULT = NODE VARIABLE
     // ----------------------
     // if the interrelationship changes a Node Property, it can only contain Node
     // variables in the associated calculation.
     // ----------------------------------------- 
     if ( application_level_ == NODE ) 
       {
         while ( nst1 != nst2 ) {
              // 1. Get the Operands
              // -------------------
              for ( typename map<string,Operand<dim> >::iterator
                    oiter=operand_list_.begin(); oiter!=operand_list_.end(); oiter++ )
                {
                   csmp::Index prop_key = (*oiter).second.Key();
                   
                   if ( prop_key.place != NODE ) {
                        throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(), 
                                              "calculations with nodal result variables can only involve node variables");
                        return;
                     }
                   switch( prop_key.type ) {
                        case SCALAR:
                             (*oiter).second = (*nst1)->Read( prop_key ); 
                          break;
                        case VECTOR:
                             (*nst1)->Read( prop_key, vc ); 
                             (*oiter).second = vc;
                          break;
                        case TENSOR:
                             (*nst1)->Read( prop_key, ts ); 
                             (*oiter).second = ts;
                          break;
                        default: 
                        throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                               "variable type could not be resolved");
                     }
              } // end looping through operands
              
            // 2. Doing the calculation of the result
            // --------------------------------------
            Calculate();
            
            // 3. Mapping the result back to the current node or the new variable container
            // ----------------------------------------------------------------------------
            // 3.1 range checking
            if ( ResultWithinRange() == false ) {
                 csmp_error.Note( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                           "calculated result out of range; nothing is done.");
                 (*result_).second.Out();
                 csmp_error.Note( INFO, "This affects the result property ", result_property_.c_str() );
              }
            // 3.2 if the calculated variable exists in database the result is stored in Model  
            else {
                 switch( (*result_).second.Type() ) {
                     case SCALAR: 
                          (*result_).second.AssignTo( sc );
                          (*nst1)->Store( (*result_).second.Key(), sc );
                       break;
                     case VECTOR:
                          (*result_).second.AssignTo( vc );
                          (*nst1)->Store( (*result_).second.Key(), vc );
                       break;
                     case TENSOR:
                          (*result_).second.AssignTo( ts );
                          (*nst1)->Store( (*result_).second.Key(), ts );
                      break;
                    default: 
                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                          "variable type could not be resolved.");
                  }
              }
             nst1++;  

           } // end of traversing NODE array
         return;
      }


     // ----------------------------------
     // RESULT = ELEMENT_INTEGRATION_POINT VARIABLE
     // ----------------------------------
     // if the interrelationship changes an segment property, averages of the node 
     // properties are used to compute this value
     // ----------------------------------------- 
     if ( application_level_ == ELEMENT_INTEGRATION_POINT ) 
       {
          for ( auto eit=gref.CellsBegin(); eit!=gref.CellsEnd(); eit++ )
            // looping over the constraint points inside the current element    
            for ( auto i{0U}; i<(*eit)->IntegrationPoints(); i++ ) 
              {
                // 1. Get the Operands
                // -------------------
                for ( typename map<string,Operand<dim> >::iterator
                      oiter=operand_list_.begin(); oiter!=operand_list_.end(); oiter++ )
                  {
                     csmp::Index prop_key = (*oiter).second.Key();
                     
                     switch( prop_key.place )
                       {
                           case NODE: // properties are interpolated to integration point
                               switch( prop_key.type ) {
                                    case SCALAR:
                                         (*oiter).second = (*eit)->PropertyValueAtIntegrationPoint( prop_key, i );
                                      break;
                                    case VECTOR:
                                         (*eit)->PropertyValueAtIntegrationPoint( prop_key, i, vc );
                                         (*oiter).second = vc;
                                      break;
                                    case TENSOR:
                                         (*eit)->PropertyValueAtIntegrationPoint( prop_key, i, ts );
                                         (*oiter).second = ts;
                                      break;
                                    default: 
                                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                                          "node variable type could not be resolved.");
                                  }
                            break;
                          case ELEMENT_INTEGRATION_POINT:
                               switch( prop_key.type ) {
                                    case SCALAR:
                                         (*oiter).second = (*eit)->Read( i, prop_key ); 
                                      break;
                                    case VECTOR:
                                         (*eit)->Read( i, prop_key, vc ); 
                                         (*oiter).second = vc;
                                      break;
                                    case TENSOR:
                                         (*eit)->Read( i, prop_key, ts ); 
                                         (*oiter).second = ts;
                                      break;
                                    default: 
                                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                                           "integration point variable; type could not be resolved.");
                                  }
                            break;
                          case ELEMENT:
                               switch( prop_key.type ) {
                                    case SCALAR:
                                         (*oiter).second = (*eit)->Read( prop_key ); 
                                      break;
                                    case VECTOR:
                                         (*eit)->Read( prop_key, vc ); 
                                         (*oiter).second = vc;
                                      break;
                                    case TENSOR:
                                         (*eit)->Read( prop_key, ts ); 
                                         (*oiter).second = ts;
                                      break;
                                    default: 
                                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                                           "element variable type could not be resolved.");
                                  }
                            break;
                          case REGION:
                               switch( prop_key.type ) {
                                    case SCALAR:
                                         (*oiter).second = dynamic_cast<Region<dim>&>(gref).Read( prop_key ); 
                                      break;
                                    case VECTOR:
                                         dynamic_cast<Region<dim>&>(gref).Read( prop_key, vc ); 
                                         (*oiter).second = vc;
                                      break;
                                    case TENSOR:
                                         dynamic_cast<Region<dim>&>(gref).Read( prop_key, ts ); 
                                         (*oiter).second = ts;
                                      break;
                                    default: 
                                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                                           "region variable type could not be resolved.");
                                  }
                            break;
                           case BOUNDARY:
                               switch( prop_key.type ) {
                                    case SCALAR:
                                         (*oiter).second = dynamic_cast<Boundary<dim>&>(gref).Read( prop_key ); 
                                      break;
                                    case VECTOR:
                                         dynamic_cast<Boundary<dim>&>(gref).Read( prop_key, vc ); 
                                         (*oiter).second = vc;
                                      break;
                                    case TENSOR:
                                         dynamic_cast<Boundary<dim>&>(gref).Read( prop_key, ts ); 
                                         (*oiter).second = ts;
                                      break;
                                    default: 
                                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                                           "boundary variable type could not be resolved.");
                                  }
                            break;
                           case SPLIT_BOUNDARY:
                               switch( prop_key.type ) {
                                    case SCALAR:
                                         (*oiter).second = dynamic_cast<SplitBoundary<dim>&>(gref).Read( prop_key ); 
                                      break;
                                    case VECTOR:
                                         dynamic_cast<SplitBoundary<dim>&>(gref).Read( prop_key, vc ); 
                                         (*oiter).second = vc;
                                      break;
                                    case TENSOR:
                                         dynamic_cast<SplitBoundary<dim>&>(gref).Read( prop_key, ts ); 
                                         (*oiter).second = ts;
                                      break;
                                    default: 
                                    throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                                           "boundary variable type could not be resolved.");
                                  }
                            break;
                           default: csmp_error.Note( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                  "so fR, calculations with 'IntegrationPoint' target cannot involve Face variables.");
                                  return;
                     } // end switch

                  } // end looping through operands
                  
                // 2. Doing the calculation of the result
                // --------------------------------------
                Calculate();
                
                // 3. Mapping the result back to the current segment or the new variable container
                // -------------------------------------------------------------------------------
                // 3.1 range checking
                if ( ResultWithinRange() == false )
                  {
                     csmp_error.Note( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                            "calculated result out of range; nothing is done");
                     (*result_).second.Out();
                     csmp_error.Note( INFO, "This affects the result property ", result_property_.c_str() );
                  }
                // 3.2 if the calculated variable exists in database the result is stored in Model  
                else 
                  {
                     switch( (*result_).second.Type() ) {
                         case SCALAR: 
                              (*result_).second.AssignTo( sc );
                              (*eit)->Store( i, (*result_).second.Key(), sc );
                           break;
                         case VECTOR:
                              (*result_).second.AssignTo( vc );
                              (*eit)->Store( i, (*result_).second.Key(), vc );
                           break;
                         case TENSOR:
                              (*result_).second.AssignTo( ts );
                              (*eit)->Store( i, (*result_).second.Key(), ts );
                            break;
                          default: 
                          throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                                 "variable type could not be resolved");
                      }
                 
                 } 
             
             } // end integration point loop     

         return;

      } // end IntegrationPoints
           
           
     // -------------------------
     // RESULT = ELEMENT VARIABLE
     // -------------------------
     // if the interrelationship changes an element property, averages of the node segment etc.
     // properties are used to compute this value
     // ----------------------------------------- 
     if (application_level_ == ELEMENT ) {
         while ( est1 != est2 )
           {
              // 1. Get the Operands
              // -------------------
              for ( typename map<string,Operand<dim> >::iterator
                    oiter=operand_list_.begin(); oiter!=operand_list_.end(); oiter++ )
                {
                   csmp::Index prop_key = (*oiter).second.Key();
                   
                   switch( prop_key.place )
                     {
                        case NODE:
                             switch( prop_key.type )
                               {
                                  case SCALAR:
                                       (*est1)->PropertyValueAtBaryCenter( prop_key, sc ); 
                                       (*oiter).second = sc;
                                    break;
                                  case VECTOR:
                                       (*est1)->PropertyValueAtBaryCenter( prop_key, vc ); 
                                       (*oiter).second = vc;
                                    break;
                                  case TENSOR:
                                       (*est1)->PropertyValueAtBaryCenter( prop_key, ts ); 
                                       (*oiter).second = ts;
                                    break;
                                  default: 
                                  throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                                         "variable type could not be resolved");
                                }
                          break;
                        case ELEMENT_INTEGRATION_POINT:
                             switch( prop_key.type )
                               {
                                  case SCALAR:
                                       (*est1)->PropertyValueAtBaryCenter( prop_key, sc ); 
                                       (*oiter).second = sc;
                                    break;
                                  case VECTOR:
                                       (*est1)->PropertyValueAtBaryCenter( prop_key, vc ); 
                                       (*oiter).second = vc;
                                    break;
                                  case TENSOR:
                                       (*est1)->PropertyValueAtBaryCenter( prop_key, ts ); 
                                       (*oiter).second = ts;
                                    break;
                                  default: 
                                  throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                                         "variable type could not be resolved");
                                }
                          break;
                        case ELEMENT:
                             switch( prop_key.type )
                               {
                                  case SCALAR:
                                       (*oiter).second = (*est1)->Read( prop_key ); 
                                    break;
                                  case VECTOR:
                                       (*est1)->Read( prop_key, vc ); 
                                       (*oiter).second = vc;
                                    break;
                                  case TENSOR:
                                       (*est1)->Read( prop_key, ts ); 
                                       (*oiter).second = ts;
                                    break;
                                  default: 
                                  throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                                         "variable type could not be resolved");
                                }
                           break;
                         case REGION:
                             switch( prop_key.type )
                               {
                                  case SCALAR:
                                       (*oiter).second = dynamic_cast<Region<dim>&>(gref).Read( prop_key ); 
                                    break;
                                  case VECTOR:
                                       dynamic_cast<Region<dim>&>(gref).Read( prop_key, vc ); 
                                       (*oiter).second = vc;
                                    break;
                                  case TENSOR:
                                       dynamic_cast<Region<dim>&>(gref).Read( prop_key, ts ); 
                                       (*oiter).second = ts;
                                    break;
                                  default: 
                                  throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                                         "variable type could not be resolved");
                              }
                           break;
                        case BOUNDARY:
                             switch( prop_key.type )
                               {
                                  case SCALAR:
                                       (*oiter).second = dynamic_cast<Boundary<dim>&>(gref).Read( prop_key ); 
                                    break;
                                  case VECTOR:
                                       dynamic_cast<Boundary<dim>&>(gref).Read( prop_key, vc ); 
                                       (*oiter).second = vc;
                                    break;
                                  case TENSOR:
                                       dynamic_cast<Boundary<dim>&>(gref).Read( prop_key, ts ); 
                                       (*oiter).second = ts;
                                    break;
                                  default: 
                                  throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                                         "variable type could not be resolved");
                               }
                             break;
                       case SPLIT_BOUNDARY:
                             switch( prop_key.type )
                               {
                                  case SCALAR:
                                       (*oiter).second = dynamic_cast<SplitBoundary<dim>&>(gref).Read( prop_key ); 
                                    break;
                                  case VECTOR:
                                       dynamic_cast<SplitBoundary<dim>&>(gref).Read( prop_key, vc ); 
                                       (*oiter).second = vc;
                                    break;
                                  case TENSOR:
                                       dynamic_cast<SplitBoundary<dim>&>(gref).Read( prop_key, ts ); 
                                       (*oiter).second = ts;
                                    break;
                                  default: 
                                  throw csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                                         "variable type could not be resolved");
                                }
                           break;
                        default: csmp::Exception( ERROR, "Interrelation<dim>::Apply", name_.c_str(), 
                                                 "calculation with 'Element' variable; input variable placement cannot be resolved");
                          return;
                   } // end switch
              } // end looping through operands
              
            // 2. Doing the calculation of the result
            // --------------------------------------
            Calculate();
                      
            // 3. Mapping the result back to the current element or the new variable container
            // -------------------------------------------------------------------------------
            // 3.1 range checking
            if ( ResultWithinRange() == false )
              {
                 csmp_error.Note( ERROR, "Interrelation<dim>::Apply", name_.c_str(),
                                        "calculated result out of range; nothing is done");
                 (*result_).second.Out();
                 csmp_error.Note( INFO, "This affects the result property ", result_property_.c_str() );
              }
            // 3.2 if the calculated variable exists in database the result is stored in Model  
            else 
              {
                 switch( (*result_).second.Type() )
                  {
                     case SCALAR: 
                          (*result_).second.AssignTo( sc );
                          (*est1)->Store( (*result_).second.Key(), sc );
                       break;
                     case VECTOR:
                          (*result_).second.AssignTo( vc );
                          (*est1)->Store( (*result_).second.Key(), vc );
                       break;
                     case TENSOR:
                          (*result_).second.AssignTo( ts );
                          (*est1)->Store( (*result_).second.Key(), ts );
                        break;
                      default: 
                      throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                                             "variable type could not be resolved.");
                  }
              }
           est1++;
         } // end of traversing element array
       return;
      }

      if ( application_level_ == FACE or application_level_ == INTER_FACE ) 
        throw csmp::Exception( ERROR, "Interrelation<dim>::Apply:", name_.c_str(),
                              "calculations with Face or InterFace as result variables are not handled yet.");
 
 } // end Apply (Region level)


template class Interrelation<1U>;
template class Interrelation<2U>;
template class Interrelation<3U>;


template void  Interrelation<1U>::Apply( ModelSubDomain<1U,Element>& );
template void  Interrelation<2U>::Apply( ModelSubDomain<2U,Element>& );
template void  Interrelation<3U>::Apply( ModelSubDomain<3U,Element>& );

template void  Interrelation<1U>::Apply( ModelSubDomain<1U,Face>& );
template void  Interrelation<2U>::Apply( ModelSubDomain<2U,Face>& );
template void  Interrelation<3U>::Apply( ModelSubDomain<3U,Face>& );

template void  Interrelation<1U>::Apply( ModelSubDomain<1U,InterFace>& );
template void  Interrelation<2U>::Apply( ModelSubDomain<2U,InterFace>& );
template void  Interrelation<3U>::Apply( ModelSubDomain<3U,InterFace>& );

} // end namespace csmp


 
  
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
