#include "Visitor.h"
#include "Exception.h"

using namespace std;

// Base class for the visitor scheme.

namespace csmp {

/**
 
Constructor of visitor object. The user can specify the application level and 
application target of the visitor by supplying enum PLACEMENT 
argument. Per default, visitors act on Element objects.  

@param level PLACEMENT that will be used to determine where in the object
hierarchy the visitation begins, i.e. Model or Region etc.

@param target for which the Visit() method will be called that actually
carries out the computation, for instance on the Node or on the Element.

*/
template<uint32_t dim>
Visitor<dim>::Visitor( PLACEMENT level, PLACEMENT target ) 
 : application_level_(level), 
   application_target_(target),
   verbose_(false)
 {
    if ( application_level_ > application_target_ )
      throw csmp::Exception( FATAL_ERROR, "Visitor<dim>(constructor):",
                            "The application level of a Visitor cannot be lower down in hierarchy than the target; read documentation." );
 }

template<uint32_t dim>
Visitor<dim>::~Visitor()
{
}

/**
 
The application level determines the type of object that the visitor
will be applied to by the Model. The visitor itself may then still
continue to migrate down the hierarchy of objects, but not upwards.
This method sets, and or (without argument) returns the application
level of the visitor on which it is called.  

@section arguments Input Arguments

The class of object that the visitor shall be applied to.

@attention If no argument is specified, the method will return the current setting
of target.

*/
template<uint32_t dim>
void  Visitor<dim>::ApplicationLevel( PLACEMENT p ) 
 { application_level_ = p; }

template<uint32_t dim>
PLACEMENT  Visitor<dim>::ApplicationLevel() const 
 { return application_level_; }




/**
 
The target of a visitation is the lowest level of the hierarchy
it will be applied to, i.e. Node, IntegrationPoint, Element or Face.  

If no argument is specified, the method will return the current setting
of target.

@section messages Messages 

If the argument supplied is not applicable, the method reports an error.  
*/
template<uint32_t dim>
void  Visitor<dim>::ApplicationTarget( PLACEMENT p )           
 { application_target_ = p; }


template<uint32_t dim>
csmp::PLACEMENT  Visitor<dim>::ApplicationTarget() const 
{ return application_target_; }


/**
 
The key method of the visitor design pattern. Visit() must be defined
for the target object (for instance Node objects, or Element objects).
Since the argument is the actual object, its entire interface is 
accessible to the visitor in order to carry out operations. Overload
this method in the derived class to gain access to target objects
and perform the actions that are desired for the visitation.  

@param f is a pointer to the target Face (forming part of a boundary).

*/
template<uint32_t dim>
void Visitor<dim>::Visit( Face<dim>* f )
 {
    cerr <<"\nVisitor::Visit(Face<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }


template<uint32_t dim>
void Visitor<dim>::Visit( Edge<dim>* f )
 {
    cerr <<"\nVisitor::Visit(Edge<3U>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }


template<uint32_t dim>
void Visitor<dim>::Visit( InterFace<dim>* )
 {
    cerr <<"\nVisitor::Visit(InterFace<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }
 
template<uint32_t dim>
void Visitor<dim>::Visit( Node<dim>* )
 {
    cerr <<"\nVisitor::Visit(Node<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }

template<uint32_t dim>
void Visitor<dim>::Visit( Element<dim>* )
 {
    cerr <<"\nVisitor::Visit(Element<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }
 
template<uint32_t dim>
void Visitor<dim>::Visit( Region<dim>* )
 {
    cerr <<"\nVisitor::Visit(Region<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }

template<uint32_t dim>
void Visitor<dim>::Visit( Boundary<dim>* )
 {
    cerr <<"\nVisitor::Visit(Boundary<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }
 
template<uint32_t dim>
void Visitor<dim>::Visit( SplitBoundary<dim>* )
 {
    cerr <<"\nVisitor::Visit(Boundary<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr << endl;
 }
 
template<uint32_t dim>
void Visitor<dim>::Visit( Model<dim>* )
 {
    cerr<<"\nVisitor::Visit(Model<dim>*): Method not implemented in the subclass of the Visitor class that you are using."<< endl;
 }

template<uint32_t dim>
void Visitor<dim>::SetInitialProperties( Model<dim>* )
 {
    cerr<<"\nVisitor::SetInitialProperties(Model<dim>*): Method not implemented in the subclass of the Visitor class that you are using.";
    cerr<<" Press <enter> to continue."<<endl;
    cin.get();
 }

template class Visitor<1U>;
template class Visitor<2U>;
template class Visitor<3U>;

} // end namespace csmp

