#ifndef CSMP_VISITOR_H
#define CSMP_VISITOR_H

#include "CSMP_definitions.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Face;
template<uint32_t> class Edge;
template<uint32_t> class InterFace;
template<uint32_t> class Element;
template<uint32_t> class Model;
template<uint32_t> class Region;
template<uint32_t> class Boundary;
template<uint32_t> class SplitBoundary;

/**

Visitors are designed to facilitate computations on any of the objects
that are contained inside a Model. Thus, one may for instance derive
an equation of state type visitor to go around the nodes looking up
pressure and temperature and use this information to calculate nodal
fluid densities. As an important difference to Interrelations, visitors
are not restricted to the computation of a single variable and they
can also be used to gain access to the full finite-element or 
finite-volume functionality of CSMP. Thus, one can use a visitor to
find a specific finite element that contains a user-specified point. This
is done by across mesh traversal. In contrast, interrelations are restricted to 
sequential node-by-node, element-by-element calculations.  
 
@author S.K. Matthai
@author Stephen G. Roberts
@author S. Geiger
@date 2001

@section motivation Motivation

To make the classical Visitor design pattern available to CSMP users
(for documentation of pattern, refer to Gamma et al. 1996, Addison
and Wesley, p. 331). 
 
@section applicability Applicability

Visitors in CSMP are designed to allow visitations of the Model
hierarchy of objects: Model > Region > Element > > Node, 
as well as interfaces. The user can derive visitors from this base class
in order to perform different operations on any of the objects
encountered during a visitation.  
 
 
@section examples Application Examples

Have a look at the AlterationVisitor of the equation of state visitors
in order to understand this design pattern.

*/
template<uint32_t dim>
class Visitor {
  public:
    explicit Visitor( PLACEMENT level=MODEL, PLACEMENT target=ELEMENT );
    virtual ~Visitor();

    // the highest level in the hierarchy (Model, Region, Element, Node)
    // where the visitor will start to make changes
    void         ApplicationLevel( PLACEMENT p );
    PLACEMENT    ApplicationLevel() const;

    // the lowest level in the hierarchy where visitor will make changes
    void         ApplicationTarget( PLACEMENT p );
    PLACEMENT    ApplicationTarget() const;

    virtual void SetInitialProperties(Model<dim>* model); //Added to allow flexibility in the creation and usage of a visitors that need specific properties
                                                          //to exist/be-calculated before they can be properly initialized.  Model is added as a parameter
                                                          //to aid the extraction of needed parameters. (J.E.M. 07.01.2016)
    virtual void Visit( Model<dim>* );   
    virtual void Visit( Region<dim>* );   
    virtual void Visit( Boundary<dim>* );   
    virtual void Visit( SplitBoundary<dim>* );
    virtual void Visit( Element<dim>* ); 
    virtual void Visit( Face<dim>* );
    virtual void Visit( Edge<dim>* );
    virtual void Visit( InterFace<dim>* );
    virtual void Visit( Node<dim>* );   

    bool Verbose(){return this->verbose_;}
    void Verbose(bool verbose){this->verbose_=verbose;}

  protected:
    PLACEMENT    application_level_;
    PLACEMENT    application_target_;
  private:
    bool verbose_;
};




} // csmp

#endif
