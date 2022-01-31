#ifndef CSMP_ELEMENT_TO_NODE_PROPERTY_VISITOR_H
#define CSMP_ELEMENT_TO_NODE_PROPERTY_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "AP_BoolVector.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Node;
template<size_t> class Element;

/// replaces value of property a with propery b
template<typename Var, size_t dim>
class ElementToNodePropertyVisitor : public Visitor<dim> {
  public:
    ElementToNodePropertyVisitor( const PropertyDatabase<dim>& p,
                                  const char* elmt_prop,
                                  const char* node_prop,
                                  size_t nodes );

    virtual ~ElementToNodePropertyVisitor();
    
    virtual void Visit( Node<dim>* );  
    virtual void Visit( Element<dim>* );  
    
    void ApplyWeightingToExtrapolatedValues();

  private:
    csmp::Index       eprop_key_, nprop_key_;
    Var               evariable_, nvariable_;
    std::vector<Var>  vars_vector_; 
    std::vector<double>   summed_weights_;
    BoolVector        weighting_completed_;
};


} // csmp


#endif







