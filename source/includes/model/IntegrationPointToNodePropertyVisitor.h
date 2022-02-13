#ifndef CSMP_INTEGRATION_POINT_TO_NODE_PROPERTY_VISITOR_H
#define CSMP_INTEGRATION_POINT_TO_NODE_PROPERTY_VISITOR_H

#include "Visitor.h"
#include "AP_BoolVector.h"

namespace csmp {

template<uint32_t> class PropertyDatabase;
template<uint32_t> class Node;
template<uint32_t> class Element;

/** replaces value of property a with propery b

@todo (1) IntegrationPointToNodePropertyVisitor is not extrapolating correctly (A)
@todo (1) SKM: consider using the element specific method ExtrapolateIntegrationPointToNodeProperty for this task

*/
template<typename Var, uint32_t dim>
class IntegrationPointToNodePropertyVisitor : public Visitor<dim> {
  public:
    IntegrationPointToNodePropertyVisitor( const PropertyDatabase<dim>& p, 
                                          const char* cpoint_prop,
                                          const char* node_prop,
                                          size_t nodes );

    virtual ~IntegrationPointToNodePropertyVisitor();
    
    virtual void Visit( Node<dim>* );  
    virtual void Visit( Element<dim>* );  
    
    void ApplyWeightingToExtrapolatedValues();

  private:
    csmp::Index       cprop_key_, nprop_key_;
    Var               variable_;
    std::vector<Var>  vars_vector_; 
    std::vector<double>   cp_vars_, nd_vars_, 
                      summed_weights_;
    BoolVector        weighting_completed_;
};


} // csmp


#endif







