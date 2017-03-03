#ifndef CSMP_FEM_DATA_OUTPUT_VISITOR_H
#define CSMP_FEM_DATA_OUTPUT_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"

namespace csmp {

struct Index;
template<size_t> class PropertyDatabase;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Region;
template<typename> class FEM_Data;

/// replaces value of property a with propery b
template<typename Var, size_t dim>
class FEM_DataOutputVisitor : public Visitor<dim> {
  public:
    FEM_DataOutputVisitor( const PropertyDatabase<dim>&, 
                           const char* input_prop,
                           FEM_Data<Var>& data );

    virtual ~FEM_DataOutputVisitor();
    
    virtual void Visit( Region<dim>* );  
    virtual void Visit( Element<dim>* );  
    virtual void Visit( Node<dim>* );  
    
    void       Reset();

  private:
    FEM_Data<Var>&  input_data_ref_;
    csmp::Index     prop_key_;
    Var             variable_;
    size_t          counter_;
};

} // csmp


#endif







