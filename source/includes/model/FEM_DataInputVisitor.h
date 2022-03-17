#ifndef CSMP_FEM_DATA_INPUT_VISITOR_H
#define CSMP_FEM_DATA_INPUT_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"

namespace csmp {

/// to be applied to entire model
template<uint32_t> class PropertyDatabase;
template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Region;
template<uint32_t> class Model;
template<typename> class FEM_Data;

template<typename Var, uint32_t dim>
class FEM_DataInputVisitor : public Visitor<dim> {
  public:
    FEM_DataInputVisitor( Model<dim>&, 
                          const char* input_prop,
                          const FEM_Data<Var>& );

    virtual ~FEM_DataInputVisitor();
    
    virtual void Visit( Region<dim>* );  
    virtual void Visit( Element<dim>* );  
    virtual void Visit( Node<dim>* );  
    
    void       Reset();

  private:
    const FEM_Data<Var>&  input_data_ref_;
    csmp::Index           prop_key_;
    Var                   variable_;
    size_t                counter_;
};

} // csmp

#endif







