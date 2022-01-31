#ifndef CSMP_DATA_INPUT_VISITOR_H
#define CSMP_DATA_INPUT_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Face;
template<size_t> class InterFace;
template<size_t> class Region;
template<size_t> class Model;

/**

   Use this visitor to input raw data into model
   which is assumed to be ordered in input vector
   by variable components (for example vector v1x, v1y, v1z, v2x...)
   visitor will renumber elements or nodes.
   
   @author S. K. Matthai
   @date 2006

*/
template<typename Var, size_t dim>
class DataInputVisitor : public Visitor<dim> {
  public:
    DataInputVisitor( Model<dim>& sg, 
                      const char* input_prop,
                      const std::vector<double>& input_data );

    virtual ~DataInputVisitor();
    
    virtual void Visit( Model<dim>* );
    virtual void Visit( Region<dim>* );
    virtual void Visit( Element<dim>* );  
    virtual void Visit( Face<dim>* );
    virtual void Visit( InterFace<dim>* );
    virtual void Visit( Node<dim>* );
    
    void  Reset();

  private:
    const std::vector<double>&  input_data_ref_;
    csmp::Index             prop_key_;
    Var                     variable_;
    size_t                  counter_;
};

} // csmp

#endif







