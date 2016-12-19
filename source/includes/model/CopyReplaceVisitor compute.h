#ifndef CSMP_COPY_REPLACE_VISITOR_H
#define CSMP_COPY_REPLACE_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"
#include "AP_BoolVector.h"

namespace csmp {

template<size_t> class PropertyDatabase;
template<size_t> class Node;
template<size_t> class Element;
template<size_t> class Region;



/** Replaces properties with the user supplied values 

    @author S.K. Matthai
    @date 2006

    @todo (2-F) Implement for new variable placements (FVIP)
*/
template<typename Var, size_t dim>
class CopyReplaceVisitor : public Visitor<dim> {
  public:
    CopyReplaceVisitor( const PropertyDatabase<dim>&, 
                        const char* prop_a, const char* prop_b, 
                        size_t nodes=0U ); // target property is not a node

    virtual ~CopyReplaceVisitor();
    
    virtual void Visit( Model<dim>* );
    virtual void Visit( Region<dim>* );
    virtual void Visit( Element<dim>* );  
    virtual void Visit( Node<dim>* );  

    void ComputeContribution( Element<dim>* );
    void ComputeContribution( Node<dim>* );
    
    void       Reset();

  private:
    BoolVector   nodes_visited_;
    csmp::Index  prop_key_a_, prop_key_b_;
    Var          variable_;
#if defined(_OPENMP ) 
    std::vector<Var> thread_variable_;
#endif
};


} // csmp


#endif







