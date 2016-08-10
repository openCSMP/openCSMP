#ifndef CSMP_MESH_RENUMBER_VISITOR_H
#define CSMP_MESH_RENUMBER_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;

/**
    Node and Element visitor that applies
    a unique Idx numbering in the range 0..n-1
    
    @author S.K. Matthai 
    @date 2006
*/
template<size_t dim>
class MeshRenumberVisitor : public Visitor<dim> {
  public:
    MeshRenumberVisitor();

    virtual ~MeshRenumberVisitor();
    
    virtual void Visit( Element<dim>* );  
    
    size_t  VisitedElements() const;
    size_t  VisitedNodes() const; 
    
    void    Reset();

  private:
    std::set<Node<dim>*>  node_ptrs_;
    size_t                e_counter_;
    size_t                n_counter_;
};

} // csmp

#endif








