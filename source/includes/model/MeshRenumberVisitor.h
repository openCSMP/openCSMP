#ifndef CSMP_MESH_RENUMBER_VISITOR_H
#define CSMP_MESH_RENUMBER_VISITOR_H

#include "CSMP_definitions.h"
#include "Visitor.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Element;

/**
    Node and Element visitor that applies
    a unique Idx numbering in the range 0..n-1
    
    @author S.K. Matthai 
    @date 2006
*/
template<uint32_t dim>
class MeshRenumberVisitor final : public Visitor<dim> {
  public:
    MeshRenumberVisitor();
    
    void Visit( Element<dim>* ) override final;
    void Visit( Model<dim>* ) override final {}
    
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








