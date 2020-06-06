#include <assert.h>
#include <iostream>
#include "Manifold.h"
#include "Node.h"
#include "Element.h"
namespace csmp
{

template<size_t Dim>
Manifold<Dim>::Manifold(std::vector<Node<Dim>*> nodes) :
nodes_(nodes) 
{
  assert(nodes_.size() > 1, "Manifold constructor: input has less than 2 nodes");
}

template<size_t Dim>
Manifold<Dim>::~Manifold()
{
  nodes_.clear();
}

template<size_t Dim>
void Manifold<Dim>::Sort(const Index& index)
{
  for (size_t i = 0; i < nodes_.size(); ++i) 
  {
    auto first = nodes_[i]->Parent(0)->Read(index);
    for (size_t j = i + 1; j < nodes_.size(); ++j)
    {
      auto second = nodes_[j]->Parent(0)->Read(index);
      if (first > second) 
      {
        std::swap(nodes_[i], nodes_[j]);
      }
    }
  }
}


template<size_t Dim>
size_t Manifold<Dim>::Nodes() const
{
  return nodes_.size();
}

template<size_t Dim>
Node<Dim>* Manifold<Dim>::N(const size_t& index)
{
  return nodes_[index];
}

template<size_t Dim>
void Manifold<Dim>::Out() 
{
  std::cout << "Manifold with nodes IDs:\t";
  for (auto const& node : nodes_) 
  {
    std::cout << node->Idx() << "\t";
  }
  std::cout << std::endl;
}

template class Manifold<1U>;
template class Manifold<2U>;
template class Manifold<3U>;

}// csmp