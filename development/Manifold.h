#ifndef MANIFOLD_H
#define MANIFOLD_H

#include <vector>

namespace csmp
{

template<size_t Dim>
class Node;

class Index;

/**
* @brief: a class contains coincident nodes at split-boundary
* @date: 2020/02/02
*/

template<size_t Dim>
class Manifold
{
public:
  /// constructor
  explicit Manifold(std::vector<Node<Dim>*>);
  
  /// copy constructor
  Manifold(const Manifold<Dim>&) = delete;
  
  /// asignment constructor 
  Manifold(Manifold<Dim>&&) = delete;
  
  /// destructor
  ~Manifold();

  /// assigment operator
  Manifold<Dim>& operator = (const Manifold<Dim>&) = delete;
  
  /// asscending sort
  void Sort(const Index&);

  /// number of nodes 
  size_t Nodes() const;

  /// access to nodes
  Node<Dim>* N(const size_t&);

  /// print out information
  void Out();

private:
  std::vector<Node<Dim>*> nodes_;
};
  
}

#endif

