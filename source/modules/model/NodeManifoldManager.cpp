#include "NodeManifoldManager.h"
#include "ErrorHandler.h"


using namespace std;

namespace csmp {

template<size_t dim>
NodeManifoldManager<dim>::NodeManifoldManager()
  {
  }

template<size_t dim>
NodeManifoldManager<dim>::~NodeManifoldManager()
  {
  }

template<size_t dim>
typename std::deque<NodeManifold<dim>* >::iterator NodeManifoldManager<dim>::NodeManifoldsBegin()
  { return node_manifolds_.begin(); }


template<size_t dim>
typename std::deque<NodeManifold<dim>* >::iterator NodeManifoldManager<dim>::NodeManifoldsEnd()
  { return node_manifolds_.end(); }

template<size_t dim>
typename std::deque<NodeManifold<dim>* >::const_iterator NodeManifoldManager<dim>::NodeManifoldsBegin() const
  { return node_manifolds_.begin(); }


template<size_t dim>
typename std::deque<NodeManifold<dim>* >::const_iterator NodeManifoldManager<dim>::NodeManifoldsEnd() const
  { return node_manifolds_.end(); }

/*
template<size_t dim>
typename std::deque< Manifold<dim>* >::iterator NodeManifoldManager<dim>::Manifold( const Manifold<dim>* md )
  { 
    for( manifoldIterator it( ManifoldsBegin() ); it != ManifoldsEnd(); ++it )
      if( &it == md )
        return it;
    return ManifoldsEnd();
  }
*/

template<size_t dim>
size_t  NodeManifoldManager<dim>::NodeManifolds() const
  { return node_manifolds_.size(); }

template<size_t dim>
bool NodeManifoldManager<dim>::InsertNodeManifold(NodeManifold<dim>* md){
    if(std::find(node_manifolds_.begin(), node_manifolds_.end(), md) == node_manifolds_.end()) {
        node_manifolds_.push_back(md);
        return true;
    } else {
        ErrorHandler&  csmp_error( ErrorHandler::Instance() );
        csmp_error.notice( INFO, "NodeManifoldManager<dim>::InsertNodeManifold:",
                                 "NodeManifold already exists. Nothing was done.");   
    }
    
    return false;
}

template<size_t dim>
bool NodeManifoldManager<dim>::RemoveNodeManifold(NodeManifold<dim>* md) 
{
    if(std::find(node_manifolds_.begin(), node_manifolds_.end(), md) != node_manifolds_.end()) {
        remove(node_manifolds_.begin(),node_manifolds_.end(),md);
        return true;
    } else {
        ErrorHandler&  csmp_error( ErrorHandler::Instance() );
        csmp_error.notice( INFO, "NodeManifoldManager<dim>::RemoveNodeManifold:",
                                 "NodeManifold does not exist. Nothing was done.");                          
    }
    
    return false;
}

template<size_t dim>
bool NodeManifoldManager<dim>::AddNodeManifold( std::vector<Node<dim>*> nodes )
{
    NodeManifold<dim>* md = new NodeManifold<dim>(nodes);
    return InsertNodeManifold(md);
}

template<size_t dim>
bool NodeManifoldManager<dim>::AddNodeManifold( std::vector<Node<dim>*> nodes,
                                               PARENT_GEOMETRIC_ENTITY parent_geometry )
{
    NodeManifold<dim>* md = new NodeManifold<dim>(nodes, parent_geometry);
    return InsertNodeManifold(md);
}

template<size_t dim>
bool NodeManifoldManager<dim>::AddNodeManifold(typename std::vector<Node<dim>*>::iterator nodesBegin,
                                               typename std::vector<Node<dim>*>::iterator nodesEnd)
{
    NodeManifold<dim>* md = new NodeManifold<dim>(nodesBegin, nodesEnd);
    return InsertNodeManifold(md);
}

template<size_t dim>
bool NodeManifoldManager<dim>::AddNodeManifold(typename std::vector<Node<dim>*>::iterator nodesBegin,
                                               typename std::vector<Node<dim>*>::iterator nodesEnd,
                                               PARENT_GEOMETRIC_ENTITY parent_geometry)
{
    NodeManifold<dim>* md = new NodeManifold<dim>(nodesBegin, nodesEnd, parent_geometry);
    return InsertNodeManifold(md);
}

template<size_t dim>
bool NodeManifoldManager<dim>::InsertNodeToManifold(NodeManifold<dim>* md, Node<dim>* nd)
{
    return md->Add(nd);
}

template<size_t dim>
bool NodeManifoldManager<dim>::RemoveNodeFromManifold(NodeManifold<dim>* md, Node<dim>* nd)
{
    return md->Delete(nd);
}



template class NodeManifoldManager<1U>;
template class NodeManifoldManager<2U>;
template class NodeManifoldManager<3U>;

} // csmp
