
#include "NodeManifold.h"
#include "Node.h"
#include "Element.h"
#include "Exception.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<size_t dim>
NodeManifold<dim>::NodeManifold()
{
}

template<size_t dim>
NodeManifold<dim>::NodeManifold( const vector<Node<dim>*>& nodes )
 : nodes_(nodes),
   parent_geometry_(NOT_SPECIFIED)
{
    assert( nodes_.size() > 1 );
    for( auto nd : nodes_ )
      nd->AssignParentManifold(this);
}

template<size_t dim>
NodeManifold<dim>::NodeManifold( const vector<Node<dim>*>& nodes, PARENT_GEOMETRIC_ENTITY parent_geometry )
 : nodes_(nodes),
   parent_geometry_(parent_geometry) 
{
    assert(nodes_.size() > 1);
    for( auto nd : nodes_ )
      nd->AssignParentManifold(this);
}

template<size_t dim>
NodeManifold<dim>::NodeManifold( typename std::vector<Node<dim>*>::iterator nodesBegin,
                                 typename std::vector<Node<dim>*>::iterator nodesEnd )
: parent_geometry_(NOT_SPECIFIED)
{
    assert(nodesBegin != nodesEnd);
    nodes_.assign(nodesBegin, nodesEnd);
    for(auto nd : nodes_)
        nd->AssignParentManifold(this);    

} 

template<size_t dim>
NodeManifold<dim>::NodeManifold( typename std::vector<Node<dim>*>::iterator nodesBegin,
                                 typename std::vector<Node<dim>*>::iterator nodesEnd,
                                 PARENT_GEOMETRIC_ENTITY parent_geometry )
: parent_geometry_(parent_geometry)
{
    assert(nodesBegin != nodesEnd);
    nodes_.assign(nodesBegin, nodesEnd);
    for(auto nd : nodes_)
        nd->AssignParentManifold(this);     

} 

template<size_t dim>
NodeManifold<dim>::~NodeManifold()
{
    for(auto nd : nodes_)
        nd->AssignParentManifold(nullptr);    
  
    nodes_.clear();
}


/**
    Copy constructor also copies the pointer assignments (!).
*/
template<size_t dim>
NodeManifold<dim>::NodeManifold( const NodeManifold& md)
 : nodes_(md.nodes_),
   parent_geometry_(md.parent_geometry_)   
{
}

/**
    Move constructor also copies the pointer assignments (!).
*/
template<size_t dim>
NodeManifold<dim>::NodeManifold( NodeManifold&& md)
  : nodes_{md.nodes_},
    parent_geometry_(md.parent_geometry_)  
{
}

template<size_t dim>
NodeManifold<dim>& NodeManifold<dim>::operator=( const NodeManifold<dim>& md )
 {
    if ( &md != this ) 
      {
         nodes_ = md.nodes_;
         parent_geometry_ = md.parent_geometry_;  
      }
    return *this;
 }

/**
    Move assignment, relying on that similar operators exist for the nodes components.
    
    @note this assumes that the supplied node is a temporary.
*/
template<size_t dim>
NodeManifold<dim>& NodeManifold<dim>::operator=( NodeManifold<dim>&& md )
 {
    assert( this != &md );
    nodes_ = md.nodes_;
    parent_geometry_ = md.parent_geometry_;  
 
    return *this;
 }



template<size_t dim>
void NodeManifold<dim>::SortByVariableIndex(const Index& index)
{
    assert(index.type==SCALAR);
    assert(index.place==NODE || index.place==ELEMENT);

    if(index.type!=SCALAR) {
        throw csmp::Exception( INFO, "NodeManifold<dim>::SortByVariableIndex",
                                     "variable type not supported, must be SCALAR");       
    }
 
    if(index.place==NODE){
        for (size_t i = 0; i < nodes_.size(); ++i) {
            auto first = nodes_[i]->Read(index);
            assert(!isnan(first));
            for (size_t j = i + 1; j < nodes_.size(); ++j) {
                auto second = nodes_[j]->Read(index);
                assert(!isnan(second));
                if (first > second) std::swap(nodes_[i], nodes_[j]);
            }
        }    
    } else if(index.place==ELEMENT) {
        for (size_t i = 0; i < nodes_.size(); ++i) {
            //double64 first = std::numeric_limits<double64>::quiet_NaN();
            double64 first_value(0.);
            size_t first_count(0);
            for(size_t e = 0; e < nodes_[i]->Parents(); e++) {
                //if( (dim==2 && nodes_[i]->Parent(e)->IsSurfaceElement()) || (dim==3 && nodes_[i]->Parent(e)->IsVolumeElement()) ) {
                if( !isnan(nodes_[i]->Parent(e)->Read(index)) ) {    
                    //first = nodes_[i]->Parent(e)->Read(index);
                    //break;
                    first_value += nodes_[i]->Parent(e)->Read(index);
                    first_count ++;
                }
            }
            //assert(!isnan(first));
            assert(first_count != 0);
            first_value /= first_count;
            assert(!isnan(first_value));

            for (size_t j = i + 1; j < nodes_.size(); ++j) {
                //double64 second = std::numeric_limits<double64>::quiet_NaN();
                double64 second_value(0.);
                size_t second_count(0);                
                for(size_t e = 0; e < nodes_[j]->Parents(); e++) {
                    //if( (dim==2 && nodes_[j]->Parent(e)->IsSurfaceElement()) || (dim==3 && nodes_[j]->Parent(e)->IsVolumeElement()) ) {
                    if( !isnan(nodes_[j]->Parent(e)->Read(index)) ) {    
                        //second = nodes_[j]->Parent(e)->Read(index);
                        //break;
                        second_value += nodes_[j]->Parent(e)->Read(index);
                        second_count ++;
                    }
                }
                //assert(!isnan(second));
                assert(second_count != 0);
                second_value /= second_count;
                assert(!isnan(second_value));                
                //if (first > second) std::swap(nodes_[i], nodes_[j]);
                if (first_value > second_value) std::swap(nodes_[i], nodes_[j]);
            }
        } 
    } else {
        throw csmp::Exception( INFO, "NodeManifold<dim>::SortByVariableIndex",
                                     "variable placement not supported");      

    }
}


template<size_t dim>
size_t NodeManifold<dim>::Nodes() const
{
  return nodes_.size();
}

template<size_t dim>
Node<dim>* NodeManifold<dim>::N(const size_t& index)
{
  return nodes_[index];
}

template<size_t dim>
bool NodeManifold<dim>::Add(Node<dim>* nd)
{
    //a node can only exist in one manifold
    if(nd->ParentManifold() != nullptr) {
        ErrorHandler&  csmp_error( ErrorHandler::Instance() );
        csmp_error.notice( INFO, "NodeManifold<dim>::Add(Node<dim>* nd)",
                                 "Node already in another manifold. Nothing was done.");         
        return false;
    }

    //the node is already in current manifold
    if(std::find(nodes_.begin(), nodes_.end(), nd) != nodes_.end()) {
        ErrorHandler&  csmp_error( ErrorHandler::Instance() );
        csmp_error.notice( INFO, "NodeManifold<dim>::Add(Node<dim>* nd)",
                                 "Node already exist. Nothing was done."); 
        return false;
    }

    nodes_.push_back(nd);  
    nd->AssignParentManifold(this);  

    return true;

} 

template<size_t dim>
bool NodeManifold<dim>::Delete(Node<dim>* nd)
{    
    if(std::find(nodes_.begin(), nodes_.end(), nd) != nodes_.end()) {
        remove(nodes_.begin(),nodes_.end(),nd);
        nd->AssignParentManifold(nullptr); 
        return true;
    } else {
        ErrorHandler&  csmp_error( ErrorHandler::Instance() );
        csmp_error.notice( INFO, "NodeManifold<dim>::Delete(Node<dim>* nd)",
                                 "Node does not exist. Nothing was done.");                          
    }
    
    return false;

} 

template<size_t dim>
void NodeManifold<dim>::Out() 
{
  std::cout << "NodeManifold with nodes IDs:\t";
  for (auto const& node : nodes_) 
  {
    std::cout << node->Idx() << "\t";
  }
  std::cout << std::endl;
}

template class NodeManifold<1U>;
template class NodeManifold<2U>;
template class NodeManifold<3U>;

}// csmp
