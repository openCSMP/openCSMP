//
//  NodeFunctions_Test.h
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 7/07/2024.
//

#ifndef CSMP_NODE_FUNCTIONS_TEST_H
#define CSMP_NODE_FUNCTIONS_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"
#include "Node.h"

namespace csmp {

template<uint32_t> class Model;
template<uint32_t> class Region;

/**

TODO: these functions still need to be testsed
             
/// if there is only one parent element expected, then use this method instead of 'parentElementsSharedByFace'
@code
template<uint32_t dim>
std::pair<Element<dim>*,size_t>  parentElement( std::vector<Node<dim>*>::const_iterator first,
                                                std::vector<Node<dim>*>::const_iterator last );
@endcode

/// for a node that lies on an internal surface, method reports which of volumetric parent elements lie on the inside and which on the outside of this surface; throws if assumptions are not met
@code
template<uint32_t dim>
std::pair<std::vector<Element<dim>*>,std::vector<Element<dim>*>>  parentElementsAdjacentTo( const Node<dim>* const low_dim_node );
@endcode

/// determines role of Node (simple mesh node vs. geometric constraint), using BOX flagging and parent element connectivity,
@code
template<uint32_t dim>
TOPOTYPE checkModelPartThatNodeBelongsTo( const Node<dim>* const );
@endcode

/// prints Idx and boundary flag values of the nodes that this node is connected with
@code
template<uint32_t dim>
void printNeighbors( const Node<dim>* const );
@endcode

/// prints current parent information and checks for duplicate parents
@code
template<uint32_t dim>
void printParents( const Node<dim>* const );
@endcode

/// calculates the size of the Node excluding the dynamic contribution to the stored variables
@code
template<uint32_t dim>
size_t sizeOf( const Node<dim>* const );
@endcode

    @author SKM
    @date 7/7/24
*/
class NodeFunctions_Test : public Test {
  public:
    NodeFunctions_Test();
    virtual ~NodeFunctions_Test() {}
    
    virtual void run();
  
  private:

    // testing whether functionality is suitable for creating a sparsity pattern
    template<uint32_t dim>
    bool Test_nodeNeighbors( const Model<dim>& model );
    
    bool Test_parentElementsSharedByFace();
    
    const bool verbose_ = true;
};

} // end csmp

#endif /* CSMP_NODE_FUNCTIONS_TEST_H */
