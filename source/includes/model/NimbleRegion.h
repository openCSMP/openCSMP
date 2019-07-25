#ifndef CSMP_NIMBLE_REGION_H
#define CSMP_NIMBLE_REGION_H

#include "CSMP_definitions.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;

/**
    Flexible region of Element objects.
 
    A light (non-unique) version of Region that can grow or shrink,
    retaining the knowledge of its perimeter.

    @attention NimbleRegion does not have a name or properties of its own and
    it is not managed by the RegionInterface, but by DES or any other parallel algorithm.
 
    @attention Thus far, Nimbleregion only works with linear elements.
 
    ===========================================================

   Design specs: What should such a quick region be able to do?

    0. Must be a template with integral parameter size_t dim
    1. Be a suitable plugin into the PDE_Integrator, meaning that it has the methods:
    2. Must have a contructor that accepts iterators to node pointers
    3. Should be able to grow or shrink at low computational cost
    4. have the interfaces that are needed by the PDE_Integrator to function:

   5. Since it is based on vectors, it must run high-mem so that these do not need to be reallocated all the time.
      This means that are not be pruned back after creation.

   6. may or may not have a name (Name() or ID).

   7. Know its perimeter nodes / have nodes sorted into interior and perimeter ranges.
      This is not critical because the PDE_Integrator finds essential conditions by looking at node flags anyway.
 
      Create method, that allows DES algorithm to convey this information.
*/
template<size_t dim>
class NimbleRegion {
  public:
    /// for flexibility with regard to application domain
    typedef Element<dim>  CellType;

    /// constructs region from supplied nodes, relying on existing node-parent-element connectivity to identify elements
    void Initialise( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );
  
    /// construction of region from nodes, relying on existing node-parent-element connectivity to identify elements
    NimbleRegion( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );

    /// adds multiple nodes and potential extra elements to region, does not remove any nodes or elements
    void Grow( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );
 
    /// removes nodes and elements that might have become disconnected from the region
    void Shrink( std::vector<Node<dim>*>& );

    /// drop all storage
    void Clear();

    // retrieving information
  
    size_t Nodes() const;
    size_t InteriorNodes() const;
    size_t PerimeterNodes() const;
    size_t Elements() const;
    size_t RenumberNodes() const;

    // accessors
  
    CellType* const E(size_t);
    Node<dim>* const N(size_t);

    const CellType* const E(size_t) const;
    const Node<dim>* const N(size_t) const;

    // iterators
  
    typename std::vector<csmp::Node<dim>*>::iterator        NodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator        PerimeterNodesBegin();
    typename std::vector<csmp::Node<dim>*>::iterator        NodesEnd();
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator  PerimeterNodesBegin() const;
    typename std::vector<csmp::Node<dim>*>::const_iterator  NodesEnd() const;

    typename std::vector<CellType*>::iterator               ElementsBegin();
    typename std::vector<CellType*>::iterator               ElementsEnd();
    typename std::vector<CellType*>::const_iterator         ElementsBegin() const;
    typename std::vector<CellType*>::const_iterator         ElementsEnd() const;
  
    /// writes the current element and node memberships to the console
    void Out() const;

  private:
    NimbleRegion() = delete;
  
    // current implementation based on idea that vectors are resized with little overhead as long as their capacity is not changed
    // sets are used to keep nodes and elements unique
    std::vector<Node<dim>*>     nodes_;             ///< sorted into interior and perimeter ranges
    size_t                      n_interior_nodes_;  ///< first perimeter node
    std::vector<Element<dim>*>  elements_;          ///< all elements, interior and exterior
    const bool                  verbose_ = true;    ///< flag for testing and reporting
};

} // csmp   

#endif /* CSMP_NIMBLE_REGION_H */
