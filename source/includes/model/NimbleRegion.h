#ifndef CSMP_NIMBLE_REGION_H
#define CSMP_NIMBLE_REGION_H

#include "CSMP_definitions.h"
#include "Exception.h"
#include "ErrorHandler.h"

namespace csmp {

template<size_t> class Node;
template<size_t> class Element;

/**
    Flexible collection of Element and Node object pointers
    to be used by the PDE_Integrator, which also requires access
    to ranges of (perimeter) nodes.
 
    NimbleRegion is a light (non-unique) version of Region that can grow or shrink,
    retaining the knowledge of its perimeter.
 
    Like Region, NimbleRegion keeps track of its perimeter, which
    always includes a halo that comes about as all elements that
    contain user-supplied nodes (and additional ones) will
    be included into the contructed region.
 
    The region is based on the STL vector and operated under a high-memory policy.
    This means that the Update() and other functions retain previous vector storage
    allocations, just resizing these vectors to zero when a reconstruction is performed,
    to avoid costly memory allocations.
    To delete both elements and the internal storage of the vectors, delete and create a new
    NimbleRegion object.
 
    @attention NimbleRegion does not have a name or properties of its own and
    it is not managed by the RegionInterface,
    create it using its constructor and manage it with DES or any other parallel algorithm.
 
    @attention Grow(), Shrink() and ChangeShape() operations rely on the connectivity
    between equidimensional elements and nodes.
 
    @attention Thus far, Nimbleregion only works with linear elements.
 
    Design specifications in detail: NimbleRegion for use in PDE_Integrators?

    0. Must be a template with integral parameter size_t dim
 
    1. Be a suitable plugin into the PDE_Integrator, meaning that it has the methods:
 
      1.1 Must have a contructor that accepts iterators to node pointers
      1.2 Should be able to change shape at low computational cost
      1.3 Have the interfaces needed by the PDE_Integrator to function:
          - Must be based on vectors,
          - Know its perimeter nodes / have nodes sorted into interior and perimeter ranges
            so that boundary conditions can be applied.
 
    2. InterFace / task share with DES algorithm to convey update information.
*/
template<size_t dim>
class NimbleRegion {
  public:
    /// for flexibility with regard to application domain
    typedef Element<dim>  CellType;
    
  public:
    /// construction of region from nodes, relying on existing node-parent-element connectivity to identify elements
    NimbleRegion( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );

    /// constructs region from supplied nodes, relying on existing node-parent-element connectivity to identify elements and perimeter nodes
    void Update( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );

    /// as above, but with different way to find perimeter (FAIL: perimeter incorrect for discontiguous patches)
    void Update2( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );
    
    void Update3( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );
  
/* CARRY OUT DIAGNOSTICS WHETHER THESE INTERFACES ARE WORTH IMPLEMENTING

    /// adjust to arbitrary yet small shape modifications, where most of the region stays the same; argument new current nodes
    void ChangeShape( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );
    void ChangeShape( typename std::set<Node<dim>*>& nodes_to_remove, typename std::set<Node<dim>*>& nodes_to_add );

    /// adds multiple nodes and potential extra elements to region, does not remove any nodes or elements
    void Grow( typename std::vector<Node<dim>*>::iterator first, typename std::vector<Node<dim>*>::iterator last );
 
    /// removes nodes and elements that might have become disconnected from the region
    void Shrink( std::vector<Node<dim>*>& );
*/
    /// resize to zero while keeping the memory
    void Clear();
  
    /// removing storage
    void Erase();

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
    
    std::string Name() const
        {throw csmp::Exception( ERROR, "NimbleRegion<dim>::Name","Method not implemented");};
    bool IsPerimeterNode( const csmp::Node<dim>* ) const 
        {throw csmp::Exception( ERROR, "NimbleRegion<dim>::IsPerimeterNode","Method not implemented");};
    void UpdateMemberIndexes() const
        {throw csmp::Exception( ERROR, "NimbleRegion<dim>::UpdateMemberIndexes","Method not implemented");};    
    bool IsPerimeterElement( const size_t eidx ) const
        {throw csmp::Exception( ERROR, "NimbleRegion<dim>::IsPerimeterElement","Method not implemented");};  
    size_t SharedPerimeterNodes(typename std::vector<csmp::Node<dim>*>::const_iterator start, 
                                typename std::vector<csmp::Node<dim>*>::const_iterator end ) const
        {throw csmp::Exception( ERROR, "NimbleRegion<dim>::SharedPerimeterNodes","Method not implemented");};  
    

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
