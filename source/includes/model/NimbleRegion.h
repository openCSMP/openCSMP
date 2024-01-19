#ifndef CSMP_NIMBLE_REGION_H
#define CSMP_NIMBLE_REGION_H

#include "ModelSubDomain.h"

namespace csmp {

template<uint32_t> class Node;
template<uint32_t> class Element;
template<uint32_t> class Region;

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

    0. Must be a template with integral parameter uint32_t dim
 
    1. Be a suitable plugin into the PDE_Integrator, meaning that it has the methods:
 
      1.1 Must have a contructor that accepts iterators to node pointers
      1.2 Should be able to change shape at low computational cost
      1.3 Have the interfaces needed by the PDE_Integrator to function:
          - Must be based on vectors,
          - Know its perimeter nodes / have nodes sorted into interior and perimeter ranges
            so that boundary conditions can be applied.
 
    2. InterFace / task share with DES algorithm to convey update information.
*/
template<uint32_t dim>
class NimbleRegion : public ModelSubDomain<dim,Element>,
                     public LocalVariableStorage<dim,Region> {
  public:
    /// construction of region from nodes, relying on existing node-parent-element connectivity to identify elements
    NimbleRegion( const PropertyDatabase<dim>&,
                  typename std::vector<Node<dim>*>::const_iterator first,
                  typename std::vector<Node<dim>*>::const_iterator last );

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

    /// Local variable storage interface
    virtual PLACEMENT Placement() const { return REGION; }

    // retrieving information
  
    size_t Nodes() const;
    size_t InteriorNodes() const;
    size_t PerimeterNodes() const;
    size_t Cells() const;
    size_t RenumberNodes() const;


    // accessors
  
    Element<dim>* const E(size_t);
    Node<dim>* const N(size_t);

    const Element<dim>* const E(size_t) const;
    const Node<dim>* const N(size_t) const;
  
    /// writes the current element and node memberships to the console
    void Out() const;

  private:
    NimbleRegion() = delete;
    
    /// constructs region from supplied nodes, relying on existing node-parent-element connectivity to identify elements and perimeter nodes
    void Update( typename std::vector<Node<dim>*>::const_iterator first, typename std::vector<Node<dim>*>::const_iterator last );

    /// as above, but with different way to find perimeter (FAIL: perimeter incorrect for discontiguous patches)
    void Update2( typename std::vector<Node<dim>*>::const_iterator first, typename std::vector<Node<dim>*>::const_iterator last );
    
    void Update3( typename std::vector<Node<dim>*>::const_iterator first, typename std::vector<Node<dim>*>::const_iterator last );
  
  private:
    // current implementation based on idea that vectors are resized with little overhead as long as their capacity is not changed
    const bool                  verbose_ = true;    ///< flag for testing and reporting
                       
    friend class NimbleRegion_Test;
};

} // csmp   

#endif /* CSMP_NIMBLE_REGION_H */
