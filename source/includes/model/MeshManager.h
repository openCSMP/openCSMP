#ifndef CSMP_MESH_MANAGER_H
#define CSMP_MESH_MANAGER_H

#include <deque>
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "PrimitiveContainer.h"

namespace csmp {

class FiniteElement;
class FiniteElementManager;
class ModelTopology;
template<size_t> class PropertyDatabase;
template<size_t> class VSet;
template<size_t> class FiniteVolumeStencilManager;

/**
 
@brief Helper class of the Model which takes care of the storage of Element, Face and InterFace objects;
internal application is hidden and may vary between models (deque-storage is default).
 
@author S.K. Matthaei
@date 2007

@remark gain access via Mesh() public interface of Model.

@attention the MeshManager takes care of the creation and destruction of Elements, Faces or Interfaces.
Region or Boundary objects merely contain pointers to these.

*/
template<size_t dim>
class MeshManager {
  public:
    MeshManager();
    MeshManager( const PropertyDatabase<dim>&, const FiniteElementManager&, VSet<dim>& );
    MeshManager&  operator=( const MeshManager& );
    MeshManager( const MeshManager& );
    ~MeshManager();

    /// sets up distributed storage for variables, finite elements, and mesh connectivity
    bool Initialize( const PropertyDatabase<dim>&,
                     const FiniteElementManager&, 
                     const VSet<dim>& );

    /// from verified mesh and property data including Face and Interface objects
    bool Reconstruct( const PropertyDatabase<dim>&,
                      const FiniteElementManager&,
                      const VSet<dim>& );
  
    /// assigns finite volume stencils to the FV pointers stored in each element
    void InitializeFiniteVolumeStencils( const PropertyDatabase<dim>&,
                                         const FiniteElementManager&, 
                                         FiniteVolumeStencilManager<dim>& );

    /// returns true if the mesh consists of multiple element types
    bool    HybridElementMesh() const;
  
    /// counts and returns current indices of elements that may give rise to problems during the assignment of boundary conditions
    size_t DetectElementsWithAllNodesOnBoundary( std::set<size_t>& ) const;

    /// returns number of nodes=vertices in the current mesh
    size_t  Nodes() const;
  
    /// returns number of elements in the current mesh
    size_t  Elements() const;
  
    /// returns number of Faces=lower-dimensional elements in current mesh
    size_t  Faces() const;
  
    /// returns number of interfaces=faces with multiplicated nodes
    size_t  InterFaces() const;

    /// the node from which the node tree is build
    csmp::Node<dim>&        RootNode();
  
    /// the element at the basis of element tree
    csmp::Element<dim>&     RootElement();
  
    csmp::Face<dim>&        RootFace();
    csmp::InterFace<dim>&   RootInterFace();
  
    // mutator iterators
    typename PrimitiveContainer<csmp::Node<dim> >::iterator        NodesBegin();
    typename PrimitiveContainer<csmp::Node<dim> >::iterator        NodesEnd();
    typename PrimitiveContainer<csmp::Element<dim> >::iterator     ElementsBegin();
    typename PrimitiveContainer<csmp::Element<dim> >::iterator     ElementsEnd();
    typename PrimitiveContainer<csmp::Face<dim> >::iterator        FacesBegin();
    typename PrimitiveContainer<csmp::Face<dim> >::iterator        FacesEnd();
    typename PrimitiveContainer<csmp::InterFace<dim> >::iterator   InterFacesBegin();
    typename PrimitiveContainer<csmp::InterFace<dim> >::iterator   InterFacesEnd();

    // const iterators
    typename PrimitiveContainer<csmp::Node<dim> >::const_iterator  NodesBegin()      const;
    typename PrimitiveContainer<csmp::Node<dim> >::const_iterator  NodesEnd()        const;
    typename PrimitiveContainer<Element<dim> >::const_iterator     ElementsBegin()   const;
    typename PrimitiveContainer<Element<dim> >::const_iterator     ElementsEnd()     const;
    typename PrimitiveContainer<Face<dim> >::const_iterator        FacesBegin()      const;
    typename PrimitiveContainer<Face<dim> >::const_iterator        FacesEnd()        const;
    typename PrimitiveContainer<InterFace<dim> >::const_iterator   InterFacesBegin() const;
    typename PrimitiveContainer<InterFace<dim> >::const_iterator   InterFacesEnd()   const;

    const csmp::Element<dim>&   RootElement()   const;
    const csmp::Node<dim>&      RootNode()      const;
    const csmp::Face<dim>&      RootFace()      const;
    const csmp::InterFace<dim>& RootInterFace() const;

    // Accessors
    /// checks whether entitities in PrimitiveContainer were deleted or added which may have upset construction
    bool NodeAtIndexIsSafe() const;
    csmp::Node<dim>& NodeAtIndex( size_t );
    const csmp::Node<dim>& NodeAtIndex( size_t ) const;

    bool ElementAtIndexIsSafe() const;
    csmp::Element<dim>& ElementAtIndex( size_t );
    const csmp::Element<dim>& ElementAtIndex( size_t ) const;

    bool FaceAtIndexIsSafe() const;
    csmp::Face<dim>& FaceAtIndex( size_t );
    const csmp::Face<dim>& FaceAtIndex( size_t ) const;

    bool InterFaceAtIndexIsSafe() const;
    csmp::InterFace<dim>& InterFaceAtIndex( size_t );
    const csmp::InterFace<dim>& InterFaceAtIndex( size_t ) const;

    /// returns a pointer to the Node storead at index i in the PrimitiveContainer; check validy before with NodeAtIndexIsSafe()
    csmp::Node<dim>*      const PointerToNodeAtIndex( size_t ) const;
    csmp::Element<dim>*   const PointerToElementAtIndex( size_t ) const;
    csmp::Face<dim>*      const PointerToFaceAtIndex( size_t ) const;
    csmp::InterFace<dim>* const PointerToInterFaceAtIndex( size_t ) const;

    /// checks whether the node pointed to by n was destructed
    bool Verify( const Node<dim>* n ) const     { return !node_collection_.OnFreeList(n); }
    bool Verify(const Element<dim>* n ) const   { return !elmt_collection_.OnFreeList(n); }
    bool Verify(const Face<dim>* n ) const      { return !face_collection_.OnFreeList(n); }
    bool Verify(const InterFace<dim>* n ) const { return !interface_collection_.OnFreeList(n); }

    // Insertion
    /// emplaces (avoiding copying by move semantics) a temporary node into the PrimitiveContainer
    csmp::Node<dim>*  PushBack( csmp::Node<dim>&& );

    /// emplaces (avoiding copying by move semantics) temporary multiple nodes into the PrimitiveContainer
    template<typename... Args>
    csmp::Node<dim>*  EmplaceNode(Args&&... args) { node_collection_.Emplace(std::forward<Args>(args)...); }

    csmp::Element<dim>*  PushBack( csmp::Element<dim>&& );

    template<typename... Args>
    csmp::Element<dim>* EmplaceElement(Args&&... args) { elmt_collection_.Emplace(std::forward<Args>(args)...); }

    csmp::Face<dim>* const PushBack( csmp::Face<dim>&& );

    template<typename... Args>
    csmp::Face<dim>*  EmplaceFace(Args&&... args) { face_collection_.Emplace(std::forward<Args>(args)...); }

    csmp::InterFace<dim>* PushBack( csmp::InterFace<dim>&& );

    template<typename... Args>
    csmp::InterFace<dim>* EmplaceInterFace(Args&&... args) { interface_collection_.Emplace(std::forward<Args>(args)...); }
 
    /// pushes back Node object if it does not already exist in the node deque
    csmp::Node<dim>*      PushBackIfUnique( csmp::Node<dim>&& );
    csmp::Element<dim>*   PushBackIfUnique( csmp::Element<dim>&& );
    csmp::Face<dim>*      PushBackIfUnique( csmp::Face<dim>&& );
    csmp::InterFace<dim>* PushBackIfUnique( csmp::InterFace<dim>&& );

    /// removes primitives
    void Erase( csmp::Node<dim>& );
    void Erase( csmp::Element<dim>& );
    void Erase( csmp::Face<dim>& );
    void Erase( csmp::InterFace<dim>& );

    /// erase all objects of the given type
    bool EraseNodes();
    bool EraseElements();
    bool EraseFaces();
    bool EraseInterFaces();

    /// (Re)number all cells; either continuous for all cells or seperate ranges for all entity types (const because idx is mutable)
    void AssignUniqueNumbers( bool in_a_single_sequence ) const;
  
    /// outputs mesh as polygonal dataset (VSet, see HDF doc of NCSA, Urbana, Champagne, Il, US)
    void OutputMeshTo( VSet<dim>& ) const;
  
    /// relying exactly on the output from OutputMeshTo()
    void InputMeshFrom( const PropertyDatabase<dim>&, FiniteElementManager&, const VSet<dim>&, bool rebuild_mesh_from_scratch );
  
    /// storing distributed variables associated with the mesh in the VSet
    void OutputStoredVariablesTo( const PropertyDatabase<dim>&, VSet<dim>& ) const;
    void InputStoredVariablesFrom( const PropertyDatabase<dim>&, const VSet<dim>& );

    /// Rebuild parent relationships, for example after a region was removed
    void RebuildParentRelationships(typename std::vector<csmp::Node<dim>*>::iterator begin, typename std::vector<csmp::Node<dim>*>::iterator end);
  
    /// prints stored objects and their connectivity to screen
    void Out() const;
    
  private:

    /// constructs Elements (no Faces or InterFaces) and initialises their property storage
    bool BuildElementsAndVariableStorage( const PropertyDatabase<dim>&,
                                          const FiniteElementManager&, 
                                          const VSet<dim>& );

    /// 'native' method for the reconstruction of the mesh (Element-, Face- & InterFace objects) from CSMP binary file
    bool ReconstructMeshAndVariableStorage( const PropertyDatabase<dim>&,
                                            const FiniteElementManager&,
                                            const VSet<dim>& );

    /// connects the elements to their nodes etc. using information provided in VSet, except for element neighbors
    bool InitializeConnectivity( const VSet<dim>& );

    /// 'native' method that does not do any connectivity testing anymore; use if model is loaded from CSMP binary
    bool InitializeVerifiedConnectivity( const VSet<dim>& );

    /// these collections are used only if adaptive remeshing is turned off, else access is via root node or element only
    PrimitiveContainer<csmp::Node<dim> >      node_collection_;
    PrimitiveContainer<csmp::Element<dim> >   elmt_collection_;
    PrimitiveContainer<csmp::Face<dim> >      face_collection_;
    PrimitiveContainer<csmp::InterFace<dim> > interface_collection_;
    bool                                      hybrid_element_mesh_;
};

} // end namespace csmp


#endif


