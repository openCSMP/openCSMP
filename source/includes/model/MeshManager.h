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
    // IMPORTANT: These methods should only be used during mesh construction fom a file.
  
    /// checks whether there are any nodes are on the freelist, i.e. have been deleted, if so index operations are not safe
    bool NodeAtIndexIsSafe() const;
  
    /// retrieve node by index; call NodeAtIndexIsSafe() to see whether this operation is safe first
    csmp::Node<dim>& NodeAtIndex(size_t i);
    /// retrieve node by index; call NodeAtIndexIsSafe() to see whether this operation is safe first
    const csmp::Node<dim>& NodeAtIndex(size_t i) const;

    /// checks whether there are any elements are on the freelist, i.e. have been deleted, if so index operations are not safe
    bool ElementAtIndexIsSafe() const;

    /// retrieve element by index; call NodeAtIndexIsSafe() to see whether this operation is safe first
    csmp::Element<dim>& ElementAtIndex(size_t i);
    /// retrieve element by index; call NodeAtIndexIsSafe() to see whether this operation is safe first
    const csmp::Element<dim>& ElementAtIndex(size_t i) const;

    csmp::Face<dim>& FaceAtIndex(size_t i);
    const csmp::Face<dim>& FaceAtIndex(size_t i) const;

    /// checks whether the given Node is on the freelist, if so it has been destructed and cannot be accessed anymore
    bool VerifyNode(const Node<dim>* n) const { return !node_collection_.OnFreeList(n); }
    /// checks whether the given Element is on the freelist, if so it has been destructed and cannot be accessed anymore
    bool VerifyElement(const Element<dim>* n) const { return !elmt_collection_.OnFreeList(n); }

    // Insertion methods
  
    /// emplaces (no copying) a node at the end of the deque in which the Node objects may be stored (when adaptive_remeshing_ false)
    csmp::Node<dim>*      PushBack( csmp::Node<dim>&& );

    /// emplaces (no copying) an element at the end of the deque in which the Element objects may be stored (when adaptive_remeshing_ false)
    csmp::Element<dim>*   PushBack( csmp::Element<dim>&& );

    /// emplaces (no copying) a face at the end of the deque in which the Face objects may be stored (when adaptive_remeshing_ false)
    csmp::Face<dim>* const PushBack( csmp::Face<dim>&& );

    /// emplaces (no copying) an interface at the end of the deque in which the InterFace objects may be stored (when adaptive_remeshing_ false)
    csmp::InterFace<dim>* PushBack( csmp::InterFace<dim>&& );

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

    /// these collections store the actual FE primitives: Nodes, Elements etc. and be used to access them, else access is via root node or element
    PrimitiveContainer<csmp::Node<dim> >      node_collection_;
    PrimitiveContainer<csmp::Element<dim> >   elmt_collection_;
    PrimitiveContainer<csmp::Face<dim> >      face_collection_;
    PrimitiveContainer<csmp::InterFace<dim> > interface_collection_;
    bool                                      hybrid_element_mesh_;   ///< true if the mesh consists of different FE types
};

} // end namespace csmp









/* ANDREW BROMAGE - code created by him but is not used and does not compile

    /// emplaces (no copying) a node at the end of the deque in which the Node objects may be stored (when adaptive_remeshing_ false)
    template<typename... Args>
    csmp::Node<dim>* EmplaceNode(Args&&... args)
    {
        node_collection_.Emplace(std::forward<Args>(args)...);
    }

    /// emplaces (no copying) an element at the end of the deque in which the Element objects may be stored (when adaptive_remeshing_ false)
    template<typename... Args>
    csmp::Element<dim>* EmplaceElement(Args&&... args)
    {
        elmt_collection_.Emplace(std::forward<Args>(args)...);
    }

    /// emplaces (no copying) a face at the end of the deque in which the Face objects may be stored (when adaptive_remeshing_ false)
    template<typename... Args>
    csmp::Face<dim>* EmplaceFace(Args&&... args)
    {
        face_collection_.Emplace(std::forward<Args>(args)...);
    }

    /// emplaces (no copying) an interface at the end of the deque in which the InterFace objects may be stored (when adaptive_remeshing_ false)
    template<typename... Args>
    csmp::InterFace<dim>* EmplaceInterFace(Args&&... args)
    {
        interface_collection_.Emplace(std::forward<Args>(args)...);
    }

*/



#endif

