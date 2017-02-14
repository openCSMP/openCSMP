#ifndef CSMP_MESH_MANAGER_H
#define CSMP_MESH_MANAGER_H

#include <deque>
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"

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
    MeshManager( const MeshManager& );
    MeshManager&  operator=( const MeshManager& );
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
    typename std::deque<csmp::Node<dim> >::iterator        NodesBegin();
    typename std::deque<csmp::Node<dim> >::iterator        NodesEnd();
    typename std::deque<csmp::Element<dim> >::iterator     ElementsBegin();
    typename std::deque<csmp::Element<dim> >::iterator     ElementsEnd();
    typename std::deque<csmp::Face<dim> >::iterator        FacesBegin();
    typename std::deque<csmp::Face<dim> >::iterator        FacesEnd();
    typename std::deque<csmp::InterFace<dim> >::iterator   InterFacesBegin();
    typename std::deque<csmp::InterFace<dim> >::iterator   InterFacesEnd();

    // const iterators
    typename std::deque<csmp::Node<dim> >::const_iterator  NodesBegin()      const;
    typename std::deque<csmp::Node<dim> >::const_iterator  NodesEnd()        const;
    typename std::deque<Element<dim> >::const_iterator     ElementsBegin()   const;
    typename std::deque<Element<dim> >::const_iterator     ElementsEnd()     const;
    typename std::deque<Face<dim> >::const_iterator        FacesBegin()      const;
    typename std::deque<Face<dim> >::const_iterator        FacesEnd()        const;
    typename std::deque<InterFace<dim> >::const_iterator   InterFacesBegin() const;
    typename std::deque<InterFace<dim> >::const_iterator   InterFacesEnd()   const;

    const csmp::Element<dim>&   RootElement()   const;
    const csmp::Node<dim>&      RootNode()      const;
    const csmp::Face<dim>&      RootFace()      const;
    const csmp::InterFace<dim>& RootInterFace() const;

    // Insertion
    /// emplaces (no copying) a node at the end of the deque in which the Node objects may be stored (when adaptive_remeshing_ false)
    csmp::Node<dim>*      PushBack( csmp::Node<dim>&& );

    /// emplaces (no copying) an element at the end of the deque in which the Element objects may be stored (when adaptive_remeshing_ false)
    csmp::Element<dim>*   PushBack( csmp::Element<dim>&& );

    /// emplaces (no copying) a face at the end of the deque in which the Face objects may be stored (when adaptive_remeshing_ false)
    csmp::Face<dim>* const PushBack( csmp::Face<dim>&& );

    /// emplaces (no copying) an interface at the end of the deque in which the InterFace objects may be stored (when adaptive_remeshing_ false)
    csmp::InterFace<dim>* PushBack( csmp::InterFace<dim>&& );

    // TODO: SKM: check these (which do not make much sense for deque containers); remove what is not needed
    /// pushes back Node object if it does not already exist in the node deque
    csmp::Node<dim>*      PushBackIfUnique( csmp::Node<dim>&& );
    csmp::Element<dim>*   PushBackIfUnique( csmp::Element<dim>&& );
    csmp::Face<dim>*      PushBackIfUnique( csmp::Face<dim>&& );
    csmp::InterFace<dim>* PushBackIfUnique( csmp::InterFace<dim>&& );

    /// removes node without invalidating any other pointers, references or values of the container unless the last element is deleted (return = false)
    bool Erase( const csmp::Node<dim>& );
    bool Erase( const csmp::Element<dim>& );
    bool Erase( const csmp::Face<dim>& );
    bool Erase( const csmp::InterFace<dim>& );

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

    /// prints stored objects and their connectivity to a stream
    void Out() const { Out(std::cout); }
    void Out(std::ostream& os) const;
    
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
  
  private:

    /// these collections are used only if adaptive remeshing is turned off, else access is via root node or element only
    bool                              adaptive_remeshing_;
    bool                              hybrid_element_mesh_;
    std::deque<csmp::Node<dim> >      node_collection_;
    std::deque<csmp::Element<dim> >   elmt_collection_;
    std::deque<csmp::Face<dim> >      face_collection_;
    std::deque<csmp::InterFace<dim> > interface_collection_;
};


template<size_t dim>
inline bool  MeshManager<dim>::HybridElementMesh() const 
 { return hybrid_element_mesh_; }

// size of containers

template<size_t dim>
inline size_t  MeshManager<dim>::Nodes() const
 { assert( !adaptive_remeshing_ ); return node_collection_.size(); }

template<size_t dim>
inline size_t  MeshManager<dim>::Elements() const
 { assert( !adaptive_remeshing_ ); return elmt_collection_.size(); }

template<size_t dim>
inline size_t  MeshManager<dim>::Faces() const
 { assert( !adaptive_remeshing_ ); return face_collection_.size(); }

template<size_t dim>
inline size_t  MeshManager<dim>::InterFaces() const
 { assert( !adaptive_remeshing_ ); return interface_collection_.size(); }


// accessors

template<size_t dim>
inline csmp::Node<dim>&    MeshManager<dim>::RootNode()
 { assert( !node_collection_.empty() ); return node_collection_[0]; }

template<size_t dim>
inline csmp::Element<dim>&  MeshManager<dim>::RootElement()
 { assert( !elmt_collection_.empty() ); return elmt_collection_[0]; }

template<size_t dim>
inline csmp::Face<dim>&  MeshManager<dim>::RootFace()
 { assert( !face_collection_.empty() ); return face_collection_[0]; }

template<size_t dim>
inline csmp::InterFace<dim>&  MeshManager<dim>::RootInterFace()
 { assert( !interface_collection_.empty() ); return interface_collection_[0]; }

template<size_t dim>
inline typename std::deque<csmp::Node<dim> >::iterator  MeshManager<dim>::NodesBegin()
 { assert( !node_collection_.empty() ); return node_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::Node<dim> >::iterator  MeshManager<dim>::NodesEnd()
 { assert( !node_collection_.empty() ); return node_collection_.end(); }

template<size_t dim>
inline typename std::deque<csmp::Element<dim> >::iterator  MeshManager<dim>::ElementsBegin()
 { assert( !elmt_collection_.empty() ); return elmt_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::Element<dim> >::iterator  MeshManager<dim>::ElementsEnd()
 { assert( !elmt_collection_.empty() ); return elmt_collection_.end(); }

template<size_t dim>
inline typename std::deque<csmp::Face<dim> >::iterator  MeshManager<dim>::FacesBegin()
 { assert( !face_collection_.empty() ); return face_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::Face<dim> >::iterator  MeshManager<dim>::FacesEnd()
 { assert( !face_collection_.empty() ); return face_collection_.end(); }

template<size_t dim>
inline typename std::deque<csmp::InterFace<dim> >::iterator  MeshManager<dim>::InterFacesBegin()
 { assert( !interface_collection_.empty() ); return interface_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::InterFace<dim> >::iterator  MeshManager<dim>::InterFacesEnd()
 { assert( !interface_collection_.empty() ); return interface_collection_.end(); }

// const accessors

template<size_t dim>
inline const csmp::Node<dim>&    MeshManager<dim>::RootNode() const
 { assert( !node_collection_.empty() ); return node_collection_[0]; }

template<size_t dim>
inline const csmp::Element<dim>&  MeshManager<dim>::RootElement() const
 { assert( !elmt_collection_.empty() ); return elmt_collection_[0]; }

template<size_t dim>
inline const csmp::Face<dim>&  MeshManager<dim>::RootFace() const
 { assert( !face_collection_.empty() ); return face_collection_[0]; }

template<size_t dim>
inline const csmp::InterFace<dim>&  MeshManager<dim>::RootInterFace() const
 { assert( !interface_collection_.empty() ); return interface_collection_[0]; }

template<size_t dim>
inline typename std::deque<csmp::Node<dim> >::const_iterator  MeshManager<dim>::NodesBegin() const
 { assert( !node_collection_.empty() ); return node_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::Node<dim> >::const_iterator  MeshManager<dim>::NodesEnd() const
 { assert( !node_collection_.empty() ); return node_collection_.end(); }

template<size_t dim>
inline typename std::deque<csmp::Element<dim> >::const_iterator  MeshManager<dim>::ElementsBegin() const
 { assert( !elmt_collection_.empty() ); return elmt_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::Element<dim> >::const_iterator  MeshManager<dim>::ElementsEnd() const
 { assert( !elmt_collection_.empty() ); return elmt_collection_.end(); }

template<size_t dim>
inline typename std::deque<csmp::Face<dim> >::const_iterator  MeshManager<dim>::FacesBegin() const
 { assert( !face_collection_.empty() ); return face_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::Face<dim> >::const_iterator  MeshManager<dim>::FacesEnd() const
 { assert( !face_collection_.empty() ); return face_collection_.end(); }

template<size_t dim>
inline typename std::deque<csmp::InterFace<dim> >::const_iterator  MeshManager<dim>::InterFacesBegin() const
 { assert( !interface_collection_.empty() ); return interface_collection_.begin(); }

template<size_t dim>
inline typename std::deque<csmp::InterFace<dim> >::const_iterator  MeshManager<dim>::InterFacesEnd() const
 { assert( !interface_collection_.empty() ); return interface_collection_.end(); }


} // end namespace csmp


#endif


