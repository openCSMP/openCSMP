#ifndef CONTROL_VOLUME_ELEMENT_H
#define CONTROL_VOLUME_ELEMENT_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {


template<uint32_t> class Element;
template<uint32_t> class Node;
template<uint32_t> class ElementFace;
template<uint32_t> class CVE_Patch;


template<uint32_t dim>
class ControlVolumeElement
{
public:
    ControlVolumeElement();
    ~ControlVolumeElement();
    ControlVolumeElement(Element<dim>* element, size_t no_faces, size_t local_property_storage_size);
    ControlVolumeElement(Node<dim>* node, size_t no_faces, size_t local_property_storage_size); // constructor for zero dimensional elements
    void                                        ConstructControlVolumeElement(Element<dim>* element, size_t no_faces, size_t local_property_storage_size);
    void                                        ConstructControlVolumeElement(Node<dim>* node, size_t no_faces, size_t local_property_storage_size); // method for zero dimensional elements
    void                                        Face(ElementFace<dim>* element_face, size_t face_number, char direction);
    ElementFace<dim>*                           Face(size_t face_number) const;
    double                                      FaceNormalDirection(size_t face_number) const;
    size_t                                      Neighbors() const;
    ControlVolumeElement<dim>*                  Neighbor(size_t neighbor_num) const;
    Element<dim>*                               E();
    Node<dim>*                                  N(size_t node_number);
    size_t                                      Nodes() const;
    Point<dim>                                  BaryCenter() const;
    double                                      Volume() const;
    void                                        Volume(double vol);
    size_t                                      Faces() const;
    size_t                                      ID() const;
    size_t                                      Dimension();
    double                                      PropertyValue(const size_t& index) const;
    void                                        PropertyValue(const size_t& index, const double& value);
    void                                        PropertyValueAdd(const size_t& index, const double& value);
    void                                        ResizePropertyStorage(const size_t& new_size);
    size_t                                      PropertyStorageSize() const;
    void                                        Active(bool active);
    bool                                        Active() const;

    void                                        ReservePatchStorage(size_t size);
    void                                        PushbackPatch(CVE_Patch<dim>* patch);
    size_t                                      Patches();
    CVE_Patch<dim>*                             Patch(size_t patch_number);
    typename std::vector<CVE_Patch<dim>*>::iterator      PatchesBegin();
    typename std::vector<CVE_Patch<dim>*>::iterator      PatchesEnd();
    
private:
    static size_t                    controlVolumeElementID_Counter_;
    size_t                           controlVolumeElementID_;
    Element<dim>*                    elementPtr_;
    Node<dim>*                       nodePtr_; // this is used for zero dimensional elements. in 1D and 2D it is zero
    std::vector<ElementFace<dim>*>   controlVolumeElementFaces_;
    std::vector<int8_t>              controlVolumeElementFacesNormalDirection_;
    std::vector<double>              localPropertyStorage_;
    double                           volume_;
    unsigned char                    elementDimension_;
    bool                             active_;
    std::vector<CVE_Patch<dim>*>     parentPatches_;
    

};




template<uint32_t dim>
inline void ControlVolumeElement<dim>::Face(ElementFace<dim>* element_face, size_t face_number, char direction)
{
    controlVolumeElementFaces_[face_number] = element_face;
    controlVolumeElementFacesNormalDirection_[face_number] = direction;
}



template<uint32_t dim>
inline ElementFace<dim>* ControlVolumeElement<dim>::Face(size_t face_number) const
{
    return controlVolumeElementFaces_[face_number];
}



template<uint32_t dim>
inline double ControlVolumeElement<dim>::FaceNormalDirection(size_t face_number) const
{
    return static_cast<double>(controlVolumeElementFacesNormalDirection_[face_number]);
}



template<uint32_t dim>
inline size_t ControlVolumeElement<dim>::Neighbors() const
{
    return controlVolumeElementFaces_.size();
}



template<uint32_t dim>
inline ControlVolumeElement<dim>* ControlVolumeElement<dim>::Neighbor(size_t neighbor_num) const
{
    if (controlVolumeElementFacesNormalDirection_[neighbor_num] == 1)
        return controlVolumeElementFaces_[neighbor_num]->OutsideCVE();
    else
        return controlVolumeElementFaces_[neighbor_num]->InsideCVE();
}



template<uint32_t dim>
inline Element<dim>* ControlVolumeElement<dim>::E()
{
    return elementPtr_;
}



template<uint32_t dim>
inline Node<dim>* ControlVolumeElement<dim>::N(size_t node_number)
{
    if (elementDimension_ == 0)
        return nodePtr_;
    else
        return elementPtr_->N(node_number);
}



template<uint32_t dim>
inline size_t ControlVolumeElement<dim>::Nodes() const
{
    if (elementDimension_ == 0)
        return 1U;
    else
        return elementPtr_->Nodes();

}



template<uint32_t dim>
inline Point<dim> ControlVolumeElement<dim>::BaryCenter() const
{
    if (elementDimension_ == 0)
        return nodePtr_->Coordinate();
    else
        return elementPtr_->BaryCenter();

}



template<uint32_t dim>
inline double ControlVolumeElement<dim>::Volume() const
{
    return volume_;
}



template<uint32_t dim>
inline void ControlVolumeElement<dim>::Volume(double vol)
{
    volume_ = vol;
}



template<uint32_t dim>
inline size_t ControlVolumeElement<dim>::Faces() const
{
    return controlVolumeElementFaces_.size();
}



template<uint32_t dim>
inline size_t ControlVolumeElement<dim>::ID() const
{
    return controlVolumeElementID_;
}



template<uint32_t dim>
inline size_t ControlVolumeElement<dim>::Dimension()
{
    return (size_t) elementDimension_;
}




template<uint32_t dim>
inline double ControlVolumeElement<dim>::PropertyValue(const size_t& index) const
{
    return localPropertyStorage_[index];
}



template<uint32_t dim>
inline void ControlVolumeElement<dim>::PropertyValue(const size_t& index, const double& value)
{
    localPropertyStorage_[index] = value;
}



template<uint32_t dim>
inline void ControlVolumeElement<dim>::PropertyValueAdd(const size_t& index, const double& value)
{
    localPropertyStorage_[index] += value;
}



template<uint32_t dim>
inline void ControlVolumeElement<dim>::ResizePropertyStorage(const size_t& new_size)
{
    localPropertyStorage_.resize(new_size);
    std::vector<double>(localPropertyStorage_).swap(localPropertyStorage_);
}



template<uint32_t dim>
inline size_t ControlVolumeElement<dim>::PropertyStorageSize() const
{
    return localPropertyStorage_.size();
}


template<uint32_t dim>
inline void ControlVolumeElement<dim>::Active(bool active)
{
    active_ = active;
}


template<uint32_t dim>
inline bool ControlVolumeElement<dim>::Active() const
{
    return active_;
}


template<uint32_t dim>
inline void ControlVolumeElement<dim>::ReservePatchStorage(size_t size)
{
    parentPatches_.reserve(size);
}



template<uint32_t dim>
inline void ControlVolumeElement<dim>::PushbackPatch(CVE_Patch<dim>* patch)
{
    parentPatches_.push_back(patch);
}



template<uint32_t dim>
inline size_t ControlVolumeElement<dim>::Patches()
{
    return parentPatches_.size();
}



template<uint32_t dim>
inline CVE_Patch<dim>* ControlVolumeElement<dim>::Patch(size_t patch_number)
{
    return parentPatches_[patch_number];
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim>*>::iterator ControlVolumeElement<dim>::PatchesBegin()
{
    return parentPatches_.begin();
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim>*>::iterator ControlVolumeElement<dim>::PatchesEnd()
{
    return parentPatches_.end();
}


} // end csmp

#endif
