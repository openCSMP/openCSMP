#ifndef ELEMENT_FACE_H
#define ELEMENT_FACE_H

#include "CSMP_definitions.h"
#include "Point.h"
#include "Box.h"

namespace csmp {

template<uint32_t> class ControlVolumeElement;

template<uint32_t dim>
class ElementFace
{
public:
    ElementFace();
    ~ElementFace();
    ElementFace(ControlVolumeElement<dim>& CVE1, 
                ControlVolumeElement<dim>& CVE2, 
                size_t local_property_storage_size, 
                BOX_BOUNDARY place);
    void                            ConstructFace(ControlVolumeElement<dim>& CVE1, 
                                                  ControlVolumeElement<dim>& CVE2, 
                                                  size_t local_property_storage_size, 
                                                  BOX_BOUNDARY place);
    BOX_BOUNDARY                    Placement();
    double                          Area();
    void                            Area(double area);
    Point<dim>                      UnitNormal();
    size_t                          ID();
    ControlVolumeElement<dim>*      InsideCVE();
    ControlVolumeElement<dim>*      OutsideCVE();
    double                          PropertyValue(const size_t& index) const;
    void                            PropertyValue(const size_t& index, const double& value);
    void                            PropertyValueAdd(const size_t& index, const double& value);
    void                            ResizePropertyStorage(const size_t& new_size);
    size_t                          PropertyStorageSize() const;
    void                            Active(bool active);
    bool                            Active() const;
    


private:
    static size_t                   elementFaceID_Counter_;
    size_t                          elementFaceID_;
    BOX_BOUNDARY                    place_;
    double                        elementFaceArea_;
    Point<dim>                      faceUnitNormal_;
    ControlVolumeElement<dim>*      insideControlVolumeElement_;
    ControlVolumeElement<dim>*      outsideControlVolumeElement_;
    std::vector<double>           localPropertyStorage_;
    bool                            active_;

};


template<uint32_t dim>
inline BOX_BOUNDARY ElementFace<dim>::Placement()
{
    return place_;
}



template<uint32_t dim>
inline double ElementFace<dim>::Area()
{
    return elementFaceArea_;
}



template<uint32_t dim>
inline void ElementFace<dim>::Area(double area)
{
    elementFaceArea_ = area;
}



template<uint32_t dim>
inline Point<dim> ElementFace<dim>::UnitNormal()
{
    return faceUnitNormal_;
}



template<uint32_t dim>
inline size_t ElementFace<dim>::ID()
{
    return elementFaceID_;
}



template<uint32_t dim>
inline ControlVolumeElement<dim>* ElementFace<dim>::InsideCVE()
{
    return insideControlVolumeElement_;
}



template<uint32_t dim>
inline ControlVolumeElement<dim>* ElementFace<dim>::OutsideCVE()
{
    return outsideControlVolumeElement_;
}



template<uint32_t dim>
inline double ElementFace<dim>::PropertyValue(const size_t& index) const
{
    return localPropertyStorage_[index];
}



template<uint32_t dim>
inline void ElementFace<dim>::PropertyValue(const size_t& index, const double& value)
{
    localPropertyStorage_[index] = value;
}



template<uint32_t dim>
inline void ElementFace<dim>::PropertyValueAdd(const size_t& index, const double& value)
{
    localPropertyStorage_[index] += value;
}



template<uint32_t dim>
inline void ElementFace<dim>::ResizePropertyStorage(const size_t& new_size)
{
    localPropertyStorage_.resize(new_size);
    std::vector<double>(localPropertyStorage_).swap(localPropertyStorage_);
}



template<uint32_t dim>
inline size_t ElementFace<dim>::PropertyStorageSize() const
{
    return localPropertyStorage_.size();
}



template<uint32_t dim>
inline void ElementFace<dim>::Active(bool active)
{
    active_ = active;
}



template<uint32_t dim>
inline bool ElementFace<dim>::Active() const
{
    return active_;
}

} // end csmp

#endif
