// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVE_SECTOR_H
#define CVE_SECTOR_H

#include "CSMP_definitions.h"
#include "Point.h"

namespace csmp {

template<uint32_t> class Element;
template<uint32_t> class Node;
template<uint32_t> class ElementFace;
template<uint32_t> class ControlVolumeElement;

template<uint32_t dim>
class CVE_Sector {

public:
    CVE_Sector(ControlVolumeElement<dim>& cvelm, Node<dim>& nd, double lower_dimensional_width);

    ~CVE_Sector();
    void                        GenerateSector(ControlVolumeElement<dim>& cvelm, Node<dim>& nd, double lower_dimensional_width);
    size_t                      Faces() const;
    double                      FaceArea(size_t face_num) const;
    Point<dim>                  FaceNormal(size_t face_num) const;
    double                      FaceNormalDirection(size_t face_num) const;
    double                      Volume() const;
    size_t                      Dimension() const;
    Element<dim>*               E() const;
    ControlVolumeElement<dim>*  CVE() const;
    ElementFace<dim>*           Face(size_t face_num) const;
    unsigned char               PatchFaceNumber(size_t face_num);
    void                        PatchFaceNumber(size_t face_num, size_t patch_face_num);
    const Point<dim>&           TotalFacetsAreaVector() const;

private:  
    ControlVolumeElement<dim>*  CVE_ptr_;
    std::vector<double>         faceArea_;
    std::vector<unsigned char>  CVE_SectorFaceConnection_;
    std::vector<unsigned char>  patchFaceNumber_;
    double                      volume_;
    Point<dim>                  totalFacetsAreaVector_;

};



template<uint32_t dim>
inline size_t CVE_Sector<dim>::Faces() const
{
    return CVE_SectorFaceConnection_.size();
}




template<uint32_t dim>
inline double CVE_Sector<dim>::FaceArea(size_t face_num) const
{
    return faceArea_[face_num];
}




template<uint32_t dim>
inline Point<dim> CVE_Sector<dim>::FaceNormal(size_t face_num) const
{
    return CVE_ptr_->Face(CVE_SectorFaceConnection_[face_num])->UnitNormal();
}




template<uint32_t dim>
inline double CVE_Sector<dim>::FaceNormalDirection(size_t face_num) const
{
    return static_cast<double>(CVE_ptr_->FaceNormalDirection(CVE_SectorFaceConnection_[face_num]));
}




template<uint32_t dim>
inline double CVE_Sector<dim>::Volume() const
{
    return volume_;
}




template<uint32_t dim>
inline size_t CVE_Sector<dim>::Dimension() const
{
    return CVE_ptr_->Dimension();;
}




template<uint32_t dim>
inline Element<dim>* CVE_Sector<dim>::E() const
{
    return CVE_ptr_->E();
}




template<uint32_t dim>
inline ControlVolumeElement<dim>* CVE_Sector<dim>::CVE() const
{
    return CVE_ptr_;
}




template<uint32_t dim>
inline ElementFace<dim>* CVE_Sector<dim>::Face(size_t face_num) const
{
    return CVE_ptr_->Face(CVE_SectorFaceConnection_[face_num]);
}




template<uint32_t dim>
inline unsigned char CVE_Sector<dim>::PatchFaceNumber(size_t face_num)
{
    return patchFaceNumber_[face_num];
}




template<uint32_t dim>
inline void CVE_Sector<dim>::PatchFaceNumber(size_t face_num, size_t patch_face_num)
{
    patchFaceNumber_[face_num] = (unsigned char) patch_face_num;
}




template<uint32_t dim>
inline const Point<dim>& CVE_Sector<dim>::TotalFacetsAreaVector() const
{
    return totalFacetsAreaVector_;
}

} // end namespace csmp

#endif
