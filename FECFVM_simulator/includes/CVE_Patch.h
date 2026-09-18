// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CVE_PATCH_H
#define CVE_PATCH_H

#include "CVE_Sector.h"

namespace csmp {


template<uint32_t> class Node;
template<uint32_t> class ControlVolumeElement;
template<uint32_t> class ElementFace;

template<uint32_t dim>
class CVE_Patch {
  public:
      CVE_Patch(Node<dim>& nd);
      ~CVE_Patch();
      
      uint32_t                                          Sectors() const;
      CVE_Sector<dim>*                                  SectorAccessor(uint32_t sec_num);
      typename std::vector<CVE_Sector<dim> >::iterator  SectorsBegin();
      typename std::vector<CVE_Sector<dim> >::iterator  SectorsEnd();
      double                                            Volume();
      Node<dim>*                                        CenterNode();
      void                                              CreateAndAddSector(ControlVolumeElement<dim>& cvelm, double lower_dimensional_width);
      void                                              ReserveSectorsStorage(size_t size);
      void                                              ReserveFacesStorage(size_t size);
      void                                              PushbackFace(ElementFace<dim>* face);
      void                                              FaceCyclicDirection(size_t face, char direction);
      char                                              FaceCyclicDirection(size_t face);
      void                                              BoundaryFaces(size_t faces);
      void                                              InteriorFaces(size_t faces);
      uint32_t                                          Faces() const;
      uint32_t                                          InteriorFaces() const;
      ElementFace<dim>*                                 FaceAccessor(uint32_t face_num);
      typename std::vector<ElementFace<dim>*>::iterator FacesBegin();
      typename std::vector<ElementFace<dim>*>::iterator PerimeterFacesBegin();
      typename std::vector<ElementFace<dim>*>::iterator FacesEnd();
      void                                              Active(bool active);
      bool                                              Active() const;

  private:
      
      std::vector<CVE_Sector<dim> >                 sectorsCollection_;
      std::vector<ElementFace<dim>* >               facesCollection_;
      std::vector<char>                             facesCyclicDirection_;
      Node<dim>*                                    centerNode_ptr_;
      unsigned char                                 interiorFaces_;
      unsigned char                                 boundaryFaces_;
      bool                                          active_;
};



template<uint32_t dim>
inline void CVE_Patch<dim>::ReserveSectorsStorage(size_t size)
{
    sectorsCollection_.reserve(size);
}




template<uint32_t dim>
inline void CVE_Patch<dim>::ReserveFacesStorage(size_t size)
{
    facesCollection_.reserve(size);
    facesCyclicDirection_.reserve(size);
}




template<uint32_t dim>
inline void CVE_Patch<dim>::PushbackFace(ElementFace<dim>* face)
{
    facesCollection_.push_back(face);
    facesCyclicDirection_.push_back(0);
}




template<uint32_t dim>
inline void CVE_Patch<dim>::FaceCyclicDirection(size_t face, char direction)
{
    facesCyclicDirection_[face] = direction;
}




template<uint32_t dim>
inline char CVE_Patch<dim>::FaceCyclicDirection(size_t face)
{
    return facesCyclicDirection_[face];
}




template<uint32_t dim>
inline uint32_t CVE_Patch<dim>::Sectors() const
{
    return static_cast<uint32_t>(sectorsCollection_.size());
}




template<uint32_t dim>
inline CVE_Sector<dim>* CVE_Patch<dim>::SectorAccessor(uint32_t sec_num)
{
    return &sectorsCollection_[sec_num];
}




template<uint32_t dim>
inline typename std::vector<CVE_Sector<dim> >::iterator CVE_Patch<dim>::SectorsBegin()
{
    return sectorsCollection_.begin();
}




template<uint32_t dim>
inline typename std::vector<CVE_Sector<dim> >::iterator CVE_Patch<dim>::SectorsEnd()
{
    return sectorsCollection_.end();
}




template<uint32_t dim>
inline Node<dim>* CVE_Patch<dim>::CenterNode()
{
    return centerNode_ptr_;
}




template<uint32_t dim>
inline uint32_t CVE_Patch<dim>::Faces() const
{
    return interiorFaces_ + boundaryFaces_;
}




template<uint32_t dim>
inline uint32_t CVE_Patch<dim>::InteriorFaces() const
{
    return interiorFaces_;
}




template<uint32_t dim>
inline void CVE_Patch<dim>::BoundaryFaces(size_t faces)
{
    boundaryFaces_ = faces;
}




template<uint32_t dim>
inline void CVE_Patch<dim>::InteriorFaces(size_t faces)
{
    interiorFaces_ = faces;
}




template<uint32_t dim>
inline ElementFace<dim>* CVE_Patch<dim>::FaceAccessor(uint32_t face_num)
{
    return facesCollection_[face_num];
}




template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator CVE_Patch<dim>::FacesBegin()
{
    return facesCollection_.begin();
}




template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator CVE_Patch<dim>::PerimeterFacesBegin()
{
    return facesCollection_.begin() + (long) interiorFaces_;
}




template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator CVE_Patch<dim>::FacesEnd()
{
    return facesCollection_.end();
}



template<uint32_t dim>
inline void CVE_Patch<dim>::Active(bool active)
{
    active_ = active;
}



template<uint32_t dim>
inline bool CVE_Patch<dim>::Active() const
{
    return active_;
}


} // end namespace csmp

#endif
