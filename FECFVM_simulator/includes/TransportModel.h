// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef TRANSPORT_MODEL_H
#define TRANSPORT_MODEL_H

#include "CSMP_definitions.h"
#include "ElementFace.h"
#include "ControlVolumeElement.h"
#include "CVE_Patch.h"



namespace csmp {


template<uint32_t> class Model;


template<uint32_t dim>
class TransportModel
{

public:
    TransportModel(Model<dim>& model, const char* region_name, const char* vol_mod,
                   size_t element_local_storage_size, size_t boundary_element_local_storage_size,
                   size_t face_local_storage_size, size_t boundary_face_local_storage_size);
    ~TransportModel();

    void                                                         CreateControlVolumeElementsAndFaces(const char* region_name,
                                                                                                     size_t element_local_storage_size, 
                                                                                                     size_t boundary_element_local_storage_size,
                                                                                                     size_t face_local_storage_size, 
                                                                                                     size_t boundary_face_local_storage_size);
    void                                                         CreatePatches(const char* vol_mod);
    void                                                         CorrectVolumes(const char* vol_mod);

    size_t                                                       CVEs() const;
    size_t                                                       InteriorCVEs() const;
    size_t                                                       SurfaceCVEs() const;
    size_t                                                       InteriorSurfaceCVEs() const;
    size_t                                                       LineCVEs() const;
    size_t                                                       InteriorLineCVEs() const;
    size_t                                                       PointCVEs() const;
    size_t                                                       InteriorPointCVEs() const;
    ControlVolumeElement<dim>*                                   CVE_Accessor(size_t CVE_num);

    size_t                                                       Faces() const;
    size_t                                                       InteriorFaces() const;
    size_t                                                       LineFaces() const;
    size_t                                                       InteriorLineFaces() const;
    size_t                                                       PointFaces() const;
    size_t                                                       InteriorPointFaces() const;
    ElementFace<dim>*                                            FaceAccessor(size_t face_num);

    size_t                                                       Patches() const;
    size_t                                                       InteriorPatches() const;
    CVE_Patch<dim>*                                              PatchAccessor(size_t patch_num);
    void                                                         AssignParentPatchesToCVEs();


    typename std::vector<ControlVolumeElement<dim> >::iterator   CVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   PerimeterCVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   CVEsEnd();
    typename std::vector<ControlVolumeElement<dim> >::iterator   SurfaceCVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   PerimeterSurfaceCVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   SurfaceCVEsEnd();
    typename std::vector<ControlVolumeElement<dim> >::iterator   LineCVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   PerimeterLineCVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   LineCVEsEnd();
    typename std::vector<ControlVolumeElement<dim> >::iterator   PointCVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   PerimeterPointCVEsBegin();
    typename std::vector<ControlVolumeElement<dim> >::iterator   PointCVEsEnd();
    
    typename std::vector<ElementFace<dim> >::iterator            FacesBegin();
    typename std::vector<ElementFace<dim> >::iterator            PerimeterFacesBegin();
    typename std::vector<ElementFace<dim> >::iterator            FacesEnd();
    typename std::vector<ElementFace<dim> >::iterator            LineFacesBegin();
    typename std::vector<ElementFace<dim> >::iterator            PerimeterLineFacesBegin();
    typename std::vector<ElementFace<dim> >::iterator            LineFacesEnd();
    typename std::vector<ElementFace<dim> >::iterator            PointFacesBegin();
    typename std::vector<ElementFace<dim> >::iterator            PerimeterPointFacesBegin();
    typename std::vector<ElementFace<dim> >::iterator            PointFacesEnd();

    typename std::vector<CVE_Patch<dim> >::iterator              PatchesBegin();
    typename std::vector<CVE_Patch<dim> >::iterator              PerimeterPatchesBegin();
    typename std::vector<CVE_Patch<dim> >::iterator              PatchesEnd();

    void                                                         ZeroCVEsProperty(size_t prop);
    void                                                         ZeroFacesProperty(size_t prop);
    
    void                                                         TransportModelReport();
    


private:
    
    Model<dim>&                                           modelRef_;
    const std::string                                     regionName_;
    // stores CVEs in the following order:
    // INTERIOR SURFACE CVE >> INTERIOR LINE CVE >> INTERIOR POINT CVE >>
    // BOUNDARY SURFACE CVE >> BOUNDARY LINE CVE >> BOUNDARY POINT CVE
    std::vector<ControlVolumeElement<dim> >               CVE_Collection_;
    // stores faces in the following order:
    // INTERIOR LINE FACE >> INTERIOR POINT FACE >> 
    // BOUNDARY LINE FACE >> BOUNDARY POINT FACE >>
    std::vector<ElementFace<dim> >                        faceCollection_;
    // stores patches. first comes interior patches and at the end of the
    // vector boundary patches
    std::vector<CVE_Patch<dim> >                          patchCollection_;
    size_t                                                interiorSurfaceCVE_;
    size_t                                                boundarySurfaceCVE_;
    size_t                                                interiorLineCVE_;
    size_t                                                boundaryLineCVE_;
    size_t                                                interiorPointCVE_;
    size_t                                                boundaryPointCVE_;
    size_t                                                interiorLineFaces_;
    size_t                                                boundaryLineFaces_;
    size_t                                                interiorPointFaces_;
    size_t                                                boundaryPointFaces_;
    size_t                                                interiorPatches_;
    size_t                                                boundaryPatches_;

};



template<uint32_t dim>
inline size_t TransportModel<dim>::CVEs() const
{
    return (interiorSurfaceCVE_ + boundarySurfaceCVE_ + interiorLineCVE_ + boundaryLineCVE_ + 
            interiorPointCVE_ + boundaryPointCVE_);
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorCVEs() const
{
    return (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_);
}



template<uint32_t dim>
inline size_t TransportModel<dim>::SurfaceCVEs() const
{
    return (interiorSurfaceCVE_ + boundarySurfaceCVE_);
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorSurfaceCVEs() const
{
    return interiorSurfaceCVE_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::LineCVEs() const
{
    return (interiorLineCVE_ + boundaryLineCVE_);
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorLineCVEs() const
{
    return interiorLineCVE_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::PointCVEs() const
{
    return (interiorPointCVE_ + boundaryPointCVE_);
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorPointCVEs() const
{
    return interiorPointCVE_;
}



template<uint32_t dim>
inline ControlVolumeElement<dim>* TransportModel<dim>::CVE_Accessor(size_t CVE_num)
{
    return &CVE_Collection_[CVE_num];
}



template<uint32_t dim>
inline size_t TransportModel<dim>::Faces() const
{
    return interiorLineFaces_ + interiorPointFaces_ + boundaryLineFaces_ + boundaryPointFaces_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorFaces() const
{
    return interiorLineFaces_ + interiorPointFaces_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::LineFaces() const
{
    return interiorLineFaces_ + boundaryLineFaces_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorLineFaces() const
{
    return interiorLineFaces_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::PointFaces() const
{
    return interiorPointFaces_ + boundaryPointFaces_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorPointFaces() const
{
    return interiorPointFaces_;
}



template<uint32_t dim>
inline ElementFace<dim>* TransportModel<dim>::FaceAccessor(size_t face_num)
{
    return &faceCollection_[face_num];
}



template<uint32_t dim>
inline size_t TransportModel<dim>::Patches() const
{
    return interiorPatches_ + boundaryPatches_;
}



template<uint32_t dim>
inline size_t TransportModel<dim>::InteriorPatches() const
{
    return interiorPatches_;
}



template<uint32_t dim>
inline CVE_Patch<dim>* TransportModel<dim>::PatchAccessor(size_t patch_num)
{
    return &patchCollection_[patch_num];
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::CVEsBegin()
{
    return CVE_Collection_.begin();
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::PerimeterCVEsBegin()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_);
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::CVEsEnd()
{
    return CVE_Collection_.end();
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::SurfaceCVEsBegin()
{
    return CVE_Collection_.begin();
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::PerimeterSurfaceCVEsBegin()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_);
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::SurfaceCVEsEnd()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_ + 
                                             boundarySurfaceCVE_);
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::LineCVEsBegin()
{
    return CVE_Collection_.begin() + (long) interiorSurfaceCVE_;
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::PerimeterLineCVEsBegin()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_ + 
                                             boundarySurfaceCVE_);
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::LineCVEsEnd()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_ + 
                                             boundarySurfaceCVE_ + boundaryLineCVE_);
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::PointCVEsBegin()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_);
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::PerimeterPointCVEsBegin()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_ + 
                                             boundarySurfaceCVE_ + boundaryLineCVE_);
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim> >::iterator TransportModel<dim>::PointCVEsEnd()
{
    return CVE_Collection_.begin() + (long) (interiorSurfaceCVE_ + interiorLineCVE_ + interiorPointCVE_ + 
                                             boundarySurfaceCVE_ + boundaryLineCVE_ + boundaryPointCVE_);
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::FacesBegin()
{
    return faceCollection_.begin();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::PerimeterFacesBegin()
{
    return faceCollection_.begin() + (long) (interiorLineFaces_ + interiorPointFaces_);
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::FacesEnd()
{
    return faceCollection_.end();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::LineFacesBegin()
{
    return faceCollection_.begin();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::PerimeterLineFacesBegin()
{
    return faceCollection_.begin() + (long) (interiorLineFaces_ + interiorPointFaces_);
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::LineFacesEnd()
{
    return faceCollection_.begin() + (long) (interiorLineFaces_ + interiorPointFaces_ + 
                                             boundaryLineFaces_);
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::PointFacesBegin()
{
    return faceCollection_.begin() + (long) interiorLineFaces_;
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::PerimeterPointFacesBegin()
{
    return faceCollection_.begin() + (long) (interiorLineFaces_ + interiorPointFaces_ + 
                                             boundaryLineFaces_);
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim> >::iterator TransportModel<dim>::PointFacesEnd()
{
    return faceCollection_.end();
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim> >::iterator TransportModel<dim>::PatchesBegin()
{
    return patchCollection_.begin();
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim> >::iterator TransportModel<dim>::PerimeterPatchesBegin()
{
    return patchCollection_.begin() + (long) interiorPatches_;
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim> >::iterator TransportModel<dim>::PatchesEnd()
{
    return patchCollection_.end();
}



template<uint32_t dim>
inline void TransportModel<dim>::ZeroCVEsProperty(size_t prop)
{
    for (typename std::vector<ControlVolumeElement<dim> >::iterator eit = CVEsBegin(); eit != CVEsEnd(); eit++)
        eit->PropertyValue(prop, 0.);
}



template<uint32_t dim>
inline void TransportModel<dim>::ZeroFacesProperty(size_t prop)
{
    for (typename std::vector<ElementFace<dim> >::iterator fit = FacesBegin(); fit != FacesEnd(); fit++)
        fit->PropertyValue(prop, 0.);
}



} // end csmp

#endif
