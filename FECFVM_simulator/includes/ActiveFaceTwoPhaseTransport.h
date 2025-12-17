#ifndef ACTIVE_FACE_TWO_PHASE_TRANSPORT
#define ACTIVE_FACE_TWO_PHASE_TRANSPORT

#include "CSMP_definitions.h"
#include "TwoPhaseElementBasedTransport.h"
#include "LU_Solver.h"
#include "VectorVariable.h"

namespace csmp {

    template<uint32_t> class CVE_Patch;

    template<uint32_t dim>
    class ActiveFaceTwoPhaseTransport : public TwoPhaseElementBasedTransport<dim>
    {
    public:
        ActiveFaceTwoPhaseTransport(Model<dim>& model,
                                    TransportModel<dim>& transport_model,
                                    const char* region_name,
                                    const char* porosity,
                                    const char* permeability,
                                    const char* saturation_oil,
                                    const char* saturation_water, // it needs this key to update it
                                    const char* divergence,       // class uses this storage to store the divergence
                                    const char* total_velocity,
                                    const char* element_oil_injection_volume_rate, // only positive values are allowed
                                    const char* element_water_injection_volume_rate, // only positive values are allowed
                                    const char* element_total_production_volume_rate, // only negative values are allowed
                                    IMPES_Setup<dim>& setup);

        ~ActiveFaceTwoPhaseTransport();

        virtual double                      Transport(double total_time_increment, TwoPhaseModel<dim>& flow_model);

    private:
        void                                  FindActiveRegions(TwoPhaseModel<dim>& flow_model, const double activation_criteria);
        double                              UpdateActiveRegions(TwoPhaseModel<dim>& flow_model, double criteria); // returns CFL based minimum time step size for new cves
        void                                  CheckForInactiveCVEs(TwoPhaseModel<dim>& flow_model, const double deactivation_criteria);
        void                                  ZeroActiveCVEsProperty(size_t property_number);
        std::pair<double, double>         CFL_Timestep_MaxDivergenceActiveRegions();


        void                                  ConstructFluxesInteriorPatch(CVE_Patch<dim>& patch);
        void                                  ConstructFluxesBoundaryPatch(CVE_Patch<dim>& patch);
        void                                  ConstructFluxesActiveRegions();

    public:
        typename std::vector<ElementFace<dim>*>::iterator                                  ActiveInteriorLineFacesBegin();
        typename std::vector<ElementFace<dim>*>::iterator                                  ActiveInteriorLineFacesEnd();
        typename std::vector<ElementFace<dim>*>::iterator                                  ActiveInteriorPointFacesBegin();
        typename std::vector<ElementFace<dim>*>::iterator                                  ActiveInteriorPointFacesEnd();
        typename std::vector<ElementFace<dim>*>::iterator                                  ActiveModelBoundaryFacesBegin();
        typename std::vector<ElementFace<dim>*>::iterator                                  ActiveModelBoundaryFacesEnd();
        typename std::map<ElementFace<dim>*, double>::iterator                           MonitoringLineFacesFluxBegin();
        typename std::map<ElementFace<dim>*, double>::iterator                           MonitoringLineFacesFluxEnd();
        typename std::map<ElementFace<dim>*, double>::iterator                           MonitoringPointFacesFluxBegin();
        typename std::map<ElementFace<dim>*, double>::iterator                           MonitoringPointFacesFluxEnd();
        typename std::vector<ControlVolumeElement<dim>*>::iterator                         Active1D2D_CVEsBegin();
        typename std::vector<ControlVolumeElement<dim>*>::iterator                         Active1D2D_CVEsEnd();
        typename std::vector<ControlVolumeElement<dim>*>::iterator                         Active0D_CVEsBegin();
        typename std::vector<ControlVolumeElement<dim>*>::iterator                         Active0D_CVEsEnd();
        typename std::vector<CVE_Patch<dim>*>::iterator                                    ActiveInteriorPatchesBegin();
        typename std::vector<CVE_Patch<dim>*>::iterator                                    ActiveInteriorPatchesEnd();
        typename std::vector<CVE_Patch<dim>*>::iterator                                    ActiveBoundaryPatchesBegin();
        typename std::vector<CVE_Patch<dim>*>::iterator                                    ActiveBoundaryPatchesEnd();


    private:
        std::vector<ElementFace<dim>*>                                    activeInteriorLineFaces_;
        std::vector<ElementFace<dim>*>                                    activeInteriorPointFaces_;
        std::vector<ElementFace<dim>*>                                    activeModelBoundaryFaces_;
        std::vector<ControlVolumeElement<dim>*>                           active1D2D_CVEs_;
        std::vector<ControlVolumeElement<dim>*>                           active0D_CVEs_;
        std::map<ElementFace<dim>*, double>                             monitoringLineFacesFlux_; // stores face pointer, base flux and face number in the active cve
        std::map<ElementFace<dim>*, double>                             monitoringPointFacesFlux_; // stores face pointer, base flux and face number in the active cve
        std::vector<CVE_Patch<dim>*>                                      activeInteriorPatches_;
        std::vector<CVE_Patch<dim>*>                                      activeBoundaryPatches_;

        LU_Solver solver_;
        DenseMatrix<DM_MIN> A_;
        std::vector<double> b_;
        std::vector<double> fluxes_;
        std::vector<double> sector_face_area_;
        VectorVariable<dim> vt_;
        VectorVariable<dim> vt_inside_;
        VectorVariable<dim> vt_outside_;

    };


template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveInteriorLineFacesBegin()
{
    return activeInteriorLineFaces_.begin();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveInteriorLineFacesEnd()
{
    return activeInteriorLineFaces_.end();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveInteriorPointFacesBegin()
{
    return activeInteriorPointFaces_.begin();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveInteriorPointFacesEnd()
{
    return activeInteriorPointFaces_.end();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveModelBoundaryFacesBegin()
{
    return activeModelBoundaryFaces_.begin();
}



template<uint32_t dim>
inline typename std::vector<ElementFace<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveModelBoundaryFacesEnd()
{
    return activeModelBoundaryFaces_.end();
}



template<uint32_t dim>
inline typename std::map<ElementFace<dim>*, double>::iterator ActiveFaceTwoPhaseTransport<dim>::MonitoringLineFacesFluxBegin()
{
    return monitoringLineFacesFlux_.begin();
}



template<uint32_t dim>
inline typename std::map<ElementFace<dim>*, double>::iterator ActiveFaceTwoPhaseTransport<dim>::MonitoringLineFacesFluxEnd()
{
    return monitoringLineFacesFlux_.end();
}



template<uint32_t dim>
inline typename std::map<ElementFace<dim>*, double>::iterator ActiveFaceTwoPhaseTransport<dim>::MonitoringPointFacesFluxBegin()
{
    return monitoringPointFacesFlux_.begin();
}



template<uint32_t dim>
inline typename std::map<ElementFace<dim>*, double>::iterator ActiveFaceTwoPhaseTransport<dim>::MonitoringPointFacesFluxEnd()
{
    return monitoringPointFacesFlux_.end();
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim>*> ::iterator ActiveFaceTwoPhaseTransport<dim>::Active1D2D_CVEsBegin()
{
    return active1D2D_CVEs_.begin();
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim>*> ::iterator ActiveFaceTwoPhaseTransport<dim>::Active1D2D_CVEsEnd()
{
    return active1D2D_CVEs_.end();
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim>*> ::iterator ActiveFaceTwoPhaseTransport<dim>::Active0D_CVEsBegin()
{
    return active0D_CVEs_.begin();
}



template<uint32_t dim>
inline typename std::vector<ControlVolumeElement<dim>*> ::iterator ActiveFaceTwoPhaseTransport<dim>::Active0D_CVEsEnd()
{
    return active0D_CVEs_.end();
}


template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveInteriorPatchesBegin()
{
    return activeInteriorPatches_.begin();
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveInteriorPatchesEnd()
{
    return activeInteriorPatches_.end();
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveBoundaryPatchesBegin()
{
    return activeBoundaryPatches_.begin();
}



template<uint32_t dim>
inline typename std::vector<CVE_Patch<dim>*>::iterator ActiveFaceTwoPhaseTransport<dim>::ActiveBoundaryPatchesEnd()
{
    return activeBoundaryPatches_.end();
}

} // end csmp

#endif
