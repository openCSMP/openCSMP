// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ActiveElementTwoPhaseTransport.h"
#include "Model.h"
#include "TransportModel.h"
#include "TwoPhaseModel.h"
#include "SignalHandler.h"
#include "IMPES_Setup.h"
#include "DenseMatrix.h"
#include "LU_Solver.h"

#include "ModelTime.h"
#include "PropertyHandle.h"
#include "VTK_Interface.h"
#include "VariableOperations.h"

namespace csmp {

template<uint32_t dim>
ActiveElementTwoPhaseTransport<dim>::ActiveElementTwoPhaseTransport(Model<dim>& model,
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
                                                                    IMPES_Setup<dim>& setup)
 : TwoPhaseElementBasedTransport<dim>(model, transport_model, region_name, porosity, 
                                      permeability, saturation_oil, saturation_water,
                                      divergence, total_velocity,
                                      element_oil_injection_volume_rate,
                                      element_water_injection_volume_rate,
                                      element_total_production_volume_rate,
                                      setup)
{
    this->transportModelRef_.ZeroCVEsProperty(0U);
}




template<uint32_t dim>
ActiveElementTwoPhaseTransport<dim>::~ActiveElementTwoPhaseTransport()
{
}




template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::ZeroActiveAndMonitoringCVEsProperty(size_t property_number)
{
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active1D2D_CVEs_.begin(); 
        cveit != active1D2D_CVEs_.end(); ++cveit)
    {
        (*cveit)->PropertyValue(property_number, 0.);
    }

    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active0D_CVEs_.begin(); 
        cveit != active0D_CVEs_.end(); ++cveit)
    {
        (*cveit)->PropertyValue(property_number, 0.);
    }

    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = monitoring1D2D_CVEs_.begin(); 
        cveit != monitoring1D2D_CVEs_.end(); ++cveit)
    {
        cveit->first->PropertyValue(property_number, 0.);
    }

    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = monitoring0D_CVEs_.begin(); 
        cveit != monitoring0D_CVEs_.end(); ++cveit)
    {
        cveit->first->PropertyValue(property_number, 0.);
    }

} // end ZeroActiveAndMonitoringCVEsProperty




template<uint32_t dim>
std::pair<double, double> ActiveElementTwoPhaseTransport<dim>::CFL_Timestep_MaxDivergenceActiveRegions()
{
    double min_cfl_timestep(std::numeric_limits<double>::max());
    double max_divergence(0.);

    // looping over active 2D and 1D cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active1D2D_CVEs_.begin();
         cveit != active1D2D_CVEs_.end(); ++cveit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); ++i)
        {
            total_flux += std::abs((*cveit)->Face(i)->PropertyValue(0U));
            flux_imbalance += (*cveit)->Face(i)->PropertyValue(0U) * (*cveit)->FaceNormalDirection(i);
        }

        flux_imbalance += ((*cveit)->E()->Read(this->totalProdRateKey_) + 
                           (*cveit)->E()->Read(this->oilInjRateKey_) + 
                           (*cveit)->E()->Read(this->waterInjRateKey_)) * 
                           (*cveit)->Volume();

        (*cveit)->E()->Store(this->divergenceKey_, ScalarVariable(PLAIN, flux_imbalance));

        const double element_timestep((*cveit)->Volume() * (*cveit)->E()->Read(this->porosityKey_) * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);

    }

    // looping over active 0D cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active0D_CVEs_.begin();
         cveit != active0D_CVEs_.end(); ++cveit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); ++i)
        {
            total_flux += std::abs((*cveit)->Face(i)->PropertyValue(0U));
            flux_imbalance += (*cveit)->Face(i)->PropertyValue(0U) * (*cveit)->FaceNormalDirection(i);
        }

        const double element_timestep((*cveit)->Volume() * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);

    }

    // looping over monitoring 2D and 1D cves
    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = monitoring1D2D_CVEs_.begin();
         cveit != monitoring1D2D_CVEs_.end(); ++cveit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (uint32_t i = 0U; i < cveit->first->Faces(); ++i)
        {
            total_flux += std::abs(cveit->first->Face(i)->PropertyValue(0U));
            flux_imbalance += cveit->first->Face(i)->PropertyValue(0U) * cveit->first->FaceNormalDirection(i);
        }

        flux_imbalance += (cveit->first->E()->Read(this->totalProdRateKey_) + 
                           cveit->first->E()->Read(this->oilInjRateKey_) + 
                           cveit->first->E()->Read(this->waterInjRateKey_)) * 
                           cveit->first->Volume();

        cveit->first->E()->Store(this->divergenceKey_, ScalarVariable(PLAIN, flux_imbalance));

        const double element_timestep(cveit->first->Volume() * cveit->first->E()->Read(this->porosityKey_) * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);

    }

    // looping over monitoring 0D cves
    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = monitoring0D_CVEs_.begin();
         cveit != monitoring0D_CVEs_.end(); ++cveit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (uint32_t i = 0; i < cveit->first->Faces(); ++i)
        {
            total_flux += std::abs(cveit->first->Face(i)->PropertyValue(0U));
            flux_imbalance += cveit->first->Face(i)->PropertyValue(0U) * cveit->first->FaceNormalDirection(i);
        }

        const double element_timestep(cveit->first->Volume() * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);

    }

    return std::make_pair(min_cfl_timestep, max_divergence);

} // end CFL_Timestep_MaxDivergenceActiveCVEs


/*

// This is the modified method which uses the square root of permeability instead of permeability itself
// to increase the flow inside low permeability zones


template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::ConstructFluxesInteriorPatch(CVE_Patch<dim>& patch)
{
    static double perm_inside;
    static double perm_outside;

    // in this method (ConstructFluxesInteriorPatch) it is assumed that the influx to the sector/element is positive
    // and outflux from the sector/element is negative

    b_.resize(patch.Faces());
    fluxes_.resize(patch.Faces());
    sector_face_area_.resize(patch.Faces());
    A_.Resize(patch.Faces(), patch.Faces());
    A_.Zero();

    double min_perm(std::numeric_limits<double>::max());

    // we construct one less equation since it is a cyclic problem
    for (size_t s = 0U; s < patch.Sectors() - 1U; ++s)
    {            
        if (patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_) < min_perm)
            min_perm = patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_);

        if (patch.SectorAccessor(s)->CVE()->Dimension() != 0U)
        {
            patch.SectorAccessor(s)->E()->Read(this->totalVelocityKey_, vt_);

            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);

            }

            b_[s] = vt_.DotProduct(patch.SectorAccessor(s)->TotalFacetsAreaVector());

            b_[s] -= (patch.SectorAccessor(s)->E()->Read(this->oilInjRateKey_) + 
                    patch.SectorAccessor(s)->E()->Read(this->waterInjRateKey_) +
                    patch.SectorAccessor(s)->E()->Read(this->totalProdRateKey_)) *
                    patch.SectorAccessor(s)->CVE()->Volume() / 
                    patch.SectorAccessor(s)->E()->Nodes();

        }
        else
        {
            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);

                b_[s] = 0.;

            }

        }

    }

    if (patch.SectorAccessor(patch.Sectors() - 1U)->CVE()->Dimension() != 0)
    {
        if (patch.SectorAccessor(patch.Sectors() - 1U)->E()->Read(this->permeabilityKey_) < min_perm)
            min_perm = patch.SectorAccessor(patch.Sectors() - 1U)->E()->Read(this->permeabilityKey_);
    }

    // constructing auxiliary equations to close the system
    double total_flux(0.);
    size_t equation_number(patch.Sectors() - 1U);
    std::vector<bool> not_visited_face(patch.Faces(), true);

    for (size_t s = 0U; s < patch.Sectors(); ++s)
    {
        // each sub-patch (a sectors collection) which is isolated by two line elements gives an
        // auxiliary equation
        if (patch.SectorAccessor(s)->CVE()->Dimension() == 2U)
        {
            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                if (not_visited_face[patch.SectorAccessor(s)->PatchFaceNumber(f)])
                {
                    not_visited_face[patch.SectorAccessor(s)->PatchFaceNumber(f)] = false;
                    perm_inside = std::sqrt(patch.SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(this->permeabilityKey_)); // --------->
                    perm_outside = std::sqrt(patch.SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(this->permeabilityKey_)); // ---------->
                    patch.SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(this->totalVelocityKey_, vt_inside_);
                    patch.SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(this->totalVelocityKey_, vt_outside_);

                    total_flux -= vt_inside_.DotProduct(patch.SectorAccessor(s)->FaceNormal(f)) * 
                        patch.SectorAccessor(s)->FaceArea(f) / perm_inside
                        * patch.FaceCyclicDirection(patch.SectorAccessor(s)->PatchFaceNumber(f));
                    total_flux -= vt_outside_.DotProduct(patch.SectorAccessor(s)->FaceNormal(f)) * 
                        patch.SectorAccessor(s)->FaceArea(f) / perm_outside
                        * patch.FaceCyclicDirection(patch.SectorAccessor(s)->PatchFaceNumber(f));

                    A_(equation_number, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (1. / perm_inside + 1. / perm_outside) * min_perm * patch.SectorAccessor(s)->FaceArea(f)
                        * patch.FaceCyclicDirection(patch.SectorAccessor(s)->PatchFaceNumber(f));
                }

            }

            if (s == patch.Sectors() - 1U)
            {
                b_[equation_number] = total_flux * min_perm;
            }

        }
        else
        {
            if (equation_number < patch.Faces()) b_[equation_number] = total_flux * min_perm;
            total_flux = 0.;
            ++equation_number;
            continue;
        }

    }

    solver_.Solve(A_, b_, fluxes_, patch.Faces());

    // assigning the results to the corresponding faces
    typename std::vector<double>::iterator rit(fluxes_.begin());
    typename std::vector<double>::iterator ait(sector_face_area_.begin());
    for (typename std::vector<ElementFace<dim>*>::iterator fit = patch.FacesBegin(); fit != patch.FacesEnd(); ++fit, ++rit, ++ait)
    {
        (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
    }

} // end ConstructFluxesInteriorPatch





template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::ConstructFluxesBoundaryPatch(CVE_Patch<dim>& patch)
{
    static double perm_inside;
    static double perm_outside;

    // construct fluxes for boundary patches
    b_.resize(patch.Faces());
    fluxes_.resize(patch.Faces());
    sector_face_area_.resize(patch.Faces());
    A_.Resize(patch.Faces(), patch.Faces());
    A_.Zero();

    double min_perm(std::numeric_limits<double>::max());

    for (size_t s = 0U; s < patch.Sectors(); ++s)
    {            
        if (patch.SectorAccessor(s)->CVE()->Dimension() != 0U)
        {
            if (patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_) < min_perm)
                min_perm = patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_);

            patch.SectorAccessor(s)->E()->Read(this->totalVelocityKey_, vt_);

            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) =  
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);
            }

            b_[s] = vt_.DotProduct(patch.SectorAccessor(s)->TotalFacetsAreaVector());

            b_[s] -= (patch.SectorAccessor(s)->E()->Read(this->oilInjRateKey_) + 
                    patch.SectorAccessor(s)->E()->Read(this->waterInjRateKey_) +
                    patch.SectorAccessor(s)->E()->Read(this->totalProdRateKey_)) *
                    patch.SectorAccessor(s)->CVE()->Volume() / 
                    patch.SectorAccessor(s)->E()->Nodes();
        }
        else
        {
            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);

                b_[s] = 0.;

            }

        }

    }

    // constructing auxiliary equations to close the system
    size_t equation_number(patch.Sectors());
    double perm_boundary_face_1(std::sqrt(patch.FaceAccessor(patch.InteriorFaces())->InsideCVE()->E()->Read(this->permeabilityKey_))); // ----------->
    double area_boundary_face_1(patch.FaceAccessor(patch.InteriorFaces())->Area() / 2.);

    for (size_t f = patch.InteriorFaces() + 1U; f < patch.Faces(); ++f)
    {
        if (patch.FaceAccessor(f)->InsideCVE()->Dimension() != 0U)
        {
            double perm_boundary_face_2(std::sqrt(patch.FaceAccessor(f)->InsideCVE()->E()->Read(this->permeabilityKey_))); // -------------------->
            double area_boundary_face_2(patch.FaceAccessor(f)->InsideCVE()->Dimension() == 2U ?
                                          patch.FaceAccessor(f)->Area() / 2. : patch.FaceAccessor(f)->Area());

            b_[equation_number] = 0.;
            A_(equation_number, patch.InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
            A_(equation_number, f) = -area_boundary_face_2 / perm_boundary_face_2;

            ++equation_number;
        }
        else
        {
            for (size_t n = 0U; n < patch.FaceAccessor(f)->InsideCVE()->Faces(); ++n)
            {
                if (patch.FaceAccessor(f)->InsideCVE()->Face(n)->Placement() == NOT)
                {
                    double perm_boundary_face_2(std::sqrt(patch.FaceAccessor(f)->InsideCVE()->Neighbor(n)->E()->Read(this->permeabilityKey_))); // ----------->
                    double area_boundary_face_2(patch.FaceAccessor(f)->InsideCVE()->Face(n)->Area());

                    // finding the face number in the patch
                    size_t face_num;
                    for (face_num = 0U; face_num < patch.InteriorFaces(); ++face_num)
                    {
                        if (patch.FaceAccessor(face_num) == patch.FaceAccessor(f)->InsideCVE()->Face(n)) break;
                    }

                    b_[equation_number] = 0.;
                    A_(equation_number, patch.InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
                    A_(equation_number, face_num) = area_boundary_face_2 / perm_boundary_face_2 * 
                                                   patch.FaceAccessor(f)->InsideCVE()->FaceNormalDirection(n);

                    ++equation_number;
                }

            }

        }

    }

    solver_.Solve(A_, b_, fluxes_, patch.Faces());

    // assigning the results to the corresponding faces
    typename std::vector<double>::iterator rit(fluxes_.begin());
    typename std::vector<double>::iterator ait(sector_face_area_.begin());
    for (typename std::vector<ElementFace<dim>*>::iterator fit = patch.FacesBegin(); fit != patch.FacesEnd(); ++fit, ++rit, ++ait)
    {
        (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
    }

} // end ConstructFluxesBoundaryPatch


*/




template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::ConstructFluxesInteriorPatch(CVE_Patch<dim>& patch)
{
    static double perm_inside;
    static double perm_outside;

    // in this method (ConstructFluxesInteriorPatch) it is assumed that the influx to the sector/element is positive
    // and outflux from the sector/element is negative

    b_.resize(patch.Faces());
    fluxes_.resize(patch.Faces());
    sector_face_area_.resize(patch.Faces());
    A_.Resize(patch.Faces(), patch.Faces());
    A_.Zero();

    double min_perm(std::numeric_limits<double>::max());

    // we construct one less equation since it is a cyclic problem
    for (uint32_t s = 0U; s < patch.Sectors() - 1U; ++s)
    {            
        if (patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_) < min_perm)
            min_perm = patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_);

        if (patch.SectorAccessor(s)->CVE()->Dimension() != 0U)
        {
            patch.SectorAccessor(s)->E()->Read(this->totalVelocityKey_, vt_);

            for (uint32_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);

            }

            b_[s] = dotProduct(vt_, patch.SectorAccessor(s)->TotalFacetsAreaVector());

            b_[s] -= (patch.SectorAccessor(s)->E()->Read(this->oilInjRateKey_) + 
                    patch.SectorAccessor(s)->E()->Read(this->waterInjRateKey_) +
                    patch.SectorAccessor(s)->E()->Read(this->totalProdRateKey_)) *
                    patch.SectorAccessor(s)->CVE()->Volume() / 
                    patch.SectorAccessor(s)->E()->Nodes();

        }
        else
        {
            for (uint32_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);

                b_[s] = 0.;

            }

        }

    }

    if (patch.SectorAccessor(patch.Sectors() - 1U)->CVE()->Dimension() != 0)
    {
        if (patch.SectorAccessor(patch.Sectors() - 1U)->E()->Read(this->permeabilityKey_) < min_perm)
            min_perm = patch.SectorAccessor(patch.Sectors() - 1U)->E()->Read(this->permeabilityKey_);
    }

    // constructing auxiliary equations to close the system
    double total_flux(0.);
    uint32_t equation_number(patch.Sectors() - 1U);
    std::vector<bool> not_visited_face(patch.Faces(), true);

    for (uint32_t s = 0U; s < patch.Sectors(); ++s)
    {
        // each sub-patch (a sectors collection) which is isolated by two line elements gives an
        // auxiliary equation
        if (patch.SectorAccessor(s)->CVE()->Dimension() == 2U)
        {
            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                if (not_visited_face[patch.SectorAccessor(s)->PatchFaceNumber(f)])
                {
                    not_visited_face[patch.SectorAccessor(s)->PatchFaceNumber(f)] = false;
                    perm_inside = patch.SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(this->permeabilityKey_);
                    perm_outside = patch.SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(this->permeabilityKey_);
                    patch.SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(this->totalVelocityKey_, vt_inside_);
                    patch.SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(this->totalVelocityKey_, vt_outside_);

                    total_flux -= dotProduct(vt_inside_,patch.SectorAccessor(s)->FaceNormal(f)) *
                        patch.SectorAccessor(s)->FaceArea(f) / perm_inside
                        * patch.FaceCyclicDirection(patch.SectorAccessor(s)->PatchFaceNumber(f));
                    total_flux -= dotProduct(vt_outside_,patch.SectorAccessor(s)->FaceNormal(f)) *
                        patch.SectorAccessor(s)->FaceArea(f) / perm_outside
                        * patch.FaceCyclicDirection(patch.SectorAccessor(s)->PatchFaceNumber(f));

                    A_(equation_number, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (1. / perm_inside + 1. / perm_outside) * min_perm * patch.SectorAccessor(s)->FaceArea(f)
                        * patch.FaceCyclicDirection(patch.SectorAccessor(s)->PatchFaceNumber(f));
                }

            }

            if (s == patch.Sectors() - 1U)
            {
                b_[equation_number] = total_flux * min_perm;
            }

        }
        else
        {
            if (equation_number < patch.Faces()) b_[equation_number] = total_flux * min_perm;
            total_flux = 0.;
            ++equation_number;
            continue;
        }

    }

    solver_.Solve(A_, b_, fluxes_, patch.Faces());

    // assigning the results to the corresponding faces
    typename std::vector<double>::iterator rit(fluxes_.begin());
    typename std::vector<double>::iterator ait(sector_face_area_.begin());
    for (typename std::vector<ElementFace<dim>*>::iterator fit = patch.FacesBegin(); fit != patch.FacesEnd(); ++fit, ++rit, ++ait)
    {
        (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
    }

} // end ConstructFluxesInteriorPatch





template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::ConstructFluxesBoundaryPatch(CVE_Patch<dim>& patch)
{
    // construct fluxes for boundary patches
    b_.resize(patch.Faces());
    fluxes_.resize(patch.Faces());
    sector_face_area_.resize(patch.Faces());
    A_.Resize(patch.Faces(), patch.Faces());
    A_.Zero();

    double min_perm(std::numeric_limits<double>::max());

    for (uint32_t s = 0U; s < patch.Sectors(); ++s)
    {            
        if (patch.SectorAccessor(s)->CVE()->Dimension() != 0U)
        {
            if (patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_) < min_perm)
                min_perm = patch.SectorAccessor(s)->E()->Read(this->permeabilityKey_);

            patch.SectorAccessor(s)->E()->Read(this->totalVelocityKey_, vt_);

            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) =  
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);
            }

            b_[s] = dotProduct(vt_,patch.SectorAccessor(s)->TotalFacetsAreaVector());

            b_[s] -= (patch.SectorAccessor(s)->E()->Read(this->oilInjRateKey_) + 
                    patch.SectorAccessor(s)->E()->Read(this->waterInjRateKey_) +
                    patch.SectorAccessor(s)->E()->Read(this->totalProdRateKey_)) *
                    patch.SectorAccessor(s)->CVE()->Volume() / 
                    patch.SectorAccessor(s)->E()->Nodes();
        }
        else
        {
            for (size_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
            {
                // storing the sector face area for later uses
                sector_face_area_[patch.SectorAccessor(s)->PatchFaceNumber(f)] = patch.SectorAccessor(s)->FaceArea(f);

                A_(s, patch.SectorAccessor(s)->PatchFaceNumber(f)) = 
                    patch.SectorAccessor(s)->FaceNormalDirection(f) * patch.SectorAccessor(s)->FaceArea(f);

                b_[s] = 0.;

            }
        }
    }

    // constructing auxiliary equations to close the system
    uint32_t equation_number(patch.Sectors());
    double perm_boundary_face_1(patch.FaceAccessor(patch.InteriorFaces())->InsideCVE()->E()->Read(this->permeabilityKey_));
    double area_boundary_face_1(patch.FaceAccessor(patch.InteriorFaces())->Area() / 2.);

    for (uint32_t f = patch.InteriorFaces() + 1U; f < patch.Faces(); ++f)
    {
        if (patch.FaceAccessor(f)->InsideCVE()->Dimension() != 0U)
        {
            double perm_boundary_face_2(patch.FaceAccessor(f)->InsideCVE()->E()->Read(this->permeabilityKey_));
            double area_boundary_face_2(patch.FaceAccessor(f)->InsideCVE()->Dimension() == 2U ?
                                          patch.FaceAccessor(f)->Area() / 2. : patch.FaceAccessor(f)->Area());

            b_[equation_number] = 0.;
            A_(equation_number, patch.InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
            A_(equation_number, f) = -area_boundary_face_2 / perm_boundary_face_2;

            ++equation_number;
        }
        else
        {
            for (uint32_t n = 0U; n < patch.FaceAccessor(f)->InsideCVE()->Faces(); ++n)
            {
                if (patch.FaceAccessor(f)->InsideCVE()->Face(n)->Placement() == NOT)
                {
                    double perm_boundary_face_2(patch.FaceAccessor(f)->InsideCVE()->Neighbor(n)->E()->Read(this->permeabilityKey_));
                    double area_boundary_face_2(patch.FaceAccessor(f)->InsideCVE()->Face(n)->Area());

                    // finding the face number in the patch
                    uint32_t face_num{0U};
                    for ( ; face_num < patch.InteriorFaces(); ++face_num)
                    {
                        if (patch.FaceAccessor(face_num) == patch.FaceAccessor(f)->InsideCVE()->Face(n)) break;
                    }

                    b_[equation_number] = 0.;
                    A_(equation_number, patch.InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
                    A_(equation_number, face_num) = area_boundary_face_2 / perm_boundary_face_2 * 
                                                   patch.FaceAccessor(f)->InsideCVE()->FaceNormalDirection(n);

                    ++equation_number;
                }
            }
        }
    }

    solver_.Solve(A_, b_, fluxes_, patch.Faces());

    // assigning the results to the corresponding faces
    typename std::vector<double>::iterator rit(fluxes_.begin());
    typename std::vector<double>::iterator ait(sector_face_area_.begin());
    for (typename std::vector<ElementFace<dim>*>::iterator fit = patch.FacesBegin(); fit != patch.FacesEnd(); ++fit, ++rit, ++ait)
    {
        (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
    }

} // end ConstructFluxesBoundaryPatch




template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::ConstructFluxesActiveRegions()
{
    std::cout << "\nActiveElementTwoPhaseTransport<" << dim << ">::ConstructFluxesActiveRegions()\n";

    this->transportModelRef_.ZeroFacesProperty(0U);

    for (typename std::vector<CVE_Patch<dim>*>::iterator pit = ActiveInteriorPatchesBegin();
         pit != ActiveInteriorPatchesEnd(); ++pit)
    {
        ConstructFluxesInteriorPatch(*(*pit));
    }

    for (typename std::vector<CVE_Patch<dim>*>::iterator pit = ActiveBoundaryPatchesBegin();
        pit != ActiveBoundaryPatchesEnd(); ++pit)
    {
        ConstructFluxesBoundaryPatch(*(*pit));
    }

} // end ConstructFluxesActiveRegions





template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::FindActiveRegions(TwoPhaseModel<dim>& flow_model, const double activation_criteria)
{
    std::cout << "\nActiveElementTwoPhaseTransport<" << dim << ">::FindActiveRegions()\n";

    activeInteriorLineFaces_.resize(0U);
    activeInteriorPointFaces_.resize(0U);
    activeModelBoundaryFaces_.resize(0U);
    active1D2D_CVEs_.resize(0U);
    active0D_CVEs_.resize(0U);
    monitoring1D2D_CVEs_.clear();
    monitoring0D_CVEs_.clear();
    monitoringLineFaces_.clear();
    monitoringPointFaces_.clear();
    activeInteriorPatches_.resize(0U);
    activeBoundaryPatches_.resize(0U);

    // changing the face active flag to false
    for (typename std::vector<ElementFace<dim> >::iterator fit = this->transportModelRef_.FacesBegin(); 
         fit != this->transportModelRef_.FacesEnd(); ++fit)
    {
        fit->Active(false);
    }

    // changing the patch active flag to false
    for (typename std::vector<CVE_Patch<dim> >::iterator pit = this->transportModelRef_.PatchesBegin(); 
        pit != this->transportModelRef_.PatchesEnd(); ++pit)
    {
        pit->Active(false);
    }

    // looping over inside 2D and 1D cves to find active cves
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = this->transportModelRef_.CVEsBegin();
         cveit != this->transportModelRef_.PointCVEsBegin(); ++cveit)
    {
        cveit->Active(false);

        double accumulation(0.);
        double divergence(0.);

        for (size_t n = 0U; n < cveit->Faces(); ++n)
        {
            if ((cveit->Face(n)->InsideCVE()->Dimension() == 0U) || 
                (cveit->Face(n)->OutsideCVE()->Dimension() == 0U))
            {
                accumulation += this->CalculateOilFluxInteriorPointFace(*(cveit->Face(n)), flow_model) * cveit->FaceNormalDirection(n);
            } 
            else
            {
                accumulation += this->CalculateOilFluxInteriorLineFace(*(cveit->Face(n)), flow_model) * cveit->FaceNormalDirection(n);
            }

            divergence += cveit->Face(n)->PropertyValue(0U) * cveit->FaceNormalDirection(n);

        } // end loop faces

        flow_model.Initialize(*(cveit->E()));
        flow_model.EffectiveSaturation();
        const double lambda_w_in(flow_model.MobilityPhase(1U));
        const double lambda_o_in(flow_model.MobilityPhase(2U));

        const double oil_source_sink(cveit->Volume() * (cveit->E()->Read(this->oilInjRateKey_) +
                                       cveit->E()->Read(this->totalProdRateKey_) * lambda_o_in / (lambda_w_in + lambda_o_in)));
        accumulation += oil_source_sink;
        accumulation -= divergence * lambda_o_in / (lambda_w_in + lambda_o_in);

        if (std::abs(accumulation) > activation_criteria)
        {
            active1D2D_CVEs_.push_back(&(*cveit));
            cveit->Active(true);

            // making the parent patches active and insert in the corresponding container
            for (typename std::vector<CVE_Patch<dim>*>::iterator pit = cveit->PatchesBegin(); 
                pit != cveit->PatchesEnd(); ++pit)
            {
                if (!((*pit)->Active()))
                {
                    (*pit)->Active(true);

                    if ((*pit)->CenterNode()->AtBoundary() == NOT)
                    {
                        activeInteriorPatches_.push_back(*pit);
                    } 
                    else
                    {
                        activeBoundaryPatches_.push_back(*pit);
                    }

                } // end if patch is not active

            } // end looping over the patches
 
        } // end if it should be an active cell

    } // end looping over inside 1D and 2D cves

    // looping over inside 0D cves to find active cves
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = this->transportModelRef_.PointCVEsBegin();
         cveit != this->transportModelRef_.PerimeterCVEsBegin(); ++cveit)
    {
        cveit->Active(false);

        double vis_w_avg(0.);
        double vis_o_avg(0.);

        for (size_t n = 0U; n < cveit->Neighbors(); ++n)
        {
            flow_model.Initialize(*(cveit->Neighbor(n)->E()));
            flow_model.EffectiveSaturation();

            vis_w_avg += flow_model.ViscosityWettingPhase();
            vis_o_avg += flow_model.ViscosityNonWettingPhase();

        }

        vis_w_avg /= (*cveit).Neighbors();
        vis_o_avg /= (*cveit).Neighbors();

        double lambda_w_in;
        double lambda_o_in;

        if (((*cveit).PropertyValue(1U) > 0.) && ((*cveit).PropertyValue(1U) < 1.))
        {
            lambda_o_in = (*cveit).PropertyValue(1U) / vis_o_avg;
            lambda_w_in = (1. - (*cveit).PropertyValue(1U)) / vis_w_avg;                    
        }
        else if ((*cveit).PropertyValue(1U) >= 1.)
        {
            lambda_o_in = 1. / vis_o_avg;
            lambda_w_in = 0.;
        }
        else
        {
            lambda_o_in = 0.;
            lambda_w_in = 1. / vis_w_avg;
        }
        
        double accumulation(0.);
        double divergence(0.);

        for (size_t n = 0U; n < cveit->Faces(); ++n)
        {
            accumulation += this->CalculateOilFluxInteriorPointFace(*(cveit->Face(n)), flow_model) * cveit->FaceNormalDirection(n);
            divergence += cveit->Face(n)->PropertyValue(0U) * cveit->FaceNormalDirection(n);

        } // end looping over neighbors

        accumulation -= divergence * lambda_o_in / (lambda_w_in + lambda_o_in);

        if (std::abs(accumulation) > activation_criteria)
        {
            active0D_CVEs_.push_back(&(*cveit));
            cveit->Active(true);

            // making the parent patches active and insert in the corresponding container
            for (typename std::vector<CVE_Patch<dim>*>::iterator pit = cveit->PatchesBegin(); 
                pit != cveit->PatchesEnd(); ++pit)
            {
                if (!((*pit)->Active()))
                {
                    (*pit)->Active(true);
                    activeInteriorPatches_.push_back(*pit);
                } // end if patch is not active

            } // end looping over the patches

        }

    } // end looping over inside 0D cves


    // looping over boundary 2D and 1D cves to find active cves
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = this->transportModelRef_.PerimeterCVEsBegin();
         cveit != this->transportModelRef_.PerimeterPointCVEsBegin(); ++cveit)
    {
        cveit->Active(false);

        double accumulation(0.);
        double divergence(0.);

        for (size_t n = 0U; n < cveit->Faces(); ++n)
        {
            if (cveit->Face(n)->Placement() != NOT)
            {
                accumulation += this->CalculateOilFluxBoundaryLineFace(*(cveit->Face(n)), flow_model) * cveit->FaceNormalDirection(n);
            } 
            else
            {
                if ((cveit->Face(n)->InsideCVE()->Dimension() == 0U) || 
                    (cveit->Face(n)->OutsideCVE()->Dimension() == 0U))
                {
                    accumulation += this->CalculateOilFluxInteriorPointFace(*(cveit->Face(n)), flow_model) * cveit->FaceNormalDirection(n);
                } 
                else
                {
                    accumulation += this->CalculateOilFluxInteriorLineFace(*(cveit->Face(n)), flow_model) * cveit->FaceNormalDirection(n);
                }

            }

            divergence += cveit->Face(n)->PropertyValue(0U) * cveit->FaceNormalDirection(n);

        } // end loop faces

        flow_model.Initialize(*(cveit->E()));
        flow_model.EffectiveSaturation();
        const double lambda_w_in(flow_model.MobilityPhase(1U));
        const double lambda_o_in(flow_model.MobilityPhase(2U));

        const double oil_source_sink(cveit->Volume() * (cveit->E()->Read(this->oilInjRateKey_) +
                                       cveit->E()->Read(this->totalProdRateKey_) * lambda_o_in / (lambda_w_in + lambda_o_in)));
        accumulation += oil_source_sink;
        accumulation -= divergence * lambda_o_in / (lambda_w_in + lambda_o_in);

        if (std::abs(accumulation) > activation_criteria)
        {
            active1D2D_CVEs_.push_back(&(*cveit));
            cveit->Active(true);

            // making the parent patches active and insert in the corresponding container
            for (typename std::vector<CVE_Patch<dim>*>::iterator pit = cveit->PatchesBegin(); 
                pit != cveit->PatchesEnd(); ++pit)
            {
                if (!((*pit)->Active()))
                {
                    (*pit)->Active(true);

                    if ((*pit)->CenterNode()->AtBoundary() == NOT)
                    {
                        activeInteriorPatches_.push_back(*pit);
                    } 
                    else
                    {
                        activeBoundaryPatches_.push_back(*pit);
                    }

                } // end if patch is not active

            } // end looping over the patches

        }

    } // end looping over inside 1D and 2D cves
    
    std::multiset<ElementFace<dim>*> monitoring_line_faces;
    std::multiset<ElementFace<dim>*> monitoring_point_faces;
    std::set<ControlVolumeElement<dim>*> monitoring_1d2d_cves;
    std::set<ControlVolumeElement<dim>*> monitoring_0d_cves;

    // finding active faces and monitoring cves (looping over 1D-2D cves)
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active1D2D_CVEs_.begin();
         cveit != active1D2D_CVEs_.end(); ++cveit)
    {
        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            if ((*cveit)->Face(n)->Placement() == NOT)
            {
                // make the face active and inserting the face in the corresponding active face container
                if (!((*cveit)->Face(n)->Active()))
                {
                    (*cveit)->Face(n)->Active(true);

                    if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                    {
                        activeInteriorPointFaces_.push_back((*cveit)->Face(n));
                    } 
                    else
                    {
                        activeInteriorLineFaces_.push_back((*cveit)->Face(n));
                    }

                }
                // if the neighbor cve is not active put it in the monitoring container
                if (!((*cveit)->Neighbor(n)->Active()))
                {
                    if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                    {
                        monitoring_0d_cves.insert((*cveit)->Neighbor(n));
                    } 
                    else
                    {
                        monitoring_1d2d_cves.insert((*cveit)->Neighbor(n));
                    }
                }
            }
            // if it is a boundary face
            else
            {
                (*cveit)->Face(n)->Active(true);
                activeModelBoundaryFaces_.push_back((*cveit)->Face(n));
            }
        
        } // end looping over cve's faces

    } // end looping over 1d2d cves

    // finding active faces and monitoring cves (looping over 0D cves)
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active0D_CVEs_.begin();
        cveit != active0D_CVEs_.end(); ++cveit)
    {
        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // make the face active and inserting the face in the corresponding active face container
            if (!((*cveit)->Face(n)->Active()))
            {
                (*cveit)->Face(n)->Active(true);
                activeInteriorPointFaces_.push_back((*cveit)->Face(n));
            }
            // if the neighbor cve is not active put it in the monitoring container
            if (!((*cveit)->Neighbor(n)->Active()))
            {
                    monitoring_1d2d_cves.insert((*cveit)->Neighbor(n));
            }

        } // end looping over cve's faces

    } // end looping over 0d cves

    // looping over monitoring 1d2d cves to find monitoring faces and make the parent cve's patches active
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = monitoring_1d2d_cves.begin();
        cveit != monitoring_1d2d_cves.end(); ++cveit)
    {
        // we set the flag for monitoring cves as active
        (*cveit)->Active(true);

        // making the parent patches active and insert in the corresponding container
        for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
            pit != (*cveit)->PatchesEnd(); ++pit)
        {
            if (!((*pit)->Active()))
            {
                (*pit)->Active(true);

                if ((*pit)->CenterNode()->AtBoundary() == NOT)
                {
                    activeInteriorPatches_.push_back(*pit);
                } 
                else
                {
                    activeBoundaryPatches_.push_back(*pit);
                }

            } // end if patch is not active

        } // end looping over the patches

        // insert cve in the monitoring container with initial oil saturation as the value of the key
        flow_model.Initialize(*((*cveit)->E()));
        monitoring1D2D_CVEs_[*cveit] = flow_model.Saturation(2U);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // find inactive faces and inserting them in the corresponding monitoring face container
            if (!((*cveit)->Face(n)->Active()))
            {
                if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                {
                    monitoring_point_faces.insert((*cveit)->Face(n));
                } 
                else
                {
                    monitoring_line_faces.insert((*cveit)->Face(n));
                }

            }

        } // end looping over the faces

    } // end looping over 1d2d cves

    // looping over monitoring 0d cves to find monitoring faces and making parent cve's patches active
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = monitoring_0d_cves.begin();
        cveit != monitoring_0d_cves.end(); ++cveit)
    {
        // we set the flag for monitoring cves as active
        (*cveit)->Active(true);

        // making the parent patches active and insert in the corresponding container
        for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
            pit != (*cveit)->PatchesEnd(); ++pit)
        {
            if (!((*pit)->Active()))
            {
                (*pit)->Active(true);
                activeInteriorPatches_.push_back(*pit);
            } // end if patch is not active

        } // end looping over the patches

        // insert cve in the monitoring container with initial oil saturation as the value of the key
        monitoring0D_CVEs_[*cveit] = (*cveit)->PropertyValue(1U);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // find inactive faces and inserting them in the monitoring point face container
            if (!((*cveit)->Face(n)->Active()))
            {
                monitoring_point_faces.insert((*cveit)->Face(n));
            }

        } // end looping over the faces

    } // end looping over 0d cves

    // putting monitoring faces in the corresponding container
    // looping over monitoring line faces
    for (typename std::multiset<ElementFace<dim>*>::iterator fit = monitoring_line_faces.begin(); 
         fit != monitoring_line_faces.end(); ++fit)
    {
        // if the face is not shared between two monitoring cves
        if (monitoring_line_faces.count(*fit) == 1U)
        {
            monitoringLineFaces_.insert(*fit);
        }
        // if the face is shared between two monitoring cves we make it active
        else
        {
            if (!((*fit)->Active()))
            {
                (*fit)->Active(true);
                activeInteriorLineFaces_.push_back(*fit);
            }

        }

    } // end looping over monitoring line faces

    // looping over monitoring point faces
    for (typename std::multiset<ElementFace<dim>*>::iterator fit = monitoring_point_faces.begin(); 
        fit != monitoring_point_faces.end(); ++fit)
    {
        // if the face is not shared between two monitoring cves
        if (monitoring_point_faces.count(*fit) == 1U)
        {
            monitoringPointFaces_.insert(*fit);
        }
        // if the face is shared between two monitoring cves we make it active
        else
        {
            if (!((*fit)->Active()))
            {
                (*fit)->Active(true);
                activeInteriorPointFaces_.push_back(*fit);
            }

        }

    } // end looping over monitoring point faces

} // end FindActiveRegions




template<uint32_t dim>
double ActiveElementTwoPhaseTransport<dim>::UpdateActiveRegions(TwoPhaseModel<dim>& flow_model, double criteria)
{
    std::vector<ControlVolumeElement<dim>*> new_active_1d2d_cves;
    std::vector<ControlVolumeElement<dim>*> new_active_0d_cves;
    std::set<ControlVolumeElement<dim>*> new_monitoring_1d2d_cves;
    std::set<ControlVolumeElement<dim>*> new_monitoring_0d_cves;
    std::vector<CVE_Patch<dim>*> new_active_interior_patches;
    std::vector<CVE_Patch<dim>*> new_active_boundary_patches;

    // finding new active 1d2d cves
    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = monitoring1D2D_CVEs_.begin();
         cveit != monitoring1D2D_CVEs_.end(); ++cveit)
    {
        flow_model.Initialize(*(cveit->first->E()));

        if (std::abs(flow_model.Saturation(2U) - cveit->second) > criteria)
        {
            new_active_1d2d_cves.push_back(cveit->first);

            for (size_t i = 0U; i < cveit->first->Faces(); ++i)
            {
                if (!(cveit->first->Face(i)->Active()))
                {
                    cveit->first->Face(i)->Active(true);

                    if (cveit->first->Face(i)->Placement() == NOT)
                    {
                        if (cveit->first->Neighbor(i)->Dimension() == 0U)
                        {
                            activeInteriorPointFaces_.push_back(cveit->first->Face(i));
                            monitoringPointFaces_.erase(cveit->first->Face(i));
                            new_monitoring_0d_cves.insert(cveit->first->Neighbor(i));
                        } 
                        else
                        {
                            activeInteriorLineFaces_.push_back(cveit->first->Face(i));
                            monitoringLineFaces_.erase(cveit->first->Face(i));
                            new_monitoring_1d2d_cves.insert(cveit->first->Neighbor(i));
                        }

                    }
                    else
                    {
                        activeModelBoundaryFaces_.push_back(cveit->first->Face(i));
                        monitoringLineFaces_.erase(cveit->first->Face(i));
                    }

                }

            }

        } // end if it becomes an active cve

    } // end looping over monitoring 1d2d cves

    // finding new active 0d cves
    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = monitoring0D_CVEs_.begin();
        cveit != monitoring0D_CVEs_.end(); ++cveit)
    {
        if (std::abs(cveit->first->PropertyValue(1U) - cveit->second) > criteria)
        {
            new_active_0d_cves.push_back(cveit->first);

            for (size_t i = 0U; i < cveit->first->Faces(); ++i)
            {
                if (!(cveit->first->Face(i)->Active()))
                {
                    cveit->first->Face(i)->Active(true);
                    activeInteriorPointFaces_.push_back(cveit->first->Face(i));
                    monitoringPointFaces_.erase(cveit->first->Face(i));
                    new_monitoring_1d2d_cves.insert(cveit->first->Neighbor(i));
                }

            }

        } // end if it becomes an active cve

    } // end looping over monitoring 0d cves

    // removing new active cves from the monitoring container and adding them to the active container
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = new_active_1d2d_cves.begin();
        cveit != new_active_1d2d_cves.end(); ++cveit)
    {
        active1D2D_CVEs_.push_back(*cveit);
        monitoring1D2D_CVEs_.erase(*cveit);
    }

    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = new_active_0d_cves.begin();
        cveit != new_active_0d_cves.end(); ++cveit)
    {
        active0D_CVEs_.push_back(*cveit);
        monitoring0D_CVEs_.erase(*cveit);
    }

    // we set the flag for monitoring cves as active and activate their parent patches
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = new_monitoring_1d2d_cves.begin();
        cveit != new_monitoring_1d2d_cves.end(); ++cveit)
    {
        (*cveit)->Active(true);

        // making the parent patches active and insert in the corresponding container
        for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
            pit != (*cveit)->PatchesEnd(); ++pit)
        {
            if (!((*pit)->Active()))
            {
                (*pit)->Active(true);

                if ((*pit)->CenterNode()->AtBoundary() == NOT)
                {
                    activeInteriorPatches_.push_back(*pit);
                    new_active_interior_patches.push_back(*pit);
                } 
                else
                {
                    activeBoundaryPatches_.push_back(*pit);
                    new_active_boundary_patches.push_back(*pit);
                }

            } // end if patch is not active

        } // end looping over the patches

    }

    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = new_monitoring_0d_cves.begin();
        cveit != new_monitoring_0d_cves.end(); ++cveit)
    {
        (*cveit)->Active(true);

        for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
            pit != (*cveit)->PatchesEnd(); ++pit)
        {
            if (!((*pit)->Active()))
            {
                (*pit)->Active(true);
                activeInteriorPatches_.push_back(*pit);
                new_active_interior_patches.push_back(*pit);
            } // end if patch is not active

        } // end looping over the patches

    }

    // constructing fluxes for new active patches
    for (typename std::vector<CVE_Patch<dim>*>::iterator pit = new_active_interior_patches.begin();
         pit != new_active_interior_patches.end(); ++pit)
    {
        ConstructFluxesInteriorPatch(*(*pit));
    }

    for (typename std::vector<CVE_Patch<dim>*>::iterator pit = new_active_boundary_patches.begin();
        pit != new_active_boundary_patches.end(); ++pit)
    {
        ConstructFluxesBoundaryPatch(*(*pit));
    }


    // looping over new monitoring 1d2d cves to put them in the container and to find new monitoring faces
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = new_monitoring_1d2d_cves.begin();
        cveit != new_monitoring_1d2d_cves.end(); ++cveit)
    {
        // insert cve in the monitoring container with initial oil saturation as the value of the key
        flow_model.Initialize(*((*cveit)->E()));
        monitoring1D2D_CVEs_[*cveit] = flow_model.Saturation(2U);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // find inactive faces and inserting them in the corresponding monitoring face container
            if (!((*cveit)->Neighbor(n)->Active()))
            {
                if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                {
                    monitoringPointFaces_.insert((*cveit)->Face(n));
                } 
                else
                {
                    monitoringLineFaces_.insert((*cveit)->Face(n));
                }

            }
            else
            {
                if (!((*cveit)->Face(n)->Active()))
                {
                    (*cveit)->Face(n)->Active(true);

                    if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                    {
                        activeInteriorPointFaces_.push_back((*cveit)->Face(n));
                        monitoringPointFaces_.erase((*cveit)->Face(n));
                    } 
                    else
                    {
                        if ((*cveit)->Face(n)->Placement() == NOT)
                        {
                            activeInteriorLineFaces_.push_back((*cveit)->Face(n));
                            monitoringLineFaces_.erase((*cveit)->Face(n));
                        } 
                        else
                        {
                            monitoringLineFaces_.insert((*cveit)->Face(n));
                            (*cveit)->Face(n)->Active(false);
                        }

                    }

                }

            }

        } // end looping over the faces

    } // end looping over 1d2d cves

    // looping over monitoring 0d cves to find monitoring faces
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = new_monitoring_0d_cves.begin();
        cveit != new_monitoring_0d_cves.end(); ++cveit)
    {
        // insert cve in the monitoring container with initial oil saturation as the value of the key
        monitoring0D_CVEs_[*cveit] = (*cveit)->PropertyValue(1U);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // find inactive faces and inserting them in the monitoring point face container
            if (!((*cveit)->Neighbor(n)->Active()))
            {
                monitoringPointFaces_.insert((*cveit)->Face(n));
            }
            else
            {
                if (!((*cveit)->Face(n)->Active()))
                {
                    (*cveit)->Face(n)->Active(true);
                    activeInteriorPointFaces_.push_back((*cveit)->Face(n));
                    monitoringPointFaces_.erase((*cveit)->Face(n));
                }

            }

        } // end looping over the faces

    } // end looping over 0d cves

    // calculating divergence and minimum CFL based time step for new active region
    double min_cfl_timestep(std::numeric_limits<double>::max());

    // looping over 2D and 1D cves
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = new_monitoring_1d2d_cves.begin();
        cveit != new_monitoring_1d2d_cves.end(); ++cveit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); ++i)
        {
            total_flux += std::abs((*cveit)->Face(i)->PropertyValue(0U));
            flux_imbalance += (*cveit)->Face(i)->PropertyValue(0U) * (*cveit)->FaceNormalDirection(i);
        }

        flux_imbalance += ((*cveit)->E()->Read(this->totalProdRateKey_) + 
                           (*cveit)->E()->Read(this->oilInjRateKey_) + 
                           (*cveit)->E()->Read(this->waterInjRateKey_)) * 
                           (*cveit)->Volume();

        (*cveit)->E()->Store(this->divergenceKey_, ScalarVariable(PLAIN, flux_imbalance));

        const double element_timestep((*cveit)->Volume() * (*cveit)->E()->Read(this->porosityKey_) * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;

    }
    // looping over 0D cves
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = new_monitoring_0d_cves.begin();
        cveit != new_monitoring_0d_cves.end(); ++cveit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); ++i)
        {
            total_flux += std::abs((*cveit)->Face(i)->PropertyValue(0U));
            flux_imbalance += (*cveit)->Face(i)->PropertyValue(0U) * (*cveit)->FaceNormalDirection(i);
        }

        const double element_timestep((*cveit)->Volume() * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;

    }

    return min_cfl_timestep;

} // end UpdateActiveElements





template<uint32_t dim>
void ActiveElementTwoPhaseTransport<dim>::CheckForInactiveCVEs(TwoPhaseModel<dim>& flow_model, const double activation_criteria)
{
    std::cout << "\nActiveElementTwoPhaseTransport<" << dim << ">::CheckForInactiveCVEs()\n";

    // create a copy of previous active and monitoring cves
    std::vector<ControlVolumeElement<dim>*> previous_active1D2D_CVEs(active1D2D_CVEs_);
    std::vector<ControlVolumeElement<dim>*> previous_active0D_CVEs(active0D_CVEs_);
    // we store the previous value for monitoring cves, so if they are still monitoring we assign that value to the map
    std::map<ControlVolumeElement<dim>*, double> previous_monitoring1D2D_CVEs(monitoring1D2D_CVEs_);
    std::map<ControlVolumeElement<dim>*, double> previous_monitoring0D_CVEs(monitoring0D_CVEs_);

    //for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = Monitoring1D2D_CVEsBegin();
    //    cveit != Monitoring1D2D_CVEsEnd(); ++cveit)
    //{
    //    previous_active1D2D_CVEs.push_back(cveit->first);
    //}

    //for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = Monitoring0D_CVEsBegin();
    //    cveit != Monitoring0D_CVEsEnd(); ++cveit)
    //{
    //    previous_active0D_CVEs.push_back(cveit->first);
    //}

    // changing the active monitoring cves flag to false
    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = Monitoring1D2D_CVEsBegin();
        cveit != Monitoring1D2D_CVEsEnd(); ++cveit)
    {
        cveit->first->Active(false);
    }

    for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = Monitoring0D_CVEsBegin();
        cveit != Monitoring0D_CVEsEnd(); ++cveit)
    {
        cveit->first->Active(false);
    }

    // changing the active face flag to false
    for (typename std::vector<ElementFace<dim>*>::iterator fit = ActiveInteriorLineFacesBegin(); 
        fit != ActiveInteriorLineFacesEnd(); ++fit)
    {
        (*fit)->Active(false);
    }

    for (typename std::vector<ElementFace<dim>*>::iterator fit = ActiveInteriorPointFacesBegin(); 
        fit != ActiveInteriorPointFacesEnd(); ++fit)
    {
        (*fit)->Active(false);
    }

    for (typename std::vector<ElementFace<dim>*>::iterator fit = ActiveModelBoundaryFacesBegin(); 
        fit != ActiveModelBoundaryFacesEnd(); ++fit)
    {
        (*fit)->Active(false);
    }

    // changing the previous active patch flag to false
    for (typename std::vector<CVE_Patch<dim>*>::iterator pit = ActiveInteriorPatchesBegin();
         pit != ActiveInteriorPatchesEnd(); ++pit)
    {
        (*pit)->Active(false);
    }

    for (typename std::vector<CVE_Patch<dim>*>::iterator pit = ActiveBoundaryPatchesBegin();
        pit != ActiveBoundaryPatchesEnd(); ++pit)
    {
        (*pit)->Active(false);
    }

    activeInteriorLineFaces_.resize(0U);
    activeInteriorPointFaces_.resize(0U);
    activeModelBoundaryFaces_.resize(0U);
    active1D2D_CVEs_.resize(0U);
    active0D_CVEs_.resize(0U);
    monitoring1D2D_CVEs_.clear();
    monitoring0D_CVEs_.clear();
    monitoringLineFaces_.clear();
    monitoringPointFaces_.clear();
    activeInteriorPatches_.resize(0U);
    activeBoundaryPatches_.resize(0U);

    // looping over previous active 2D and 1D cves to find new active cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = previous_active1D2D_CVEs.begin();
        cveit != previous_active1D2D_CVEs.end(); ++cveit)
    {
        (*cveit)->Active(false);

        double accumulation(0.);
        double divergence(0.);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            if (((*cveit)->Face(n)->InsideCVE()->Dimension() == 0U) || 
                ((*cveit)->Face(n)->OutsideCVE()->Dimension() == 0U))
            {
                accumulation += this->CalculateOilFluxInteriorPointFace(*((*cveit)->Face(n)), flow_model) * (*cveit)->FaceNormalDirection(n);
            } 
            else
            {
                accumulation += this->CalculateOilFluxInteriorLineFace(*((*cveit)->Face(n)), flow_model) * (*cveit)->FaceNormalDirection(n);
            }

            divergence += (*cveit)->Face(n)->PropertyValue(0U) * (*cveit)->FaceNormalDirection(n);

        } // end loop faces

        flow_model.Initialize(*((*cveit)->E()));
        flow_model.EffectiveSaturation();
        const double lambda_w_in(flow_model.MobilityPhase(1U));
        const double lambda_o_in(flow_model.MobilityPhase(2U));

        const double oil_source_sink((*cveit)->Volume() * ((*cveit)->E()->Read(this->oilInjRateKey_) +
            (*cveit)->E()->Read(this->totalProdRateKey_) * lambda_o_in / (lambda_w_in + lambda_o_in)));
        accumulation += oil_source_sink;
        accumulation -= divergence * lambda_o_in / (lambda_w_in + lambda_o_in);

        if (std::abs(accumulation) > activation_criteria)
        {
            active1D2D_CVEs_.push_back(*cveit);
            (*cveit)->Active(true);

            // making the parent patches active and insert in the corresponding container
            for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
                pit != (*cveit)->PatchesEnd(); ++pit)
            {
                if (!((*pit)->Active()))
                {
                    (*pit)->Active(true);

                    if ((*pit)->CenterNode()->AtBoundary() == NOT)
                    {
                        activeInteriorPatches_.push_back(*pit);
                    } 
                    else
                    {
                        activeBoundaryPatches_.push_back(*pit);
                    }

                } // end if patch is not active

            } // end looping over the patches

        }

    } // end looping over inside 1D and 2D cves

    // looping over previous active 0D cves to find new active cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = previous_active0D_CVEs.begin();
        cveit != previous_active0D_CVEs.end(); ++cveit)
    {
        (*cveit)->Active(false);

        double vis_w_avg(0.);
        double vis_o_avg(0.);

        for (size_t n = 0U; n < (*cveit)->Neighbors(); ++n)
        {
            flow_model.Initialize(*((*cveit)->Neighbor(n)->E()));
            flow_model.EffectiveSaturation();

            vis_w_avg += flow_model.ViscosityWettingPhase();
            vis_o_avg += flow_model.ViscosityNonWettingPhase();
        }

        vis_w_avg /= (*cveit)->Neighbors();
        vis_o_avg /= (*cveit)->Neighbors();

        double lambda_w_in;
        double lambda_o_in;

        if (((*cveit)->PropertyValue(1U) > 0.) && ((*cveit)->PropertyValue(1U) < 1.))
        {
            lambda_o_in = (*cveit)->PropertyValue(1U) / vis_o_avg;
            lambda_w_in = (1. - (*cveit)->PropertyValue(1U)) / vis_w_avg;                    
        }
        else if ((*cveit)->PropertyValue(1U) >= 1.)
        {
            lambda_o_in = 1. / vis_o_avg;
            lambda_w_in = 0.;
        }
        else
        {
            lambda_o_in = 0.;
            lambda_w_in = 1. / vis_w_avg;
        }

        double accumulation(0.);
        double divergence(0.);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            accumulation += this->CalculateOilFluxInteriorPointFace(*((*cveit)->Face(n)), flow_model) * (*cveit)->FaceNormalDirection(n);
            divergence += (*cveit)->Face(n)->PropertyValue(0U) * (*cveit)->FaceNormalDirection(n);
        } // end looping over neighbors

        accumulation -= divergence * lambda_o_in / (lambda_w_in + lambda_o_in);

        if (std::abs(accumulation) > activation_criteria)
        {
            active0D_CVEs_.push_back(*cveit);
            (*cveit)->Active(true);

            // making the parent patches active and insert in the corresponding container
            for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
                pit != (*cveit)->PatchesEnd(); ++pit)
            {
                if (!((*pit)->Active()))
                {
                    (*pit)->Active(true);
                    activeInteriorPatches_.push_back(*pit);
                } // end if patch is not active

            } // end looping over the patches

        }

    } // end looping over inside 0D cves

    std::multiset<ElementFace<dim>*> monitoring_line_faces;
    std::multiset<ElementFace<dim>*> monitoring_point_faces;
    std::set<ControlVolumeElement<dim>*> monitoring_1d2d_cves;
    std::set<ControlVolumeElement<dim>*> monitoring_0d_cves;

    // finding active faces and monitoring cves (looping over 1D-2D cves)
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active1D2D_CVEs_.begin();
        cveit != active1D2D_CVEs_.end(); ++cveit)
    {
        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            if ((*cveit)->Face(n)->Placement() == NOT)
            {
                // make the face active and inserting the face in the corresponding active face container
                if (!((*cveit)->Face(n)->Active()))
                {
                    (*cveit)->Face(n)->Active(true);

                    if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                    {
                        activeInteriorPointFaces_.push_back((*cveit)->Face(n));
                    } 
                    else
                    {
                        activeInteriorLineFaces_.push_back((*cveit)->Face(n));
                    }

                }
                // if the neighbor cve is not active put it in the monitoring container
                if (!((*cveit)->Neighbor(n)->Active()))
                {
                    if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                    {
                        monitoring_0d_cves.insert((*cveit)->Neighbor(n));
                    } 
                    else
                    {
                        monitoring_1d2d_cves.insert((*cveit)->Neighbor(n));
                    }

                }

            }
            // if it is a boundary face
            else
            {
                (*cveit)->Face(n)->Active(true);
                activeModelBoundaryFaces_.push_back((*cveit)->Face(n));
            }

        } // end looping over cve's faces

    } // end looping over 1d2d cves

    // finding active faces and monitoring cves (looping over 0D cves)
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active0D_CVEs_.begin();
        cveit != active0D_CVEs_.end(); ++cveit)
    {
        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // make the face active and inserting the face in the corresponding active face container
            if (!((*cveit)->Face(n)->Active()))
            {
                (*cveit)->Face(n)->Active(true);
                activeInteriorPointFaces_.push_back((*cveit)->Face(n));
            }
            // if the neighbor cve is not active put it in the monitoring container
            if (!((*cveit)->Neighbor(n)->Active()))
            {
                monitoring_1d2d_cves.insert((*cveit)->Neighbor(n));
            }

        } // end looping over cve's faces

    } // end looping over 0d cves

    // looping over monitoring 1d2d cves to find monitoring faces and make cve's parent patches active
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = monitoring_1d2d_cves.begin();
        cveit != monitoring_1d2d_cves.end(); ++cveit)
    {
        // we set the flag for monitoring cves as active
        (*cveit)->Active(true);

        // making the parent patches active and insert in the corresponding container
        for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
            pit != (*cveit)->PatchesEnd(); ++pit)
        {
            if (!((*pit)->Active()))
            {
                (*pit)->Active(true);

                if ((*pit)->CenterNode()->AtBoundary() == NOT)
                {
                    activeInteriorPatches_.push_back(*pit);
                } 
                else
                {
                    activeBoundaryPatches_.push_back(*pit);
                }

            } // end if patch is not active

        } // end looping over the patches

        // if cves was already monitoring cve we use the previous value as the key value, 
        // otherwise we use current saturation
        if (previous_monitoring1D2D_CVEs.find(*cveit) != previous_monitoring1D2D_CVEs.end())
        {
            monitoring1D2D_CVEs_[*cveit] = previous_monitoring1D2D_CVEs[*cveit];
        } 
        else
        {
            flow_model.Initialize(*((*cveit)->E()));
            monitoring1D2D_CVEs_[*cveit] = flow_model.Saturation(2U);
        }

        //// insert cve in the monitoring container with initial oil saturation as the value of the key
        //flow_model.Initialize(*((*cveit)->E()));
        //monitoring1D2D_CVEs_[*cveit] = flow_model.Saturation(2U);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // find inactive faces and inserting them in the corresponding monitoring face container
            if (!((*cveit)->Face(n)->Active()))
            {
                if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                {
                    monitoring_point_faces.insert((*cveit)->Face(n));
                } 
                else
                {
                    monitoring_line_faces.insert((*cveit)->Face(n));
                }

            }

        } // end looping over the faces

    } // end looping over 1d2d cves

    // looping over monitoring 0d cves to find monitoring faces and make cve's parent patches active
    for (typename std::set<ControlVolumeElement<dim>*>::iterator cveit = monitoring_0d_cves.begin();
        cveit != monitoring_0d_cves.end(); ++cveit)
    {
        // we set the flag for monitoring cves as active
        (*cveit)->Active(true);

        // making the parent patches active and insert in the corresponding container
        for (typename std::vector<CVE_Patch<dim>*>::iterator pit = (*cveit)->PatchesBegin(); 
            pit != (*cveit)->PatchesEnd(); ++pit)
        {
            if (!((*pit)->Active()))
            {
                (*pit)->Active(true);
                activeInteriorPatches_.push_back(*pit);
            } // end if patch is not active

        } // end looping over the patches

        // if cves was already monitoring cve we use the previous value as the key value, 
        // otherwise we use current saturation
        if (previous_monitoring0D_CVEs.find(*cveit) != previous_monitoring0D_CVEs.end())
        {
            monitoring0D_CVEs_[*cveit] = previous_monitoring0D_CVEs[*cveit];
        } 
        else
        {
            monitoring0D_CVEs_[*cveit] = (*cveit)->PropertyValue(1U);
        }

        //// insert cve in the monitoring container with initial oil saturation as the value of the key
        //monitoring0D_CVEs_[*cveit] = (*cveit)->PropertyValue(1U);

        for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
        {
            // find inactive faces and inserting them in the monitoring point face container
            if (!((*cveit)->Face(n)->Active()))
            {
                monitoring_point_faces.insert((*cveit)->Face(n));
            }

        } // end looping over the faces

    } // end looping over 0d cves

    // putting monitoring faces in the corresponding container
    // looping over monitoring line faces
    for (typename std::multiset<ElementFace<dim>*>::iterator fit = monitoring_line_faces.begin(); 
        fit != monitoring_line_faces.end(); ++fit)
    {
        // if the face is not shared between two monitoring cves
        if (monitoring_line_faces.count(*fit) == 1U)
        {
            monitoringLineFaces_.insert(*fit);
        }
        // if the face is shared between two monitoring cves we make it active
        else
        {
            if (!((*fit)->Active()))
            {
                (*fit)->Active(true);
                activeInteriorLineFaces_.push_back(*fit);
            }

        }

    } // end looping over monitoring line faces

    // looping over monitoring point faces
    for (typename std::multiset<ElementFace<dim>*>::iterator fit = monitoring_point_faces.begin(); 
        fit != monitoring_point_faces.end(); ++fit)
    {
        // if the face is not shared between two monitoring cves
        if (monitoring_point_faces.count(*fit) == 1U)
        {
            monitoringPointFaces_.insert(*fit);
        }
        // if the face is shared between two monitoring cves we make it active
        else
        {
            if (!((*fit)->Active()))
            {
                (*fit)->Active(true);
                activeInteriorPointFaces_.push_back(*fit);
            }

        }

    } // end looping over monitoring point faces

    //// looping over previous monitoring cves. if they are not included as monitoring cves in the new setup we add them
    //std::vector<ControlVolumeElement<dim>*> previos_monitoring1D2D_CVEs_needToBeAdded;
    //std::vector<ControlVolumeElement<dim>*> previos_monitoring0D_CVEs_needToBeAdded;
    //
    //// looping over previous monitoring 1d2d cves
    //for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = previous_monitoring1D2D_CVEs.begin();
    //     cveit != previous_monitoring1D2D_CVEs.end(); ++cveit)
    //{
    //    if (!(cveit->first->Active()))
    //    {
    //        // first we check if there is saturation changes in the cve. if there is no saturation changes
    //        // we don't need to include the cve as a monitoring cve
    //        double accumulation(0.);
    //        double divergence(0.);

    //        for (size_t n = 0U; n < cveit->first->Faces(); ++n)
    //        {
    //            if ((cveit->first->Face(n)->InsideCVE()->Dimension() == 0U) || 
    //                (cveit->first->Face(n)->OutsideCVE()->Dimension() == 0U))
    //            {
    //                accumulation += this->CalculateOilFluxInteriorPointFace(*(cveit->first->Face(n)), flow_model) * cveit->first->FaceNormalDirection(n);
    //            } 
    //            else
    //            {
    //                accumulation += this->CalculateOilFluxInteriorLineFace(*(cveit->first->Face(n)), flow_model) * cveit->first->FaceNormalDirection(n);
    //            }

    //            divergence += cveit->first->Face(n)->PropertyValue(0U) * cveit->first->FaceNormalDirection(n);

    //        } // end loop faces

    //        flow_model.Initialize(*(cveit->first->E()));
    //        flow_model.EffectiveSaturation();
    //        const double lambda_w_in(flow_model.MobilityPhase(1U));
    //        const double lambda_o_in(flow_model.MobilityPhase(2U));

    //        const double oil_source_sink(cveit->first->Volume() * (cveit->first->E()->Read(this->oilInjRateKey_) +
    //            cveit->first->E()->Read(this->totalProdRateKey_) * lambda_o_in / (lambda_w_in + lambda_o_in)));
    //        accumulation += oil_source_sink;
    //        accumulation -= divergence * lambda_o_in / (lambda_w_in + lambda_o_in);

    //        if (std::abs(accumulation) > activation_criteria)
    //        {
    //            cveit->first->Active(true);
    //            previos_monitoring1D2D_CVEs_needToBeAdded.push_back(cveit->first);
    //            monitoring1D2D_CVEs_[cveit->first] = cveit->second;

    //            // making the parent patches active and insert in the corresponding container
    //            for (typename std::vector<CVE_Patch<dim>*>::iterator pit = cveit->first->PatchesBegin(); 
    //                pit != cveit->first->PatchesEnd(); ++pit)
    //            {
    //                if (!((*pit)->Active()))
    //                {
    //                    (*pit)->Active(true);

    //                    if ((*pit)->CenterNode()->AtBoundary() == NOT)
    //                    {
    //                        activeInteriorPatches_.push_back(*pit);
    //                    } 
    //                    else
    //                    {
    //                        activeBoundaryPatches_.push_back(*pit);
    //                    }

    //                } // end if patch is not active

    //            } // end looping over the patches

    //        } // if there is saturation changes

    //    } // end if previous monitoring cve is not active

    //} // end looping over previous monitoring cves

    //// looping over previous monitoring 0d cves
    //for (typename std::map<ControlVolumeElement<dim>*, double>::iterator cveit = previous_monitoring0D_CVEs.begin();
    //    cveit != previous_monitoring0D_CVEs.end(); ++cveit)
    //{
    //    if (!(cveit->first->Active()))
    //    {
    //        // first we check if there is saturation changes in the cve. if there is no saturation changes
    //        // we don't need to include the cve as a monitoring cve
    //        double vis_w_avg(0.);
    //        double vis_o_avg(0.);

    //        for (size_t n = 0U; n < cveit->first->Neighbors(); ++n)
    //        {
    //            flow_model.Initialize(*(cveit->first->Neighbor(n)->E()));
    //            flow_model.EffectiveSaturation();

    //            vis_w_avg += flow_model.ViscosityWettingPhase();
    //            vis_o_avg += flow_model.ViscosityNonWettingPhase();
    //        }

    //        vis_w_avg /= cveit->first->Neighbors();
    //        vis_o_avg /= cveit->first->Neighbors();

    //        double lambda_w_in;
    //        double lambda_o_in;

    //        if ((cveit->first->PropertyValue(1U) > 0.) && (cveit->first->PropertyValue(1U) < 1.))
    //        {
    //            lambda_o_in = cveit->first->PropertyValue(1U) / vis_o_avg;
    //            lambda_w_in = (1. - cveit->first->PropertyValue(1U)) / vis_w_avg;                    
    //        }
    //        else if (cveit->first->PropertyValue(1U) >= 1.)
    //        {
    //            lambda_o_in = 1. / vis_o_avg;
    //            lambda_w_in = 0.;
    //        }
    //        else
    //        {
    //            lambda_o_in = 0.;
    //            lambda_w_in = 1. / vis_w_avg;
    //        }

    //        double accumulation(0.);
    //        double divergence(0.);

    //        for (size_t n = 0U; n < cveit->first->Faces(); ++n)
    //        {
    //            accumulation += this->CalculateOilFluxInteriorPointFace(*(cveit->first->Face(n)), flow_model) * cveit->first->FaceNormalDirection(n);
    //            divergence += cveit->first->Face(n)->PropertyValue(0U) * cveit->first->FaceNormalDirection(n);
    //        } // end looping over neighbors

    //        accumulation -= divergence * lambda_o_in / (lambda_w_in + lambda_o_in);

    //        if (std::abs(accumulation) > activation_criteria)
    //        {
    //            cveit->first->Active(true);
    //            previos_monitoring0D_CVEs_needToBeAdded.push_back(cveit->first);
    //            monitoring0D_CVEs_[cveit->first] = cveit->second;

    //            // making the parent patches active and insert in the corresponding container
    //            for (typename std::vector<CVE_Patch<dim>*>::iterator pit = cveit->first->PatchesBegin(); 
    //                pit != cveit->first->PatchesEnd(); ++pit)
    //            {
    //                if (!((*pit)->Active()))
    //                {
    //                    (*pit)->Active(true);
    //                    activeInteriorPatches_.push_back(*pit);
    //                } // end if patch is not active

    //            } // end looping over the patches

    //        } // end if there is saturation changes 

    //    } // end if previous monitoring cve is not active

    //} // end looping over previous monitoring cves

    //// looping over previous monitoring 1d2d cves those are need to be added to find new active and monitoring faces
    //for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = previos_monitoring1D2D_CVEs_needToBeAdded.begin();
    //     cveit != previos_monitoring1D2D_CVEs_needToBeAdded.end(); ++cveit)
    //{
    //    for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
    //    {
    //        if ((*cveit)->Face(n)->Placement() == NOT)
    //        {
    //            if ((*cveit)->Neighbor(n)->Active())
    //            {
    //                if (!((*cveit)->Face(n)->Active()))
    //                {
    //                    (*cveit)->Face(n)->Active(true);

    //                    if ((*cveit)->Neighbor(n)->Dimension() == 0U)
    //                    {
    //                        activeInteriorPointFaces_.push_back((*cveit)->Face(n));
    //                        monitoringPointFaces_.erase((*cveit)->Face(n));
    //                    } 
    //                    else
    //                    {
    //                        activeInteriorLineFaces_.push_back((*cveit)->Face(n));
    //                        monitoringLineFaces_.erase((*cveit)->Face(n));
    //                    }

    //                }

    //            } // end if neighbor is active cve
    //            // if neighbor is not active
    //            else
    //            {
    //                if ((*cveit)->Neighbor(n)->Dimension() == 0U)
    //                {
    //                    monitoringPointFaces_.insert((*cveit)->Face(n));
    //                } 
    //                else
    //                {
    //                    monitoringLineFaces_.insert((*cveit)->Face(n));
    //                }

    //            } // end if neighbor is not active

    //        } // end if it is an interior face
    //        else
    //        {
    //            if (!((*cveit)->Face(n)->Active()))
    //            {
    //                monitoringLineFaces_.insert((*cveit)->Face(n));
    //            }

    //        } // end if it is a boundary face

    //    } // end looping over the faces

    //} // end looping over previous monitoring 1d2d cves those are need to be added

    //// looping over previous monitoring 0d cves those are need to be added to find new active and monitoring faces
    //for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = previos_monitoring0D_CVEs_needToBeAdded.begin();
    //     cveit != previos_monitoring0D_CVEs_needToBeAdded.end(); ++cveit)
    //{
    //    for (size_t n = 0U; n < (*cveit)->Faces(); ++n)
    //    {
    //        if ((*cveit)->Neighbor(n)->Active())
    //        {
    //            if (!((*cveit)->Face(n)->Active()))
    //            {
    //                (*cveit)->Face(n)->Active(true);
    //                activeInteriorPointFaces_.push_back((*cveit)->Face(n));
    //                monitoringPointFaces_.erase((*cveit)->Face(n));
    //            }

    //        } // end if neighbor is active cve
    //        // if neighbor is not active
    //        else
    //        {
    //                monitoringPointFaces_.insert((*cveit)->Face(n));
    //        } // end if neighbor is not active

    //    } // end looping over cve's faces

    //} // end looping over previous monitoring 1d2d cves those are need to be added


} // end CheckForInactiveCVEs





template<uint32_t dim>
double ActiveElementTwoPhaseTransport<dim>::Transport(double total_time_increment, TwoPhaseModel<dim>& flow_model)
{
    std::cout << "\nActiveElementTwoPhaseTransport<" << dim << ">::Transport()\n";

    typedef typename std::vector<ElementFace<dim>*>::iterator face_ptr_itr;
    typedef typename std::vector<ControlVolumeElement<dim>*>::iterator CVE_ptr_itr;
    typedef typename std::map<ControlVolumeElement<dim>*, double>::iterator m_cve_ptr_itr;
    typedef typename std::set<ElementFace<dim>*>::iterator m_face_ptr_itr;

    static bool first_time(true);
    static size_t deactivate_counter(1U);

    if (first_time)
    {
        this->ConstructFluxes();
        FindActiveRegions(flow_model, this->setup_.ActivationCriteria());
        ConstructFluxesActiveRegions();
        first_time = false;
    } 
    else
    {
        if (this->setup_.Deactivate() != 0U)
        {
            if (deactivate_counter % this->setup_.Deactivate() == 0U)
            {
                ConstructFluxesActiveRegions();
                CheckForInactiveCVEs(flow_model, this->setup_.ActivationCriteria());
                ConstructFluxesActiveRegions();
                deactivate_counter = 1U;
            }
            else
            {
                ConstructFluxesActiveRegions();
                ++deactivate_counter;
            }

        }
        else
        {
            ConstructFluxesActiveRegions();
        }
        
    }

    double normal_changes(this->setup_.NormalSaturationChanges());
    double max_changes(this->setup_.MaximumSaturationChanges());
    double max_internal_time_step_size(this->setup_.MaximumTransportTimestepSize());
    size_t maximum_number_of_transport_steps(this->setup_.MaximumTransportTimesteps());
    double CFL_multiplier(this->setup_.CFL_Multiplier());
    bool with_divergence_correction(this->setup_.WithDivergenceCorrection());

    std::pair<double, double> CFL_timestep_divergence(CFL_Timestep_MaxDivergenceActiveRegions());
    double CFL_timestep(CFL_timestep_divergence.first * CFL_multiplier);
    std::cout << std::scientific;
    std::cout.precision(4);
    std::cout << std::fixed;
    std::cout << "\nCFL Time step Criteria: " << CFL_timestep << " sec\n";
    std::cout << "Max Divergence: " << CFL_timestep_divergence.second << "\n";
    std::cout << std::fixed;
    std::cout << "Transport Interval: " << total_time_increment << " sec\n";
    std::cout << "Transport Maximum Number of Steps: " << maximum_number_of_transport_steps << "\n";
    std::cout << "Active Percentage: " << std:: fixed << (double) (active1D2D_CVEs_.size() + active0D_CVEs_.size() + monitoring1D2D_CVEs_.size() + monitoring0D_CVEs_.size()) /
                                           (double) this->transportModelRef_.CVEs() * 100. << " %\n";
    std::cout << "\n";

    static double internal_time_increment(total_time_increment);
    if (internal_time_increment > CFL_timestep) internal_time_increment = CFL_timestep;
    if (internal_time_increment > total_time_increment) internal_time_increment = total_time_increment;
    double total_time(0.);
    double previous_time_step(0.);
    double maximum_change_in_timestep(0.);
    size_t transport_step_number(0U);

    SignalHandler sig;

    bool exiting_signal_received(false);

    // outer time loop for achieving total_time_increment
    do
    {
        // inner time loop for controlling the stability based on normal_changes and max_changes
        do
        {
            ZeroActiveAndMonitoringCVEsProperty(0U);

            // transport for inner line faces
            for (face_ptr_itr fit = ActiveInteriorLineFacesBegin(); fit != ActiveInteriorLineFacesEnd(); ++fit)
            {                
                const double oil_face_flux(this->CalculateOilFluxInteriorLineFace(*(*fit), flow_model));
                (*fit)->InsideCVE()->PropertyValueAdd(0U, oil_face_flux);
                (*fit)->OutsideCVE()->PropertyValueAdd(0U, -oil_face_flux);

            } // end for loop transport for inner line faces

            // transport for boundary line faces
            for (face_ptr_itr fit = ActiveModelBoundaryFacesBegin(); fit != ActiveModelBoundaryFacesEnd(); ++fit)
            {                
                const double oil_face_flux(this->CalculateOilFluxBoundaryLineFace(*(*fit), flow_model));
                (*fit)->InsideCVE()->PropertyValueAdd(0U, oil_face_flux);

            } // end for loop transport for boundary line faces

            //// looping over interior point CVEs
            //for (CVE_ptr_itr cveit = Active0D_CVEsBegin(); cveit != Active0D_CVEsEnd(); ++cveit)
            //{
            //    this->TransportAcrossZeroDimensionalCVE_Faces(*(*cveit), flow_model);
            //} // end looping over interior point CVEs


            // transport through point faces
            // looping over interior active point faces
            for (face_ptr_itr fit = ActiveInteriorPointFacesBegin(); fit != ActiveInteriorPointFacesEnd(); fit++)
            {                
                const double oil_face_flux(this->CalculateOilFluxInteriorPointFace(*(*fit), flow_model));
                (*fit)->InsideCVE()->PropertyValueAdd(0U, oil_face_flux);
                (*fit)->OutsideCVE()->PropertyValueAdd(0U, -oil_face_flux);

            } // end for loop transport for inner line faces

            // we replace the current saturation of monitoring cves with initial value to do transport over monitoring faces
            // for 1d2d cves
            for (m_cve_ptr_itr cveit = Monitoring1D2D_CVEsBegin(); cveit != Monitoring1D2D_CVEsEnd(); ++cveit)
            {
                const double oil_saturation(cveit->first->E()->Read(this->oilSaturationKey_));
                cveit->first->E()->Store(this->oilSaturationKey_, ScalarVariable(cveit->first->E()->Status(this->oilSaturationKey_), 
                    cveit->second));
                cveit->first->E()->Store(this->waterSaturationKey_, ScalarVariable(cveit->first->E()->Status(this->waterSaturationKey_), 
                    1. - cveit->second));

                monitoring1D2D_CVEs_[cveit->first] = oil_saturation;
            }

            // for 0d cves
            for (m_cve_ptr_itr cveit = Monitoring0D_CVEsBegin(); cveit != Monitoring0D_CVEsEnd(); ++cveit)
            {
                const double oil_saturation(cveit->first->PropertyValue(1U));
                cveit->first->PropertyValue(1U, cveit->second);
                monitoring0D_CVEs_[cveit->first] = oil_saturation;
            }


            // transport over monitoring line faces
            for (m_face_ptr_itr fit = MonitoringLineFacesBegin(); fit != MonitoringLineFacesEnd(); ++fit)
            {
                const double oil_flux((*fit)->Placement() == NOT ? this->CalculateOilFluxInteriorLineFace(*(*fit), flow_model) : 
                                        this->CalculateOilFluxBoundaryLineFace(*(*fit), flow_model));

                if ((*fit)->InsideCVE()->Active())
                {
                    (*fit)->InsideCVE()->PropertyValueAdd(0U, oil_flux);
                } 
                else
                {
                    (*fit)->OutsideCVE()->PropertyValueAdd(0U, -oil_flux);
                }
            }

            // transport over monitoring point faces
            for (m_face_ptr_itr fit = MonitoringPointFacesBegin(); fit != MonitoringPointFacesEnd(); ++fit)
            {
                const double oil_flux(this->CalculateOilFluxInteriorPointFace(*(*fit), flow_model));

                if ((*fit)->InsideCVE()->Active())
                {
                    (*fit)->InsideCVE()->PropertyValueAdd(0U, oil_flux);
                } 
                else
                {
                    (*fit)->OutsideCVE()->PropertyValueAdd(0U, -oil_flux);
                }
            }

            // now we replace back the saturation of monitoring cves with current value
            // for 1d2d cves
            for (m_cve_ptr_itr cveit = Monitoring1D2D_CVEsBegin(); cveit != Monitoring1D2D_CVEsEnd(); ++cveit)
            {
                const double oil_saturation(cveit->first->E()->Read(this->oilSaturationKey_));
                cveit->first->E()->Store(this->oilSaturationKey_, ScalarVariable(cveit->first->E()->Status(this->oilSaturationKey_), 
                    cveit->second));
                cveit->first->E()->Store(this->waterSaturationKey_, ScalarVariable(cveit->first->E()->Status(this->waterSaturationKey_), 
                    1. - cveit->second));

                monitoring1D2D_CVEs_[cveit->first] = oil_saturation;
            }

            // for 0d cves
            for (m_cve_ptr_itr cveit = Monitoring0D_CVEsBegin(); cveit != Monitoring0D_CVEsEnd(); ++cveit)
            {
                const double oil_saturation(cveit->first->PropertyValue(1U));
                cveit->first->PropertyValue(1U, cveit->second);
                monitoring0D_CVEs_[cveit->first] = oil_saturation;
            }

            maximum_change_in_timestep = 0.;

            // looping over active inside 2D and 1D cves
            for (CVE_ptr_itr cveit = Active1D2D_CVEsBegin(); cveit != Active1D2D_CVEsEnd(); ++cveit)
            {
                flow_model.Initialize(*((*cveit)->E()));
                flow_model.EffectiveSaturation();
                const double lambda_w(flow_model.MobilityPhase(1U));
                const double lambda_o(flow_model.MobilityPhase(2U));

                const double oil_source_sink((*cveit)->Volume() * ((*cveit)->E()->Read(this->oilInjRateKey_) +
                                               (*cveit)->E()->Read(this->totalProdRateKey_) * 
                                               lambda_o / (lambda_w + lambda_o)));

                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes((*cveit)->PropertyValue(0U) + oil_source_sink);
                // adding flux mismatch
                if (with_divergence_correction)
                    changes -= (*cveit)->E()->Read(this->divergenceKey_) * lambda_o / (lambda_w + lambda_o);
                // calculating element saturation changes for the timestep
                changes *= (internal_time_increment / (*cveit)->Volume() / (*cveit)->E()->Read(this->porosityKey_));

                (*cveit)->PropertyValue(0U, changes);

                maximum_change_in_timestep = ((std::abs(changes) > maximum_change_in_timestep)
                    ? std::abs(changes) : maximum_change_in_timestep);

            }

            // looping over active inside 0D cves
            for (CVE_ptr_itr cveit = Active0D_CVEsBegin(); cveit != Active0D_CVEsEnd(); ++cveit)
            {
                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes((*cveit)->PropertyValue(0U));
                // assuming porosity is 1.
                changes *= (internal_time_increment / (*cveit)->Volume());

                (*cveit)->PropertyValue(0U, changes);

                maximum_change_in_timestep = ((std::abs(changes) > maximum_change_in_timestep)
                    ? std::abs(changes) : maximum_change_in_timestep);

            }

            // looping over monitoring 2D and 1D cves
            for (m_cve_ptr_itr cveit = Monitoring1D2D_CVEsBegin(); cveit != Monitoring1D2D_CVEsEnd(); ++cveit)
            {
                flow_model.Initialize(*(cveit->first->E()));
                flow_model.EffectiveSaturation();
                const double lambda_w(flow_model.MobilityPhase(1U));
                const double lambda_o(flow_model.MobilityPhase(2U));

                const double oil_source_sink(cveit->first->Volume() * (cveit->first->E()->Read(this->oilInjRateKey_) +
                                               cveit->first->E()->Read(this->totalProdRateKey_) * 
                                               lambda_o / (lambda_w + lambda_o)));

                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes(cveit->first->PropertyValue(0U) + oil_source_sink);
                // adding flux mismatch
                if (with_divergence_correction)
                    changes -= cveit->first->E()->Read(this->divergenceKey_) * lambda_o / (lambda_w + lambda_o);
                // calculating element saturation changes for the timestep
                changes *= (internal_time_increment / cveit->first->Volume() / cveit->first->E()->Read(this->porosityKey_));

                cveit->first->PropertyValue(0U, changes);

                maximum_change_in_timestep = ((std::abs(changes) > maximum_change_in_timestep)
                    ? std::abs(changes) : maximum_change_in_timestep);

            }

            // looping over active inside 0D cves
            for (m_cve_ptr_itr cveit = Monitoring0D_CVEsBegin(); cveit != Monitoring0D_CVEsEnd(); ++cveit)
            {
                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes(cveit->first->PropertyValue(0U));
                // assuming porosity is 1.
                changes *= (internal_time_increment / cveit->first->Volume());

                cveit->first->PropertyValue(0U, changes);

                maximum_change_in_timestep = ((std::abs(changes) > maximum_change_in_timestep)
                    ? std::abs(changes) : maximum_change_in_timestep);

            }

            // store the time step size for calculating total time if do loop criteria satisfied
            previous_time_step = internal_time_increment;
            // calculating internal time step size based on maximum changes in the time step and normal changes
            internal_time_increment = internal_time_increment * normal_changes / maximum_change_in_timestep;
            // comparing to the CFL condition time step size
            internal_time_increment = (internal_time_increment > CFL_timestep ? CFL_timestep : internal_time_increment);

            internal_time_increment = ((internal_time_increment > max_internal_time_step_size)
                ? max_internal_time_step_size : internal_time_increment);

            std::cout << "Internal time step = " << previous_time_step << "\t\t";
            std::cout << "Max changes = " << maximum_change_in_timestep;
            if (maximum_change_in_timestep > max_changes)
                std::cout << "\t(CUT)\n";
            else
                std::cout << "\n";

        } while (maximum_change_in_timestep > max_changes); // for internal time do loop


        // calculating interior active 1D-2D cves saturation
        for (CVE_ptr_itr cveit = Active1D2D_CVEsBegin(); cveit != Active1D2D_CVEsEnd(); ++cveit)
        {
            if ((*cveit)->E()->Status(this->oilSaturationKey_) != DIRICH)
            {
                double saturation((*cveit)->E()->Read(this->oilSaturationKey_));
                saturation += (*cveit)->PropertyValue(0U);
                (*cveit)->E()->Store(this->oilSaturationKey_, ScalarVariable(PLAIN, saturation));
                (*cveit)->E()->Store(this->waterSaturationKey_, ScalarVariable(PLAIN, 1. - saturation));
            }

        }

        // calculating interior active 0D cves saturation
        for (CVE_ptr_itr cveit = Active0D_CVEsBegin(); cveit != Active0D_CVEsEnd(); ++cveit)
        {
            (*cveit)->PropertyValue(1U, (*cveit)->PropertyValue(1U) + (*cveit)->PropertyValue(0U));
        }

        // calculating monitoring 1D-2D cves saturation
        for (m_cve_ptr_itr cveit = Monitoring1D2D_CVEsBegin(); cveit != Monitoring1D2D_CVEsEnd(); ++cveit)
        {
            if (cveit->first->E()->Status(this->oilSaturationKey_) != DIRICH)
            {
                double saturation(cveit->first->E()->Read(this->oilSaturationKey_));
                saturation += cveit->first->PropertyValue(0U);
                cveit->first->E()->Store(this->oilSaturationKey_, ScalarVariable(PLAIN, saturation));
                cveit->first->E()->Store(this->waterSaturationKey_, ScalarVariable(PLAIN, 1. - saturation));
            }

        }

        // calculating monitoring 0D cves saturation
        for (m_cve_ptr_itr cveit = Monitoring0D_CVEsBegin(); cveit != Monitoring0D_CVEsEnd(); ++cveit)
        {
            cveit->first->PropertyValue(1U, cveit->first->PropertyValue(1U) + cveit->first->PropertyValue(0U));
        }
        const double new_cfl_timestep(UpdateActiveRegions(flow_model, this->setup_.Limiter()) * CFL_multiplier);
        if (CFL_timestep > new_cfl_timestep) CFL_timestep = new_cfl_timestep;

        total_time += previous_time_step;
        ++transport_step_number;



/*

PropertyHandle<2U> active_element(this->modelRef_, "active element", SCALAR, ELEMENT);
Index active_element_key(this->modelRef_.Database().StorageKey("active element"));


for (std::vector<ControlVolumeElement<2U> >::iterator cveit = this->transportModelRef_.CVEsBegin();
    cveit != this->transportModelRef_.PointCVEsBegin(); cveit++)
{
    if (cveit->Active())
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 1.));
    }
    else
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 0.));
    }
}

for (std::vector<ControlVolumeElement<2U> >::iterator cveit = this->transportModelRef_.PerimeterCVEsBegin();
    cveit != this->transportModelRef_.PerimeterPointCVEsBegin(); cveit++)
{
    if (cveit->Active())
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 1.));
    }
    else
    {
        cveit->E()->Store(active_element_key, ScalarVariable(PLAIN, 0.));
    }
}

VTK_Interface<2U> vtkOutput;

vtkOutput.OutputDataToVTK(this->modelRef_, "active element", "active element", ModelTime::Instance().modelTime + total_time, false);
vtkOutput.OutputDataToVTK(this->modelRef_, "saturation oil", "saturation oil", ModelTime::Instance().modelTime + total_time, false);


*/



        if ((internal_time_increment > (total_time_increment - total_time)) && ((total_time_increment - total_time) > 0.))
            internal_time_increment = (total_time_increment - total_time);

        if (internal_time_increment > max_internal_time_step_size) internal_time_increment = max_internal_time_step_size;
        if (internal_time_increment > CFL_timestep) internal_time_increment = CFL_timestep;

        if (sig.SignalRaised()) std::raise(SIGINT);

        if (sig.OutputSignal() || sig.RestartSignal() || sig.QuitSignal()) exiting_signal_received = true;

    } while ((total_time < total_time_increment) && (transport_step_number < maximum_number_of_transport_steps) && !exiting_signal_received); // end of outer do loop to achieve total_time_increment 

    std::cout << std::scientific;

    return total_time;

} // end Transport






template class ActiveElementTwoPhaseTransport<2U>;

}
