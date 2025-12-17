#include "ActiveFaceTwoPhaseTransport.h"
#include "Model.h"
#include "TransportModel.h"
#include "TwoPhaseModel.h"
#include "SignalHandler.h"
#include "IMPES_Setup.h"
#include "VariableOperations.h"

namespace csmp {



template<uint32_t dim>
ActiveFaceTwoPhaseTransport<dim>::ActiveFaceTwoPhaseTransport(Model<dim>& model,
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
ActiveFaceTwoPhaseTransport<dim>::~ActiveFaceTwoPhaseTransport()
{
}



template<uint32_t dim>
std::pair<double, double> ActiveFaceTwoPhaseTransport<dim>::CFL_Timestep_MaxDivergenceActiveRegions()
{
    double min_cfl_timestep(std::numeric_limits<double>::max());
    double max_divergence(0.);

    // looping over 2D and 1D cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active1D2D_CVEs_.begin();
        cveit != active1D2D_CVEs_.end(); cveit++)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); i++)
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

    // looping over 0D cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active0D_CVEs_.begin();
        cveit != active0D_CVEs_.end(); cveit++)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); i++)
        {
            total_flux += std::abs((*cveit)->Face(i)->PropertyValue(0U));
            flux_imbalance += (*cveit)->Face(i)->PropertyValue(0U) * (*cveit)->FaceNormalDirection(i);
        }

        const double element_timestep((*cveit)->Volume() * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);

    }

    return std::make_pair(min_cfl_timestep, max_divergence);

} // end CFL_Timestep_MaxDivergenceActiveRegions


/*

// This is the modified method which uses the square root of permeability instead of permeability itself
// to increase the flow inside low permeability zones


template<uint32_t dim>
void ActiveFaceTwoPhaseTransport<dim>::ConstructFluxesInteriorPatch(CVE_Patch<dim>& patch)
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
                    perm_inside = std::sqrt(patch.SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(this->permeabilityKey_)); // ---------------->
                    perm_outside = std::sqrt(patch.SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(this->permeabilityKey_)); // --------------->
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
void ActiveFaceTwoPhaseTransport<dim>::ConstructFluxesBoundaryPatch(CVE_Patch<dim>& patch)
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
    double perm_boundary_face_1(std::sqrt(patch.FaceAccessor(patch.InteriorFaces())->InsideCVE()->E()->Read(this->permeabilityKey_))); // ------------->
    double area_boundary_face_1(patch.FaceAccessor(patch.InteriorFaces())->Area() / 2.);

    for (size_t f = patch.InteriorFaces() + 1U; f < patch.Faces(); ++f)
    {
        if (patch.FaceAccessor(f)->InsideCVE()->Dimension() != 0U)
        {
            double perm_boundary_face_2(std::sqrt(patch.FaceAccessor(f)->InsideCVE()->E()->Read(this->permeabilityKey_))); // --------------->
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
void ActiveFaceTwoPhaseTransport<dim>::ConstructFluxesInteriorPatch(CVE_Patch<dim>& patch)
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
            for (uint32_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
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
void ActiveFaceTwoPhaseTransport<dim>::ConstructFluxesBoundaryPatch(CVE_Patch<dim>& patch)
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

            for (uint32_t f = 0U; f < patch.SectorAccessor(s)->Faces(); ++f)
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
void ActiveFaceTwoPhaseTransport<dim>::ConstructFluxesActiveRegions()
{
    std::cout << "\nActiveFaceTwoPhaseTransport<" << dim << ">::ConstructFluxesActiveRegions()\n";

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
void ActiveFaceTwoPhaseTransport<dim>::ZeroActiveCVEsProperty(size_t property_number)
{
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active1D2D_CVEs_.begin(); cveit != active1D2D_CVEs_.end(); cveit++)
    {
        (*cveit)->PropertyValue(property_number, 0.);
    }

    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active0D_CVEs_.begin(); cveit != active0D_CVEs_.end(); cveit++)
    {
        (*cveit)->PropertyValue(property_number, 0.);
    }

} // end ZeroActiveCVEsProperty





template<uint32_t dim>
void ActiveFaceTwoPhaseTransport<dim>::FindActiveRegions(TwoPhaseModel<dim>& flow_model, const double activation_criteria)
{
    std::cout << "\nActiveFaceTwoPhaseTransport<" << dim << ">::FindActiveRegions()\n";

    activeInteriorLineFaces_.resize(0U);
    activeInteriorPointFaces_.resize(0U);
    activeModelBoundaryFaces_.resize(0U);
    active1D2D_CVEs_.resize(0U);
    active0D_CVEs_.resize(0U);
    monitoringLineFacesFlux_.clear();
    monitoringPointFacesFlux_.clear();
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

        }

    } // end looping over inside 1D and 2D cves

    // looping over inside 0D cves to find active cves
    for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = this->transportModelRef_.PointCVEsBegin();
        cveit != this->transportModelRef_.PerimeterCVEsBegin(); ++cveit)
    {
        cveit->Active(false);

        double vis_w_avg(0.);
        double vis_o_avg(0.);

        for (size_t n = 0U; n < cveit->Neighbors(); n++)
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

        for (size_t n = 0U; n < cveit->Faces(); n++)
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
    
    // finding active and monitoring faces (looping over 1D-2D cves)
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active1D2D_CVEs_.begin();
         cveit != active1D2D_CVEs_.end(); cveit++)
    {
        for (size_t n = 0U; n < (*cveit)->Neighbors(); n++)
        {
            if ((*cveit)->Neighbor(n) != (*cveit))
            {
                if ((*cveit)->Neighbor(n)->Active())
                {
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

                }
                // if neighbor is not an active cve
                else
                {
                    if ((*cveit)->Neighbor(n)->Dimension() == 0U)
                    {
                        const double oil_flux(this->CalculateOilFluxInteriorPointFace(*((*cveit)->Face(n)), flow_model));
                        const double water_flux((*cveit)->Face(n)->PropertyValue(0U) - oil_flux);
                        const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                                                      oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

                        monitoringPointFacesFlux_[(*cveit)->Face(n)] = oil_flux_ratio;

                    } // end if neighbor is a 0D cve
                    // if neighbor is 1D or 2D cve
                    else
                    {
                        const double oil_flux(this->CalculateOilFluxInteriorLineFace(*((*cveit)->Face(n)), flow_model));
                        const double water_flux((*cveit)->Face(n)->PropertyValue(0U) - oil_flux);
                        const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                                                      oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

                        monitoringLineFacesFlux_[(*cveit)->Face(n)] = oil_flux_ratio;

                    } // end if neighbor is not a 0D cve

                } // end if neighbor is not an active cve

            }
            // if it is a boundary face
            else
            {
                (*cveit)->Face(n)->Active(true);
                activeModelBoundaryFaces_.push_back((*cveit)->Face(n));
            }
        }
    } // end looping over 1D 2D cves

    // finding active and monitoring faces (looping over 0D cves)
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = active0D_CVEs_.begin();
        cveit != active0D_CVEs_.end(); cveit++)
    {
        for (size_t n = 0U; n < (*cveit)->Neighbors(); n++)
        {
            if ((*cveit)->Neighbor(n)->Active())
            {
                if (!((*cveit)->Face(n)->Active()))
                {
                    (*cveit)->Face(n)->Active(true);
                    activeInteriorPointFaces_.push_back((*cveit)->Face(n));
                }
            }
            else
            {
                const double oil_flux(this->CalculateOilFluxInteriorPointFace(*((*cveit)->Face(n)), flow_model));
                const double water_flux((*cveit)->Face(n)->PropertyValue(0U) - oil_flux);
                const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                                               oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

                monitoringPointFacesFlux_[(*cveit)->Face(n)] = oil_flux_ratio;

            }

        }

    }

} // end FindActiveRegions




template<uint32_t dim>
double ActiveFaceTwoPhaseTransport<dim>::UpdateActiveRegions(TwoPhaseModel<dim>& flow_model, double criteria)
{
    std::vector<ControlVolumeElement<dim>*> new_active_1d2d_cves;
    std::vector<ControlVolumeElement<dim>*> new_active_0d_cves;
    std::vector<ElementFace<dim>*> new_active_line_faces;
    std::vector<ElementFace<dim>*> new_active_boundary_faces;
    std::vector<ElementFace<dim>*> new_active_point_faces;
    std::vector<CVE_Patch<dim>*> new_active_interior_patches;
    std::vector<CVE_Patch<dim>*> new_active_boundary_patches;

    // finding new active line faces
    for (typename std::map<ElementFace<dim>*, double>::iterator fit = monitoringLineFacesFlux_.begin();
        fit != monitoringLineFacesFlux_.end(); ++fit)
    {
        const double oil_flux(this->CalculateOilFluxInteriorLineFace(*(fit->first), flow_model));
        const double water_flux(fit->first->PropertyValue(0U) - oil_flux);
        const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                                      oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

        if (std::abs(oil_flux_ratio - fit->second) > criteria)
        {
            fit->first->Active(true);

            // adding the face to the active faces vector
            // since we never monitor boundary faces all faces in the monitoring region are interior faces
            // (if the parent element is active, they will be active as well)
            activeInteriorLineFaces_.push_back(fit->first);
            new_active_line_faces.push_back(fit->first);

            // finding new active cves
            if (!(fit->first->InsideCVE()->Active()))
            {
                fit->first->InsideCVE()->Active(true);
                active1D2D_CVEs_.push_back(fit->first->InsideCVE());
                new_active_1d2d_cves.push_back(fit->first->InsideCVE());
            }
            if (!(fit->first->OutsideCVE()->Active()))
            {
                fit->first->OutsideCVE()->Active(true);
                active1D2D_CVEs_.push_back(fit->first->OutsideCVE());
                new_active_1d2d_cves.push_back(fit->first->OutsideCVE());
            }

        }

    } // end for loop over active region boundary line faces

    // removing new active line faces from monitoring active region boundary line faces
    for (size_t i = 0U; i < new_active_line_faces.size(); i++)
    {
        monitoringLineFacesFlux_.erase(new_active_line_faces[i]);
    }

    for (size_t i = 0U; i < new_active_boundary_faces.size(); i++)
    {
        monitoringLineFacesFlux_.erase(new_active_boundary_faces[i]);
    }

    // finding new active point faces
    for (typename std::map<ElementFace<dim>*, double>::iterator fit = monitoringPointFacesFlux_.begin(); 
        fit != monitoringPointFacesFlux_.end(); fit++)
    {
        const double oil_flux(this->CalculateOilFluxInteriorPointFace(*(fit->first), flow_model));
        const double water_flux(fit->first->PropertyValue(0U) - oil_flux);
        const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                                       oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

        if (std::abs(oil_flux_ratio - fit->second) > criteria)
        {
            // adding the face to the active faces vector
            activeInteriorPointFaces_.push_back(fit->first);
            new_active_point_faces.push_back(fit->first);
            fit->first->Active(true);

            // finding new active cves
            if (!(fit->first->InsideCVE()->Active()))
            {
                fit->first->InsideCVE()->Active(true);

                if (fit->first->InsideCVE()->Dimension() == 0U)
                {
                    active0D_CVEs_.push_back(fit->first->InsideCVE());
                    new_active_0d_cves.push_back(fit->first->InsideCVE());
                }
                else
                {
                    active1D2D_CVEs_.push_back(fit->first->InsideCVE());
                    new_active_1d2d_cves.push_back(fit->first->InsideCVE());
                }
            }
            if (!(fit->first->OutsideCVE()->Active()))
            {
                fit->first->OutsideCVE()->Active(true);

                if (fit->first->OutsideCVE()->Dimension() == 0U)
                {
                    active0D_CVEs_.push_back(fit->first->OutsideCVE());
                    new_active_0d_cves.push_back(fit->first->OutsideCVE());
                }
                else
                {
                    active1D2D_CVEs_.push_back(fit->first->OutsideCVE());
                    new_active_1d2d_cves.push_back(fit->first->OutsideCVE());
                }

            }

        }

    } // end for loop over active region boundary point faces

    // removing new active point faces from monitoring active region boundary point faces
    for (size_t i = 0U; i < new_active_point_faces.size(); i++)
    {
        monitoringPointFacesFlux_.erase(new_active_point_faces[i]);
    }

    // finding new monitoring boundary faces and adding to the corresponding map
    // first looping over new active 1D-2D cves
    for (size_t i = 0U; i < new_active_1d2d_cves.size(); i++)
    {
        for (size_t n = 0U; n < new_active_1d2d_cves[i]->Neighbors(); n++)
        {
            // if it is a boundary face we activate it and don't monitor it
            if (new_active_1d2d_cves[i]->Face(n)->Placement() != NOT)
            {
                if (!(new_active_1d2d_cves[i]->Face(n)->Active()))
                {
                    new_active_1d2d_cves[i]->Face(n)->Active(true);
                    activeModelBoundaryFaces_.push_back(new_active_1d2d_cves[i]->Face(n));
                }
            } 
            else
            {
                if (!(new_active_1d2d_cves[i]->Neighbor(n)->Active()))
                {
                    // 0D cve neighbor
                    if (new_active_1d2d_cves[i]->Neighbor(n)->Dimension() == 0U)
                    {
                        const double oil_flux(this->CalculateOilFluxInteriorPointFace(*(new_active_1d2d_cves[i]->Face(n)), flow_model));
                        const double water_flux(new_active_1d2d_cves[i]->Face(n)->PropertyValue(0U) - oil_flux);
                        const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                                                      oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

                        monitoringPointFacesFlux_[new_active_1d2d_cves[i]->Face(n)] = oil_flux_ratio;

                    } 
                    // 1D-2D cve neighbor
                    else
                    {
                        const double oil_flux(this->CalculateOilFluxInteriorLineFace(*(new_active_1d2d_cves[i]->Face(n)), flow_model));
                        const double water_flux(new_active_1d2d_cves[i]->Face(n)->PropertyValue(0U) - oil_flux);
                        const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                            oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

                        monitoringLineFacesFlux_[new_active_1d2d_cves[i]->Face(n)] = oil_flux_ratio;

                    }

                }
                // if neighbor is active but the shared face is not active, we activate the shared face
                else
                {
                    if (!(new_active_1d2d_cves[i]->Face(n)->Active()))
                    {
                        new_active_1d2d_cves[i]->Face(n)->Active(true);

                        if (new_active_1d2d_cves[i]->Neighbor(n)->Dimension() == 0U)
                        {
                            activeInteriorPointFaces_.push_back(new_active_1d2d_cves[i]->Face(n));
                            monitoringPointFacesFlux_.erase(new_active_1d2d_cves[i]->Face(n));
                        }
                        else
                        {
                            activeInteriorLineFaces_.push_back(new_active_1d2d_cves[i]->Face(n));
                            monitoringLineFacesFlux_.erase(new_active_1d2d_cves[i]->Face(n));
                        }

                    }

                }

            }

        } // end looping over neighbors

    } // end looping over new active 1D-2D cves

    // looping over new active 0D cves to add new monitoring point faces
    for (size_t i = 0U; i < new_active_0d_cves.size(); i++)
    {
        for (size_t n = 0U; n < new_active_0d_cves[i]->Neighbors(); n++)
        {
            if (!(new_active_0d_cves[i]->Neighbor(n)->Active()))
            {
                const double oil_flux(this->CalculateOilFluxInteriorPointFace(*(new_active_0d_cves[i]->Face(n)), flow_model));
                const double water_flux(new_active_0d_cves[i]->Face(n)->PropertyValue(0U) - oil_flux);
                const double oil_flux_ratio((std::abs(oil_flux) + std::abs(water_flux)) != 0. ?
                                              oil_flux / (std::abs(oil_flux) + std::abs(water_flux)) : 0.);

                monitoringPointFacesFlux_[new_active_0d_cves[i]->Face(n)] = oil_flux_ratio;
            }
            // if neighbor is active but the shared face is not active, we activate the shared face
            else
            {
                if (!(new_active_0d_cves[i]->Face(n)->Active()))
                {
                    new_active_0d_cves[i]->Face(n)->Active(true);
                    activeInteriorPointFaces_.push_back(new_active_0d_cves[i]->Face(n));
                    monitoringPointFacesFlux_.erase(new_active_0d_cves[i]->Face(n));
                }

            }

        }

    } // end looping over new active 0D cves

    // activating the parent patches of new active 1d2d cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = new_active_1d2d_cves.begin();
        cveit != new_active_1d2d_cves.end(); ++cveit)
    {
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

    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = new_active_0d_cves.begin();
        cveit != new_active_0d_cves.end(); ++cveit)
    {
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

    // calculating divergence and minimum CFL based time step for new active region
    double min_cfl_timestep(std::numeric_limits<double>::max());

    // looping over 2D and 1D cves
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = new_active_1d2d_cves.begin();
        cveit != new_active_1d2d_cves.end(); cveit++)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); i++)
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
    for (typename std::vector<ControlVolumeElement<dim>*>::iterator cveit = new_active_0d_cves.begin();
        cveit != new_active_0d_cves.end(); cveit++)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for (size_t i = 0; i < (*cveit)->Faces(); i++)
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
double ActiveFaceTwoPhaseTransport<dim>::Transport(double total_time_increment, TwoPhaseModel<dim>& flow_model)
{
    std::cout << "\nActiveFaceTwoPhaseTransport<" << dim << ">::Transport()\n";

    typedef typename std::vector<ElementFace<dim>*>::iterator face_ptr_itr;
    typedef typename std::vector<ControlVolumeElement<dim>*>::iterator CVE_ptr_itr;

    static bool first_time(true);
    if (first_time)
    {
        this->ConstructFluxes();
        FindActiveRegions(flow_model, this->setup_.ActivationCriteria());
        ConstructFluxesActiveRegions();
        first_time = false;
    } 
    else
    {
        ConstructFluxesActiveRegions();
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
    std::cout << "Active Percentage: " << std:: fixed << (double) (active1D2D_CVEs_.size() + active0D_CVEs_.size()) /
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
            ZeroActiveCVEsProperty(0U);
            //transportModelRef_.ZeroCVEsProperty(0U);

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

            // transport through point faces
            // looping over interior active point faces
            for (face_ptr_itr fit = ActiveInteriorPointFacesBegin(); fit != ActiveInteriorPointFacesEnd(); fit++)
            {                
                const double oil_face_flux(this->CalculateOilFluxInteriorPointFace(*(*fit), flow_model));
                (*fit)->InsideCVE()->PropertyValueAdd(0U, oil_face_flux);
                (*fit)->OutsideCVE()->PropertyValueAdd(0U, -oil_face_flux);

            } // end for loop transport for inner line faces

            // transport over monitoring faces with constant initial values
            typedef typename std::map<ElementFace<dim>*, double>::iterator face_flux_ptr_itr;
            for (face_flux_ptr_itr fit = MonitoringLineFacesFluxBegin(); fit != MonitoringLineFacesFluxEnd(); fit++)
            {
                const double oil_flux(this->CalculateOilFluxInteriorLineFace(*(fit->first), flow_model));
                const double water_flux(fit->first->PropertyValue(0U) - oil_flux);
                const double total_flux(std::abs(oil_flux) + std::abs(water_flux));

                if (fit->first->InsideCVE()->Active())
                {
                    fit->first->InsideCVE()->PropertyValueAdd(0U, fit->second * total_flux);
                } 
                else
                {
                    fit->first->OutsideCVE()->PropertyValueAdd(0U, -fit->second * total_flux);
                }
            }

            for (face_flux_ptr_itr fit = MonitoringPointFacesFluxBegin(); fit != MonitoringPointFacesFluxEnd(); fit++)
            {
                const double oil_flux(this->CalculateOilFluxInteriorPointFace(*(fit->first), flow_model));
                const double water_flux(fit->first->PropertyValue(0U) - oil_flux);
                const double total_flux(std::abs(oil_flux) + std::abs(water_flux));

                if (fit->first->InsideCVE()->Active())
                {
                    fit->first->InsideCVE()->PropertyValueAdd(0U, fit->second * total_flux);
                } 
                else
                {
                    fit->first->OutsideCVE()->PropertyValueAdd(0U, -fit->second * total_flux);
                }
            }

            maximum_change_in_timestep = 0.;

            // looping over active inside 2D and 1D cves
            for (CVE_ptr_itr cveit = Active1D2D_CVEsBegin(); cveit != Active1D2D_CVEsEnd(); cveit++)
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
            for (CVE_ptr_itr cveit = Active0D_CVEsBegin(); cveit != Active0D_CVEsEnd(); cveit++)
            {
                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes((*cveit)->PropertyValue(0U));
                // assuming porosity is 1.
                changes *= (internal_time_increment / (*cveit)->Volume());

                (*cveit)->PropertyValue(0U, changes);

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
        for (CVE_ptr_itr cveit = Active1D2D_CVEsBegin(); cveit != Active1D2D_CVEsEnd(); cveit++)
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
        for (CVE_ptr_itr cveit = Active0D_CVEsBegin(); cveit != Active0D_CVEsEnd(); cveit++)
        {
            (*cveit)->PropertyValue(1U, (*cveit)->PropertyValue(1U) + (*cveit)->PropertyValue(0U));
        }

        const double new_cfl_timestep(UpdateActiveRegions(flow_model, this->setup_.Limiter()) * CFL_multiplier);
        if (CFL_timestep > new_cfl_timestep) CFL_timestep = new_cfl_timestep;

        total_time += previous_time_step;
        transport_step_number++;

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





template class ActiveFaceTwoPhaseTransport<2U>;


} // end namespace csmp
