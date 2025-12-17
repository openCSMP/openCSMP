#include "TwoPhaseElementBasedTransport.h"
#include "IMPES_Setup.h"
#include "Model.h"
#include "TransportModel.h"
#include "TwoPhaseModel.h"
#include "PropertyHandle.h"
#include "VTK_Interface.h"
#include "DenseMatrix.h"
#include "LU_Solver.h"
#include "SignalHandler.h"
#include "VariableOperations.h"


namespace csmp{

template<uint32_t dim>
TwoPhaseElementBasedTransport<dim>::TwoPhaseElementBasedTransport(Model<dim>& model,
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
  : modelRef_(model),
    regionRef_(model.Region(region_name)),
    transportModelRef_(transport_model),
    porosityKey_(model.Database().StorageKey(porosity)),
    permeabilityKey_(model.Database().StorageKey(permeability)),
    oilSaturationKey_(model.Database().StorageKey(saturation_oil)),
    waterSaturationKey_(model.Database().StorageKey(saturation_water)),
    divergenceKey_(model.Database().StorageKey(divergence)),
    totalVelocityKey_(model.Database().StorageKey(total_velocity)),
    totalProdRateKey_(model.Database().StorageKey(element_total_production_volume_rate)),
    oilInjRateKey_(model.Database().StorageKey(element_oil_injection_volume_rate)),
    waterInjRateKey_(model.Database().StorageKey(element_water_injection_volume_rate)),
    setup_(setup),
    withGravity_(setup.WithGravitationalForces()),
    withCapillary_(setup.WithCapillaryForces()),
    G_(setup.WithGravitationalForces() ? 9.8066 : 0.),
    EPSILON_(1e-18)

{
    if ( porosityKey_.place != ELEMENT || porosityKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'porosity' variable must be a scalar placed on the element" );

    if ( oilSaturationKey_.place != ELEMENT || oilSaturationKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'saturation oil' variable must be a scalar placed on the element" );

    if ( waterSaturationKey_.place != ELEMENT || waterSaturationKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'saturation water' variable must be a scalar placed on the element" );

    if ( permeabilityKey_.place != ELEMENT || permeabilityKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'permeability' variable must be a scalar placed on the element" );

    if ( divergenceKey_.place != ELEMENT || divergenceKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'divergence' variable must be a scalar placed on the element" );

    if ( totalVelocityKey_.place != ELEMENT || totalVelocityKey_.type != VECTOR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'total velocity' variable must be a vector placed on the element" );

    if ( totalProdRateKey_.place != ELEMENT || totalProdRateKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'total production volume rate' variable must be a scalar placed on the element" );

    if ( oilInjRateKey_.place != ELEMENT || oilInjRateKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'oil injection volume rate' variable must be a scalar placed on the element" );

    if ( waterInjRateKey_.place != ELEMENT || waterInjRateKey_.type != SCALAR )
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport(constructor)",
        "The 'water injection volume rate' variable must be a scalar placed on the element" );


    if (withCapillary_) CalculateBarycentersDistanceDerivativeTerm();

}



template<uint32_t dim>
TwoPhaseElementBasedTransport<dim>::~TwoPhaseElementBasedTransport()
{
}



template<uint32_t dim>
void TwoPhaseElementBasedTransport<dim>::AssignTransportBoundaryCondition(const char* boundary_node_oil_saturation)
{
    std::cout << "\nTwoPhaseElementBasedTransport<" << dim << "> Assigning boundary condition for transport.\n";

    csmp::Index oil_saturation_key((modelRef_.Database().StorageKey(boundary_node_oil_saturation)));

    if (oil_saturation_key.place != NODE || oil_saturation_key.type != SCALAR)
        throw csmp::Exception( FATAL_ERROR, "TwoPhaseElementBasedTransport::AssignTransportBoundaryCondition",
        "The 'boundary oil saturation' variable must be a scalar placed on the node" );

    // cleaning the previous boundary conditions. the second local storage is used for boundary condition
    for ( auto fit = transportModelRef_.PerimeterFacesBegin();
        fit != transportModelRef_.FacesEnd(); ++fit)
    {
        fit->PropertyValue(1U, 0.);
    }

    for ( auto pit = transportModelRef_.PerimeterPatchesBegin(); pit != transportModelRef_.PatchesEnd(); ++pit)
    {
        // we do this check to assign value to the faces from non-corner nodes
        if (((*pit).CenterNode()->AtBoundary() == BOTTOM) || ((*pit).CenterNode()->AtBoundary() == RIGHT) ||
            ((*pit).CenterNode()->AtBoundary() == TOP) || ((*pit).CenterNode()->AtBoundary() == LEFT))
        {
            for ( auto fit = pit->PerimeterFacesBegin(); fit != pit->FacesEnd(); ++fit)
            {
                (*fit)->PropertyValue(1U, (*pit).CenterNode()->Read(oil_saturation_key));
            }
        }

    } // end looping over boundary patches

} // end AssignTransportBoundaryCondition





// This method should be called only when capillary forces are exists
template<uint32_t dim>
void TwoPhaseElementBasedTransport<dim>::CalculateBarycentersDistanceDerivativeTerm()
{

    std::cout << "\nTwoPhaseElementBasedTransport<" << dim << ">::CalculateBarycentersDistanceDerivativeTerm()\n";

    for ( auto fit = transportModelRef_.FacesBegin(); fit != transportModelRef_.PerimeterFacesBegin(); ++fit)
      {
          const Point<dim> bc_in(fit->InsideCVE()->BaryCenter());
          const Point<dim> bc_out(fit->OutsideCVE()->BaryCenter());

          Point<dim> bc_vector(bc_out - bc_in);
          double normal_distance_derivative_term(dotProduct(bc_vector, fit->UnitNormal()));

          normal_distance_derivative_term /= (bc_in.DistanceTo(bc_out) * bc_in.DistanceTo(bc_out));

          fit->PropertyValue(1U, normal_distance_derivative_term);

      }

} // end CalculateBarycentersDistanceDerivativeTerm




template<uint32_t dim>
std::pair<double, double> TwoPhaseElementBasedTransport<dim>::CFL_Timestep_MaxDivergence()
{
    double min_cfl_timestep(std::numeric_limits<double>::max());
    double max_divergence(0.);

    // looping over interior 2D and 1D elements
    for (typename std::vector<ControlVolumeElement<dim> >::iterator eit = transportModelRef_.CVEsBegin();
        eit != transportModelRef_.PointCVEsBegin(); ++eit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for ( uint32_t i{0U}; i < eit->Faces(); ++i)
          {
              total_flux += std::abs(eit->Face(i)->PropertyValue(0U));
              flux_imbalance += eit->Face(i)->PropertyValue(0U) * eit->FaceNormalDirection(i);
          }

        flux_imbalance += (eit->E()->Read(totalProdRateKey_) + 
                           eit->E()->Read(oilInjRateKey_) + 
                           eit->E()->Read(waterInjRateKey_)) * 
                           eit->Volume();

        eit->E()->Store(divergenceKey_, ScalarVariable(PLAIN, flux_imbalance));

        const double element_timestep(eit->Volume() * eit->E()->Read(porosityKey_) * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);

    }

    // looping over boundary 2D and 1D elements
    for (typename std::vector<ControlVolumeElement<dim> >::iterator eit = transportModelRef_.PerimeterCVEsBegin();
        eit != transportModelRef_.PerimeterPointCVEsBegin(); ++eit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for ( uint32_t i{0U}; i < eit->Faces(); ++i)
          {
             total_flux += std::abs(eit->Face(i)->PropertyValue(0U));
             flux_imbalance += eit->Face(i)->PropertyValue(0U) * eit->FaceNormalDirection(i);
          }

        flux_imbalance += (eit->E()->Read(totalProdRateKey_) + 
                           eit->E()->Read(oilInjRateKey_) + 
                           eit->E()->Read(waterInjRateKey_)) * 
                           eit->Volume();

        eit->E()->Store(divergenceKey_, ScalarVariable(PLAIN, flux_imbalance));

        const double element_timestep(eit->Volume() * eit->E()->Read(porosityKey_) * 2. / total_flux);      

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);

    }

    // looping over 0D elements
    for (typename std::vector<ControlVolumeElement<dim> >::iterator eit = transportModelRef_.PointCVEsBegin();
        eit != transportModelRef_.PerimeterCVEsBegin(); ++eit)
    {
        double total_flux(0.);
        double flux_imbalance(0);

        for ( uint32_t i{0U}; i < eit->Faces(); ++i)
          {
             total_flux += std::abs(eit->Face(i)->PropertyValue(0U));
             flux_imbalance += eit->Face(i)->PropertyValue(0U) * eit->FaceNormalDirection(i);
          }

        const double element_timestep(eit->Volume() * 2. / total_flux);

        if (element_timestep < min_cfl_timestep) min_cfl_timestep = element_timestep;
        if (max_divergence < std::abs(flux_imbalance)) max_divergence = std::abs(flux_imbalance);
    }

    return std::make_pair(min_cfl_timestep, max_divergence);

} // end CFL_Timestep


/*

// This is the modified method which uses the square root of permeability instead of permeability itself
// to increase the flow inside low permeability zones

template<uint32_t dim>
void TwoPhaseElementBasedTransport<dim>::ConstructFluxes()
{

    std::cout << "\nTwoPhaseElementBasedTransport<" << dim << ">::ConstructFluxes()\n";

    VectorVariable<dim> vt;
    VectorVariable<dim> vt_inside;
    VectorVariable<dim> vt_outside;
    double perm_inside;
    double perm_outside;
    LU_Solver solver;
    DenseMatrix<DM_MIN> A;
    std::vector<double> b;
    std::vector<double> fluxes;
    std::vector<double> sector_face_area; // storing the sector face area for later uses

    typedef typename std::vector<CVE_Patch<dim> >::iterator patch_itr;

    // make previous fluxes zero
    transportModelRef_.ZeroFacesProperty(0U);

    // in this method (ConstructFluxes) it is assumed that the influx to the sector/element is positive
    // and outflux from the sector/element is negative
    // construct fluxes for inside patches
    for (patch_itr pit = transportModelRef_.PatchesBegin(); pit != transportModelRef_.PerimeterPatchesBegin(); ++pit)
    {
        b.resize((*pit).Faces());
        fluxes.resize((*pit).Faces());
        sector_face_area.resize((*pit).Faces());
        A.Resize((*pit).Faces(), (*pit).Faces());
        A.Zero();

        double min_perm(std::numeric_limits<double>::max());

        // we construct one less equation since it is a cyclic problem
        for (size_t s = 0U; s < (*pit).Sectors() - 1U; ++s)
        {            
            if ((*pit).SectorAccessor(s)->E()->Read(permeabilityKey_) < min_perm)
                min_perm = (*pit).SectorAccessor(s)->E()->Read(permeabilityKey_);

            if ((*pit).SectorAccessor(s)->CVE()->Dimension() != 0U)
            {
                (*pit).SectorAccessor(s)->E()->Read(totalVelocityKey_, vt);

                for (size_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);

                }

                b[s] = vt.DotProduct((*pit).SectorAccessor(s)->TotalFacetsAreaVector());

                b[s] -= ((*pit).SectorAccessor(s)->E()->Read(oilInjRateKey_) + 
                    (*pit).SectorAccessor(s)->E()->Read(waterInjRateKey_) +
                    (*pit).SectorAccessor(s)->E()->Read(totalProdRateKey_)) *
                    (*pit).SectorAccessor(s)->CVE()->Volume() / 
                    (*pit).SectorAccessor(s)->E()->Nodes();

            }
            else
            {
                for (size_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);

                    b[s] = 0.;

                }

            }

        }

        if ((*pit).SectorAccessor((*pit).Sectors() - 1U)->CVE()->Dimension() != 0)
        {
            if ((*pit).SectorAccessor((*pit).Sectors() - 1U)->E()->Read(permeabilityKey_) < min_perm)
                min_perm = (*pit).SectorAccessor((*pit).Sectors() - 1U)->E()->Read(permeabilityKey_);
        }

        // constructing auxiliary equations to close the system
        double total_flux(0.);
        size_t equation_number((*pit).Sectors() - 1U);
        std::vector<bool> not_visited_face((*pit).Faces(), true);

        for (size_t s = 0U; s < (*pit).Sectors(); ++s)
        {
            // each sub-patch (a sectors collection) which is isolated by two line elements gives an
            // auxiliary equation
            if ((*pit).SectorAccessor(s)->CVE()->Dimension() == 2U)
            {
                for (size_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    if (not_visited_face[(*pit).SectorAccessor(s)->PatchFaceNumber(f)])
                    {
                        not_visited_face[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = false;
                        perm_inside = std::sqrt((*pit).SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(permeabilityKey_)); // ------------------->
                        perm_outside = std::sqrt((*pit).SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(permeabilityKey_)); // ------------------>
                        (*pit).SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(totalVelocityKey_, vt_inside);
                        (*pit).SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(totalVelocityKey_, vt_outside);

                        total_flux -= vt_inside.DotProduct((*pit).SectorAccessor(s)->FaceNormal(f)) * 
                            (*pit).SectorAccessor(s)->FaceArea(f) / perm_inside
                            * (*pit).FaceCyclicDirection((*pit).SectorAccessor(s)->PatchFaceNumber(f));
                        total_flux -= vt_outside.DotProduct((*pit).SectorAccessor(s)->FaceNormal(f)) * 
                            (*pit).SectorAccessor(s)->FaceArea(f) / perm_outside
                            * (*pit).FaceCyclicDirection((*pit).SectorAccessor(s)->PatchFaceNumber(f));

                        A(equation_number, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                            (1. / perm_inside + 1. / perm_outside) * min_perm * (*pit).SectorAccessor(s)->FaceArea(f)
                            * (*pit).FaceCyclicDirection((*pit).SectorAccessor(s)->PatchFaceNumber(f));
                    }

                }

                if (s == (*pit).Sectors() - 1U)
                {
                    b[equation_number] = total_flux * min_perm;
                }

            }
            else
            {
                if (equation_number < (*pit).Faces()) b[equation_number] = total_flux * min_perm;
                total_flux = 0.;
                ++equation_number;
                continue;
            }

        }

        solver.Solve(A, b, fluxes, (*pit).Faces());

        // assigning the results to the corresponding faces
        typename std::vector<double>::iterator rit(fluxes.begin());
        typename std::vector<double>::iterator ait(sector_face_area.begin());
        for (typename std::vector<ElementFace<dim>*>::iterator fit = (*pit).FacesBegin(); fit != (*pit).FacesEnd(); ++fit, ++rit, ++ait)
        {
            (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
        }

    } // end interior patches for loop


    // construct fluxes for boundary patches
    for (patch_itr pit = transportModelRef_.PerimeterPatchesBegin(); pit != transportModelRef_.PatchesEnd(); ++pit)
    {
        b.resize((*pit).Faces());
        fluxes.resize((*pit).Faces());
        sector_face_area.resize((*pit).Faces());
        A.Resize((*pit).Faces(), (*pit).Faces());
        A.Zero();

        double min_perm(std::numeric_limits<double>::max());

        for (size_t s = 0U; s < (*pit).Sectors(); ++s)
        {            
            if ((*pit).SectorAccessor(s)->CVE()->Dimension() != 0U)
            {
                if ((*pit).SectorAccessor(s)->E()->Read(permeabilityKey_) < min_perm)
                    min_perm = (*pit).SectorAccessor(s)->E()->Read(permeabilityKey_);

                (*pit).SectorAccessor(s)->E()->Read(totalVelocityKey_, vt);

                for (size_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) =  
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);
                }

                b[s] = vt.DotProduct((*pit).SectorAccessor(s)->TotalFacetsAreaVector());

                b[s] -= ((*pit).SectorAccessor(s)->E()->Read(oilInjRateKey_) + 
                    (*pit).SectorAccessor(s)->E()->Read(waterInjRateKey_) +
                    (*pit).SectorAccessor(s)->E()->Read(totalProdRateKey_)) *
                    (*pit).SectorAccessor(s)->CVE()->Volume() / 
                    (*pit).SectorAccessor(s)->E()->Nodes();
            }
            else
            {
                for (size_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);

                    b[s] = 0.;

                }

            }

        }

        // constructing auxiliary equations to close the system
        size_t equation_number((*pit).Sectors());
        double perm_boundary_face_1(std::sqrt((*pit).FaceAccessor((*pit).InteriorFaces())->InsideCVE()->E()->Read(permeabilityKey_))); // ------------>
        double area_boundary_face_1((*pit).FaceAccessor((*pit).InteriorFaces())->Area() / 2.);

        for (size_t f = (*pit).InteriorFaces() + 1U; f < (*pit).Faces(); ++f)
        {
            if ((*pit).FaceAccessor(f)->InsideCVE()->Dimension() != 0U)
            {
                double perm_boundary_face_2(std::sqrt((*pit).FaceAccessor(f)->InsideCVE()->E()->Read(permeabilityKey_))); // -------------------->
                double area_boundary_face_2((*pit).FaceAccessor(f)->InsideCVE()->Dimension() == 2U ?
                    (*pit).FaceAccessor(f)->Area() / 2. : (*pit).FaceAccessor(f)->Area());

                b[equation_number] = 0.;
                A(equation_number, (*pit).InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
                A(equation_number, f) = -area_boundary_face_2 / perm_boundary_face_2;

                ++equation_number;
            }
            else
            {
                for (size_t n = 0U; n < (*pit).FaceAccessor(f)->InsideCVE()->Faces(); ++n)
                {
                    if ((*pit).FaceAccessor(f)->InsideCVE()->Face(n)->Placement() == NOT)
                    {
                        double perm_boundary_face_2(std::sqrt((*pit).FaceAccessor(f)->InsideCVE()->Neighbor(n)->E()->Read(permeabilityKey_))); // ------------->
                        double area_boundary_face_2((*pit).FaceAccessor(f)->InsideCVE()->Face(n)->Area());

                        // finding the face number in the patch
                        size_t face_num;
                        for (face_num = 0U; face_num < (*pit).InteriorFaces(); ++face_num)
                        {
                            if ((*pit).FaceAccessor(face_num) == (*pit).FaceAccessor(f)->InsideCVE()->Face(n)) break;
                        }

                        b[equation_number] = 0.;
                        A(equation_number, (*pit).InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
                        A(equation_number, face_num) = area_boundary_face_2 / perm_boundary_face_2 * 
                            (*pit).FaceAccessor(f)->InsideCVE()->FaceNormalDirection(n);

                        ++equation_number;
                    }

                }

            }

        }

        solver.Solve(A, b, fluxes, (*pit).Faces());

        // assigning the results to the corresponding faces
        typename std::vector<double>::iterator rit(fluxes.begin());
        typename std::vector<double>::iterator ait(sector_face_area.begin());
        for (typename std::vector<ElementFace<dim>*>::iterator fit = (*pit).FacesBegin(); fit != (*pit).FacesEnd(); ++fit, ++rit, ++ait)
        {
            (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
        }

    } // end boundary patches for loop

} // end ConstructFluxes

*/


template<uint32_t dim>
void TwoPhaseElementBasedTransport<dim>::ConstructFluxes()
{

    std::cout << "\nTwoPhaseElementBasedTransport<" << dim << ">::ConstructFluxes()\n";

    VectorVariable<dim> vt;
    VectorVariable<dim> vt_inside;
    VectorVariable<dim> vt_outside;
    double perm_inside;
    double perm_outside;
    LU_Solver solver;
    DenseMatrix<DM_MIN> A;
    std::vector<double> b;
    std::vector<double> fluxes;
    std::vector<double> sector_face_area; // storing the sector face area for later uses

    typedef typename std::vector<CVE_Patch<dim> >::iterator patch_itr;

    // make previous fluxes zero
    transportModelRef_.ZeroFacesProperty(0U);

    // in this method (ConstructFluxes) it is assumed that the influx to the sector/element is positive
    // and outflux from the sector/element is negative
    // construct fluxes for inside patches
    for (patch_itr pit = transportModelRef_.PatchesBegin(); pit != transportModelRef_.PerimeterPatchesBegin(); ++pit)
    {
        b.resize((*pit).Faces());
        fluxes.resize((*pit).Faces());
        sector_face_area.resize((*pit).Faces());
        A.Resize((*pit).Faces(), (*pit).Faces());
        A.Zero();

        double min_perm(std::numeric_limits<double>::max());

        // we construct one less equation since it is a cyclic problem
        for (uint32_t s = 0U; s < (*pit).Sectors() - 1U; ++s)
        {            
            if ((*pit).SectorAccessor(s)->E()->Read(permeabilityKey_) < min_perm)
                min_perm = (*pit).SectorAccessor(s)->E()->Read(permeabilityKey_);

            if ((*pit).SectorAccessor(s)->CVE()->Dimension() != 0U)
            {
                (*pit).SectorAccessor(s)->E()->Read(totalVelocityKey_, vt);

                for (uint32_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);

                }

                b[s] = dotProduct(vt,(*pit).SectorAccessor(s)->TotalFacetsAreaVector());

                b[s] -= ((*pit).SectorAccessor(s)->E()->Read(oilInjRateKey_) + 
                    (*pit).SectorAccessor(s)->E()->Read(waterInjRateKey_) +
                    (*pit).SectorAccessor(s)->E()->Read(totalProdRateKey_)) *
                    (*pit).SectorAccessor(s)->CVE()->Volume() / 
                    (*pit).SectorAccessor(s)->E()->Nodes();

            }
            else
            {
                for (uint32_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);

                    b[s] = 0.;

                }

            }

        }

        if ((*pit).SectorAccessor((*pit).Sectors() - 1U)->CVE()->Dimension() != 0)
        {
            if ((*pit).SectorAccessor((*pit).Sectors() - 1U)->E()->Read(permeabilityKey_) < min_perm)
                min_perm = (*pit).SectorAccessor((*pit).Sectors() - 1U)->E()->Read(permeabilityKey_);
        }

        // constructing auxiliary equations to close the system
        double total_flux(0.);
        uint32_t equation_number((*pit).Sectors() - 1U);
        std::vector<bool> not_visited_face((*pit).Faces(), true);

        for (uint32_t s = 0U; s < (*pit).Sectors(); ++s)
        {
            // each sub-patch (a sectors collection) which is isolated by two line elements gives an
            // auxiliary equation
            if ((*pit).SectorAccessor(s)->CVE()->Dimension() == 2U)
            {
                for (uint32_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    if (not_visited_face[(*pit).SectorAccessor(s)->PatchFaceNumber(f)])
                    {
                        not_visited_face[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = false;
                        perm_inside = (*pit).SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(permeabilityKey_);
                        perm_outside = (*pit).SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(permeabilityKey_);
                        (*pit).SectorAccessor(s)->Face(f)->InsideCVE()->E()->Read(totalVelocityKey_, vt_inside);
                        (*pit).SectorAccessor(s)->Face(f)->OutsideCVE()->E()->Read(totalVelocityKey_, vt_outside);

                        total_flux -= dotProduct(vt_inside,(*pit).SectorAccessor(s)->FaceNormal(f)) *
                            (*pit).SectorAccessor(s)->FaceArea(f) / perm_inside
                            * (*pit).FaceCyclicDirection((*pit).SectorAccessor(s)->PatchFaceNumber(f));
                        total_flux -= dotProduct(vt_outside,(*pit).SectorAccessor(s)->FaceNormal(f)) *
                            (*pit).SectorAccessor(s)->FaceArea(f) / perm_outside
                            * (*pit).FaceCyclicDirection((*pit).SectorAccessor(s)->PatchFaceNumber(f));

                        A(equation_number, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                            (1. / perm_inside + 1. / perm_outside) * min_perm * (*pit).SectorAccessor(s)->FaceArea(f)
                            * (*pit).FaceCyclicDirection((*pit).SectorAccessor(s)->PatchFaceNumber(f));
                    }

                }

                if (s == (*pit).Sectors() - 1U)
                {
                    b[equation_number] = total_flux * min_perm;
                }

            }
            else
            {
                if (equation_number < (*pit).Faces()) b[equation_number] = total_flux * min_perm;
                total_flux = 0.;
                ++equation_number;
                continue;
            }

        }

        solver.Solve(A, b, fluxes, (*pit).Faces());

        // assigning the results to the corresponding faces
        typename std::vector<double>::iterator rit(fluxes.begin());
        typename std::vector<double>::iterator ait(sector_face_area.begin());
        for (typename std::vector<ElementFace<dim>*>::iterator fit = (*pit).FacesBegin(); fit != (*pit).FacesEnd(); ++fit, ++rit, ++ait)
        {
            (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
        }

    } // end interior patches for loop


    // construct fluxes for boundary patches
    for (patch_itr pit = transportModelRef_.PerimeterPatchesBegin(); pit != transportModelRef_.PatchesEnd(); ++pit)
    {
        b.resize((*pit).Faces());
        fluxes.resize((*pit).Faces());
        sector_face_area.resize((*pit).Faces());
        A.Resize((*pit).Faces(), (*pit).Faces());
        A.Zero();

        double min_perm(std::numeric_limits<double>::max());

        for (uint32_t s = 0U; s < (*pit).Sectors(); ++s)
        {            
            if ((*pit).SectorAccessor(s)->CVE()->Dimension() != 0U)
            {
                if ((*pit).SectorAccessor(s)->E()->Read(permeabilityKey_) < min_perm)
                    min_perm = (*pit).SectorAccessor(s)->E()->Read(permeabilityKey_);

                (*pit).SectorAccessor(s)->E()->Read(totalVelocityKey_, vt);

                for (uint32_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) =  
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);
                }

                b[s] = dotProduct(vt,(*pit).SectorAccessor(s)->TotalFacetsAreaVector());

                b[s] -= ((*pit).SectorAccessor(s)->E()->Read(oilInjRateKey_) + 
                    (*pit).SectorAccessor(s)->E()->Read(waterInjRateKey_) +
                    (*pit).SectorAccessor(s)->E()->Read(totalProdRateKey_)) *
                    (*pit).SectorAccessor(s)->CVE()->Volume() / 
                    (*pit).SectorAccessor(s)->E()->Nodes();
            }
            else
            {
                for (uint32_t f = 0U; f < (*pit).SectorAccessor(s)->Faces(); ++f)
                {
                    // storing the sector face area for later uses
                    sector_face_area[(*pit).SectorAccessor(s)->PatchFaceNumber(f)] = (*pit).SectorAccessor(s)->FaceArea(f);

                    A(s, (*pit).SectorAccessor(s)->PatchFaceNumber(f)) = 
                        (*pit).SectorAccessor(s)->FaceNormalDirection(f) * (*pit).SectorAccessor(s)->FaceArea(f);

                    b[s] = 0.;

                }

            }

        }

        // constructing auxiliary equations to close the system
        uint32_t equation_number((*pit).Sectors());
        double perm_boundary_face_1((*pit).FaceAccessor((*pit).InteriorFaces())->InsideCVE()->E()->Read(permeabilityKey_));
        double area_boundary_face_1((*pit).FaceAccessor((*pit).InteriorFaces())->Area() / 2.);

        for (uint32_t f = (*pit).InteriorFaces() + 1U; f < (*pit).Faces(); ++f)
        {
            if ((*pit).FaceAccessor(f)->InsideCVE()->Dimension() != 0U)
            {
                double perm_boundary_face_2((*pit).FaceAccessor(f)->InsideCVE()->E()->Read(permeabilityKey_));
                double area_boundary_face_2((*pit).FaceAccessor(f)->InsideCVE()->Dimension() == 2U ?
                    (*pit).FaceAccessor(f)->Area() / 2. : (*pit).FaceAccessor(f)->Area());

                b[equation_number] = 0.;
                A(equation_number, (*pit).InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
                A(equation_number, f) = -area_boundary_face_2 / perm_boundary_face_2;

                ++equation_number;
            }
            else
            {
                for (uint32_t n = 0U; n < (*pit).FaceAccessor(f)->InsideCVE()->Faces(); ++n)
                {
                    if ((*pit).FaceAccessor(f)->InsideCVE()->Face(n)->Placement() == NOT)
                    {
                        double perm_boundary_face_2((*pit).FaceAccessor(f)->InsideCVE()->Neighbor(n)->E()->Read(permeabilityKey_));
                        double area_boundary_face_2((*pit).FaceAccessor(f)->InsideCVE()->Face(n)->Area());

                        // finding the face number in the patch
                        uint32_t face_num;
                        for (face_num = 0U; face_num < (*pit).InteriorFaces(); ++face_num)
                        {
                            if ((*pit).FaceAccessor(face_num) == (*pit).FaceAccessor(f)->InsideCVE()->Face(n)) break;
                        }

                        b[equation_number] = 0.;
                        A(equation_number, (*pit).InteriorFaces()) = area_boundary_face_1 / perm_boundary_face_1;
                        A(equation_number, face_num) = area_boundary_face_2 / perm_boundary_face_2 * 
                            (*pit).FaceAccessor(f)->InsideCVE()->FaceNormalDirection(n);

                        ++equation_number;
                    }

                }

            }

        }

        solver.Solve(A, b, fluxes, (*pit).Faces());

        // assigning the results to the corresponding faces
        typename std::vector<double>::iterator rit(fluxes.begin());
        typename std::vector<double>::iterator ait(sector_face_area.begin());
        for (typename std::vector<ElementFace<dim>*>::iterator fit = (*pit).FacesBegin(); fit != (*pit).FacesEnd(); ++fit, ++rit, ++ait)
        {
            (*fit)->PropertyValueAdd(0U, (*rit) * (*ait));
        }

    } // end boundary patches for loop

} // end ConstructFluxes



template<uint32_t dim>
void TwoPhaseElementBasedTransport<dim>::FindUpstreamProperties(const double& face_flux, const double& area, 
                                                                const Point<dim>& face_normal, const double& perm,
                                                                const double& lambda_w_out, const double& lambda_o_out,
                                                                const double& lambda_w_in, const double& lambda_o_in,
                                                                const double& ro_w_out, const double& ro_o_out,
                                                                const double& ro_w_in, const double& ro_o_in,
                                                                const double& normal_capillary_grad,
                                                                double& lambda_w_upstream, double& lambda_o_upstream,
                                                                double& ro_w_upstream, double& ro_o_upstream) const
{
    if (!withGravity_ && !withCapillary_)
    {
        if (face_flux >= 0.)
        {
            lambda_w_upstream = lambda_w_out;
            lambda_o_upstream = lambda_o_out;
            ro_w_upstream = ro_w_out;
            ro_o_upstream = ro_o_out;
            return;
        }
        else
        {
            lambda_w_upstream = lambda_w_in;
            lambda_o_upstream = lambda_o_in;
            ro_w_upstream = ro_w_in;
            ro_o_upstream = ro_o_in;
            return;
        }
    }
    else
    {
        const double face_normal_y(face_normal[1U]);

        if (face_flux >= 0.)
        {
            if (face_normal_y >= 0.)
            {
                if (ro_w_out * face_normal_y * G_ - normal_capillary_grad > ro_o_out * face_normal_y * G_)
                {
                    // lambda_w = lambda_w_left
                    double dp;
                    const double lambda_o(0.5 * (lambda_o_out + lambda_o_in));
                    const double ro_o(0.5 * (ro_o_out + ro_o_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w_out * ro_w_out + lambda_o * ro_o) * G_ * face_normal_y;
                    dp += lambda_w_out * normal_capillary_grad;
                    dp /= (lambda_w_out + lambda_o + EPSILON_);

                    lambda_w_upstream = lambda_w_out;
                    lambda_o_upstream = (dp + ro_o * G_ * face_normal_y >= 0.? lambda_o_out : lambda_o_in);
                    ro_w_upstream = ro_w_out;
                    ro_o_upstream = (dp + ro_o * G_ * face_normal_y >= 0.? ro_o_out : ro_o_in);

                    return;

                } // end if (ro_w_left * face_normal_y * g - normal_capillary_grad > ro_o_left * face_normal_y * g)
                else
                {
                    // lambda_o = lambda_o_left
                    double dp;
                    const double lambda_w(0.5 * (lambda_w_out + lambda_w_in));
                    const double ro_w(0.5 * (ro_w_out + ro_w_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w * ro_w + lambda_o_out * ro_o_out) * G_ * face_normal_y;
                    dp += lambda_w * normal_capillary_grad;
                    dp /= (lambda_w + lambda_o_out + EPSILON_);

                    lambda_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad >= 0.? lambda_w_out : lambda_w_in);
                    lambda_o_upstream = lambda_o_out;
                    ro_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad >= 0.? ro_w_out : ro_w_in);
                    ro_o_upstream = ro_o_out;

                    return;

                } // end else (ro_w_left * face_normal_y * g - normal_capillary_grad > ro_o_left * face_normal_y * g)

            } // end if (*pit).SectorAccessor(s)->LeftFaceNormal()[1U] >= 0.
            else
            {
                if (ro_w_out * face_normal_y * G_ - normal_capillary_grad < ro_o_out * face_normal_y * G_)
                {
                    // lambda_o = lambda_o_left
                    double dp;
                    const double lambda_w(0.5 * (lambda_w_out + lambda_w_in));
                    const double ro_w(0.5 * (ro_w_out + ro_w_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w * ro_w + lambda_o_out * ro_o_out) * G_ * face_normal_y;
                    dp += lambda_w * normal_capillary_grad;
                    dp /= (lambda_w + lambda_o_out + EPSILON_);

                    lambda_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad >= 0.? lambda_w_out : lambda_w_in);
                    lambda_o_upstream = lambda_o_out;
                    ro_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad >= 0.? ro_w_out : ro_w_in);
                    ro_o_upstream = ro_o_out;

                    return;

                } // end if (ro_w_left * face_normal_y * g - normal_capillary_grad < ro_o_left * face_normal_y * g)
                else
                {
                    // lambda_w = lambda_w_left
                    double dp;
                    const double lambda_o(0.5 * (lambda_o_out + lambda_o_in));
                    const double ro_o(0.5 * (ro_o_out + ro_o_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w_out * ro_w_out + lambda_o * ro_o) * G_ * face_normal_y;
                    dp += lambda_w_out * normal_capillary_grad;
                    dp /= (lambda_w_out + lambda_o + EPSILON_);

                    lambda_w_upstream = lambda_w_out;
                    lambda_o_upstream = (dp + ro_o * G_ * face_normal_y >= 0.? lambda_o_out : lambda_o_in);
                    ro_w_upstream = ro_w_out;
                    ro_o_upstream = (dp + ro_o * G_ * face_normal_y >= 0.? ro_o_out : ro_o_in);

                    return;

                } // end else (ro_w_left * face_normal_y * g - normal_capillary_grad < ro_o_left * face_normal_y * g)

            } // end else (*pit).SectorAccessor(s)->LeftFaceNormal()[1U] >= 0. 

        } // end if face_flux_[counter][s] >= 0.
        else
        {
            if (face_normal_y >= 0.)
            {
                if (ro_w_out * face_normal_y * G_ - normal_capillary_grad > ro_o_out * face_normal_y * G_)
                {
                    // lambda_o = lambda_o_right
                    double dp;
                    const double lambda_w(0.5 * (lambda_w_out + lambda_w_in));
                    const double ro_w(0.5 * (ro_w_out + ro_w_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w * ro_w + lambda_o_in * ro_o_in) * G_ * face_normal_y;
                    dp += lambda_w * normal_capillary_grad;
                    dp /= (lambda_w + lambda_o_in + EPSILON_);

                    lambda_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad <= 0.? lambda_w_in : lambda_w_out);
                    lambda_o_upstream = lambda_o_in;
                    ro_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad <= 0.? ro_w_in : ro_w_out);
                    ro_o_upstream = ro_o_in;

                    return;

                } // end if (ro_w_left * face_normal_y * g - normal_capillary_grad > ro_o_left * face_normal_y * g)
                else
                {
                    // lambda_w = lambda_w_right
                    double dp;
                    const double lambda_o(0.5 * (lambda_o_out + lambda_o_in));
                    const double ro_o(0.5 * (ro_o_out + ro_o_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w_in * ro_w_in + lambda_o * ro_o) * G_ * face_normal_y;
                    dp += lambda_w_in * normal_capillary_grad;
                    dp /= (lambda_w_in + lambda_o + EPSILON_);

                    lambda_w_upstream = lambda_w_in;
                    lambda_o_upstream = (dp + ro_o * G_ * face_normal_y <= 0.? lambda_o_in : lambda_o_out);
                    ro_w_upstream = ro_w_in;
                    ro_o_upstream = (dp + ro_o * G_ * face_normal_y <= 0.? ro_o_in : ro_o_out);

                    return;

                } // end else (ro_w_left * face_normal_y * g - normal_capillary_grad > ro_o_left * face_normal_y * g)

            } // end if (*pit).SectorAccessor(s)->LeftFaceNormal()[1U] >= 0.
            else
            {
                if (ro_w_out * face_normal_y * G_ - normal_capillary_grad < ro_o_out * face_normal_y * G_)
                {
                    // lambda_w = lambda_w_right
                    double dp;
                    const double lambda_o(0.5 * (lambda_o_out + lambda_o_in));
                    const double ro_o(0.5 * (ro_o_out + ro_o_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w_in * ro_w_in + lambda_o * ro_o) * G_ * face_normal_y;
                    dp += lambda_w_in * normal_capillary_grad;
                    dp /= (lambda_w_in + lambda_o + EPSILON_);

                    lambda_w_upstream = lambda_w_in;
                    lambda_o_upstream = (dp + ro_o * G_ * face_normal_y <= 0.? lambda_o_in : lambda_o_out);
                    ro_w_upstream = ro_w_in;
                    ro_o_upstream = (dp + ro_o * G_ * face_normal_y <= 0.? ro_o_in : ro_o_out);

                    return;

                } // end if (ro_w_left * face_normal_y * g - normal_capillary_grad < ro_o_left * face_normal_y * g)
                else
                {
                    // lambda_o = lambda_o_right
                    double dp;
                    const double lambda_w(0.5 * (lambda_w_out + lambda_w_in));
                    const double ro_w(0.5 * (ro_w_out + ro_w_in));

                    dp = face_flux / area / perm;

                    dp -= (lambda_w * ro_w + lambda_o_in * ro_o_in) * G_ * face_normal_y;
                    dp += lambda_w * normal_capillary_grad;
                    dp /= (lambda_w + lambda_o_in + EPSILON_);

                    lambda_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad <= 0.? lambda_w_in : lambda_w_out);
                    lambda_o_upstream = lambda_o_in;
                    ro_w_upstream = (dp + ro_w * G_ * face_normal_y - normal_capillary_grad <= 0.? ro_w_in : ro_w_out);
                    ro_o_upstream = ro_o_in;

                    return;

                } // end else (ro_w_left > ro_o_left)

            } // end else (*pit).SectorAccessor(s)->LeftFaceNormal()[1U] >= 0. 

        } // end else face_flux_[counter][s] >= 0.

    }

} // end FindUpstreamProperties



template<uint32_t dim>
double TwoPhaseElementBasedTransport<dim>::CalculateOilFluxInteriorLineFace(ElementFace<dim>& face, TwoPhaseModel<dim>& flow_model)
{

    flow_model.Initialize(*(face.InsideCVE()->E()));
    flow_model.EffectiveSaturation();
    const double perm_in(flow_model.Permeability());
    const double ro_w_in(flow_model.DensityWettingPhase());
    const double ro_o_in(flow_model.DensityNonWettingPhase());
    const double lambda_w_in(flow_model.MobilityPhase(1U));
    const double lambda_o_in(flow_model.MobilityPhase(2U));
    const double pc_in(flow_model.pc_Phase());

    flow_model.Initialize(*(face.OutsideCVE()->E()));
    flow_model.EffectiveSaturation();
    const double perm_out(flow_model.Permeability());
    const double ro_w_out(flow_model.DensityWettingPhase());
    const double ro_o_out(flow_model.DensityNonWettingPhase());
    const double lambda_w_out(flow_model.MobilityPhase(1U));
    const double lambda_o_out(flow_model.MobilityPhase(2U));
    const double pc_out(flow_model.pc_Phase());

    const double capillary_grad(withCapillary_ ? face.PropertyValue(1U) * (pc_out - pc_in) : 0.);
    const double perm(2. / (1. / perm_out + 1. / perm_in));

    double lambda_w, lambda_o;
    double ro_w, ro_o;

    FindUpstreamProperties(face.PropertyValue(0U), face.Area(), face.UnitNormal(),
                           perm, lambda_w_out, lambda_o_out, lambda_w_in, lambda_o_in, 
                           ro_w_out, ro_o_out, ro_w_in, ro_o_in, capillary_grad,
                           lambda_w, lambda_o, ro_w, ro_o);

    // face vertical normal area
    const double fvna(face.UnitNormal()[1U] * face.Area());

    double oil_face_flux(face.PropertyValue(0U) / perm);

    oil_face_flux += lambda_w * capillary_grad * face.Area();
    oil_face_flux -= (lambda_w * ro_w + lambda_o * ro_o) * G_ * fvna;
    oil_face_flux /= (lambda_w + lambda_o + EPSILON_);
    oil_face_flux += ro_o * G_ * fvna;
    oil_face_flux *= lambda_o * perm;

    return oil_face_flux;

} // end CalculateOilFluxInteriorLineFace




template<uint32_t dim>
double TwoPhaseElementBasedTransport<dim>::CalculateOilFluxBoundaryLineFace(ElementFace<dim>& face, TwoPhaseModel<dim>& flow_model)
{
/*
    if (face.PropertyValue(0U) > 0.)
    {
        return face.PropertyValue(0U) * face.PropertyValue(1U);
    }
    else
    {
        flow_model.Initialize(*(face.InsideCVE()->E()));
        flow_model.EffectiveSaturation();
        const double lambda_w(flow_model.MobilityPhase(1U));
        const double lambda_o(flow_model.MobilityPhase(2U));

        return face.PropertyValue(0U) * lambda_o / (lambda_w + lambda_o + EPSILON_);
    }
*/

    if (face.PropertyValue(0U) > 0.)
    {
        return face.PropertyValue(0U) * face.PropertyValue(1U);
    }
    else
    {
/*
        flow_model.Initialize(*(face.InsideCVE()->E()));
        flow_model.EffectiveSaturation();
        const double lambda_w_in(flow_model.MobilityPhase(1U));
        const double lambda_o_in(flow_model.MobilityPhase(2U));

        return face.PropertyValue(0U) * lambda_o_in / (lambda_w_in + lambda_o_in + EPSILON_);
*/
        
        flow_model.Initialize(*(face.InsideCVE()->E()));
        flow_model.EffectiveSaturation();
        const double perm_in(flow_model.Permeability());
        const double ro_w_in(flow_model.DensityWettingPhase());
        const double ro_o_in(flow_model.DensityNonWettingPhase());
        const double lambda_w_in(flow_model.MobilityPhase(1U));
        const double lambda_o_in(flow_model.MobilityPhase(2U));
//        const double pc_in(flow_model.pc_Phase());

        double lambda_w, lambda_o;
        double ro_w, ro_o;

        FindUpstreamProperties(face.PropertyValue(0U), face.Area(), face.UnitNormal(),
                               perm_in, lambda_w_in, lambda_o_in, lambda_w_in, lambda_o_in, 
                               ro_w_in, ro_o_in, ro_w_in, ro_o_in, 0,
                               lambda_w, lambda_o, ro_w, ro_o);

        // face vertical normal area
        const double fvna(face.UnitNormal()[1U] * face.Area());

        double oil_face_flux(face.PropertyValue(0U) / perm_in);

        oil_face_flux += lambda_w * 0. * face.Area();
        oil_face_flux -= (lambda_w * ro_w + lambda_o * ro_o) * G_ * fvna;
        oil_face_flux /= (lambda_w + lambda_o + EPSILON_);
        oil_face_flux += ro_o * G_ * fvna;
        oil_face_flux *= lambda_o * perm_in;

        return oil_face_flux;
    
    }
    
} // end CalculateOilFluxBoundaryLineFace




template<uint32_t dim>
double TwoPhaseElementBasedTransport<dim>::CalculateOilFluxInteriorPointFace(ElementFace<dim>& face, TwoPhaseModel<dim>& flow_model)
{
    double perm_in(0.);
    double ro_w_in(0.);
    double ro_o_in(0.);
    double lambda_w_in(0.);
    double lambda_o_in(0.);
    double pc_in(0.);

    double perm_out(0.);
    double ro_w_out(0.);
    double ro_o_out(0.);
    double lambda_w_out(0.);
    double lambda_o_out(0.);
    double pc_out(0.);

    if (face.InsideCVE()->Dimension() == 0U)
    {
        double vis_w_avg(0.);
        double vis_o_avg(0.);

        for (size_t n = 0U; n < face.InsideCVE()->Neighbors(); ++n)
        {
            flow_model.Initialize(*(face.InsideCVE()->Neighbor(n)->E()));
            flow_model.EffectiveSaturation();

            vis_w_avg += flow_model.ViscosityWettingPhase();
            vis_o_avg += flow_model.ViscosityNonWettingPhase();
            ro_w_in += flow_model.DensityWettingPhase();
            ro_o_in += flow_model.DensityNonWettingPhase();
            pc_in += flow_model.pc_Phase() * flow_model.Permeability();
            perm_in += flow_model.Permeability();

        }

        vis_w_avg /= face.InsideCVE()->Neighbors();
        vis_o_avg /= face.InsideCVE()->Neighbors();
        ro_w_in /= face.InsideCVE()->Neighbors();
        ro_o_in /= face.InsideCVE()->Neighbors();
        pc_in /= perm_in;
        perm_in /= face.InsideCVE()->Neighbors();

        if ((face.InsideCVE()->PropertyValue(1U) > 0.) && (face.InsideCVE()->PropertyValue(1U) < 1.))
        {
            lambda_o_in = face.InsideCVE()->PropertyValue(1U) / vis_o_avg;
            lambda_w_in = (1. - face.InsideCVE()->PropertyValue(1U)) / vis_w_avg;                    
        }
        else if (face.InsideCVE()->PropertyValue(1U) >= 1.)
        {
            lambda_o_in = 1. / vis_o_avg;
            lambda_w_in = 0.;
        }
        else
        {
            lambda_o_in = 0.;
            lambda_w_in = 1. / vis_w_avg;
        }

        flow_model.Initialize(*(face.OutsideCVE()->E()));
        flow_model.EffectiveSaturation();

        perm_out = flow_model.Permeability();
        ro_w_out = flow_model.DensityWettingPhase();
        ro_o_out = flow_model.DensityNonWettingPhase();
        lambda_w_out = flow_model.MobilityPhase(1U);
        lambda_o_out = flow_model.MobilityPhase(2U);
        pc_out = flow_model.pc_Phase();

    }
    else
    {
        flow_model.Initialize(*(face.InsideCVE()->E()));
        flow_model.EffectiveSaturation();

        perm_in = flow_model.Permeability();
        ro_w_in = flow_model.DensityWettingPhase();
        ro_o_in = flow_model.DensityNonWettingPhase();
        lambda_w_in = flow_model.MobilityPhase(1U);
        lambda_o_in = flow_model.MobilityPhase(2U);
        pc_in = flow_model.pc_Phase();

        double vis_w_avg(0.);
        double vis_o_avg(0.);

        for (size_t n = 0U; n < face.OutsideCVE()->Neighbors(); ++n)
        {
            flow_model.Initialize(*(face.OutsideCVE()->Neighbor(n)->E()));
            flow_model.EffectiveSaturation();

            vis_w_avg += flow_model.ViscosityWettingPhase();
            vis_o_avg += flow_model.ViscosityNonWettingPhase();
            ro_w_out += flow_model.DensityWettingPhase();
            ro_o_out += flow_model.DensityNonWettingPhase();
            pc_out += flow_model.pc_Phase() * flow_model.Permeability();
            perm_out += flow_model.Permeability();

        }

        vis_w_avg /= face.OutsideCVE()->Neighbors();
        vis_o_avg /= face.OutsideCVE()->Neighbors();
        ro_w_out /= face.OutsideCVE()->Neighbors();
        ro_o_out /= face.OutsideCVE()->Neighbors();
        pc_out/= perm_out;
        perm_out /= face.OutsideCVE()->Neighbors();

        if ((face.OutsideCVE()->PropertyValue(1U) > 0.) && (face.OutsideCVE()->PropertyValue(1U) < 1.))
        {
            lambda_o_out = face.OutsideCVE()->PropertyValue(1U) / vis_o_avg;
            lambda_w_out = (1. - face.OutsideCVE()->PropertyValue(1U)) / vis_w_avg;                    
        }
        else if (face.OutsideCVE()->PropertyValue(1U) >= 1.)
        {
            lambda_o_out = 1. / vis_o_avg;
            lambda_w_out = 0.;
        }
        else
        {
            lambda_o_out = 0.;
            lambda_w_out = 1. / vis_w_avg;
        }

    }

    const double capillary_grad(withCapillary_ ? face.PropertyValue(1U) * (pc_out - pc_in) : 0.);
    const double perm(2. / (1. / perm_out + 1. / perm_in));

    double lambda_w, lambda_o;
    double ro_w, ro_o;

    FindUpstreamProperties(face.PropertyValue(0U), face.Area(), face.UnitNormal(),
                           perm, lambda_w_out, lambda_o_out, lambda_w_in, lambda_o_in, 
                           ro_w_out, ro_o_out, ro_w_in, ro_o_in, capillary_grad,
                           lambda_w, lambda_o, ro_w, ro_o);

    // face vertical normal area
    const double fvna(face.UnitNormal()[1U] * face.Area());

    double oil_face_flux(face.PropertyValue(0U) / perm);

    oil_face_flux += lambda_w * capillary_grad * face.Area();
    oil_face_flux -= (lambda_w * ro_w + lambda_o * ro_o) * G_ * fvna;
    oil_face_flux /= (lambda_w + lambda_o + EPSILON_);
    oil_face_flux += ro_o * G_ * fvna;
    oil_face_flux *= lambda_o * perm;

    return oil_face_flux;

} // end CalculateOilFluxInteriorPointFace






template<uint32_t dim>
void TwoPhaseElementBasedTransport<dim>::TransportAcrossZeroDimensionalCVE_Faces(ControlVolumeElement<dim>& cve, TwoPhaseModel<dim>& flow_model)
{
    double vis_w_avg(0.);
    double vis_o_avg(0.);
    double ro_w_avg(0.);
    double ro_o_avg(0.);
    double pc_avg(0.);
    double perm_avg(0.);

    for (size_t n = 0U; n < cve.Neighbors(); ++n)
    {
        flow_model.Initialize(*(cve.Neighbor(n)->E()));
        flow_model.EffectiveSaturation();

        vis_w_avg += flow_model.ViscosityWettingPhase();
        vis_o_avg += flow_model.ViscosityNonWettingPhase();
        ro_w_avg += flow_model.DensityWettingPhase();
        ro_o_avg += flow_model.DensityNonWettingPhase();
        pc_avg += flow_model.pc_Phase() * flow_model.Permeability();
        perm_avg += flow_model.Permeability();
    }

    vis_w_avg /= cve.Neighbors();
    vis_o_avg /= cve.Neighbors();
    ro_w_avg /= cve.Neighbors();
    ro_o_avg /= cve.Neighbors();
    pc_avg /= perm_avg;
    perm_avg /= cve.Neighbors();

    double lambda_w_in;
    double lambda_o_in;

    if ((cve.PropertyValue(1U) > 0.) && (cve.PropertyValue(1U) < 1.))
    {
        lambda_o_in = cve.PropertyValue(1U) / vis_o_avg;
        lambda_w_in = (1. - cve.PropertyValue(1U)) / vis_w_avg;                    
    }
    else if (cve.PropertyValue(1U) >= 1.)
    {
        lambda_o_in = 1. / vis_o_avg;
        lambda_w_in = 0.;
    }
    else
    {
        lambda_o_in = 0.;
        lambda_w_in = 1. / vis_w_avg;
    }

    for (size_t n = 0U; n < cve.Neighbors(); ++n)
    {
        flow_model.Initialize(*(cve.Neighbor(n)->E()));
        flow_model.EffectiveSaturation();
        const double perm_out(flow_model.Permeability());
        const double ro_w_out(flow_model.DensityWettingPhase());
        const double ro_o_out(flow_model.DensityNonWettingPhase());
        const double lambda_w_out(flow_model.MobilityPhase(1U));
        const double lambda_o_out(flow_model.MobilityPhase(2U));
        const double pc_out(flow_model.pc_Phase());

        const double capillary_grad(withCapillary_ ? cve.Face(n)->PropertyValue(1U) * cve.FaceNormalDirection(n) * (pc_out - pc_avg) : 0.);
        const double perm(2. / (1. / perm_out + 1. / perm_avg));

        double lambda_w, lambda_o;
        double ro_w, ro_o;

        FindUpstreamProperties(cve.Face(n)->PropertyValue(0U) * cve.FaceNormalDirection(n), 
                               cve.Face(n)->Area(), cve.Face(n)->UnitNormal() * cve.FaceNormalDirection(n),
                               perm, lambda_w_out, lambda_o_out, lambda_w_in, lambda_o_in, 
                               ro_w_out, ro_o_out, ro_w_avg, ro_o_avg, capillary_grad,
                               lambda_w, lambda_o, ro_w, ro_o);

        // face vertical normal area
        const double fvna(cve.Face(n)->UnitNormal()[1U] * cve.FaceNormalDirection(n) * cve.Face(n)->Area());

        double oil_face_flux(cve.Face(n)->PropertyValue(0U) * cve.FaceNormalDirection(n) / perm);

        oil_face_flux += lambda_w * capillary_grad * cve.Face(n)->Area();
        oil_face_flux -= (lambda_w * ro_w + lambda_o * ro_o) * G_ * fvna;
        oil_face_flux /= (lambda_w + lambda_o + EPSILON_);
        oil_face_flux += ro_o * G_ * fvna;
        oil_face_flux *= lambda_o * perm;

        cve.PropertyValueAdd(0U, oil_face_flux);
        cve.Neighbor(n)->PropertyValueAdd(0U, -oil_face_flux);

    } // end looping over neighbors

} // end TransportAcrossZeroDimensionalCVE_Faces




template<uint32_t dim>
double TwoPhaseElementBasedTransport<dim>::Transport(double total_time_increment,
                                                       TwoPhaseModel<dim>& flow_model)
{
    std::cout << "\nTwoPhaseElementBasedTransport<" << dim << ">::Transport()\n";

    double normal_changes(setup_.NormalSaturationChanges());
    double max_changes(setup_.MaximumSaturationChanges());
    double max_internal_time_step_size(setup_.MaximumTransportTimestepSize());
    size_t maximum_number_of_transport_steps(setup_.MaximumTransportTimesteps());
    double CFL_multiplier(setup_.CFL_Multiplier());
    bool with_divergence_correction(setup_.WithDivergenceCorrection());

    ConstructFluxes();

    typedef typename std::vector<ElementFace<dim> >::iterator face_itr;
    typedef typename std::vector<ControlVolumeElement<dim> >::iterator CVE_itr;

    std::pair<double, double> CFL_timestep_divergence(CFL_Timestep_MaxDivergence());
    const double CFL_timestep(CFL_timestep_divergence.first * CFL_multiplier);
    std::cout.precision(4);
    std::cout << std::scientific;
    std::cout << std::fixed;
    std::cout << "\nCFL Time step Criteria: " << CFL_timestep << " sec\n";
    std::cout << "Max Divergence: " << CFL_timestep_divergence.second << "\n";
    std::cout << std::fixed;
    std::cout << "Transport Interval: " << total_time_increment << " sec\n";
    std::cout << "Transport Maximum Number of Steps: " << maximum_number_of_transport_steps << "\n";
    std::cout << "\n";

    static double internal_time_increment(total_time_increment);
    internal_time_increment = internal_time_increment < CFL_timestep ? internal_time_increment : CFL_timestep;
    internal_time_increment = internal_time_increment > total_time_increment ? total_time_increment : internal_time_increment;
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
            transportModelRef_.ZeroCVEsProperty(0U);

            // transport for inner line faces
            for (face_itr fit = transportModelRef_.LineFacesBegin(); fit != transportModelRef_.PointFacesBegin(); ++fit)
            {                
                const double oil_face_flux(CalculateOilFluxInteriorLineFace(*fit, flow_model));
                fit->InsideCVE()->PropertyValueAdd(0U, oil_face_flux);
                fit->OutsideCVE()->PropertyValueAdd(0U, -oil_face_flux);

            } // end for loop transport for inner line faces

            // transport for boundary line faces
            for (face_itr fit = transportModelRef_.PerimeterLineFacesBegin(); fit != transportModelRef_.PerimeterPointFacesBegin(); ++fit)
            {                
                const double oil_face_flux(CalculateOilFluxBoundaryLineFace(*fit, flow_model));
                fit->InsideCVE()->PropertyValueAdd(0U, oil_face_flux);

            } // end for loop transport for boundary line faces

            // looping over interior point CVEs
            for (CVE_itr cveit = transportModelRef_.PointCVEsBegin(); cveit != transportModelRef_.PerimeterCVEsBegin(); ++cveit)
            {
                TransportAcrossZeroDimensionalCVE_Faces(*cveit, flow_model);
            } // end looping over interior point CVEs

            maximum_change_in_timestep = 0.;

            // looping over inside 2D and 1D cves
            for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = transportModelRef_.CVEsBegin();
                cveit != transportModelRef_.PointCVEsBegin(); ++cveit)
            {
                flow_model.Initialize(*(cveit->E()));
                flow_model.EffectiveSaturation();
                const double lambda_w(flow_model.MobilityPhase(1U));
                const double lambda_o(flow_model.MobilityPhase(2U));

                const double oil_source_sink(cveit->Volume() * (cveit->E()->Read(oilInjRateKey_) +
                                               cveit->E()->Read(totalProdRateKey_) * 
                                               lambda_o / (lambda_w + lambda_o)));

                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes(cveit->PropertyValue(0U) + oil_source_sink);
                // adding flux mismatch
                if (with_divergence_correction)
                    changes -= cveit->E()->Read(divergenceKey_) * lambda_o / (lambda_w + lambda_o);
                // calculating element saturation changes for the timestep
                changes *= (internal_time_increment / cveit->Volume() / cveit->E()->Read(porosityKey_));

                cveit->PropertyValue(0U, changes);

                maximum_change_in_timestep = ((std::abs(changes) > maximum_change_in_timestep)
                    ? std::abs(changes) : maximum_change_in_timestep);
            }

            // looping over inside 0D cves
            for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = transportModelRef_.PointCVEsBegin();
                cveit != transportModelRef_.PerimeterCVEsBegin(); ++cveit)
            {
                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes(cveit->PropertyValue(0U));
                // assuming porosity is 1.
                changes *= (internal_time_increment / cveit->Volume());

                cveit->PropertyValue(0U, changes);

                maximum_change_in_timestep = ((std::abs(changes) > maximum_change_in_timestep)
                    ? std::abs(changes) : maximum_change_in_timestep);
            }

            // looping over boundary 2D and 1D cves
            for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = transportModelRef_.PerimeterCVEsBegin();
                cveit != transportModelRef_.PerimeterPointCVEsBegin(); ++cveit)
            {
                flow_model.Initialize(*(cveit->E()));
                flow_model.EffectiveSaturation();
                const double lambda_w(flow_model.MobilityPhase(1U));
                const double lambda_o(flow_model.MobilityPhase(2U));

                const double oil_source_sink(cveit->Volume() * (cveit->E()->Read(oilInjRateKey_) +
                    cveit->E()->Read(totalProdRateKey_) * 
                    lambda_o / (lambda_w + lambda_o)));

                // calculating oil saturation changes by integrating oil fluxes over the element
                double changes(cveit->PropertyValue(0U) + oil_source_sink);
                // adding flux mismatch
                if (with_divergence_correction)
                    changes -= cveit->E()->Read(divergenceKey_) * lambda_o / (lambda_w + lambda_o);
                // calculating element saturation changes for the timestep
                changes *= (internal_time_increment / cveit->Volume() / cveit->E()->Read(porosityKey_));

                cveit->PropertyValue(0U, changes);

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


        // calculating interior element saturation
        for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = transportModelRef_.CVEsBegin();
            cveit != transportModelRef_.PointCVEsBegin(); ++cveit)
        {
            if (cveit->E()->Status(oilSaturationKey_) != DIRICH)
            {
                double saturation(cveit->E()->Read(oilSaturationKey_));
                saturation += cveit->PropertyValue(0U);
                cveit->E()->Store(oilSaturationKey_, ScalarVariable(PLAIN, saturation));
                cveit->E()->Store(waterSaturationKey_, ScalarVariable(PLAIN, 1. - saturation));
            }

        }

        // calculating interior 0D cves saturation
        for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = transportModelRef_.PointCVEsBegin();
            cveit != transportModelRef_.PerimeterCVEsBegin(); ++cveit)
        {
            cveit->PropertyValue(1U, cveit->PropertyValue(1U) + cveit->PropertyValue(0U));
        }

        // calculating boundary element saturation
        for (typename std::vector<ControlVolumeElement<dim> >::iterator cveit = transportModelRef_.PerimeterCVEsBegin();
            cveit != transportModelRef_.PerimeterPointCVEsBegin(); ++cveit)
        {
            if (cveit->E()->Status(oilSaturationKey_) != DIRICH)
            {
                double saturation(cveit->E()->Read(oilSaturationKey_));
                saturation += cveit->PropertyValue(0U);
                cveit->E()->Store(oilSaturationKey_, ScalarVariable(PLAIN, saturation));
                cveit->E()->Store(waterSaturationKey_, ScalarVariable(PLAIN, 1. - saturation));
            }

        }

        total_time += previous_time_step;
        ++transport_step_number;

        if ((internal_time_increment > (total_time_increment - total_time)) && ((total_time_increment - total_time) > 0.))
            internal_time_increment = (total_time_increment - total_time);

        internal_time_increment = ((internal_time_increment > max_internal_time_step_size)
            ? max_internal_time_step_size : internal_time_increment);

        if (sig.SignalRaised()) std::raise(SIGINT);

        if (sig.OutputSignal() || sig.RestartSignal() || sig.QuitSignal()) exiting_signal_received = true;

    } while ((total_time < total_time_increment) && (transport_step_number < maximum_number_of_transport_steps) && !exiting_signal_received); // end of outer do loop to achieve total_time_increment 

    std::cout << std::scientific;

    return total_time;

} // end Transport





template<uint32_t dim>
void TwoPhaseElementBasedTransport<dim>::Divergence(double time)
{
    std::cout << "\nTwoPhaseElementBasedTransport<" << dim << "> Calculating the divergence of total velocity field.\n";

    ConstructFluxes();

    // looping over interior 2D and 1D elements
    for (typename std::vector<ControlVolumeElement<dim> >::iterator eit = transportModelRef_.CVEsBegin();
        eit != transportModelRef_.PointCVEsBegin(); ++eit)
    {
        double flux_imbalance(0);

        for (size_t i = 0; i < eit->Faces(); ++i)
        {
            flux_imbalance += eit->Face(i)->PropertyValue(0U) * eit->FaceNormalDirection(i);
        }

        flux_imbalance += (eit->E()->Read(totalProdRateKey_) + 
                           eit->E()->Read(oilInjRateKey_) + 
                           eit->E()->Read(waterInjRateKey_)) * 
                           eit->Volume();

        eit->E()->Store(divergenceKey_, ScalarVariable(PLAIN, flux_imbalance));

    }

    // looping over boundary 2D and 1D elements
    for (typename std::vector<ControlVolumeElement<dim> >::iterator eit = transportModelRef_.PerimeterCVEsBegin();
        eit != transportModelRef_.PerimeterPointCVEsBegin(); ++eit)
    {
        double flux_imbalance(0);

        for (size_t i = 0; i < eit->Faces(); ++i)
        {
            flux_imbalance += eit->Face(i)->PropertyValue(0U) * eit->FaceNormalDirection(i);
        }

        flux_imbalance += (eit->E()->Read(totalProdRateKey_) + 
                           eit->E()->Read(oilInjRateKey_) + 
                           eit->E()->Read(waterInjRateKey_)) * 
                           eit->Volume();

        eit->E()->Store(divergenceKey_, ScalarVariable(PLAIN, flux_imbalance));

    }

    VTK_Interface<dim> vtk_output;

    vtk_output.OutputDataToVTK(modelRef_, "divergence", "divergence", static_cast<size_t>(time), false);


} // end Divergence


template class TwoPhaseElementBasedTransport<2U>;

} // end namespace csmp
