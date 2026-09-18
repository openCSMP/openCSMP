// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PorosityPermeabilityCouplingVisitor.h"

using namespace std;

namespace csmp
{
template<uint32_t dim>
PorosityPermeabilityCouplingVisitor<dim>::PorosityPermeabilityCouplingVisitor(
    Model<dim> &model, bool with_pore_volume_fields)
    : Visitor<dim>( ELEMENT ),
    model                     ( model ),
    permeability_key          ( model.Database().StorageKey("permeability") ),
    porosity_key              ( model.Database().StorageKey("porosity") ),
    // The next five are only used by CalculatePoreVolumeChangeFactor().
    // When with_pore_volume_fields is false, they are left as default-
    // constructed Index objects and re-looked-up if that method is called.
    porosity_node_key         ( with_pore_volume_fields
                          ? model.Database().StorageKey("nodal porosity") : Index() ),
    bulk_volume_key           ( with_pore_volume_fields
                        ? model.Database().StorageKey("bulk volume") : Index() ),
    pore_volume_before_key    ( with_pore_volume_fields
                               ? model.Database().StorageKey("pore volume before") : Index() ),
    pore_volume_key           ( with_pore_volume_fields
                        ? model.Database().StorageKey("pore volume") : Index() ),
    volume_change_factor_key  ( with_pore_volume_fields
                                 ? model.Database().StorageKey("pore volume change factor") : Index() ),

    depth_dependent   (true),
    coupling          ("VermaPruess"),

    phi_min           (0.03), // lower limit porosity
    phi_max           (0.1)   // upper limit porosity


{
    with_pore_volume_fields_ = with_pore_volume_fields;
}

template<uint32_t dim>
PorosityPermeabilityCouplingVisitor<dim>::~PorosityPermeabilityCouplingVisitor()
{

}

template<uint32_t dim>
void PorosityPermeabilityCouplingVisitor<dim>::Visit( Model<dim> *m )
{

}

template<uint32_t dim>
void PorosityPermeabilityCouplingVisitor<dim>::Visit( Region<dim> *r )
{

}

template<uint32_t dim>
void PorosityPermeabilityCouplingVisitor<dim>::Visit(Element<dim> *e)
{
    phi_D = 0;

    // Calculate depth-dependent permeability-porosity coupling
    if (depth_dependent)
    {
        if (coupling == "Costa")
        {
            phi_D = CalculateDepthDependentPorosityPermeabilityCoupling_Costa(e);
        }
        else if (coupling == "VermaPruess")
        {
            phi_D = CalculateDepthDependentPorosityPermeabilityCoupling_VermaPruess(e);
        }
    }
    // include more permeability-porosity coupling here (e.g., temperature-dependent)

    // select 'dominating' coupling to store on element
    phi = phi_D;

    // apply upper and lower limits
    phi = std::min(phi(), phi_max);
    phi = std::max(phi(), phi_min);

    // Store data
    e->Store(porosity_key, phi);

} // end visit

// Coupling of porosity to given depth-dependent permeability following Eq. 5 (Fig. 3) of Costa 2006 (doi.org/10.1029/2005GL025134)
template<uint32_t dim>
double PorosityPermeabilityCouplingVisitor<dim>::CalculateDepthDependentPorosityPermeabilityCoupling_Costa( Element<dim> *e )
{
    double a, b;
    e->Read(permeability_key, k);

    // Coupling of porosity to given depth-dependent permeability following Eq. 5 (Fig. 3) of Costa 2006 (doi.org/10.1029/2005GL025134)
    // We apply a log-log fit to their data where k=1E-15 <=> phi=0.1 and k=3E-17 <=> phi=0.03
    b = (log(0.1) - log(0.03)) / (log(1e-15) - log(3e-17));
    a = log(0.03) - b * log(3e-17);

    return exp(a + b * log(k()));
}

// Coupling or porosity using a Verma-Pruess relationship
template<uint32_t dim>
double PorosityPermeabilityCouplingVisitor<dim>::CalculateDepthDependentPorosityPermeabilityCoupling_VermaPruess( Element<dim> *e )
{
    // old version
    // double k_ref(1.E-14), phi_ref(0.25), phi_crit(0.03), n(2.);

    // This version fits k = [1E-16, 1E-15, 1E-14, 1E-13] to phi = [0.030, 0.05, 0.1, 0.215]
    double k_ref(1.E-14), phi_ref(0.0992), phi_crit(0.0166), n(2.633);
    e->Read(permeability_key, k);

    return phi_crit + (phi_ref-phi_crit) * std::pow( k()/k_ref, 1.0/n );
}

template<uint32_t dim>
void PorosityPermeabilityCouplingVisitor<dim>::CalculatePoreVolumeChangeFactor()
{
    if (!with_pore_volume_fields_) {
        throw csmp::Exception(
            FATAL_ERROR,
            "PorosityPermeabilityCouplingVisitor::CalculatePoreVolumeChangeFactor",
            "this visitor was constructed with with_pore_volume_fields=false; the "
            "pore-volume fields were not looked up. Construct with the default (true) "
            "to use this method.");
    }

    ScalarVariable volume_change_factor, bulk_volume, pore_volume_before, phi_node,
        pore_volume_, pore_volume;

    typename vector<Node<dim>*>::const_iterator node;

    for (node = model.Region("Model").NodesBegin();
         node != model.Region("Model").NodesEnd();
         node++)
    {
        (*node)->Read(bulk_volume_key, bulk_volume);
        (*node)->Read(pore_volume_before_key, pore_volume_before);
        (*node)->Read(porosity_node_key, phi_node);

        pore_volume = bulk_volume() * phi_node();

        volume_change_factor() = pore_volume() / pore_volume_before();

        (*node)->Store(volume_change_factor_key, volume_change_factor);
        // (*node)->Store(pore_volume_key, pore_volume); // Calculated using the PoreVolumeVisitor
    }
}

template<uint32_t dim>
void PorosityPermeabilityCouplingVisitor<dim>::DepthDependent( bool d_dep )
{
    depth_dependent = d_dep;
}

template<uint32_t dim>
void PorosityPermeabilityCouplingVisitor<dim>::CouplingOption( std::string coupling_option )
{
    coupling = coupling_option;

    if (coupling!="Costa" && coupling!="VermaPruess")
    {
        throw csmp::Exception( FATAL_ERROR, "PorosityPermeabilityCouplingVisitor::Visit",
                              "Invalid permeability-porosity-coupling option: " + coupling +
                                  "\nAvailable options: 'Costa' or 'VermaPruess'. " );
    }
}

template<uint32_t dim>
void PorosityPermeabilityCouplingVisitor<dim>::ChangeUpperLowerPorosityLimits(
    double lower_limit, double upper_limit  )
{
    phi_min = lower_limit;
    phi_max = upper_limit;
}


template class PorosityPermeabilityCouplingVisitor<1U>;
template class PorosityPermeabilityCouplingVisitor<2U>;
template class PorosityPermeabilityCouplingVisitor<3U>;

} // end namespace csmp
