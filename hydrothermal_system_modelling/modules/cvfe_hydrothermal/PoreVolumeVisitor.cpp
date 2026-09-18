// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "PoreVolumeVisitor.h"
#include "Exception.h"
#include "Model.h"

using namespace std;

namespace csmp {

/** custom constructor */
template<uint32_t dim>
PoreVolumeVisitor<dim>::PoreVolumeVisitor(Model<dim>& model,
                                          const char* porosity,
                                          const char* volume,
                                          const char* pore_volume,
                                          const char* thickness)
{

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);
    property_key_ = model.Database().StorageKey(pore_volume);
    phi_key_      = model.Database().StorageKey(porosity);
    vol_key_      = model.Database().StorageKey(volume);
    thickness_key_= model.Database().StorageKey(thickness);
    output_solid_volume = false;

    if ( property_key_.type != SCALAR || property_key_.place != NODE )
        throw csmp::Exception( ERROR, "PoreVolumeVisitor::(constructor)",
                              pore_volume, " must be a nodal scalar property." );

    if ( phi_key_.type != SCALAR || phi_key_.place != NODE )
        throw csmp::Exception( ERROR, "PoreVolumeVisitor::(constructor)",
                              porosity, " must be a scalar property." );

    if ( vol_key_.type != SCALAR || vol_key_.place != NODE )
        throw csmp::Exception( ERROR, "PoreVolumeVisitor::(constructor)",
                              volume, " must be a scalar property." );
}

template<uint32_t dim>
PoreVolumeVisitor<dim>::PoreVolumeVisitor(Model<dim>& model,
                                          const char* porosity,
                                          const char* volume,
                                          const char* pore_volume,
                                          const char* solid_volume,
                                          const char* thickness)
{

    this->ApplicationLevel(MODEL);
    this->ApplicationTarget(NODE);
    property_key_ = model.Database().StorageKey(pore_volume);
    phi_key_      = model.Database().StorageKey(porosity);
    vol_key_      = model.Database().StorageKey(volume);
    solid_vol_key_= model.Database().StorageKey(solid_volume);
    thickness_key_= model.Database().StorageKey(thickness);
    output_solid_volume = true;

    if ( property_key_.type != SCALAR || property_key_.place != NODE )
        throw csmp::Exception( ERROR, "PoreVolumeVisitor::(constructor)",
                              pore_volume, " must be a nodal scalar property." );

    if ( phi_key_.type != SCALAR || phi_key_.place != NODE )
        throw csmp::Exception( ERROR, "PoreVolumeVisitor::(constructor)",
                              porosity, " must be a scalar property." );

    if ( vol_key_.type != SCALAR || vol_key_.place != NODE )
        throw csmp::Exception( ERROR, "PoreVolumeVisitor::(constructor)",
                              volume, " must be a scalar property." );

    if ( solid_vol_key_.type != SCALAR || solid_vol_key_.place != NODE )
        throw csmp::Exception( ERROR, "PoreVolumeVisitor::(constructor)",
                              solid_volume, " must be a scalar property." );
}

/** deconstructor */
template<uint32_t dim>
PoreVolumeVisitor<dim>::~PoreVolumeVisitor()
{
}

/** visit function for Node */
template<uint32_t dim>
void PoreVolumeVisitor<dim>::Visit(Node<dim>* n)
{

    n->Read( phi_key_, phi );

    vol()       = 0.;
    pore_vol()  = 0.;
    solid_vol() = 0.;

    for (unsigned int i = 0; i < n->Parents(); i++)
    {
        volume = n->Parent(i)->Volume();

        // Volume() returns the element's own measure: a volume for a volumetric
        // element, an AREA for a surface element, a LENGTH for a line element.
        // Multiplying by `thickness` turns a surface element's area into a proper
        // volume, which is the point for lower-dimensional fracture/fault regions.
        //
        // TODO (Benoit DD/MM/YYYY): line elements go through the same path, so a
        // well node picks up length x thickness / 2 per well line element — an m^2
        // quantity added to an m^3 accumulator. The error is small (tens of m^2
        // against control volumes of 1e3-1e4 m^3) and it reaches the well index
        // only through sqrt() and then log(), so it is negligible in practice.
        // It is still wrong, and a well node's bulk volume should probably count
        // only its reservoir parents. Skipping IsLine() parents here would fix it
        // — but check first what else reads "bulk volume" on well nodes before
        // changing a quantity this widely used.
        volume *= n->Parent(i)->Read(thickness_key_);

        volume /= n->Parent(i)->Nodes();
        vol() += volume;
    }

    pore_vol() = vol()*phi();
    n->Store( property_key_, pore_vol );
    n->Store( vol_key_, vol );

    if (output_solid_volume )
    {
        solid_vol()=vol()-pore_vol();
        n->Store( solid_vol_key_, solid_vol );
    }
}

/** visit function for Region */
template<uint32_t dim>
void PoreVolumeVisitor<dim>::Visit(Model<dim>* n)
{
    // no calcution for the region
}

template class PoreVolumeVisitor<1U>;
template class PoreVolumeVisitor<2U>;
template class PoreVolumeVisitor<3U>;

} // csmp