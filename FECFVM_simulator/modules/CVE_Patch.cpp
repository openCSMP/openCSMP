// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "CVE_Patch.h"
#include "Node.h"


namespace csmp {




template<uint32_t dim>
CVE_Patch<dim>::CVE_Patch(Node<dim>& nd)
 : centerNode_ptr_(&nd),
   interiorFaces_(0),
   boundaryFaces_(0),
   active_(false)
{
}




template<uint32_t dim>
CVE_Patch<dim>::~CVE_Patch()
{
}




template<uint32_t dim>
void CVE_Patch<dim>::CreateAndAddSector(ControlVolumeElement<dim>& cvelm, double lower_dimensional_width)
{
    sectorsCollection_.push_back(CVE_Sector<dim> (cvelm, *centerNode_ptr_, lower_dimensional_width));
}




template<uint32_t dim>
double CVE_Patch<dim>::Volume()
{
    double volume(0.);
    
    for (typename std::vector<CVE_Sector<dim> >::iterator csi = SectorsBegin(); csi != SectorsEnd(); csi++)
        volume += csi->Volume();

    return volume;
}




template class CVE_Patch<2U>;

} // end namespace csmp
