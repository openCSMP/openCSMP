// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "compute2PhaseMobilityAtBaryCenter.h"
#include "Model.h"
#include "Region.h"
#include "TwoPhaseModel.h"

using namespace std;
namespace csmp
{
template<uint32_t dim,template<uint32_t> class CELL >
void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<dim,CELL>& sg,
                                        TwoPhaseModel<dim>& relperm)
{
  ScalarVariable  sc, mob_t;

  //Index saturationOilKey=relperm.NonWettingPhaseSaturationKey();
  //Index saturationWaterKey=relperm.WettingPhaseSaturationKey();
  Index totalMobilityKey=relperm.TotalMobilityKey();

  // 0. Computing the saturation of water = 1 - So
  // ---------------------------------------------
  /*const typename vector<Node<dim>*>::iterator nodesend=sg.NodesEnd();
  for ( typename vector<Node<dim>*>::iterator
        it=sg.NodesBegin(); it!=nodesend; ++it )
  {
    sc = 1. - (*it)->Read( saturationOilKey );
    (*it)->Store(  saturationWaterKey, sc );
  }
  */

  const auto elementsend=sg.CellsEnd();
  // 2. Computing the multiphase flow properties
  for ( auto it=sg.CellsBegin(); it!=elementsend; ++it )
  {
    // 0. setting up the relative permeability model
    // ---------------------------------------------
    relperm.Initialize( *(*it) );
    relperm.InitializeForBaryCenter( *(*it) );
    relperm.EffectiveSaturation();

    // 1. k * total mobility
    // -------------------------------------------------------------------------------------------
    mob_t = relperm.TotalMobility();
    (*it)->Store(  totalMobilityKey, mob_t );


    // relperm.Out();
  }

} // end compute2PhaseFlowPropertiesAtBaryCenter



template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<1U,Element>&, TwoPhaseModel<1U>&);
template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<2U,Element>&, TwoPhaseModel<2U>&);
template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<3U,Element>&, TwoPhaseModel<3U>&);

///Roman, 2014: should TwoPhaseModel have methods for Faces and InterFaces?
//template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<1U,Face>&, TwoPhaseModel<1U>&);
//template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<2U,Face>&, TwoPhaseModel<2U>&);
//template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<3U,Face>&, TwoPhaseModel<3U>&);

//template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<1U,InterFace>&, TwoPhaseModel<1U>&);
//template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<2U,InterFace>&, TwoPhaseModel<2U>&);
//template void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<3U,InterFace>&, TwoPhaseModel<3U>&);
}
