#include "compute2PhaseMobilityAtBaryCenter.h"
#include "Model.h"
#include "TwoPhaseModel.h"

using namespace std;
namespace csmp
{
template<size_t dim,template<size_t > class SIMPLEX >
void compute2PhaseMobilityAtBaryCenter( ModelSubDomain<dim,SIMPLEX>& sg,
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

  const typename vector<SIMPLEX<dim>*>::iterator elementsend=sg.ElementsEnd();
  // 2. Computing the multiphase flow properties
  for ( typename vector<SIMPLEX<dim>*>::iterator
        it=sg.ElementsBegin(); it!=elementsend; ++it )
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
