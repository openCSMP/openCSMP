#include "FlowFunctions1.h"
#include "ErrorHandler.h"
#include "FiniteElementPlacement.h"
#include "FiniteVolumePlacement.h"

using namespace std;

namespace csmp {

template<size_t dim>
FlowFunctions1<dim>::FlowFunctions1( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
 {
 }

template<size_t dim>
FlowFunctions2<dim>::FlowFunctions2( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
//   ExperimentalSaturationFunctions<dim,FlowFunctions2>("rocktypes")
 {
 }

template class FlowFunctions1<1U>;
template class FlowFunctions1<2U>;
template class FlowFunctions1<3U>;

template class FlowFunctions2<1U>;
template class FlowFunctions2<2U>;
template class FlowFunctions2<3U>;


} // end namespace csmp





