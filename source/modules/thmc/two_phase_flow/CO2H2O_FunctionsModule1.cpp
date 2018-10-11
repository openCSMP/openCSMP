#include "CO2H2O_FunctionsModule1.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<size_t dim>
CO2H2O_FunctionsModule0<dim>::CO2H2O_FunctionsModule0( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db),
   Fluid<dim,CO2H2O_FunctionsModule0> (false)
 {
 }

template<size_t dim>
CO2H2O_FunctionsModule1<dim>::CO2H2O_FunctionsModule1( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db),
   Fluid<dim,CO2H2O_FunctionsModule1> (false)
 {
 }

template<size_t dim>
CO2H2O_FunctionsModule2<dim>::CO2H2O_FunctionsModule2( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db),
   Fluid<dim,CO2H2O_FunctionsModule2> (false) 
 {
 }

template class CO2H2O_FunctionsModule0<1U>;
template class CO2H2O_FunctionsModule0<2U>;
template class CO2H2O_FunctionsModule0<3U>;

template class CO2H2O_FunctionsModule1<1U>;
template class CO2H2O_FunctionsModule1<2U>;
template class CO2H2O_FunctionsModule1<3U>;

template class CO2H2O_FunctionsModule2<1U>;
template class CO2H2O_FunctionsModule2<2U>;
template class CO2H2O_FunctionsModule2<3U>;


} // end namespace csmp





