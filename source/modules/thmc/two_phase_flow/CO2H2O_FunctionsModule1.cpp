#include "CO2H2O_FunctionsModule1.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<size_t dim>
CO2H2O_FunctionsModule0<dim>::CO2H2O_FunctionsModule0( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
 {
 }

template<size_t dim>
CO2H2O_FunctionsModule1<dim>::CO2H2O_FunctionsModule1( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
 {
 }

template<size_t dim>
CO2H2O_FunctionsModule2<dim>::CO2H2O_FunctionsModule2( PropertyDatabase<dim>& db, const char* model_name )
 : variables::VariableSet_CO2GeoSequestration(db), ExperimentalSaturationFunctions<dim, csmp::CO2H2O_FunctionsModule2>( (string(model_name) + "-rock_types.txt").c_str() )
 {
    // setting the maximum number of rock types in the database to the number of records red from the -rock_types.txt file
    db.SetRangeOf( "rocktype", 0., this->RockTypes() );
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





