#include "FlowFunctionsModule.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

template<size_t dim>
FlowFunctionsModule1<dim>::FlowFunctionsModule1( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
 {
 }

template<size_t dim>
FlowFunctionsModule2<dim>::FlowFunctionsModule2( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
 {
 }

template<size_t dim>
FlowFunctionsModule3<dim>::FlowFunctionsModule3( PropertyDatabase<dim>& db, const char* model_name )
 : variables::VariableSet_CO2GeoSequestration(db), ExperimentalSaturationFunctions<dim, csmp::FlowFunctionsModule3>( (string(model_name) + "-rock_types.txt").c_str() )
 {
    // setting the maximum number of rock types in the database to the number of records red from the -rock_types.txt file
    db.SetRangeOf( "rocktype", 0., this->RockTypes() );
 }
 
 template<size_t dim>
FlowFunctionsModule4<dim>::FlowFunctionsModule4( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
 {
 }

template<size_t dim>
FlowFunctionsModule5<dim>::FlowFunctionsModule5( const PropertyDatabase<dim>& db )
 : variables::VariableSet_CO2GeoSequestration(db)
 {
 }

template<size_t dim>
FlowFunctionsModule6<dim>::FlowFunctionsModule6( PropertyDatabase<dim>& db, const char* model_name )
 : variables::VariableSet_CO2GeoSequestration(db), ExperimentalSaturationFunctions<dim, csmp::FlowFunctionsModule6>( (string(model_name) + "-rock_types.txt").c_str() )
 {
    // setting the maximum number of rock types in the database to the number of records red from the -rock_types.txt file
    db.SetRangeOf( "rocktype", 0., this->RockTypes() );
 }


template class FlowFunctionsModule1<1U>;
template class FlowFunctionsModule1<2U>;
template class FlowFunctionsModule1<3U>;

template class FlowFunctionsModule2<1U>;
template class FlowFunctionsModule2<2U>;
template class FlowFunctionsModule2<3U>;

template class FlowFunctionsModule3<1U>;
template class FlowFunctionsModule3<2U>;
template class FlowFunctionsModule3<3U>;

template class FlowFunctionsModule4<1U>;
template class FlowFunctionsModule4<2U>;
template class FlowFunctionsModule4<3U>;

template class FlowFunctionsModule5<1U>;
template class FlowFunctionsModule5<2U>;
template class FlowFunctionsModule5<3U>;

template class FlowFunctionsModule6<1U>;
template class FlowFunctionsModule6<2U>;
template class FlowFunctionsModule6<3U>;


} // end namespace csmp





