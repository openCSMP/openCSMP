#include "IsnanIsinf_Test.h"

using namespace std;

namespace csmp
{

IsnanIsinf_Test::IsnanIsinf_Test()
{
}
	
IsnanIsinf_Test::~IsnanIsinf_Test()
{
}
	
void IsnanIsinf_Test::run()
{
	double zero = 0.0;

	// isnan tests
	_info("isnan(0.0) = "	<< isnan(0.0));
	_test( !isnan(0.0) );

	_info("isnan( 1.0/0.0 ) = "<< isnan(1.0/zero));
	_test( !isnan(1.0/zero) );

	_info("isnan( -1.0/0.0 ) = "<< isnan(-1.0/zero));
	_test( !isnan(-1.0/zero) );

	_info("isnan( sqrt(-1.0) ) = "<< isnan(sqrt(-1.0)));
	_test( isnan(sqrt(-1.0)) );

	// isinf tests

	_info("isinf( 0.0 ) = "<< isinf(0.0));
	_test( !isinf(0.0) );

	_info("isinf( 1.0/0.0 ) = "<< isinf(1.0/zero));
	_test( isinf(1.0/zero) );

	_info("isinf( -1.0/0.0 ) = "<< isinf(-1.0/zero));
	_test( isinf(-1.0/zero) );

	_info("isinf( sqrt(-1.0) ) = "<< isinf(sqrt(-1.0)));
	_test( !isinf(sqrt(-1.0)) );

	// isfinite tests
	_info("isfinite( 0.0 ) = "<< isfinite(0.0));
	_test( isfinite(0.0) );

	_info("isfinite( 1.0/0.0 ) = "<< isfinite(1.0/zero));
	_test( !isfinite(1.0/zero) );

	_info("isfinite( -1.0/0.0 ) = "<< isfinite(-1.0/zero));
	_test( !isfinite(-1.0/zero) );

	_info("isfinite( sqrt(-1.0) ) = "	<< isfinite(sqrt(-1.0)));
	_test( !isfinite(sqrt(-1.0)) );
}



} //end namespace csmp
