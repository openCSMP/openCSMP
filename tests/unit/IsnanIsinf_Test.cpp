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
	double64 zero = 0.0;

	// isnan tests
	getInfoStream() << "isnan(0.0) = "	<< isnan(0.0);
	_test( !isnan(0.0) );

	getInfoStream() << "isnan( 1.0/0.0 ) = "<< isnan(1.0/zero);
	_test( !isnan(1.0/zero) );

	getInfoStream() << "isnan( -1.0/0.0 ) = "<< isnan(-1.0/zero);
	_test( !isnan(-1.0/zero) );

	getInfoStream() << "isnan( sqrt(-1.0) ) = "<< isnan(sqrt(-1.0));
	_test( isnan(sqrt(-1.0)) );

	// isinf tests

	getInfoStream() << "isinf( 0.0 ) = "<< isinf(0.0);
	_test( !isinf(0.0) );

	getInfoStream() << "isinf( 1.0/0.0 ) = "<< isinf(1.0/zero);
	_test( isinf(1.0/zero) );

	getInfoStream() << "isinf( -1.0/0.0 ) = "<< isinf(-1.0/zero);
	_test( isinf(-1.0/zero) );

	getInfoStream() << "isinf( sqrt(-1.0) ) = "<< isinf(sqrt(-1.0));
	_test( !isinf(sqrt(-1.0)) );

	// isfinite tests
	getInfoStream() << "isfinite( 0.0 ) = "<< isfinite(0.0);
	_test( isfinite(0.0) );

	getInfoStream() << "isfinite( 1.0/0.0 ) = "<< isfinite(1.0/zero);
	_test( !isfinite(1.0/zero) );

	getInfoStream() << "isfinite( -1.0/0.0 ) = "<< isfinite(-1.0/zero);
	_test( !isfinite(-1.0/zero) );

	getInfoStream() << "isfinite( sqrt(-1.0) ) = "	<< isfinite(sqrt(-1.0));
	_test( !isfinite(sqrt(-1.0)) );
}



} //end namespace csmp
