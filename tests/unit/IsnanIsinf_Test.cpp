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
	cout << "isnan(0.0) = "	<< isnan(0.0) << endl;
	_test( isnan(0.0) == 0 );

	cout << "isnan( 1.0/0.0 ) = "<< isnan(1.0/zero) << endl;
	_test( isnan(1.0/zero) == 0 );

	cout << "isnan( -1.0/0.0 ) = "<< isnan(-1.0/zero) << endl;
	_test( isnan(-1.0/zero) == 0 );

	cout << "isnan( sqrt(-1.0) ) = "<< isnan(sqrt(-1.0)) << endl;
	_test( isnan(sqrt(-1.0)) == 1 );

	// isinf tests

	cout << "isinf( 0.0 ) = "<< isinf(0.0) << endl;
	_test( isinf(0.0) == 0 );

	cout << "isinf( 1.0/0.0 ) = "<< isinf(1.0/zero) << endl;
	_test( isinf(1.0/zero) == 1 );

	cout << "isinf( -1.0/0.0 ) = "<< isinf(-1.0/zero) << endl;
	_test( isinf(-1.0/zero) == 1 );

	cout << "isinf( sqrt(-1.0) ) = "<< isinf(sqrt(-1.0)) << endl;
	_test( isinf(sqrt(-1.0)) == 0 );

	// isfinite tests
	cout << "isfinite( 0.0 ) = "<< isfinite(0.0) << endl;
	_test( isfinite(0.0) == 1 );

	cout << "isfinite( 1.0/0.0 ) = "<< isfinite(1.0/zero) << endl;
	_test( isfinite(1.0/zero) == 0 );

	cout << "isfinite( -1.0/0.0 ) = "<< isfinite(-1.0/zero) << endl;
	_test( isfinite(-1.0/zero) == 0 );

	cout << "isfinite( sqrt(-1.0) ) = "	<< isfinite(sqrt(-1.0)) << endl;
	_test( isfinite(sqrt(-1.0)) == 0 );

	
		
}



} //end namespace csmp
