// unit test
// Created by : Siroos
// Date of Modification: 02.11.2010

#include "CubicSpline_Test.h"
#include "vector"
#include <iomanip>
#include <iostream>

using namespace std;

namespace csmp
{
    CubicSpline_Test::CubicSpline_Test()
        : tolerance_( 1.E-2 )
    {
        datafile = "cspline_test_data";
        getInfoStream() << "input data file = " << datafile;

        // x values
        x.push_back(0.0);
        x.push_back(0.5);
        x.push_back(1.0);
        x.push_back(1.5);
        x.push_back(2.0);
        x.push_back(2.5);
        x.push_back(3.0);
        x.push_back(3.5);
        x.push_back(4.0);
        x.push_back(4.5);
        x.push_back(5.0);

        // f(x) = sin(x) values
        fx.push_back(0.0);
        fx.push_back(0.4794);
        fx.push_back(0.8415);
        fx.push_back(0.9975);
        fx.push_back(0.9093);
        fx.push_back(0.5985);
        fx.push_back(0.1411);
        fx.push_back(-0.3508);
        fx.push_back(-0.7568);
        fx.push_back(-0.9775);
        fx.push_back(-0.9589);

        // derivatives of f(x) values f'(x) = cos(x)
        dfx.push_back(1.0);
        dfx.push_back(0.8776);
        dfx.push_back(0.5403);
        dfx.push_back(0.0707);
        dfx.push_back(-0.4161);
        dfx.push_back(-0.8011);
        dfx.push_back(-0.9900);
        dfx.push_back(-0.9365);
        dfx.push_back(-0.6536);
        dfx.push_back(-0.2108);
        dfx.push_back(0.2837);

        df1 = 1.0;
        dfn = 0.2837;
    }

	CubicSpline_Test::~CubicSpline_Test()
    {
		if(cspline2_copy!=NULL)
			delete cspline2_copy;
	}

    void CubicSpline_Test::run()
    {
        InitializeTest( datafile.c_str() );
        InitializeTest( x, fx, df1, dfn);
        OutTest();
        ValueTest();
        DerivativeTest();
        MaxDerivativeTest();
        Range_xTest();
        Range_fxTest();
    }

    void CubicSpline_Test::InitializeTest( const char* file)
    {
        cspline1.Initialize( file );
        cspline1_copy = cspline1;
    }

    void CubicSpline_Test::InitializeTest( const std::vector<double64>& rx, const std::vector<double64>& rfx,
                                           const double64 df1, const double64 dfn )
    {
        cspline2.Initialize( rx, rfx, df1, dfn);
        cspline2_copy = new CubicSpline( cspline2 );
    }

    void CubicSpline_Test::OutTest()
    {
        cspline1.Out(getInfoStream());
    }

    void CubicSpline_Test::ValueTest()          
    {
        getInfoStream() << "\n\n----------------------------------------------" << endl;
        getInfoStream() << "  x\tinput f(x) \t interpolation value" << endl;
        getInfoStream() << "----------------------------------------------" << endl;
        for ( int i = 0; i < x.size(); i++ )
        {
           getInfoStream() << "\n" << setprecision(2) << x[i] << fixed
                << "\t" << setprecision(4) << fx[i]  << fixed
                << "\t\t" << setprecision(20) << cspline1.Value(x[i]) << fixed << endl;
           _equal( cspline1.Value(x[i]), fx[i], tolerance_ );
           _equal( cspline1_copy.Value(x[i]), fx[i], tolerance_ );
           _equal( cspline2.Value(x[i]), fx[i], tolerance_ );
           _equal( cspline2_copy->Value(x[i]), fx[i], tolerance_ );
        }
    }

    void CubicSpline_Test::DerivativeTest()
    {
        getInfoStream() << "\n----------------------------------------------" << endl;
        getInfoStream() << "  x\tinput df(x) \t derevitives" << endl;
        getInfoStream() << "----------------------------------------------" << endl;
        for ( int i = 0; i < x.size(); i++ )
        {
           getInfoStream() << "\n" << setprecision(2) << x[i] << fixed
                << "\t" << setprecision(4) << dfx[i]  << fixed
                << "\t\t" << setprecision(10) << cspline1.Derivative( x[i] ) << fixed << endl;
           _equal( cspline1.Derivative(x[i]), dfx[i], tolerance_ );
           _equal( cspline1_copy.Derivative(x[i]), dfx[i], tolerance_ );
           _equal( cspline2.Derivative(x[i]), dfx[i], tolerance_ );
           _equal( cspline2_copy->Derivative(x[i]), dfx[i], tolerance_ );

        }
    }

    void CubicSpline_Test::MaxDerivativeTest()
    {
        double64 mdv1 = cspline1.MaxDerivative();
        getInfoStream() << "\nMaximum Derivative of input data f'(x)";
        getInfoStream() << "\n" << setprecision(10) << mdv1 << fixed << endl;
        _equal( cspline1.MaxDerivative(), dfx[0], tolerance_ );
        _equal( cspline1_copy.MaxDerivative(), dfx[0], tolerance_ );
        _equal( cspline2.MaxDerivative(), dfx[0], tolerance_ );
        _equal( cspline2_copy->MaxDerivative(), dfx[0], tolerance_ );
    }

    void CubicSpline_Test::Range_xTest()
    {
        double64 x_range = cspline1.Range_x();
        getInfoStream() << "\nRange of x data";
        getInfoStream() << "\n" << setprecision(10) << x_range << fixed << endl;
        _test( cspline1.Range_x() == 5. );
        _test( cspline1_copy.Range_x() == 5. );
        _test( cspline2.Range_x() == 5. );
        _test( cspline2_copy->Range_x() == 5. );
    }

    void CubicSpline_Test::Range_fxTest()
    {
        double64 fx_range = cspline1.Range_fx();
        getInfoStream() << "\nRange of f(x) data";
        getInfoStream() << "\n" << setprecision(20) << fx_range << fixed << endl;
        _equal( cspline1.Range_fx(), ( fx[3]-fx[9] ), tolerance_ );
        _equal( cspline1_copy.Range_fx(), ( fx[3]-fx[9] ), tolerance_ );
        _equal( cspline2.Range_fx(), ( fx[3]-fx[9] ), tolerance_ );
        _equal( cspline2_copy->Range_fx(), ( fx[3]-fx[9] ), tolerance_ );
    }

} // end csmp

