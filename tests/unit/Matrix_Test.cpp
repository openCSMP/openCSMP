// unit test
// Created by : Siroos and Georg
// Date of Modification: 21.11.2010

#include "Matrix_Test.h"

using namespace std;

namespace csmp {

Matrix_Test::Matrix_Test( bool verbose )
 :  A(4,4),
    ResA(4,4),
    B(4,4),
    ResAB(4,4),
    C(4,4),
    ResC(4,4),
    D(4,4),
    ResD(4,4),
    E(4,4),
    F(4,4),
    ResF(4,4),
    G(4,4),
    ResG(4,4),
    H(4,4),
    ResH(4,4),
    Matrix2x2(2,2),
    Matrix3x3(3,3),
   verbose_(verbose)
{
}

Matrix_Test::~Matrix_Test()
{
}

void Matrix_Test::run()
{
   if ( verbose_ ) {
        cout << "\n==============";
        cout << "\nTesting Matrix" << endl;
        cout << "==============" << endl;
        cout << "\nTesting Identity function" << endl;
        cout << "=========================" << endl;
     }

    B( 0, 1 ) = 0.;
    B( 0, 2 ) = 2.;
    B( 0, 3 ) = 1.;

    B( 1, 0 ) = 3.;
    B( 1, 1 ) = 1.;
    B( 1, 2 ) = 0.;
    B( 1, 3 ) = 2.;

    B( 2, 0 ) = 1.;
    B( 2, 1 ) = 0.;
    B( 2, 2 ) = 2.;
    B( 2, 3 ) = 4.;

    B( 3, 0 ) = 0.;
    B( 3, 1 ) = 1.;
    B( 3, 2 ) = 3.;
    B( 3, 3 ) = 5.;

    if ( verbose_ ) cout << "\nB.Identity()";

    B.Identity();
    if ( verbose_ ) B.Out();

    for ( int i = 0; i < 4; i++ ) {
        for ( int j = 0; j < 4; j++ )
        {
            if ( i==j ) _equal( B(i,j), 1., 1E-6 );
            else        _equal( B(i,j), 0., 1E-6 );
        }
    }

//Testing + and - operator
    if ( verbose_ ) {
        cout << "\nTesting + and += operator" << endl;
        cout << "======================" << endl;
      }

    A( 0, 0 ) = 1.;
    A( 0, 1 ) = 0.;
    A( 0, 2 ) = 2.;
    A( 0, 3 ) = 1.;

    A( 1, 0 ) = 3.;
    A( 1, 1 ) = 1.;
    A( 1, 2 ) = 0.;
    A( 1, 3 ) = 2.;

    A( 2, 0 ) = 1.;
    A( 2, 1 ) = 0.;
    A( 2, 2 ) = 2.;
    A( 2, 3 ) = 4.;

    A( 3, 0 ) = 0.;
    A( 3, 1 ) = 1.;
    A( 3, 2 ) = 3.;
    A( 3, 3 ) = 5.;

    if ( verbose_ ) {
    cout << "\nTest Matrix A";
    A.Out();

    cout << "Test Matrix B";
      }
    B.Fill(1.);
    if ( verbose_ ) {
    B.Out();

    cout << "Matrix C = A + B";
      }
    C = A + B;
    C.Out();
    E = C;
  
    if ( verbose_ ) {
    cout << "Matrix A += B =! C";
      }
    A += B;
    if ( verbose_ ) A.Out();

    ResC( 0, 0 ) = 2.;
    ResC( 0, 1 ) = 1.;
    ResC( 0, 2 ) = 3.;
    ResC( 0, 3 ) = 2.;

    ResC( 1, 0 ) = 4.;
    ResC( 1, 1 ) = 2.;
    ResC( 1, 2 ) = 1.;
    ResC( 1, 3 ) = 3.;

    ResC( 2, 0 ) = 2.;
    ResC( 2, 1 ) = 1.;
    ResC( 2, 2 ) = 3.;
    ResC( 2, 3 ) = 5.;

    ResC( 3, 0 ) = 1.;
    ResC( 3, 1 ) = 2.;
    ResC( 3, 2 ) = 4.;
    ResC( 3, 3 ) = 6.;

    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            //cout << C(i, j) << "  ";
            //cout << A(i, j) << endl;
            _equal( A(i,j ), ResC(i,j), 1E-6 );
            _equal( C(i,j ), ResC(i,j), 1E-6 );
        }
    }

    if ( verbose_ ) cout << "Matrix D = C - B";
    D = C - B;
    if ( verbose_ ) D.Out();

    if ( verbose_ ) cout << "Matrix A -= B =! D";
    A -= B;
    if ( verbose_ ) A.Out();

    ResA( 0, 0 ) = 1.;
    ResA( 0, 1 ) = 0.;
    ResA( 0, 2 ) = 2.;
    ResA( 0, 3 ) = 1.;

    ResA( 1, 0 ) = 3.;
    ResA( 1, 1 ) = 1.;
    ResA( 1, 2 ) = 0.;
    ResA( 1, 3 ) = 2.;

    ResA( 2, 0 ) = 1.;
    ResA( 2, 1 ) = 0.;
    ResA( 2, 2 ) = 2.;
    ResA( 2, 3 ) = 4.;

    ResA( 3, 0 ) = 0.;
    ResA( 3, 1 ) = 1.;
    ResA( 3, 2 ) = 3.;
    ResA( 3, 3 ) = 5.;


    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            //cout << D(i, j) << "  ";
            //cout << A(i, j) << endl;
            _equal( D(i,j), ResA(i,j), 1E-6 );
            _equal( A(i,j), ResA(i,j), 1E-6 );
        }
    }


    //Testing * operators
    if ( verbose_ ) {
    cout << "\nTesting * operators" << endl;
    cout << "===================" << endl;

        //Matrix-scalar multiplication
        cout << "\nMatrix-scalar multiplication method" << endl;
        cout << "-----------------------------------" << endl;

        cout << "\nTest matrix A";
        A.Out();

        cout << "Matrix A *= ( 3. )";
    }
    C = A;
    C *= ( 3. );
    if ( verbose_ ) C.Out();

    ResC( 0, 0 ) = 3.;
    ResC( 0, 1 ) = 0.;
    ResC( 0, 2 ) = 6.;
    ResC( 0, 3 ) = 3.;

    ResC( 1, 0 ) = 9.;
    ResC( 1, 1 ) = 3.;
    ResC( 1, 2 ) = 0.;
    ResC( 1, 3 ) = 6.;

    ResC( 2, 0 ) = 3.;
    ResC( 2, 1 ) = 0.;
    ResC( 2, 2 ) = 6.;
    ResC( 2, 3 ) = 12.;

    ResC( 3, 0 ) = 0.;
    ResC( 3, 1 ) = 3.;
    ResC( 3, 2 ) = 9.;
    ResC( 3, 3 ) = 15.;


    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            _equal( C(i,j), ResC(i,j), 1E-6 );
        }
    }


    //Matrix multiplication by vector
    //First vector-matrix multiplication method
    if ( verbose_ ) {
        cout << "\nFirst vector-matrix multiplication method: y = A * x" << endl;
        cout << "----------------------------------------------------" << endl;

        cout << "\nTest matrix A";

        A.Out();
      }
    x.push_back(1.);
    x.push_back(2.);
    x.push_back(0.);
    x.push_back(1.5);

    if ( verbose_ ) cout << "Test vector x =   ";

    vector<double>::const_iterator it( x.begin() );
    for ( it = x.begin(); it != x.end(); it++ )
    {
        if ( verbose_ ) cout << *it << setw(5);
    }

    if ( verbose_ ) cout << "\n\ny = A * x";
    y = A * x;

    vector<double>::const_iterator it1(y.begin());

    cout << "\n\ny = ";

    for ( it1 = y.begin(); it1 != y.end(); it1++ )
    {
        if ( verbose_ ) cout << *it1 << setw(5);
    }

    sol_y.push_back(2.5);
    sol_y.push_back(8.);
    sol_y.push_back(7.);
    sol_y.push_back(9.5);

    vector<double>::const_iterator itsol_y(sol_y.begin());

    for ( itsol_y=sol_y.begin(), it1=y.begin(); itsol_y != sol_y.end(); itsol_y++, it1++ )
        {
            //cout << "\n" << *itsol_y << setw(7) << *it1;
            _equal( *itsol_y, *it1, 1E-6);
        }


    //Second matrix-vector multiplication method
    if ( verbose_ ) {
        cout << "\n\nSecond matrix-vector multiplication method: A *= x" << endl;
        cout << "--------------------------------------------------" << endl;
      }
    B = A;
    if ( verbose_ ) cout << "A *= x";
    B.operator *= (x);

    if ( verbose_ ) cout << "\nA = ";

    for ( int i = 0; i < 4; i++ )
    {
        if ( verbose_ ) cout << B(i,0) << setw(5);

        //Retrieving values of matrix A for comparison
        solB.push_back(B(i,0));
    }

    vector<double>::const_iterator itB( solB.begin() );

    //cout << endl << "\nA  " << setw(5) << "y";

    for ( itB=solB.begin(), itsol_y=sol_y.begin(); itB != solB.end(); itB++, itsol_y++ )
    {
        //cout << "\n" << *itB << setw(7) << *it1;
        _equal( *itB, *itsol_y, 1E-6 );
    }

    if ( verbose_ ) {
    cout << endl;


    cout << "\nThird matrix-vector multiplication method: A *= ( C array )" << endl;
    cout << "--------------------------------------------------" << endl;
     }
  
    D = A;

    if ( verbose_ ) {
        cout << "\nTest matrix A";
        D.Out();
      }

    double Ca[4] = { 1., 2., 0., 1.5 };

    if ( verbose_ ) cout << "Test array Ca";

    for ( int i = 0; i < 4; i++ )
      if ( verbose_ ) cout << "\n" << Ca[i];


    if ( verbose_ ) cout << "\n\nMatrix A *= ( C array )";
    D.operator *=( Ca );
    D.Out();

    for ( int i = 0; i < 4; i++ )
    {
        if ( verbose_ ) cout << D(i,0) << setw(5);

        //Retrieving values of matrix A for comparison
        solD.push_back(D(i,0));
    }


    //cout << endl << "\nA  " << setw(5) << "y";

    vector<double>::const_iterator itD( solD.begin() );
    for ( auto sol = sol_y.begin(); sol != sol_y.end(); sol++, itD++ )
        {
            //cout << "\n" << *itsol_y << setw(7) << *itD;
            _equal( *sol, *itD, 1E-6);
        }

    if ( verbose_ ) cout << endl;


    //Matrix-by-matrix multiplication
    if ( verbose_ ) {
        cout << "\n\nFirst matrix-matrix multiplication method: A * B" << endl;
        cout << "------------------------------------------------" << endl;

        cout << "\nTest matrix A";
        A.Out();

        cout << "Test matrix B";
      }
    B = E;
    if ( verbose_ ) B.Out();

    if ( verbose_ ) cout << "C = A * B";
    C = A * B;
    if ( verbose_ ) C.Out();

    ResAB( 0, 0 ) = 7.;
    ResAB( 0, 1 ) = 5.;
    ResAB( 0, 2 ) = 13.;
    ResAB( 0, 3 ) = 18.;

    ResAB( 1, 0 ) = 12.;
    ResAB( 1, 1 ) = 9.;
    ResAB( 1, 2 ) = 18.;
    ResAB( 1, 3 ) = 21.;

    ResAB( 2, 0 ) = 10.;
    ResAB( 2, 1 ) = 11.;
    ResAB( 2, 2 ) = 25.;
    ResAB( 2, 3 ) = 36.;

    ResAB( 3, 0 ) = 15.;
    ResAB( 3, 1 ) = 15.;
    ResAB( 3, 2 ) = 30.;
    ResAB( 3, 3 ) = 48.;

    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            _equal( C(i,j), ResAB(i,j), 1E-6 );
        }
    }

    if ( verbose_ ) {
        cout << "\nSecond matrix-matrix multiplication method: A *= B" << endl;
        cout << "--------------------------------------------------" << endl;
      }
  
    D = A;

    if ( verbose_ ) {
        cout << "\nTest matrix A";
        D.Out();

        cout << "\nTest matrix B";
        B.Out();

        cout << "\nMatrix A *= B";
      }
//D.operator *=(B);
    D*=B;
    if ( verbose_ ) D.Out();

    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            _equal( D(i,j), ResAB(i,j), 1E-6 );
        }
    }


//Testing fill operators
    if ( verbose_ ) {
        cout << "\nTesting fill functions" << endl;
        cout << "======================" << endl;

        //fill
        cout << "\nMatrix fill()" << endl;
        cout << "-------------" << endl;

        cout << "\nB.fill(1.)";
      }
    B.Fill(1.);
    if ( verbose_ ) B.Out();

    C = 1.;

    for ( int i = 0; i < 4; i++ )
    {
        for ( int j = 0; j < 4; j++ )
        {
            _equal( B(i,j), C(i,j), 1E-6 );
        }
    }

    //fillCol
    if ( verbose_ ) {
        cout << "\nMatrix fillCol()" << endl;
        cout << "---------------" << endl;

        cout << "\nB.fillCol(2,5.5)";
      }
    B.FillCol(2,5.5);
    B.Out();

    for ( int i = 0; i < 4; i++ )
    {
        _equal( B(i,2), 5.5, 1E-6 );
    }

    //fillRow
    if ( verbose_ ) {
        cout << "\nMatrix fillRow()" << endl;
        cout << "----------------" << endl;

        cout << "\nB.fillRow(2,3.3)";
      }
    B.FillRow(2,3.3);
    if ( verbose_ ) B.Out();

    for ( int j = 0; j < 4; j++ )
    {
        _equal( B(2,j), 3.3, 1E-6 );
    }

    C = B;

//Testing summation functions
    if ( verbose_ ) {
        cout << "\nTesting summation function" << endl;
        cout << "==========================" << endl;

        //ColSum()
        cout << "\nMatrix ColSum()" << endl;
        cout << "---------------" << endl;

        cout << "\nB.ColSum(2) = " << B.ColSum(2) << endl;
      }
  
    double ColSumB2 = B.ColSum(2);
    double CheckColSumB2 = 0.;

    for ( int i = 0; i < 4; i++ )
    {
        CheckColSumB2 += B(i,2);
    }

    if ( verbose_ ) cout << "\nB.CheckColSum(2) = " << CheckColSumB2 << endl;
    _equal( CheckColSumB2, ColSumB2, 1E-6 );


    //RowSum
    if ( verbose_ ) {
        cout << "\nMatrix RowSum()" << endl;
        cout << "---------------" << endl;

        cout << "\nB.RowSum(2) = " << B.RowSum(2) << endl;
      }
    double RowSumB2 = B.RowSum(2);
    double CheckRowSumB2 = 0.;

    for ( int j = 0; j < 4; j++ )
    {
        CheckRowSumB2 += B(2,j);
    }

    if ( verbose_ ) cout << "\nB.CheckRowSum(2) = " << CheckRowSumB2 << endl;
    _equal( CheckRowSumB2, RowSumB2, 1E-6 );

    //Testing Zero functions
    if ( verbose_ ) {
        cout << "\nTesting Zero functions" << endl;
        cout << "======================" << endl;

        //ZeroCol
        cout << "\nMatrix ZeroCol()" << endl;
        cout << "----------------" << endl;

        cout << "\nB.ZeroCol(2)";
      }
        B.ZeroCol(2);
        if ( verbose_ ) B.Out();

        for ( int i = 0; i < 4; i++ )
        {
            _equal( B(i,2), 0., 1E-6 );
        }

        //ZeroRow
        if ( verbose_ ) {
            cout << "\nMatrix ZeroRow()" << endl;
            cout << "----------------" << endl;

            cout << "\nB.ZeroRow(2)";
          }
        B.ZeroRow(2);
        if ( verbose_ ) B.Out();

        for ( int j = 0; j < 4; j++ )
        {
            _equal( B(2,j), 0., 1E-6 );
        }

        //Zero
        if ( verbose_ ) {
            cout << "\nMatrix Zero()" << endl;
            cout << "-------------" << endl;

            cout << "\nB.Zero()";
          }
        B.Zero();
        if ( verbose_ ) B.Out();

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(i,j), 0., 1E-6 );
            }
        }


    //Testing operator= ()
    if ( verbose_ ) {
        cout << "\nTesting operator= ()" << endl;
        cout << "====================" << endl;

        //Testing operator= (const Matrix& )
        cout << "\nTesting operator= ( const Matrix& )" << endl;
        cout << "---------------------------------------" << endl;

        cout << "\nMatrix B";
        B.Out();

        cout << "\nMatrix B = C";
      }
        B = C;
        if ( verbose_ ) B.Out();

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(i,j), C(i,j), 1E-6 );
            }
        }


    //Testing transposed functions
    if ( verbose_ ) {
        cout << "\nTesting transposed matrix functions" << endl;
        cout << "===================================" << endl;

        //Transposed matrix
        cout << "\nMatrix Transposed()" << endl;
        cout << "-------------------" << endl;

        cout << "\nTest matrix A";
        A.Out();
      }
        A.Transposed(F);
      if ( verbose_ ) {
          cout << "\nF = A^T = A.Transposed(F)";
          F.Out();
        }
  
        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( F(i,j), A(j,i), 1E-6 );
            }
        }

        //Transposed product
        if ( verbose_ ) {
            cout << "\nMatrix TransposedProduct()" << endl;
            cout << "--------------------------" << endl;

            cout << "\nTest matrix A";
            A.Out();
          }
        A.TransposedProduct(H);
  
        if ( verbose_ ) {
            cout << "\nH = A^T * A = A.TransposedProduct(F)";
            H.Out();
          }

        ResH( 0, 0 ) = 11.;
        ResH( 0, 1 ) = 3.;
        ResH( 0, 2 ) = 4;
        ResH( 0, 3 ) = 11.;

        ResH( 1, 0 ) = 3.;
        ResH( 1, 1 ) = 2.;
        ResH( 1, 2 ) = 3.;
        ResH( 1, 3 ) = 7.;

        ResH( 2, 0 ) = 4.;
        ResH( 2, 1 ) = 3.;
        ResH( 2, 2 ) = 17.;
        ResH( 2, 3 ) = 25.;

        ResH( 3, 0 ) = 11.;
        ResH( 3, 1 ) = 7.;
        ResH( 3, 2 ) = 25.;
        ResH( 3, 3 ) = 46.;

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( H(i,j), ResH(j,i), 1E-6 );
            }
        }


        //MultiplyWithTransposedOf
        if ( verbose_ ) {
            cout << "\nMatrix MultiplyWithTransposedOf()" << endl;
            cout << "---------------------------------" << endl;

            cout << "\nTest matrix A";
            A.Out();

            cout << "\nTest matrix B";
            B.Out();

            cout << "\nF = A B^T = A.MultiplyWithTransposedOf( B, F )";
          }
        A.MultiplyWithTransposedOf( B, F );
        if ( verbose_ ) F.Out();

        ResF( 0, 0 ) = 13.;
        ResF( 0, 1 ) = 13.;
        ResF( 0, 2 ) = 13.2;
        ResF( 0, 3 ) = 13.;

        ResF( 1, 0 ) = 6.;
        ResF( 1, 1 ) = 6.;
        ResF( 1, 2 ) = 19.8;
        ResF( 1, 3 ) = 6.;

        ResF( 2, 0 ) = 16.;
        ResF( 2, 1 ) = 16.;
        ResF( 2, 2 ) = 23.1;
        ResF( 2, 3 ) = 16.;

        ResF( 3, 0 ) = 22.5;
        ResF( 3, 1 ) = 22.5;
        ResF( 3, 2 ) = 29.7;
        ResF( 3, 3 ) = 22.5;

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( ResF(i,j), F(i,j), 1E-6 );
            }
        }

        //MultiplyTransposedOfWith
        if ( verbose_ ) {
            cout << "\nMatrix MultiplyTransposedOfWith()" << endl;
            cout << "---------------------------------" << endl;

            cout << "\nTest matrix A";
            A.Out();

            cout << "\nTest matrix B";
            B.Out();

            cout << "\nG = A^T B = A.MultiplyTransposedOfWith( B, G )";
          }
        A.MultiplyTransposedOfWith( B, G );
        if ( verbose_ ) G.Out();

        ResG( 0, 0 ) = 7.3;
        ResG( 0, 1 ) = 7.3;
        ResG( 0, 2 ) = 25.3;
        ResG( 0, 3 ) = 7.3;

        ResG( 1, 0 ) = 2.;
        ResG( 1, 1 ) = 2.;
        ResG( 1, 2 ) = 11.;
        ResG( 1, 3 ) = 2.;

        ResG( 2, 0 ) = 11.6;
        ResG( 2, 1 ) = 11.6;
        ResG( 2, 2 ) = 34.1;
        ResG( 2, 3 ) = 11.6;

        ResG( 3, 0 ) = 21.2;
        ResG( 3, 1 ) = 21.2;
        ResG( 3, 2 ) = 57.2;
        ResG( 3, 3 ) = 21.2;

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( ResG(i,j), G(i,j), 1E-6 );
            }
        }


    //Testing operator= ()
    if ( verbose_ ) {
        cout << "\nTesting operator= ()" << endl;
        cout << "====================================" << endl;

        //Operator= ( TensorVariable 2U )
        cout << "\nOperator= ( TensorVariable 2U )" << endl;
        cout << "-------------------------------" << endl;

        cout << "\n2U test matrix Matrix2x2";
        Matrix2x2.Out();
      }
  
        double val = 16.;

        if ( verbose_ ) cout << "\nMatrix2x2 = val";
        Matrix2x2 = val;
        if ( verbose_ ) Matrix2x2.Out();

        for ( int i = 0; i < 2; i++ )
        {
            for ( int j = 0; j < 2; j++ )
            {
                _equal( Matrix2x2(i,j), 16., 1E-6 );
            }
        }


        //Operator= ( TensorVariable 2U )
        if ( verbose_ ) {
            cout << "\nOperator= ( TensorVariable 2U )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n2U test matrix Matrix2x2";
          }
        Matrix2x2.Fill(5.);
        if ( verbose_ ) {
            Matrix2x2.Out();

            cout << "\n2U tensor variable TensorVariable2U";
          }
        TensorVariable2U.Identity();
        TensorVariable2U += 1.;
        if ( verbose_ ) TensorVariable2U.Out();

        if ( verbose_ ) cout << "\nMatrix2x2 = TensorVariable2U";
        Matrix2x2 = TensorVariable2U;
        if ( verbose_ ) Matrix2x2.Out();

        for ( int i = 0; i < 2; i++ )
        {
            for ( int j = 0; j < 2; j++ )
            {
                if ( i==j ) _equal( Matrix2x2(i,j), 2., 1E-6 );
                else        _equal( Matrix2x2(i,j), 1., 1E-6 );
            }
        }


        //Operator= ( TensorVariable 3U )
        if ( verbose_ ) {
            cout << "\nOperator= ( TensorVariable 3U )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n3U test matrix Matrix3x3";
          }
        Matrix3x3.Fill(5.);
        if ( verbose_ ) Matrix3x3.Out();

        if ( verbose_ ) cout << "\n3U tensor variable TensorVariable3U";
        TensorVariable3U.Identity();
        TensorVariable3U += 2.;
        if ( verbose_ ) TensorVariable3U.Out();

        if ( verbose_ ) cout << "\nMatrix3x3 = TensorVariable3U";
        Matrix3x3 = TensorVariable3U;
        if ( verbose_ ) Matrix3x3.Out();

        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                if ( i==j ) _equal( Matrix3x3(i,j), 3., 1E-6 );
                else        _equal( Matrix3x3(i,j), 2., 1E-6 );
            }
        }

        //Testing operator*= ( TensorVariable )
        if ( verbose_ ) {
            cout << "\nTesting operator*= ( TensorVariable )" << endl;
            cout << "====================================" << endl;

            //Operator *= ( double )
            cout << "\nOperator *= ( double )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n3U test matrix Matrix3x3";
          }
            Matrix3x3.Resize(3,3);
            Matrix3x3.Fill(5.);
            if ( verbose_ ) Matrix3x3.Out();

            if ( verbose_ ) cout << "\nMatrix3x3 *= ( 5 )";
            Matrix3x3 *= ( 5 );
            if ( verbose_ ) Matrix3x3.Out();

            _equal( Matrix3x3(0,0), 25., 1E-6 );

            //Operator *= ( TensorVariable 2U )
            if ( verbose_ ) {
                cout << "\nOperator *= ( TensorVariable 2U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix Matrix3x3";
              }
            Matrix3x3.Resize(2,2);
            Matrix3x3.Fill(5.);
            if ( verbose_ ) Matrix3x3.Out();

            if ( verbose_ ) cout << "\n2U tensor variable TensorVariable2U";
            TensorVariable2U.Identity();
            TensorVariable2U += 1.;
            TensorVariable2U.Out();

            if ( verbose_ ) cout << "\nMatrix3x3 *= ( TensorVariable2U )";
            Matrix3x3 *= ( TensorVariable2U );
            Matrix3x3.Out();

            for ( int i = 0; i < 2; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    _equal( Matrix3x3(i,j), 15., 1E-6 );
                }
            }

            //Operator*= ( TensorVariable 3U )
            if ( verbose_ ) {
                cout << "\nOperator*= ( TensorVariable 3U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix Matrix3x3";
              }
            Matrix3x3.Resize(3,3);
            Matrix3x3.Fill(5.);
            if ( verbose_ ) Matrix3x3.Out();

            if ( verbose_ ) cout << "\n3U tensor variable TensorVariable3U";
            TensorVariable3U.Identity();
            TensorVariable3U += 2.;
            if ( verbose_ ) TensorVariable3U.Out();

            if ( verbose_ ) cout << "\nMatrix3x3 *= ( TensorVariable3U )";
            Matrix3x3 *= ( TensorVariable3U );
            if ( verbose_ ) Matrix3x3.Out();

            for ( int i = 0; i < 3; i++ )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    _equal( Matrix3x3(i,j), 35., 1E-6 );
                }
            }


        //Testing operator*= ( VectorVariable )
        if ( verbose_ ) {
            cout << "\nTesting operator*= ( VectorVariable )" << endl;
            cout << "====================================" << endl;

            //Operator *= ( VectorVariable 2U )
            cout << "\nOperator *= ( VectorVariable 2U )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n2U test matrix Matrix2x2";
          }
            Matrix3x3.Resize(2,2);
            Matrix3x3.Fill(5.);
            if ( verbose_ ) Matrix3x3.Out();

            if ( verbose_ ) cout << "\n2U Vector variable VectorVariable2U";
            VectorVariable2U = 5;
            VectorVariable2U += 1.;
            if ( verbose_ ) VectorVariable2U.Out();

            if ( verbose_ ) cout << "\nMatrix2x2 *= ( VectorVariable2U )";
            Matrix3x3 *= ( VectorVariable2U );
            if ( verbose_ ) Matrix3x3.Out();

            for ( int i = 0; i < 2; i++ )
            {
                _equal( Matrix3x3(i,0), 60., 1E-6 );

            }


            //Operator*= ( VectorVariable 3U )
            if ( verbose_ ) {
                cout << "\nOperator*= ( VectorVariable 3U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix Matrix3x3";
              }
            Matrix3x3.Resize(3,3);
            Matrix3x3.Fill(5.);
            if ( verbose_ ) Matrix3x3.Out();

            if ( verbose_ ) cout << "\n3U Vector variable VectorVariable3U";
            VectorVariable3U = 5;
            VectorVariable3U += 2.;
            if ( verbose_ ) VectorVariable3U.Out();

            if ( verbose_ ) cout << "\nMatrix3x3 *= ( VectorVariable3U )";
            Matrix3x3 *= ( VectorVariable3U );
            if ( verbose_ ) Matrix3x3.Out();

            for ( int i = 0; i < 3; i++ )
            {
                _equal( Matrix3x3(i,0), 105., 1E-6 );
            }


        //RowCondenseTo( vector<double> )
        if ( verbose_ ) {
            cout << "\nTesting RowCondenseTo( vector<double> )" << endl;
            cout << "===========================================" << endl;
            //RowCondenseTo ( vector<double> )
            cout << "\nRowCondenseTo ( vector<double> )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n3U test matrix Matrix3x3";
        }
            Matrix3x3.Resize(3,3);
            Matrix3x3.Fill(9.);
            if ( verbose_ ) Matrix3x3.Out();

            vector<double>::const_iterator itx(x.begin());

            if ( verbose_ ) cout << "Test vector x";

            for ( itx = x.begin(); itx != x.end(); itx++ )
            {
                if ( verbose_ ) cout << "\n" << *itx;
            }

            if ( verbose_ ) cout << "\n\nMatrix3x3.RowCondenseTo( x )";
            Matrix3x3.RowCondenseTo( x );

            for ( itx = x.begin(); itx != x.end(); itx++)
            {
                if ( verbose_ ) cout << "\n" << *itx ;
                _equal( *itx, 27., 1E-6 );
            }

            if ( verbose_ ) cout << endl;


        //ExportTo ( TensorVariable )
        if ( verbose_ ) {
            cout << "\nTesting ExportTo ( TensorVariable )" << endl;
            cout << "===========================================" << endl;

            //ExportTo ( TensorVariable 2U )
            cout << "\nExportTo ( TensorVariable 2U )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n2U test matrix Matrix2x2";
          }
            Matrix3x3.Resize(2,2);
            Matrix3x3.Identity();
            if ( verbose_ ) Matrix3x3.Out();

            if ( verbose_ ) cout << "\n2U tensor variable TensorVariable2U";
            TensorVariable2U=0.;
            if ( verbose_ ) TensorVariable2U.Out();

            if ( verbose_ ) cout << "\nMatrix2x2 *= ( TensorVariable2U )";
            Matrix3x3.ExportTo(TensorVariable2U);
            if ( verbose_ ) TensorVariable3U.Out();

            for ( int i = 0; i < 2; i++ )
            {
                _equal( Matrix3x3(i,i), 1., 1E-6 );
            }


            //ExportTo ( TensorVariable 3U )
            if ( verbose_ ) {
                cout << "\n\nExportTo ( TensorVariable 3U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix Matrix3x3";
              }
            Matrix3x3.Resize(3,3);
            Matrix3x3.Identity();
            if ( verbose_ ) Matrix3x3.Out();

            if ( verbose_ ) cout << "\n3U tensor variable TensorVariable3U";
            TensorVariable3U=0.;
            if ( verbose_ ) TensorVariable3U.Out();

            if ( verbose_ ) cout << "\nMatrix3x3 *= ( TensorVariable3U )";
            Matrix3x3.ExportTo(TensorVariable3U);
            if ( verbose_ ) TensorVariable3U.Out();

            for ( int i = 0; i < 3; i++ )
            {
                _equal( Matrix3x3(i,i), 1., 1E-6 );
            }

            //Testing L1 matrix norm()
            if ( verbose_ ) {
                cout << "\nTesting NormL1()" << endl;
                cout << "================" << endl;

                    cout << "\n3U test matrix A";
              }
                    A(2,2) = -99.;
                    A(2,1) = -199.;
                    if ( verbose_ ) A.Out();
                    double normL1 = A.NormL1();

                    if ( verbose_ ) cout << "\nThe maximum absolute column sum norm of test matrix A is " << normL1 << endl;
                    _equal( normL1, 201., 1E-6 );


            //Testing L infinity matrix norm()
            if ( verbose_ ) {
                cout << "\nTesting NormL_Infinity()" << endl;
                cout << "================" << endl;

                    cout << "\n3U test matrix A";
              }
                    A(2,2) = -99.;
                    A(2,1) = -199.;
                    if ( verbose_ ) A.Out();
                    double normInf = A.NormL_Infinity();

                    if ( verbose_ ) cout << "\nThe maximum absolute row sum norm of test matrix A is " << normInf << endl;
                    _equal( normInf, 303., 1E-6 );


            //Testing In() has been done manually and confirmed
            //  Otherwise script based In() call must be answered



}


} // end csmp

