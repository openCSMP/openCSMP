// unit test
// Created by : Siroos and Georg
// Date of Modification: 21.11.2010

#include "DenseMatrix_Test.h"

using namespace std;

namespace csmp
{
    DenseMatrix_Test::DenseMatrix_Test()
    {       
    }

    DenseMatrix_Test::~DenseMatrix_Test()
    {
    }

    void DenseMatrix_Test::run()
    {

    cout << "\n===================";
    cout << "\nTesting DenseMatrix" << endl;
    cout << "===================" << endl;

    //Testing Identity functions
        cout << "\nTesting Identity function" << endl;
        cout << "=========================" << endl;

        cout << "\nB.Identity()";
        B.Identity();
        B.Out();

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                if ( i==j ) _equal( B(i,j), 1., 1E-6 );
                else        _equal( B(i,j), 0., 1E-6 );
            }
        }

    //Testing + and - operator

        cout << "\nTesting + and += operator" << endl;
        cout << "======================" << endl;

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

        cout << "\nTest Matrix A";
        A.Out();

        cout << "Test Matrix B";
        B.Fill(1.);
        B.Out();

        cout << "Matrix C = A + B";
        C = A + B;
        C.Out();
        E = C;

        cout << "Matrix A += B =! C";
        A += B;
        A.Out();

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

        cout << "Matrix D = C - B";
        D = C - B;
        D.Out();

        cout << "Matrix A -= B =! D";
        A -= B;
        A.Out();

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
    cout << "\nTesting * operators" << endl;
    cout << "===================" << endl;

        //Matrix-scalar multiplication
        cout << "\nMatrix-scalar multiplication method" << endl;
        cout << "-----------------------------------" << endl;

        cout << "\nTest matrix A";
        A.Out();

        cout << "Matrix A *= ( 3 )";
        C = A;
        C.operator *= ( 3. );
        C.Out();

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


        //Matrix-scalar multiplication
        cout << "\nMatrix-scalar multiplication method" << endl;
        cout << "-----------------------------------" << endl;

        cout << "\nTest matrix A";
        A.Out();

        cout << "Matrix A *= ( 3 )";
        C = A;
        C.operator *= ( 3. );
        C.Out();


        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( C(i,j), ResC(i,j), 1E-6 );
            }
        }


        //Matrix multiplication by vector
        //First vector-matrix multiplication method
        cout << "\nFirst vector-matrix multiplication method: y = A * x" << endl;
        cout << "----------------------------------------------------" << endl;

        cout << "\nTest matrix A";

        A.Out();

        x.push_back(1.);
        x.push_back(2.);
        x.push_back(0.);
        x.push_back(1.5);

        cout << "Test vector x =   ";

        vector<double64>::const_iterator it( x.begin() );
        for ( vector<double64>::const_iterator
             it = x.begin(); it != x.end(); it++ )
        {
            cout << *it << setw(5);
        }

        cout << "\n\ny = A * x";
        y = A * x;

        vector<double64>::const_iterator it1(y.begin());

        cout << "\n\ny = ";

        for (vector<double64>::const_iterator
             it1 = y.begin(); it1 != y.end(); it1++ )
        {
            cout << *it1 << setw(5);
        }

        sol_y.push_back(2.5);
        sol_y.push_back(8.);
        sol_y.push_back(7.);
        sol_y.push_back(9.5);

        vector<double64>::const_iterator itsol_y(sol_y.begin());

        for (vector<double64>::const_iterator
             itsol_y = sol_y.begin(); itsol_y != sol_y.end(); itsol_y++, it1++ )
            {
                //cout << "\n" << *itsol_y << setw(7) << *it1;
                _equal( *itsol_y, *it1, 1E-6);
            }


        //\nSecond matrix-vector multiplication method
        cout << "\n\nSecond matrix-vector multiplication method: A *= x" << endl;
        cout << "--------------------------------------------------" << endl;

        B = A;
        cout << "A *= x";
        B.operator *= (x);

        cout << "\nA = ";

        for ( int i = 0; i < 4; i++ )
        {
            cout << B(i,0) << setw(5);

            //Retrieving values of matrix A for comparison
            solB.push_back(B(i,0));
        }

        vector<double64>::const_iterator itB( solB.begin() );

        //cout << endl << "\nA  " << setw(5) << "y";

        for (vector<double64>::const_iterator
             itB = solB.begin(); itB != solB.end(); itB++, itsol_y++ )
        {
            //cout << "\n" << *itB << setw(7) << *it1;
            _equal( *itB, *itsol_y, 1E-6 );
        }

        cout << endl;


        cout << "\nThird matrix-matrix multiplication method: A *= ( C array )" << endl;
        cout << "--------------------------------------------------" << endl;

        D = A;

        cout << "\nTest matrix A";
        D.Out();


        double64 Ca[4] = { 1., 2., 0., 1.5 };

        cout << "Test array Ca";

        for ( int i = 0; i < 4; i++ )
            cout << "\n" << Ca[i];


        cout << "\n\nMatrix A *= ( C array )";
        D.operator *=( Ca );
        D.Out();

        for ( int i = 0; i < 4; i++ )
        {
            cout << D(i,0) << setw(5);

            //Retrieving values of matrix A for comparison
            solD.push_back(D(i,0));
        }

        vector<double64>::const_iterator itD( solD.begin() );

        //cout << endl << "\nA  " << setw(5) << "y";

        for (vector<double64>::const_iterator
             itsol_y = sol_y.begin(); itsol_y != sol_y.end(); itsol_y++, itD++ )
            {
                //cout << "\n" << *itsol_y << setw(7) << *itD;
                _equal( *itsol_y, *itD, 1E-6);
            }

        cout << endl;


        //Matrix-by-matrix multiplication
        cout << "\n\nFirst matrix-matrix multiplication method: A * B" << endl;
        cout << "------------------------------------------------" << endl;

        cout << "\nTest matrix A";
        A.Out();

        cout << "Test matrix B";
        B = E;
        B.Out();

        cout << "C = A * B";
        C = A * B;
        C.Out();

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

        cout << "\nSecond matrix-matrix multiplication method: A *= B" << endl;
        cout << "--------------------------------------------------" << endl;

        D = A;

        cout << "\nTest matrix A";
        D.Out();

        cout << "\nMatrix A *= B";
        D.operator *=(B);
        D.Out();

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( D(i,j), ResAB(i,j), 1E-6 );
            }
        }


    //Testing fill operators
        cout << "\nTesting fill functions" << endl;
        cout << "======================" << endl;

        //fill
        cout << "\nMatrix fill()" << endl;
        cout << "-------------" << endl;

        cout << "\nB.fill(1.)";
        B.Fill(1.);
        B.Out();

        C = 1.;

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(i,j), C(i,j), 1E-6 );
            }
        }

        //fillCol
        cout << "\nMatrix fillCol()" << endl;
        cout << "---------------" << endl;

        cout << "\nB.fillCol(2,5.5)";
        B.FillCol(2,5.5);
        B.Out();

        for ( int i = 0; i < 4; i++ )
        {
            _equal( B(i,2), 5.5, 1E-6 );
        }

        //fillRow
        cout << "\nMatrix fillRow()" << endl;
        cout << "----------------" << endl;

        cout << "\nB.fillRow(2,3.3)";
        B.FillRow(2,3.3);
        B.Out();

        for ( int j = 0; j < 4; j++ )
        {
            _equal( B(2,j), 3.3, 1E-6 );
        }

        C = B;

    //Testing summation functions
        cout << "\nTesting summation function" << endl;
        cout << "==========================" << endl;

        //ColSum()
        cout << "\nMatrix ColSum()" << endl;
        cout << "---------------" << endl;

        cout << "\nB.ColSum(2) = " << B.ColSum(2) << endl;

        double64 ColSumB2 = B.ColSum(2);
        double64 CheckColSumB2 = 0.;

        for ( int i = 0; i < 4; i++ )
        {
            CheckColSumB2 += B(i,2);
        }

        cout << "\nB.CheckColSum(2) = " << CheckColSumB2 << endl;
        _equal( CheckColSumB2, ColSumB2, 1E-6 );


        //RowSum
        cout << "\nMatrix RowSum()" << endl;
        cout << "---------------" << endl;

        cout << "\nB.RowSum(2) = " << B.RowSum(2) << endl;

        double64 RowSumB2 = B.RowSum(2);
        double64 CheckRowSumB2 = 0.;

        for ( int j = 0; j < 4; j++ )
        {
            CheckRowSumB2 += B(2,j);
        }

        cout << "\nB.CheckRowSum(2) = " << CheckRowSumB2 << endl;
        _equal( CheckRowSumB2, RowSumB2, 1E-6 );

        //Testing Zero functions
            cout << "\nTesting Zero functions" << endl;
            cout << "======================" << endl;

            //ZeroCol
            cout << "\nMatrix ZeroCol()" << endl;
            cout << "----------------" << endl;

            cout << "\nB.ZeroCol(2)";
            B.ZeroCol(2);
            B.Out();

            for ( int i = 0; i < 4; i++ )
            {
                _equal( B(i,2), 0., 1E-6 );
            }

            //ZeroRow
            cout << "\nMatrix ZeroRow()" << endl;
            cout << "----------------" << endl;

            cout << "\nB.ZeroRow(2)";
            B.ZeroRow(2);
            B.Out();

            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(2,j), 0., 1E-6 );
            }

            //Zero
            cout << "\nMatrix Zero()" << endl;
            cout << "-------------" << endl;

            cout << "\nB.Zero()";
            B.Zero();
            B.Out();

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( B(i,j), 0., 1E-6 );
                }
            }


        //Testing operator= ()
            cout << "\nTesting operator= ()" << endl;
            cout << "====================" << endl;

            //Testing operator= (const DenseMatrix& )
            cout << "\nTesting operator= ( const DenseMatrix& )" << endl;
            cout << "---------------------------------------" << endl;

            cout << "\nMatrix B";
            B.Out();

            cout << "\nMatrix B = C";
            B = C;
            B.Out();

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( B(i,j), C(i,j), 1E-6 );
                }
            }


        //Testing transposed functions
            cout << "\nTesting transposed matrix functions" << endl;
            cout << "===================================" << endl;

            //Transposed matrix
            cout << "\nMatrix Transposed()" << endl;
            cout << "-------------------" << endl;

            cout << "\nTest matrix A";
            A.Out();

            A.Transposed(F);
            cout << "\nF = A^T = A.Transposed(F)";
            F.Out();

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( F(i,j), A(j,i), 1E-6 );
                }
            }

            //MultiplyWithTransposedOf
            cout << "\nMatrix MultiplyWithTransposedOf()" << endl;
            cout << "---------------------------------" << endl;

            cout << "\nTest matrix A";
            A.Out();


            cout << "\nF = A B^T = A.MultiplyWithTransposedOf( B, F )";
            A.MultiplyWithTransposedOf( B, F );
            F.Out();

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
            cout << "\nMatrix MultiplyTransposedOfWith()" << endl;
            cout << "---------------------------------" << endl;

            cout << "\nTest matrix A";
            A.Out();

            cout << "\nTest matrix B";
            B.Out();

            cout << "\nG = A^T B = A.MultiplyTransposedOfWith( B, G )";
            A.MultiplyTransposedOfWith( B, G );
            G.Out();

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


        //Testing operator= ( TensorVariable )
            cout << "\nTesting operator= ( TensorVariable )" << endl;
            cout << "====================================" << endl;

            //Operator= ( TensorVariable 1U )
            cout << "\nOperator= ( TensorVariable 1U )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n3U test matrix DenseMatrix3x3";
            DenseMatrix3x3.Fill(5.);
            DenseMatrix3x3.Out();

            cout << "\n1U tensor variable TensorVariable1U";
            TensorVariable1U.Identity();
            TensorVariable1U.Out();

            cout << "\nDenseMatrix3x3 = TensorVariable1U";
            DenseMatrix3x3 = TensorVariable1U;
            DenseMatrix3x3.Out();

            _equal( DenseMatrix3x3(0,0), 1., 1E-6 );

            //Operator= ( TensorVariable 2U )
            cout << "\nOperator= ( TensorVariable 2U )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n3U test matrix DenseMatrix3x3";
            DenseMatrix3x3.Fill(5.);
            DenseMatrix3x3.Out();

            cout << "\n2U tensor variable TensorVariable2U";
            TensorVariable2U.Identity();
            TensorVariable2U += 1.;
            TensorVariable2U.Out();

            cout << "\nDenseMatrix3x3 = TensorVariable2U";
            DenseMatrix3x3 = TensorVariable2U;
            DenseMatrix3x3.Out();

            for ( int i = 0; i < 2; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    if ( i==j ) _equal( DenseMatrix3x3(i,j), 2., 1E-6 );
                    else        _equal( DenseMatrix3x3(i,j), 1., 1E-6 );
                }
            }


            //Operator= ( TensorVariable 3U )
            cout << "\nOperator= ( TensorVariable 3U )" << endl;
            cout << "-------------------------------" << endl;

            cout << "\n3U test matrix DenseMatrix3x3";
            DenseMatrix3x3.Fill(5.);
            DenseMatrix3x3.Out();

            cout << "\n3U tensor variable TensorVariable3U";
            TensorVariable3U.Identity();
            TensorVariable3U += 2.;
            TensorVariable3U.Out();

            cout << "\nDenseMatrix3x3 = TensorVariable3U";
            DenseMatrix3x3 = TensorVariable3U;
            DenseMatrix3x3.Out();

            for ( int i = 0; i < 3; i++ )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    if ( i==j ) _equal( DenseMatrix3x3(i,j), 3., 1E-6 );
                    else        _equal( DenseMatrix3x3(i,j), 2., 1E-6 );
                }
            }

            //Testing operator*= ( TensorVariable )
                cout << "\nTesting operator*= ( TensorVariable )" << endl;
                cout << "====================================" << endl;

                //Operator *= ( double64 )
                cout << "\nOperator *= ( double64 )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\nDenseMatrix3x3 *= ( 5 )";
                DenseMatrix3x3 *= ( 5 );
                DenseMatrix3x3.Out();

                _equal( DenseMatrix3x3(0,0), 25., 1E-6 );

                //Operator *= ( TensorVariable 2U )
                cout << "\nOperator *= ( TensorVariable 2U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\n2U tensor variable TensorVariable2U";
                TensorVariable2U.Identity();
                TensorVariable2U += 1.;
                TensorVariable2U.Out();

                cout << "\nDenseMatrix3x3 *= ( TensorVariable2U )";
                DenseMatrix3x3 *= ( TensorVariable2U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 2; i++ )
                {
                    for ( int j = 0; j < 2; j++ )
                    {
                        _equal( DenseMatrix3x3(i,j), 15., 1E-6 );
                    }
                }


                //Operator*= ( TensorVariable 3U )
                cout << "\nOperator*= ( TensorVariable 3U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\n3U tensor variable TensorVariable3U";
                TensorVariable3U.Identity();
                TensorVariable3U += 2.;
                TensorVariable3U.Out();

                cout << "\nDenseMatrix3x3 *= ( TensorVariable3U )";
                DenseMatrix3x3 *= ( TensorVariable3U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 3; i++ )
                {
                    for ( int j = 0; j < 3; j++ )
                    {
                        _equal( DenseMatrix3x3(i,j), 35., 1E-6 );
                    }
                }

            //Testing operator*= ( VectorVariable )
                cout << "\nTesting operator*= ( VectorVariable )" << endl;
                cout << "====================================" << endl;

                //Operator *= ( VectorVariable 2U )
                cout << "\nOperator *= ( VectorVariable 2U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\n2U Vector variable VectorVariable2U";
                VectorVariable2U = 5;
                VectorVariable2U += 1.;
                VectorVariable2U.Out();

                cout << "\nDenseMatrix2x2 *= ( VectorVariable2U )";
                DenseMatrix3x3 *= ( VectorVariable2U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(i,0), 60., 1E-6 );

                }


                //Operator*= ( VectorVariable 3U )
                cout << "\nOperator*= ( VectorVariable 3U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\n3U Vector variable VectorVariable3U";
                VectorVariable3U = 5;
                VectorVariable3U += 2.;
                VectorVariable3U.Out();

                cout << "\nDenseMatrix3x3 *= ( VectorVariable3U )";
                DenseMatrix3x3 *= ( VectorVariable3U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(i,0), 105., 1E-6 );
                }


            //AssignToDiagonal ( TensorVariable )
                cout << "\nTesting AssignToDiagonal ( VectorVariable )" << endl;
                cout << "===========================================" << endl;

                //AssignToDiagonal ( diag_elmts, ScalarVariable )
                cout << "\nAssignToDiagonal ( diag_elmts, ScalarVariable )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\nDenseMatrix3x3.AssignToDiagonal ( 2, 6. )";
                DenseMatrix3x3.AssignToDiagonal(2, 6.);
                DenseMatrix3x3.Out();

                for ( int i = 0, j = 0; i < 2; i++, j++ )
                {
                   _equal( DenseMatrix3x3(i,j), 6., 1E-6 );
                }

                //AssignToDiagonal ( VectorVariable1U )
                cout << "\nAssignToDiagonal ( VectorVariable 1U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n1U test matrix DenseMatrix1x1";
                DenseMatrix3x3.Resize(1,1);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\n1U Vector variable VectorVariable1U";
                VectorVariable1U = 5;
                VectorVariable1U += 6.;
                VectorVariable1U.Out();

                cout << "\nDenseMatrix1x1.AssignToDiagonal ( VectorVariable1U )";
                DenseMatrix3x3.AssignToDiagonal ( VectorVariable1U );
                DenseMatrix3x3.Out();

                _equal( DenseMatrix3x3(0,0), 11., 1E-6 );

                //AssignToDiagonal ( VectorVariable2U )
                cout << "\nAssignToDiagonal ( VectorVariable 2U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\n2U Vector variable VectorVariable2U";
                VectorVariable2U = 9;
                VectorVariable2U += 1.;
                VectorVariable2U.Out();

                cout << "\nDenseMatrix2x2.AssignToDiagonal ( VectorVariable2U )";
                DenseMatrix3x3.AssignToDiagonal ( VectorVariable2U );
                DenseMatrix3x3.Out();

                for ( int i = 0, j = 0; i < 2; i++, j++ )
                {
                   _equal( DenseMatrix3x3(i,j), 10., 1E-6 );
                }


                //AssignToDiagonal ( VectorVariable3U )
                cout << "\nAssignToDiagonal ( VectorVariable3U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\n3U Vector variable VectorVariable3U";
                VectorVariable3U = 8;
                VectorVariable3U += 2.;
                VectorVariable3U.Out();

                cout << "\nDenseMatrix3x3.AssignToDiagonal ( VectorVariable3U )";
                DenseMatrix3x3.AssignToDiagonal ( VectorVariable3U );
                DenseMatrix3x3.Out();

                for ( int i = 0, j = 0 ; i < 3; i++, j++ )
                {
                    _equal( DenseMatrix3x3(i,j), 10., 1E-6 );
                }


                //AssignToDiagonal ( size_t diag elmnts, double64 )
                cout << "\nAssignToDiagonal ( size_t diag elmnts, double64 )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out();

                cout << "\nDenseMatrix3x3.AssignToDiagonal ( size_t diag elmnts, double64 )";
                DenseMatrix3x3.AssignToDiagonal(3, 10.);
                DenseMatrix3x3.Out();

                for ( int i = 0, j = 0 ; i < 3; i++, j++ )
                {
                    _equal( DenseMatrix3x3(i,j), 10., 1E-6 );
                }


            //RowCondenseTo( vector<double64> )
                cout << "\nTesting RowCondenseTo( vector<double64> )" << endl;
                cout << "===========================================" << endl;
                //RowCondenseTo ( vector<double64> )
                cout << "\nRowCondenseTo ( vector<double64> )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(9.);
                DenseMatrix3x3.Out();

                vector<double64>::const_iterator itx(x.begin());

                cout << "Test vector x";

                for (vector<double64>::const_iterator
                     itx = x.begin(); itx != x.end(); itx++)
                {
                    cout << "\n" << *itx;
                }

                cout << "\n\nDenseMatrix3x3.RowCondenseTo( x )";
                DenseMatrix3x3.RowCondenseTo( x );

                for (vector<double64>::const_iterator
                     itx = x.begin(); itx != x.end(); itx++)
                {
                    cout << "\n" << *itx ;
                    _equal( *itx, 27., 1E-6 );
                }

                cout << endl;


            //ExportTo ( TensorVariable )
                cout << "\nTesting ExportTo ( TensorVariable )" << endl;
                cout << "===========================================" << endl;

                //ExportTo ( TensorVariable 2U )
                cout << "\nExportTo ( TensorVariable 2U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out();

                cout << "\n2U tensor variable TensorVariable2U";
                TensorVariable2U = 0.;
                TensorVariable2U.Out();

                cout << "\nDenseMatrix2x2 *= ( TensorVariable2U )";
                DenseMatrix3x3.ExportTo(TensorVariable2U);
                TensorVariable3U.Out();

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(i,i), 1., 1E-6 );
                }

                //ExportTo ( TensorVariable 3U )
                cout << "\n\nExportTo ( TensorVariable 3U )" << endl;
                cout << "-------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out();

                cout << "\n3U tensor variable TensorVariable3U";
                TensorVariable3U = 0.;
                TensorVariable3U.Out();

                cout << "\nDenseMatrix3x3 *= ( TensorVariable3U )";
                DenseMatrix3x3.ExportTo(TensorVariable3U);
                TensorVariable3U.Out();

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(i,i), 1., 1E-6 );
                }

            //Assign ( size_t , Point )
                cout << "\n\nTesting Assign ( size_t , Point )" << endl;
                cout << "===========================================" << endl;

                //AssignRow ( size_t i, Point<1U> )
                cout << "\nAssignRow ( size_t i, Point<1U> )" << endl;
                cout << "---------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out();

                cout << "Point Point<1U>";
                Point1U= 13.;
                Point1U.Out();

                cout << "\nAssignRow( 1, Point1U )";
                DenseMatrix3x3.AssignRow( 1, Point1U );
                DenseMatrix3x3.Out();

                _equal( DenseMatrix3x3(1,0), 13., 1E-6 );

                //AssignRow ( size_t i, Point<2U> )
                cout << "AssignRow ( size_t i, Point<2U> )" << endl;
                cout << "---------------------------------" << endl;

                cout << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Zero();
                DenseMatrix3x3.Out();

                cout << "Point Point<2U>";
                Point2U = 11.;
                Point2U.Out();

                cout << "\nAssignRow( 1, Point2U )";
                DenseMatrix3x3.AssignRow( 1, Point2U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(1,i), 11., 1E-6 );
                }

                //AssignRow ( size_t i, Point<3U> )
                cout << "\nAssignRow ( size_t i, Point<3U> )" << endl;
                cout << "---------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out();

                cout << "Point Point<3U>";
                Point3U.Set (13., 12., 11.);
                Point3U.Out();

                cout << "\nAssignRow( 1, Point2U )";
                DenseMatrix3x3.AssignRow( 1, Point3U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(1,i), 13.-i, 1E-6 );
                }

                //AssignCol ( size_t i, Point<1U> )
                cout << "\nAssignCol ( size_t i, Point<1U> )" << endl;
                cout << "---------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out();

                cout << "Point Point<1U>";
                Point1U= 13.;
                Point1U.Out();

                cout << "\nAssignCol( 1, Point1U )";
                DenseMatrix3x3.AssignCol( 1, Point1U );
                DenseMatrix3x3.Out();

                _equal( DenseMatrix3x3(0,1), 13., 1E-6 );

                //AssignCol ( size_t i, Point<2U> )
                cout << "\nAssignCol ( size_t i, Point<2U> )" << endl;
                cout << "---------------------------------" << endl;

                cout << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Zero();
                DenseMatrix3x3.Out();

                cout << "Point Point<2U>";
                Point2U = 11.;
                Point2U.Out();

                cout << "\nAssignCol( 1, Point2U )";
                DenseMatrix3x3.AssignCol( 1, Point2U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(i,1), 11., 1E-6 );
                }

                //AssignCol ( size_t i, Point<3U> )
                cout << "\nAssignCol ( size_t i, Point<3U> )" << endl;
                cout << "---------------------------------" << endl;

                cout << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out();

                cout << "Point Point<3U>";
                Point3U.Set (13., 12., 11.);
                Point3U.Out();

                cout << "\nAssignCol( 1, Point2U )";
                DenseMatrix3x3.AssignCol( 2, Point3U );
                DenseMatrix3x3.Out();

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(i,2), 13.-i, 1E-6 );
                }


                //Testing operator*= ( Point<1U> )
                    cout << "\nTesting operator*= ( Point<1U> )" << endl;
                    cout << "====================================" << endl;

                    //Operator *= ( Point<1U> )
                    cout << "\nOperator *= ( Point<1U> )" << endl;
                    cout << "-------------------------------" << endl;

                    cout << "\n3U test matrix DenseMatrix3x3";
                    DenseMatrix3x3.Resize(1,1);
                    DenseMatrix3x3.Fill(5.);
                    DenseMatrix3x3.Out();

                    cout << "\nPoint<1U> Point1U";
                    Point1U.Out();

                    cout << "\nDenseMatrix3x3 *= ( Point1U )";
                    DenseMatrix3x3 *= ( Point1U );
                    DenseMatrix3x3.Out();

                    _equal( DenseMatrix3x3(0,0), 65., 1E-6 );

                    //Operator *= ( Point<2U> )
                    cout << "\nOperator *= ( Point<2U> )" << endl;
                    cout << "-------------------------------" << endl;

                    cout << "\n3U test matrix DenseMatrix3x3";
                    DenseMatrix3x3.Resize(2,2);
                    DenseMatrix3x3.Fill(5.);
                    DenseMatrix3x3.Out();

                    cout << "\nPoint<2U> Point2U";
                    Point2U.Out();

                    cout << "\nDenseMatrix3x3 *= ( Point2U )";
                    DenseMatrix3x3 *= ( Point2U );
                    DenseMatrix3x3.Out();

                    for ( int i = 0; i < 2; i++ )
                        _equal( DenseMatrix3x3(i,0), 110., 1E-6 );


                    //Operator*= ( Point<3U> )
                    cout << "\nOperator*= ( Point<3U> )" << endl;
                    cout << "-------------------------------" << endl;

                    cout << "\n3U test matrix DenseMatrix3x3";
                    DenseMatrix3x3.Resize(3,3);
                    DenseMatrix3x3.Fill(5.);
                    DenseMatrix3x3.Out();

                    cout << "\nPoint<3U> Point3U";
                    Point3U.Out();

                    cout << "\nDenseMatrix3x3 *= ( Point3U )";
                    DenseMatrix3x3 *= ( Point3U );
                    DenseMatrix3x3.Out();

                    for ( int i = 0; i < 3; i++ )
                        _equal( DenseMatrix3x3(i,0), 180., 1E-6 );


                //Testing L1 matrix norm()
                    cout << "\nTesting NormL1()" << endl;
                    cout << "================" << endl;

                        cout << "\n3U test matrix A";
                        A(2,2) = -99.;
                        A(2,1) = -199.;
                        A.Out();
                        double64 normL1 = A.NormL1();

                        cout << "\nThe maximum absolute column sum norm of test matrix A is " << normL1 << endl;
                        _equal( normL1, 201., 1E-6 );


                //Testing L infinity matrix norm()
                    cout << "\nTesting NormL_Infinity()" << endl;
                    cout << "================" << endl;

                        cout << "\n3U test matrix A";
                        A(2,2) = -99.;
                        A(2,1) = -199.;
                        A.Out();
                        double64 normInf = A.NormL_Infinity();

                        cout << "\nThe maximum absolute row sum norm of test matrix A is " << normInf << endl;
                        _equal( normInf, 303., 1E-6 );

/*
                //Testing In()
                    cout << "\nTesting In()" << endl;
                    cout << "============" << endl;


                        cout << "\n3U test matrix DenseMatrix3x3";
                        DenseMatrix3x3.Resize(3,3);
                        int rows = 2;
                        int cols = 2;
                        DenseMatrix3x3.In();
*/

    }


} // end csmp

