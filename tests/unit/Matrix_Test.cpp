// unit test
// Created by : Siroos and Georg
// Date of Modification: 21.11.2010

#include "Matrix_Test.h"

using namespace std;

namespace csmp
{
    Matrix_Test::Matrix_Test()
        :
        A(4,4),
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
        Matrix3x3(3,3)

    {
    }

    Matrix_Test::~Matrix_Test()
    {
    }

    void Matrix_Test::run()
    {
        std::ostream& os = getInfoStream();

    os << "\n==============";
    os << "\nTesting Matrix" << endl;
    os << "==============" << endl;

    //Testing Identity functions
        os << "\nTesting Identity function" << endl;
        os << "=========================" << endl;


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

        os << "\nB.Identity()";

        B.Identity();
        B.Out(os);

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                if ( i==j ) _equal( B(i,j), 1., 1E-6 );
                else        _equal( B(i,j), 0., 1E-6 );
            }
        }

    //Testing + and - operator

        os << "\nTesting + and += operator" << endl;
        os << "======================" << endl;

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

        os << "\nTest Matrix A";
        A.Out(os);

        os << "Test Matrix B";
        B.Fill(1.);
        B.Out(os);

        os << "Matrix C = A + B";
        C = A + B;
        C.Out(os);
        E = C;

        os << "Matrix A += B =! C";
        A += B;
        A.Out(os);

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
                //os << C(i, j) << "  ";
                //os << A(i, j) << endl;
                _equal( A(i,j ), ResC(i,j), 1E-6 );
                _equal( C(i,j ), ResC(i,j), 1E-6 );
            }
        }

        os << "Matrix D = C - B";
        D = C - B;
        D.Out(os);

        os << "Matrix A -= B =! D";
        A -= B;
        A.Out(os);

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
                //os << D(i, j) << "  ";
                //os << A(i, j) << endl;
                _equal( D(i,j), ResA(i,j), 1E-6 );
                _equal( A(i,j), ResA(i,j), 1E-6 );
            }
        }


//Testing * operators
    os << "\nTesting * operators" << endl;
    os << "===================" << endl;

        //Matrix-scalar multiplication
        os << "\nMatrix-scalar multiplication method" << endl;
        os << "-----------------------------------" << endl;

        os << "\nTest matrix A";
        A.Out(os);

        os << "Matrix A *= ( 3. )";
        C = A;
        C *= ( 3. );
        C.Out(os);

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
        os << "\nFirst vector-matrix multiplication method: y = A * x" << endl;
        os << "----------------------------------------------------" << endl;

        os << "\nTest matrix A";

        A.Out(os);

        x.push_back(1.);
        x.push_back(2.);
        x.push_back(0.);
        x.push_back(1.5);

        os << "Test vector x =   ";

        vector<double64>::const_iterator it( x.begin() );
        for ( vector<double64>::const_iterator
             it = x.begin(); it != x.end(); it++ )
        {
            os << *it << setw(5);
        }

        os << "\n\ny = A * x";
        y = A * x;

        vector<double64>::const_iterator it1(y.begin());

        os << "\n\ny = ";

        for (vector<double64>::const_iterator
             it1 = y.begin(); it1 != y.end(); it1++ )
        {
            os << *it1 << setw(5);
        }

        sol_y.push_back(2.5);
        sol_y.push_back(8.);
        sol_y.push_back(7.);
        sol_y.push_back(9.5);

        vector<double64>::const_iterator itsol_y(sol_y.begin());

        for (vector<double64>::const_iterator
             itsol_y = sol_y.begin(); itsol_y != sol_y.end(); itsol_y++, it1++ )
            {
                //os << "\n" << *itsol_y << setw(7) << *it1;
                _equal( *itsol_y, *it1, 1E-6);
            }


        //\nSecond matrix-vector multiplication method
        os << "\n\nSecond matrix-vector multiplication method: A *= x" << endl;
        os << "--------------------------------------------------" << endl;

        B = A;
        os << "A *= x";
        B.operator *= (x);

        os << "\nA = ";

        for ( int i = 0; i < 4; i++ )
        {
            os << B(i,0) << setw(5);

            //Retrieving values of matrix A for comparison
            solB.push_back(B(i,0));
        }

        vector<double64>::const_iterator itB( solB.begin() );

        //os << endl << "\nA  " << setw(5) << "y";

        for (vector<double64>::const_iterator
             itB = solB.begin(); itB != solB.end(); itB++, itsol_y++ )
        {
            //os << "\n" << *itB << setw(7) << *it1;
            _equal( *itB, *itsol_y, 1E-6 );
        }

        os << endl;


        os << "\nThird matrix-vector multiplication method: A *= ( C array )" << endl;
        os << "--------------------------------------------------" << endl;

        D = A;

        os << "\nTest matrix A";
        D.Out(os);


        double64 Ca[4] = { 1., 2., 0., 1.5 };

        os << "Test array Ca";

        for ( int i = 0; i < 4; i++ )
            os << "\n" << Ca[i];


        os << "\n\nMatrix A *= ( C array )";
        D.operator *=( Ca );
        D.Out(os);

        for ( int i = 0; i < 4; i++ )
        {
            os << D(i,0) << setw(5);

            //Retrieving values of matrix A for comparison
            solD.push_back(D(i,0));
        }

        vector<double64>::const_iterator itD( solD.begin() );

        //os << endl << "\nA  " << setw(5) << "y";

        for (vector<double64>::const_iterator
             itsol_y = sol_y.begin(); itsol_y != sol_y.end(); itsol_y++, itD++ )
            {
                //os << "\n" << *itsol_y << setw(7) << *itD;
                _equal( *itsol_y, *itD, 1E-6);
            }

        os << endl;


        //Matrix-by-matrix multiplication
        os << "\n\nFirst matrix-matrix multiplication method: A * B" << endl;
        os << "------------------------------------------------" << endl;

        os << "\nTest matrix A";
        A.Out(os);

        os << "Test matrix B";
        B = E;
        B.Out(os);

        os << "C = A * B";
        C = A * B;
        C.Out(os);

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

        os << "\nSecond matrix-matrix multiplication method: A *= B" << endl;
        os << "--------------------------------------------------" << endl;

        D = A;

        os << "\nTest matrix A";
        D.Out(os);

        os << "\nTest matrix B";
        B.Out(os);

        os << "\nMatrix A *= B";
 //D.operator *=(B);
        D*=B;
        D.Out(os);

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( D(i,j), ResAB(i,j), 1E-6 );
            }
        }


    //Testing fill operators
        os << "\nTesting fill functions" << endl;
        os << "======================" << endl;

        //fill
        os << "\nMatrix fill()" << endl;
        os << "-------------" << endl;

        os << "\nB.fill(1.)";
        B.Fill(1.);
        B.Out(os);

        C = 1.;

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(i,j), C(i,j), 1E-6 );
            }
        }

        //fillCol
        os << "\nMatrix fillCol()" << endl;
        os << "---------------" << endl;

        os << "\nB.fillCol(2,5.5)";
        B.FillCol(2,5.5);
        B.Out(os);

        for ( int i = 0; i < 4; i++ )
        {
            _equal( B(i,2), 5.5, 1E-6 );
        }

        //fillRow
        os << "\nMatrix fillRow()" << endl;
        os << "----------------" << endl;

        os << "\nB.fillRow(2,3.3)";
        B.FillRow(2,3.3);
        B.Out(os);

        for ( int j = 0; j < 4; j++ )
        {
            _equal( B(2,j), 3.3, 1E-6 );
        }

        C = B;

    //Testing summation functions
        os << "\nTesting summation function" << endl;
        os << "==========================" << endl;

        //ColSum()
        os << "\nMatrix ColSum()" << endl;
        os << "---------------" << endl;

        os << "\nB.ColSum(2) = " << B.ColSum(2) << endl;

        double64 ColSumB2 = B.ColSum(2);
        double64 CheckColSumB2 = 0.;

        for ( int i = 0; i < 4; i++ )
        {
            CheckColSumB2 += B(i,2);
        }

        os << "\nB.CheckColSum(2) = " << CheckColSumB2 << endl;
        _equal( CheckColSumB2, ColSumB2, 1E-6 );


        //RowSum
        os << "\nMatrix RowSum()" << endl;
        os << "---------------" << endl;

        os << "\nB.RowSum(2) = " << B.RowSum(2) << endl;

        double64 RowSumB2 = B.RowSum(2);
        double64 CheckRowSumB2 = 0.;

        for ( int j = 0; j < 4; j++ )
        {
            CheckRowSumB2 += B(2,j);
        }

        os << "\nB.CheckRowSum(2) = " << CheckRowSumB2 << endl;
        _equal( CheckRowSumB2, RowSumB2, 1E-6 );

        //Testing Zero functions
            os << "\nTesting Zero functions" << endl;
            os << "======================" << endl;

            //ZeroCol
            os << "\nMatrix ZeroCol()" << endl;
            os << "----------------" << endl;

            os << "\nB.ZeroCol(2)";
            B.ZeroCol(2);
            B.Out(os);

            for ( int i = 0; i < 4; i++ )
            {
                _equal( B(i,2), 0., 1E-6 );
            }

            //ZeroRow
            os << "\nMatrix ZeroRow()" << endl;
            os << "----------------" << endl;

            os << "\nB.ZeroRow(2)";
            B.ZeroRow(2);
            B.Out(os);

            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(2,j), 0., 1E-6 );
            }

            //Zero
            os << "\nMatrix Zero()" << endl;
            os << "-------------" << endl;

            os << "\nB.Zero()";
            B.Zero();
            B.Out(os);

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( B(i,j), 0., 1E-6 );
                }
            }


        //Testing operator= ()
            os << "\nTesting operator= ()" << endl;
            os << "====================" << endl;

            //Testing operator= (const Matrix& )
            os << "\nTesting operator= ( const Matrix& )" << endl;
            os << "---------------------------------------" << endl;

            os << "\nMatrix B";
            B.Out(os);

            os << "\nMatrix B = C";
            B = C;
            B.Out(os);

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( B(i,j), C(i,j), 1E-6 );
                }
            }


        //Testing transposed functions
            os << "\nTesting transposed matrix functions" << endl;
            os << "===================================" << endl;

            //Transposed matrix
            os << "\nMatrix Transposed()" << endl;
            os << "-------------------" << endl;

            os << "\nTest matrix A";
            A.Out(os);

            A.Transposed(F);
            os << "\nF = A^T = A.Transposed(F)";
            F.Out(os);

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( F(i,j), A(j,i), 1E-6 );
                }
            }

            //Transposed product
            os << "\nMatrix TransposedProduct()" << endl;
            os << "--------------------------" << endl;

            os << "\nTest matrix A";
            A.Out(os);

            A.TransposedProduct(H);
            os << "\nH = A^T * A = A.TransposedProduct(F)";
            H.Out(os);

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
            os << "\nMatrix MultiplyWithTransposedOf()" << endl;
            os << "---------------------------------" << endl;

            os << "\nTest matrix A";
            A.Out(os);

            os << "\nTest matrix B";
            B.Out(os);

            os << "\nF = A B^T = A.MultiplyWithTransposedOf( B, F )";
            A.MultiplyWithTransposedOf( B, F );
            F.Out(os);

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
            os << "\nMatrix MultiplyTransposedOfWith()" << endl;
            os << "---------------------------------" << endl;

            os << "\nTest matrix A";
            A.Out(os);

            os << "\nTest matrix B";
            B.Out(os);

            os << "\nG = A^T B = A.MultiplyTransposedOfWith( B, G )";
            A.MultiplyTransposedOfWith( B, G );
            G.Out(os);

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
            os << "\nTesting operator= ()" << endl;
            os << "====================================" << endl;

            //Operator= ( TensorVariable 2U )
            os << "\nOperator= ( TensorVariable 2U )" << endl;
            os << "-------------------------------" << endl;

            os << "\n2U test matrix Matrix2x2";
            Matrix2x2.Out(os);

            double64 val = 16.;

            os << "\nMatrix2x2 = val";
            Matrix2x2 = val;
            Matrix2x2.Out(os);

            for ( int i = 0; i < 2; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    _equal( Matrix2x2(i,j), 16., 1E-6 );
                }
            }


            //Operator= ( TensorVariable 2U )
            os << "\nOperator= ( TensorVariable 2U )" << endl;
            os << "-------------------------------" << endl;

            os << "\n2U test matrix Matrix2x2";
            Matrix2x2.Fill(5.);
            Matrix2x2.Out(os);

            os << "\n2U tensor variable TensorVariable2U";
            TensorVariable2U.Identity();
            TensorVariable2U += 1.;
            TensorVariable2U.Out(os);

            os << "\nMatrix2x2 = TensorVariable2U";
            Matrix2x2 = TensorVariable2U;
            Matrix2x2.Out(os);

            for ( int i = 0; i < 2; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    if ( i==j ) _equal( Matrix2x2(i,j), 2., 1E-6 );
                    else        _equal( Matrix2x2(i,j), 1., 1E-6 );
                }
            }


            //Operator= ( TensorVariable 3U )
            os << "\nOperator= ( TensorVariable 3U )" << endl;
            os << "-------------------------------" << endl;

            os << "\n3U test matrix Matrix3x3";
            Matrix3x3.Fill(5.);
            Matrix3x3.Out(os);

            os << "\n3U tensor variable TensorVariable3U";
            TensorVariable3U.Identity();
            TensorVariable3U += 2.;
            TensorVariable3U.Out(os);

            os << "\nMatrix3x3 = TensorVariable3U";
            Matrix3x3 = TensorVariable3U;
            Matrix3x3.Out(os);

            for ( int i = 0; i < 3; i++ )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    if ( i==j ) _equal( Matrix3x3(i,j), 3., 1E-6 );
                    else        _equal( Matrix3x3(i,j), 2., 1E-6 );
                }
            }

            //Testing operator*= ( TensorVariable )
                os << "\nTesting operator*= ( TensorVariable )" << endl;
                os << "====================================" << endl;

                //Operator *= ( double64 )
                os << "\nOperator *= ( double64 )" << endl;
                os << "-------------------------------" << endl;

                os << "\n3U test matrix Matrix3x3";
                Matrix3x3.Resize(3,3);
                Matrix3x3.Fill(5.);
                Matrix3x3.Out(os);

                os << "\nMatrix3x3 *= ( 5 )";
                Matrix3x3 *= ( 5 );
                Matrix3x3.Out(os);

                _equal( Matrix3x3(0,0), 25., 1E-6 );

                //Operator *= ( TensorVariable 2U )
                os << "\nOperator *= ( TensorVariable 2U )" << endl;
                os << "-------------------------------" << endl;

                os << "\n3U test matrix Matrix3x3";
                Matrix3x3.Resize(2,2);
                Matrix3x3.Fill(5.);
                Matrix3x3.Out(os);

                os << "\n2U tensor variable TensorVariable2U";
                TensorVariable2U.Identity();
                TensorVariable2U += 1.;
                TensorVariable2U.Out(os);

                os << "\nMatrix3x3 *= ( TensorVariable2U )";
                Matrix3x3 *= ( TensorVariable2U );
                Matrix3x3.Out(os);

                for ( int i = 0; i < 2; i++ )
                {
                    for ( int j = 0; j < 2; j++ )
                    {
                        _equal( Matrix3x3(i,j), 15., 1E-6 );
                    }
                }

                //Operator*= ( TensorVariable 3U )
                os << "\nOperator*= ( TensorVariable 3U )" << endl;
                os << "-------------------------------" << endl;

                os << "\n3U test matrix Matrix3x3";
                Matrix3x3.Resize(3,3);
                Matrix3x3.Fill(5.);
                Matrix3x3.Out(os);

                os << "\n3U tensor variable TensorVariable3U";
                TensorVariable3U.Identity();
                TensorVariable3U += 2.;
                TensorVariable3U.Out(os);

                os << "\nMatrix3x3 *= ( TensorVariable3U )";
                Matrix3x3 *= ( TensorVariable3U );
                Matrix3x3.Out(os);

                for ( int i = 0; i < 3; i++ )
                {
                    for ( int j = 0; j < 3; j++ )
                    {
                        _equal( Matrix3x3(i,j), 35., 1E-6 );
                    }
                }


            //Testing operator*= ( VectorVariable )
                os << "\nTesting operator*= ( VectorVariable )" << endl;
                os << "====================================" << endl;

                //Operator *= ( VectorVariable 2U )
                os << "\nOperator *= ( VectorVariable 2U )" << endl;
                os << "-------------------------------" << endl;

                os << "\n2U test matrix Matrix2x2";
                Matrix3x3.Resize(2,2);
                Matrix3x3.Fill(5.);
                Matrix3x3.Out(os);

                os << "\n2U Vector variable VectorVariable2U";
                VectorVariable2U = 5;
                VectorVariable2U += 1.;
                VectorVariable2U.Out(os);

                os << "\nMatrix2x2 *= ( VectorVariable2U )";
                Matrix3x3 *= ( VectorVariable2U );
                Matrix3x3.Out(os);

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( Matrix3x3(i,0), 60., 1E-6 );

                }


                //Operator*= ( VectorVariable 3U )
                os << "\nOperator*= ( VectorVariable 3U )" << endl;
                os << "-------------------------------" << endl;

                os << "\n3U test matrix Matrix3x3";
                Matrix3x3.Resize(3,3);
                Matrix3x3.Fill(5.);
                Matrix3x3.Out(os);

                os << "\n3U Vector variable VectorVariable3U";
                VectorVariable3U = 5;
                VectorVariable3U += 2.;
                VectorVariable3U.Out(os);

                os << "\nMatrix3x3 *= ( VectorVariable3U )";
                Matrix3x3 *= ( VectorVariable3U );
                Matrix3x3.Out(os);

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( Matrix3x3(i,0), 105., 1E-6 );
                }


            //RowCondenseTo( vector<double64> )
                os << "\nTesting RowCondenseTo( vector<double64> )" << endl;
                os << "===========================================" << endl;
                //RowCondenseTo ( vector<double64> )
                os << "\nRowCondenseTo ( vector<double64> )" << endl;
                os << "-------------------------------" << endl;

                os << "\n3U test matrix Matrix3x3";
                Matrix3x3.Resize(3,3);
                Matrix3x3.Fill(9.);
                Matrix3x3.Out(os);

                vector<double64>::const_iterator itx(x.begin());

                os << "Test vector x";

                for (vector<double64>::const_iterator
                     itx = x.begin(); itx != x.end(); itx++)
                {
                    os << "\n" << *itx;
                }

                os << "\n\nMatrix3x3.RowCondenseTo( x )";
                Matrix3x3.RowCondenseTo( x );

                for (vector<double64>::const_iterator
                     itx = x.begin(); itx != x.end(); itx++)
                {
                    os << "\n" << *itx ;
                    _equal( *itx, 27., 1E-6 );
                }

                os << endl;


            //ExportTo ( TensorVariable )
                os << "\nTesting ExportTo ( TensorVariable )" << endl;
                os << "===========================================" << endl;

                //ExportTo ( TensorVariable 2U )
                os << "\nExportTo ( TensorVariable 2U )" << endl;
                os << "-------------------------------" << endl;

                os << "\n2U test matrix Matrix2x2";
                Matrix3x3.Resize(2,2);
                Matrix3x3.Identity();
                Matrix3x3.Out(os);

<<<<<<< HEAD
                cout << "\n2U tensor variable TensorVariable2U";
                TensorVariable2U=0.;
                TensorVariable2U.Out();
=======
                os << "\n2U tensor variable TensorVariable2U";
                TensorVariable2U.Zero();
                TensorVariable2U.Out(os);
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3

                os << "\nMatrix2x2 *= ( TensorVariable2U )";
                Matrix3x3.ExportTo(TensorVariable2U);
                TensorVariable3U.Out(os);

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( Matrix3x3(i,i), 1., 1E-6 );
                }


                //ExportTo ( TensorVariable 3U )
                os << "\n\nExportTo ( TensorVariable 3U )" << endl;
                os << "-------------------------------" << endl;

                os << "\n3U test matrix Matrix3x3";
                Matrix3x3.Resize(3,3);
                Matrix3x3.Identity();
                Matrix3x3.Out(os);

<<<<<<< HEAD
                cout << "\n3U tensor variable TensorVariable3U";
                TensorVariable3U=0.;
                TensorVariable3U.Out();
=======
                os << "\n3U tensor variable TensorVariable3U";
                TensorVariable3U.Zero();
                TensorVariable3U.Out(os);
>>>>>>> b71117cb444b1539e747fa5855ab341062d3c4b3

                os << "\nMatrix3x3 *= ( TensorVariable3U )";
                Matrix3x3.ExportTo(TensorVariable3U);
                TensorVariable3U.Out(os);

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( Matrix3x3(i,i), 1., 1E-6 );
                }

                //Testing L1 matrix norm()
                    os << "\nTesting NormL1()" << endl;
                    os << "================" << endl;

                        os << "\n3U test matrix A";
                        A(2,2) = -99.;
                        A(2,1) = -199.;
                        A.Out(os);
                        double64 normL1 = A.NormL1();

                        os << "\nThe maximum absolute column sum norm of test matrix A is " << normL1 << endl;
                        _equal( normL1, 201., 1E-6 );


                //Testing L infinity matrix norm()
                    os << "\nTesting NormL_Infinity()" << endl;
                    os << "================" << endl;

                        os << "\n3U test matrix A";
                        A(2,2) = -99.;
                        A(2,1) = -199.;
                        A.Out(os);
                        double64 normInf = A.NormL_Infinity();

                        os << "\nThe maximum absolute row sum norm of test matrix A is " << normInf << endl;
                        _equal( normInf, 303., 1E-6 );


                //Testing In() has been done manually and confirmed
                //  Otherwise script based In() call must be answered



    }


} // end csmp

