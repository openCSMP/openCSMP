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

    getInfoStream() << "\n===================";
    getInfoStream() << "\nTesting DenseMatrix" << endl;
    getInfoStream() << "===================" << endl;

    //Testing Identity functions
        getInfoStream() << "\nTesting Identity function" << endl;
        getInfoStream() << "=========================" << endl;

        getInfoStream() << "\nB.Identity()";
        B.Identity();
        B.Out(getInfoStream());

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                if ( i==j ) _equal( B(i,j), 1., 1E-6 );
                else        _equal( B(i,j), 0., 1E-6 );
            }
        }

    //Testing + and - operator

        getInfoStream() << "\nTesting + and += operator" << endl;
        getInfoStream() << "======================" << endl;

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

        getInfoStream() << "\nTest Matrix A";
        A.Out(getInfoStream());

        getInfoStream() << "Test Matrix B";
        B.Fill(1.);
        B.Out(getInfoStream());

        getInfoStream() << "Matrix C = A + B";
        C = A + B;
        C.Out(getInfoStream());
        E = C;

        getInfoStream() << "Matrix A += B =! C";
        A += B;
        A.Out(getInfoStream());

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
                //getInfoStream() << C(i, j) << "  ";
                //getInfoStream() << A(i, j) << endl;
                _equal( A(i,j ), ResC(i,j), 1E-6 );
                _equal( C(i,j ), ResC(i,j), 1E-6 );
            }
        }

        getInfoStream() << "Matrix D = C - B";
        D = C - B;
        D.Out(getInfoStream());

        getInfoStream() << "Matrix A -= B =! D";
        A -= B;
        A.Out(getInfoStream());

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
                //getInfoStream() << D(i, j) << "  ";
                //getInfoStream() << A(i, j) << endl;
                _equal( D(i,j), ResA(i,j), 1E-6 );
                _equal( A(i,j), ResA(i,j), 1E-6 );
            }
        }


//Testing * operators
    getInfoStream() << "\nTesting * operators" << endl;
    getInfoStream() << "===================" << endl;

        //Matrix-scalar multiplication
        getInfoStream() << "\nMatrix-scalar multiplication method" << endl;
        getInfoStream() << "-----------------------------------" << endl;

        getInfoStream() << "\nTest matrix A";
        A.Out(getInfoStream());

        getInfoStream() << "Matrix A *= ( 3 )";
        C = A;
        C.operator *= ( 3. );
        C.Out(getInfoStream());

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
        getInfoStream() << "\nMatrix-scalar multiplication method" << endl;
        getInfoStream() << "-----------------------------------" << endl;

        getInfoStream() << "\nTest matrix A";
        A.Out(getInfoStream());

        getInfoStream() << "Matrix A *= ( 3 )";
        C = A;
        C.operator *= ( 3. );
        C.Out(getInfoStream());


        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( C(i,j), ResC(i,j), 1E-6 );
            }
        }


        //Matrix multiplication by vector
        //First vector-matrix multiplication method
        getInfoStream() << "\nFirst vector-matrix multiplication method: y = A * x" << endl;
        getInfoStream() << "----------------------------------------------------" << endl;

        getInfoStream() << "\nTest matrix A";

        A.Out(getInfoStream());

        x.push_back(1.);
        x.push_back(2.);
        x.push_back(0.);
        x.push_back(1.5);

        getInfoStream() << "Test vector x =   ";

        vector<double64>::const_iterator it( x.begin() );
        for ( vector<double64>::const_iterator
             it = x.begin(); it != x.end(); it++ )
        {
            getInfoStream() << *it << setw(5);
        }

        getInfoStream() << "\n\ny = A * x";
        y = A * x;

        vector<double64>::const_iterator it1(y.begin());

        getInfoStream() << "\n\ny = ";

        for (vector<double64>::const_iterator
             it1 = y.begin(); it1 != y.end(); it1++ )
        {
            getInfoStream() << *it1 << setw(5);
        }

        sol_y.push_back(2.5);
        sol_y.push_back(8.);
        sol_y.push_back(7.);
        sol_y.push_back(9.5);

        vector<double64>::const_iterator itsol_y(sol_y.begin());

        for (vector<double64>::const_iterator
             itsol_y = sol_y.begin(); itsol_y != sol_y.end(); itsol_y++, it1++ )
            {
                //getInfoStream() << "\n" << *itsol_y << setw(7) << *it1;
                _equal( *itsol_y, *it1, 1E-6);
            }


        //\nSecond matrix-vector multiplication method
        getInfoStream() << "\n\nSecond matrix-vector multiplication method: A *= x" << endl;
        getInfoStream() << "--------------------------------------------------" << endl;

        B = A;
        getInfoStream() << "A *= x";
        B.operator *= (x);

        getInfoStream() << "\nA = ";

        for ( int i = 0; i < 4; i++ )
        {
            getInfoStream() << B(i,0) << setw(5);

            //Retrieving values of matrix A for comparison
            solB.push_back(B(i,0));
        }

        vector<double64>::const_iterator itB( solB.begin() );

        //getInfoStream() << endl << "\nA  " << setw(5) << "y";

        for (vector<double64>::const_iterator
             itB = solB.begin(); itB != solB.end(); itB++, itsol_y++ )
        {
            //getInfoStream() << "\n" << *itB << setw(7) << *it1;
            _equal( *itB, *itsol_y, 1E-6 );
        }

        getInfoStream() << endl;


        getInfoStream() << "\nThird matrix-matrix multiplication method: A *= ( C array )" << endl;
        getInfoStream() << "--------------------------------------------------" << endl;

        D = A;

        getInfoStream() << "\nTest matrix A";
        D.Out(getInfoStream());


        double64 Ca[4] = { 1., 2., 0., 1.5 };

        getInfoStream() << "Test array Ca";

        for ( int i = 0; i < 4; i++ )
            getInfoStream() << "\n" << Ca[i];


        getInfoStream() << "\n\nMatrix A *= ( C array )";
        D.operator *=( Ca );
        D.Out(getInfoStream());

        for ( int i = 0; i < 4; i++ )
        {
            getInfoStream() << D(i,0) << setw(5);

            //Retrieving values of matrix A for comparison
            solD.push_back(D(i,0));
        }

        vector<double64>::const_iterator itD( solD.begin() );

        //getInfoStream() << endl << "\nA  " << setw(5) << "y";

        for (vector<double64>::const_iterator
             itsol_y = sol_y.begin(); itsol_y != sol_y.end(); itsol_y++, itD++ )
            {
                //getInfoStream() << "\n" << *itsol_y << setw(7) << *itD;
                _equal( *itsol_y, *itD, 1E-6);
            }

        getInfoStream() << endl;


        //Matrix-by-matrix multiplication
        getInfoStream() << "\n\nFirst matrix-matrix multiplication method: A * B" << endl;
        getInfoStream() << "------------------------------------------------" << endl;

        getInfoStream() << "\nTest matrix A";
        A.Out(getInfoStream());

        getInfoStream() << "Test matrix B";
        B = E;
        B.Out(getInfoStream());

        getInfoStream() << "C = A * B";
        C = A * B;
        C.Out(getInfoStream());

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

        getInfoStream() << "\nSecond matrix-matrix multiplication method: A *= B" << endl;
        getInfoStream() << "--------------------------------------------------" << endl;

        D = A;

        getInfoStream() << "\nTest matrix A";
        D.Out(getInfoStream());

        getInfoStream() << "\nMatrix A *= B";
        D.operator *=(B);
        D.Out(getInfoStream());

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( D(i,j), ResAB(i,j), 1E-6 );
            }
        }


    //Testing fill operators
        getInfoStream() << "\nTesting fill functions" << endl;
        getInfoStream() << "======================" << endl;

        //fill
        getInfoStream() << "\nMatrix fill()" << endl;
        getInfoStream() << "-------------" << endl;

        getInfoStream() << "\nB.fill(1.)";
        B.Fill(1.);
        B.Out(getInfoStream());

        C = 1.;

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(i,j), C(i,j), 1E-6 );
            }
        }

        //fillCol
        getInfoStream() << "\nMatrix fillCol()" << endl;
        getInfoStream() << "---------------" << endl;

        getInfoStream() << "\nB.fillCol(2,5.5)";
        B.FillCol(2,5.5);
        B.Out(getInfoStream());

        for ( int i = 0; i < 4; i++ )
        {
            _equal( B(i,2), 5.5, 1E-6 );
        }

        //fillRow
        getInfoStream() << "\nMatrix fillRow()" << endl;
        getInfoStream() << "----------------" << endl;

        getInfoStream() << "\nB.fillRow(2,3.3)";
        B.FillRow(2,3.3);
        B.Out(getInfoStream());

        for ( int j = 0; j < 4; j++ )
        {
            _equal( B(2,j), 3.3, 1E-6 );
        }

        C = B;

    //Testing summation functions
        getInfoStream() << "\nTesting summation function" << endl;
        getInfoStream() << "==========================" << endl;

        //ColSum()
        getInfoStream() << "\nMatrix ColSum()" << endl;
        getInfoStream() << "---------------" << endl;

        getInfoStream() << "\nB.ColSum(2) = " << B.ColSum(2) << endl;

        double64 ColSumB2 = B.ColSum(2);
        double64 CheckColSumB2 = 0.;

        for ( int i = 0; i < 4; i++ )
        {
            CheckColSumB2 += B(i,2);
        }

        getInfoStream() << "\nB.CheckColSum(2) = " << CheckColSumB2 << endl;
        _equal( CheckColSumB2, ColSumB2, 1E-6 );


        //RowSum
        getInfoStream() << "\nMatrix RowSum()" << endl;
        getInfoStream() << "---------------" << endl;

        getInfoStream() << "\nB.RowSum(2) = " << B.RowSum(2) << endl;

        double64 RowSumB2 = B.RowSum(2);
        double64 CheckRowSumB2 = 0.;

        for ( int j = 0; j < 4; j++ )
        {
            CheckRowSumB2 += B(2,j);
        }

        getInfoStream() << "\nB.CheckRowSum(2) = " << CheckRowSumB2 << endl;
        _equal( CheckRowSumB2, RowSumB2, 1E-6 );

        //Testing Zero functions
            getInfoStream() << "\nTesting Zero functions" << endl;
            getInfoStream() << "======================" << endl;

            //ZeroCol
            getInfoStream() << "\nMatrix ZeroCol()" << endl;
            getInfoStream() << "----------------" << endl;

            getInfoStream() << "\nB.ZeroCol(2)";
            B.ZeroCol(2);
            B.Out(getInfoStream());

            for ( int i = 0; i < 4; i++ )
            {
                _equal( B(i,2), 0., 1E-6 );
            }

            //ZeroRow
            getInfoStream() << "\nMatrix ZeroRow()" << endl;
            getInfoStream() << "----------------" << endl;

            getInfoStream() << "\nB.ZeroRow(2)";
            B.ZeroRow(2);
            B.Out(getInfoStream());

            for ( int j = 0; j < 4; j++ )
            {
                _equal( B(2,j), 0., 1E-6 );
            }

            //Zero
            getInfoStream() << "\nMatrix Zero()" << endl;
            getInfoStream() << "-------------" << endl;

            getInfoStream() << "\nB.Zero()";
            B.Zero();
            B.Out(getInfoStream());

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( B(i,j), 0., 1E-6 );
                }
            }


        //Testing operator= ()
            getInfoStream() << "\nTesting operator= ()" << endl;
            getInfoStream() << "====================" << endl;

            //Testing operator= (const DenseMatrix& )
            getInfoStream() << "\nTesting operator= ( const DenseMatrix& )" << endl;
            getInfoStream() << "---------------------------------------" << endl;

            getInfoStream() << "\nMatrix B";
            B.Out(getInfoStream());

            getInfoStream() << "\nMatrix B = C";
            B = C;
            B.Out(getInfoStream());

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( B(i,j), C(i,j), 1E-6 );
                }
            }


        //Testing transposed functions
            getInfoStream() << "\nTesting transposed matrix functions" << endl;
            getInfoStream() << "===================================" << endl;

            //Transposed matrix
            getInfoStream() << "\nMatrix Transposed()" << endl;
            getInfoStream() << "-------------------" << endl;

            getInfoStream() << "\nTest matrix A";
            A.Out(getInfoStream());

            A.Transposed(F);
            getInfoStream() << "\nF = A^T = A.Transposed(F)";
            F.Out(getInfoStream());

            for ( int i = 0; i < 4; i++ )
            {
                for ( int j = 0; j < 4; j++ )
                {
                    _equal( F(i,j), A(j,i), 1E-6 );
                }
            }

            //MultiplyWithTransposedOf
            getInfoStream() << "\nMatrix MultiplyWithTransposedOf()" << endl;
            getInfoStream() << "---------------------------------" << endl;

            getInfoStream() << "\nTest matrix A";
            A.Out(getInfoStream());


            getInfoStream() << "\nF = A B^T = A.MultiplyWithTransposedOf( B, F )";
            A.MultiplyWithTransposedOf( B, F );
            F.Out(getInfoStream());

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
            getInfoStream() << "\nMatrix MultiplyTransposedOfWith()" << endl;
            getInfoStream() << "---------------------------------" << endl;

            getInfoStream() << "\nTest matrix A";
            A.Out(getInfoStream());

            getInfoStream() << "\nTest matrix B";
            B.Out(getInfoStream());

            getInfoStream() << "\nG = A^T B = A.MultiplyTransposedOfWith( B, G )";
            A.MultiplyTransposedOfWith( B, G );
            G.Out(getInfoStream());

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
            getInfoStream() << "\nTesting operator= ( TensorVariable )" << endl;
            getInfoStream() << "====================================" << endl;

            //Operator= ( TensorVariable 1U )
            getInfoStream() << "\nOperator= ( TensorVariable 1U )" << endl;
            getInfoStream() << "-------------------------------" << endl;

            getInfoStream() << "\n3U test matrix DenseMatrix3x3";
            DenseMatrix3x3.Fill(5.);
            DenseMatrix3x3.Out(getInfoStream());

            getInfoStream() << "\n1U tensor variable TensorVariable1U";
            TensorVariable1U.Identity();
            TensorVariable1U.Out(getInfoStream());

            getInfoStream() << "\nDenseMatrix3x3 = TensorVariable1U";
            DenseMatrix3x3 = TensorVariable1U;
            DenseMatrix3x3.Out(getInfoStream());

            _equal( DenseMatrix3x3(0,0), 1., 1E-6 );

            //Operator= ( TensorVariable 2U )
            getInfoStream() << "\nOperator= ( TensorVariable 2U )" << endl;
            getInfoStream() << "-------------------------------" << endl;

            getInfoStream() << "\n3U test matrix DenseMatrix3x3";
            DenseMatrix3x3.Fill(5.);
            DenseMatrix3x3.Out(getInfoStream());

            getInfoStream() << "\n2U tensor variable TensorVariable2U";
            TensorVariable2U.Identity();
            TensorVariable2U += 1.;
            TensorVariable2U.Out(getInfoStream());

            getInfoStream() << "\nDenseMatrix3x3 = TensorVariable2U";
            DenseMatrix3x3 = TensorVariable2U;
            DenseMatrix3x3.Out(getInfoStream());

            for ( int i = 0; i < 2; i++ )
            {
                for ( int j = 0; j < 2; j++ )
                {
                    if ( i==j ) _equal( DenseMatrix3x3(i,j), 2., 1E-6 );
                    else        _equal( DenseMatrix3x3(i,j), 1., 1E-6 );
                }
            }


            //Operator= ( TensorVariable 3U )
            getInfoStream() << "\nOperator= ( TensorVariable 3U )" << endl;
            getInfoStream() << "-------------------------------" << endl;

            getInfoStream() << "\n3U test matrix DenseMatrix3x3";
            DenseMatrix3x3.Fill(5.);
            DenseMatrix3x3.Out(getInfoStream());

            getInfoStream() << "\n3U tensor variable TensorVariable3U";
            TensorVariable3U.Identity();
            TensorVariable3U += 2.;
            TensorVariable3U.Out(getInfoStream());

            getInfoStream() << "\nDenseMatrix3x3 = TensorVariable3U";
            DenseMatrix3x3 = TensorVariable3U;
            DenseMatrix3x3.Out(getInfoStream());

            for ( int i = 0; i < 3; i++ )
            {
                for ( int j = 0; j < 3; j++ )
                {
                    if ( i==j ) _equal( DenseMatrix3x3(i,j), 3., 1E-6 );
                    else        _equal( DenseMatrix3x3(i,j), 2., 1E-6 );
                }
            }

            //Testing operator*= ( TensorVariable )
                getInfoStream() << "\nTesting operator*= ( TensorVariable )" << endl;
                getInfoStream() << "====================================" << endl;

                //Operator *= ( double64 )
                getInfoStream() << "\nOperator *= ( double64 )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3 *= ( 5 )";
                DenseMatrix3x3 *= ( 5 );
                DenseMatrix3x3.Out(getInfoStream());

                _equal( DenseMatrix3x3(0,0), 25., 1E-6 );

                //Operator *= ( TensorVariable 2U )
                getInfoStream() << "\nOperator *= ( TensorVariable 2U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n2U tensor variable TensorVariable2U";
                TensorVariable2U.Identity();
                TensorVariable2U += 1.;
                TensorVariable2U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3 *= ( TensorVariable2U )";
                DenseMatrix3x3 *= ( TensorVariable2U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 2; i++ )
                {
                    for ( int j = 0; j < 2; j++ )
                    {
                        _equal( DenseMatrix3x3(i,j), 15., 1E-6 );
                    }
                }


                //Operator*= ( TensorVariable 3U )
                getInfoStream() << "\nOperator*= ( TensorVariable 3U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n3U tensor variable TensorVariable3U";
                TensorVariable3U.Identity();
                TensorVariable3U += 2.;
                TensorVariable3U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3 *= ( TensorVariable3U )";
                DenseMatrix3x3 *= ( TensorVariable3U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 3; i++ )
                {
                    for ( int j = 0; j < 3; j++ )
                    {
                        _equal( DenseMatrix3x3(i,j), 35., 1E-6 );
                    }
                }

            //Testing operator*= ( VectorVariable )
                getInfoStream() << "\nTesting operator*= ( VectorVariable )" << endl;
                getInfoStream() << "====================================" << endl;

                //Operator *= ( VectorVariable 2U )
                getInfoStream() << "\nOperator *= ( VectorVariable 2U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n2U Vector variable VectorVariable2U";
                VectorVariable2U = 5;
                VectorVariable2U += 1.;
                VectorVariable2U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix2x2 *= ( VectorVariable2U )";
                DenseMatrix3x3 *= ( VectorVariable2U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(i,0), 60., 1E-6 );

                }


                //Operator*= ( VectorVariable 3U )
                getInfoStream() << "\nOperator*= ( VectorVariable 3U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n3U Vector variable VectorVariable3U";
                VectorVariable3U = 5;
                VectorVariable3U += 2.;
                VectorVariable3U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3 *= ( VectorVariable3U )";
                DenseMatrix3x3 *= ( VectorVariable3U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(i,0), 105., 1E-6 );
                }


            //AssignToDiagonal ( TensorVariable )
                getInfoStream() << "\nTesting AssignToDiagonal ( VectorVariable )" << endl;
                getInfoStream() << "===========================================" << endl;

                //AssignToDiagonal ( diag_elmts, ScalarVariable )
                getInfoStream() << "\nAssignToDiagonal ( diag_elmts, ScalarVariable )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3.AssignToDiagonal ( 2, 6. )";
                DenseMatrix3x3.AssignToDiagonal(2, 6.);
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0, j = 0; i < 2; i++, j++ )
                {
                   _equal( DenseMatrix3x3(i,j), 6., 1E-6 );
                }

                //AssignToDiagonal ( VectorVariable1U )
                getInfoStream() << "\nAssignToDiagonal ( VectorVariable 1U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n1U test matrix DenseMatrix1x1";
                DenseMatrix3x3.Resize(1,1);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n1U Vector variable VectorVariable1U";
                VectorVariable1U = 5;
                VectorVariable1U += 6.;
                VectorVariable1U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix1x1.AssignToDiagonal ( VectorVariable1U )";
                DenseMatrix3x3.AssignToDiagonal ( VectorVariable1U );
                DenseMatrix3x3.Out(getInfoStream());

                _equal( DenseMatrix3x3(0,0), 11., 1E-6 );

                //AssignToDiagonal ( VectorVariable2U )
                getInfoStream() << "\nAssignToDiagonal ( VectorVariable 2U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n2U Vector variable VectorVariable2U";
                VectorVariable2U = 9;
                VectorVariable2U += 1.;
                VectorVariable2U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix2x2.AssignToDiagonal ( VectorVariable2U )";
                DenseMatrix3x3.AssignToDiagonal ( VectorVariable2U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0, j = 0; i < 2; i++, j++ )
                {
                   _equal( DenseMatrix3x3(i,j), 10., 1E-6 );
                }


                //AssignToDiagonal ( VectorVariable3U )
                getInfoStream() << "\nAssignToDiagonal ( VectorVariable3U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n3U Vector variable VectorVariable3U";
                VectorVariable3U = 8;
                VectorVariable3U += 2.;
                VectorVariable3U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3.AssignToDiagonal ( VectorVariable3U )";
                DenseMatrix3x3.AssignToDiagonal ( VectorVariable3U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0, j = 0 ; i < 3; i++, j++ )
                {
                    _equal( DenseMatrix3x3(i,j), 10., 1E-6 );
                }


                //AssignToDiagonal ( size_t diag elmnts, double64 )
                getInfoStream() << "\nAssignToDiagonal ( size_t diag elmnts, double64 )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(5.);
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3.AssignToDiagonal ( size_t diag elmnts, double64 )";
                DenseMatrix3x3.AssignToDiagonal(3, 10.);
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0, j = 0 ; i < 3; i++, j++ )
                {
                    _equal( DenseMatrix3x3(i,j), 10., 1E-6 );
                }


            //RowCondenseTo( vector<double64> )
                getInfoStream() << "\nTesting RowCondenseTo( vector<double64> )" << endl;
                getInfoStream() << "===========================================" << endl;
                //RowCondenseTo ( vector<double64> )
                getInfoStream() << "\nRowCondenseTo ( vector<double64> )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Fill(9.);
                DenseMatrix3x3.Out(getInfoStream());

                vector<double64>::const_iterator itx(x.begin());

                getInfoStream() << "Test vector x";

                for (vector<double64>::const_iterator
                     itx = x.begin(); itx != x.end(); itx++)
                {
                    getInfoStream() << "\n" << *itx;
                }

                getInfoStream() << "\n\nDenseMatrix3x3.RowCondenseTo( x )";
                DenseMatrix3x3.RowCondenseTo( x );

                for (vector<double64>::const_iterator
                     itx = x.begin(); itx != x.end(); itx++)
                {
                    getInfoStream() << "\n" << *itx ;
                    _equal( *itx, 27., 1E-6 );
                }

                getInfoStream() << endl;


            //ExportTo ( TensorVariable )
                getInfoStream() << "\nTesting ExportTo ( TensorVariable )" << endl;
                getInfoStream() << "===========================================" << endl;

                //ExportTo ( TensorVariable 2U )
                getInfoStream() << "\nExportTo ( TensorVariable 2U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n2U tensor variable TensorVariable2U";
                TensorVariable2U.Zero();
                TensorVariable2U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix2x2 *= ( TensorVariable2U )";
                DenseMatrix3x3.ExportTo(TensorVariable2U);
                TensorVariable3U.Out(getInfoStream());

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(i,i), 1., 1E-6 );
                }

                //ExportTo ( TensorVariable 3U )
                getInfoStream() << "\n\nExportTo ( TensorVariable 3U )" << endl;
                getInfoStream() << "-------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "\n3U tensor variable TensorVariable3U";
                TensorVariable3U.Zero();
                TensorVariable3U.Out(getInfoStream());

                getInfoStream() << "\nDenseMatrix3x3 *= ( TensorVariable3U )";
                DenseMatrix3x3.ExportTo(TensorVariable3U);
                TensorVariable3U.Out(getInfoStream());

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(i,i), 1., 1E-6 );
                }

            //Assign ( size_t , Point )
                getInfoStream() << "\n\nTesting Assign ( size_t , Point )" << endl;
                getInfoStream() << "===========================================" << endl;

                //AssignRow ( size_t i, Point<1U> )
                getInfoStream() << "\nAssignRow ( size_t i, Point<1U> )" << endl;
                getInfoStream() << "---------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "Point Point<1U>";
                Point1U= 13.;
                Point1U.Out(getInfoStream());

                getInfoStream() << "\nAssignRow( 1, Point1U )";
                DenseMatrix3x3.AssignRow( 1, Point1U );
                DenseMatrix3x3.Out(getInfoStream());

                _equal( DenseMatrix3x3(1,0), 13., 1E-6 );

                //AssignRow ( size_t i, Point<2U> )
                getInfoStream() << "AssignRow ( size_t i, Point<2U> )" << endl;
                getInfoStream() << "---------------------------------" << endl;

                getInfoStream() << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Zero();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "Point Point<2U>";
                Point2U = 11.;
                Point2U.Out(getInfoStream());

                getInfoStream() << "\nAssignRow( 1, Point2U )";
                DenseMatrix3x3.AssignRow( 1, Point2U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(1,i), 11., 1E-6 );
                }

                //AssignRow ( size_t i, Point<3U> )
                getInfoStream() << "\nAssignRow ( size_t i, Point<3U> )" << endl;
                getInfoStream() << "---------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "Point Point<3U>";
                Point3U.Set (13., 12., 11.);
                Point3U.Out(getInfoStream());

                getInfoStream() << "\nAssignRow( 1, Point2U )";
                DenseMatrix3x3.AssignRow( 1, Point3U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(1,i), 13.-i, 1E-6 );
                }

                //AssignCol ( size_t i, Point<1U> )
                getInfoStream() << "\nAssignCol ( size_t i, Point<1U> )" << endl;
                getInfoStream() << "---------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "Point Point<1U>";
                Point1U= 13.;
                Point1U.Out(getInfoStream());

                getInfoStream() << "\nAssignCol( 1, Point1U )";
                DenseMatrix3x3.AssignCol( 1, Point1U );
                DenseMatrix3x3.Out(getInfoStream());

                _equal( DenseMatrix3x3(0,1), 13., 1E-6 );

                //AssignCol ( size_t i, Point<2U> )
                getInfoStream() << "\nAssignCol ( size_t i, Point<2U> )" << endl;
                getInfoStream() << "---------------------------------" << endl;

                getInfoStream() << "\n2U test matrix DenseMatrix2x2";
                DenseMatrix3x3.Resize(2,2);
                DenseMatrix3x3.Zero();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "Point Point<2U>";
                Point2U = 11.;
                Point2U.Out(getInfoStream());

                getInfoStream() << "\nAssignCol( 1, Point2U )";
                DenseMatrix3x3.AssignCol( 1, Point2U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 2; i++ )
                {
                    _equal( DenseMatrix3x3(i,1), 11., 1E-6 );
                }

                //AssignCol ( size_t i, Point<3U> )
                getInfoStream() << "\nAssignCol ( size_t i, Point<3U> )" << endl;
                getInfoStream() << "---------------------------------" << endl;

                getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                DenseMatrix3x3.Resize(3,3);
                DenseMatrix3x3.Identity();
                DenseMatrix3x3.Out(getInfoStream());

                getInfoStream() << "Point Point<3U>";
                Point3U.Set (13., 12., 11.);
                Point3U.Out(getInfoStream());

                getInfoStream() << "\nAssignCol( 1, Point2U )";
                DenseMatrix3x3.AssignCol( 2, Point3U );
                DenseMatrix3x3.Out(getInfoStream());

                for ( int i = 0; i < 3; i++ )
                {
                    _equal( DenseMatrix3x3(i,2), 13.-i, 1E-6 );
                }


                //Testing operator*= ( Point<1U> )
                    getInfoStream() << "\nTesting operator*= ( Point<1U> )" << endl;
                    getInfoStream() << "====================================" << endl;

                    //Operator *= ( Point<1U> )
                    getInfoStream() << "\nOperator *= ( Point<1U> )" << endl;
                    getInfoStream() << "-------------------------------" << endl;

                    getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                    DenseMatrix3x3.Resize(1,1);
                    DenseMatrix3x3.Fill(5.);
                    DenseMatrix3x3.Out(getInfoStream());

                    getInfoStream() << "\nPoint<1U> Point1U";
                    Point1U.Out(getInfoStream());

                    getInfoStream() << "\nDenseMatrix3x3 *= ( Point1U )";
                    DenseMatrix3x3 *= ( Point1U );
                    DenseMatrix3x3.Out(getInfoStream());

                    _equal( DenseMatrix3x3(0,0), 65., 1E-6 );

                    //Operator *= ( Point<2U> )
                    getInfoStream() << "\nOperator *= ( Point<2U> )" << endl;
                    getInfoStream() << "-------------------------------" << endl;

                    getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                    DenseMatrix3x3.Resize(2,2);
                    DenseMatrix3x3.Fill(5.);
                    DenseMatrix3x3.Out(getInfoStream());

                    getInfoStream() << "\nPoint<2U> Point2U";
                    Point2U.Out(getInfoStream());

                    getInfoStream() << "\nDenseMatrix3x3 *= ( Point2U )";
                    DenseMatrix3x3 *= ( Point2U );
                    DenseMatrix3x3.Out(getInfoStream());

                    for ( int i = 0; i < 2; i++ )
                        _equal( DenseMatrix3x3(i,0), 110., 1E-6 );


                    //Operator*= ( Point<3U> )
                    getInfoStream() << "\nOperator*= ( Point<3U> )" << endl;
                    getInfoStream() << "-------------------------------" << endl;

                    getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                    DenseMatrix3x3.Resize(3,3);
                    DenseMatrix3x3.Fill(5.);
                    DenseMatrix3x3.Out(getInfoStream());

                    getInfoStream() << "\nPoint<3U> Point3U";
                    Point3U.Out(getInfoStream());

                    getInfoStream() << "\nDenseMatrix3x3 *= ( Point3U )";
                    DenseMatrix3x3 *= ( Point3U );
                    DenseMatrix3x3.Out(getInfoStream());

                    for ( int i = 0; i < 3; i++ )
                        _equal( DenseMatrix3x3(i,0), 180., 1E-6 );


                //Testing L1 matrix norm()
                    getInfoStream() << "\nTesting NormL1()" << endl;
                    getInfoStream() << "================" << endl;

                        getInfoStream() << "\n3U test matrix A";
                        A(2,2) = -99.;
                        A(2,1) = -199.;
                        A.Out(getInfoStream());
                        double64 normL1 = A.NormL1();

                        getInfoStream() << "\nThe maximum absolute column sum norm of test matrix A is " << normL1 << endl;
                        _equal( normL1, 201., 1E-6 );


                //Testing L infinity matrix norm()
                    getInfoStream() << "\nTesting NormL_Infinity()" << endl;
                    getInfoStream() << "================" << endl;

                        getInfoStream() << "\n3U test matrix A";
                        A(2,2) = -99.;
                        A(2,1) = -199.;
                        A.Out(getInfoStream());
                        double64 normInf = A.NormL_Infinity();

                        getInfoStream() << "\nThe maximum absolute row sum norm of test matrix A is " << normInf << endl;
                        _equal( normInf, 303., 1E-6 );

/*
                //Testing In()
                    getInfoStream() << "\nTesting In()" << endl;
                    getInfoStream() << "============" << endl;


                        getInfoStream() << "\n3U test matrix DenseMatrix3x3";
                        DenseMatrix3x3.Resize(3,3);
                        int rows = 2;
                        int cols = 2;
                        DenseMatrix3x3.In();
*/

    }


} // end csmp

