// unit test
// Created by : Siroos, Georg and Alina
// Date of Modification: 21.11.2010

#include "SparseMatrix_Test.h"
#include "PL_Utilities.h"
#include "DenseMatrix.h"
#include <iomanip>
#include <iostream>

using namespace std;

namespace csmp
{
    SparseMatrix_Test::SparseMatrix_Test()
        :
        A(4),
        ResA(4),
        B(6),
        BB(6),
        BBB(6),
        ResAB(4),
        C(4),
        ResC(4),
        D(4),
        ResD(4),
        E(4),
        F(4),
        ResF(4),
        G(4),
        ResG(4),
        H(4),
        SparseMatrix2x2(2),
        SparseMatrix3x3(2),
        fTolerance(1.e-6)
    {
    }

    SparseMatrix_Test::~SparseMatrix_Test()
    {
    }


    void SparseMatrix_Test::run()
    {

    cout << "\n====================";
    cout << "\nTesting SparseMatrix" << endl;
    cout << "====================" << endl;

    //Testing Assign function
        cout << "\nTesting Assign function" << endl;
        cout << "=======================" << endl;

        A.Assign( 0, 0, 1. );
        A.Assign( 0, 1, 2. );
        A.Assign( 0, 2, 2. );
        A.Assign( 0, 3, 1. );

        A.Assign( 1, 0, 3. );
        A.Assign( 1, 1, 1. );
        A.Assign( 1, 2, 5. );
        A.Assign( 1, 3, 2. );

        A.Assign( 2, 0, 1. );
        A.Assign( 2, 1, 5. );
        A.Assign( 2, 2, 2. );
        A.Assign( 2, 3, 4. );

        A.Assign( 3, 0, 4. );
        A.Assign( 3, 1, 1. );
        A.Assign( 3, 2, 3. );
        A.Assign( 3, 3, 5. );

        cout << "SparsityMatrix A = " << endl;
        A.Out();

        _test( A( 0,1 ) == 2. );
        _test( A( 3,3 ) == 5. );

        DenseMatrix<4> DenseMatrix2x2;
        DenseMatrix2x2( 0, 0 ) = 2;
        DenseMatrix2x2( 0, 1 ) = 3;
        DenseMatrix2x2( 1, 0 ) = 4;
        DenseMatrix2x2( 1, 1 ) = 5;

        Matrix Matrix2x2Val(4,4);
        Matrix2x2Val( 0, 0 ) = 40.;
        Matrix2x2Val( 0, 1 ) = 30.;
        Matrix2x2Val( 1, 0 ) = 20.;
        Matrix2x2Val( 1, 1 ) = 10.;

        B.Assign( DenseMatrix2x2, Matrix2x2Val );
        cout << "SparsityMatrix B = " << endl;
        B.Out();
        BB = B;
        D=A;

        _test( B( 2, 2 ) == 40. );
        _test( B( 3, 3 ) == 30. );
        _test( B( 4, 4 ) == 20. );
        _test( B( 5, 5 ) == 10. );


    //Testing MultiplyWith function
        cout << "\nTesting MultiplyWith function" << endl;
        cout << "=============================" << endl;

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
        A.MultiplyWith( x, y );

        vector<double64>::const_iterator it1(y.begin());

        cout << "\n\ny = ";

        for (vector<double64>::const_iterator
             it1 = y.begin(); it1 != y.end(); it1++ )
        {
            cout << *it1 << setw(5);
        }

        sol_y.push_back(6.5);
        sol_y.push_back(8.);
        sol_y.push_back(17.);
        sol_y.push_back(13.5);

        vector<double64>::const_iterator itsol_y(sol_y.begin());

        for (vector<double64>::const_iterator
             itsol_y = sol_y.begin(); itsol_y != sol_y.end(); itsol_y++, it1++ )
            {
                //cout << "\n" << *itsol_y << setw(7) << *it1;
                _equal( *itsol_y, *it1, 1E-6);
            }


    //Testing Operator= ( SparseMatrix )
    cout << "\n Testing Operator= ( SparseMatrix )" << endl;
    cout << "==========================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();

        cout << "SparseMatrix B = SparseMatrix A";
        B = A;
        B.Out();

        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                _test( B( i, j ) == A( i, j ) );

            }
        }


    //Testing MultiplyEntryWith function
        cout << "\nTesting MultiplyEntryWith function" << endl;
        cout << "==================================" << endl;

        cout << "\nSparseMatrix B ( 0, 0 ) * 16. = 1. * 16.";
        B.MultiplyEntryWith( 0, 0, 16. );
        B.Out();
        B.RemoveEntry( 0, 1 );

        cout << "\nSparseMatrix B ( 0, 1 ) * 99. = 0. * 99.";
        B.MultiplyEntryWith( 0, 1, 99. );
        B.Out();

        _equal( B( 0, 0 ) , 16. ,fTolerance );
        _equal( B( 0, 1 ) , 0. ,fTolerance );


    //Testing Add function
        cout << "\nTesting Add function" << endl;
        cout << "====================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();
        B.Add( 0, 0, 16. );
        B.Add( 0, 1, 15. );
        cout << "SparseMatrix B ( 0, 0 ) + 16. = 32.";
        cout << "\nSparseMatrix B ( 0, 1 ) + 15. = 15.";
        B.Out();

        _test( B( 0, 0 ) == 32. );
        _test( B( 0, 1 ) == 15. );


    //Testing Rows function
        cout << "\nTesting Rows function" << endl;
        cout << "=====================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();

        _test( B.Rows() == 4. );


    //Testing Cols function
        cout << "\nTesting Cols function" << endl;
        cout << "=====================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();

        _test( B.Cols() == 4. );


    //Testing ZeroRow function
        cout << "\nTesting ZeroRow function" << endl;
        cout << "========================" << endl;

        cout << "\nSparseMatrix B before ZeroRow(3)";
        B.Out();
        B.ZeroRow(3);
        cout << "\nSparseMatrix B after ZeroRow(3)";
        B.Out();

        _test( B.Entries() == 12 );

        for ( int i = 3; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _test( B( i, j ) == 0. );
                cout << "B( "<< i <<", " << j << " ) = " << B( i, j ) << endl;
            }
        }

    //Testing ZeroColumn function
        cout << "\nTesting ZeroColumn function" << endl;
        cout << "===========================" << endl;

        cout << "\nSparseMatrix B before ZeroColumn( 3 )";
        B.Out();
        B.ZeroColumn(3);
        cout << "\nSparseMatrix B after ZeroColumn( 3 )";
        B.Out();

        _test( B.Entries() == 9 );

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 3; j < 4; j++ )
            {
                _test( B( i, j ) == 0. );
                cout << "B( "<< i <<", " << j << " ) = " << B( i, j ) << endl;
            }
        }
        cout << endl;

    //Testing RemoveEntry function
        cout << "\nTesting RemoveEntry function" << endl;
        cout << "============================" << endl;

        cout << "\nSparseMatrix B entries before RemoveEntry( 0, 0) = ";
        cout << B.Entries() << endl;
        B.RemoveEntry( 0, 0 );
        cout << "\nSparseMatrix B entries after RemoveEntry( 0, 0 ) = ";
        cout << B.Entries() << endl;
        cout << "\nEntry B( 0, 0 ) = " <<  B( 0, 0 ) << endl;

        _test( B.Entries() == 8. );
        _test( B( 0, 0 ) == 0. );


    //Testing Zero function
        cout << "\nTesting Zero function" << endl;
        cout << "=====================" << endl;

        C = B;
        cout << "\nSparseMatrix B before Zero()";
        B.Out();
        B.Zero();
        cout << "SparseMatrix B after Zero()";
        B.Out();

        _test( B.Entries() == 0. );


    //Testing Erase function
        cout << "\nTesting Erase function" << endl;
        cout << "======================" << endl;

        B = C;
        cout << "\nSparseMatrix B before Erase()";
        B.Out();
        B.Erase();
        cout << "SparseMatrix B after Erase()";
        B.Out();

        _test( B.Entries() == 0. );


    //Testing ColumnIndices function
        cout << "\nTesting ColumnIndices function" << endl;
        cout << "==============================" << endl;

        cout << "\nSparseMatrix B";
        BB.Assign( 3, 4, 5. );
        BBB = BB;
        BB.Out();
        BB.ColumnIndices( 3, sizetVector );
        cout << "Vector sizetVector resulting from B.ColumnIndices( 3, sizetVector ) = (" << sizetVector[0] << "," << sizetVector[1] << ")" << endl;

        _test( sizetVector.size() == 2. );
        _test( sizetVector[0] == 3. );
        _test( sizetVector[1] == 4. );


    //Testing Resize function
        cout << "\nTesting Resize function" << endl;
        cout << "=======================" << endl;

        cout << "\nSparseMatrix A before Resize()";
        A.Out();

        cout << "SparseMatrix A after Resize(2)";
        A.Resize(2);
        A.Out();

        _test( A.Rows() == 2. );
        _test( A.Cols() == 2. );

        cout << "SparseMatrix A after Resize(5)";
        A.Resize(5);
        A.Out();

       _test( A.Rows() == 5. );
       _test( A.Cols() == 5. );


   //Testing At function
       cout << "\nTesting At function" << endl;
       cout << "===================" << endl;

       cout << "\nSparseMatrix B";
       BB.Out();

       cout << "SparseMatrix B.At( 3, 4 ) = " << BB.At( 3, 4 ) << endl;

      _test( BB.At( 3, 4 ) == 5. );


   //Testing Operator() ( size_t, size_t )
      cout << "\nTesting Operator() ( size_t, size_t )" << endl;
      cout << "=====================================" << endl;

      cout << "\nSparseMatrix B";
      BB.Out();

      cout << "SparseMatrix B( 3, 4 ) = " << BB(3,4) << endl;

     _test( BB( 3, 4 ) == 5. );


      //Testing Symmetric function
         cout << "\nTesting Symmetric function" << endl;
         cout << "==========================" << endl;

         BB.Assign( 0, 0, 10. );
         BB.Assign( 1, 1, 20. );
         cout << "\nSparseMatrix B";
         BB.Out();

         cout << "SparseMatrix B is symmetric = " << BB.Symmetric();

         _test( BB.Symmetric() == 0 );

         BB.RemoveEntry( 3, 4 );
         cout << "\nSparseMatrix B entries after RemoveEntry( 3, 4 )";
         BB.Out();

         cout << "SparseMatrix B is symmetric = " << BB.Symmetric();

         _test( BB.Symmetric() == 1 );


     //Testing ZeroesInDiagonal function
         cout << "\nTesting ZeroesInDiagonal function" << endl;
         cout << "=================================" << endl;

         cout << "\nSparseMatrix B";
         BB.Out();
         cout << "Does SparseMatrix B has zeros in its diagonal? " << BB.ZeroesInDiagonal() << endl;
         _test( BB.ZeroesInDiagonal() == 0 );

         BB.RemoveEntry( 0, 0 );
         //BB.Assign( 0, 0, 0. );
         cout << "\nSparseMatrix B entries after Assign( 0, 0, 0. )";
         BB.Out();
         cout << "Does SparseMatrix B has zeros in its diagonal? " << BB.ZeroesInDiagonal() << endl;
         _test( BB.ZeroesInDiagonal() == 1 );

         BB.RemoveEntry( 0, 0 );
         cout << "\nSparseMatrix B entries after RemoveEntry( 0, 0 )";
         BB.Out();
         cout << "Does SparseMatrix B has zeros in its diagonal? " << BB.ZeroesInDiagonal();
         _test( BB.ZeroesInDiagonal() == 1 );


     //Testing DiagonallyPositive function
         cout << "\nTesting DiagonallyPositive function" << endl;
         cout << "===================================" << endl;

         cout << "\nSparseMatrix B";
         BB.Out();
         cout << "Is SparseMatrix B diagonally positive? " << BB.DiagonallyPositive() << endl;
         _test( BB.DiagonallyPositive() == 1 );

         BB.Assign( 0, 0, -5. );
         BB.Assign( 2, 2, -10. );
         cout << "\nSparseMatrix B entries after Assign( 0, 0, -5. ) and Assign( 2, 2, -10. )";
         BB.Out();
         cout << "Is SparseMatrix B diagonally positive? " << BB.DiagonallyPositive() << endl;
         _test( BB.DiagonallyPositive() == 0 );


     //Testing RecountEntries function
         cout << "\nTesting RecountEntries function" << endl;
         cout << "===============================" << endl;

         cout << "\nSparseMatrix B";
         BB.Out();
         cout << "Recounted Entries of SparseMatrix B = " << BB.RecountEntries() << endl;
         _test( BB.RecountEntries() == BB.Entries() );

         BB.Assign( 0, 1, 11. );
         cout << "Recounted Entries of SparseMatrix B = " << BB.RecountEntries() << endl;
         _test( BB.RecountEntries() == BB.Entries() );


     //Testing RemoveHalo function
         cout << "\nTesting RemoveHalo function" << endl;
         cout << "===========================" << endl;

         cout << "\nSparseMatrix B";
         BBB.Out();

         BBB.RemoveHalo(2);
         cout << "SparseMatrix B entries after RemoveHalo = " << BBB.Entries() << endl;
         BBB.Out();
         _test( BBB.Entries() == 3 );


     //Testing SparsityPattern function
         cout << "\nTesting SparsityPattern function" << endl;
         cout << "================================" << endl;

         cout << "\nSparseMatrix A";
         BBB.Out();

         BBB.SparsityPattern( "Sparsity" );
         string datafile = "cspline_test_data";
         cout << "Data file Sparsity has been created" << endl;

     //Testing InfinityNorm function
         cout << "\nTesting InfinityNorm function" << endl;
         cout << "================================" << endl;
         cout << "\nSparseMatrix D";
         D.MultiplyEntryWith(0,1,-1.);
         D.MultiplyEntryWith(1,2,-1.);
         D.MultiplyEntryWith(2,3,-1.);
         D.MultiplyEntryWith(3,0,-1);
         D.Out();
         _test( D.InfinityNorm() == 13 );
         cout << "SparseMatrix D infinity norm = " << D.InfinityNorm() << endl;

  }


} // end csmp

