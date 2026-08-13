// unit test
// Created by : Siroos, Georg and Alina
// Date of Modification: 21.11.2010

#include "SparseMatrix_Test.h"
#include "DenseMatrix.h"
#include "Matrix.h"
#include "compareFloats.h"

using namespace std;

namespace csmp {

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
    SparseMatrix3x3(2)
{
}


// checked: SKM 8/5/2022
void SparseMatrix_Test::run()
 {
    if ( verbose_ ) {
        cout << "\n====================";
        cout << "\nTesting SparseMatrix" << endl;
        cout << "====================" << endl;

        //Testing Assign function
        cout << "\nTesting Assign function" << endl;
        cout << "=======================" << endl;
     }
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
   
    if ( verbose_ ) {
        cout << "SparsityMatrix A = " << endl;
        A.Out();
      }
        _test( essentiallyEqual( A( 0,1 ), 2.) );
        _test( essentiallyEqual( A( 3,3 ), 5.) );

        DenseMatrix<4> DenseMatrix2x2;
        DenseMatrix2x2( 0, 0 ) = 2.;
        DenseMatrix2x2( 0, 1 ) = 3.;
        DenseMatrix2x2( 1, 0 ) = 4.;
        DenseMatrix2x2( 1, 1 ) = 5.;

        Matrix Matrix2x2Val(4,4);
        Matrix2x2Val( 0, 0 ) = 40.;
        Matrix2x2Val( 0, 1 ) = 30.;
        Matrix2x2Val( 1, 0 ) = 20.;
        Matrix2x2Val( 1, 1 ) = 10.;

        B.Assign( DenseMatrix2x2, Matrix2x2Val );
        if ( verbose_ ) {
            cout << "SparsityMatrix B = " << endl;
            B.Out();
          }
        BB = B;
        D  = A;

        _test( essentiallyEqual( B( 2, 2 ), 40.) );
        _test( essentiallyEqual( B( 3, 3 ), 30.) );
        _test( essentiallyEqual( B( 4, 4 ), 20.) );
        _test( essentiallyEqual( B( 5, 5 ), 10.) );


    //Testing MultiplyWith function
    if ( verbose_ ) {
        cout << "\nTesting MultiplyWith function" << endl;
        cout << "=============================" << endl;
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
        A.MultiplyWith( x, y );

        vector<double>::const_iterator it1(y.begin());

        if ( verbose_ ) cout << "\n\ny = ";

        for ( it1 = y.begin(); it1 != y.end(); it1++ )
        {
            if ( verbose_ ) cout << *it1 << setw(5);
        }

        sol_y.push_back(6.5);
        sol_y.push_back(8.);
        sol_y.push_back(17.);
        sol_y.push_back(13.5);

        vector<double>::const_iterator itsol_y(sol_y.begin());

        for ( itsol_y = sol_y.begin(), it1=y.begin(); itsol_y != sol_y.end(); itsol_y++, it1++ )
            {
                //cout << "\n" << *itsol_y << setw(7) << *it1;
                _equal( *itsol_y, *it1, 1E-6);
            }


    //Testing Operator= ( SparseMatrix )
    if ( verbose_ ) {
        cout << "\n Testing Operator= ( SparseMatrix )" << endl;
        cout << "==========================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();

        cout << "SparseMatrix B = SparseMatrix A";
     }
        B = A;
        if ( verbose_ ) B.Out();

        for ( int i = 0; i < 3; i++ )
        {
            for ( int j = 0; j < 3; j++ )
            {
                _test( essentiallyEqual( B( i, j ), A( i, j ) ) );

            }
        }


    //Testing MultiplyEntryWith function
    if ( verbose_ ) {
        cout << "\nTesting MultiplyEntryWith function" << endl;
        cout << "==================================" << endl;

        cout << "\nSparseMatrix B ( 0, 0 ) * 16. = 1. * 16.";
     }
        B.MultiplyEntryWith( 0, 0, 16. );
        if ( verbose_ ) B.Out();
        B.RemoveEntry( 0, 1 );

        if ( verbose_ ) cout << "\nSparseMatrix B ( 0, 1 ) * 99. = 0. * 99.";
        B.MultiplyEntryWith( 0, 1, 99. );
        if ( verbose_ ) B.Out();

        _equal( B( 0, 0 ) , 16. ,fTolerance );
        _equal( B( 0, 1 ) , 0. ,fTolerance );


    //Testing Add function
    if ( verbose_ ) {
        cout << "\nTesting Add function" << endl;
        cout << "====================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();
      }
        B.Add( 0, 0, 16. );
        B.Add( 0, 1, 15. );
      if ( verbose_ ) {
          cout << "SparseMatrix B ( 0, 0 ) + 16. = 32.";
          cout << "\nSparseMatrix B ( 0, 1 ) + 15. = 15.";
          B.Out();
        }
   
        _test( essentiallyEqual( B( 0, 0 ), 32.) );
        _test( essentiallyEqual( B( 0, 1 ), 15.) );


    //Testing Rows function
    if ( verbose_ ) {
        cout << "\nTesting Rows function" << endl;
        cout << "=====================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();
      }
        _test( B.Rows() == 4. );


    //Testing Cols function
    if ( verbose_ ) {
        cout << "\nTesting Cols function" << endl;
        cout << "=====================" << endl;

        cout << "\nSparseMatrix B";
        B.Out();
      }
        _test( B.Cols() == 4. );


    //Testing ZeroRow function
    if ( verbose_ ) {
        cout << "\nTesting ZeroRow function" << endl;
        cout << "========================" << endl;

        cout << "\nSparseMatrix B before ZeroRow(3)";
        B.Out();
      }
        B.ZeroRow(3);
    if ( verbose_ ) {
        cout << "\nSparseMatrix B after ZeroRow(3)";
        B.Out();
      }
        _test( B.Entries() == 12 );

        for ( int i = 3; i < 4; i++ )
        {
            for ( int j = 0; j < 4; j++ )
            {
                _test( B( i, j ) == 0. );
                if ( verbose_ ) cout << "B( "<< i <<", " << j << " ) = " << B( i, j ) << endl;
            }
        }

    //Testing ZeroColumn function
    if ( verbose_ ) {
        cout << "\nTesting ZeroColumn function" << endl;
        cout << "===========================" << endl;

        cout << "\nSparseMatrix B before ZeroColumn( 3 )";
        B.Out();
      }
        B.ZeroColumn(3);
    if ( verbose_ ) {
        cout << "\nSparseMatrix B after ZeroColumn( 3 )";
        B.Out();
      }

        _test( B.Entries() == 9 );

        for ( int i = 0; i < 4; i++ )
        {
            for ( int j = 3; j < 4; j++ )
            {
                _test( essentiallyEqual( B( i, j ), 0.) );
                if ( verbose_ ) cout << "B( "<< i <<", " << j << " ) = " << B( i, j ) << endl;
            }
        }
     if ( verbose_ )  cout << endl;

    //Testing RemoveEntry function
    if ( verbose_ ) {
        cout << "\nTesting RemoveEntry function" << endl;
        cout << "============================" << endl;

        cout << "\nSparseMatrix B entries before RemoveEntry( 0, 0) = ";
        cout << B.Entries() << endl;
      }
        B.RemoveEntry( 0, 0 );
    if ( verbose_ ) {
        cout << "\nSparseMatrix B entries after RemoveEntry( 0, 0 ) = ";
        cout << B.Entries() << endl;
        cout << "\nEntry B( 0, 0 ) = " <<  B( 0, 0 ) << endl;
      }
        _test( B.Entries() == 8. );
        _test( essentiallyEqual( B( 0, 0 ), 0.) );


    //Testing Zero function
    if ( verbose_ ) {
        cout << "\nTesting Zero function" << endl;
        cout << "=====================" << endl;
      }
        C = B;
    if ( verbose_ ) {
        cout << "\nSparseMatrix B before Zero()";
        B.Out();
      }
        B.Zero();
    if ( verbose_ ) {
        cout << "SparseMatrix B after Zero()";
        B.Out();
      }
        _test( B.Entries() == 0. );


    //Testing Erase function
    if ( verbose_ ) {
        cout << "\nTesting Erase function" << endl;
        cout << "======================" << endl;
      }
        B = C;
    if ( verbose_ ) {
        cout << "\nSparseMatrix B before Erase()";
        B.Out();
      }
        B.Erase();
    if ( verbose_ ) {
        cout << "SparseMatrix B after Erase()";
        B.Out();
      }
        _test( B.Entries() == 0. );


    //Testing ColumnIndices function
    if ( verbose_ ) {
        cout << "\nTesting ColumnIndices function" << endl;
        cout << "==============================" << endl;

        cout << "\nSparseMatrix B";
      }
        BB.Assign( 3, 4, 5. );
        BBB = BB;
        if ( verbose_ ) BB.Out();
        BB.ColumnIndices( 3, sizetVector );
        if ( verbose_ )
          cout << "Vector sizetVector resulting from B.ColumnIndices( 3, sizetVector ) = (" << sizetVector[0] << "," << sizetVector[1] << ")" << endl;

        _test( sizetVector.size() == 2. );
        _test( sizetVector[0] == 3. );
        _test( sizetVector[1] == 4. );


    //Testing Resize function
    if ( verbose_ ) {
        cout << "\nTesting Resize function" << endl;
        cout << "=======================" << endl;

        cout << "\nSparseMatrix A before Resize()";
        A.Out();

        cout << "SparseMatrix A after Resize(2)";
      }
        A.Resize(2);
        if ( verbose_ ) A.Out();

        _test( A.Rows() == 2. );
        _test( A.Cols() == 2. );

        if ( verbose_ ) cout << "SparseMatrix A after Resize(5)";
        A.Resize(5);
        A.Out();

       _test( A.Rows() == 5. );
       _test( A.Cols() == 5. );


   //Testing At function
   if ( verbose_ ) {
       cout << "\nTesting At function" << endl;
       cout << "===================" << endl;

       cout << "\nSparseMatrix B";
       BB.Out();

       cout << "SparseMatrix B.At( 3, 4 ) = " << BB.At( 3, 4 ) << endl;
     }
      _test( essentiallyEqual( BB.At( 3, 4 ), 5.) );


   //Testing Operator() ( size_t, size_t )
   if ( verbose_ ) {
      cout << "\nTesting Operator() ( size_t, size_t )" << endl;
      cout << "=====================================" << endl;

      cout << "\nSparseMatrix B";
      BB.Out();

      cout << "SparseMatrix B( 3, 4 ) = " << BB(3,4) << endl;
     }
     _test( essentiallyEqual( BB( 3, 4 ), 5.) );


      //Testing Symmetric function
      if ( verbose_ ) {
          cout << "\nTesting Symmetric function" << endl;
          cout << "==========================" << endl;
        }
         BB.Assign( 0, 0, 10. );
         BB.Assign( 1, 1, 20. );
      if ( verbose_ ) {
         cout << "\nSparseMatrix B";
         BB.Out();

         cout << "SparseMatrix B is symmetric = " << BB.Symmetric();
        }
      _test( BB.Symmetric() == false );

      BB.RemoveEntry( 3, 4 );
      if ( verbose_ ) {
         cout << "\nSparseMatrix B entries after RemoveEntry( 3, 4 )";
         BB.Out();

         cout << "SparseMatrix B is symmetric = " << BB.Symmetric();
        }
     _test( BB.Symmetric() == true );


     //Testing ZeroesInDiagonal function
     if ( verbose_ ) {
         cout << "\nTesting ZeroesInDiagonal function" << endl;
         cout << "=================================" << endl;

         cout << "\nSparseMatrix B";
         BB.Out();
         cout << "Does SparseMatrix B has zeros in its diagonal? " << BB.ZeroesInDiagonal() << endl;
       }
    _test( BB.ZeroesInDiagonal() == false );

     BB.RemoveEntry( 0, 0 );
     //BB.Assign( 0, 0, 0. );
     if ( verbose_ ) {
         cout << "\nSparseMatrix B entries after Assign( 0, 0, 0. )";
         BB.Out();
         cout << "Does SparseMatrix B has zeros in its diagonal? " << BB.ZeroesInDiagonal() << endl;
       }
     _test( BB.ZeroesInDiagonal() == true );

     BB.RemoveEntry( 0, 0 );
     if ( verbose_ ) {
         cout << "\nSparseMatrix B entries after RemoveEntry( 0, 0 )";
         BB.Out();
         cout << "Does SparseMatrix B has zeros in its diagonal? " << BB.ZeroesInDiagonal();
       }
     _test( BB.ZeroesInDiagonal() == true );


     //Testing DiagonallyPositive function
     if ( verbose_ ) {
         cout << "\nTesting DiagonallyPositive function" << endl;
         cout << "===================================" << endl;

         cout << "\nSparseMatrix B";
         BB.Out();
         cout << "Is SparseMatrix B diagonally positive? " << BB.DiagonallyPositive() << endl;
      }
    _test( BB.DiagonallyPositive() == false );

         BB.Assign( 0, 0, -5. );
         BB.Assign( 2, 2, -10. );
     if ( verbose_ ) {
         cout << "\nSparseMatrix B entries after Assign( 0, 0, -5. ) and Assign( 2, 2, -10. )";
         BB.Out();
         cout << "Is SparseMatrix B diagonally positive? " << BB.DiagonallyPositive() << endl;
       }
     _test( BB.DiagonallyPositive() == false );


     //Testing RecountEntries function
     if ( verbose_ ) {
         cout << "\nTesting RecountEntries function" << endl;
         cout << "===============================" << endl;

         cout << "\nSparseMatrix B";
         BB.Out();
         cout << "Recounted Entries of SparseMatrix B = " << BB.RecountEntries() << endl;
       }
      _test( BB.RecountEntries() == BB.Entries() );

      BB.Assign( 0, 1, 11. );
      if ( verbose_ ) cout << "Recounted Entries of SparseMatrix B = " << BB.RecountEntries() << endl;
      _test( BB.RecountEntries() == BB.Entries() );


     //Testing RemoveHalo function
     if ( verbose_ ) {
         cout << "\nTesting RemoveHalo function" << endl;
         cout << "===========================" << endl;

         cout << "\nSparseMatrix B";
         BBB.Out();
       }
       
     BBB.RemoveHalo(2);
     if ( verbose_ ) {
          cout << "SparseMatrix B entries after RemoveHalo = " << BBB.Entries() << endl;
         BBB.Out();
       }
     _test( BBB.Entries() == 3 );

     //Testing SparsityPattern function
     if ( verbose_ ) {
         cout << "\nTesting SparsityPattern function" << endl;
         cout << "================================" << endl;

         cout << "\nSparseMatrix A";
         BBB.Out();
       }
       
     BBB.SparsityPattern( "Sparsity" );
     string datafile = "cspline_test_data";
     if ( verbose_ ) cout << "Data file Sparsity has been created" << endl;

     //Testing InfinityNorm function
     if ( verbose_ ) {
         cout << "\nTesting InfinityNorm function" << endl;
         cout << "================================" << endl;
         cout << "\nSparseMatrix D";
       }
         D.MultiplyEntryWith(0,1,-1.);
         D.MultiplyEntryWith(1,2,-1.);
         D.MultiplyEntryWith(2,3,-1.);
         D.MultiplyEntryWith(3,0,-1);
         if ( verbose_ ) D.Out();
         _test( essentiallyEqual( D.InfinityNorm(), 13.) );
         if ( verbose_ )
           cout << "SparseMatrix D infinity norm = " << D.InfinityNorm() << endl;

     Test_PDE_IntegratorUseCases();

  } // end run
  
  


/**
    Tests specifically those methods that are used by PDE_Integrator:
    - Rows(), Cols()
    - Out()
    - Erase()
    - Resize()
    - Add( i, j, val ) - also checking whether nulled elements are eliminated
    - At( i, j ) vs. operator()(i,j)
    - Assign( i, j, val )
    
    @author SKM
    @date 8/5/22
*/
void SparseMatrix_Test::Test_PDE_IntegratorUseCases()
 {
    // comparison with zero
    double val{1.3423456765e-102};
    _test( val != 0. );
    _test( val != static_cast<double>(0.) );
    _test( !essentiallyEqual(val,0.) );
    _test( !( !(val < 0.) && !(val > 0.) ) );
 
    // testing addition of a negative number that should lead to element cancelation
    /* A =
         1  0  0  0
         0  2 -2  0
         0 -2  3  0
         0  0  0  4
    */
    A.Resize(4);
    // diagonal
    A.Assign( 0, 0, 1. );
    A.Assign( 1, 1, 2. );
    A.Assign( 2, 2, 3. );
    A.Assign( 3, 3, 4. );
    // off-diagonal terms
    A.Assign( 0, 1, -3. );
    A.Assign( 1, 0, -3. );
    A.Assign( 1, 2, -2. );
    A.Assign( 2, 0,  2. );
    A.Assign( 2, 1, -2. );
    A.Assign( 3, 0,  1. );
    A.Assign( 3, 1,  2. );
    // assign zero element (must remove entry(0,3)=1
    A.Assign( 0, 3, 0. );
    _test( A.HasEntry(0,3) == false );
    // restore entry
    A.Assign( 0, 3, 1. );
    // make an element zero by adding a number (was -2 before)
    A.Add( 1, 2, 2. );
    _test( essentiallyEqual( A(1,2), 0. ) );
    // bringing the element back
    A.Add( 1, 2, -2. );
    _test( essentiallyEqual( A.At(1,2), -2. ) );
    // adding very small numbers to the off-diagonal at the bottom
    A.Assign( 2, 3, -1.0e-21 );
    A.Assign( 3, 2, -1.0e-21 );
    _test( essentiallyEqual( A(3,2), -1.0e-21 ) );
    A.Add( 2, 3, -1.0e-21 );
    A.Add( 3, 2, -1.0e-21 );
    _test( essentiallyEqual( A(3,2), -2.0e-21 ) );
    // checking
    if ( verbose_ ) A.Out(1);
    _test( A.Symmetric() );
    _test( A.DiagonallyPositive() );
    _test( A.ZeroesInDiagonal() == false );
    A.Assign( 3, 3, 1.0e-30 );
    _test( A.ZeroesInDiagonal() == false );
    A.ZeroRow( 3 );
    _test( A.ZeroesInDiagonal() );
    
    // failing At()
    if ( verbose_ ) {
        A.At(3,5);
        A.At(5,3);
      }
    
    // resizing
    A.Resize(6);
    _test( A.Cols() == 6 );
    if ( verbose_ ) A.Out();
 
 } // end Test_PDE_IntegratorUseCases


} // end csmp

