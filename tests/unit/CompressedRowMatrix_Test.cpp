#include "CompressedRowMatrix_Test.h"
#include "CompressedRowMatrix.h"

using namespace std;
namespace csmp
{
  void CompressedRowMatrix_Test::run()
  {
    cout << "\n===========================";
    cout << "\nTesting CompressedRowMatrix" << endl;
    cout << "===========================" << endl;


    SparseMatrix A(4);
    A.Assign( 0, 0, 1. );
    A.Assign( 0, 2, 2. );
    A.Assign( 0, 3, 1. );

    A.Assign( 1, 0, 3. );
    A.Assign( 1, 1, 1. );
    A.Assign( 1, 3, 2. );

    A.Assign( 2, 0, 1. );
    A.Assign( 2, 2, 2. );
    A.Assign( 2, 3, 4. );

    A.Assign( 3, 1, 1. );
    A.Assign( 3, 2, 3. );
    A.Assign( 3, 3, 5. );

    std::cout<<"Matrix A set up"<<endl;

	CompressedRowMatrix testMatrix1(A);
	CompressedRowMatrix testMatrix2;
	testMatrix2.InitializePointBased(A,4);
	CompressedRowMatrix testMatrix3(testMatrix1);
	    
	cout<< "\nPrinting the A matrix:\n";
	A.Out();

	cout << "Printing the testMatrix1 initialized by A using the '()' operator:\n";
	cout << testMatrix1(0, 0) << " " << testMatrix1(0, 1) << " " << testMatrix1(0, 2) << " " << testMatrix1(0, 3) << endl;;
	cout << testMatrix1(1, 0) << " " << testMatrix1(1, 1) << " " << testMatrix1(1, 2) << " " << testMatrix1(1, 3) << endl;;
	cout << testMatrix1(2, 0) << " " << testMatrix1(2, 1) << " " << testMatrix1(2, 2) << " " << testMatrix1(2, 3) << endl;;
	cout << testMatrix1(3, 0) << " " << testMatrix1(3, 1) << " " << testMatrix1(3, 2) << " " << testMatrix1(3, 3) << endl;;

	cout << "Printing the testMatrix2 initialized for a PointBased approach\n";
	cout << testMatrix1(0, 0) << " " << testMatrix1(0, 1) << " " << testMatrix1(0, 2) << " " << testMatrix1(0, 3) << endl;;
	cout << testMatrix1(1, 0) << " " << testMatrix1(1, 1) << " " << testMatrix1(1, 2) << " " << testMatrix1(1, 3) << endl;;
	cout << testMatrix1(2, 0) << " " << testMatrix1(2, 1) << " " << testMatrix1(2, 2) << " " << testMatrix1(2, 3) << endl;;
	cout << testMatrix1(3, 0) << " " << testMatrix1(3, 1) << " " << testMatrix1(3, 2) << " " << testMatrix1(3, 3) << endl;;

	cout << "Printing the testMatrix3 initialized by testMatrix1 using the '=' operator:\n";
	cout << testMatrix1(0, 0) << " " << testMatrix1(0, 1) << " " << testMatrix1(0, 2) << " " << testMatrix1(0, 3) << endl;;
	cout << testMatrix1(1, 0) << " " << testMatrix1(1, 1) << " " << testMatrix1(1, 2) << " " << testMatrix1(1, 3) << endl;;
	cout << testMatrix1(2, 0) << " " << testMatrix1(2, 1) << " " << testMatrix1(2, 2) << " " << testMatrix1(2, 3) << endl;;
	cout << testMatrix1(3, 0) << " " << testMatrix1(3, 1) << " " << testMatrix1(3, 2) << " " << testMatrix1(3, 3) << endl;;



	testMatrix1.Out();

    //Testing the () operator
    _equal( testMatrix1( 0, 0 ), 1., tolerance_ );
    _equal( testMatrix1( 3, 3 ), 5., tolerance_ );
    _equal( testMatrix1( 3, 2 ), 3., tolerance_ );

    _equal( testMatrix2( 0, 0 ), 1., tolerance_ );
    _equal( testMatrix2( 3, 3 ), 5., tolerance_ );
    _equal( testMatrix2( 3, 2 ), 3., tolerance_ );
    
    // trying to recuperate an element that is not stored (should be zero)
    _equal( testMatrix2( 2, 1 ), 0., tolerance_ );
    
	//Testing the = operator
	_equal(testMatrix3(0, 0), 1., tolerance_);
	_equal(testMatrix3(3, 3), 5., tolerance_);
	_equal(testMatrix3(3, 2), 3., tolerance_);

    testMatrix1.Out();

    std::cout <<"End of testMatrix1.out()"<<endl;

    //Testing methods to extract total amount of "existing" entries (note: they could be zero any way!)
    _equal( testMatrix1.TotalExistingEntries(), A.Entries(), tolerance_ );
	_equal(testMatrix2.TotalExistingEntries(), A.Entries(), tolerance_);
	_equal(testMatrix3.TotalExistingEntries(), A.Entries(), tolerance_);

    //Testing the size of ia vector
    _equal( testMatrix1.ia.size()-1U, A.Rows(), tolerance_ );
	_equal(testMatrix2.ia.size() - 1U, A.Rows(), tolerance_);
	_equal(testMatrix3.ia.size() - 1U, A.Rows(), tolerance_);

	//Testing the last element of "ia" vector
	_equal(testMatrix1.ia[testMatrix1.Rows() ], A.Entries()+1U, tolerance_);
	_equal(testMatrix2.ia[testMatrix2.Rows() ], A.Entries() + 1U, tolerance_);
	_equal(testMatrix3.ia[testMatrix3.Rows() ], A.Entries() + 1U, tolerance_);

	//Testing the size of "ja" vector
	_equal(testMatrix1.ja.size() , A.Entries(), tolerance_);
	_equal(testMatrix2.ja.size(), A.Entries(), tolerance_);
	_equal(testMatrix3.ja.size(), A.Entries(), tolerance_);

	//Testing the size of "a' vector
	_equal(testMatrix1.a.size(), A.Entries(), tolerance_);
	_equal(testMatrix2.a.size(), A.Entries(), tolerance_);
	_equal(testMatrix3.a.size(), A.Entries(), tolerance_);
  }

}// end namespace csmp
