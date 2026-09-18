#include "CompressedRowMatrix_Test.h"
#include "CompressedRowMatrix.h"
#include "VSet.h"
#include "Model.h"
#include "Region.h"
#include "Element.h"
#include "Node.h"
#include "vsetMakers.h"
#include "Exception.h"

using namespace std;

namespace csmp {

void CompressedRowMatrix_Test::run()
  {
    if ( verbose_ ) {
        cout << "\n===========================";
        cout << "\nTesting CompressedRowMatrix" << endl;
        cout << "===========================" << endl;
      }

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

    if ( verbose_ ) std::cout<<"Matrix A set up"<<endl;

  const bool SAMG_formatting{ true };
	CompressedRowMatrix testMatrix1( A, SAMG_formatting );
  _test( testMatrix1.IsFormattedForSAMG() == true );
  _test( testMatrix1.VerifySparsityPatternSAMG() == true );
 
 
	CompressedRowMatrix testMatrix2;
	testMatrix2.InitializePointBasedSAMG(A,4);
  _test( testMatrix1.IsFormattedForSAMG() == true );
  _test( testMatrix1.VerifySparsityPatternSAMG() == true );

	CompressedRowMatrix testMatrix3(testMatrix1);
 
  // keep in zero-based format
  testMatrix1.Initialize( A );
  _test( testMatrix1.IsFormattedForSAMG() == false );
  _test( testMatrix1.VerifySparsityPattern() == true );
  // transfer to Fortran
  testMatrix1.ConvertToSAMGFormat();
	
  if ( verbose_ ) {
      cout<< "\nPrinting the A matrix:\n";
      A.Out(); // prints only the non-zero elements

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
    }

    //Testing the () operator
    _equal( testMatrix1( 0, 0 ), 1., tolerance_ );
    _equal( testMatrix1( 3, 3 ), 5., tolerance_ );
    _equal( testMatrix1( 3, 2 ), 3., tolerance_ );

    _equal( testMatrix2( 0, 0 ), 1., tolerance_ );
    _equal( testMatrix2( 3, 3 ), 5., tolerance_ );
    _equal( testMatrix2( 3, 2 ), 3., tolerance_ );
    
    // trying to recuperate an element that is not stored (should be zero)
    _equal( testMatrix2( 2, 1 ), 0., tolerance_ );
    
    //Testing operator= and member function At()
    _equal(testMatrix3(0, 0), 1., tolerance_);
    _equal(testMatrix3(3, 3), 5., tolerance_);
    _equal(testMatrix3(3, 2), 3., tolerance_);

    if ( verbose_ ) std::cout <<"End of testMatrix1.out()"<<endl;

    //Testing methods to extract total amount of "existing" entries (note: they could be zero any way!)
    _test(testMatrix1.Entries() == A.Entries() );
	  _test(testMatrix2.Entries() == A.Entries() );
	  _test(testMatrix3.Entries() == A.Entries() );

    //Testing the size of ia vector
    _test( testMatrix1.ia.size() - 1U == A.Rows() );
	  _test( testMatrix2.ia.size() - 1U == A.Rows() );
	  _test( testMatrix3.ia.size() - 1U == A.Rows() );

    //Testing the last element of "ia" vector
    _test(testMatrix1.ia[ testMatrix1.Rows() ] == static_cast<int32_t>(A.Entries() + 1) );
    _test(testMatrix2.ia[ testMatrix2.Rows() ] == static_cast<int32_t>(A.Entries() + 1) );
    _test(testMatrix3.ia[ testMatrix3.Rows() ] == static_cast<int32_t>(A.Entries() + 1) );

    //Testing the size of "ja" vector
    _test(testMatrix1.ja.size() == A.Entries() );
    _test(testMatrix2.ja.size() == A.Entries() );
    _test(testMatrix3.ja.size() == A.Entries() );

    //Testing the size of "a' vector
    _test(testMatrix1.a.size() == A.Entries() );
    _test(testMatrix2.a.size() == A.Entries() );
    _test(testMatrix3.a.size() == A.Entries() );
  
  // testing other members of compressed row matrix
  testMatrix1.Zero();
  // is the storage preserved
  _test( testMatrix1.Entries() == A.Entries() );
  
  TestCRM_Multiplication();
  
  Test_generateSparsityPattern();
  
} // end run
  
  
  
  /**
 * Test for CompressedRowMatrix::MultiplyWith
 * Uses Assign() to build the test matrix and _equal() for verification.
 */
void CompressedRowMatrix_Test::TestCRM_Multiplication()
{
    tolerance_ = 1e-12;
    const size_t nnu = 3;
    
    // 1. Initialize the structure
    // (setting up the sparsity pattern/size first)
    vector<int32_t> ia{0, 2, 5, 7}; /*ia*/
    vector<int32_t> ja{0, 1,  1, 0, 2,  2, 1}; /*ja*/
    vector<double>  a{0.,0.,0.,0.,0.,0.,0.}; /*a*/
    CompressedRowMatrix A( ia, ja, a );
    // ... logic to size A to 3x3 if needed ...

    // 2. Build the Matrix A using Assign()
    // Matrix:
    // [ 10.0   2.0   0.0 ]
    // [  3.0  20.0   4.0 ]
    // [  0.0   5.0  30.0 ]
    
    // CRITICAL: Assign diagonal FIRST for each row to satisfy SAMG/CRM convention
    A.Assign(0, 0, 10.0); // Row 0 Diagonal
    A.Assign(0, 1,  2.0); 

    A.Assign(1, 1, 20.0); // Row 1 Diagonal
    A.Assign(1, 0,  3.0); 
    A.Assign(1, 2,  4.0);

    A.Assign(2, 2, 30.0); // Row 2 Diagonal
    A.Assign(2, 1,  5.0);

    // 3. Setup vectors for y = Ax
    std::vector<double> x = {1.0, 2.0, 3.0};
    std::vector<double> y(nnu, 0.0);

    // 4. Perform Multiplication
    A.MultiplyWith(x, y);

    // 5. Automated Verification
    // y[0] = 10(1) + 2(2) = 14
    // y[1] = 3(1) + 20(2) + 4(3) = 55
    // y[2] = 5(2) + 30(3) = 100

    _equal( y[0],  14.0, tolerance_ );
    _equal( y[1],  55.0, tolerance_ );
    _equal( y[2], 100.0, tolerance_ );

    // 6. Final check: Verify Diagonal-First storage via Out() or direct access
    // This confirms Assign() placed the diagonal at the correct physical offset
    _equal( A.At(0, 0), 10.0, tolerance_ );
    _equal( A.At(1, 1), 20.0, tolerance_ );
    _equal( A.At(2, 2), 30.0, tolerance_ );
    
} // end TestCRM_Multiplication


// helper function: build degree of freedom mapping
static void  buildDOF_Mapping( const map<Parameter,size_t> test_operands,
                               const ModelSubDomain<3,Element>& gref,
                               vector<size_t>& DOF_indexes )
  {
    const uint32_t dim{3};
    gref.RenumberNodes();
    size_t DOF{0U};
    
    for ( auto& it : test_operands )
      {
        auto        niter(gref.NodesBegin());
        csmp::Index prop_key = it.first.key;
        size_t      offset = it.second;

        if (prop_key.place != NODE)
          throw csmp::Exception(ERROR, 
            "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::ReduceSystemSizeEliminatingEssentialConditions",
            "So far no conditions are assigned to elements, faces, segments");

        size_t  position{0U};
        switch (prop_key.type) {
            case SCALAR:
              while (niter != gref.NodesEnd()) {
                  position = (*niter)->Idx() + offset;
                  assert( position < DOF_indexes.size() );
                  if ((*niter)->Status(prop_key) == DIRICH) 
                    DOF_indexes[position] = NULL_IDX;
                  else {
                      DOF_indexes[position] = DOF;
                      DOF = DOF + 1U;
                    }
                  niter++;
                }
              break;
            case VECTOR:
              while ( niter != gref.NodesEnd()) {
                     for ( uint32_t i{0U}; i < dim; ++i ) {
                          position = (*niter)->Idx() * dim + i + offset;
                          assert( position < DOF_indexes.size() );
                          if ( (*niter)->Status(prop_key,i) == DIRICH ) 
                            DOF_indexes[position] = NULL_IDX;
                          else {
                               DOF_indexes[position] = DOF;
                               DOF = DOF + 1U;
                            }
                       }
                    niter++;
                 }
              break;
            default:
              throw csmp::Exception(FATAL_ERROR,
                "PDE_Integrator<dim,CELLTYPE,MATRIXTYPE>::ReduceSystemSizeEliminatingEssentialConditions",
                "Variable type not recognised by this method");
            }
          
      } // end for (all Dirichlet flagged variables)

  } // end build DOF mapping




void CompressedRowMatrix_Test::Test_generateSparsityPattern()
  {
     VSet<3U> vset;
     create_Prism_Hexa_VSet( vset, true /* Skewed */ ); // includes Box boundary flagging
     Model<3> model( vset, "CSMP-1phase-variables.txt" );
     // boundary conditions
     model.InputBoundaryValue( LEFT, "fluid pressure", makeScalar(DIRICH,2e5) );
     model.InputBoundaryValue( RIGHT, "fluid pressure", makeScalar(DIRICH,1e5) );
     
     // 1. single SCALAR variable test case
     // -----------------------------------
     // build PDE_Integrator internals
     map<Parameter,size_t> test_operands;
     size_t                total_degrees_of_freedom{ vset.Vertices() }; // scalar "fluid pressure" placed on Node
     vector<size_t>        DOF_indexes(total_degrees_of_freedom,0);
     
     test_operands.insert( make_pair(model.Database().Parameter("fluid pressure"),0) );
     
     buildDOF_Mapping( test_operands, model.Region("Model"), DOF_indexes );
     
     // deprecated: generateSparsityPatternEliminatingEssentialConditions( CRM, test_operands, DOF_indexes, model.Region("Model") );
     SparsityPattern sp = generateSparsityPattern( test_operands, DOF_indexes, model.Region("Model") );

     CompressedRowMatrix  CRM( sp.ia, sp.ja, sp.a );

     _test( CRM.VerifySparsityPattern() == true );
     
  } // end Test_generateSparsityPattern



} // end namespace csmp
