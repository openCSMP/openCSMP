//
//  AccumulationSpeedProfiling_Test.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 30/11/2024.
//

#include "AccumulationSpeedProfiling_Test.h"
#include "Model.h"
#include "ANSYS_Model3D.h"
#include "NumIntegral_dNT_op_dN_dV.h"
#include "NumIntegral_NT_op_N_dV.h"
#include "VelocityAndVolumeFlux.h"

#include <Eigen/Dense>

using namespace std;

namespace csmp {

AccumulationSpeedProfiling_Test::~AccumulationSpeedProfiling_Test()
 {
    delete model_ptr_;
 }



void AccumulationSpeedProfiling_Test::run()
 {
    BuildModel3D();
    AccumulateIntegralOnPolyhedralMesh();
    MatricesAsFunctionVersusGlobalVariables();
    MatrixMultiplication_DenseMatrix_vs_Eigen_matrix();
    EigenMatrixTemplateTest();
    EigenMatrixDynamicArgumentsTemplateTest();
    EigenDynamicMatrixTest();
    MatrixResizeTest();
    
 } // end run
 
 
Model<3U>* AccumulationSpeedProfiling_Test::BuildModel3D()
 {
    ANSYS_Model3D m0( "prism_test", "prism_test", "CSMP-1phase-variables.txt", true );
    m0.OutputToBinaryFile("PrismTest_speed_test");
    
    model_ptr_ = new Model<3U>( string("PrismTest_speed_test") );
   
    // setting up necessary variables
    model_ptr_->InputPropertyValue( "permeability", makeScalar( PLAIN, 1.0e-12 ) );
    model_ptr_->InputPropertyValue( "conductivity", makeScalar( PLAIN, 1.0e-9 ) );
    model_ptr_->InputPropertyValue( "porosity", makeScalar( PLAIN, 0.25 ) );
    model_ptr_->InputPropertyValue( "fluid volume source", makeScalar( PLAIN, 0. ) );
    model_ptr_->InputPropertyValue( "velocity", makeVector( ANY, ANY, ANY, 0., 0., 0. ) );
    model_ptr_->InputPropertyValue( "pore velocity", makeVector( ANY, ANY, ANY, 0., 0., 0. )  );
    model_ptr_->InputPropertyValue( "volume flux", makeScalar( PLAIN, 0. ) );
    
    return model_ptr_;

 } // end BuildModel3D
 
 
 
/*
   Learnings: reading CSMP's distributed storage is sooo slow. It even eclipses the matrix multiplications
*/
void AccumulationSpeedProfiling_Test::AccumulateIntegralOnPolyhedralMesh()
 {
    if ( verbose_ )
      cout <<"\n"<<"AccumulationSpeedProfiling_Test::AccumulateIntegralOnPolyhedralMesh: testing..."<< endl;

    // create some common PDE operators for testing
    // --------------------------------------------
    NumIntegral_dNT_op_dN_dV<3U> lap( model_ptr_->Database(), "conductivity", "fluid pressure",  "fluid pressure" );

    // 1. accumulating into the sparse matrix 'G'
    // ------------------------------------------
    Region<3U>&  model_domain = model_ptr_->Region("Model");
    SparseMatrix G( model_domain.Cells() );

    vector<size_t>   DOF_indexes( model_domain.Nodes() ); ///< indices of DOFs, but only of the non-Dirichlet dofs, size enumerated 0 - DOF-1 (including Dirich DOF)
    vector<double>   pivotVector( model_domain.Nodes() ); ///< mapping from DOFs to actual node numbers

    double           dt{1.7};

    iota( DOF_indexes.begin(), DOF_indexes.end(), 0 ); //all indices are used
    fill( pivotVector.begin(), pivotVector.end(), 0. );

    // accumulating the Laplacian
    auto t0 = chrono::high_resolution_clock::now();
    for ( const auto& it : model_domain.CellVector() )
     {
       lap.GetOperands( *it );
       lap.ComputeContribution( *it );
       lap.MultiplyWithTimeFactor( dt );
       lap.AssignToGlobal( *it, G, pivotVector, DOF_indexes );
     }
    auto t1 = chrono::high_resolution_clock::now();
	  cout <<"\n\t"<<"1. Time to accumulate Laplacian: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
       
              
    // 2. accumulating into the righhand vector 'rhs'
    // --------------------------------------------------------------
    NumIntegral_NT_op_N_dV<3U>   src( model_ptr_->Database(),  "fluid volume source", "fluid pressure" );
    vector<double> rhs( model_domain.Nodes(), 0. );

    // accumulating of mass matrix
    t0 = chrono::high_resolution_clock::now();
    for ( const auto& it : model_domain.CellVector() )
       {
         src.GetOperands( *it );
         src.ComputeContribution( *it );
         src.MultiplyWithTimeFactor( dt );
         src.AssignToGlobal( *it, rhs, DOF_indexes );
       }
    t1 = chrono::high_resolution_clock::now();
	  cout <<"\n\t"<<"2. Time to accumulate mass matrix into 'rhs' vector: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;


    // 3. postprocessing velocities etc.
    // --------------------------------------------------------------
    VelocityAndVolumeFlux<3U>  vel( *model_ptr_,  "conductivity", "porosity", "fluid pressure", false );
    // needs to get to postprocessing
    vel.ApplicationCycle( 1 );
    // generating a fluid pressure gradient field from node positions
    model_ptr_->AssignNodeCoordinatesTo( "fluid pressure", 'y' ); // x, y, z
//    printRangeOfVariable( *model_ptr_, "fluid pressure" );

    // accumulating of mass matrix
    t0 = chrono::high_resolution_clock::now();
    for ( const auto& it : model_domain.CellVector() )
       {
         vel.GetOperands( *it );
         vel.ComputeContribution( *it );
         vel.WriteOperands( *it );
       }
    t1 = chrono::high_resolution_clock::now();
	  cout <<"\n\t"<<"3. Time to postprocess fluid velocity: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
   
    printRangeOfVariable( *model_ptr_, "velocity" );
//    printRangeOfVariable( *model_ptr_, "pore velocity" );
//    printRangeOfVariable( *model_ptr_, "volume flux" );

 } // end AccumulateIntegralOnPolyhedralMesh




// same speed as when C is preallocated
static double matrixMultiplication( const DenseMatrix<DM6>& A,  const DenseMatrix<DM6>& B )
 {
     DenseMatrix<DM6> C = A * B;
     return ( C(0,0) + C(1,1) + C(2,2) + C(3,3) + C(4,4) );
 }

// speed reduced to 2/3 relative to matrixMultiplication
static double matrixMultiplicationInsideFunction()
 {
     DenseMatrix<DM6> A = { {1,2,3,4,5,6}, {6,5,4,3,2,1}, {1,2,3,4,5,6}, {6,5,4,3,2,1}, {6,7,8,9,10,11}, {6,7,8,9,10,11} };
     DenseMatrix<DM6> B = { {1,1,1,1,1,1}, {2,2,2,2,2,2}, {3,3,3,3,3,3}, {4,4,4,4,4,4}, {5,5,5,5,5,5}, {6,6,6,6,6,6} };
     DenseMatrix<DM6> C = A * B;
     return ( C(0,0) + C(1,1) + C(2,2) + C(3,3) + C(4,4) + C(5,5) );
 }

/*
    Learnings: storing a matrix locally has very little impact on speed
*/
void AccumulationSpeedProfiling_Test::MatricesAsFunctionVersusGlobalVariables()
 {
    cout <<"\n\n"<<"AccumulationSpeedProfiling_Test::MatricesAsFunctionVersusGlobalVariables: testing..."<< endl;
  
     DenseMatrix<DM6> A = { {1,2,3,4,5,6}, {6,5,4,3,2,1}, {1,2,3,4,5,6}, {6,5,4,3,2,1}, {6,7,8,9,10,11}, {6,7,8,9,10,11} };
     DenseMatrix<DM6> B = { {1,1,1,1,1,1}, {2,2,2,2,2,2}, {3,3,3,3,3,3}, {4,4,4,4,4,4}, {5,5,5,5,5,5}, {6,6,6,6,6,6} };
     DenseMatrix<DM6> C;
     
     //A.Out();
     //B.Out();
     
     // 0. base case: global matrix multiplication using operator overloading
     double tracker{0.};
     auto t0 = chrono::high_resolution_clock::now();
     for ( int i{0}; i<1e6; ++i ) {
          C = A * B;
          tracker += ( C(0,0) + C(1,1) + C(2,2) + C(3,3) + C(4,4) + C(5,5) );
       }
     auto t1 = chrono::high_resolution_clock::now();
	   cout <<"\n\t"<<"1. Time for 10^6 multiplications of 5x5 matrices: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
     cout << tracker << endl;
     //C.Out();

     // 1. using a function where the C matrix is a local variable
     tracker = 0.;
     t0 = chrono::high_resolution_clock::now();
     for ( int i{0}; i<1e6; ++i ) {
          tracker += matrixMultiplication( A, B );
       }
     t1 = chrono::high_resolution_clock::now();
	   cout <<"\n\t"<<"1. Time for 10^6 multiplications of 5x5 matrices; matrix C within function: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
     cout << tracker << endl;
     //C.Out();

     // 2. using a function in which the matrices are generated and multiplied
     tracker = 0.;
     t0 = chrono::high_resolution_clock::now();
     for ( int i{0}; i<1e6; ++i ) {
          tracker += matrixMultiplicationInsideFunction();
       }
     t1 = chrono::high_resolution_clock::now();
	   cout <<"\n\t"<<"1. Time for 10^6 multiplications of 5x5 matrices; A, B, C inside function: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
     cout << tracker << endl;
     //C.Out();

 } // end MatricesAsFunctionVersusGlobalVariables



/*
    Learnings: for C = A* B, csmp::DenseMatrix is faster than Eigen's dynamc MatrixXd; both are significantly slower than Eigen's full template
    For chains of multiplications, the expression templates kick in and even MatrixXd is significantly faster
    Eigen's matrices have much more functionality and shine when operations are strung together
*/
void AccumulationSpeedProfiling_Test::MatrixMultiplication_DenseMatrix_vs_Eigen_matrix()
 {
     cout <<"\n\n"<<"AccumulationSpeedProfiling_Test::MatrixMultiplication_DenseMatrix_vs_Eigen_matrix: testing..."<< endl;
     // Create a random number generator
     std::random_device rd;  // Seed for the random number generator
     std::mt19937 gen(rd()); // Standard Mersenne Twister generator
     std::uniform_real_distribution<> dis(1.0, 100.0); // Range [1, 100]
  
     // CSMP++
     DenseMatrix<DM6> A;
     DenseMatrix<DM6> B;
     DenseMatrix<DM6> C;
     
     // EIGEN (dynamic matrix) using expression templated multiplication
     Eigen::MatrixXd EA(6,6);
     Eigen::MatrixXd EB(6,6);
     Eigen::MatrixXd EC(6,6);
     // templatised version
     Eigen::Matrix<double,6,6> ETA;
     Eigen::Matrix<double,6,6> ETB;
     Eigen::Matrix<double,6,6> ETC;
     
     // 0. matrix multiplication using csmp::DenseMatrix
     double tracker{0.};
     auto t0 = chrono::high_resolution_clock::now();
     for ( int n{0}; n<1e6; ++n ) {
          // matrix initialisation
          for (int i = 0; i < A.Rows(); ++i)
            for (int j = 0; j < A.Cols(); ++j) {
                 A(i, j) = dis(gen);
                 B(i, j) = dis(gen);
              }
          // calculation
          C  = A * B;
          C += A;
          C += B;
          tracker += ( C(0,0) + C(1,1) + C(2,2) + C(3,3) + C(4,4) + C(5,5) );
       }
     auto t1 = chrono::high_resolution_clock::now();
	   cout <<"\n\t"<<"1. Time for 10^6 multiplications of 5x5 csmp::DenseMatrix objects: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
     cout << tracker << endl;
     //C.Out();
     

     // 1. matrix multiplication using Eigen::MatrixXd
     tracker = 0.;
     t0 = chrono::high_resolution_clock::now();
     for ( int n{0}; n<1e6; ++n ) {
          // matrix initialisation
          for (int i = 0; i < EA.rows(); ++i)
            for (int j = 0; j < EA.cols(); ++j) {
                 EA(i, j) = dis(gen);
                 EB(i, j) = dis(gen);
              }
          // calculation
          EC = EA * EB + EA + EB;
          tracker += ( EC(0,0) + EC(1,1) + EC(2,2) + EC(3,3) + EC(4,4) + EC(5,5) );
       }
     t1 = chrono::high_resolution_clock::now();
	   cout <<"\n\t"<<"1. Time for 10^6 multiplications of 5x5 Eigen::MatrixXd objects: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
     cout << tracker << endl;
     //C.Out();


     // 2. matrix multiplication using Eigen::Matrix<>
     tracker = 0.;
     t0 = chrono::high_resolution_clock::now();
     for ( int n{0}; n<1e6; ++n ) {
          // matrix initialisation
          for (int i = 0; i < ETA.rows(); ++i)
            for (int j = 0; j < ETA.cols(); ++j) {
                 ETA(i, j) = dis(gen);
                 ETB(i, j) = dis(gen);
              }
          // calculation
          ETC = ETA * ETB + ETA + ETB;
          tracker += ( ETC(0,0) + ETC(1,1) + ETC(2,2) + ETC(3,3) + ETC(4,4) + ETC(5,5) );
       }
     t1 = chrono::high_resolution_clock::now();
	   cout <<"\n\t"<<"1. Time for 10^6 multiplications of 5x5 Eigen::Matrix<> objects: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
     cout << tracker << endl;
     //C.Out();
 }



// TEST OF EIGEN Dense matrix class Eigen::Matrix<double, int, int>

// Function with local matrix creation
static double computeResultLocal( const Eigen::Matrix<double, 3, 2>& B,
                                  const Eigen::Matrix<double, 2, 2>& P,
                                  const Eigen::Matrix<double, 2, 3>& A,
                                  double scalar ) {
    Eigen::Matrix<double, 3, 3> result = (B * (P * A)) * scalar;
    return result(1,1);
}

// Function with pass-by-reference
static double computeResultReference( Eigen::Matrix<double, 3, 3>& result,
                                      const Eigen::Matrix<double, 3, 2>& B,
                                      const Eigen::Matrix<double, 2, 2>& P,
                                      const Eigen::Matrix<double, 2, 3>& A,
                                      double scalar ) {
    result = (B * (P * A)) * scalar;
    return result(1,1);
}

/* Testing Eigen::Matrix<double, int, int>
    
   Learnings: no speed difference due to presence of local variable, 3x speed as compared with Eigen::MatrixXd
*/
void AccumulationSpeedProfiling_Test::EigenMatrixTemplateTest()
 {
    Eigen::Matrix<double, 3, 2> B;
    Eigen::Matrix<double, 2, 2> P;
    Eigen::Matrix<double, 2, 3> A;
    B.setRandom();
    P.setRandom();
    A.setRandom();
    double scalar{5.0}, tracker{0.};

    // Timing local creation
    auto startLocal = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1e6; ++i) {
        P.setRandom();
        tracker += computeResultLocal(B, P, A, scalar);
    }
    auto endLocal = std::chrono::high_resolution_clock::now();

    // Timing pass-by-reference
    Eigen::Matrix<double, 3, 3> result;
    auto startRef = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1e6; ++i) {
        P.setRandom();
        tracker += computeResultReference(result, B, P, A, scalar);
    }
    auto endRef = std::chrono::high_resolution_clock::now();

    std::cout << "\nAccumulationSpeedProfiling_Test::EigenMatrixTemplateTest: Local creation time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endLocal - startLocal).count()
              << " ms\n";

    std::cout << "\nAccumulationSpeedProfiling_Test::EigenMatrixTemplateTest: Pass-by-reference time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endRef - startRef).count()
              << " ms\n";
  
    cout << tracker << endl;
              
 } // end EigenMatrixTemplateTest




// TEST OF EIGEN Dense matrix class Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>

// Function with local matrix creation
static double computeResultLocal1( const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& B,
                                  const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& P,
                                  const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& A,
                                  double scalar ) {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> result = (B * (P * A)) * scalar;
    return result(1,1);
}

// Function with pass-by-reference
static double computeResultReference1( Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& result,
                                      const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& B,
                                      const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& P,
                                      const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>& A,
                                      double scalar ) {
    result = (B * (P * A)) * scalar;
    return result(1,1);
}

/* Testing Eigen::Matrix<double, int, int>
    
   Learnings: no speed difference due to presence of local variable, 3x speed as compared with Eigen::MatrixXd
*/
void AccumulationSpeedProfiling_Test::EigenMatrixDynamicArgumentsTemplateTest()
 {
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> B(3,2);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> P(2,2);
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> A(2,3);
    B.setRandom();
    P.setRandom();
    A.setRandom();
    double scalar{5.0}, tracker{0.};

    // Timing local creation
    auto startLocal = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1e6; ++i) {
        P.setRandom();
        tracker += computeResultLocal1(B, P, A, scalar);
    }
    auto endLocal = std::chrono::high_resolution_clock::now();

    // Timing pass-by-reference
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> result;
    auto startRef = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1e6; ++i) {
        P.setRandom();
        tracker += computeResultReference1(result, B, P, A, scalar);
    }
    auto endRef = std::chrono::high_resolution_clock::now();

    std::cout << "\nAccumulationSpeedProfiling_Test::EigenMatrixDynamicArgumentsTemplateTest: Local creation time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endLocal - startLocal).count()
              << " ms\n";

    std::cout << "\nAccumulationSpeedProfiling_Test::EigenMatrixDynamicArgumentsTemplateTest: Pass-by-reference time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endRef - startRef).count()
              << " ms\n";
  
    cout << tracker << endl;
              
 } // end EigenMatrixDynamicArgumentsTemplateTest






// TESTS OF EIGEN Dense matrix class Eigen::MatrixXd( int, int )

// Function with local matrix creation
static double computeDynamicResultLocal( const Eigen::MatrixXd& BT,
                                         const Eigen::MatrixXd& MTRL,
                                         const Eigen::MatrixXd& B,
                                         double scalar )
{
    Eigen::MatrixXd LHS = (BT * (MTRL * B)) * scalar;
    return LHS(1,1);
}

// Function with pass-by-reference
static double computeDynamicResultReference( Eigen::MatrixXd& LHS,
                                             const Eigen::MatrixXd& BT,
                                             const Eigen::MatrixXd& MTRL,
                                             const Eigen::MatrixXd& B,
                                             double scalar ) {
    LHS = (BT * (MTRL * B)) * scalar;
    return LHS(1,1);
}

/* Testing Eigen::Matrix<double, int, int>
    
   Learnings: computeDynamicResultLocal() is 30% slower than computeDynamicResultReference
   - overall speed is 4x less than for Eigen::Matrix<double,int,int>
*/
void AccumulationSpeedProfiling_Test::EigenDynamicMatrixTest()
 {
    Eigen::MatrixXd BT(3,2);
    Eigen::MatrixXd MTRL(2,2);
    Eigen::MatrixXd B(2,3);
    BT.setRandom();
    MTRL.setRandom();
    B.setRandom();
    double scalar{5.0}, tracker{0.};

    // Timing local creation
    auto startLocal = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1e6; ++i) {
        MTRL.setRandom();
        tracker += computeDynamicResultLocal( BT, MTRL, B, scalar );
    }
    auto endLocal = std::chrono::high_resolution_clock::now();

    // Timing pass-by-reference
    Eigen::MatrixXd LHS(3,3);
    auto startRef = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 1e6; ++i) {
        MTRL.setRandom();
        tracker += computeDynamicResultReference( LHS, BT, MTRL, B, scalar );
    }
    auto endRef = std::chrono::high_resolution_clock::now();

    std::cout << "\nAccumulationSpeedProfiling_Test::EigenDynamicMatrixTest: Local creation time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endLocal - startLocal).count()
              << " ms\n";

    std::cout << "\nAccumulationSpeedProfiling_Test::EigenDynamicMatrixTest: Pass-by-reference time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endRef - startRef).count()
              << " ms\n";
  
    cout << tracker << endl;
              
 } // end EigenDynamicMatrixTest


/*
    Compares resizing of csmp::DenseMatrix with Eigen::MatrixXd and EigenMatrix<double,Eigen::Dynamic,Eigen::Dynamic>
*/
void AccumulationSpeedProfiling_Test::MatrixResizeTest()
 {
     cout <<"\n\n"<<"AccumulationSpeedProfiling_Test::MatrixResizeTest: testing..."<< endl;
     // Create a random number generator
     std::random_device rd;  // Seed for the random number generator
     std::mt19937 gen(rd()); // Standard Mersenne Twister generator
     std::uniform_int_distribution<> dis(1,36); // Range [1, 100]
  
     // CSMP++
     {
       DenseMatrix<DM_MIN> CSMPMAT;
       CSMPMAT.Zero();
       // resizing = just changing rows_ and cols_ (high-mem policy)
       double tracker{0.};
       auto t0 = chrono::high_resolution_clock::now();
       for ( int n{0}; n<1e7; ++n ) {
            // matrix resize
            CSMPMAT.Resize( dis(gen), dis(gen) );
            tracker += CSMPMAT(0,0);
         }
       auto t1 = chrono::high_resolution_clock::now();
       cout <<"\n\t"<<"1. Time for 10^7 resize operations of csmp::DenseMatrix: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
       cout << tracker << endl;
     }
     
     // EIGEN (dynamic matrix) using expression templated multiplication
     {
       Eigen::MatrixXd EM(DM_MIN,DM_MIN);
       double tracker{0.};
       auto t0 = chrono::high_resolution_clock::now();
       for ( int n{0}; n<1e7; ++n ) {
            // matrix resize
            EM.resize( dis(gen), dis(gen) );
            tracker += EM(0,0);
         }
       auto t1 = chrono::high_resolution_clock::now();
       cout <<"\n\t"<<"2. Time for 10^7 resize operations of Eigen::MatrixXd: "<< chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
       cout << tracker << endl;
     }
 
     // templatised version
     {
       Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> EDYN(DM_MIN, DM_MIN);
       double tracker{0.};
       auto t0 = chrono::high_resolution_clock::now();
       for ( int n{0}; n<1e7; ++n ) {
            // matrix resize
            EDYN.resize( dis(gen), dis(gen) );
            tracker += EDYN(0,0);
         }
       auto t1 = chrono::high_resolution_clock::now();
       cout <<"\n\t"<<"3. Time for 10^7 resize operations of Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic>: ";
       cout << chrono::duration_cast<chrono::milliseconds>(t1-t0).count() << " milliseconds." << endl;
       cout << tracker << endl;
     }
     
    // checking whether fixed-size Eigen matrices can be assigned to dynamic ones (interoperability)
    // (without resizing the result matrix this fails and may even corrupt memory)
    {
      cout <<"\n"<<"Matrix addition test 1:";
      Eigen::MatrixXd EM(4,4);
      EM.setConstant(5.);
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> EDYN(7,7);
      EDYN.resize(4,4); // ALWAYS BEFORE VALUE ASSIGNMENT
      EDYN.setZero();
      EDYN += EM;
      cout <<"\n\t"<<"Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> after addition of MatrixXd (does not work): "<< EDYN << endl;
    }
    {
      cout <<"\n"<<"Matrix addition test 2:";
      Eigen::Matrix<double,4,4> EM;
      EM.setConstant(5.);
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> EDYN(7,7);
      EDYN.resize(4,4);
      EDYN.setZero();
      EDYN += EM;
      cout <<"\n\t"<<"Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> after addition of Eigen::Matrix<double,4,4>: "<< EDYN << endl;
    }
    {
      cout <<"\n"<<"Matrix addition test 3 (block insertion):";
      Eigen::Matrix<double,4,4> EM;
      EM.setConstant(5.);
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> EDYN(7,7);
      // top left corner (EDYN gets zero'd out before this operation)
      EDYN.block(0,0,4,4) += EM;
      cout <<"\n\t"<<"Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> after addition of Eigen::Matrix<double,4,4>: "<< EDYN << endl;
    }
   
    
 } // end MatrixResizeTest()



} // end csmp
