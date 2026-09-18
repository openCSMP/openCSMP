//
//  AccumulationSpeedProfiling_Test.h
//  Open CSMP++
//
//  Created by Stephan Matthai on 30/11/2024.
//

#ifndef CSMP_ACCUMULATION_SPEED_PROFILING_TEST_H
#define CSMP_ACCUMULATION_SPEED_PROFILING_TEST_H

#include "CSMP_definitions.h"
#include "Test.h"

namespace csmp {

template<uint32_t> class Model;

class AccumulationSpeedProfiling_Test : public Test {
  public:
    AccumulationSpeedProfiling_Test() = default;
    ~AccumulationSpeedProfiling_Test();
    
    void run() override;
    
    const static bool verbose_ = true;
    
  private:
    Model<3U>* BuildModel3D();
    void AccumulateIntegralOnPolyhedralMesh();
    void MatricesAsFunctionVersusGlobalVariables();
    void MatrixMultiplication_DenseMatrix_vs_Eigen_matrix();
    void EigenMatrixTemplateTest();
    void EigenMatrixDynamicArgumentsTemplateTest();
    void EigenDynamicMatrixTest();
    void MatrixResizeTest();
    
    Model<3U>* model_ptr_ = nullptr;
};

} // end csmp

#endif /* CSMP_ACCUMULATION_SPEED_PROFILING_TEST_H */
