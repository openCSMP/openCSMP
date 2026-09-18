// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef RANDOM_FIELD_GENERATOR_H
#define RANDOM_FIELD_GENERATOR_H

#include "CSMP_mathUtilities.h"
#include "Model.h"
#include "Matrix.h"


namespace csmp {

template<uint32_t dim>
class RandomFieldGenerator  {
  public:
    RandomFieldGenerator();
    // generating the random field
    void RandomElementField2D( Model<dim>& sg, const char* variable, double mean, double sigma, double xlength, double ylength, 
                               bool logarithmic, const char* region="Model", size_t iterations=50 );
    void RandomNodeField2D( Model<dim>& sg, const char* variable, double mean, double sigma, double xlength, double ylength, 
		                    bool logarithmic, const char* region="Model", size_t iterations=50 );
    // saving the random field
    void OutputRandomNodeField( Model<dim>& sg);
    void OutputRandomElementField( Model<dim>& sg );
    // reading in a saved random field
    void InputRandomNodeField( Model<dim>& sg, const char* variable );
    void InputRandomElementField( Model<dim>& sg, const char* variable );

  private:
    csmp::Matrix  m1_, m2_, m3_;
    const double pi_;
    std::vector<double> k_;
    std::string fname;
    Matrix& UniformRandomMatrix(size_t m, size_t n);
    Matrix& NormalRandomMatrix(size_t m, size_t n);
    double TheoreticalStandardDeviation(size_t m, double xl, double yl, double Lx, double Ly);
};

} // csp

#endif
