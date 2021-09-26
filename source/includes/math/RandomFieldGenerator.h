#ifndef RANDOM_FIELD_GENERATOR_H
#define RANDOM_FIELD_GENERATOR_H

#include "CSMP_mathUtilities.h"
#include "Model.h"
#include "Matrix.h"


namespace csmp {

template<size_t dim>
class RandomFieldGenerator  {
  public:
    RandomFieldGenerator();
    // generating the random field
    void RandomElementField2D( Model<dim>& sg, const char* variable, double64 mean, double64 sigma, double64 xlength, double64 ylength, 
                               bool logarithmic, const char* region="Model", size_t iterations=50 );
    void RandomNodeField2D( Model<dim>& sg, const char* variable, double64 mean, double64 sigma, double64 xlength, double64 ylength, 
		                    bool logarithmic, const char* region="Model", size_t iterations=50 );
    // saving the random field
    void OutputRandomNodeField( Model<dim>& sg);
    void OutputRandomElementField( Model<dim>& sg );
    // reading in a saved random field
    void InputRandomNodeField( Model<dim>& sg, const char* variable );
    void InputRandomElementField( Model<dim>& sg, const char* variable );

  private:
    csmp::Matrix  m1_, m2_, m3_;
    const double64 pi_;
    std::vector<double64> k_;
    std::string fname;
    Matrix& UniformRandomMatrix(size_t m, size_t n);
    Matrix& NormalRandomMatrix(size_t m, size_t n);
    double64 TheoreticalStandardDeviation(size_t m, double64 xl, double64 yl, double64 Lx, double64 Ly);
};

} // csp

#endif
