#ifndef CFL_TESTCASE_H
#define CFL_TESTCASE_H


#include <iostream>
#include <iomanip>
#include <limits>
#include "Test.h"


namespace csmp{

template<size_t> class Model;
struct CFL_TestData;

class CFL_TestCase : public Test
{
public:
    CFL_TestCase();
    ~CFL_TestCase();

    virtual void run();

    /// Set model
    void SetModel( Model<3>& model,
                   double mobilityRatio,
                   double velX,
                   double phi );

    void OutputVTU( Model<3>& model,
                 const char* fileName,
                 long time );

    void SetFront( Model<3>& model,
                     double satOil,
                     double xFront );

    void Test( CFL_TestData& data, double satOil, double tolFront );

};


} // csmp

#endif // CFL_TESTCASE_H
