//
//  SlopeMechanics_Example.h
//  examples
//
//  Created by Stephan Matthai on 12/10/18.
//

#ifndef CSMP_SLOPE_MECHANICS_EXAMPLE_H
#define CSMP_SLOPE_MECHANICS_EXAMPLE_H

#include "Example.h"

#define DIM 2

namespace csmp {

template<size_t> class Model;

class SlopeMechanics_Example : public Example {
  public:
    virtual ~SlopeMechanics_Example() {}
    virtual void Specifications();
    virtual void Run();
};

void porePressureBiotAlphaProduct( Model<DIM>& model, const char* target_region );

void dilatationInducedChangeInPorePressure( Model<DIM>& model );

void permeabilityPorosityCorrelation( Model<DIM>& model );

void gravityForce( Model<DIM>& model, double acc_gravity );


} // end csmp

#endif /* CSMP_SLOPE_MECHANICS_EXAMPLE_H */
