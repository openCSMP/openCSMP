#ifndef MAXIMUM_DIFFERENCE_H
#define MAXIMUM_DIFFERENCE_H

#include "Interrelation.h"

namespace csmp {

template<uint32_t dim>
class MaximumDifference : public Interrelation<dim> {
    Operand<dim>&   Var1; /// < total fluid pressure
    Operand<dim>&   Var2; /// < conductivity multiplier
    ScalarVariable  var1, var2;
    double        max_difference;

  public:
    explicit MaximumDifference( const PropertyDatabase<dim>& p,
                                const char* variable1, const char* variable2 );

    void Calculate() override final;

    void   Reset();
    double Value() const; /// < outputs the computed value
};

/// compare distributed variable values: calculates maximum difference between distributed variable values of the same type and placement; difference can be normalised
template<uint32_t dim>
double maximumDifference( const Model<dim>&,
                            const char* new_property, const char* old_property,
                            bool normalise );



} // csmp

#endif

