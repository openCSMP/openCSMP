#ifndef MAXIMUM_DIFFERENCE_H
#define MAXIMUM_DIFFERENCE_H

#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class MaximumDifference : public Interrelation<dim> {
    Operand<dim>&   Var1; /// < total fluid pressure
    Operand<dim>&   Var2; /// < conductivity multiplier
    ScalarVariable  var1, var2;
    double64        max_difference;

  public:
    explicit MaximumDifference( const PropertyDatabase<dim>& p,
                                const char* variable1, const char* variable2 );
    ~MaximumDifference() {};
    void Calculate();
    void Reset();
    double64   Value() const; /// < outputs the computed value
};




template<size_t dim>
inline void MaximumDifference<dim>::Calculate()
 {
    Var1.AssignTo( var1 );
    Var2.AssignTo( var2 );
    
    if ( std::fabs(var1() - var2()) > max_difference )
      max_difference = std::fabs(var1() - var2());
     
 } // end




} // csmp

#endif

