#ifndef __ModifyPropertyWhere_h__
#define __ModifyPropertyWhere_h__

#include "Interrelation.h"

namespace csmp {

template<size_t dim>
class ModifyPropertyWhere : public Interrelation<dim> {
    Operand<dim>&  k;
    double64                if_value, to_value;
    
  public:
    ModifyPropertyWhere( const PropertyDatabase<dim>& p, 
                         const char* prop, double64 where, double64 to );
                         
    ~ModifyPropertyWhere() {};
    void Calculate();
};


template<size_t dim>
inline void ModifyPropertyWhere<dim>::Calculate()
 {
    if ( k == if_value ) k = to_value;
 } 

} // csmp

#endif

