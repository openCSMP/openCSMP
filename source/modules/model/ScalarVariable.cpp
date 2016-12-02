#include "ScalarVariable.h"

namespace csmp {

void  ScalarVariable::Out() const
 {
    std::cout <<"\nStatus: "<< parseStatus(flag_);
    if ( isnan(data_) )
      std::cout <<", value: NAN\n";
    else
      std::cout <<", value: " << data_ << std::endl;
 }

 

std::ostream&  operator<<( std::ostream& stream, const ScalarVariable& o )
 {
    stream << o.Value() <<" ("<< parseStatus(o.Flag()) <<")";
    return stream;
 }


} // end namespace csmp
