#include "ScalarVariable.h"

namespace csmp {

void  ScalarVariable::Out(std::ostream& os) const
 {
    os <<"\nStatus: "<< parseStatus(flag_);
    if ( isnan(data_) )
      os <<", value: NAN\n";
    else
      os <<", value: " << data_ << std::endl;
 }

 

std::ostream&  operator<<( std::ostream& stream, const ScalarVariable& o )
 {
    stream << o.Value() <<" ("<< parseStatus(o.Flag()) <<")";
    return stream;
 }


} // end namespace csmp
