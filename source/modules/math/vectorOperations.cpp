//
//  vectorOperations.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 8/02/2017.
//  Copyright © 2017 Stephan Matthai. All rights reserved.
//

#include "vectorOperations.h"

using namespace std;

namespace csmp {

template<size_t dim>
double64 valueAverage( const VectorVariable<dim>& vc )
 {
     double64 sum(0.);
     for ( size_t i=0U; i<dim; ++i )
       sum += vc[i];
   
     return sum / static_cast<double64>(dim);
 }

template double64 valueAverage( const VectorVariable<1U>& );
template double64 valueAverage( const VectorVariable<2U>& );
template double64 valueAverage( const VectorVariable<3U>& );

} // end csmp
