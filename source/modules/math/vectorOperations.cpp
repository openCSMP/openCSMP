// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  vectorOperations.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 8/02/2017.
//

#include "vectorOperations.h"

using namespace std;

namespace csmp {

template<uint32_t dim>
double valueAverage( const VectorVariable<dim>& vc )
 {
     double sum(0.);
     for ( auto i{0U}; i<dim; ++i )
       sum += vc[i];
   
     return sum / static_cast<double>(dim);
 }

template double valueAverage( const VectorVariable<1U>& );
template double valueAverage( const VectorVariable<2U>& );
template double valueAverage( const VectorVariable<3U>& );

} // end csmp
