//
//  LinearAlgebraicSystem.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 6/8/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "LinearAlgebraicSystem.h"

using namespace std;

namespace csmp {


LinearAlgebraicSystem::LinearAlgebraicSystem( size_t m_x_n )
 : LHS( m_x_n ),
   RHS( m_x_n ),
   X( m_x_n )
 {
 }



void LinearAlgebraicSystem::Resize( size_t m_x_n )
 { 
    LHS.Resize( m_x_n );
    RHS.resize( m_x_n );
    X.resize( m_x_n );
 }



} // end csmp
