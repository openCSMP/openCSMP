// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  PostProcessor.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 6/8/20.
//

#include "PostProcessor.h"
#include "Region.h"
#include "ImplicitTransport.h"

using namespace std;

namespace csmp {

template<uint32_t dim, template<uint32_t> class USER>
void PostProcessor<dim,USER>::PostProcess()
 {
    for ( typename GoverningEquation<dim>::PostProcessingOperatorConstIterator 
          it=User()->PostProcessingOperatorsBegin(); it!=User()->PostProcessingOperatorsEnd(); ++it )
      {
         // loop over the domain, apply the operator, store the results
         for ( auto ot=User()->ComputationDomain().CellsBegin(); ot!=User()->ComputationDomain().CellsEnd(); ++ot )
           (*ot)->Out();
      }
          
 }

template class PostProcessor<3U,ImplicitTransport>;

} // end csmp
