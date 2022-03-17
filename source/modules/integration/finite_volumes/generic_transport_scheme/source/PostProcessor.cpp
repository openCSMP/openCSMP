//
//  PostProcessor.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 6/8/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
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
         for ( auto ot=User()->ComputationDomain().ElementsBegin(); ot!=User()->ComputationDomain().ElementsEnd(); ++ot )
           (*ot)->Out();
      }
          
 }

template class PostProcessor<3U,ImplicitTransport>;

} // end csmp
