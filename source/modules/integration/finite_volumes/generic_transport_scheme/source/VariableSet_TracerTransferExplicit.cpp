//
//  VariableSet_TracerTransferExplicit.cpp
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#include "VariableSet_TracerTransferExplicit.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp {

template<size_t dim>
VariableSet_TracerTransferExplicit::VariableSet_TracerTransferExplicit( const PropertyDatabase<dim>& p )
 : VariableSet_TracerTransfer(p)
 {
   // potentially these properties could be created here from scratch to have them only
   // when there is a need for a transport calculation.
 }



void VariableSet_TracerTransferExplicit::CheckVariables() const
 {
    VariableSet_TracerTransfer::CheckVariables();
// TODO: complete missing checks

} // end CheckVariables

 
template VariableSet_TracerTransferExplicit::VariableSet_TracerTransferExplicit( const PropertyDatabase<1U>& );
template VariableSet_TracerTransferExplicit::VariableSet_TracerTransferExplicit( const PropertyDatabase<2U>& );
template VariableSet_TracerTransferExplicit::VariableSet_TracerTransferExplicit( const PropertyDatabase<3U>& );
 
} // end csmp
