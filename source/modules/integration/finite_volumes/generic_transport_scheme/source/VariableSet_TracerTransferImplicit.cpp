//
//  VariableSet_TracerTransferImplicit.cpp
//
//  Created by Stephan Matthai on 19/06/2015.
//  Copyright (c) 2015 Stephan K. Matthai. All rights reserved.
//

#include "VariableSet_TracerTransferImplicit.h"
#include "Exception.h"
#include "PropertyDatabase.h"

namespace csmp {

template<size_t dim>
VariableSet_TracerTransferImplicit::VariableSet_TracerTransferImplicit( const PropertyDatabase<dim>& p )
    : VariableSet_TracerTransfer(p)
{
}


 
void VariableSet_TracerTransferImplicit::CheckVariables() const
 {
     VariableSet_TracerTransfer::CheckVariables();
} // end CheckVariables


 
template VariableSet_TracerTransferImplicit::VariableSet_TracerTransferImplicit( const PropertyDatabase<1U>& );
template VariableSet_TracerTransferImplicit::VariableSet_TracerTransferImplicit( const PropertyDatabase<2U>& );
template VariableSet_TracerTransferImplicit::VariableSet_TracerTransferImplicit( const PropertyDatabase<3U>& );
 
} // end csmp
