// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "SimulatorSetupParameter.h"

namespace csmp {

// comparison, not case sensitive. (from http://www.cplusplus.com/reference/list/list/sort/)
bool compare_setup_parameter_nocase (const csmp::SimulatorSetupParameter& first, const csmp::SimulatorSetupParameter& second)
{
    unsigned int i=0;
    while ( (i<first.name.length()) && (i<second.name.length()) )
    {
        if (tolower(first.name[i])<tolower(second.name[i])) return true;
        else if (tolower(first.name[i])>tolower(second.name[i])) return false;
        ++i;
    }
    return ( first.name.length() < second.name.length() );
}

} // end csmp

