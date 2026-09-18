// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  convertANSYS_to_CSMP_model.h
//  CSMP_DFM_Upscaling
//
//

#ifndef convertANSYS_to_CSMP_model_H
#define convertANSYS_to_CSMP_model_H

#include <iostream>

namespace csmp {

  /// read in the *.asc, *.dat, *-regions.txt and *-configuration.txt files and build native CSMP model from them
  bool convertANSYS3D_to_CSMP_model( const char* ansys_model_input_deck, const char* var_file_name );
  bool convertANSYS2D_to_CSMP_model( const char* ansys_model_input_deck, const char* var_file_name );

}



#endif 
