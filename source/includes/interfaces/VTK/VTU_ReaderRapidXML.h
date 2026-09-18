// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

//
//  read_VTU_File() function.
//
//  Open-CSMP
//
//  Created by Stephan Matthai on 26/1/2024.
//

#ifndef XML_VTU_RAPID_XML_BASED_READER_H
#define XML_VTU_RAPID_XML_BASED_READER_H

#include "CSMP_definitions.h"

namespace csmp {

class ScalarVariable;
class ArrayVariable;
class FlaggedArrayVariable;
template<uint32_t> class VectorVariable;
template<uint32_t> class TensorVariable;
template<uint32_t> class VSet;
class ModelTopology;


 /**
     Reads the supplied VTK (XML=.vtu) file returning the mesh and associated properties into the VSet.
     Model regions are inferred from 'attribute' or 'region id' named scalar cell variables.
     The region information is returned into a ModelTopology object.
     Method returns error codes or 0 if everything went OK.
     
     @param fname name of input file excluding the extension '.vtu.
     
     @param verbose prints entire data records retrieved from file.
     
     @return returns error code or zero if all data records were read successfully.
 */
 template<uint32_t dim>
 int read_VTU_File( const char* fname, VSet<dim>&, ModelTopology&, bool verbose=false );
     




// HELPER FUNCTIONS

/// Compile-time version of member Size(), given space dimension of the model, returns variable values that are stored in the template parameter variable.
template<uint32_t dim,template<uint32_t> class csmp_var>
int valuesInType( int components = -1 ) {
     if constexpr      ( std::is_same<csmp_var<dim>,csmp::ScalarVariable>::value_type ) return 1;
     else if constexpr ( std::is_same<csmp_var<dim>,csmp::VectorVariable<dim>>::value_type ) return dim;
     else if constexpr ( std::is_same<csmp_var<dim>,csmp::TensorVariable<dim>>::value_type ) return dim * dim;
     // Array and flagged array need a different function because they cannot be evaluated at compile time
     else if ( components >= dim ) return components;
     return csmp::UNSPECIFIED;
  }
  
/// Compile-time version of member Size() returns number of values in ScalarVariable and ArrayVariable types.
template<typename csmp_var>
int valuesInType( const csmp_var& var = csmp_var{} ) {
     if constexpr      ( std::is_same<csmp_var,csmp::ScalarVariable>::value_type ) return 1;
     else if constexpr ( std::is_same<csmp_var,csmp::ArrayVariable>::value_type ) return var.Size();
     else if constexpr ( std::is_same<csmp_var,csmp::FlaggedArrayVariable>::value_type ) return var.Size();
     else if constexpr ( std::is_same<csmp_var,double>::value_type ) return 1;
     else if constexpr ( std::is_same<csmp_var,int>::value_type ) return 1;
     return csmp::UNSPECIFIED;
  }
  
/// Variable type deduction from number of components; functions does does return flagged array type.
template<uint32_t dim>
csmp::VARIABLE_TYPE inferTypeFromNumberOfComponents( int components ) {
     static_assert( dim != 1U, "inferFromNumberOfComponents: impossible task in 1D");
     if      ( components == 1 ) return csmp::SCALAR;
     else if ( components == dim ) return csmp::VECTOR;
     else if ( components == dim * dim ) return csmp::TENSOR;
     else if ( components > dim ) return csmp::ARRAY;
     std::cerr <<"\n\n"<<"ERROR, inferTypeFromNumberOfComponents("<< components <<"): failed to identify csmp::VARIABLE_TYPE, returning ARRAY"<< std::endl;
     return csmp::ARRAY;
  }

} // end CSMP



#endif /* XML_VTU_RAPID_XML_BASED_READER_H */
