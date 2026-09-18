// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef LOCAL_VARIABLES_H
#define LOCAL_VARIABLES_H

namespace csmp {

/// Data class to store physical variable count at given placement
struct LocalVariables {

  // TODO: investigate which type is best for performance
  using int_type = uint32_t; ///< unsigned integer type that is big enough to hold 'totalDataDepth'

  LocalVariables()
    : scalars           (0),
      vectors           (0),
      tensors           (0),
      arrayCount        (0),
      arrayLength       (0),
      flaggedArrayCount (0),
      flaggedArrayLength(0),
      totalDataDepth    (0),
      totalFlagDepth    (0)
  {}

  LocalVariables( int_type scalarsVars,
                  int_type vectorVars,
                  int_type tensorVars,
                  int_type array_count,
                  int_type array_length,
                  int_type flag_array_count,
                  int_type flag_array_length,
                  int_type total_data_depth,
                  int_type total_flag_depth )

    : scalars           (scalarsVars),
      vectors           (vectorVars),
      tensors           (tensorVars),
      arrayCount        (array_count),
      arrayLength       (array_length),
      flaggedArrayCount (flag_array_count),
      flaggedArrayLength(flag_array_length),
      totalDataDepth    (total_data_depth),
      totalFlagDepth    (total_flag_depth)
  {}

  LocalVariables( const LocalVariables& lvs ) 
    : scalars           (lvs.scalars),
      vectors           (lvs.vectors),
      tensors           (lvs.tensors),
      arrayCount        (lvs.arrayCount),
      arrayLength       (lvs.arrayLength),
      flaggedArrayCount (lvs.flaggedArrayCount),
      flaggedArrayLength(lvs.flaggedArrayLength),
      totalDataDepth    (lvs.totalDataDepth),
      totalFlagDepth    (lvs.totalFlagDepth)
  {}

  LocalVariables& operator=( const LocalVariables& lvs )
  { 
    if( this != &lvs ) 
    {
        scalars             = lvs.scalars;
        vectors             = lvs.vectors;
        tensors             = lvs.tensors;
        arrayCount          = lvs.arrayCount;
        arrayLength         = lvs.arrayLength;
        flaggedArrayCount   = lvs.flaggedArrayCount;
        flaggedArrayLength  = lvs.flaggedArrayLength;
        totalDataDepth      = lvs.totalDataDepth;
        totalFlagDepth      = lvs.totalFlagDepth;
    }
    return *this; 
  }

  bool Empty() const { return ( scalars==0U && vectors==0U && tensors==0U && arrayCount==0U && flaggedArrayCount==0U ); }

  int_type  scalars,
            vectors,
            tensors,
            arrayCount,
            arrayLength,
            flaggedArrayCount,
            flaggedArrayLength,
            totalDataDepth,
            totalFlagDepth;
};


} // csmp

#endif
