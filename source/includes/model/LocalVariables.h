#ifndef LOCAL_VARIABLES_H
#define LOCAL_VARIABLES_H

namespace csmp {

/// Data class to store physical variable count at given placement
struct LocalVariables {
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

  LocalVariables( uint32_t scalarsVars,
                  uint32_t vectorVars,
                  uint32_t tensorVars,
                  uint32_t array_count,
                  uint32_t array_length,
                  uint32_t flag_array_count,
                  uint32_t flag_array_length,
                  uint32_t total_data_depth,
                  uint32_t total_flag_depth )

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

  uint32_t  scalars,
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
