#ifndef BOOL_VECTOR_H
#define BOOL_VECTOR_H

// include declarations

#include "CSMP_definitions.h"

/** Provides functionality to handle a bit vector.
@author A Paluszny
@date 05 2007
 */
namespace csmp {

class BoolVector
{
public:
  BoolVector();
  explicit BoolVector(size_t uSize);
  BoolVector(const BoolVector & bv);
  BoolVector& operator=( const BoolVector& );
  ~BoolVector();
  void Resize(size_t uSize);
  size_t Size() const;

  /// sets the value of the bit with the specified offset. (uOffset: 0-> uSize_-1)
  void SetBit( size_t uOffset, bool bValue );

  /// sets the value of the bit with the specified offset. (uOffset: 0-> uSize_-1)
  bool GetBit( size_t uOffset ) const;

  /// sets all bits to a value
  void SetAll( bool value );
  
  void Out() const { Out(std::cout); }
  void Out(std::ostream& os) const;

private:
  /// bit vector
  uint8* pBits_;

  /// size of the vector in bytes
  size_t uSize_;
};

}

#endif //BoolVector_h


 

 

