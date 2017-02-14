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



/** Default constructor.

*/
inline BoolVector::BoolVector() :
 pBits_ ( NULL ),
 uSize_ ( 0U )
{
}

/** destructor.
 */
inline BoolVector::~BoolVector() 
{
  delete[] pBits_;
}



inline size_t BoolVector::Size() const
{
  return uSize_;
}

/** Sets the value of the bit with the specified offset.

@param uOffset offset of the bit in the bit vector.
@param bValue  logical value to set the bit.
 */
inline void BoolVector::SetBit(const size_t uOffset, const bool bValue)
{
  (bValue != 0) ? *(pBits_ + (uOffset >> 3)) |=  (1 << (uOffset & 0x7)) : 
    *(pBits_ + (uOffset >> 3)) &= ~(1 << (uOffset & 0x7));                 
}

/** Gets the value of the bit with the specified offset.

@param uOffset offset of the bit in the bit vector.
@return logical value of the bit.
 */
inline bool BoolVector::GetBit(const size_t uOffset) const
{
  return ((*(pBits_ + (uOffset >> 3)) & (1 << (uOffset & 0x7))) != 0);
}

} 

#endif //BoolVector_h


 

 

