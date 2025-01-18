// include declarations
#include "AP_BoolVector.h"
#include <memory>

namespace csmp {

/** Size constructor.

@param uSize initial size of the bit vector.
*/
BoolVector::BoolVector(size_t uSize) :
 pBits_ (NULL),
 uSize_ (uSize) 
{
  if( uSize != 0 ) //if not, it will allocate 1 byte although empty
  {
   pBits_ = new uint8_t[uSize/8+1];
  
   // uSize/8+1 is the byte size
   // initialize everything to false   
   memset(pBits_,0, uSize/8+1);
  } 
}

/** Copy constructor.

 */
BoolVector::BoolVector(const BoolVector& bv) :
 pBits_ (NULL),
 uSize_ (bv.uSize_) 
{
  if( uSize_ != 0 ) //special case
  {
   pBits_ = new uint8_t[uSize_/8+1];
  
   // copy bits
   memcpy(pBits_,bv.pBits_,uSize_/8+1);
  }   
}

/** Assign operator

 */
BoolVector& BoolVector::operator=( const BoolVector& bv )
{
  if ( &bv != this ) {
      uSize_ = bv.uSize_; 
      if( pBits_ != NULL ) delete [] pBits_;
      if( uSize_ == 0 ) pBits_ = NULL;
      else   
        {
           const size_t uHowManyBytes ( bv.uSize_/8+1 );
           pBits_ = new uint8_t[uHowManyBytes];
           memcpy(pBits_,bv.pBits_,uHowManyBytes);
        }
    }
  return *this;
}



/** Default constructor.

*/
BoolVector::BoolVector() :
 pBits_ ( NULL ),
 uSize_ ( 0U )
{
}

/** destructor.
 */
BoolVector::~BoolVector()
{
  delete[] pBits_;
}



size_t BoolVector::Size() const
{
  return uSize_;
}

/** Sets the value of the bit with the specified offset.

@param uOffset offset of the bit in the bit vector.
@param bValue  logical value to set the bit.
 */
void BoolVector::SetBit(const size_t uOffset, const bool bValue)
{
  (bValue != 0) ? *(pBits_ + (uOffset >> 3)) |=  (1 << (uOffset & 0x7)) : 
    *(pBits_ + (uOffset >> 3)) &= ~(1 << (uOffset & 0x7));                 
}

/** Gets the value of the bit with the specified offset.

@param uOffset offset of the bit in the bit vector.
@return logical value of the bit.
 */
bool BoolVector::GetBit(const size_t uOffset) const
{
  return ((*(pBits_ + (uOffset >> 3)) & (1 << (uOffset & 0x7))) != 0);
}


/** 

@param uSize initial size of the bit vector.
*/

void BoolVector::Resize(size_t uSize)
{
  if( uSize == 0 ) //special case
  {
    uSize_ = 0;
    
    if(pBits_ != NULL)
     delete[] pBits_;
    
    pBits_ = NULL;
    
    return;
  }  
  
  const size_t uByteSize(uSize/8+1); 
  const size_t uCurrentByteSize(uSize_/8+1); 
  
  //remember old bit container ptr
  uint8_t * pBits_Tmp = pBits_;
  
  //compute how many to copy
  const size_t uHowMany(std::min(uCurrentByteSize,uByteSize));
  
  //resize bit container
  pBits_ = new uint8_t[uByteSize];
  
  //initiallize to false
  memset(pBits_,0,uByteSize);

  // copy old info that fits & delete temp
  if (pBits_Tmp != NULL)
  {
    memcpy(pBits_,pBits_Tmp,uHowMany);
    delete[] pBits_Tmp;
  }
  
  //set new size
  uSize_=uSize;
}



/** Sets all bits to a value.
*/
void BoolVector::SetAll(const bool value)
{
  if (pBits_!=NULL)
  {
    const size_t uByteSize(uSize_/8+1);
    memset(pBits_,(value?255:0),uByteSize);
  }
}



void BoolVector::Out() const
 {
    std::cout <<"\nBoolVector::Out: not defined yet."<< std::endl;
 }

} //end of namespace csmp
