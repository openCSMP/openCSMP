#ifndef FIXED_ARRAY_H
#define FIXED_ARRAY_H

#include <iostream>
#include <vector>

namespace csmp {

#define FIXED_ARRAY_debug

template<unsigned int dim>
class FixedArray {
    public:
      FixedArray() {}
      ~FixedArray() {}
      explicit FixedArray( const double64& a ) {
           for ( unsigned int i=0U; i<dim; i++ ) v_[i] = a;
        }
      FixedArray( const FixedArray& a ) {
           for ( unsigned int i=0U; i<dim; i++ ) v_[i] = a.v_[i];
        }
      FixedArray& operator=( const std::vector<double64>& vec ) { 
#ifdef FIXED_ARRAY_debug
if ( dim != vec.size() ) std::cout<<"\nFixedArray=vector: size mismatch."<< std::endl;
#endif      
           for ( unsigned int i=0U; i<dim; i++ ) v_[i]=vec[i];
           return *this;   
        }
      FixedArray& operator=( const double64& val ) {
           for ( unsigned int i=0U; i<dim; i++ ) v_[i]=val;
           return *this;   
        }  
      double64  operator[]( unsigned int i ) const { return v_[i]; }
      double64& operator[]( unsigned int i ) { return v_[i]; }
      unsigned int size() const { return dim; }
      void Out( std::vector<double64>& vec ) const {
           vec.resize(dim); for ( unsigned int i=0U; i<dim; i++ ) vec[i]=v_[i]; 
        }  
      void Out() const { Out(std::cout); }
      void Out(std::ostream& os) const { 
           std::cout << "\ncsp::FixedArray<"<< dim <<" ";
           for ( unsigned int i=0U; i<dim; i++ ) 
             std::cout <<v_[i]; std::cout << std::endl; 
        }
      // substitute call "fT* ary" with  &ary[0]
    private:
      double64 v_[dim];
 };
 
 
 
 } // end namespace csmp

#endif
