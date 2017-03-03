#ifndef CSMP_SUM_OF_PRODUCTS_H
#define CSMP_SUM_OF_PRODUCTS_H

#include "CSMP_definitions.h"

namespace csmp {

class SumOfProducts {
 public:  
   SumOfProducts();
   ~SumOfProducts();
   SumOfProducts( const SumOfProducts& s );
   SumOfProducts& operator=(const SumOfProducts& s );
   
// negative indexed terms are regarded as constants 
   void       AddProduct( const std::map<int32,double64>& product );
   double64  Sum( const double64 a[] ) const;
   void       Erase();

   void       InitializeFrom( const char* s, 
                              std::map<std::string,std::pair<int,double64> >& withval,
                              std::map<std::string,int>& withoutval );

   void       TokenizeString( char* s, const char* delim, std::list<std::string>& li ) const;

   std::vector<std::map<int,double64> >::iterator  Begin() { return sumproducts.begin(); };
   std::vector<std::map<int,double64> >::iterator  End()   { return sumproducts.end(); };

   void       Out() const;
 
 private:
   // terms that make up sum, e.g.,
   // a1 c2 - b2^2 c1 ...
   // terms:  spec.idx  product term   
   std::vector<std::map<int,double64> >  sumproducts;
};

} // end namespace csp

#endif


