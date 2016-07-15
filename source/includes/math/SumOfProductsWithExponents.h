#ifndef CSMP_SUM_OF_PRODUCTS_WITH_EXPONENTS_H
#define CSMP_SUM_OF_PRODUCTS_WITH_EXPONENTS_H

#include "CSMP_definitions.h"

/**

proxy for chemical reactions (mass action laws) with stochiometric
 exponents

CONVENTION: indepent species terms have negative indices and they are not differentiated
            for in derivative*/

namespace csmp {

//          c.idx/coefficient    e.idx/exponent
typedef std::pair<std::pair<int32,double64>,std::map<int32,std::pair<double64,double64> > >  sumofproducts_pair;

class SumOfProductsWithExponents {
   // terms that make up sum, e.g.,
   // a1 c2 - b2^2 c1 ...
   // terms    coeff   spec.idx  exponents   
   std::vector<sumofproducts_pair>  sumproducts;
   
 public:  
   SumOfProductsWithExponents();
   ~SumOfProductsWithExponents();
   SumOfProductsWithExponents( const SumOfProductsWithExponents& s );
   SumOfProductsWithExponents& operator=(const SumOfProductsWithExponents& s );
   
   void       AddProduct( int32 coeff_idx, double64 coeff, const std::map<int32,std::pair<double64,double64> >& product );
   double64  Sum( double64 a[] ) const;
   double64  DerivativeWithRespectTo( int32 spec_idx, double64 a[] ) const;
   void       Erase();
   void       Out() const;
   double64  Sign( double64 a ) const { return a >= 0.0 ? 1.0 : -1.0; };

   std::vector<sumofproducts_pair>::iterator  Begin() { return sumproducts.begin(); };
   std::vector<sumofproducts_pair>::iterator  End()   { return sumproducts.end(); };
};


/**

@class SumOfProductsWithExponents SumOfProductsWithExponents "applied_math/SumOfProductsWithExponents.h"

- Sum of products must has storage for constant product terms and their exponents. The indices
  for such constant terms should be set to <0, since they are not dependent variables.
const term, exponent

@code
idea: pair<pair<int,double64>,map<int,pair<double64,double64>,less<int> > >  sumofproducts_pair;
@endcode

when derivatives are taken and the spec.index, idx < 0, the constants are used rather than
the values supplied through a[].
*/

} // end namespace csp


#endif




