#include "SumOfProductsWithExponents.h"

#include <cmath>

using namespace std;

namespace csmp {

SumOfProductsWithExponents::SumOfProductsWithExponents()  
{}



SumOfProductsWithExponents::~SumOfProductsWithExponents() 
{}



SumOfProductsWithExponents::SumOfProductsWithExponents( const SumOfProductsWithExponents& s )  
{ *this=s; }



SumOfProductsWithExponents& SumOfProductsWithExponents::operator=(const SumOfProductsWithExponents& s )  
 { 
    if ( &s==this ) return *this; 
    sumproducts=s.sumproducts;
    return *this; 
 }


// Methods 

double64 SumOfProductsWithExponents::Sum( double64 a[] ) const
 {
    std::vector<sumofproducts_pair>::const_iterator                it;
    std::map<int,std::pair<double64,double64> >::const_iterator  sit;
    double64  sum, product;
    
    for ( sum=0.0, it=sumproducts.begin(); it!=sumproducts.end(); it++ )
      {
         for ( product=1.0, sit=(*it).second.begin(); sit!=(*it).second.end(); sit++ )
            {
               if ( (*sit).first >= 0 ) 
                 {
                    if ( (*sit).second.second == 1.0 ) product *= a[ (*sit).first+1 ] * Sign((*sit).second.first);
                    else product *= std::pow( a[ (*sit).first+1 ]*Sign((*sit).second.first), (*sit).second.second );
                 }
               else product *= std::pow( (*sit).second.first, (*sit).second.second );
            }
         sum += (*it).first.second * product;
      }
    return sum;      
 }



/**

go through all products and enlist all terms that contain requested species index while
omitting respective species from product. 'a' is the vector that contains the multiplication factors for
each component of the product */
double64 SumOfProductsWithExponents::DerivativeWithRespectTo( int32 spec_idx, double64 a[] ) const
 {
    std::vector<sumofproducts_pair>::const_iterator                it;
    std::map<int32,std::pair<double64,double64> >::const_iterator  sit;
    double64  derivative, product;
 
    assert( spec_idx >= 0 );
 
    // calculating the non-zero contributions
    for ( derivative=0.0, it=sumproducts.begin(); it!=sumproducts.end(); it++ )
      {
         // if the desired species is part of the product
         // it is calculated omitting the species term unless its exponent is not equal to 1
         // in which case it is differentiated.
         if ( ((*it).second.find( spec_idx )) != (*it).second.end() )
           {
              for ( product=1.0, sit=(*it).second.begin(); sit!=(*it).second.end(); sit++ )
                if ( (*sit).first != spec_idx )
                  {
                     // dependent variable have idx >=0 and are taken from supplied array
                     if ( (*sit).first >= 0 ) product *= std::pow( a[ (*sit).first+1 ]*Sign((*sit).second.first), (*sit).second.second );
                     // for independent variables the stored value is taken
                     else                     product *= std::pow( (*sit).second.first, (*sit).second.second );
                  }
                // if (*sit).first == spec_idx and species to be differentiated after has an 
                // exponent other than 1 it is differentiated and summed (f(x^2)=2x^2-1).
                // (derivatives can only be taken of dependent variables!)   
                else 
                if ( (*sit).first >= 0 && (*sit).second.second != 1.0 )
                  product *= ((*sit).second.second * std::pow( a[ (*sit).first+1 ]*Sign((*sit).second.first), (*sit).second.second-1.0 ));
              //                   prod.coeff    product 
              derivative += (*it).first.second * product;
           }
      }
    return derivative;

 } // end DerivativeWithRespectTo






/**

const valued product terms must have negative exponents, such that when derivatives are taken they are not
input from the supplied a[] vector. */
void SumOfProductsWithExponents::AddProduct( int32 coeff_idx, double64 coeff, const map<int32,pair<double64,double64> >& product )
 {
     sumproducts.reserve( sumproducts.size()+1 );
     sumproducts.push_back( sumofproducts_pair(make_pair(coeff_idx,coeff),product) );
 }



void SumOfProductsWithExponents::Erase()
 {
     sumproducts.erase( sumproducts.begin(), sumproducts.end() );
 }



void SumOfProductsWithExponents::Out(std::ostream& os) const
 {
    size_t n(1U);
    
    os <<"\nSumOfProductsWithExponents:: = ";
    for ( vector<sumofproducts_pair>::const_iterator
          it=sumproducts.begin(); it!=sumproducts.end(); it++, n++ )
      {
         os << (*it).first.second <<"(c"<< (*it).first.first <<") ";
         for ( map<int32,pair<double64,double64> >::const_iterator
               sit=(*it).second.begin(); sit!=(*it).second.end(); sit++ )
           {
              if ( (*sit).first >= 0 )
                {
                   if ( (*sit).second.first < 0.0 ) os <<"-y"<< (*sit).first;
                   else                             os <<"y"<< (*sit).first;
                   if ( (*sit).second.second != 1.0 ) os <<"^"<< (*sit).second.second;
                   os <<" ";
                }
              else 
                {
                   os << (*sit).second.first;
                   if ( (*sit).second.second != 1.0 ) os <<"^"<< (*sit).second.second;
                   os <<" ";
                }
           }
         if ( n < sumproducts.size() ) os <<"+ ";
      }
    os << endl;
    
 } // end out

} // end namespace csmp







