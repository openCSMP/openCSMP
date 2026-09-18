// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

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

double SumOfProductsWithExponents::Sum( double a[] ) const
 {
    std::vector<sumofproducts_pair>::const_iterator                it;
    std::map<int,std::pair<double,double> >::const_iterator  sit;
    double  sum, product;
    
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
double SumOfProductsWithExponents::DerivativeWithRespectTo( int32_t spec_idx, double a[] ) const
 {
    std::vector<sumofproducts_pair>::const_iterator                it;
    std::map<int32_t,std::pair<double,double> >::const_iterator  sit;
    double  derivative, product;
 
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
void SumOfProductsWithExponents::AddProduct( int32_t coeff_idx, double coeff, const map<int32_t,pair<double,double> >& product )
 {
     sumproducts.reserve( sumproducts.size()+1 );
     sumproducts.push_back( sumofproducts_pair(make_pair(coeff_idx,coeff),product) );
 }



void SumOfProductsWithExponents::Erase()
 {
     sumproducts.erase( sumproducts.begin(), sumproducts.end() );
 }



void SumOfProductsWithExponents::Out() const
 {
    size_t n(1U);
    
    cout <<"\nSumOfProductsWithExponents:: = ";
    for ( vector<sumofproducts_pair>::const_iterator
          it=sumproducts.begin(); it!=sumproducts.end(); it++, n++ )
      {
         cout << (*it).first.second <<"(c"<< (*it).first.first <<") ";
         for ( map<int32_t,pair<double,double> >::const_iterator
               sit=(*it).second.begin(); sit!=(*it).second.end(); sit++ )
           {
              if ( (*sit).first >= 0 )
                {
                   if ( (*sit).second.first < 0.0 ) cout <<"-y"<< (*sit).first;
                   else                             cout <<"y"<< (*sit).first;
                   if ( (*sit).second.second != 1.0 ) cout <<"^"<< (*sit).second.second;
                   cout <<" ";
                }
              else 
                {
                   cout << (*sit).second.first;
                   if ( (*sit).second.second != 1.0 ) cout <<"^"<< (*sit).second.second;
                   cout <<" ";
                }
           }
         if ( n < sumproducts.size() ) cout <<"+ ";
      }
    cout << endl;
    
 } // end out

} // end namespace csmp







