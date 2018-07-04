#include "SumOfProducts.h"

using namespace std;

namespace csmp {

SumOfProducts::SumOfProducts()  
{}


SumOfProducts::~SumOfProducts() 
{}

SumOfProducts::SumOfProducts( const SumOfProducts& s )  
{ *this=s; }



SumOfProducts& SumOfProducts::operator=(const SumOfProducts& s )  
 { 
    if ( &s==this ) return *this; 
    sumproducts=s.sumproducts;
    return *this; 
 }


// Methods 


void SumOfProducts::AddProduct( const map<int,double64>& product )
 {
     sumproducts.reserve( sumproducts.size()+1 );
     sumproducts.push_back( product );
 }




void SumOfProducts::TokenizeString( char* s, const char* delim, list<string>& li ) const
  {
     li.erase( li.begin(), li.end() );
     char* tok;
     tok = strtok( s, delim );
     do {
           li.push_back( string(tok) );
        }
     while ( (tok=strtok( NULL, delim )) != NULL ); 
  }




/**

takes character string containing the righthandside of the Algebraic equation and builds
itself from this string.*/
void SumOfProducts::InitializeFrom( const char* s, map<string,pair<int32,double64> >& withval,
                                                   map<string,int32>&              withoutval )
 {
    list<string>                                        products, product_terms;
    list<string>::iterator                              tik, tik1;
    map<string,int32>::const_iterator                   it;
    map<string,pair<int32,double64> >::const_iterator  vit;
    map<int32,double64>                                    product;
    char                        rhs[250], temp[100], temp1[50], temp2[50], *sub;
    int32                         idx, idxminus=-500; // large enough to avoid overwrite by independent vars.
    double64                      val;
    strcpy( rhs, s );
    TokenizeString( rhs, "+", products );

     // reading products
     // ----------------
     for ( tik=products.begin(); tik!=products.end(); tik++ )  
       {
          // for each product create subtokens corresponding to terms
          // --------------------------------------------------------
          strcpy( temp, (*tik).c_str() );
          TokenizeString( temp, " ", product_terms );

          // reading product terms (term by term)
          // ---------------------------------------------------------
          for ( tik1=product_terms.begin(); tik1!=product_terms.end(); tik1++ ) 
            {
               // parsing product term (if an identifier was used instead of a number
               if ( (vit=withval.find( *tik1 )) != withval.end() )
                 {
                    idx = (*vit).second.first;
                    val = (*vit).second.second;
                 }
               else if ( (it=withoutval.find( *tik1 )) != withoutval.end() )
                 {
                    idx = (*it).second;
                    val = 1.0;
                 }
               // if product term is not found in the species list, 
               // the term is checked for whether it is preceded by
               // a minus sign which needs to be separated
               else
                 {
                    // return pointer to first occurrence of "-" in String
                    strcpy( temp1, (*tik1).c_str() );
                    sub = strchr( temp1, '-' );
                    if ( sub != NULL )
                      {
                         // reading species name
                         for ( size_t i=1U; i<=strlen(temp1); i++ ) temp2[i-1U] = temp1[i];
                         if ( (vit=withval.find( *tik1 )) != withval.end() )
                           {
                              idx = (*vit).second.first;
                              val = (*vit).second.second;
                            }
                         else if ( (it=withoutval.find(string(temp2))) != withoutval.end() )
                           {
                              idx = (*it).second;
                              val = -1.0;
                           }
                       }
                     else // if there is just a number given
                       {
                          val = atof( temp1 );
                          assert( val > -1.0e+100 && val < 1.0e+100 );
                          idx = idxminus--; 
                       }
                  }
               // recording product terms
               // -----------------------
               product[ idx ] = val;

            } // end for 

          AddProduct( product );
          product.erase( product.begin(), product.end() );
       }
   
 } // end InitializeFrom







/**

Add all the terms in the sum, using the supplied values in 'a' if the product term index
is greater than zero and the stored terms (constants) where the index is less than
zero.*/
double64 SumOfProducts::Sum( const double64 a[] ) const
 {
    vector<map<int32,double64> >::const_iterator  it;
    map<int32,double64>::const_iterator           sit;
    double64  sum, product;
    
    for ( sum=0.0, it=sumproducts.begin(); it!=sumproducts.end(); it++ )
      {
         for ( sit=(*it).begin(), product=1.0; sit!=(*it).end(); sit++ )
            {
               //                           a is indexed from 1...n      sign is used
               if ( (*sit).first >= 0 ) product *= a[ (*sit).first+1 ] * (*sit).second;
               else                     product *= (*sit).second;
            }
         sum += product;
      }
    return sum;      
 }




void SumOfProducts::Erase()
 {
     sumproducts.erase( sumproducts.begin(), sumproducts.end() );
 }








void SumOfProducts::Out() const
 {
    size_t n=1U;
    cout <<"\nSumOfProducts:: = ";
    for ( vector<map<int32,double64> >::const_iterator
          it1=sumproducts.begin(); it1!=sumproducts.end(); it1++, n++ )
      {
         for ( map<int32,double64>::const_iterator sit=(*it1).begin(); sit!=(*it1).end(); sit++ )
           {
              if ( (*sit).first >= 0 ) 
                {
                   if ( (*sit).second > 0. ) cout <<"a"<< (*sit).first <<" ";
                   else                      cout <<"-a"<< (*sit).first <<" ";
                }
              else cout << (*sit).second <<" ";
           }
         if ( n < sumproducts.size() ) cout <<"+ ";
      }
    cout << endl;
    
 } // end out


} // end namespace csmp






