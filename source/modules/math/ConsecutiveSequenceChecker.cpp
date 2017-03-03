//
//  ConsecutiveSequenceChecker.cpp
//  CSMP_GitHub
//
//  Created by Stephan Matthai on 10/12/2016.
//  Copyright © 2016 Stephan Matthai. All rights reserved.
//

#include <algorithm>
#include <unordered_set>
#include <type_traits>
#include <cassert>
#include <ciso646>
#include <string>
#include "ConsecutiveSequenceChecker.h"

using namespace std;

/** checks elements of an unsigned integer set

1. Finds the Max­i­mum and min­i­mum ele­ments in array (Say the array is arrA)
2. Check if array length   = max-min+1
3. Sub­tract the min from every ele­ment of the array.
4. Check if array doesn’t have duplicates

*/
template<>
bool ConsecutiveSequenceChecker::IsRangeOfUnsignedIntConsecutive( const std::set<size_t>& iset,
                                                                  bool check_whether_max_value_is_size_minus1 )
 {
    if ( iset.empty() ) throw logic_error(" ConsecutiveSequenceChecker::IsRangeOfUnsignedIntConsecutive: set is empty.");
    // set is already sorted !
    // 1.
    // 2.
    if( (*iset.rbegin()) - (*iset.begin()) + 1U != iset.size() ) return false;
   
    if ( check_whether_max_value_is_size_minus1 )
      if( (*iset.rbegin()) + 1U != iset.size() ) return false;
   
    // SET DOES NOT CONTAIN DUPLICATES
    return true;
   
 } // end IsRangeOfUnsignedIntConsecutive



/** checks elements of an unsigned integer vector

1. Finds the max­i­mum and min­i­mum ele­ments in the array
2. Checks if array length = max-min+1
3. Sub­tract the min from every ele­ment of the array.
4. Checks if array doesn’t have duplicates

*/
template<>
bool ConsecutiveSequenceChecker::IsRangeOfUnsignedIntConsecutive( const std::vector<size_t>& ivec,
                                                                  bool check_whether_max_value_is_size_minus1 )
 {
    if ( ivec.empty() ) throw logic_error(" ConsecutiveSequenceChecker::IsRangeOfUnsignedIntConsecutive: vector is empty.");
 
    // 1. value ranges
    size_t imin = (*min_element( ivec.begin(), ivec.end() ));
    size_t imax = (*max_element( ivec.begin(), ivec.end() ));
   
    // 2. are the extrema consistent with a consecutive numbering
    if( imax - imin + 1U != ivec.size() ) return false;
   
    if ( check_whether_max_value_is_size_minus1 )
      if ( imax > 0U and imax - 1U != ivec.size() ) return false;
   
    // 3. checking for duplicates
    unordered_set<size_t> iset;
    for ( auto i=ivec.begin(); i!=ivec.end(); ++i ) {
        pair<unordered_set<size_t>::iterator,bool> it=iset.insert( (*i) );
        // fail on duplicates
        if ( !it.second ) return false;
      }
   
    // 4. checking whether numbering is consecutive
    vector<size_t>::const_iterator it1(ivec.begin()), it2(ivec.begin()); it2++;
    while ( it2 != ivec.end() ) {
         if ( (*it1) + 1U != (*it2) ) return false;
         it2++;
         it1++;
      }

    return true;
   
 } // end IsRangeOfUnsignedIntConsecutive




/**
   checks whether the keys of the map are numbered in consecutive order.
*/
template<typename K, typename V>
bool ConsecutiveSequenceChecker::IsKeyRangeOfUnsignedIntConsecutive( const std::map<K,V>& imap,
                                                                    bool check_whether_max_value_is_size_minus1 )
 {
    static_assert(is_integral<K>::value, "ConsecutiveSequenceChecker::IsKeyRangeOfUnsignedIntConsecutive: template argument on parameter must be an integer.");

    if ( imap.empty() ) throw logic_error(" ConsecutiveSequenceChecker::IsKeyRangeOfUnsignedIntConsecutive: map is empty.");
    // 1.
    // 2.
    if( imap.size() != (*imap.rbegin()).first - (*imap.begin()).first + 1U ) return false;

    if ( check_whether_max_value_is_size_minus1 )
      if( (*imap.rbegin()).first + 1U != imap.size() ) return false;
   
    return true;
   
 } // end IsKeyRangeOfUnsignedIntConsecutive

template bool ConsecutiveSequenceChecker::IsKeyRangeOfUnsignedIntConsecutive( const std::map<size_t,size_t>&, bool );



/**
    allows to check the values of the map (knowing that the key already is unique).
*/
template<>
bool ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive( const std::map<size_t,size_t>& imap,
                                                                       bool check_whether_max_value_is_size_minus1 )
 {
    //static_assert(is_integral<V>::value, "ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive: template argument on parameter must be an integer.");

    if ( imap.empty() ) throw logic_error(" ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive: map is empty.");
    // checking uniqueness first
    set<size_t> iset;
    for ( auto i=imap.begin(); i!=imap.end(); ++i ) {
        pair<set<size_t>::iterator,bool> it=iset.insert( (*i).second );
        // fail on duplicates
        if ( !it.second ) return false;
      }
   
    // checking whether range is consecutive
    if ( (*iset.rbegin()) - (*iset.begin()) + 1U != iset.size() ) return false;
   
    // checking whether largest element = number-of-elements - 1
    if ( check_whether_max_value_is_size_minus1 )
      if ( (*iset.rbegin()) + 1U != iset.size() ) return false;
   
    // checking whether value range is consecutive upon iteration over range
    map<size_t,size_t>::const_iterator it1(imap.begin()), it2(imap.begin()); it2++;
    while ( it2 != imap.end() ) {
         if ( (*it1).second + 1U != (*it2).second ) {
              cerr <<"\nConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive: detected problem (n vs. n+1): ";
              cerr << (*it1).second <<" vs "<< (*it2).second <<"\n";
              return false;
           }
         it2++;
         it1++;
      }   
   
    return true;
 
 } // end IsValueRangeOfUnsignedIntConsecutive




/**
    Checks whether the values in the supplied map are continuous, unique and bounded as 
    max value < #-elements in the map.
*/
template<>
bool ConsecutiveSequenceChecker::IsValueRangeUniqueAndBounded( const std::map<size_t,size_t>& imap,
                                                               bool check_whether_max_value_is_size_minus1 )
 {
    //static_assert(is_integral<V>::value, "ConsecutiveSequenceChecker::IsValueRangeUniqueAndBounded: template argument on parameter must be an integer.");

    if ( imap.empty() ) throw logic_error(" ConsecutiveSequenceChecker::IsValueRangeUniqueAndBounded: map is empty.");
    // checking uniqueness first
    set<size_t> iset;
    for ( auto i=imap.begin(); i!=imap.end(); ++i ) {
        pair<set<size_t>::iterator,bool> it=iset.insert( (*i).second );
        // fail on duplicates
        if ( !it.second ) return false;
      }
   
    // checking whether range is consecutive
    if ( (*iset.rbegin()) - (*iset.begin()) + 1U != iset.size() ) return false;
   
    // checking whether largest element = number-of-elements - 1
    if ( check_whether_max_value_is_size_minus1 )
      if ( (*iset.rbegin()) + 1U != iset.size() ) return false;
   
    // checking whether value range is consecutive upon iteration over range
    set<size_t>::const_iterator it1(iset.begin()), it2(iset.begin()); it2++;
    while ( it2 != iset.end() ) {
         if ( (*it1) + 1U != (*it2) ) {
              cerr <<"\nConsecutiveSequenceChecker::IsValueRangeUniqueAndBounded: detected problem (n vs. n+1): ";
              cerr << (*it1) <<" vs "<< (*it2) <<"\n";
              return false;
           }
         it2++;
         it1++;
      }
   
    return true;
 
 } // end IsValueRangeUniqueAndBounded






/*
    // checking whether value range is consecutive upon iteration over range
    map<size_t,size_t>::const_iterator it1(imap.begin()), it2(imap.begin()); it2++;
    while ( it2 != imap.end() ) {
         if ( (*it1).second + 1U != (*it2).second ) {
              cerr <<"\nConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive: detected problem (n vs. n+1): ";
              cerr << (*it1).second <<" vs "<< (*it2).second <<"\n";
              return false;
           }
         it2++;
         it1++;
      }
*/



/*
	bool WithoutAuxArray(int [] arrA)
  {
		//this method with work if numbers are non negative
		int max = findMax(arrA);
		int min = findMin(arrA);
		if(arrA.length!=max-min+1) return false;
		for(int i = 0;i<arrA.length;i++){
			arrA[i] = arrA[i]-min+1;
		}
		for(int i = 0;i<arrA.length;i++){
			int x  = Math.abs(arrA[i]);
			if(arrA[x-1]>0){
				arrA[x-1] = arrA[x-1]*-1;
			}
			else{
				return false;
			}
		}
		return true;
	}
	
  bool WithAuxArray(int [] arrA)
   {
		// this method with work even if numbers are negative
		int []  aux = new int [arrA.length];
		int max = findMax(arrA);
		int min = findMin(arrA);
		if(arrA.length!=max-min+1) return false;
		for(int i = 0;i<arrA.length;i++){
			arrA[i] = arrA[i]-min;
			aux[i] = 0;
		}
		for(int i = 0;i<arrA.length;i++){
			if(aux[arrA[i]]==0){
				aux[arrA[i]]=1;
			}
			else{
				return false;
			}
		}
		//If we have reached till here means , we satisfied all the requirements
		return true;
	}
*/


/** checks whether the values of a map are in consecutive order

    1. checks the individual element ranges for duplicates.
    
    2. combines these ranges into a single one and tests that it is consecutive.
*/
template<>
bool ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive( const multimap<string,vector<size_t> >& immap,
                                                                       bool check_whether_max_value_is_size_minus1 )
 {
    // 1. uniqueness check on all ranges
    for ( typename multimap<string,vector<size_t> >::const_iterator
          it=immap.begin(); it!=immap.end(); ++it ) {
         // for each integer vector the uniqueness of entries is checked
         unordered_set<size_t> iset;
         for ( typename vector<size_t>::const_iterator i=(*it).second.begin(); i!=(*it).second.end(); ++i ) {
              pair<typename unordered_set<size_t>::iterator,bool> set_iterator = iset.insert( (*i) );
              if ( set_iterator.second == false ) {
                   cerr <<"\nConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive(multimap): ";
                   cerr <<" integer range of '"<< (*it).first <<"' contains duplicates.\n";
                   return false;
                }
           }
      }
   
    // 2. combined range check on entire integer set
    set<size_t> iset;
    for ( typename multimap<string,vector<size_t> >::const_iterator
          it=immap.begin(); it!=immap.end(); ++it ) {
         // for each integer vector the uniqueness of entries is checked
        for ( typename vector<size_t>::const_iterator i=(*it).second.begin(); i!=(*it).second.end(); ++i ) {
              pair<typename set<size_t>::iterator,bool> set_iterator = iset.insert( (*i) );
              if ( set_iterator.second == false ) {
                   cerr <<"\nConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive(multimap): ";
                   cerr <<" integer range of '"<< (*it).first <<"' overlaps with one of the previous ones.\n";
                }
           }

      }
    // 3. checking the combined iset
    if ( !IsRangeOfUnsignedIntConsecutive( iset, check_whether_max_value_is_size_minus1 ) ) return false;
        
    return true;
 
 } // end IsValueRangeOfUnsignedIntConsecutive

//template bool ConsecutiveSequenceChecker::IsValueRangeOfUnsignedIntConsecutive( const multimap<string,vector<size_t> >&, bool );



bool ConsecutiveSequenceChecker::Test_ConsecutiveSequenceChecker()
 {
    int failure_count(0);
   
		vector<size_t> arrA = {21,22,24,26,23,25}; // unique, but non consecutive
		vector<size_t> arrB = {21,22,23,24,25,26}; // consecutive and unique
		vector<size_t> arrC = {21,22,23,24,25,25}; // non consecutive and non unique
   
    cerr <<"\nConsecutiveSequenceChecker::Test_ConsecutiveSequenceChecker: 8 tests of the member functions:\n";
   
    if ( !ConsecutiveSequenceChecker().IsRangeOfUnsignedIntConsecutive( arrA, false ) )
      cerr <<"\n1. vector A is non-consecutive as it should.\n";
    else failure_count++;
   
    if ( ConsecutiveSequenceChecker().IsRangeOfUnsignedIntConsecutive( arrB, false ) )
      cerr <<"\n2. vector B is consecutive as it should.\n";
    else failure_count++;
   
    if ( !ConsecutiveSequenceChecker().IsRangeOfUnsignedIntConsecutive( arrC, false ) )
      cerr <<"\n3. vector C is non-consecutive as it should.\n";
    else failure_count++;

    map<size_t,size_t>  mapA = {{0,1},{1,2},{2,3},{3,4},{4,5}};

    if ( ConsecutiveSequenceChecker().IsKeyRangeOfUnsignedIntConsecutive( mapA, false ) )
      cerr <<"\n4. map A keys are consecutive as they should.\n";
    else failure_count++;
   
    if ( ConsecutiveSequenceChecker().IsValueRangeOfUnsignedIntConsecutive( mapA, false ) )
      cerr <<"\n5. map A value range is consecutive as it should.\n";
    else failure_count++;

    if ( !ConsecutiveSequenceChecker().IsValueRangeOfUnsignedIntConsecutive( mapA, true ) )
      cerr <<"\n6. map A value range is larger than size range of values in the map.\n";
    else failure_count++;
   
    map<size_t,size_t>  mapB = {{0,1},{1,2},{3,3},{4,4},{5,5}};

    if ( !ConsecutiveSequenceChecker().IsKeyRangeOfUnsignedIntConsecutive( mapB, false ) )
      cerr <<"\n7. map B Key range is detected as non consecutive as it should.\n";
    else failure_count++;
   
    map<size_t,size_t>  mapC = {{0,4},{1,2},{3,3},{4,0},{5,1}};

    if ( !ConsecutiveSequenceChecker().IsValueRangeOfUnsignedIntConsecutive( mapC, false ) )
      cerr <<"\n8. map C value range is detected as non consecutive as it should.\n";
    else failure_count++;

    return (failure_count == 0 );
   
 } // end Test_ConsecutiveSequenceChecker
 
