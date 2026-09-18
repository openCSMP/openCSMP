// Copyright © 2016 Stephan Matthai. All rights reserved.
// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef CONSECUTIVE_INTEGER_CHECKER_H
#define CONSECUTIVE_INTEGER_CHECKER_H

#include <iostream>
#include <map>
#include <set>
#include <vector>


/**

@class ConsecutiveSequenceChecker - checks integer sequences for whether they are consecutive.

@author SKM 
@date 10/12/2016

1. Finds the max­i­mum and min­i­mum ele­ments in sequences
2. Checks if the sequence length = max-min+1
3. Checks if the sequence contains duplicates
4. checks whether the elements of the sequence are consecutuve

@note current implementation is restricted to size_t  unsigned integers.

@note use Test_ConsecutiveSequenceChecker() to see verification of methods.

*/
class ConsecutiveSequenceChecker {
  public:
    /// checks elements of an integer set
    template<typename T>
    static bool IsRangeOfUnsignedIntConsecutive( const std::set<T>&, bool check_whether_max_value_is_size_minus1 );

    /// checks elements of an integer vector
    template<typename T>
    static bool IsRangeOfUnsignedIntConsecutive( const std::vector<T>&, bool check_whether_max_value_is_size_minus1 );

    /// checks whether the keys of a map are in consecutive order
    template<typename K, typename V>
    static bool IsKeyRangeOfUnsignedIntConsecutive( const std::map<K,V>&, bool check_whether_max_value_is_size_minus1 );

    /// checks whether the values of a map are in consecutive order
    template<typename K, typename V>
    static bool IsValueRangeOfUnsignedIntConsecutive( const std::map<K,V>&, bool check_whether_max_value_is_size_minus1 );
    
    /// checks whether the values in the map are unique, form a continuous sequence and (optionally) are bounded by max-val == n-values-1
    template<typename K, typename V>
    static bool IsValueRangeUniqueAndBounded( const std::map<K,V>&, bool check_whether_max_value_is_size_minus1 );
  
    /// checks whether the values of a map are in consecutive order
    template<typename K, typename V>
    static bool IsValueRangeOfUnsignedIntConsecutive( const std::multimap<K,std::vector<V> >&, bool check_whether_max_value_is_size_minus1 );

    /// verifies the interfaces of this class
    static bool Test_ConsecutiveSequenceChecker();
};

#endif /* CONSECUTIVE_SEQUENCE_CHECKER_H */
