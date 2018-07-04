#ifndef CSMP_SMALLSET_H
#define CSMP_SMALLSET_H

#include "CSMP_definitions.h"

namespace csmp {

 /**
  A small set of trivial types, represented as a sorted vector.
 */

  template<typename T, unsigned N = 32>
  class SmallSet {
  public:
    static_assert(std::is_trivial<T>::value, "Type must be trivial");
    static_assert(N <= 32, "The optimised size of this set is probably too large");

    typedef T value_type;
    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;

    SmallSet()
    {
      storage_.reserve(N);
    }

    std::pair<iterator,bool> insert( const value_type& value )
    {
      auto it = std::lower_bound(storage_.begin(), storage_.end(), value);
      if (it != storage_.end() && *it == value) {
        return std::make_pair(it,false);
      }
      else {
        auto newit = storage_.insert(it, value);
        return std::make_pair(newit, true);
      }

    }

    std::pair<iterator,bool> insert( value_type&& value )
    {
      auto it = std::lower_bound(storage_.begin(), storage_.end(), value);
      if (it != storage_.end() && *it == value) {
        return std::make_pair(it,false);
      }
      else {
        auto newit = storage_.insert(it, value);
        return std::make_pair(newit, true);
      }
    }

    bool empty() const
    {
      return storage_.empty();
    }

    size_t size() const
    {
      return storage_.size();
    }

    iterator begin()
    {
      return storage_.begin();
    }

    const_iterator begin() const
    {
      return storage_.begin();
    }

    iterator end()
    {
      return storage_.end();
    }

    const_iterator end() const
    {
      return storage_.end();
    }

    iterator find(const value_type& value)
    {
      auto it = std::lower_bound(storage_.begin(), storage_.end(), value);
      if (it == storage_.end() || *it != value) {
        return storage_.end();
      }
      else {
        return it;
      }
    }

    const_iterator find(const value_type& value) const
    {
      auto it = std::lower_bound(storage_.begin(), storage_.end(), value);
      if (it == storage_.end() || *it != value) {
        return storage_.end();
      }
      else {
        return it;
      }
    }

    void erase(const value_type& value)
    {
      auto it = std::lower_bound(storage_.begin(), storage_.end(), value);
      if (it != storage_.end() && *it == value) {
        storage_.erase(it);
      }
    }

    size_t count(const value_type& value) const
    {
      auto it = std::lower_bound(storage_.begin(), storage_.end(), value);
      return (it == storage_.end() || *it != value) ? 0 : 1;
    }


  private:
    std::vector<T> storage_;
  };


} // csmp

#endif




