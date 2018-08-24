#ifndef CSMP_PRIMITIVECONTAINER_H
#define CSMP_PRIMITIVECONTAINER_H

#include "CSMP_definitions.h"

namespace csmp {


template<typename T,class Traits> class PrimitiveContainerIterator;
template<typename T,class Traits> class PrimitiveContainerConstIterator;


/** 
    @brief Traits class for primitives
 
    A memory-efficient container for storing primitives (elements, faces, interfaces, nodes) in a mesh.

    Pointers are never invalidated, even when primitives are deleted.

    Recycles storage when primitives are erased.

    The main trait required is Clear(), which frees up dynamic storage inside
    the primitive. This is called when a primitive is erased.

    @author A.J. Bromage
    @date 2017

    @todo AJB - The dereference operation Index() should be safer.
    @todo AJB - The iterators should conform to OutputIterator as well as ForwardIterator.
    @todo AJB - Bad things will happen if you erase the Root element. Don't do that.
*/
template<typename T>
struct PrimitiveTraits {
    /// Clear a primitive so that it releases any dynamic storage
    void Clear(T* prim) { }
};

template<typename T, class Traits = PrimitiveTraits<T>>
class PrimitiveContainer
{
      friend class PrimitiveContainerIterator<T,Traits>;
      friend class PrimitiveContainerConstIterator<T,Traits>;
  public:
      /// STL boilerplate
      typedef T value_type;
      typedef T& reference;
      typedef const T& const_reference;
      typedef std::ptrdiff_t difference_type;
      typedef std::size_t size_type;
      typedef PrimitiveContainerIterator<T,Traits> iterator;
      typedef PrimitiveContainerConstIterator<T,Traits> const_iterator;
  
      /// Default copy-type operations are not supported, because pointers are important. Use Assign() below.
      PrimitiveContainer(const PrimitiveContainer&) = delete;
      PrimitiveContainer(PrimitiveContainer&&) = delete;
      PrimitiveContainer& operator=(const PrimitiveContainer&) = delete;

      PrimitiveContainer() { }
      ~PrimitiveContainer() { }

      /// Allocate a primitive
      T* Alloc();

      /// Delete a primitive
      void Free(T* prim);

      /// Copy a container.
      // This does not copy free primitives.
      void Assign(const PrimitiveContainer& rhs);

      /// Root primitive
      const T& Root() const
      {
        auto ii = begin();
        assert(ii != end());
        return *ii;
      }

      T& Root()
      {
        auto ii = begin();
        assert(ii != end());
        return *ii;
      }

      /// The number of primitives in the container.
      size_type size() const
      {
          return storage_.size() - free_list_.size();
      }

      /// True if there are no primitives.
      bool empty() const { return storage_.size() == free_list_.size(); }

      /// Iterator support.
      iterator begin();
      iterator end();
      const_iterator begin() const;
      const_iterator end() const;

      /// True if it's safe to call Index()
      bool IndexOperationIsSafe() const { return free_list_.empty(); }

      /// Dereference operator. Only usable during mesh construction or after Assign().
      T& Index(size_t i);

      /// Dereference operator. Only usable during mesh construction or after Assign().
      const T& Index(size_t i) const;

      /// Emplace a primitive.
      // This returns a pointer to the new primitive.
      template<typename... Args>
      T* Emplace(Args&&... args)
      {
        if (free_list_.empty()) {
          storage_.emplace_back(std::forward<Args>(args)...);
          return &storage_.back();
        }

        auto flit = free_list_.begin();
        T* prim = *flit;
        free_list_.erase(flit);
        *prim = T(std::forward<Args>(args)...);
        return prim;
      }

      // TODO: AJB: Delete this
      bool OnFreeList(const T* item) const {
          return free_list_.count(const_cast<T*>(item)) > 0;
      }

  private:
      std::deque<T> storage_;
      std::unordered_set<T*> free_list_;
};


template<typename T, class Traits>
T*
PrimitiveContainer<T,Traits>::Alloc()
{
  if (free_list_.empty()) {
      storage_.emplace_back();
      return &storage_.back();
  }

  T* prim = free_list_.back();
  free_list_.pop_back();
  return prim;
}

template<typename T, class Traits>
void
PrimitiveContainer<T,Traits>::Free(T* prim)
{
  Traits traits;
  traits.Clear(prim);
  free_list_.insert(prim);
}


template<typename T, class Traits>
inline T&
PrimitiveContainer<T,Traits>::Index(size_t i)
{
  assert(i < storage_.size());
  T& prim = storage_[i];
  assert(!free_list_.count(&prim));
  return prim;
}


template<typename T, class Traits>
inline const T&
PrimitiveContainer<T,Traits>::Index(size_t i) const
{
  assert(i < storage_.size());
  const T& prim = storage_[i];
  assert(!free_list_.count(const_cast<T*>(&prim)));
  return prim;
}


/// The PrimitiveContainer's iterator.
template<typename T, class Traits>
class PrimitiveContainerIterator {
public:
  typedef typename std::deque<T>::iterator underlying_it;

    PrimitiveContainerIterator()
        : container_(0), it_()
    {
    }

    PrimitiveContainerIterator(
            const PrimitiveContainer<T,Traits>* container, underlying_it it)
        : container_(container), it_(it)
    {
        AdvanceToNextValidItem();
    }

    bool operator==(const PrimitiveContainerIterator& rhs) const {
        return container_ == rhs.container_ && it_ == rhs.it_;
    }

    bool operator!=(const PrimitiveContainerIterator& rhs) const {
        return !(*this == rhs);
    }

    T& operator*() const {
        return *it_;
    }

    T* operator->() const {
        return &*it_;
    }

    const PrimitiveContainerIterator& operator++() {
        ++it_;
        AdvanceToNextValidItem();
        return *this;
    }

    PrimitiveContainerIterator operator++(int) {
        auto result = *this;
        ++(*this);
        return result;
    }
  
private:
  // To allow conversion from iterator to const_iterator
  friend class PrimitiveContainerConstIterator<T,Traits>;

  void AdvanceToNextValidItem() {
    while (it_ != container_->storage_.end() && container_->free_list_.count(&*it_)) {
        ++it_;
    }
  }
  
  const PrimitiveContainer<T,Traits>* container_;
  underlying_it it_;

};

/// The PrimitiveContainer's const_iterator.
template<typename T, class Traits>
class PrimitiveContainerConstIterator {

public:
  typedef typename std::deque<T>::const_iterator underlying_it;

    PrimitiveContainerConstIterator()
        : container_(0), it_()
    {
    }

    PrimitiveContainerConstIterator(const PrimitiveContainerIterator<T,Traits>& it)
        : container_(it.container_), it_(it.it_)
    {
    }

    PrimitiveContainerConstIterator(
            const PrimitiveContainer<T,Traits>* container, underlying_it it)
        : container_(container), it_(it)
    {
        AdvanceToNextValidItem();
    }

    bool operator==(const PrimitiveContainerConstIterator& rhs) const {
        return container_ == rhs.container_ && it_ == rhs.it_;
    }

    bool operator!=(const PrimitiveContainerConstIterator& rhs) const {
        return !(*this == rhs);
    }

    const T& operator*() const {
        return *it_;
    }

    const T* operator->() const {
        return &*it_;
    }

    const PrimitiveContainerConstIterator& operator++() {
        ++it_;
        AdvanceToNextValidItem();
        return *this;
    }

    PrimitiveContainerConstIterator operator++(int) {
        auto result = *this;
        ++(*this);
        return result;
    }
  
private:
  void AdvanceToNextValidItem() {
    while (it_ != container_->storage_.end() &&
           container_->free_list_.count(const_cast<T*>(&*it_))) {
        ++it_;
    }
  }
  
  const PrimitiveContainer<T,Traits>* container_;
  underlying_it it_;
};


template<typename T,class Traits>
inline
typename PrimitiveContainer<T,Traits>::iterator
PrimitiveContainer<T,Traits>::begin()
{
    return iterator(this, storage_.begin());
}

template<typename T,class Traits>
inline
typename PrimitiveContainer<T,Traits>::iterator
PrimitiveContainer<T,Traits>::end()
{
    return iterator(this, storage_.end());
}

template<typename T,class Traits>
inline
typename PrimitiveContainer<T,Traits>::const_iterator
PrimitiveContainer<T,Traits>::begin() const
{
    return const_iterator(this, storage_.begin());
}

template<typename T,class Traits>
inline
typename PrimitiveContainer<T,Traits>::const_iterator
PrimitiveContainer<T,Traits>::end() const
{
    return const_iterator(this, storage_.end());
}


template<typename T,class Traits>
void
PrimitiveContainer<T,Traits>::Assign(const PrimitiveContainer<T,Traits>& rhs)
{
    for (auto& p : rhs) {
      storage_.push_back(p);
    }
}

} // csmp


/// Iterator traits so PrimitiveContainer will work with STL algorithms.
namespace std {

  template<typename T, class Traits>
  struct iterator_traits< csmp::PrimitiveContainerIterator<T,Traits> >
  {
    typedef ptrdiff_t difference_type;
    typedef T value_type;
    typedef T* pointer;
    typedef T& reference;
    typedef std::forward_iterator_tag iterator_category;
  };

  template<typename T, class Traits>
  struct iterator_traits< csmp::PrimitiveContainerConstIterator<T,Traits> >
  {
    typedef ptrdiff_t difference_type;
    typedef const T value_type;
    typedef const T* pointer;
    typedef const T& reference;
    typedef std::forward_iterator_tag iterator_category;
  };

} // std

#endif


