#ifndef CSMP_UNIONFIND_H
#define CSMP_UNIONFIND_H

#include "CSMP_definitions.h"

namespace csmp {


/**
  Union-find algorithm. This data structure takes connections between two
  elements, and finds all of the connected components and their sizes. A
  "component" is represented by a single representative item.

  See:

    Tarjan, Robert Endre (1975). "Efficiency of a Good But Not Linear
      Set Union Algorithm". JACM. 22 (2): 215–225.
 */
template<typename Item>
class UnionFind {
public:

    /// Say that two items are in the same subset.
    void SameComponent(Item x, Item y);

    /// Extract the components and their sizes.
    template<typename Container>
    void Components(Container& container) const
    {
        for (auto component : components_) {
            auto& rec = records_[component];
            container.push_back(std::make_pair(rec.size_, rec.item_));
        }
    }

    /// Resolve an item to its representative
    Item resolve(Item item) {
        return records_[find_root(ensure(item))].item_;
    }

private:
  struct Record {
    Item item_;
    size_t parent_;
    size_t size_;

    Record(Item item, size_t n)
      : item_(item), parent_(n), size_(1)
    {
    }
  };

  std::map<Item,size_t> item_map_;
  std::deque<Record> records_;
  std::set<size_t> components_;

  // Find the representative node of an item's component by traversing
  // the tree to the root, compressing the path by halves as we traverse.
  size_t find_root(size_t i) {
    while (records_[i].parent_ != i) {
        size_t& parent = records_[i].parent_;
        parent = records_[parent].parent_;
        i = parent;
    }
    return i;
  }

  // Ensure that an item is in the set.
  size_t ensure(Item item) {
    size_t id = records_.size();
    auto result = item_map_.insert(std::make_pair(item, records_.size()));
    if (result.second) {
        records_.emplace_back(item, id);
        components_.insert(id);
        return id;
    }
    else {
        return result.first->second;
    }
  }

};


template<typename Item>
void
UnionFind<Item>::SameComponent(Item x, Item y)
{
   size_t xr = find_root(ensure(x));
   size_t yr = find_root(ensure(y));

   if (xr == yr) {
     // The two iterms are already in the same component.
     return;
   }

   // Add the smaller component to the larger component.
   if (records_[xr].size_ < records_[yr].size_) {
       records_[xr].parent_ = yr;
       records_[yr].size_ += records_[xr].size_;
       components_.erase(xr);
   }
   else {
       records_[yr].parent_ = xr;
       records_[xr].size_ += records_[yr].size_;
       components_.erase(yr);
   }
}

} // csmp

#endif


