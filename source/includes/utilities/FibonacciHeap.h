// Copyright (c) 2015-2017, Andrew J. Bromage
// All rights reserved.
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation files
// (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
//   The above copyright notice and this permission notice shall be
//   included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
// ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#ifndef FIBONACCI_HEAP_H
#define FIBONACCI_HEAP_H

#include <vector>
#include <limits>
#include <deque>
#include <iostream>
#include <cstdint>

namespace ajb {

namespace detail {

    template<typename K, typename V>
    struct FibonacciHeap_Node
    {
        FibonacciHeap_Node(const FibonacciHeap_Node&) = delete;
        FibonacciHeap_Node(FibonacciHeap_Node&&) = delete;
        FibonacciHeap_Node& operator=(const FibonacciHeap_Node&) = delete;

        K key_;
        V value_;
        uint64_t degree_;
        bool mark_;

        FibonacciHeap_Node* prev_;
        FibonacciHeap_Node* next_;
        FibonacciHeap_Node* child_;
        FibonacciHeap_Node* parent_;

        FibonacciHeap_Node(K key, V value)
            : key_(std::move(key)),
              value_(std::move(value)),
              degree_(0),
              mark_(false),
              prev_(this),
              next_(this),
              child_(0),
              parent_(0)
        {
        }

        bool singleton() const noexcept {
            return next_ == this;
        }

        void link(FibonacciHeap_Node* n)  noexcept {
            next_->prev_ = n->prev_;
            n->prev_->next_ = next_;
            next_ = n;
            n->prev_ = this;
        }

        void unlink() noexcept {
            prev_->next_ = next_;
            next_->prev_ = prev_;
            next_ = prev_ = this;
        }

        void add_child(FibonacciHeap_Node* n) noexcept {
            if (!child_) {
                child_ = n;
            }
            else {
                child_->link(n);
            }
            n->parent_ = this;
            n->mark_ = false;
            ++degree_;
        }

        void remove_child(FibonacciHeap_Node* n) noexcept
        {
            if (n->singleton()) {
                child_ = 0;
            }
            else {
                if (child_ == n) {
                    child_ = n->next_;
                }
                n->unlink();
            }
            n->parent_ = 0;
            n->mark_ = false;
            --degree_;
        }
        
        V getV() {return value_;}
        K getK() {return key_;}
    };
}


/// A Fibonacci min-heap
/**
 * See:
 *   Fredman & Tarjan, "Fibonacci heaps and their uses in improved network
 *   optimization algorithms", JACM 34 (3): 596–615, Jul 1987.
 *
 *   doi:10.1145/28869.28874
 */
template<typename K, typename V>
class FibonacciHeap
{
    FibonacciHeap(const FibonacciHeap&) = delete;
    FibonacciHeap(FibonacciHeap&&) = delete;
    FibonacciHeap& operator=(const FibonacciHeap&) = delete;

    typedef K key_type; 
    typedef V value_type; 
    typedef detail::FibonacciHeap_Node<K,V>* node_ptr;
    typedef const detail::FibonacciHeap_Node<K,V>* const_node_ptr;

    node_ptr root_;
    size_t count_;
    size_t max_degree_;
    node_ptr consolidation_[std::numeric_limits<uint32_t>::digits+1];


    void promote_children_of_root()
    {
        if (!root_->child_) {
            return;
        }

        node_ptr child = root_->child_;
        do {
            child->parent_ = 0;
            if (child->degree_ > max_degree_) {
                max_degree_ = child->degree_;
            }
            child = child->next_;
        } while (child != root_->child_);
        root_->child_ = 0;
        root_->link(child);
    }


    void cascading_cut(node_ptr node)
    {
        node_ptr parent = node->parent_;
        for (;;) {
            parent->remove_child(node);
            insert_node(node);
            if (!parent->parent_) {
                return;
            }
            else if (!parent->mark_) {
                parent->mark_ = true;
                return;
            }
            else {
                node = parent;
                parent = parent->parent_;
            }
        }
    }

 
    void remove_current_root()
    {
        if (!root_) {
            return;
        }

        --count_;
        promote_children_of_root();

        if (root_->next_ == root_) {
            node_ptr root = root_;
            root_ = 0;
            delete root;
            return;
        }

        const size_t consolidation_size = sizeof(consolidation_) / sizeof(consolidation_[0]);
        for (auto i = 0; i < consolidation_size; ++i) {
            consolidation_[i] = 0;
        }
        node_ptr curr = root_->next_;
        max_degree_ = 0; 
        do {
            uint64_t d = curr->degree_;
            node_ptr current = curr;
            curr = curr->next_;
            while (consolidation_[d]) {
                node_ptr other = consolidation_[d];
                if (current->key_ > other->key_) {
                    std::swap(other, current);
                }
                other->unlink();
                current->add_child(other);
                consolidation_[d] = 0;
                ++d;
            }
            consolidation_[d] = current;
        } while (curr != root_);
         
        {
            node_ptr root = root_;
            root_ = 0;
            delete root;
        }

        size_t newMaxDegree = 0;

        for (size_t i = 0; i < consolidation_size; ++i) {
            node_ptr newRoot = consolidation_[i];
            if (newRoot) {
                newRoot->next_ = newRoot->prev_ = newRoot;
                insert_node(newRoot);
                if (i > newMaxDegree) {
                    newMaxDegree = i;
                }
            }
        }

        max_degree_ = newMaxDegree;
    }
    
    node_ptr insert_node(node_ptr new_node)
    {
        if (!root_) {
            root_ = new_node;
        }
        else {
            root_->link(new_node);
            if (new_node->key_ < root_->key_)
            {
                root_ = new_node;
            }
        }
        return new_node;
    }    
    

public:
    typedef node_ptr finger;
    typedef const_node_ptr const_finger;

    /// Constructor
    FibonacciHeap()
        : root_(0), count_(0), max_degree_(0)
    {    
    }    

    /// Remove all elements from the heap
    void clear()
    {
        if (root_) {
            while (root_->next_ != root_ && !root_->singleton()) {
                promote_children_of_root();
                node_ptr root = root_;
                root_ = root->next_;
                root->unlink();
                delete root;
            }
           node_ptr root = root_;
           root_ = 0;
           delete root;
        }
        count_ = 0;
        max_degree_ = 0;
    }

    /// Destructor
    ~FibonacciHeap()
    {
        clear();
    }

    /// Number of elements in the heap
    uint64_t size() const
    {
        return count_;
    }

    /// Is the heap empty
    bool empty() const
    {
        return count_ == 0;
    }


    /// Insert an element in the heap (move semantics)
    finger insert(K key, V value)
    {
        ++count_;
        return insert_node(new detail::FibonacciHeap_Node<K,V>(key, value));
    }
    
    
    finger insert(node_ptr new_node)
    {
        ++count_;
        return insert_node(new_node);
    }    


    /// Get the minimum element in the heap
    finger minimum()
    {
        return root_;
    }


    /// Get the minimum element in the heap
    const_finger minimum() const
    {
        return root_;
    }
 
 
    /// Remove the minimum element in the heap
    void remove_minimum()
    {
        remove_current_root();
    }


    /// Decrease the key of an element in the heap
    void decrease_key(finger node, K new_key)
    {
        if (new_key > node->key_) {
            throw "Trying to decrease key to a greater key.";
        }
        node->key_ = std::move(new_key);
        node_ptr parent = node->parent_;
        if (!parent) {
            if (node->key_ < root_->key_) {
                root_ = node;
            }
            return;
        }
        else if (parent->key_ <= new_key) {
            return;
        }

        cascading_cut(node);
    }


    /// Remove an element from the heap
    void remove(node_ptr node)
    {
        if (node->parent_) {
            cascading_cut(node);
        }
        root_ = node;
        remove_current_root();
    }

};

} // csmp


#endif // FIBONACCI_HEAP_H

