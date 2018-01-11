#include "FibonacciHeap_Test.h"

using namespace std;
using namespace ajb;

namespace csmp {



template<typename K, typename V>
void FibonacciHeap_Test::doBasicTest( deque<pair<K,V>>& values )
{
    FibonacciHeap<K,V> h;
    size_t s = 0;

    _test(h.empty());
    _test(h.size() == s);
    for (auto& v : values) {
        h.insert(v.first, v.second);
        _test(!h.empty());
        ++s;
        _test(h.size() == s);
    }

    sort(values.begin(), values.end());
    while (!values.empty())
    {
        _test(!h.empty());
        auto min = h.minimum();
        _test(min->key_ == values.front().first);
        _test(min->value_ == values.front().second);
        h.remove_minimum();
        values.pop_front();
        --s;
        _test(h.size() == s);
    }

    _test(h.empty());
    _test(h.size() == 0);
}


void FibonacciHeap_Test::run()
{
    {
        deque<pair<unsigned,unsigned>> vs;
        vs.push_back(make_pair(4, 0));
        vs.push_back(make_pair(2, 1));
        vs.push_back(make_pair(7, 2));
        vs.push_back(make_pair(5, 3));
        vs.push_back(make_pair(1, 4));
        vs.push_back(make_pair(8, 5));
        doBasicTest(vs);
    }
}
  
}
