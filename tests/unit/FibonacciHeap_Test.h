// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#ifndef FIBONACCI_HEAP_TEST_H
#define FIBONACCI_HEAP_TEST_H

#include "Test.h"
#include "FibonacciHeap.h"

namespace csmp {

  class FibonacciHeap_Test : public Test {
    public:
      virtual void run();

    private:
        template<typename K, typename V>
        void doBasicTest( std::deque<std::pair<K,V>>& values );
  };

}

#endif // FIBONACCI_HEAP_TEST_H
