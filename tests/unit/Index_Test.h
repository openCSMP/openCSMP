#ifndef INDEX_TEST_H
#define INDEX_TEST_H

#include "Test.h"


namespace csmp
{
    class Index;
    class Index_Test: public Test
    {
    public:
        Index_Test();
        ~Index_Test();
        void IndexOperatorsTest(Index&);
        void run();
    private:
    };
}
#endif // INDEX_TEST_H
