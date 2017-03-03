#ifndef MESH_TEST_H
#define MESH_TEST_H

#include "Test.h"
#include "MeshManager.h"

#include "VSet.h"

namespace csmp
{
    template<size_t> class PropertyDatabase;
    class FiniteElementManager;

    class Mesh_Test: public Test
    {
    public:
        Mesh_Test();
        ~Mesh_Test();
        void run();
    private:
        MeshManager<3U>* mesh_;
        PropertyDatabase<3>* database_;
        FiniteElementManager* fem_manager_;
        VSet<3U>* vset_;
    };
}
#endif // MESH_TEST_H
