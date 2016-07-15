#include "Mesh_Test.h"

#include "PropertyDatabase.h"
#include "FiniteElementManager.h"

namespace csmp
{
    Mesh_Test::Mesh_Test()
    {

    }

    Mesh_Test::~Mesh_Test()
    {
        if( database_!=NULL )
            delete database_;
        if(fem_manager_!=NULL)
            delete fem_manager_;
        if(vset_!=NULL)
            delete vset_;
        if(mesh_!=NULL)
            delete mesh_;
    }

    void Mesh_Test::run()
    {
        this->database_=new PropertyDatabase<3>("CSMP_Variables.txt");
        this->fem_manager_=new FiniteElementManager(3,2,true);

        this->vset_=new VSet<3U>();

        this->mesh_=new MeshManager<3U>(*database_,*fem_manager_,*vset_);

    }
}
