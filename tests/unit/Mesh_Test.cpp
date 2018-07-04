#include "Mesh_Test.h"

#include "PropertyDatabase.h"
#include "FiniteElementManager.h"

namespace csmp {


Mesh_Test::Mesh_Test()
  : database_(nullptr),
    fem_manager_(nullptr),
    vset_(nullptr),
    mesh_(nullptr)
  {
  }


Mesh_Test::~Mesh_Test()
  {
     delete database_;
     delete fem_manager_;
     delete vset_;
     delete mesh_;
  }


void Mesh_Test::run()
  {
      this->database_=new PropertyDatabase<3>("CSMP_Variables.txt");
      this->fem_manager_=new FiniteElementManager(3,2,true);

      this->vset_=new VSet<3U>();

      this->mesh_=new MeshManager<3U>(*database_,*fem_manager_,*vset_);
  }

} // end csmp
