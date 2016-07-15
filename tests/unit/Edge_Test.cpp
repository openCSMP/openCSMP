#include "Edge_Test.h"
#include "ANSYS_Model3D.h"
#include "VTU_Interface.h"
#include "EdgeFactory.h"

namespace csmp{

  Edge_Test::Edge_Test()
    {
    
    }
  
  Edge_Test::~Edge_Test()
    {
    
    }

  void Edge_Test::run()
    {
        ANSYS_Model3D model( "FracBox", "CSMP-variables.txt", true);
        model.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 1. ) );
        Boundary<3>& boundary1(model.Boundary("BOUNDARY1"));
        Boundary<3>& boundary3(model.Boundary("BOUNDARY3"));
        Edge<3> edge( model.Database() );
        EdgeFactory<3>::SetupFrom2Boundaries edgeSetup (boundary1, boundary3);
        EdgeFactory<3>::Create( edge, edgeSetup );
        edge.InputPropertyValue( "nodal variable", makeScalar( PLAIN, 2. ) );
        _test( edge.size() == 5 );
        Index idx( model.Database().StorageKey("nodal variable") );
        for ( std::set<Node<3>*>::const_iterator it = edge.begin(); it != edge.end(); ++it)
            _test( (*it)->Read(idx) == 2. );
        VTU_Interface<3> vtu (model);
        vtu.OutputDataToVTU("EdgeTest", "nodal variable", boundary1, static_cast<int>(0) );
        vtu.OutputDataToVTU("EdgeTest", "nodal variable", boundary3, static_cast<int>(0) );
    }


  
}// csmp
