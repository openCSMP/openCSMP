#include "IntegrationPointToNodePropertyVisitor_Test.h"
#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "IntegrationPointToNodePropertyVisitor.h"

using namespace std;

namespace csmp{

  
IntegrationPointToNodePropertyVisitor_Test::IntegrationPointToNodePropertyVisitor_Test( Model<3U>& model )
  : model_( model )
{

}

void IntegrationPointToNodePropertyVisitor_Test::run()
{
  enum{DIM=3U};
  // Assuring nodes numbered from 0..n-1
  Region<DIM>& model_domain = model_.Region( "Model" );
  model_domain.RenumberNodes();

  // Creating test properties, initializing values
  PropertyHandle<DIM> cpointProperty( model_, "cpoint property", SCALAR, ELEMENT_INTEGRATION_POINT );
  PropertyHandle<DIM> nodeProperty( model_, "node property", SCALAR, NODE );
  cpointProperty = 1.;
  nodeProperty = 0.;

  // Applying visitor to be tested
  IntegrationPointToNodePropertyVisitor<ScalarVariable,3U> test( model_.Database(), "cpoint property",
                                                    "node property", model_.Region("Model").Nodes());

  test.ApplyWeightingToExtrapolatedValues();
  model_domain.ExtrapolateIntegrationPointToNodeProperty("cpoint property","node property");
  model_domain.Accept(test);
  //test.ApplyWeightingToExtrapolatedValues();

  double nodevalue(0.);
  for (vector<Node<DIM>*>::const_iterator it=model_domain.NodesBegin();it!=model_domain.NodesEnd();it++)
  {
    //std::cout<<"Value: "<<(*it)->Read(model_.Database().StorageKey("node property"))<<std::endl;
    nodevalue =(*it)->Read(model_.Database().StorageKey("node property"));
    _equal(nodevalue,1.,10e-14);
  }

} // run

}
