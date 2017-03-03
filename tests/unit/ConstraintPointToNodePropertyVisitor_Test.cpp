#include "ConstraintPointToNodePropertyVisitor_Test.h"
#include "Model.h"
#include "Region.h"
#include "PropertyHandle.h"
#include "IntegrationPointToNodePropertyVisitor.h"

using namespace std;

namespace csmp{

  
ConstraintPointToNodePropertyVisitor_Test::ConstraintPointToNodePropertyVisitor_Test( Model<3U>& model )
  : model_( model )
{

}

void ConstraintPointToNodePropertyVisitor_Test::run()
{
  enum{DIM=3U};
  // Assuring nodes numbered from 0..n-1
  //Region<dim> region = model_.Region( "Model" );
  model_.Region("Model").RenumberNodes();

  // Creating test properties, initializing values
  PropertyHandle<DIM> cpointProperty( model_, "cpoint property", SCALAR, ELEMENT_INTEGRATION_POINT );
  PropertyHandle<DIM> nodeProperty( model_, "node property", SCALAR, NODE );
  cpointProperty = 1.;
  nodeProperty = 0.;

  // Applying visitor to be tested
  IntegrationPointToNodePropertyVisitor<ScalarVariable,3U> test( model_.Database(), "cpoint property",
                                                    "node property", model_.Region("Model").Nodes());

  test.ApplyWeightingToExtrapolatedValues();
  //model_.Region("Model").ExtrapolateConstraintPointToNodeProperty("cpoint property","node property");
  model_.Region("Model").Accept(test);
  //test.ApplyWeightingToExtrapolatedValues();

  double64 nodevalue(0.);
  for (vector<Node<DIM>*>::iterator
       it=  model_.Region("Model").NodesBegin();it!=model_.Region("Model").NodesEnd();it++)
  {
    //std::cout<<"Value: "<<(*it)->Read(model_.Database().StorageKey("node property"))<<std::endl;
    nodevalue =(*it)->Read(model_.Database().StorageKey("node property"));
    _equal(nodevalue,1.,10e-14);
  }
/*
  for (vector<Element<DIM>*>::iterator
       it=  model_.Region("Model").ElementsBegin();it!=model_.Region("Model").ElementsEnd();it++)
  {

    //std::cout<<"Value: "<<(*it)->Read(model_.Database().StorageKey("node property"))<<std::endl;
    nodevalue =(*it)->Read(model_.Database().StorageKey("node property"));
    _equal(nodevalue,1.,10e-14);
  }
*/
  //_test( false );

} // run

}