#include "Vset_TestCase.h"
#include "Boundary.h"
#include "Region.h"
#include "CSMP_highLevelUtilities.h"
#include "ScalarVariable.h"
#include "ArrayVariable.h"
#include "NodeCenteredFiniteVolumeTransport.h"

// File I/O and Initialization
#include "ANSYS_Model3D.h"

using namespace std;

namespace csmp {
Vset_TestCase::Vset_TestCase(const char* prefix)
{
    this->setName("Vset_TestCase");
    prefix_=prefix;
}


Vset_TestCase::~Vset_TestCase()
{
}


/**
    TODO: @todo the VSet test should not be dependent on the Model functionality ! - refactor
*/
void Vset_TestCase::run()
{
    //------------------------------------
    // Parameters section
    string input_file_name(prefix_);
    enum{DIM=3U};
    //end Parameters section
    //------------------------------------
    cout <<"\nStart  of - "<<this->getName()<<endl<<endl;

    //------------------------------------
    // Model Output Regions only test
    //------------------------------------
    
    ScalarVariable diff( ANY, 1. );
    VectorVariable<3> vv( DIRICH, 2. );
    VectorVariable<3> vvPlain;
    TensorVariable<3> tv( ANY, 3. );
    TensorVariable<3> tvPlain;
    
    // Testing model without boundaries, variable&topology tests
    cout <<"Building ModelOutput..."<<endl;
    ANSYS_Model3D modelOutput1( input_file_name.data(),(this->getName()+".txt").c_str(),true,true,true,false);
    ArrayVariable na( "nodal array", modelOutput1.Database(), 2., ROBIN );
    const size_t elementCount1( modelOutput1.Region("Model").Elements() );
    const size_t nodeCount1( modelOutput1.Region("Model").Nodes() );
    const size_t regionCount1( modelOutput1.Regions() );
    modelOutput1.Region("Model").InputPropertyValue( "diffusivity", diff );
    modelOutput1.Region("Model").InputPropertyValue( "nodal array", na );
    modelOutput1.OutputToBinaryFile("model1");
    cout<<"The Output of Model without boundaries done..."<<endl;
    cout <<"Building ModelInput..."<<endl;
    Model<3U> modelInput1("model1");
    Index dKey1( modelOutput1.Database().StorageKey("diffusivity") );
    Index naKey1( modelOutput1.Database().StorageKey("nodal array") );
    _test( elementCount1 == modelInput1.Region("Model").Elements() );
    _test( nodeCount1 == modelInput1.Region("Model").Nodes() );
    _test( regionCount1 == modelInput1.Regions() );
    ScalarVariable scalVal;
    ArrayVariable aVal( "nodal array", modelInput1.Database() );
    for( vector<Element<3>*>::const_iterator it( modelInput1.Region("Model").ElementsBegin() ); it != modelInput1.Region("Model").ElementsEnd(); ++it )
      {
        (*it)->Read( dKey1, scalVal );
        _test( scalVal == diff );
      }
    for( vector<Node<3>*>::const_iterator it( modelInput1.Region("Model").NodesBegin() ); it != modelInput1.Region("Model").NodesEnd(); ++it )
      {
      (*it)->Read( naKey1, aVal );
      _test( aVal == na );
      }
    
    // Testing model without boundaries, variable&topology tests
    cout <<"Building ModelOutput..."<<endl;
    ANSYS_Model3D modelOutput2( input_file_name.data(),(this->getName()+".txt").c_str(),true,true,true,true);
    Index boundaryScalarKey = modelOutput2.Database().StorageKey("boundary scalar");
    Index boundaryArrayKey = modelOutput2.Database().StorageKey("boundary array");
    Index regionVectorKey = modelOutput2.Database().StorageKey("region vector");
    Index modelTensorKey = modelOutput2.Database().StorageKey("model tensor");
    ArrayVariable ba( "boundary array",  modelOutput2.Database() );
    ArrayVariable baPlain( "boundary array",  modelOutput2.Database() );
    ba = 99.;
    const size_t elementCount2( modelOutput2.Region("Model").Elements() );
    const size_t nodeCount2( modelOutput2.Region("Model").Nodes() );
    const size_t regionCount2( modelOutput2.Regions() );
    const size_t boundaryCount2( modelOutput2.Boundaries() );
    modelOutput2.Region("Model").InputPropertyValue( "diffusivity", diff );
    modelOutput2.Region("Model").Store( regionVectorKey, vv );
    modelOutput2.Region("Model").InputPropertyValue( "nodal array", na );
    modelOutput2.Store( modelTensorKey, tv );
    modelOutput2.Boundary("BOUNDARY1").InputPropertyValue("boundary scalar", makeScalar( PLAIN, 1. ) );
    modelOutput2.Boundary("BOUNDARY2").InputPropertyValue("boundary array", ba );
    modelOutput2.OutputToBinaryFile("model2");
    cout<<"The Output of Model with boundaries done..."<<endl;
    cout <<"Building ModelInput with boundaries..."<<endl;
    Model<3U> modelInput2("model2");
    _test( modelInput2.Boundary("BOUNDARY1").Read(boundaryScalarKey) == 1. );
    modelInput2.Region("Model").Read( regionVectorKey, vvPlain );
    _test( vvPlain == vv );
    modelInput2.Read( modelTensorKey, tvPlain );
    modelInput2.Boundary("BOUNDARY2").Read( boundaryArrayKey, baPlain );
    _test( tvPlain == tv );
    _test( baPlain == ba );
    _test( boundaryCount2 == modelInput2.Boundaries() );
    _test( regionCount2 == modelInput2.Regions() );    
    _test( nodeCount2 == modelInput2.Region("Model").Nodes() );  
    

    // testing model with finite volume variables
    ANSYS_Model3D modelOutput3( input_file_name.data(),(this->getName()+".txt").c_str(),true,true,true,true);
    NodeCenteredFiniteVolumeTransport<3> fvModule1( "Model", modelOutput3, "diffusivity", "nodal variable", "element vector", "nodal variable", false, false );
    Index faipVectorKey( modelOutput3.Database().StorageKey("faip vector") );
    Index seipTensorKey( modelOutput3.Database().StorageKey("seip tensor") );
    
    modelOutput3.InputPropertyValue( "faip vector", vv );
    modelOutput3.InputPropertyValue( "seip tensor", tv );
    Element<3>* ePtr = *modelOutput3.Region("Model").ElementsBegin();
    for( size_t f(0); f < ePtr->Facets(); ++f )
      for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
        }
    modelOutput3.OutputToBinaryFile("Vset_TestCase_modelOutput3");
    
    Model<3> modelInput3("Vset_TestCase_modelOutput3");
    
    ePtr = *modelInput3.Region("Model").ElementsBegin();
    size_t ctr(0);
    for( size_t f(0); f < ePtr->Facets(); ++f )
      for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
      {
        ePtr->Read( f, fip, faipVectorKey, vvPlain );
        _test( vvPlain == vv );
        ++ctr;
      }
    for( size_t s(0); s < ePtr->Sectors(); ++s )
      for( size_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
      {
        ePtr->Read( s, sip, seipTensorKey, tvPlain );
        _test( tvPlain == tv );
        ++ctr;
      }
      _test( ctr == 10 );

      NodeCenteredFiniteVolumeTransport<3> fvModule2( "Model", modelInput3, "diffusivity", "nodal variable", "element vector", "nodal variable", false, false );
      for( size_t f(0); f < ePtr->Facets(); ++f )
        for( size_t fip(0); fip < ePtr->IntegrationPointsPerFacet(); ++fip )
        {
          ePtr->Read( f, fip, faipVectorKey, vvPlain );
          _test( vvPlain == vv );
          ++ctr;
        }
    _test( ctr == 16 );
    for( size_t s(0); s < ePtr->Sectors(); ++s )
      for( size_t sip(0); sip < ePtr->IntegrationPointsPerSector(); ++sip )
      {
        ePtr->Read( s, sip, seipTensorKey, tvPlain );
        _test( tvPlain == tv );
        ++ctr;
      }
    
    cout <<"\n\n"<<this->getName()<<" FINISHED!!!"<<endl;
    
} // end Vset_TestCase


} // namespace csmp
