#include "Placement_Test.h"

using namespace std;

namespace csmp {


  Placement_Test::~Placement_Test() {
    delete mockModel;
  }

  Placement_Test::Placement_Test() :nodeScalar(1.),
    elementScalar(2.),
    elementIPScalar(3.),
    facetIPScalar(4.),
    sectorIPScalar(5.),
    vx(1.), vy(2.),
    xx(1.), xy(2.), yx(3.), yy(4.),
    tx(1.), ty(2.), tz(3.),
    arx(1.), ary(2.), arz(3.),
    floatTolerance(1e-6)
  {
    
    this->mockModel = new  ANSYS_Model2D("mockModel", "mockModelVariable.txt", true, false);
    //this->mockModel->InstantiateFiniteVolumes();
    const Region<2U>& region = mockModel->Region("Model");

    // initialize the value for 
    // node
    csmp::INDEX<SCALAR, NODE> keyNodeScalar(mockModel->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyNodeVector(mockModel->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyNodeTensor(mockModel->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE>  keyNodeArray(mockModel->Database().StorageKey("nodeArray")); 
    csmp::INDEX<FLAGGEDARRAY, NODE> keyNodeFlaggedArray(mockModel->Database().StorageKey("nodeFlaggedArray"));
    
    // element
    csmp::INDEX<SCALAR, ELEMENT> keyElementScalar(mockModel->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyElementVector(mockModel->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyElementTensor(mockModel->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyElementArray(mockModel->Database().StorageKey("elementArray"));
    csmp::INDEX<FLAGGEDARRAY, ELEMENT>keyElementFlaggedArray(mockModel->Database().StorageKey("elementFlaggedArray"));

    // elementIP
    csmp::INDEX<SCALAR, ELEMENT_INTEGRATION_POINT> keyElementIPScalar(mockModel->Database().StorageKey("elementIPScalar"));
    csmp::INDEX<VECTOR, ELEMENT_INTEGRATION_POINT> keyElementIPVector(mockModel->Database().StorageKey("elementIPVector"));
    csmp::INDEX<TENSOR, ELEMENT_INTEGRATION_POINT> keyElementIPTensor(mockModel->Database().StorageKey("elementIPTensor"));
    csmp::INDEX<ARRAY, ELEMENT_INTEGRATION_POINT>  keyElementIPArray(mockModel->Database().StorageKey("elementIPArray"));
    csmp::INDEX<FLAGGEDARRAY, ELEMENT_INTEGRATION_POINT> keyElementIPFlaggedArray(mockModel->Database().StorageKey("elementIPFlaggedArray"));

    // sectorIP
    csmp::INDEX<SCALAR, SECTOR_INTEGRATION_POINT> keySectorScalar(mockModel->Database().StorageKey("sectorIPScalar"));
    csmp::INDEX<VECTOR, SECTOR_INTEGRATION_POINT> keySectorVector(mockModel->Database().StorageKey("sectorIPVector"));
    csmp::INDEX<TENSOR, SECTOR_INTEGRATION_POINT> keySectorTensor(mockModel->Database().StorageKey("sectorIPTensor"));
    csmp::INDEX<ARRAY, SECTOR_INTEGRATION_POINT>  keySectorArray(mockModel->Database().StorageKey("sectorIPArray"));
    csmp::INDEX<FLAGGEDARRAY, SECTOR_INTEGRATION_POINT> keySectorFlaggedArray(mockModel->Database().StorageKey("sectorIPFlaggedArray"));

    // facetIP
    csmp::INDEX<SCALAR, FACET_INTEGRATION_POINT> keyFacetScalar(mockModel->Database().StorageKey("facetIPScalar"));
    csmp::INDEX<VECTOR, FACET_INTEGRATION_POINT> keyFacetVector(mockModel->Database().StorageKey("facetIPVector"));
    csmp::INDEX<TENSOR, FACET_INTEGRATION_POINT> keyFacetTensor(mockModel->Database().StorageKey("facetIPTensor"));
    csmp::INDEX<ARRAY, FACET_INTEGRATION_POINT>  keyFacetArray(mockModel->Database().StorageKey("facetIPArray"));
    csmp::INDEX<FLAGGEDARRAY, FACET_INTEGRATION_POINT> keyFacetFlaggedArray(mockModel->Database().StorageKey("facetIPFlaggedArray"));
    
    
    ArrayVariable arrayVar(3, 0., ANY);
    arrayVar.Component(0, arx);
    arrayVar.Component(1, ary);
    arrayVar.Component(2, arz);

    vector<Node<2U>*>::const_iterator nodeIterator;
    for (nodeIterator = region.NodesBegin(); nodeIterator != region.NodesEnd(); nodeIterator++) {
      (*nodeIterator)->Store(keyNodeScalar, makeScalar(ANY, nodeScalar));
      (*nodeIterator)->Store(keyNodeVector, makeVector(ANY, ANY, vx, vy));
      (*nodeIterator)->Store(keyNodeTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
      (*nodeIterator)->Store(keyNodeArray, arrayVar); 
    }

    // element

    vector<Element<2U>*>::const_iterator eleIter;
    for (eleIter = region.ElementsBegin(); eleIter != region.ElementsEnd(); eleIter++) {
      (*eleIter)->Store(keyElementScalar, makeScalar(ANY, elementScalar));
      (*eleIter)->Store(keyElementVector, makeVector(ANY, ANY, vx, vy));
      (*eleIter)->Store(keyElementTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
      (*eleIter)->Store(keyElementArray, arrayVar);

      for (size_t i(0); i < (*eleIter)->IntegrationPoints(); ++i) {
        (*eleIter)->Store(i, keyElementIPScalar, makeScalar(ANY, elementIPScalar));
        (*eleIter)->Store(i, keyElementIPVector, makeVector(ANY, ANY, vx, vy));
        (*eleIter)->Store(i, keyElementIPTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
        (*eleIter)->Store(i, keyElementIPArray, arrayVar);
      }

      for (size_t i(0); i < (*eleIter)->Sectors(); ++i) {
        (*eleIter)->Store(i, 0, keySectorScalar, makeScalar(ANY, sectorIPScalar));
        (*eleIter)->Store(i, 0, keySectorVector, makeVector(ANY, ANY, vx, vy));
        (*eleIter)->Store(i, 0, keySectorTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
        (*eleIter)->Store(i, 0, keySectorArray, arrayVar);

      }

      for (size_t i(0); i < (*eleIter)->Facets(); ++i) {
        (*eleIter)->Store(i, 0, keyFacetScalar, makeScalar(ANY, facetIPScalar));
        (*eleIter)->Store(i, 0, keyFacetVector, makeVector(ANY, ANY, vx, vy));
        (*eleIter)->Store(i, 0, keyFacetTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
        (*eleIter)->Store(i, 0, keyFacetArray, arrayVar);
      }
    }
  }


  void Placement_Test::run() {

    NodeToElementTest();
    ElementToElementTest();
    ElementIPToElementTest();
    SectorIPToElementTest();
    FacetIPToElementTest();

    // target elementIP
    NodeToElementIPTest();
    ElementToElementIPTest();
    //SectorIPToElementIPTest();  
    //FacetIPToElementIPTest();   

    // target facetIP
    NodeToFacetIPTest();
    ElementToFacetIPTest();
    //SectorIPToFacetIPTest();   


    // target sectorIP
    NodeToSectorIPTest();
    ElementToSectorIPTest();

    // target node
    // ElementIPToNodeTest(); 
    // SectorIPToNodeTest();  
    // FacetIPToNodeTest();  

  }


  // target element

  void Placement_Test::NodeToElementTest()
  {
    //scalar, vector and tensor
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(mockModel->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(mockModel->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(mockModel->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE>  keyArray(mockModel->Database().StorageKey("nodeArray"));
   
    // FE_NODE_TO_BCTR
    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++ ) {
      
      auto ePlacement =(* eit)->AtBarycenter();

      double64 bctrSca(ePlacement.Obtain(keyScalar));
      ScalarVariable bctrSca_;
      ePlacement.Obtain(keyScalar, bctrSca_);
      VectorVariable<2U> bctrVector;
      ePlacement.Obtain(keyVector, bctrVector);
      TensorVariable<2U> bctrTensor;
      ePlacement.Obtain(keyTensor, bctrTensor);
      ArrayVariable bctrArray(3, 0., ANY);
      ePlacement.Obtain(keyArray, bctrArray);

      ScalarVariable scalVal(ANY, 0.);
      ScalarVariable scalTemp;
      VectorVariable<2U> vecVal;
      VectorVariable<2U> vecTemp;
      TensorVariable<2U> tenVal;
      TensorVariable<2U> tenTemp;
      ArrayVariable arrVal(3, 0., ANY);
      ArrayVariable arrTemp;


      size_t nodeNum((*eit)->Nodes());
      vector<double64> NN;
      (*eit)->N_AtBaryCenter(NN);

      (*eit)->PropertyValueAtBaryCenter(keyScalar, scalVal);
      (*eit)->PropertyValueAtBaryCenter(keyVector, vecVal);
      (*eit)->PropertyValueAtBaryCenter(keyTensor, tenVal);
      (*eit)->PropertyValueAtBaryCenter(keyArray, arrVal);


      //scalar
      _equal(scalVal(), bctrSca, floatTolerance);
      _equal(scalVal(), bctrSca_(), floatTolerance);

      //vector
      _equal(vecVal[0], bctrVector[0], floatTolerance);
      _equal(vecVal[1], bctrVector[1], floatTolerance);

      //tensor
      _equal(tenVal(0, 0), bctrTensor(0, 0), floatTolerance);
      _equal(tenVal(0, 1), bctrTensor(0, 1), floatTolerance);
      _equal(tenVal(1, 0), bctrTensor(1, 0), floatTolerance);
      _equal(tenVal(1, 1), bctrTensor(1, 1), floatTolerance);

      //array
      _equal(arrVal[0], bctrArray[0], floatTolerance);
      _equal(arrVal[1], bctrArray[1], floatTolerance);
      _equal(arrVal[2], bctrArray[2], floatTolerance);
    }

  }


  //
  void Placement_Test::ElementToElementTest()
  {
    //scalar and vector
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(mockModel->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(mockModel->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(mockModel->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyArray(mockModel->Database().StorageKey("elementArray"));

    // FE_ELE_TO_BCTR
    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {

      auto ePlacement = (*eit)->AtBarycenter();

      double64 bctrSca(ePlacement.Obtain(keyScalar));
      ScalarVariable bctrSca_;
      ePlacement.Obtain(keyScalar, bctrSca_);
      VectorVariable<2U> bctrVector;
      ePlacement.Obtain(keyVector, bctrVector);
      TensorVariable<2U> bctrTensor;
      ePlacement.Obtain(keyTensor, bctrTensor);
      ArrayVariable bctrArray(3, 0., ANY);
      ePlacement.Obtain(keyArray, bctrArray);

      ScalarVariable scalVal(ANY, 0.);
      VectorVariable<2U> vecVal;
      TensorVariable<2U> tenVal;
      ArrayVariable arrVal(3, 0., ANY);


      (*eit)->PropertyValueAtBaryCenter(keyScalar, scalVal);
      (*eit)->PropertyValueAtBaryCenter(keyVector, vecVal);
      (*eit)->PropertyValueAtBaryCenter(keyTensor, tenVal);
      (*eit)->PropertyValueAtBaryCenter(keyArray, arrVal);

      //scalar
      _equal(scalVal(), bctrSca, floatTolerance);
      _equal(scalVal(), bctrSca_(), floatTolerance);

      //vector
      _equal(vecVal[0], bctrVector[0], floatTolerance);
      _equal(vecVal[1], bctrVector[1], floatTolerance);

      //tensor
      _equal(tenVal(0, 0), bctrTensor(0, 0), floatTolerance);
      _equal(tenVal(0, 1), bctrTensor(0, 1), floatTolerance);
      _equal(tenVal(1, 0), bctrTensor(1, 0), floatTolerance);
      _equal(tenVal(1, 1), bctrTensor(1, 1), floatTolerance);

      //array
      _equal(arrVal[0], bctrArray[0], floatTolerance);
      _equal(arrVal[1], bctrArray[1], floatTolerance);
      _equal(arrVal[2], bctrArray[2], floatTolerance);
    }

  }

  void Placement_Test::ElementIPToElementTest()
  {

    //scalar, vector and tensor
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT_INTEGRATION_POINT> keyScalar(mockModel->Database().StorageKey("elementIPScalar"));
    csmp::INDEX<VECTOR, ELEMENT_INTEGRATION_POINT> keyVector(mockModel->Database().StorageKey("elementIPVector"));
    csmp::INDEX<TENSOR, ELEMENT_INTEGRATION_POINT> keyTensor(mockModel->Database().StorageKey("elementIPTensor"));
    csmp::INDEX<ARRAY, ELEMENT_INTEGRATION_POINT>  keyArray(mockModel->Database().StorageKey("elementIPArray"));

    // FE_EIP_TO_BCTR

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {

      auto ePlacement = (*eit)->AtBarycenter();
      
      double64 bctrSca(ePlacement.Obtain(keyScalar));
      ScalarVariable bctrSca_;
      ePlacement.Obtain(keyScalar, bctrSca_);

      VectorVariable<2U> bctrVector;
      ePlacement.Obtain(keyVector, bctrVector);
      TensorVariable<2U> bctrTensor;
      ePlacement.Obtain(keyTensor, bctrTensor);
      ArrayVariable bctrArray(3, 0., ANY);
      ePlacement.Obtain(keyArray, bctrArray);
      vector<double64> NN;
      (*eit)->N_AtBaryCenter(NN);

      ScalarVariable scalVal;
      (*eit)->PropertyValueAtBaryCenter(keyScalar, scalVal);
      VectorVariable<2U> vecVal(ANY, ANY, 0., 0.);
      (*eit)->PropertyValueAtBaryCenter(keyVector, vecVal);
      TensorVariable<2U> tenVal(ANY, ANY, 0., 0., 0., 0.);
      (*eit)->PropertyValueAtBaryCenter(keyTensor, tenVal);
      ArrayVariable arrVal(3, 0., ANY);
      (*eit)->PropertyValueAtBaryCenter(keyArray, arrVal);

      //scalar
      _equal(scalVal(), bctrSca, floatTolerance);
      _equal(scalVal(), bctrSca_(), floatTolerance); //TODO: ScalarVariable.Component(size_t arbitrary value?)

      //vector
      _equal(vecVal[0], bctrVector[0], floatTolerance);
      _equal(vecVal[1], bctrVector[1], floatTolerance);

      //tensor
      _equal(tenVal(0, 0), bctrTensor(0, 0), floatTolerance);
      _equal(tenVal(0, 1), bctrTensor(0, 1), floatTolerance);
      _equal(tenVal(1, 0), bctrTensor(1, 0), floatTolerance);
      _equal(tenVal(1, 1), bctrTensor(1, 1), floatTolerance);

      //array
      _equal(arrVal[0], bctrArray[0], floatTolerance);
      _equal(arrVal[1], bctrArray[1], floatTolerance);
      _equal(arrVal[2], bctrArray[2], floatTolerance);
    }


  }


  void Placement_Test::SectorIPToElementTest()
  {

    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, SECTOR_INTEGRATION_POINT> keyScalar(mockModel->Database().StorageKey("sectorIPScalar"));
    csmp::INDEX<VECTOR, SECTOR_INTEGRATION_POINT> keyVector(mockModel->Database().StorageKey("sectorIPVector"));
    csmp::INDEX<TENSOR, SECTOR_INTEGRATION_POINT> keyTensor(mockModel->Database().StorageKey("sectorIPTensor"));
    csmp::INDEX<ARRAY, SECTOR_INTEGRATION_POINT>  keyArray(mockModel->Database().StorageKey("sectorIPArray"));

    // FE_SIP_TO_BCTR

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {

      auto ePlacement = (*eit)->AtBarycenter();

      double64 bctrSca(ePlacement.Obtain(keyScalar));
      ScalarVariable bctrSca_;
      ePlacement.Obtain(keyScalar, bctrSca_);

      VectorVariable<2U> bctrVector;
      ePlacement.Obtain(keyVector, bctrVector);
      TensorVariable<2U> bctrTensor;
      ePlacement.Obtain(keyTensor, bctrTensor);
      ArrayVariable bctrArray(3, 0., ANY);
      ePlacement.Obtain(keyArray, bctrArray);

      vector<double64> NN;
      (*eit)->N_AtBaryCenter(NN);


      ScalarVariable scalVal;
      (*eit)->PropertyValueAtBaryCenter(keyScalar, scalVal);
      VectorVariable<2U> vecVal(ANY, ANY, 0., 0.);
      (*eit)->PropertyValueAtBaryCenter(keyVector, vecVal);
      TensorVariable<2U> tenVal(ANY, ANY, 0., 0., 0., 0.);
      (*eit)->PropertyValueAtBaryCenter(keyTensor, tenVal);
      ArrayVariable arrVal(3, 0., ANY);
      (*eit)->PropertyValueAtBaryCenter(keyArray, arrVal);

      //scalar
      _equal(scalVal(), bctrSca, floatTolerance);
      _equal(scalVal(), bctrSca_(), floatTolerance); 

      _equal(vecVal[0], bctrVector[0], floatTolerance);
      _equal(vecVal[1], bctrVector[1], floatTolerance);

      //tensor
      _equal(tenVal(0, 0), bctrTensor(0, 0), floatTolerance);
      _equal(tenVal(0, 1), bctrTensor(0, 1), floatTolerance);
      _equal(tenVal(1, 0), bctrTensor(1, 0), floatTolerance);
      _equal(tenVal(1, 1), bctrTensor(1, 1), floatTolerance);

      //array
      _equal(arrVal[0], bctrArray[0], floatTolerance);
      _equal(arrVal[1], bctrArray[1], floatTolerance);
      _equal(arrVal[2], bctrArray[2], floatTolerance);
    }
  }

  void Placement_Test::FacetIPToElementTest()
  {
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, FACET_INTEGRATION_POINT> keyScalar(mockModel->Database().StorageKey("facetIPScalar"));
    csmp::INDEX<VECTOR, FACET_INTEGRATION_POINT> keyVector(mockModel->Database().StorageKey("facetIPVector"));
    csmp::INDEX<TENSOR, FACET_INTEGRATION_POINT> keyTensor(mockModel->Database().StorageKey("facetIPTensor"));
    csmp::INDEX<ARRAY, FACET_INTEGRATION_POINT> keyArray(mockModel->Database().StorageKey("facetIPArray"));
    // FE_FIP_TO_BCTR

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {

      auto ePlacement = (*eit)->AtBarycenter();

      double64 bctrSca(ePlacement.Obtain(keyScalar));
      ScalarVariable bctrSca_;
      ePlacement.Obtain(keyScalar, bctrSca_);
      VectorVariable<2U> bctrVector;
      ePlacement.Obtain(keyVector, bctrVector);
      TensorVariable<2U> bctrTensor;
      ePlacement.Obtain(keyTensor, bctrTensor);
      ArrayVariable bctrArray(3, 0., ANY);
      ePlacement.Obtain(keyArray, bctrArray);

      ScalarVariable scalVal(ANY, 0.);
      ScalarVariable scalTemp;
      VectorVariable<2U> vecVal(ANY, ANY, 0., 0.);;
      VectorVariable<2U> vecTemp;
      TensorVariable<2U> tenVal(ANY, ANY, 0., 0., 0., 0.);
      TensorVariable<2U> tenTemp;
      ArrayVariable arrVal(3, 0., ANY);
      ArrayVariable arrTemp(3, 0., ANY);

      size_t facetNum((*eit)->Facets());
      size_t IPNum(0);

      for (size_t iFacet(0); iFacet < facetNum; iFacet++) {
        // FacetIP -> ELBCTR: mean value
        for (size_t jIP(0); jIP < (*eit)->IntegrationPointsPerFacet(); jIP++) {
          // scalar
          (*eit)->Read(iFacet, jIP, keyScalar, scalTemp);
          scalVal += scalTemp;

          //vector
          (*eit)->Read(iFacet, jIP, keyVector, vecTemp);
          vecVal += vecTemp;

          //Tensor
          (*eit)->Read(iFacet, jIP, keyTensor, tenTemp);
          tenVal += tenTemp;

          //Array
          (*eit)->Read(iFacet, jIP, keyArray, arrTemp);
          arrVal += arrTemp;

          IPNum++;
        }
      }

      scalVal /= IPNum;
      vecVal  /= IPNum;
      tenVal  /= IPNum;
      arrVal /= IPNum;

      //scalar
      _equal(scalVal(), bctrSca, floatTolerance);
      _equal(scalVal(), bctrSca_(), floatTolerance);

      //vector
      _equal(vecVal[0], bctrVector[0], floatTolerance);
      _equal(vecVal[1], bctrVector[1], floatTolerance);

      //tensor
      _equal(tenVal(0, 0), bctrTensor(0, 0), floatTolerance);
      _equal(tenVal(0, 1), bctrTensor(0, 1), floatTolerance);
      _equal(tenVal(1, 0), bctrTensor(1, 0), floatTolerance);
      _equal(tenVal(1, 1), bctrTensor(1, 1), floatTolerance);

      //array
      _equal(arrVal[0], bctrArray[0], floatTolerance);
      _equal(arrVal[1], bctrArray[1], floatTolerance);
      _equal(arrVal[2], bctrArray[2], floatTolerance);
    }
  }


  // target elementIP
  void Placement_Test::NodeToElementIPTest()
  {
    //FE_NODE_TO_EIP

    //scalar, vector and tensor
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(mockModel->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(mockModel->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(mockModel->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE>  keyArray(mockModel->Database().StorageKey("nodeArray"));

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto ipPlacementCollection = (*eit)->AllElementIntegrationPoints(); // integration point placement 

      size_t ipit(0);
      for (auto ip = ipPlacementCollection.begin(); ip != ipPlacementCollection.end(); ip++) {
    
        double64 ipSca((*ip).Obtain(keyScalar));
        ScalarVariable ipSca_;
        (*ip).Obtain(keyScalar, ipSca_);
        VectorVariable<2U> ipVector;
        (*ip).Obtain(keyVector, ipVector);
        TensorVariable<2U> ipTensor;
        (*ip).Obtain(keyTensor, ipTensor);
        ArrayVariable ipArray(3, 0., ANY);
        (*ip).Obtain(keyArray, ipArray);

        ScalarVariable scalVal(ANY, 0.);
        ScalarVariable scalTemp;
        VectorVariable<2U> vecVal;
        VectorVariable<2U> vecTemp;
        TensorVariable<2U> tenVal;
        TensorVariable<2U> tenTemp;
        ArrayVariable arrVal(3, 0., ANY);
        ArrayVariable arrTemp(3, 0., ANY);
     
        (*eit)->PropertyValueAtIntegrationPoint(keyScalar, ipit, scalVal);
        (*eit)->PropertyValueAtIntegrationPoint(keyVector, ipit, vecVal);
        (*eit)->PropertyValueAtIntegrationPoint(keyTensor, ipit, tenVal);
        (*eit)->PropertyValueAtIntegrationPoint(keyArray, ipit, arrVal);
        ipit++; 
        //scalar
        _equal(scalVal(), ipSca, floatTolerance);
        _equal(scalVal(), ipSca_(), floatTolerance);

        //vector
        _equal(vecVal[0], ipVector[0], floatTolerance);
        _equal(vecVal[1], ipVector[1], floatTolerance);

        //tensor
        _equal(tenVal(0, 0), ipTensor(0, 0), floatTolerance);
        _equal(tenVal(0, 1), ipTensor(0, 1), floatTolerance);
        _equal(tenVal(1, 0), ipTensor(1, 0), floatTolerance);
        _equal(tenVal(1, 1), ipTensor(1, 1), floatTolerance);

        //array
        _equal(arrVal[0], ipArray[0], floatTolerance);
        _equal(arrVal[1], ipArray[1], floatTolerance);
        _equal(arrVal[2], ipArray[2], floatTolerance);
      }
    }
  }




  void Placement_Test::ElementToElementIPTest()
  {
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(mockModel->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(mockModel->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(mockModel->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyArray(mockModel->Database().StorageKey("elementArray"));

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto ipElementPlacement = (*eit)->AllElementIntegrationPoints(); // integration point placement 

      size_t ipit(0);
      for (auto ip = ipElementPlacement.begin(); ip != ipElementPlacement.end(); ip++) {

        double64 ipSca((*ip).Obtain(keyScalar));
        ScalarVariable ipSca_;
        (*ip).Obtain(keyScalar, ipSca_);
        VectorVariable<2U> ipVector;
        (*ip).Obtain(keyVector, ipVector);
        TensorVariable<2U> ipTensor;
        (*ip).Obtain(keyTensor, ipTensor);
        ArrayVariable ipArray(3, 0., ANY);
        (*ip).Obtain(keyArray, ipArray);

        ScalarVariable scalVal(ANY, 0.);
        ScalarVariable scalTemp;
        VectorVariable<2U> vecVal;
        VectorVariable<2U> vecTemp;

        TensorVariable<2U> tenVal;
        TensorVariable<2U> tenTemp;

        ArrayVariable arrVal(3, 0., ANY);
        ArrayVariable arrTemp(3, 0., ANY);

        (*eit)->PropertyValueAtIntegrationPoint(keyScalar, ipit, scalVal);
        (*eit)->PropertyValueAtIntegrationPoint(keyVector, ipit, vecVal);
        (*eit)->PropertyValueAtIntegrationPoint(keyTensor, ipit, tenVal);
        (*eit)->PropertyValueAtIntegrationPoint(keyArray, ipit, arrVal);

        ipit++; 
        //scalar
        _equal(scalVal(), ipSca, floatTolerance);
        _equal(scalVal(), ipSca_(), floatTolerance);

        //vector
        _equal(vecVal[0], ipVector[0], floatTolerance);
        _equal(vecVal[1], ipVector[1], floatTolerance);

        //tensor
        _equal(tenVal(0, 0), ipTensor(0, 0), floatTolerance);
        _equal(tenVal(0, 1), ipTensor(0, 1), floatTolerance);
        _equal(tenVal(1, 0), ipTensor(1, 0), floatTolerance);
        _equal(tenVal(1, 1), ipTensor(1, 1), floatTolerance);

        //array
        _equal(arrVal[0], ipArray[0], floatTolerance);
        _equal(arrVal[1], ipArray[1], floatTolerance);
        _equal(arrVal[2], ipArray[2], floatTolerance);

      }
    }
  }




#if 0

  void Placement_Test::SectorIPToElementIPTest()
  {
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, SECTOR_INTEGRATION_POINT> keyScalar(mockModel->Database().StorageKey("sectorIPScalar"));
    csmp::INDEX<VECTOR, SECTOR_INTEGRATION_POINT> keyVector(mockModel->Database().StorageKey("sectorIPVector"));
    csmp::INDEX<TENSOR, SECTOR_INTEGRATION_POINT> keyTensor(mockModel->Database().StorageKey("sectorIPTensor"));
    csmp::INDEX<ARRAY, SECTOR_INTEGRATION_POINT>  keyArray(mockModel->Database().StorageKey("sectorIPArray"));

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) 
    {
      auto ipElementPlacement = (*eit)->AllElementIntegrationPoints(); // integration point placement 
      auto ipSectorPlacement  = (*eit)->AllSectorIntegrationPoints();
      size_t ipit(0);
      for (auto ip = ipElementPlacement.cbegin(); ip != ipElementPlacement.cend(); ip++) 
      {
        
          double64 ipSca((*ip).Obtain(keyScalar));
          ScalarVariable ipSca_;
          (*ip).Obtain(keyScalar, ipSca_);
          VectorVariable<2U> ipVector;
          (*ip).Obtain(keyVector, ipVector);
          TensorVariable<2U> ipTensor;
          (*ip).Obtain(keyTensor, ipTensor);
          ArrayVariable ipArray;
          (*ip).Obtain(keyArray, ipArray);

          ScalarVariable scalVal(ANY, 0.);
          ScalarVariable scalTemp;
          VectorVariable<2U> vecVal;
          VectorVariable<2U> vecTemp;

          TensorVariable<2U> tenVal;
          TensorVariable<2U> tenTemp;

          ArrayVariable arrVal(2, 0., ANY);
          ArrayVariable arrTemp;

          for (size_t ipSector(0); ipSector < (*eit)->IntegrationPointsPerSector(); ipSector++) {
              //(*eit)->PropertyValueAtIntegrationPoint(keyScalar, ipit, scalVal);
              //(*eit)->PropertyValueAtIntegrationPoint(keyVector, ipit, vecVal);
              //(*eit)->PropertyValueAtIntegrationPoint(keyTensor, ipit, tenVal);
          }

          ipit++; // std::distance

          //scalar
          _equal(scalVal(), ipSca, floatTolerance);
          _equal(scalVal(), ipSca_(), floatTolerance);

          //vector
          _equal(vecVal[0], ipVector[0], floatTolerance);
          _equal(vecVal[1], ipVector[1], floatTolerance);

          //tensor
          _equal(tenVal(0, 0), ipTensor(0, 0), floatTolerance);
          _equal(tenVal(0, 1), ipTensor(0, 1), floatTolerance);
          _equal(tenVal(1, 0), ipTensor(1, 0), floatTolerance);
          _equal(tenVal(1, 1), ipTensor(1, 1), floatTolerance);

          //array
          // _equal(arrVal[0], bctrArray[0], floatTolerance);
          // _equal(arrVal[1], bctrArray[1], floatTolerance);
      }
    }
 }



  void Placement_Test::FacetIPToElementIPTest()
  {
  
    
  }



#endif


  // target sectorIP
  void Placement_Test::NodeToSectorIPTest()
  {
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(mockModel->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(mockModel->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(mockModel->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE> keyArray(mockModel->Database().StorageKey("nodeArray"));

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto fit = (*eit)->AllSectorIntegrationPoints();
      size_t nNode((*eit)->Nodes());
      size_t iSector(0);

      for (auto sPlacement = fit.begin(); sPlacement != fit.end(); sPlacement++)
      {
        for (size_t ip(0); ip < (*eit)->IntegrationPointsPerFacet(); ip++)
        {
          double64 secSca((*sPlacement).Obtain(keyScalar));
          ScalarVariable secSca_;
          (*sPlacement).Obtain(keyScalar, secSca_);

          VectorVariable<2U> secVector;
          (*sPlacement).Obtain(keyVector, secVector);
          TensorVariable<2U> secTensor;
          (*sPlacement).Obtain(keyTensor, secTensor);
          ArrayVariable secArray(3, 0., ANY);
          (*sPlacement).Obtain(keyArray, secArray);


          double64 scalVal(0.);
          VectorVariable<2U> vecVal(ANY, ANY, 0., 0.);;
          VectorVariable<2U> vecTemp;

          TensorVariable<2U> tenVal(ANY, ANY, 0., 0., 0., 0.);
          TensorVariable<2U> tenTemp;

          ArrayVariable arrVal(3, 0., ANY);
          ArrayVariable arrTemp(3, 0., ANY);
          vector<double64> NN(nNode);

          Point<2U> ipCoord((*eit)->FV()->SectorIntegrationPoint(iSector, ip));
          (*eit)->N_At(ipCoord, NN); 
          
          //scalar
          scalVal = (*eit)->PropertyValueAtSectorIntegrationPoint(iSector, ip, keyScalar);
          for (size_t i = 0U; i < (*eit)->Nodes(); i++)
          {
            //vector
            (*eit)->N(i)->Read(keyVector, vecTemp);
            vecVal += vecTemp * NN[i];
            //tensor
            (*eit)->N(i)->Read(keyTensor, tenTemp);
            tenVal += tenTemp * NN[i];
            //array
            (*eit)->N(i)->Read(keyArray, arrTemp);
            arrVal += arrTemp * NN[i];

          }

          //scalar
          _equal(scalVal, secSca, floatTolerance);
          _equal(scalVal, secSca_(), floatTolerance);

          //vector
          _equal(vecVal[0], secVector[0], floatTolerance);
          _equal(vecVal[1], secVector[1], floatTolerance);

          //tensor
          _equal(tenVal(0, 0), secTensor(0, 0), floatTolerance);
          _equal(tenVal(0, 1), secTensor(0, 1), floatTolerance);
          _equal(tenVal(1, 0), secTensor(1, 0), floatTolerance);
          _equal(tenVal(1, 1), secTensor(1, 1), floatTolerance);

          //array
          _equal(arrVal[0], secArray[0], floatTolerance);
          _equal(arrVal[1], secArray[1], floatTolerance);
          _equal(arrVal[2], secArray[2], floatTolerance);

        }
        iSector++;
      }
    }

  }
  


  void Placement_Test::ElementToSectorIPTest()
  {
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(mockModel->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(mockModel->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(mockModel->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyArray(mockModel->Database().StorageKey("elementArray"));

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto ipSectorPlacement = (*eit)->AllSectorIntegrationPoints(); // integration point placement 
      for (auto ip = ipSectorPlacement.begin(); ip != ipSectorPlacement.end(); ip++) {

        double64 ipSca((*ip).Obtain(keyScalar));
        ScalarVariable ipSca_;
        (*ip).Obtain(keyScalar, ipSca_);
        VectorVariable<2U> ipVector;
        (*ip).Obtain(keyVector, ipVector);
        TensorVariable<2U> ipTensor;
        (*ip).Obtain(keyTensor, ipTensor);
        ArrayVariable ipArray(3, 0., ANY);
        (*ip).Obtain(keyArray, ipArray);

        ScalarVariable scalVal(ANY, 0.);
        VectorVariable<2U> vecVal;
        TensorVariable<2U> tenVal;
        ArrayVariable arrVal(3, 0., ANY);


        (*eit)->Read(keyScalar, scalVal);
        (*eit)->Read(keyVector, vecVal);
        (*eit)->Read(keyTensor, tenVal);
        (*eit)->Read(keyArray, arrVal);

        
        //scalar
        _equal(scalVal(), ipSca, floatTolerance);
        _equal(scalVal(), ipSca_(), floatTolerance);

        //vector
        _equal(vecVal[0], ipVector[0], floatTolerance);
        _equal(vecVal[1], ipVector[1], floatTolerance);

        //tensor
        _equal(tenVal(0, 0), ipTensor(0, 0), floatTolerance);
        _equal(tenVal(0, 1), ipTensor(0, 1), floatTolerance);
        _equal(tenVal(1, 0), ipTensor(1, 0), floatTolerance);
        _equal(tenVal(1, 1), ipTensor(1, 1), floatTolerance);

        //array
        _equal(arrVal[0], ipArray[0], floatTolerance);
        _equal(arrVal[1], ipArray[1], floatTolerance);
        _equal(arrVal[2], ipArray[2], floatTolerance);
      }
    }
  }


#if 0
  void Placement_Test::ElementIPToNodeTest() {
  

  }

  void Placement_Test::SectorIPToNodeTest() {
   

  }

  void Placement_Test::FacetIPToNodeTest() {
  
  }

#endif



  // target facetIP
  void Placement_Test::NodeToFacetIPTest()
  {
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(mockModel->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(mockModel->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(mockModel->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE> keyArray(mockModel->Database().StorageKey("nodeArray"));
    // FE_NODE_TO_FCIP

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto fit = (*eit)->AllFacetIntegrationPoints();
      size_t nNode((*eit)->Nodes());
      size_t iFacet(0);

      for (auto fPlacement = fit.begin(); fPlacement != fit.end(); fPlacement++)
      {
        for (size_t ip(0); ip < (*eit)->IntegrationPointsPerFacet(); ip++) 
        {
          double64 facSca((*fPlacement).Obtain(keyScalar));
          ScalarVariable facSca_;
          (*fPlacement).Obtain(keyScalar, facSca_);

          VectorVariable<2U> facVector;
          (*fPlacement).Obtain(keyVector, facVector);
          TensorVariable<2U> facTensor;
          (*fPlacement).Obtain(keyTensor, facTensor);
          ArrayVariable facArray(3, 0., ANY);
          (*fPlacement).Obtain(keyArray, facArray);
          double64 scalVal(0.);
          VectorVariable<2U> vecVal(ANY, ANY, 0., 0.);;
          VectorVariable<2U> vecTemp;

          TensorVariable<2U> tenVal(ANY, ANY, 0., 0., 0., 0.);
          TensorVariable<2U> tenTemp;

          ArrayVariable arrVal(3, 0., ANY);
          ArrayVariable arrTemp(3, 0., ANY);
          vector<double64> NN(nNode);

          Point<2U> ipCoord((*eit)->FV()->FacetIntegrationPoint(iFacet, ip));
          (*eit)->N_At(ipCoord, NN); //TODO: N_AtFacetIntegrationPoint -> void?
          //scalar
          scalVal = (*eit)->PropertyValueAtFacetIntegrationPoint(iFacet, ip, keyScalar);
          for (size_t i = 0U; i < (*eit)->Nodes(); i++)
          {
            
            //vector
            (*eit)->N(i)->Read(keyVector, vecTemp);
            vecVal += vecTemp * NN[i];
            //tensor
            (*eit)->N(i)->Read(keyTensor, tenTemp);
            tenVal += tenTemp * NN[i];
            //array
            (*eit)->N(i)->Read(keyArray, arrTemp);
            arrVal += arrTemp * NN[i];
          }
         
          //scalar
          _equal(scalVal, facSca, floatTolerance);
          _equal(scalVal, facSca_(), floatTolerance);

          //vector
          _equal(vecVal[0], facVector[0], floatTolerance);
          _equal(vecVal[1], facVector[1], floatTolerance);

          //tensor
          _equal(tenVal(0, 0), facTensor(0, 0), floatTolerance);
          _equal(tenVal(0, 1), facTensor(0, 1), floatTolerance);
          _equal(tenVal(1, 0), facTensor(1, 0), floatTolerance);
          _equal(tenVal(1, 1), facTensor(1, 1), floatTolerance);

          //array
          _equal(arrVal[0], facArray[0], floatTolerance);
          _equal(arrVal[1], facArray[1], floatTolerance);
          _equal(arrVal[2], facArray[2], floatTolerance);
        }

        iFacet++;

      }
    }
  }


  void Placement_Test::ElementToFacetIPTest()
  {
    const Region<2U>& region = mockModel->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(mockModel->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(mockModel->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(mockModel->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT> keyArray(mockModel->Database().StorageKey("elementArray"));
    // FE_NODE_TO_FCIP

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto fit = (*eit)->AllFacetIntegrationPoints();
      size_t iFacet(0);

      for (auto fPlacement = fit.begin(); fPlacement != fit.end(); fPlacement++)
      {
        for (size_t ip(0); ip < (*eit)->IntegrationPointsPerFacet(); ip++)
        {
          double64 facSca((*fPlacement).Obtain(keyScalar));
          ScalarVariable facSca_;
          (*fPlacement).Obtain(keyScalar, facSca_);

          VectorVariable<2U> facVector;
          (*fPlacement).Obtain(keyVector, facVector);
          TensorVariable<2U> facTensor;
          (*fPlacement).Obtain(keyTensor, facTensor);
          ArrayVariable facArray(3, 0., ANY);
          (*fPlacement).Obtain(keyArray, facArray);

          ScalarVariable scalVal(ANY, 0.);
          VectorVariable<2U> vecVal(ANY, ANY, 0., 0.);;
          TensorVariable<2U> tenVal(ANY, ANY, 0., 0., 0., 0.);
          ArrayVariable arrVal(3, 0., ANY);

          (*eit)->PropertyValueAtFacetIntegrationPoint(keyScalar, iFacet, ip, scalVal);
          (*eit)->PropertyValueAtFacetIntegrationPoint(keyVector, iFacet, ip, vecVal);
          (*eit)->PropertyValueAtFacetIntegrationPoint(keyTensor, iFacet, ip, tenVal);
          (*eit)->PropertyValueAtFacetIntegrationPoint(keyArray, iFacet, ip, arrVal);

          //scalar
          _equal(scalVal(), facSca, floatTolerance);
          _equal(scalVal(), facSca_(), floatTolerance);

          //vector
          _equal(vecVal[0], facVector[0], floatTolerance);
          _equal(vecVal[1], facVector[1], floatTolerance);

          //tensor
          _equal(tenVal(0, 0), facTensor(0, 0), floatTolerance);
          _equal(tenVal(0, 1), facTensor(0, 1), floatTolerance);
          _equal(tenVal(1, 0), facTensor(1, 0), floatTolerance);
          _equal(tenVal(1, 1), facTensor(1, 1), floatTolerance);

          //array
          _equal(arrVal[0], facArray[0], floatTolerance);
          _equal(arrVal[1], facArray[1], floatTolerance);
          _equal(arrVal[2], facArray[2], floatTolerance);
        }

        iFacet++;

      }
    }

  }

  

  // others

} //csmp