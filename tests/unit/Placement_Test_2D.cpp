#include "Placement_Test_2D.h"
#include "DenseMatrix.h"

using namespace std;

namespace csmp {


  Placement_Test_2D::~Placement_Test_2D() {
    delete _model;
  }

  Placement_Test_2D::Placement_Test_2D() :
    /*
    nodeScalar(1.),
    elementScalar(2.),
    elementIPScalar(3.),
    facetIPScalar(4.),
    sectorIPScalar(5.),
    vx(1.), vy(2.),
    xx(1.), xy(2.), yx(3.), yy(4.),
    tx(1.), ty(2.), tz(3.),
    arx(1.), ary(2.), arz(3.),
    */
    floatTolerance(1e-6)
  {
    
    //this->_model = new  ANSYS_Model2D("mockModel", "mockModelVariable.txt", true, false);
    //this->_model = new  ANSYS_Model2D("WellBore2D", "mockModelVariable.txt");
    //this->_model = new  ANSYS_Model2D("HorFracs2D", "mockModelVariable.txt");
    this->_model = new  ANSYS_Model2D( "square1x1_quad_struct", "PlacementTest_2D-variables.txt");
    
    const Region<2U>& region = _model->Region("Model");

    // initialize the value for 
    // node
    csmp::INDEX<SCALAR, NODE> keyNodeScalar(_model->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyNodeVector(_model->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyNodeTensor(_model->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE>  keyNodeArray(_model->Database().StorageKey("nodeArray")); 
    csmp::INDEX<FLAGGEDARRAY, NODE> keyNodeFlaggedArray(_model->Database().StorageKey("nodeFlaggedArray"));
    
    // element
    csmp::INDEX<SCALAR, ELEMENT> keyElementScalar(_model->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyElementVector(_model->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyElementTensor(_model->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyElementArray(_model->Database().StorageKey("elementArray"));
    csmp::INDEX<FLAGGEDARRAY, ELEMENT>keyElementFlaggedArray(_model->Database().StorageKey("elementFlaggedArray"));

    // elementIP
    csmp::INDEX<SCALAR, ELEMENT_INTEGRATION_POINT> keyElementIPScalar(_model->Database().StorageKey("elementIPScalar"));
    csmp::INDEX<VECTOR, ELEMENT_INTEGRATION_POINT> keyElementIPVector(_model->Database().StorageKey("elementIPVector"));
    csmp::INDEX<TENSOR, ELEMENT_INTEGRATION_POINT> keyElementIPTensor(_model->Database().StorageKey("elementIPTensor"));
    csmp::INDEX<ARRAY, ELEMENT_INTEGRATION_POINT>  keyElementIPArray(_model->Database().StorageKey("elementIPArray"));
    csmp::INDEX<FLAGGEDARRAY, ELEMENT_INTEGRATION_POINT> keyElementIPFlaggedArray(_model->Database().StorageKey("elementIPFlaggedArray"));

    // sectorIP
    csmp::INDEX<SCALAR, SECTOR_INTEGRATION_POINT> keySectorScalar(_model->Database().StorageKey("sectorIPScalar"));
    csmp::INDEX<VECTOR, SECTOR_INTEGRATION_POINT> keySectorVector(_model->Database().StorageKey("sectorIPVector"));
    csmp::INDEX<TENSOR, SECTOR_INTEGRATION_POINT> keySectorTensor(_model->Database().StorageKey("sectorIPTensor"));
    csmp::INDEX<ARRAY, SECTOR_INTEGRATION_POINT>  keySectorArray(_model->Database().StorageKey("sectorIPArray"));
    csmp::INDEX<FLAGGEDARRAY, SECTOR_INTEGRATION_POINT> keySectorFlaggedArray(_model->Database().StorageKey("sectorIPFlaggedArray"));

    // facetIP
    csmp::INDEX<SCALAR, FACET_INTEGRATION_POINT> keyFacetScalar(_model->Database().StorageKey("facetIPScalar"));
    csmp::INDEX<VECTOR, FACET_INTEGRATION_POINT> keyFacetVector(_model->Database().StorageKey("facetIPVector"));
    csmp::INDEX<TENSOR, FACET_INTEGRATION_POINT> keyFacetTensor(_model->Database().StorageKey("facetIPTensor"));
    csmp::INDEX<ARRAY, FACET_INTEGRATION_POINT>  keyFacetArray(_model->Database().StorageKey("facetIPArray"));
    csmp::INDEX<FLAGGEDARRAY, FACET_INTEGRATION_POINT> keyFacetFlaggedArray(_model->Database().StorageKey("facetIPFlaggedArray"));
    
    
    ArrayVariable arrayVar(3, 0., ANY);
    /*
    arrayVar.Component(0, arx);
    arrayVar.Component(1, ary);
    arrayVar.Component(2, arz);
    */
    
    // node
    vector<Node<2U>*>::const_iterator nodeIterator;
    for (nodeIterator = region.NodesBegin(); nodeIterator != region.NodesEnd(); nodeIterator++) {
      Point<2U> coord = (*nodeIterator)->Coordinate();
      nodeScalar = coord.Length();
      vx = xx = coord[0]; 
      vy = yy = coord[1];
      xy = yx = 0.;
      arrayVar.Component(0, coord[0]);
      arrayVar.Component(1, coord[1]); 
      arrayVar.Component(2, coord.Length());      
      (*nodeIterator)->Store(keyNodeScalar, makeScalar(ANY, nodeScalar));
      (*nodeIterator)->Store(keyNodeVector, makeVector(ANY, ANY, vx, vy));
      (*nodeIterator)->Store(keyNodeTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
      (*nodeIterator)->Store(keyNodeArray, arrayVar); 
    }

    // element
    vector<Element<2U>*>::const_iterator eleIter;
    for (eleIter = region.ElementsBegin(); eleIter != region.ElementsEnd(); eleIter++) {
      //barycenter
      Point<2U> byct = (*eleIter)->BaryCenter();
      elementScalar = byct.Length();
      vx = xx = byct[0]; 
      vy = yy = byct[1];
      xy = yx = 0.;
      arrayVar.Component(0, byct[0]);
      arrayVar.Component(1, byct[1]); 
      arrayVar.Component(2, byct.Length());     
      (*eleIter)->Store(keyElementScalar, makeScalar(ANY, elementScalar));
      (*eleIter)->Store(keyElementVector, makeVector(ANY, ANY, vx, vy));
      (*eleIter)->Store(keyElementTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
      (*eleIter)->Store(keyElementArray, arrayVar);
      //integration points
      for (size_t i(0); i < (*eleIter)->IntegrationPoints(); ++i) {
        Point<2U> igpt = (*eleIter)->IntegrationPoint(i);
        elementIPScalar = igpt.Length();
        vx = xx = igpt[0]; 
        vy = yy = igpt[1];
        xy = yx = 0.;
        arrayVar.Component(0, igpt[0]);
        arrayVar.Component(1, igpt[1]); 
        arrayVar.Component(2, igpt.Length());       
        (*eleIter)->Store(i, keyElementIPScalar, makeScalar(ANY, elementIPScalar));
        (*eleIter)->Store(i, keyElementIPVector, makeVector(ANY, ANY, vx, vy));
        (*eleIter)->Store(i, keyElementIPTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
        (*eleIter)->Store(i, keyElementIPArray, arrayVar);
      }
      //sector integration points
      for (size_t i(0); i < (*eleIter)->Sectors(); ++i) {
        for (size_t ip(0); ip < (*eleIter)->IntegrationPointsPerSector(); ip++) {
          Point<2U> ipCoord = (*eleIter)->FV()->SectorIntegrationPoint(i, ip);
          sectorIPScalar = ipCoord.Length();
          vx = xx = ipCoord[0]; 
          vy = yy = ipCoord[1];
          xy = yx = 0.;
          arrayVar.Component(0, ipCoord[0]);
          arrayVar.Component(1, ipCoord[1]); 
          arrayVar.Component(2, ipCoord.Length());      
          (*eleIter)->Store(i, ip, keySectorScalar, makeScalar(ANY, sectorIPScalar));
          (*eleIter)->Store(i, ip, keySectorVector, makeVector(ANY, ANY, vx, vy));
          (*eleIter)->Store(i, ip, keySectorTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
          (*eleIter)->Store(i, ip, keySectorArray, arrayVar);
        }
      }
      //facet integration points
      for (size_t i(0); i < (*eleIter)->Facets(); ++i) {
        for (size_t ip(0); ip < (*eleIter)->IntegrationPointsPerFacet(); ip++) {
          Point<2U> ipCoord = (*eleIter)->FV()->FacetIntegrationPoint(i, ip);
          facetIPScalar = ipCoord.Length();
          vx = xx = ipCoord[0]; 
          vy = yy = ipCoord[1];
          xy = yx = 0.;
          arrayVar.Component(0, ipCoord[0]);
          arrayVar.Component(1, ipCoord[1]); 
          arrayVar.Component(2, ipCoord.Length());        
          (*eleIter)->Store(i, ip, keyFacetScalar, makeScalar(ANY, facetIPScalar));
          (*eleIter)->Store(i, ip, keyFacetVector, makeVector(ANY, ANY, vx, vy));
          (*eleIter)->Store(i, ip, keyFacetTensor, makeTensor(ANY, ANY, xx, xy, yx, yy));
          (*eleIter)->Store(i, ip, keyFacetArray, arrayVar);
        }
      }
    }
  }


  void Placement_Test_2D::run() {
    //FiniteElementPlacement tests
    NodeToElementTest(); 
    std::cout<<"Passed NodeToElementTest()"<<std::endl;
    ElementToElementTest(); 
    std::cout<<"Passed ElementToElementTest()"<<std::endl;
    ElementIPToElementTest(); 
    std::cout<<"Passed ElementIPToElementTest()"<<std::endl;
    SectorIPToElementTest(); 
    std::cout<<"Passed SectorIPToElementTest()"<<std::endl;
    FacetIPToElementTest(); 
    std::cout<<"Passed FacetIPToElementTest()"<<std::endl;
    NodeToElementIPTest(); 
    std::cout<<"Passed NodeToElementIPTest()"<<std::endl;
    ElementToElementIPTest();
    std::cout<<"Passed ElementToElementIPTest())"<<std::endl;
    NodeToFacetIPTest();
    std::cout<<"Passed NodeToFacetIPTest()"<<std::endl;
    ElementToFacetIPTest();
    std::cout<<"Passed ElementToFacetIPTest()"<<std::endl;
    NodeToSectorIPTest();
    std::cout<<"Passed NodeToSectorIPTest()"<<std::endl;
    ElementToSectorIPTest();
    std::cout<<"Passed ElementToSectorIPTest()"<<std::endl;
    ElementIPToElementIPTest(); 
    std::cout<<"Passed ElementIPToElementIPTest()"<<std::endl;
    FacetIPToFacetIPTest(); 
    std::cout<<"Passed FacetIPToFacetIPTest()"<<std::endl;
    SectorIPToSectorIPTest();
    std::cout<<"Passed SectorIPToSectorIPTest()"<<std::endl;
    NodeToNodeTest(); 
    std::cout<<"Passed NodeToNodeTest()"<<std::endl;
    GradientTest(); 
    std::cout<<"Passed GradientTest()"<<std::endl;
    SectorVolumeTest();
    std::cout<<"Passed SectorVolumeTest()"<<std::endl;
    FacetNormalTest();
    std::cout<<"Passed FacetNormalTest()"<<std::endl;    
    
    //FiniteVolumePlacement tests
    FV_ElementToElementTest();
    std::cout<<"Passed FV_ElementToElementTest()"<<std::endl;
    FV_FacetAreaAndNormalTest();  
    std::cout<<"Passed FV_FacetAreaAndNormalTest()"<<std::endl; 
    FV_ProjectOntoFacetNormalTest();
    std::cout<<"Passed FV_ProjectOntoFacetNormalTest()"<<std::endl;
    FV_SectorVolumeTest();
    std::cout<<"Passed FV_SectorVolumeTest()"<<std::endl;
    FV_InsideOutsideNodeTest();
    std::cout<<"Passed FV_InsideOutsideNodeTest()"<<std::endl;
    FV_GradientTest();
    std::cout<<"Passed FV_GradientTest()"<<std::endl;    
  }


  // target element

  void Placement_Test_2D::NodeToElementTest()
  {
    //scalar, vector and tensor
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(_model->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(_model->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE>  keyArray(_model->Database().StorageKey("nodeArray"));
   
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
  void Placement_Test_2D::ElementToElementTest()
  {
    //scalar and vector
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(_model->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(_model->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(_model->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyArray(_model->Database().StorageKey("elementArray"));

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

  void Placement_Test_2D::ElementIPToElementTest()
  {

    //scalar, vector and tensor
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT_INTEGRATION_POINT> keyScalar(_model->Database().StorageKey("elementIPScalar"));
    csmp::INDEX<VECTOR, ELEMENT_INTEGRATION_POINT> keyVector(_model->Database().StorageKey("elementIPVector"));
    csmp::INDEX<TENSOR, ELEMENT_INTEGRATION_POINT> keyTensor(_model->Database().StorageKey("elementIPTensor"));
    csmp::INDEX<ARRAY, ELEMENT_INTEGRATION_POINT>  keyArray(_model->Database().StorageKey("elementIPArray"));

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


  void Placement_Test_2D::SectorIPToElementTest()
  {

    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, SECTOR_INTEGRATION_POINT> keyScalar(_model->Database().StorageKey("sectorIPScalar"));
    csmp::INDEX<VECTOR, SECTOR_INTEGRATION_POINT> keyVector(_model->Database().StorageKey("sectorIPVector"));
    csmp::INDEX<TENSOR, SECTOR_INTEGRATION_POINT> keyTensor(_model->Database().StorageKey("sectorIPTensor"));
    csmp::INDEX<ARRAY, SECTOR_INTEGRATION_POINT>  keyArray(_model->Database().StorageKey("sectorIPArray"));

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

  void Placement_Test_2D::FacetIPToElementTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, FACET_INTEGRATION_POINT> keyScalar(_model->Database().StorageKey("facetIPScalar"));
    csmp::INDEX<VECTOR, FACET_INTEGRATION_POINT> keyVector(_model->Database().StorageKey("facetIPVector"));
    csmp::INDEX<TENSOR, FACET_INTEGRATION_POINT> keyTensor(_model->Database().StorageKey("facetIPTensor"));
    csmp::INDEX<ARRAY, FACET_INTEGRATION_POINT> keyArray(_model->Database().StorageKey("facetIPArray"));
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
  void Placement_Test_2D::NodeToElementIPTest()
  {
    //FE_NODE_TO_EIP

    //scalar, vector and tensor
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(_model->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(_model->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE>  keyArray(_model->Database().StorageKey("nodeArray"));

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto ipPlacementCollection = (*eit)->AllElementIntegrationPoints(); // integration point placement 
      size_t ipit(0);
      for (auto ip = ipPlacementCollection.begin(); ip != ipPlacementCollection.end(); ip++) {
        //cout<<"ipit = "<<ipit<<endl;
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



  void Placement_Test_2D::ElementToElementIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(_model->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(_model->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(_model->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyArray(_model->Database().StorageKey("elementArray"));

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

  void Placement_Test_2D::SectorIPToElementIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, SECTOR_INTEGRATION_POINT> keyScalar(_model->Database().StorageKey("sectorIPScalar"));
    csmp::INDEX<VECTOR, SECTOR_INTEGRATION_POINT> keyVector(_model->Database().StorageKey("sectorIPVector"));
    csmp::INDEX<TENSOR, SECTOR_INTEGRATION_POINT> keyTensor(_model->Database().StorageKey("sectorIPTensor"));
    csmp::INDEX<ARRAY, SECTOR_INTEGRATION_POINT>  keyArray(_model->Database().StorageKey("sectorIPArray"));

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



  void Placement_Test_2D::FacetIPToElementIPTest()
  {
  
    
  }



#endif


  // target sectorIP
  void Placement_Test_2D::NodeToSectorIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(_model->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(_model->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE> keyArray(_model->Database().StorageKey("nodeArray"));

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
  


  void Placement_Test_2D::ElementToSectorIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(_model->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(_model->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(_model->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyArray(_model->Database().StorageKey("elementArray"));

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
  void Placement_Test_2D::ElementIPToNodeTest() {
  

  }

  void Placement_Test_2D::SectorIPToNodeTest() {
   

  }

  void Placement_Test_2D::FacetIPToNodeTest() {
  
  }

#endif



  // target facetIP
  void Placement_Test_2D::NodeToFacetIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(_model->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(_model->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE> keyArray(_model->Database().StorageKey("nodeArray"));
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


  void Placement_Test_2D::ElementToFacetIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(_model->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(_model->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(_model->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT> keyArray(_model->Database().StorageKey("elementArray"));
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
  
  
  
  void Placement_Test_2D::ElementIPToElementIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT_INTEGRATION_POINT> keyScalar(_model->Database().StorageKey("elementIPScalar"));
    csmp::INDEX<VECTOR, ELEMENT_INTEGRATION_POINT> keyVector(_model->Database().StorageKey("elementIPVector"));
    csmp::INDEX<TENSOR, ELEMENT_INTEGRATION_POINT> keyTensor(_model->Database().StorageKey("elementIPTensor"));
    csmp::INDEX<ARRAY, ELEMENT_INTEGRATION_POINT>  keyArray(_model->Database().StorageKey("elementIPArray"));
    
    //FE_READ_EIP
    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      
      std::vector<ScalarVariable> scalValVector;
      std::vector<VectorVariable<2U>> vecValVector;
      std::vector<TensorVariable<2U>> tenValVector;
      std::vector<ArrayVariable> arrValVector;
      
      (*eit)->IntegrationPointPropertyVector(keyScalar, scalValVector );
      (*eit)->IntegrationPointPropertyVector(keyVector, vecValVector );
      (*eit)->IntegrationPointPropertyVector(keyTensor, tenValVector );
      (*eit)->IntegrationPointPropertyVector(keyArray, arrValVector );
    
      
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
        
        ScalarVariable scalVal = scalValVector[ipit];
        VectorVariable<2U> vecVal = vecValVector[ipit];
        TensorVariable<2U> tenVal = tenValVector[ipit];
        ArrayVariable arrVal = arrValVector[ipit];
        
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
  


  void Placement_Test_2D::FacetIPToFacetIPTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, FACET_INTEGRATION_POINT> keyScalar(_model->Database().StorageKey("facetIPScalar"));
    csmp::INDEX<VECTOR, FACET_INTEGRATION_POINT> keyVector(_model->Database().StorageKey("facetIPVector"));
    csmp::INDEX<TENSOR, FACET_INTEGRATION_POINT> keyTensor(_model->Database().StorageKey("facetIPTensor"));
    csmp::INDEX<ARRAY, FACET_INTEGRATION_POINT> keyArray(_model->Database().StorageKey("facetIPArray"));
    
    // FE_READ_FIP
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
          
          
          ScalarVariable scalVal;
          VectorVariable<2U> vecVal;
          TensorVariable<2U> tenVal;
          ArrayVariable arrVal(3, 0., ANY);
          
          (*eit)->Read(iFacet, ip, keyScalar, scalVal);
          (*eit)->Read(iFacet, ip, keyVector, vecVal);
          (*eit)->Read(iFacet, ip, keyTensor, tenVal);
          (*eit)->Read(iFacet, ip, keyArray, arrVal);
         
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
  



  void Placement_Test_2D::SectorIPToSectorIPTest()
  {

    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, SECTOR_INTEGRATION_POINT> keyScalar(_model->Database().StorageKey("sectorIPScalar"));
    csmp::INDEX<VECTOR, SECTOR_INTEGRATION_POINT> keyVector(_model->Database().StorageKey("sectorIPVector"));
    csmp::INDEX<TENSOR, SECTOR_INTEGRATION_POINT> keyTensor(_model->Database().StorageKey("sectorIPTensor"));
    csmp::INDEX<ARRAY, SECTOR_INTEGRATION_POINT>  keyArray(_model->Database().StorageKey("sectorIPArray"));

    // FE_READ_SIP
    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      auto sit = (*eit)->AllSectorIntegrationPoints();

      size_t iSector(0);
      for (auto sPlacement = sit.begin(); sPlacement != sit.end(); sPlacement++)
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


          ScalarVariable scalVal;
          VectorVariable<2U> vecVal;
          TensorVariable<2U> tenVal;
          ArrayVariable arrVal(3, 0., ANY);


          (*eit)->Read(iSector, ip, keyScalar, scalVal);
          (*eit)->Read(iSector, ip, keyVector, vecVal);
          (*eit)->Read(iSector, ip, keyTensor, tenVal);
          (*eit)->Read(iSector, ip, keyArray, arrVal);
          

          //scalar
          _equal(scalVal(), secSca, floatTolerance);
          _equal(scalVal(), secSca_(), floatTolerance);

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



  void Placement_Test_2D::NodeToNodeTest()
  {
    //scalar, vector and tensor
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    csmp::INDEX<VECTOR, NODE> keyVector(_model->Database().StorageKey("nodeVector"));
    csmp::INDEX<TENSOR, NODE> keyTensor(_model->Database().StorageKey("nodeTensor"));
    csmp::INDEX<ARRAY, NODE>  keyArray(_model->Database().StorageKey("nodeArray"));
   
    // FE_READ_NODE
    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++ ) {
    
      std::vector<ScalarVariable> scalValVector;
      std::vector<VectorVariable<2U>> vecValVector;
      std::vector<TensorVariable<2U>> tenValVector;
      std::vector<ArrayVariable> arrValVector;
      
      (*eit)->NodePropertyVector(keyScalar, scalValVector );
      (*eit)->NodePropertyVector(keyVector, vecValVector );
      (*eit)->NodePropertyVector(keyTensor, tenValVector );
      (*eit)->NodePropertyVector(keyArray, arrValVector );   
    
    
      auto nit = (*eit)->AllNodes();
      size_t iNode(0);      
      for (auto nPlacement = nit.begin(); nPlacement != nit.end(); nPlacement++)
      {      
        double64 ndSca((*nPlacement).Obtain(keyScalar));
        ScalarVariable ndSca_;
        (*nPlacement).Obtain(keyScalar, ndSca_);

        VectorVariable<2U> ndVector;
        (*nPlacement).Obtain(keyVector, ndVector);
        TensorVariable<2U> ndTensor;
        (*nPlacement).Obtain(keyTensor, ndTensor);
        ArrayVariable ndArray(3, 0., ANY);
        (*nPlacement).Obtain(keyArray, ndArray);


        ScalarVariable scalVal = scalValVector[iNode];
        VectorVariable<2U> vecVal = vecValVector[iNode];
        TensorVariable<2U> tenVal = tenValVector[iNode];
        ArrayVariable arrVal = arrValVector[iNode];
        
        iNode++; 
        
        //scalar
        _equal(scalVal(), ndSca, floatTolerance);
        _equal(scalVal(), ndSca_(), floatTolerance);

        //vector
        _equal(vecVal[0], ndVector[0], floatTolerance);
        _equal(vecVal[1], ndVector[1], floatTolerance);

        //tensor
        _equal(tenVal(0, 0), ndTensor(0, 0), floatTolerance);
        _equal(tenVal(0, 1), ndTensor(0, 1), floatTolerance);
        _equal(tenVal(1, 0), ndTensor(1, 0), floatTolerance);
        _equal(tenVal(1, 1), ndTensor(1, 1), floatTolerance);


        //array
        _equal(arrVal[0], ndArray[0], floatTolerance);
        _equal(arrVal[1], ndArray[1], floatTolerance);
        _equal(arrVal[2], ndArray[2], floatTolerance);

      }
    }
  }



  void Placement_Test_2D::GradientTest()
  {
    //scalar
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++ ) { 
    
      auto ePlacement =(* eit)->AtBarycenter();
      Point<2U> gradP = ePlacement.Gradient(keyScalar);
      
      std::vector<double64> gradV(2U, 0.);
      DenseMatrix<DM_MIN>    DN;
      (* eit)->dN_AtBaryCenter( DN );
        
      for ( size_t j=0U; j<(* eit)->Nodes(); j++ ) {
        double64 v = (* eit)->N(j)->Read(keyScalar);
        for ( size_t k=0U; k<2U; k++ ) gradV[k] += DN(k,j) * v;
      }
        
      _equal(gradP[0], gradV[0], floatTolerance);
      _equal(gradP[1], gradV[1], floatTolerance);

    }
  }        



  void Placement_Test_2D::SectorVolumeTest()
  {
    //scalar
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    
    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {

      size_t iSector(0);
      for (auto sPlacement : (*eit)->AllSectorIntegrationPoints()) 
      {    
        double64 secVol = sPlacement.SectorVolume();
        double64 sector_volume = (*eit)->SectorVolume(iSector);
        iSector++;
        
        _equal(secVol, sector_volume, floatTolerance);
      }
    }
  } 



  void Placement_Test_2D::FacetNormalTest()
  {
    const Region<2U>& region = _model->Region("Model");

    for (auto eit = region.ElementsBegin(); eit != region.ElementsEnd(); eit++) {
      size_t iFacet(0);
      for (auto fip : (*eit)->AllFacetIntegrationPoints())
      {
        Point<2U> fn = fip.FacetNormal();
        Point<2U> facet_normal = (*eit)->FacetNormal(iFacet);
        _equal(fn[0], facet_normal[0], floatTolerance);
        _equal(fn[1], facet_normal[1], floatTolerance);
        iFacet++;
       }
     }
   }



  void Placement_Test_2D::FV_ElementToElementTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, ELEMENT> keyScalar(_model->Database().StorageKey("elementScalar"));
    csmp::INDEX<VECTOR, ELEMENT> keyVector(_model->Database().StorageKey("elementVector"));
    csmp::INDEX<TENSOR, ELEMENT> keyTensor(_model->Database().StorageKey("elementTensor"));
    csmp::INDEX<ARRAY, ELEMENT>  keyArray(_model->Database().StorageKey("elementArray"));

    // FV_READ_ELMT
    for (auto nit = region.NodesBegin(); nit != region.NodesEnd(); nit++) {
      
      size_t iElement (0);
      for (auto ePlacement : (*nit)->AllElements()) {
     
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
        VectorVariable<2U> vecVal(ANY, ANY, 0., 0.);;
        TensorVariable<2U> tenVal(ANY, ANY, 0., 0., 0., 0.);
        ArrayVariable arrVal(3, 0., ANY);

        (*nit)->Parent(iElement)->Read(keyScalar, scalVal);
        (*nit)->Parent(iElement)->Read(keyVector, vecVal);
        (*nit)->Parent(iElement)->Read(keyTensor, tenVal);
        (*nit)->Parent(iElement)->Read(keyArray, arrVal);
        iElement++;
       
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
  }  


  void Placement_Test_2D::FV_FacetAreaAndNormalTest()
  {
    const Region<2U>& region = _model->Region("Model");

    for (auto nit = region.NodesBegin(); nit != region.NodesEnd(); nit++) {
      
      std::vector<double64> faVec;
      std::vector<Point<2U>> fnVec;
      for (auto fip : (*nit)->AllFacetIntegrationPoints()) {
        double64 fa = fip.FacetArea();
        faVec.push_back(fa);
        Point<2U> fn = fip.FacetNormal();
        fnVec.push_back(fn);
      }
      
      std::vector<double64> facetAreaVec;
      std::vector<Point<2U>> facetNormalVec;
      const size_t parent_elements((*nit)->Parents());      
      for ( size_t i=0U; i<parent_elements; ++i ) {
         const Element<2U>* eptr = (*nit)->Parent(i);
         const size_t sector_node = (*nit)->ParentNodeNumber(i);
         for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
           const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
           double64 facetArea = eptr->FacetArea(facet);
           facetAreaVec.push_back(facetArea);
           Point<2U> facetNormal = eptr->FacetNormal(facet);
           facetNormalVec.push_back(facetNormal);
         }
      }

      const size_t vecSize = faVec.size();
      for (size_t n=0U;n<vecSize;n++) {
        double64 fa = faVec[n];
        double64 facet_area = facetAreaVec[n];
        Point<2U> fn = fnVec[n];
        Point<2U> facet_normal = facetNormalVec[n];
        _equal(fa, facet_area, floatTolerance);
        _equal(fn[0], facet_normal[0], floatTolerance);
        _equal(fn[1], facet_normal[1], floatTolerance);
      }
      
    }
  }  




  void Placement_Test_2D::FV_ProjectOntoFacetNormalTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<VECTOR, ELEMENT> keyVector(_model->Database().StorageKey("elementVector"));

    for (auto nit = region.NodesBegin(); nit != region.NodesEnd(); nit++) {
      
      std::vector<double64> proVec;
      for (auto fip : (*nit)->AllFacetIntegrationPoints()) {
        VectorVariable<2U> vecVal;
        fip.Obtain(keyVector, vecVal);
        double64 pro = fip.ProjectOntoFacetNormal(vecVal);
        proVec.push_back(pro);
      }
      
      std::vector<double64> projectVec;
      const size_t parent_elements((*nit)->Parents());      
      for ( size_t i=0U; i<parent_elements; ++i ) {
         const Element<2U>* eptr = (*nit)->Parent(i);
         const size_t sector_node = (*nit)->ParentNodeNumber(i);
         for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
           const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
           double64 project = eptr->ProjectionOnFacetNormal(facet,keyVector);
           projectVec.push_back(project);
         }
      }
 
      const size_t vecSize = proVec.size();
      for (size_t n=0U;n<vecSize;n++) {
        double64 pro = proVec[n];
        double64 project = projectVec[n];
        _equal(pro, project, floatTolerance);
      }
      
    }
  }  
  

  void Placement_Test_2D::FV_SectorVolumeTest()
  {
    const Region<2U>& region = _model->Region("Model");
    
    for (auto nit = region.NodesBegin(); nit != region.NodesEnd(); nit++) {

      std::vector<double64> secVolVec;
      for (auto sPlacement : (*nit)->AllSectorIntegrationPoints()) 
      {    
        double64 secVol = sPlacement.SectorVolume();
        secVolVec.push_back(secVol);               
      }
     
      std::vector<double64> sectorVolumeVec;
      const size_t parent_elements((*nit)->Parents());      
      for ( size_t i=0U; i<parent_elements; ++i ) {
         const Element<2U>* eptr = (*nit)->Parent(i);
         const size_t sector_node = (*nit)->ParentNodeNumber(i);
         double64 sectorVolume = eptr->SectorVolume(sector_node);  
         sectorVolumeVec.push_back(sectorVolume);  
      }
      
      const size_t vecSize = secVolVec.size();
      for (size_t n=0U;n<vecSize;n++) {
        double64 secVol = secVolVec[n];
        double64 sectorVolume = sectorVolumeVec[n];
        _equal(secVol, sectorVolume, floatTolerance);
      }   
    }   
  } 
  

  void Placement_Test_2D::FV_InsideOutsideNodeTest()
  {
    const Region<2U>& region = _model->Region("Model");
    for (auto nit = region.NodesBegin(); nit != region.NodesEnd(); nit++) {
      
      std::vector<size_t> inVec, outVec, signVec;
      for (auto fip : (*nit)->AllFacetIntegrationPoints()) {
        size_t inside = fip.InsideNode().NodeIdx();
        size_t outside = fip.OutsideNode().NodeIdx();
        size_t sign = fip.FromInside() ? 1 : -1;
        inVec.push_back(inside);
        outVec.push_back(outside);
        signVec.push_back(sign);
      }
      
      std::vector<size_t> insideNodeVec, outsideNodeVec, nodeSignVec;
      const size_t parent_elements((*nit)->Parents());      
      for ( size_t i=0U; i<parent_elements; ++i ) {
         const Element<2U>* eptr = (*nit)->Parent(i);
         const size_t sector_node = (*nit)->ParentNodeNumber(i);
         for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
           const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
           size_t inside_nid(eptr->FV()->InsideNode(facet));
           size_t outside_nid(eptr->FV()->OutsideNode(facet));
           size_t insideNode = eptr->N(inside_nid)->Idx();
           size_t outsideNode = eptr->N(outside_nid)->Idx();
           size_t nodeSign =( sector_node == inside_nid) ? 1 : -1;
           insideNodeVec.push_back(insideNode);
           outsideNodeVec.push_back(outsideNode);
           nodeSignVec.push_back(nodeSign);
         }
      }
 
      const size_t vecSize = inVec.size();
      for (size_t n=0U;n<vecSize;n++) {
        size_t inside = inVec[n];
        size_t outside = outVec[n];
        size_t insideNode = insideNodeVec[n];
        size_t outsideNode = outsideNodeVec[n];  
        size_t sign = signVec[n];
        size_t nodeSign = nodeSignVec[n];      
        _equal(inside, insideNode, floatTolerance);
        _equal(outside, outsideNode, floatTolerance);
        _equal(sign, nodeSign, floatTolerance);
      }
      
    }
  }


  void Placement_Test_2D::FV_GradientTest()
  {
    const Region<2U>& region = _model->Region("Model");
    csmp::INDEX<SCALAR, NODE> keyScalar(_model->Database().StorageKey("nodeScalar"));
    
    for (auto nit = region.NodesBegin(); nit != region.NodesEnd(); nit++ ) { 
      std::vector<Point<2U>> gradVec;
      for (auto fip : (*nit)->AllFacetIntegrationPoints()) {
        Element<2U>& e = fip.TheElement();
        Point<2U> grad = fip.Gradient(keyScalar);
        gradVec.push_back(grad);
      }
      
      std::vector<Point<2U>> gradientVec;
      const size_t parent_elements((*nit)->Parents());      
      for ( size_t i=0U; i<parent_elements; ++i ) {
         const Element<2U>* eptr = (*nit)->Parent(i);
         const size_t sector_node = (*nit)->ParentNodeNumber(i);
         for ( size_t j=0U; j<eptr->FV()->FacetsPerSector(sector_node); ++j ) {
           const size_t facet = eptr->FV()->FacetSurroundingSector( sector_node, j );
           for (size_t ip(0); ip < eptr->IntegrationPointsPerFacet(); ip++) {
             DenseMatrix<DM_MIN> DN;
             if (eptr->IsLineElement()) {
               eptr->dN_AtBaryCenter(DN); //dN_At function not support line element yet
             } else {
               Point<2U> coord = eptr->FV()->FacetIntegrationPoint(facet, ip);
               eptr->dN_At (coord, DN);
             }
             
             Point<2U> gradient (0.);
             for ( size_t j=0U; j<eptr->Nodes(); j++ ) {
               double64 v = eptr->N(j)->Read(keyScalar);
               for ( size_t k=0U; k<2U; k++ ) gradient[k] += DN(k,j) * v;
             }
             gradientVec.push_back(gradient);
           }
         }
      }             
      
      const size_t vecSize = gradVec.size();
      for (size_t n=0U;n<vecSize;n++) {      
        Point<2U> grad = gradVec[n];
        Point<2U> gradient = gradientVec[n];
        _equal(grad[0], gradient[0], floatTolerance);
        _equal(grad[1], gradient[1], floatTolerance);

      }
    }
  }  
  

  // others

} //csmp
