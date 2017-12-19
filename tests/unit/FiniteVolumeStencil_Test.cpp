#include "FiniteVolumeStencil_Test.h"
#include "FiniteElement.h"

#include "IsoparametricLinearLineElement.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPrism.h"
#include "AP_algebraUtilities.h"

using namespace std;

namespace csmp {

/**  Method:
 
 
FiniteVolumeStencil_Test::FiniteVolumeStencil_Test( ) 
 
 
Description: 

Constructor. This function initializes the objects which will be used during the test.

@section arguments Input Arguments 
None.
@section application Application

The idea of this function is to specify the characteristics of the objects to be tested. 

 
tested: is a test function*/
FiniteVolumeStencil_Test::FiniteVolumeStencil_Test( bool verbose )
  : verbose_(verbose)
 {
    const size_t dim(3U);
    fvs_.reserve(7U);
    vecFEs_.reserve(7U);
    
    // initialize the vector of stencils
    //linear elements
    // 1,2,3D
    fvs_.push_back( FiniteVolumeStencil<3>( "ISOPARAMETRIC_LINEAR_BAR" ) );
    // 2D and 3D 
    fvs_.push_back( FiniteVolumeStencil<3>( "ISOPARAMETRIC_LINEAR_TRIANGLE" ) );
    fvs_.push_back( FiniteVolumeStencil<3>( "ISOPARAMETRIC_LINEAR_QUADRILATERAL" ) );
    // 3D only
    fvs_.push_back( FiniteVolumeStencil<3>( "ISOPARAMETRIC_LINEAR_TETRAHEDRON" ) );
    fvs_.push_back( FiniteVolumeStencil<3>( "ISOPARAMETRIC_LINEAR_PYRAMID" ) );
    fvs_.push_back( FiniteVolumeStencil<3>( "ISOPARAMETRIC_LINEAR_HEXAHEDRON" ) );
    fvs_.push_back( FiniteVolumeStencil<3>( "ISOPARAMETRIC_LINEAR_PRISM" ) );
     
    /*finite volume stencil 0: ISOPARAMETRIC_LINEAR_BAR*/
    vecFEs_.push_back(new IsoparametricLinearLineElement(dim));
    /*finite volume stencil 1: ISOPARAMETRIC_LINEAR_TRIANGLE*/
    vecFEs_.push_back(new IsoparametricLinearTriangle(dim));
    /*finite volume stencil 2: ISOPARAMETRIC_LINEAR_QUADRILATERAL*/
    vecFEs_.push_back(new IsoparametricLinearQuadrilateral(dim));
    /*finite volume stencil 3: ISOPARAMETRIC_LINEAR_TETRAHEDRON*/
    vecFEs_.push_back(new IsoparametricLinearTetrahedron());
    /*finite volume stencil 4: ISOPARAMETRIC_LINEAR_PYRAMID*/
    vecFEs_.push_back(new IsoparametricLinearPyramid());
    /*finite volume stencil 5: ISOPARAMETRIC_LINEAR_HEXAHEDRON*/
    vecFEs_.push_back(new IsoparametricLinearHexahedron());
	  /*finite volume stencil 6: ISOPARAMETRIC_LINEAR_PRISM*/
  	vecFEs_.push_back(new IsoparametricLinearPrism());
  	
  	//load reference coords as physical coords.
    DenseMatrix<DM_MIN> matCoords;
    for(vector<FiniteElement*>::const_iterator vIterFEs = vecFEs_.begin(); vIterFEs != vecFEs_.end(); vIterFEs++)
    {
      (*vIterFEs)->ReferenceCoordinates( matCoords );
	    for( size_t i = 0; i < (*vIterFEs)->Nodes(); i++ )
	     for( size_t d = 0; d < dim; d++ )
	     {
	       const double64 v(matCoords(i,d)); 
	       (*vIterFEs)->XYZ(i, d, isnan(v)?0.:v);
	     }
	  }	 
 }
 

FiniteVolumeStencil_Test::~FiniteVolumeStencil_Test()
 {
 	//cleanup finite element objects
  	vector<FiniteElement*>::iterator vDelIterFEs;
  	for( vDelIterFEs = vecFEs_.begin(); vDelIterFEs != vecFEs_.end(); vDelIterFEs++)
  	{
  		delete *vDelIterFEs;
  	}
 }
 
 
void FiniteVolumeStencil_Test::run() // runs all the tests for the class (register other methods)
 {
  if ( verbose_ ) cout << "Starting...FiniteVolumeStencil_Test::run()" << endl;
 	
 	displayReferenceCoordinates();
	
	facetAndSectorNumbersTest(); // are they right for all element types
	
	orientationAndLengthOfNormalsTest();
     
    normalTransformationTest();
	
	weightsAndFacetIntegrationPointsTest(); //does not apply for pyramids
     
	sectorFacetConnectivityTest(); // facets that delimit a certain sector for all the element types
	
	weightsOfSectorIntegrationPointsTest(); // argument list could deal with multiple integration points
	
	shapeFunctionsTest();
	
	shapeFunctionDerivativesTest();
	
	if ( verbose_ ) cout << "Done...FiniteVolumeStencil_Test::run()" << endl;
 }
 
 
 
 
 
/** void FiniteVolumeStencil_Test::displayReferenceCoordinates()

    Description: 
    This test function displays the reference coordinates of each type of element.

    @section application Application

    The idea of this function is to provide insight of how elements are being constructed.
     
    tested: is a test function.
*/
void FiniteVolumeStencil_Test::displayReferenceCoordinates()
{
  if ( verbose_ ) cout << "TESTING: displayReferenceCoordinates()" << endl;
 
	for(vector<FiniteElement*>::const_iterator vIterFEs = vecFEs_.begin(); vIterFEs != vecFEs_.end(); vIterFEs++)
  {
  	//print type
    if ( verbose_ ) cout << "Type: " << parseFiniteElementType( (*vIterFEs)->ElementType() ) << endl;
  	DenseMatrix<DM_MIN> matCoords;
		(*vIterFEs)->ReferenceCoordinates(matCoords);
	  if ( verbose_ ) matCoords.Out();
	}	
}

 
/**  Method:
 
 
void FiniteVolumeStencil_Test::facetAndSectorNumbersTest()
 
 

Description: 
Tests for each finite element type supported if the basic query functionality of the class is correct.
@section arguments Input Arguments 


@section application Application
The idea is to test the following member functions:
  	size_t  Facets() const;
    size_t  Sectors() const;
    size_t  IntegrationPointsPerFacet() const;
    size_t  IntegrationPointsPerSector() const;
    void       FacetEdgeNodes( size_t facet, size_t& inside_node, size_t& outside_node ) const;
  	FV_FACET_TYPE FacetType() const;

 
tested: is a test funtion*/
void FiniteVolumeStencil_Test::facetAndSectorNumbersTest() // are they right for all element types
 {
    if ( verbose_ ) {
        cout << "TESTING: facetAndSectorNumbersTest()" << endl;
        cout << vecFEs_.size() << " vs " << fvs_.size();
 	    }
  	assert(vecFEs_.size() == fvs_.size());
  	
  	vector<size_t> vecNodeIdsFEs;
  	
  	//for each finite element type vs finite volume stencil type
  	vector<FiniteElement*>::const_iterator vIterFEs(vecFEs_.begin());
  	vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVS;
  	const vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVSEnd(fvs_.end());
  	for(vIterFVS = fvs_.begin(); vIterFVS != vIterFVSEnd; vIterFVS++, vIterFEs++)
  	{
  	 const CSMP_FEM_TYPE elType((*vIterFEs)->ElementType());
     if ( verbose_ ) cout << "Type: " << parseFiniteElementType( elType ) << endl;
		 
		 //check number of sectors , should be = to number of nodes
     const bool bCheckNodeNr(vIterFVS->Sectors() == (*vIterFEs)->Nodes());
     _test(bCheckNodeNr);
     
  	 //check number of integration points per facet, currently one per facet
  	 _test(vIterFVS->IntegrationPointsPerFacet() == 1);
  	 //check number of integration points per sector, currently one per sector
  	 _test(vIterFVS->IntegrationPointsPerSector() == 1);
  	 
  	 //check total number of facets, should be = to number of segments 
  	 const size_t iNrOfFacets(vIterFVS->Facets()); 
  	 
#ifdef PYRAMID_TRIANGULAR_FACETS
  	   if(vIterFVS->ParentElement() == "ISOPARAMETRIC_LINEAR_PYRAMID")
  	     _test(iNrOfFacets == (*vIterFEs)->Segments()+4);
  	   else
#endif
  	     _test(iNrOfFacets == (*vIterFEs)->Segments());
  	    
  	 //for each facet: check edge nodes per facet, should be = to edge nodes per segment
#ifdef PYRAMID_TRIANGULAR_FACETS
  	 //this is not valid for the pyramid finite element with 12 triangular facets 
  	 if(vIterFVS->ParentElement() != "ISOPARAMETRIC_LINEAR_PYRAMID")
#endif
     {
  	     size_t iInsideNode(0U), iOutsideNode(0U);
	     for(size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++)
  	     {	
	   	    //get facet nodes from stencil
	  	    vIterFVS->FacetEdgeNodes(iFacet, iInsideNode, iOutsideNode);
	  	    	
	  	    //get facet nodes from finite element
  		    (*vIterFEs)->NodesOfSegment(iFacet, vecNodeIdsFEs);
	  	    assert(vecNodeIdsFEs.size() >= 2);
	  	    	
			_test(vecNodeIdsFEs[0] == iInsideNode);
	  	    _test(vecNodeIdsFEs[1] == iOutsideNode);
  	     }
      }

  	  // tests which should work for all element types
	  for(size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++)
  	  {	
        size_t expectedFacetPoints = ~(size_t)0;
        switch (vIterFVS->FacetType(iFacet))
        {
            case POINT_FACET:
                expectedFacetPoints = 1;
                break;
            case UNIT_LINEAR_FACET:
                expectedFacetPoints = 2;
                break;
            case TRIANGULAR_FACET:
                expectedFacetPoints = 3;
                break;
            case QUADRILATERAL_FACET:
                expectedFacetPoints = 4;
                break;
        }
        _test(vIterFVS->FacetPoints(iFacet) == expectedFacetPoints);
      }
   }	
 }
 
 
 
 
 
/**  void FiniteVolumeStencil_Test::sectorFacetConnectivityTest()
 
@section Description:

This test checks for each element type, that for each sector the corresponding surrounding facets are correct. 
It checks that the facets surrounding a sector correspond to the segments surrounding the node 
@section arguments Input Arguments 


@section application Application

The idea is to test the following member functions:

 	size_t  FacetsPerSector() const;
    void       Facet(s)SurroundingSector(const size_t sector, vector<size_t>& surf ) const;
    fT         SectorIntegrationPoint( size_t sector, size_t ip, size_t r_or_s_or_t ) const;
    void       SectorIntegrationPoint( size_t sector, size_t ip, vector<fT>& rst ) const;
 
*/
void FiniteVolumeStencil_Test::sectorFacetConnectivityTest() // facets that delimit a certain sector for all the element types
 {
  if ( verbose_ ) cout << "TESTING: sectorFacetConnectivityTest()" << endl;
 	//temp vector to get facets
 	vector<size_t> vecFacets;
 	set<size_t> setFacets;
 	
 	//neighbors of node relationship container
 	map<size_t/*node*/, set<size_t>/*neighbor nodes*/> mapNodeToSegmentNeighbors;
 	vector<size_t> vecNodesOfSegment; 
 	
  //for each finite element type vs finite volume stencil type
  vector<FiniteElement*>::const_iterator vIterFEs(vecFEs_.begin());
	vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVS;
	const vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVSEnd(fvs_.end());
	for(vIterFVS = fvs_.begin(); vIterFVS != vIterFVSEnd; vIterFVS++, vIterFEs++)
	{
	  if(  vIterFVS->ParentElement() == "ISOPARAMETRIC_LINEAR_PYRAMID" ) continue;
  	if(  vIterFVS->ParentElement() == "ISOPARAMETRIC_LINEAR_TETRAHEDRON" ) continue;
  	
  	const CSMP_FEM_TYPE elType((*vIterFEs)->ElementType());
    if ( verbose_ ) cout << "Type: " << parseFiniteElementType( elType ) << endl;
		
		//build neighbors of node relationship - initialize map
 		mapNodeToSegmentNeighbors.clear();
 		const size_t iNumberOfNodes((*vIterFEs)->Nodes());
 		for(size_t iNode = 0U; iNode < iNumberOfNodes; iNode++)
 		{
 			mapNodeToSegmentNeighbors.insert(make_pair(iNode, set<size_t>()));
 		}
 	
 		//fill map
 		const size_t iNumberOfSegments((*vIterFEs)->Segments());
 		for(size_t iSegment = 0U; iSegment < iNumberOfSegments; iSegment++)
 		{
 			vecNodesOfSegment.clear();
 			(*vIterFEs)->NodesOfSegment( iSegment, vecNodesOfSegment );
 			assert(vecNodesOfSegment.size() >= 2);
 			mapNodeToSegmentNeighbors[ vecNodesOfSegment[0] ].insert(iSegment);
 			mapNodeToSegmentNeighbors[ vecNodesOfSegment[1] ].insert(iSegment);
 		} 		
  	  	
  	//get number of sectors
  	const size_t iNrOfSectors(vIterFVS->Sectors()); 
  	    
  	//test size_t  FacetsPerSector() const;
  	//check if the number of facets per sector matches the reconstructed map
  	size_t iMaxNodes(0U);
 		for(size_t iSector = 0U; iSector < iNrOfSectors; iSector++)
  	{
  		if(mapNodeToSegmentNeighbors[iSector].size() > iMaxNodes)
  			iMaxNodes = mapNodeToSegmentNeighbors[iSector].size(); 
   	}
   		
    //TRACE cout << "iMaxNodes: " << iMaxNodes << " FacetsPerSector: " << vIterFVS->FacetsPerSector() << endl;
  	_test(iMaxNodes == vIterFVS->FacetsPerSector(0));
   		
   		//get reference coordinates
		DenseMatrix<DM_MIN> matCoords;
		(*vIterFEs)->ReferenceCoordinates(matCoords);
		const size_t iDim(matCoords.Cols());
  		
  	assert(mapNodeToSegmentNeighbors.size() == iNrOfSectors);
  		
   	//test integration point location for sectors
   	for(size_t iSector = 0U; iSector < iNrOfSectors; iSector++)
  	{
  		vecFacets.clear();
  		
  		copy(mapNodeToSegmentNeighbors[iSector].begin(), mapNodeToSegmentNeighbors[iSector].end(), back_insert_iterator< vector<size_t> >(vecFacets));
  	
  	  //1. Calculate sector centroid
  	  Point<3> vecCentroidGenerated;
  		//load the actual coordinates of the nodes which belong to this face
			vector< Point<3> > vecPointsOfSector;
			//insert node corresponding to this sector 
			Point<3> vecPt;
			vecPt=0.;
			for(size_t i = 0U; i < iDim; i++)
				vecPt[i] = matCoords(/*this node*/iSector,i); 
			
			vecPointsOfSector.push_back(vecPt);
				
			vector<size_t>::iterator vIterSegment;
			vector<Point<3> > vecPointsOfFacets;
			for(vIterSegment = vecFacets.begin(); vIterSegment != vecFacets.end(); vIterSegment++ )
			{
				vector< Point<3> > vecOfPointsOfTheFacet;
				vecOfPointsOfTheFacet = GetPointsOfFacet(*vIterSegment, elType, *vIterFVS, **vIterFEs);
		    
		    for(vector< Point<3> >::const_iterator vIterToInsert = vecOfPointsOfTheFacet.begin(); vIterToInsert != vecOfPointsOfTheFacet.end(); vIterToInsert++)
  		  {
  		   bool bInsert(true);
  		   for(vector< Point<3> >::const_iterator vIterPts = vecPointsOfSector.begin(); vIterPts != vecPointsOfSector.end(); vIterPts++)
  		    if(*vIterPts == *vIterToInsert) { bInsert=false; break; }
  		   if(bInsert)
  		    vecPointsOfSector.push_back(*vIterToInsert);
  		  }
  		}
  		
  		//2. Get centroid from sector
  		if(vecPointsOfSector.size() == 4)
  			areaCenterOfMass(vecPointsOfSector, vecCentroidGenerated);
  		else
  			vertexCenterOfMass(vecPointsOfSector, vecCentroidGenerated);
  		
      //TRACE
      if ( verbose_ ) {
          cout << "PTS(facet)::: sector:" << iSector << ":" <<endl;
          for(size_t i = 0U; i < vecPointsOfSector.size(); i++)
          {
            cout << "P"<<i<<": ";
            vecPointsOfSector[i].Out();
          }	
          cout << "Sector " << iSector << " centroid generated:";
          vecCentroidGenerated.Out();
  		  }
  		//get sector integration point location
  		Point<3> vecCentroidOfSector_I1, vecCentroidOfSector_I2;
  		vecCentroidOfSector_I1 = 0.;
  		vecCentroidOfSector_I2 = 0.;
  		
  		//get centroid through first interface
  		for ( size_t j = 0U; j< vIterFVS->IntegrationPointsPerSector(); j++ ) 
	     	vecCentroidOfSector_I1 += vIterFVS->SectorIntegrationPoint(iSector, j); 
			
  		//get centroid through second interface
  		for(size_t iD = 0U; iD < iDim; iD++)
			 for ( size_t j = 0U; j< vIterFVS->IntegrationPointsPerSector(); j++ ) 
	    	vecCentroidOfSector_I2[iD] += vIterFVS->SectorIntegrationPoint(iSector, j, iD); 
			
      if ( verbose_ ) {
          cout << "Sector " << iSector << " centroid original (I1):";
          vecCentroidOfSector_I1.Out();
          cout << "Sector " << iSector << " centroid original (I2):";
          vecCentroidOfSector_I2.Out();
        }
  		
  		//3. Compare calculated vs accessed sector centroids.
  		
  		//test first interface
  		for(size_t iD = 0U; iD < iDim; iD++)
			 if(!isnan(vecCentroidOfSector_I1[iD]))
			  _equal(vecCentroidGenerated[iD], vecCentroidOfSector_I1[iD], 1.e-7);
			//test second interface
			for(size_t iD = 0U; iD < iDim; iD++)
			 if(!isnan(vecCentroidOfSector_I2[iD]))
			    _equal(vecCentroidGenerated[iD], vecCentroidOfSector_I2[iD], 1.e-7);
			
  		if(mapNodeToSegmentNeighbors[iSector].size() > iMaxNodes)
  	  	iMaxNodes = mapNodeToSegmentNeighbors[iSector].size(); 
   	}
   			
  	//clear facets
  	setFacets.clear();
  	//for each sector: check that the facets surrounding a sector correspond to the segments surrounding the node 
  	for(size_t iSector = 0U; iSector < iNrOfSectors; iSector++)
  	{
      if ( verbose_ ) cout << "SECTOR "<< iSector <<": ";
  		// the facets surrounding a sector is to the segments surrounding the node
  		// as the sector is to the node  
  		setFacets.clear();
  		for(size_t iFacet = 0; iFacet < vIterFVS->FacetsPerSector(iSector); iFacet++)
  		{
  		  setFacets.insert(vIterFVS->FacetSurroundingSector(iSector, iFacet));
  		}
  		//compare sets
  		_test(setFacets == mapNodeToSegmentNeighbors[iSector]);
  		
  	}	
  } 
}
 



/**   void FiniteVolumeStencil_Test::orientationAndLengthOfNormalsTest()

Description: 

This test checks for each element type, that the corresponding normal vector is in fact normalized.
@section arguments Input Arguments 

@section application Application
The idea is to test the following member functions:

	void   UnitParametricNormalTo( size_t facet, vector<fT>& nrml ) const;
    fT   UnitParametricNormalComponent( size_t facet, size_t x_or_y_or_z ) const;
 
*/
void FiniteVolumeStencil_Test::orientationAndLengthOfNormalsTest()
 {
  if ( verbose_ ) cout << "TESTING: orientationAndLengthOfNormalsTest()" << endl;

	//for each element type, for each given facet, check if the normal is correct
 	
 	Point<3> vecNormal;
 	
 	vector<FiniteElement*>::const_iterator vIterFEs(vecFEs_.begin());
  vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVS;
  const vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVSEnd(fvs_.end());
  for(vIterFVS = fvs_.begin(); vIterFVS != vIterFVSEnd; vIterFVS++, vIterFEs++)
  {
  		//get finite element type
  		const CSMP_FEM_TYPE elType((*vIterFEs)->ElementType());

  		//print type
      if ( verbose_ ) {
          cout << "\n******************************************************************";
          cout << "\nType: " << parseFiniteElementType(elType);
          cout << "\nNumber of Facets: " << vIterFVS->Facets();
        }
  		
  		//step through facets
  		const size_t iNrOfFacets(vIterFVS->Facets());
  		for(size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++)
  		{
          if ( verbose_ ) {
              //get facet
              cout << "\n******************************************************************";
              cout << "\nFACET " << iFacet << ": ";
            }
	  	    vecNormal = vIterFVS->UnitParametricNormalTo( iFacet );	  	    
          
          //TRACE
          if ( verbose_ ) cout << "\nNormal: " << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2];
			   
			    //1. Test Length
			    //length of normal, should be 1		
			    const double64 fLength(vecNormal.Length());
			    //test if normals are unit vectors
					cout << setprecision(15) << "\nLength of Normal: " << fLength << endl;
			    _equal( fLength, 1., 1.e-7 );
				
			    //2. Generate normal on-the-fly
			    vector< Point< 3> > vecOfPointsOfTheFacet;
			    
			    size_t iNrOfFacetPoints(0);
			    if     ( (*vIterFEs)->IsLineElement() )    iNrOfFacetPoints = 1;
			    else if( (*vIterFEs)->IsSurfaceElement() ) iNrOfFacetPoints = 2;
			    else if( (*vIterFEs)->IsVolumeElement() )  iNrOfFacetPoints = 4;
			    
			  #ifdef PYRAMID_TRIANGULAR_FACETS
			    if(vIterFVS->ParentElement() == "ISOPARAMETRIC_LINEAR_PYRAMID" && iFacet > 3 ) iNrOfFacetPoints = 3;
			  #endif
			  
			    for(size_t iPoint = 0; iPoint < iNrOfFacetPoints; iPoint++)
			    {
			      vecOfPointsOfTheFacet.push_back(vIterFVS->FacetPoint(iFacet,iPoint));
			      vecOfPointsOfTheFacet[iPoint].Out();
			    }
			    
			    //generate normal
			    //calculate normal at the integration point	
			    Point<3U> vecGeneratedNormal, vecTemp;
			    vecTemp = vIterFVS->FacetIntegrationPoint(iFacet, 0U);
			    
					if ( verbose_ ) cout << "\nNumber of facet points2: " << vecOfPointsOfTheFacet.size();
			    normalOfPolygon(vecOfPointsOfTheFacet, vecTemp, vecGeneratedNormal);
			
    			//print normals
          if ( verbose_ ) cout << "\nGenerated: " << vecGeneratedNormal[0] << ", " << vecGeneratedNormal[1] << ", " << vecGeneratedNormal[2] << endl;
	  	    
			    //3. Test the orientation of the normal in each direction
			    bool bEquivalent1(true), bEquivalent2(true);
			
			    for(size_t iVal = 0; iVal < 3; iVal++)
				    bEquivalent1 &= ( fabs(vecNormal[iVal] - vecGeneratedNormal[iVal]) < 1.e-7 );
				
			    for(size_t iVal = 0; iVal < 3; iVal++)
				    bEquivalent2 &= ( fabs(vecNormal[iVal] + vecGeneratedNormal[iVal]) < 1.e-7 );
				
		  	  _test(bEquivalent1 || bEquivalent2);
		  }
  	}
 }
    
    struct FE_Transformation {
        FE_Transformation(csmp::FiniteElement* element)
        : element_(element)
        {
        }

        FiniteElement* element_;
    };
    
    /**   void FiniteVolumeStencil_Test::normalTransformationTest()
     
     Description:
     
     This test checks for each facet normal transformation, and ensures that
     it matches the computed facet normal.

     @section arguments Input Arguments
     
     @section application Application
     The idea is to test the following member functions:
     
     wp   FacetNormalTransformationNodeWeights( size_t facet, size_t node ) const;
     
     */
    void FiniteVolumeStencil_Test::normalTransformationTest()
    {
        if ( verbose_ ) cout << "TESTING: normalTransformationTest()" << endl;
        
        //for each element type, for each given facet, check if the normal is correct
        
        vector<FiniteElement*>::const_iterator vIterFEs(vecFEs_.begin());
        vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVS;
        const vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVSEnd(fvs_.end());
        DenseMatrix<DM_MIN> matCoords;

        for(vIterFVS = fvs_.begin(); vIterFVS != vIterFVSEnd; vIterFVS++, vIterFEs++)
        {
            //get finite element type
            const CSMP_FEM_TYPE elType((*vIterFEs)->ElementType());
            
#ifndef PYRAMID_TRIANGULAR_FACETS
            if((*vIterFEs)->ElementType() == ISOPARAMETRIC_LINEAR_PYRAMID ) {
                std::cerr<< "XXX WARNING: ISOPARAMETRIC_LINEAR_PYRAMID with quadrilateral apex facets NYI\n";
            }
#endif

            // Get coordinates
            (*vIterFEs)->ReferenceCoordinates( matCoords );
            
            //print type
            if ( verbose_ ) {
                cout << "\n******************************************************************";
                cout << "\nType: " << parseFiniteElementType(elType);
                cout << "\nNumber of Facets: " << vIterFVS->Facets();
            }
            
            //step through facets
            const size_t iNrOfFacets(vIterFVS->Facets());
            const size_t iNrOfNodes((*vIterFEs)->Nodes());

            for(size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++)
            {
                if ( verbose_ ) {
                    //get facet
                    cout << "\n******************************************************************";
                    cout << "\nFACET " << iFacet << ": ";
                }
                const Point<3> parametricNormal(vIterFVS->UnitParametricNormalTo( iFacet ));

                //TRACE
                if ( verbose_ ) cout << "\nParametric normal: " << parametricNormal[0] << ", " << parametricNormal[1] << ", " << parametricNormal[2];
                
                Point<3> computedNormal;

                //1. Compute normal by transformation
                if ((*vIterFEs)->IsVolumeElement()) {
                    Point<3> v0(0,0,0);
                    Point<3> v1(0,0,0);

                    for (size_t iNode = 0U; iNode < iNrOfNodes; ++iNode) {
                        const Point<3u> n(matCoords(iNode,0),matCoords(iNode,1),matCoords(iNode,2));
                        auto weights = vIterFVS->FacetNormalTransformationNodeWeights(iFacet, iNode);
                        v0 += weights.first * n;
                        v1 += weights.second * n;
                    }
                    computedNormal = crossProduct(v1,v0);
                }
                else if ((*vIterFEs)->IsSurfaceElement()) {
                    Point<3u> tangent(0,0,0);
                    Point<3u> bitangent(0,0,0);

                    for (size_t iNode = 0U; iNode < iNrOfNodes; ++iNode) {
                        const Point<3u> n(matCoords(iNode,0),matCoords(iNode,1),matCoords(iNode,2));
                        auto weights = vIterFVS->FacetNormalTransformationNodeWeights(iFacet, iNode);
                        tangent += weights.first * n;
                        bitangent += weights.second * n;
                    }
                    double64 length = exteriorProductLength(tangent, bitangent);
                    tangent.NormalizeLengthTo(1.0);
                    computedNormal = bitangent - dotProduct(tangent,bitangent) * tangent;
                    computedNormal.NormalizeLengthTo(length);
                }
                else if ((*vIterFEs)->IsLineElement()) {
                    Point<3u> v0(0,0,0);
                    for (size_t iNode = 0U; iNode < iNrOfNodes; ++iNode) {
                        const Point<3u> n(matCoords(iNode,0),matCoords(iNode,1),matCoords(iNode,2));
                        auto weights = vIterFVS->FacetNormalTransformationNodeWeights(iFacet, iNode);
                        v0 += weights.first * n;
                    }
                    computedNormal = v0;
                }


                double64 computedNormalLength(computedNormal.Length());
                computedNormal.NormalizeLengthTo(1.0);
                if ( verbose_ ) cout << "\nComputed normal: " << computedNormal[0] << ", " << computedNormal[1] << ", " << computedNormal[2];
                if ( verbose_ ) cout << "\nComputed normal length: " << computedNormalLength;

                //2. Find normal and area from the stencil
                
                const Point<3U> stencilNormal = vIterFVS->UnitParametricNormalTo(iFacet);
                const double64 stencilArea = vIterFVS->FacetIntegrationWeight(iFacet, 0);
                
                //print normals
                if ( verbose_ ) cout << "\nFrom stencil: " << stencilNormal[0] << ", " << stencilNormal[1] << ", " << stencilNormal[2] << endl;
                
                //3. Test the orientation of the normal and calcualted area
                bool bEquivalent(true);
                
                for(size_t iVal = 0; iVal < 3; iVal++)
                    bEquivalent &= ( fabs(computedNormal[iVal] - stencilNormal[iVal]) < 1.e-7 );
                
                _test(bEquivalent);
                std::cerr << "\nElement type " << parseFiniteElementType(elType) << " facet " << iFacet << '\n';
                _equal(stencilArea, computedNormalLength, 1.0e-8);
            }
        }
    }
    
    
    



/**  void FiniteVolumeStencil_Test::weightsAndFacetIntegrationPointsTest()
 
@section Description

Check if the area of the facet corresponds to the weight of the integration point. These should match.
Check if the centroid of the facet corresponds to the integration point.

@section application Application

The idea is to test the following member functions:
	fT         FacetIntegrationWeight( size_t facet, size_t ip ) const; // = facet area
    fT         FacetProjectionWeight( size_t facet, size_t ip ) const; // equal to scaling of normal length
	fT         FacetIntegrationPoint( size_t facet, size_t ip, size_t r_or_s_or_t ) const;
    void       FacetIntegrationPoint( size_t facet, size_t ip, vector<fT>& rst ) const;
 
*/
void FiniteVolumeStencil_Test::weightsAndFacetIntegrationPointsTest()
{
  if ( verbose_ ) cout << "TESTING: weightsAndFacetIntegrationPointsTest()" << endl;
  
  //for each finite element type
 	vector<FiniteElement*>::const_iterator vIterFEs(vecFEs_.begin());
  vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVS;
  const vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVSEnd(fvs_.end());
  for(vIterFVS = fvs_.begin(); vIterFVS != vIterFVSEnd; vIterFVS++, vIterFEs++)
  {
  	const CSMP_FEM_TYPE elType((*vIterFEs)->ElementType());
		
		//print type
		static int nr(0);
		if ( verbose_ ) cout << "\nType: " << parseFiniteElementType( elType ) << ": " << nr++ << endl;
  	
		//get reference coordinates
		DenseMatrix<DM_MIN> matCoords;
		(*vIterFEs)->ReferenceCoordinates(matCoords);
		
		//the dimension is given by the number of columns
		const size_t iDim(matCoords.Cols());
		
		double64 fFacetWeight(0.), fFacetArea(0.);		
  		
		vector< Point<3U> > vecOfPointsOfTheFacet;
		//for each facet, get the area, and compare it to the weight of the integration point
    //step through each facet to get weight
		const size_t iNrOfFacets(vIterFVS->Facets());
    for(size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++)
  	{
  		if(  vIterFVS->ParentElement() == "ISOPARAMETRIC_LINEAR_PYRAMID" ) continue;
  		
      if ( verbose_ ) cout << "--->FACET: " << iFacet << endl;
  		fFacetWeight = 0.; fFacetArea = 0.;
			vecOfPointsOfTheFacet = GetPointsOfFacet(iFacet, elType, *vIterFVS, **vIterFEs);
		  
		  //TRACE
      if ( verbose_ ) {
          cout << "PTS::: facet:" << iFacet << ":" <<endl;
          for(size_t i = 0U; i < vecOfPointsOfTheFacet.size(); i++)
          {
            cout << "P"<<i<<": ";
            vecOfPointsOfTheFacet[i].Out();
          }
       }
  			
			//get facet weight
			//sum all facet weights (there should be one integration point per facet)
			for ( size_t j = 0U; j< vIterFVS->IntegrationPointsPerFacet(); j++ ) 
			 fFacetWeight += vIterFVS->FacetIntegrationWeight( iFacet, j );
		
      //TEST WEIGHT
      //get the area/distance
       size_t iNrOfFacetPoints(0);
       
       if     ( (*vIterFEs)->IsLineElement() )    iNrOfFacetPoints = 1;
			 else if( (*vIterFEs)->IsSurfaceElement() ) iNrOfFacetPoints = 2;
			 else if( (*vIterFEs)->IsVolumeElement() )  iNrOfFacetPoints = 4;
			 
			   
			 if(iNrOfFacetPoints == 1)
        fFacetArea = 1.;
       else if(iNrOfFacetPoints == 2) //compute only distance
        fFacetArea = vecOfPointsOfTheFacet[0].DistanceTo(vecOfPointsOfTheFacet[1]);
       else if (iNrOfFacetPoints == 4)
        _test(areaOfPolygon(vecOfPointsOfTheFacet, iDim, fFacetArea));
      
      if ( verbose_ ) cout << "\nFacet Area" << fFacetArea << " vs " << fFacetWeight;
      _equal( fFacetArea, fFacetWeight, 1.e-7 );
      //END TEST WEIGHT
      
      //TRACE
      if ( verbose_ ) cout << "\nFacetArea (computed): " << fFacetArea << " vs Facet Weight (stored):" << fFacetWeight << endl;
        	
      //get facet integration point -first interface
  		Point<3> vecCentroidFacet_I1, vecCentroidFacet_I2;
  		vecCentroidFacet_I1 = 0.;
  		vecCentroidFacet_I2 = 0.;
  		for ( size_t j = 0U; j< vIterFVS->IntegrationPointsPerFacet(); j++ ) 
      	vecCentroidFacet_I1 += vIterFVS->FacetIntegrationPoint(iFacet, j); 
				
			//get facet integration point -second interface
  		vecCentroidFacet_I2 =	vIterFVS->FacetIntegrationPoint(iFacet, 0U); 
			
			//generate centroid
			Point<3> vecCentroidGeneratedPointsFacets;
			
			if(vecOfPointsOfTheFacet.size() == 4)
        areaCenterOfMass(vecOfPointsOfTheFacet, vecCentroidGeneratedPointsFacets);
			else
			  vertexCenterOfMass(vecOfPointsOfTheFacet, vecCentroidGeneratedPointsFacets);
			
			//TRACE
      if ( verbose_ ) {
          cout << "\nOriginal centroid (I1):";
          vecCentroidFacet_I1.Out();
          cout << "\nOriginal centroid (I2):";
          vecCentroidFacet_I2.Out();
        }
			
      //test first interface
      for(size_t iD = 0U; iD < iDim; iD++)
			  if(!isnan(vecCentroidFacet_I1[iD]))
			   _equal(vecCentroidFacet_I1[iD], vecCentroidGeneratedPointsFacets[iD],1.e-7);
			//test second interface
      for(size_t iD = 0U; iD < iDim; iD++)
			  if(!isnan(vecCentroidFacet_I2[iD]))
			   _equal(vecCentroidFacet_I2[iD], vecCentroidGeneratedPointsFacets[iD],1.e-7);
		
			//END TEST LOCATION
      }
  	}
 }
 
 
 
 

/**  void FiniteVolumeStencil_Test::weightsOfSectorIntegrationPointsTest()

Verifies if the sum of weights corresponds to the volume of the finite element.

@section application Application

The idea is to test the following member functions:

	fT  SectorIntegrationWeight( size_t sector, size_t ip ) const;
 
*/
void FiniteVolumeStencil_Test::weightsOfSectorIntegrationPointsTest() // argument list could deal with multiple integration points
 { 
  if ( verbose_ ) cout << "TESTING: weightsOfSectorIntegrationPointsTest()" << endl;
    //calculate total weight of sectors, which should correspond to the total volume of the finite element
 	double64 fSumOfWeights(0.);
 	
 	//map types to total area
 	map<CSMP_FEM_TYPE, double64> mapFEMTypeToArea;
 	mapFEMTypeToArea.insert(make_pair(ISOPARAMETRIC_LINEAR_BAR, 2.));
 	mapFEMTypeToArea.insert(make_pair(ISOPARAMETRIC_LINEAR_TRIANGLE, 1./2.));
 	mapFEMTypeToArea.insert(make_pair(ISOPARAMETRIC_LINEAR_QUADRILATERAL, 4.));
 	mapFEMTypeToArea.insert(make_pair(ISOPARAMETRIC_LINEAR_TETRAHEDRON, 1./6.));
 	mapFEMTypeToArea.insert(make_pair(ISOPARAMETRIC_LINEAR_PYRAMID, 4./3.));
 	mapFEMTypeToArea.insert(make_pair(ISOPARAMETRIC_LINEAR_HEXAHEDRON, 8.));
 	mapFEMTypeToArea.insert(make_pair(ISOPARAMETRIC_LINEAR_PRISM, 1.));
 	
 	//for each element type
 	vector<FiniteElement*>::const_iterator vIterFEs(vecFEs_.begin());
  	vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVS;
  	const vector<FiniteVolumeStencil< 3> >::const_iterator vIterFVSEnd(fvs_.end());
  	for(vIterFVS = fvs_.begin(); vIterFVS != vIterFVSEnd; vIterFVS++, vIterFEs++)
  	{
  		//print type
      if ( verbose_ )
        cout << "Type: " << parseFiniteElementType( (*vIterFEs)->ElementType() ) << endl;
  		
 		  fSumOfWeights = 0.;
  		
  		//step through each sector to get weight
		  const size_t iNrOfSectors(vIterFVS->Sectors());
  		for(size_t iSector = 0U; iSector < iNrOfSectors; iSector++)
  			for ( size_t j = 0U; j< vIterFVS->IntegrationPointsPerSector(); j++ ) 
        		fSumOfWeights += vIterFVS->SectorIntegrationWeight( iSector, j );
        
      if ( verbose_ )
        cout << "SUM OF WEIGHTS:" << setprecision(15) << fSumOfWeights << "vs." << mapFEMTypeToArea[(*vIterFEs)->ElementType()] << endl;
		
		  _equal( fSumOfWeights, mapFEMTypeToArea[(*vIterFEs)->ElementType()], 1.e-9);
  	}
 }
 
 
 
 
template<size_t dim>
vector< Point<dim> > FiniteVolumeStencil_Test::GetPointsOfFacet(const size_t& iFacet, const CSMP_FEM_TYPE& elType, 
                                                                const FiniteVolumeStencil<dim>& fvs, const FiniteElement& fe)
{
	DenseMatrix<DM_MIN> matCoords;
	fe.ReferenceCoordinates(matCoords);
		
	//the dimension is given by the number of columns
	vector< Point<dim> > vecOfPointsOfTheFacet;
	
	vector<size_t>::iterator vIterNode;
	
	//build neighbor face VertexCenterOfMasss for each edge
	multimap<size_t /*iEdge*/, Point<dim> /*vecVertexCenterOfMass*/> mmapEdgeToNeighborVertexCenterOfMasss;
	typename multimap<size_t /*iEdge*/, Point<dim> >::const_iterator iterMmap;		
	{
		//no need to build map for the following types:
		if( elType != ISOPARAMETRIC_LINEAR_BAR && 
		    elType != ISOPARAMETRIC_LINEAR_TRIANGLE &&
			elType != ISOPARAMETRIC_LINEAR_QUADRILATERAL)
		{
			//build map
			
			//step through edges (edge = segment)
			const size_t iNrOfEdges(fe.Segments());
			for(size_t iEdge = 0U; iEdge < iNrOfEdges; iEdge++)
			{
				//get nodes of segment
				vector<size_t> vecNodesOfEdge;
				fe.NodesOfSegment( iEdge, vecNodesOfEdge );
				
				//step through faces
				const size_t iNrOfFaces(fe.Faces());
				for(size_t iFace = 0U; iFace < iNrOfFaces; iFace++)
				{
					//get nodes of face
			    	vector<size_t> vecNodesOfFace;
					fe.NodesOfFace( iFace, vecNodesOfFace );
				
					//is the face adjacent to the edge?
					//check if each node of the edge is contained in the face
					bool bAdjacent = true;
					for(vIterNode = vecNodesOfEdge.begin(); vIterNode != vecNodesOfEdge.end(); vIterNode++ )
					{
						if(find(vecNodesOfFace.begin(), vecNodesOfFace.end(), *vIterNode) == vecNodesOfFace.end())
						{
							//we did NOT find it!
							bAdjacent = false;
						}
					}
					
					if(bAdjacent)
					{
						//load the actual coordinates of the nodes which belong to this face
						vector< Point<dim> > vecPointsOfFace;
						for(vIterNode = vecNodesOfFace.begin(); vIterNode != vecNodesOfFace.end(); vIterNode++ )
						{
							Point<dim> vecPt;
							for(size_t i = 0U; i < dim; i++)
							{
								vecPt[i] = matCoords(*vIterNode,i); 
							}
							vecPointsOfFace.push_back(vecPt);
						}
						
						//add to map
						Point<dim> vecTemp;
						vertexCenterOfMass(vecPointsOfFace, vecTemp);
 						mmapEdgeToNeighborVertexCenterOfMasss.insert(make_pair(iEdge, vecTemp));
					}	
				} //end for - faces	 
			} //end for - edges 				
		}//end if
	} 	
	pair< typename multimap<size_t, Point<dim> >::const_iterator, typename multimap<size_t, Point<dim> >::const_iterator> iterPair;
	
	//the points which belong to a facet are: 
	//vecCenterOfEdge := center of the edge which identifies the facet, 
	//vecCenterOfContiguous_Faces := center of the contiguous faces to the edge -> stored in the map: mmapEdgeToNeighborVertexCenterOfMasss
	//vecVertexCenterOfMass := center of the polyhedron
	Point<dim> vecVertexCenterOfMass, vecCenterOfEdge;

	//get element's VertexCenterOfMass
	vecVertexCenterOfMass = fvs.Barycenter();
	
	if(elType == ISOPARAMETRIC_LINEAR_BAR)
	{
		vecOfPointsOfTheFacet.push_back(vecVertexCenterOfMass);
		return vecOfPointsOfTheFacet;
	}
	
	//fill it with zeros
	vecCenterOfEdge = 0.;
	 
	//calculate facet area
	//get vecCenterOfEdge
	size_t iInsideNode(0U), iOutsideNode(0U);
	fvs.FacetEdgeNodes(iFacet, iInsideNode, iOutsideNode);
	
	for(size_t i = 0U; i < matCoords.Cols(); i++)
	{
		//dimension of the vector corresponds to the number of columns of the reference matrix
		vecCenterOfEdge[i] += matCoords(iInsideNode, i);
		vecCenterOfEdge[i] += matCoords(iOutsideNode, i);
		vecCenterOfEdge[i] /= 2.;
	}
	
	vecOfPointsOfTheFacet.push_back(vecCenterOfEdge);
	
	//get center of adjacent faces
	//loop through adjacent faces, add VertexCenterOfMasss alternatively back and front
	bool bAlternate(true);
	iterPair = mmapEdgeToNeighborVertexCenterOfMasss.equal_range(iFacet/*edge number should correspond to facet number*/);
	for(iterMmap = iterPair.first; iterMmap != iterPair.second; iterMmap++)
	{
		if(bAlternate)
			vecOfPointsOfTheFacet.push_back(iterMmap->second);
		else
			vecOfPointsOfTheFacet.insert(vecOfPointsOfTheFacet.begin(), iterMmap->second);
		bAlternate = !bAlternate;
	}
	
	//add VertexCenterOfMass
	vecOfPointsOfTheFacet.push_back(vecVertexCenterOfMass);
	
	return vecOfPointsOfTheFacet;       	
}





/**  template<size_t dim>
      void IsoparametricQuadraticElements_Test<dim>::shapeFunctionsTest()

Tests if the shape functions sum to one at the integration points.

@section application Application
Runs the test for all Isoparametric Quadratic elements, for any dimensions.
 
*/
void FiniteVolumeStencil_Test::shapeFunctionsTest()
{
  //shape functions must sum to 1
  if ( verbose_ ) cout << "\nTESTING: shapeFunctionsTest()" << endl;
 
  for(vector<FiniteElement*>::const_iterator vIterFEs = vecFEs_.begin(); vIterFEs != vecFEs_.end(); vIterFEs++)
  {
  	//print type
    if ( verbose_ ) cout << "\nType: " << parseFiniteElementType( (*vIterFEs)->ElementType() ) << endl;
  	
  	std::vector<double64> sf;
  	double64 sum(0.);
  	    
  	//check that they sum up at the integration points
  	for(size_t i = 0; i < (*vIterFEs)->IntegrationPoints(); i++)
  	{
  	  (*vIterFEs)->N_AtIntegrationPoint( i, sf );
  	  sum = 0.;
      if ( verbose_ ) cout << "\n";
  	  for(size_t j = 0; j < sf.size(); j++)
  	  {  
  	    sum += sf[j];
        //TRACE cout << "(" << shapeFuncsAtIP[j] <<")";
  	  }
      //TRACE cout << "\nShape functions should be 1. sum: " << sum;
  	  _equal( sum, 1., 1.e-15 );
  	}	
  	//check that they sum up at the barycenters
  	try
  	{
  	  (*vIterFEs)->N_AtBaryCenter( sf );
  	  sum = 0.;
  	  for(size_t j = 0; j < sf.size(); j++)  
  	   sum += sf[j];
      if ( verbose_ ) cout << "\nShape functions should be 1. sum: " << sum;
  	  _equal( sum, 1., 1.e-14 ); 
  	}
  	catch( invalid_argument e ) //its not defined
    {
       if ( verbose_ ) cout << "Exception: " << e.what(); //do not fail test
    } 
	}	
}





/**  template<size_t dim>
void IsoparametricQuadraticElements_Test<dim>::shapeFunctionDerivativesTest()
 
Tests if the shape functions sum to zero at the integration points.

@section application Application
Runs the test for all Isoparametric Quadratic elements, for any dimensions.
 
*/
void FiniteVolumeStencil_Test::shapeFunctionDerivativesTest()
{
  //shape functions must sum to 1
  if ( verbose_ ) cout << "\nTESTING: shapeFunctionDerivativesTest()" << endl;
  for(vector<FiniteElement*>::const_iterator vIterFEs = vecFEs_.begin(); vIterFEs != vecFEs_.end(); vIterFEs++)
  {
  	//print type
    if ( verbose_ ) cout << "\nType: " << parseFiniteElementType( (*vIterFEs)->ElementType() ) << endl;
  	
  	double64 sum(0.);
  	DenseMatrix<DM_MIN> sf;
  	  
  	//check that they sum up at the integration points
  	for(size_t i = 0; i < (*vIterFEs)->IntegrationPoints(); i++)
  	{
  	  try
  	  {
  	    (*vIterFEs)->dN_AtIntegrationPoint( sf, i );
  	    sum = 0.;
        if ( verbose_ ) cout << "\n Integration Point: " << i;
  	    for(size_t j = 0; j < sf.Rows(); j++)  
  	     sum += sf.RowSum(j);
        if ( verbose_ ) cout << "\nShape functions should be 0. sum: " << sum;
  	    _equal( sum, 0., 1.e-14 );
  	  }
  	  catch( std::range_error e )
  	  {
  	    //fail the test but do not die
        if ( verbose_ ) cout << "\nThere was an exception: " << e.what() << endl;
        _test(false);
  	  }
  	}	  
    //check that they sum up at the nodes
  	for(size_t i = 0; i < (*vIterFEs)->Nodes(); i++)
  	{
  	  try
  	  {  	  
    	  (*vIterFEs)->dN_AtNode( sf, i );
    	  sum = 0.;
        if ( verbose_ ) cout << "\n Node: " << i;
    	  for(size_t j = 0; j < sf.Rows(); j++)  
    	    sum += sf.RowSum(j);
        if ( verbose_ ) cout << "\nShape functions should be 0. sum: " << sum;
    	  _equal( sum, 0., 1.e-14 ); 
  	  }
  	  catch( std::range_error e )
  	  {
  	    //fail the test but do not die
        if ( verbose_ ) cout << "\nThere was an exception: " << e.what() << endl;
        _test(false);
  	  }
 	  
  	}
  	//check that they sum up at the barycenters
  	(*vIterFEs)->dN_AtBarycenter( sf );
  	sum = 0.;
  	for(size_t j = 0; j < sf.Rows(); j++)  
  	  sum += sf.RowSum(j);
    if ( verbose_ ) cout << "\nShape functions should be 0. sum: " << sum;
  	_equal( sum, 0., 1.e-14 );
  }	
}

} // end namespace 
