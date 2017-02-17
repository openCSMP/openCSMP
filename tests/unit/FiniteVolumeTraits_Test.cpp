#include "FiniteVolumeTraits_Test.h"
#include "IsoparametricLinearLineElement.h"
#include "IsoparametricLinearLineElement.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearQuadrilateral.h"
#include "IsoparametricLinearTetrahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPrism.h"
#include "Element.h"
#include "Point.h"
#include "VSet.h"
#include "Test.h"
#include "Model.h"

#define TIMES 1

using namespace std;

//This test only verifies values using single precision results.
//Rhino was used to measure: facet area, facet normal, sector volume,
//sector integration point (rst to xyz test), facet integration point (rst to xyz test),
//and total volume.
//Rhino settings used for this test were:
//Precision:
//   - Absolute Tolerance: 1.e-7
//   - Relative Tolerance: 1.e-7
//   - Angle Tolerance: 0.001
//Grid Spacing: 10
//Units: milimeters
//Display precision: 1.0000000

namespace csmp {


namespace {
    template<size_t dim>
    void dumpVector(const char* text, const Point<dim>& v)
    {
		if ( dim==3 )
		{
		    _info(text << "[" << v[0] << ", " << v[1] << ", " << v[2] << "]");
		}
		else
        {
		    _info(text << "[" << v[0] << ", " << v[1] << "]");
        }
    }
}

/**  Method:


FiniteVolumeTraits_Test<dim>::FiniteVolumeTraits_Test()



Description:
Constructor. This function specifies which specific modalities of the test are on.
@section arguments Input Arguments


@section application Application
The idea of this function is to control the test functionality in one place, in an easy, fast way.

tested: is a test function*/
FiniteVolumeTraits_Test::FiniteVolumeTraits_Test()
 : m_bTestFacetAreas(true),
 	m_bTestParametricFacetArea(true),
	m_bTestFacetNormals(true),
	m_bTestParametricFacetNormals(true),
	m_bTestSectorVolumes(true),
	m_bProjectionOnFacetNormal(false),
	/*this test does not pass, because the parametric specification of the integral does not work as expected
	New Comment by J.E.M.: Projecting a vector native of one reference frame onto a vector native of another
	is not correct.  First the vector vVariable should have been transformed to parametric coordinates
	and then the test would have been valid. The test was inconsistent.
	*/
	m_bPropertyValueAt(true),
	m_bRSTToXYZ(true)
{
 }

/**  Method:


void FiniteVolumeTraits_Test<dim>::Test_CreateVSet()



Description:
This function tests the creation of a VSet (not related to the  class).

@section arguments Input Arguments


tested: is a test function*/
void FiniteVolumeTraits_Test::Test_CreateVSet()
{
    const size_t iNodes(9);
    const size_t iNrOfElements(4);

  	IsoparametricLinearQuadrilateral  iso_quadrilateral;
  	VSet<2U>  vset( 4U, 4U, ISOPARAMETRIC_LINEAR_QUADRILATERAL, iNodes, iNrOfElements );

  	//define nodes
  	std::deque<double64> px(iNodes);
  	std::deque<double64> py(iNodes);
  	std::deque<double64> pz(iNodes);

  	px[0]=0.;py[0]=0.;pz[0]=0.;
  	px[1]=1.;py[1]=0.;pz[1]=0.;
  	px[2]=2.;py[2]=0.;pz[2]=0.;
  	px[3]=0.;py[3]=1.;pz[3]=0.;
  	px[4]=1.;py[4]=1.;pz[4]=0.;
  	px[5]=2.;py[5]=1.;pz[5]=0.;
  	px[6]=0.;py[6]=2.;pz[6]=0.;
  	px[7]=1.;py[7]=2.;pz[7]=0.;
  	px[8]=2.;py[8]=2.;pz[8]=0.;

  	//load nodes
  	vset.AddXYZ( px, py, pz );

    //define elements
    std::deque<std::vector<size_t> > deqElements(iNrOfElements);
    std::vector<size_t> vecNodes(4);
    vecNodes[0]=0;
    vecNodes[1]=1;
    vecNodes[2]=4;
    vecNodes[3]=3;
    deqElements[0]=vecNodes;
    vecNodes[0]=1;
    vecNodes[1]=2;
    vecNodes[2]=5;
    vecNodes[3]=4;
    deqElements[1]=vecNodes;
    vecNodes[0]=4;
    vecNodes[1]=5;
    vecNodes[2]=8;
    vecNodes[3]=7;
    deqElements[2]=vecNodes;
    vecNodes[0]=3;
    vecNodes[1]=4;
    vecNodes[2]=7;
    vecNodes[3]=6;
    deqElements[3]=vecNodes;

    vset.AddPlist( deqElements.begin(),deqElements.end());

    //define neighbors
    std::deque<std::vector<long64> > deqElementNeighbors(iNrOfElements);
    std::vector<long64> vecNeighbors(4);
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=2;
    vecNeighbors[2]=4;
    vecNeighbors[3]=LEFT_OUTSIDE;
    deqElementNeighbors[0]=vecNeighbors;
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=3;
    vecNeighbors[3]=1;
    deqElementNeighbors[1]=vecNeighbors;
    vecNeighbors[0]=2;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=TOP_OUTSIDE;
    vecNeighbors[3]=4;
    deqElementNeighbors[2]=vecNeighbors;
    vecNeighbors[0]=1;
    vecNeighbors[1]=3;
    vecNeighbors[2]=TOP_OUTSIDE;
    vecNeighbors[3]=LEFT_OUTSIDE;
    deqElementNeighbors[3]=vecNeighbors;
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());

  	vset.AddBFlag( 0, CNR1 );
    vset.AddBFlag( 1, BOTTOM_OUTSIDE );
    vset.AddBFlag( 2, CNR2 );
    vset.AddBFlag( 3, LEFT_OUTSIDE );
    vset.AddBFlag( 5, RIGHT_OUTSIDE );
    vset.AddBFlag( 6, CNR4 );
	  vset.AddBFlag( 7, TOP_OUTSIDE );
    vset.AddBFlag( 8, CNR3 );

    vset.CheckFix();
    vset.Out();

    //create the super group
   Model<2>  superGroup( vset, "fe_test_variables.txt", true );

 	 VectorVariable<2> vVariable(PLAIN,PLAIN, sqrt(2.)/2., sqrt(2.)/2.);
   _info(vVariable.Length());

  superGroup.InputPropertyValue( "velocity", vVariable );
}


/**  Method:


FiniteVolumeTraits_Test<dim>::~FiniteVolumeTraits_Test()



Description:
Destructor.
@section arguments Input Arguments
None.
@section application Application
Nothing to destroy.

tested: is a test function*/
FiniteVolumeTraits_Test::~FiniteVolumeTraits_Test()
 {
 }

/**  Method:


void FiniteVolumeTraits_Test<dim>::run()



Description:
Specifies the procedure of the test.
@section arguments Input Arguments
None.
@section application Application
Runs the tests.

tested: is a test function*/
void FiniteVolumeTraits_Test::run() // runs all the tests for the class (register other methods)
 {
  _info("---------------------------------------------------------------------");
  _info("Starting...FiniteVolumeTraits_Test::run()");

	IsoparametricLinearLineElement_Test(1.e-7,1.e-7);
	IsoparametricLinearTriangle_Test<2U>(1.e-7,1.e-7);
	IsoparametricLinearTriangle_Test<3U>(1.e-7,1.e-7);
	IsoparametricLinearQuadrilateral_Test<2U>(1.e-7,1.e-7);
	IsoparametricLinearQuadrilateral_Test<3U>(1.e-7,1.e-7);
	IsoparametricLinearTetrahedron_Test(1.e-7,1.e-7);
	IsoparametricLinearPyramid_Test(1.e-7,1.e-7);
	IsoparametricLinearPrism_Test(1.e-7,1.e-7);

	//Reviewed by J.E.M 08-09-2010.
	Test_UnitaryIsoparametricLinearHexahedron(1.e-7,1.e-7);
	Test_IsoparametricLinearHexahedron1(1.e-7,1.e-7);
	Test_IsoparametricLinearHexahedron2(1.e-7,1.e-7);



  _info("Done...FiniteVolumeTraits_Test::run()");
  _info("---------------------------------------------------------------------");

 }


/**  Method:


void FiniteVolumeTraits_Test<dim>::IsoparametricLinearLineElement_Test(double64 fTolerance, double64 fToleranceInternal)



Description:
Specifies the procedure of the test.

The tested element is supposed to be a uniform transformation of the original object by:
Translate(-3,0,0)
Rotate(N0, 40degrees) -2D
Shear(N0, N1, 35degrees)

This test does the following:
-construct the element.
-test the facet area, to the one measured in rhino
-test the parametric facet area
-test the volume of the sector, and check if it corresponds to the measured one
-test the sum of the volumes
-test the projection on the facet normal of a given property
-convert the facet integration points to physical space and check if they correspond to the measured values
-convert the sector integration points to physical space and check if they correspond to the measured values

@section arguments Input Arguments
double64 fTolerance - specifies the tolerance of the test
@section application Application
Runs the test for the Isoparametric Linear Bar element, for all dimensions.

tested: is a test function*/
void FiniteVolumeTraits_Test::IsoparametricLinearLineElement_Test(double64 fTolerance, double64)
 {
	FiniteElement* feptr			= new IsoparametricLinearLineElement(1U);
	FiniteVolumeStencil<1U>* fvptr	= new FiniteVolumeStencil<1U>("ISOPARAMETRIC_LINEAR_BAR");

	Element<1U>  elmt_( feptr );
 	elmt_.Assign( fvptr );

	Node<1U>  node1, node2, node3;

	node1.Idx( 1 );
	node1.x( -3. );

	node2.Idx( 2 );
	node2.x( -2.3680806 );

    const size_t dim(1);

	if(dim > 1)
	{
		node1.y( 0. );
		node2.y( 2.3583554 );
	}
	if(dim == 3)
	{
		node1.z( 0. );
		node2.z( 0. );
	}

	elmt_.Idx( 1 );
	elmt_.Assign( 0U, &node1 );
	elmt_.Assign( 1U, &node2 );

	//test facet area
	//testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestFacetAreas )
	  {

    double64 fArea = ( elmt_ ).FacetArea( 0U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

	  }

	if ( m_bTestParametricFacetArea )
	  {
    double64 fArea = ( elmt_ ).ParametricFacetArea( 0U );
		_info("Area is: " << fArea);
		_equal( fArea, elmt_.FV_Stencil()->FacetIntegrationWeight(0U,0U), fTolerance );
	  }

	//test sector volumes
	//testing: fT   SectorVolume( size_t sector ) const;
    double64 fVolSum(0.f);

    if ( m_bTestSectorVolumes )
	  {
    double64 fSectorVolume = ( elmt_ ).SectorVolume( 0U );
		_info("Volume is: " << setprecision(15) << fSectorVolume);
		if(dim == 1)
			_equal( fSectorVolume, 0.3159597, fTolerance );
		else
			_equal( fSectorVolume, 1.2207746, fTolerance );
		fVolSum += fSectorVolume;

    fSectorVolume = ( elmt_ ).SectorVolume( 1U );
		_info("Volume is: " << fSectorVolume);
		if(dim == 1)
			_equal( fSectorVolume, 0.3159597, fTolerance );
		else
			_equal( fSectorVolume, 1.2207746, fTolerance );
		fVolSum += fSectorVolume;

		//test if sector volumes add up to volume of element
		_info("Volume of the Element: " << elmt_.Volume());
		_info("vs Volume Sum: " << fVolSum);
		_equal(fVolSum, elmt_.Volume(), fTolerance);
	  }

	if (m_bProjectionOnFacetNormal)
      {
	   const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
	   for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
	   {
	     //const CSPINDEX& prop_key
			 VectorVariable<1U> vVariable;
	     for(size_t iD= 0U; iD < dim; iD++) {
	     	vVariable.Flag(iD)=PLAIN;
	     	vVariable(iD)=3.;
	     }

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct( ( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
			_equal(fProjectionVal, fProjectionValP, fTolerance);
	   }
     }

   if(m_bRSTToXYZ)
   {
   	Point<dim> vecRST, vecXYZ;

	//check facet integration points
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(0U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.6840403);
    _equal(vecXYZ[0], -2.6840403, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.1791777);
      _equal(vecXYZ[1], 1.1791777, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //check sector integration points
    //sector 0
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(0U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.8420201);
    _equal(vecXYZ[0], -2.8420201, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 0.5895889);
      _equal(vecXYZ[1], 0.5895889, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //sector 1
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(1U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.5260604);
    _equal(vecXYZ[0], -2.5260604, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.7687666);
      _equal(vecXYZ[1], 1.7687666, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }
   }

   delete feptr;
   delete fvptr;

 } // end test


/**  Method:


void FiniteVolumeTraits_Test<dim>::IsoparametricLinearTriangle_Test(
                     double64 fTolerance, double64 fToleranceInternal)



Description:

Specifies the procedure of the test.

The tested element is supposed to be a uniform transformation of the original object by:
Translate(-3,0,0)
Rotate(N0, 40degrees) -2D
Shear(N0, N1, 35degrees)

This test does the following:
-construct the element.
-test the facet area, to the one measured in rhino
-test the facet normals
-test the parametric facet area
-test the volume of the sector, and check if it corresponds to the measured one
-test the sum of the volumes
-test the projection on the facet normal of a given property
-convert the facet integration points to physical space and check if they correspond to the measured values
-convert the sector integration points to physical space and check if they correspond to the measured values

@section arguments Input Arguments

double64 fTolerance - specifies the tolerance of the test

@section application Application

Runs the test for the Isoparametric Linear Triangle element, for all dimensions != 1.


tested: is a test function*/
template<size_t dim>
void FiniteVolumeTraits_Test::IsoparametricLinearTriangle_Test(double64 fTolerance, double64 fToleranceInternal)
{
	FiniteElement* feptr			= new IsoparametricLinearTriangle(dim);
	FiniteVolumeStencil<dim>* fvptr	= new FiniteVolumeStencil<dim>("ISOPARAMETRIC_LINEAR_TRIANGLE");

	Element<dim>  elmt_( feptr );
	elmt_.Assign( fvptr );

	Node<dim>  node1, node2, node3;
	{
		node1.Idx( 1 );
		node1.x( -3. );
		node1.y( 0. );

		node2.Idx( 2 );
		node2.x( -2.6840403 );
		node2.y( 1.1791777 );

		node3.Idx( 3 );
		node3.x( -3.6427876 );
		node3.y( 0.7660444 );

	}

	elmt_.Idx( 1 );
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );

	//test facet area
	//testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestFacetAreas )
	  {
    double64 fArea = ( elmt_ ).FacetArea( 0U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.2733262, fTolerance );

    fArea = ( elmt_ ).FacetArea( 1U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.3287479, fTolerance );

    fArea = ( elmt_ ).FacetArea( 2U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.3399482, fTolerance );

		//test mapped vs original
    double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.2733262, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.3287479, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.3399482, fTolerance );

		//test mapped vs computed
		_info("Test mapped area vs computed.");
	  for(size_t i = 0; i < 3; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearTriangle_Test: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 3; i++)
        j += ( elmt_ ).FacetArea( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 3; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		 file.close();
	  }

	//test facet area
	//fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
	if ( m_bTestFacetNormals )
	  {
		csmp::Point<dim> vecNormal, vecNormalMapped;

		//ignore return parameter, jacobian

		//facet 0
    vecNormal = ( elmt_ ).FacetNormal(0U);
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);

        dumpVector<dim>("Normal is: ", vecNormal);

		_equal( vecNormal[0], 0.2151953, fTolerance );
		_equal( vecNormal[1], 0.9765710, fTolerance );
		_equal( vecNormalMapped[0], 0.2151953, fTolerance );
		_equal( vecNormalMapped[1], 0.9765710, fTolerance );
		if ( dim==3 )
		  _equal( vecNormal[2], 0., fTolerance );
		if ( dim==3 )
		  _equal( vecNormalMapped[2], 0., fTolerance );

        dumpVector<dim>("Mapped Normal is: ", vecNormalMapped);

		//facet 1
    vecNormal = ( elmt_ ).FacetNormal(1U);
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);

        dumpVector<dim>("Normal is: ", vecNormal);

		_equal( vecNormal[0], -0.9861773, fTolerance );
		_equal( vecNormal[1], -0.1656933, fTolerance );
		_equal( vecNormalMapped[0], -0.9861773, fTolerance );
		_equal( vecNormalMapped[1], -0.1656933, fTolerance );
		if ( dim==3 )
		  _equal( vecNormal[2], 0, fTolerance );
		if ( dim==3 )
		  _equal( vecNormalMapped[2], 0, fTolerance );

        dumpVector<dim>("Mapped Normal is: ", vecNormalMapped);

		//facet 2
    vecNormal = ( elmt_ ).FacetNormal(2U);
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);

        dumpVector<dim>("Normal is: ", vecNormal);

		_equal( vecNormal[0], 0.7806635, fTolerance );
		_equal( vecNormal[1], -0.6249516, fTolerance );
		_equal( vecNormalMapped[0], 0.7806635, fTolerance );
		_equal( vecNormalMapped[1], -0.6249516, fTolerance );
		if ( dim==3 )
		  _equal( vecNormal[2], 0, fTolerance );
	  if ( dim==3 )
		  _equal( vecNormalMapped[2], 0, fTolerance );

        dumpVector<dim>("Mapped Normal is: ", vecNormalMapped);
      }

    //test mapped vs computed
    _info("Test mapped normal vs computed.");
	  for(size_t i = 0; i < 3; i++)
	    for(size_t j = 0; j < 3; j++)
      _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearTriangle_Test<3>: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 3; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 3; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		 file.close();

	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	  {
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
        double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
				_info("Area is: " << fArea);
		    _equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
      }
	  }

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	  {
		Point<dim>  vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
		    //ignore return parameter, jacobian
        vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);
        dumpVector<dim>("Normal is: ", vecNormal);

		    _equal( vecNormal[0],  elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
		    _equal( vecNormal[1],  elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
		    if ( dim==3 )
		      _equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		  }
	  }

	//test sector volumes
	//testing: fT   SectorVolume( size_t sector ) const;
	if ( m_bTestSectorVolumes )
	  {
		double64 fVolSum(0.);

    double64 fSectorVolume = ( elmt_ ).SectorVolume( 0U );
		_info("Volume is: " << fSectorVolume);
		_equal( fSectorVolume, 1./6., fTolerance );
		fVolSum+=fSectorVolume;

    fSectorVolume = ( elmt_ ).SectorVolume( 1U );
		_info("Volume is: " << fSectorVolume);
		_equal( fSectorVolume, 1./6., fTolerance );
        fVolSum+=fSectorVolume;

    fSectorVolume = ( elmt_ ).SectorVolume( 2U );
		_info("Volume is: " << fSectorVolume);
		_equal( fSectorVolume, 1./6., fTolerance );
	    fVolSum+=fSectorVolume;

		_info("Volume of the Element: " << elmt_.Volume() << " vs Volume Sum: " << fVolSum);
		_equal(fVolSum, elmt_.Volume(), fTolerance);

	  }

   if(m_bProjectionOnFacetNormal)
   {
	   const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
	   for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
	   {
	     VectorVariable<dim> vVariable;
	     for(size_t iD= 0U; iD < dim; iD++) {
	     	vVariable.Flag(iD)=PLAIN;
	     	vVariable(iD)=3.;
	     }

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct(( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
		_equal(fProjectionVal, fProjectionValP, fTolerance);
	}
   }

   if(m_bRSTToXYZ)
   {
   	Point<dim> vecRST, vecXYZ;

	  //check facet integration points
   	//facet 0
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(0U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.9754814);
    _equal(vecXYZ[0], -2.9754814, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.6189981);
      _equal(vecXYZ[1], 0.6189981, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }


    //facet 1
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(1U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.1361783);
    _equal(vecXYZ[0], -3.1361783, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.8105092);
      _equal(vecXYZ[1], 0.8105092, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //facet 2
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(2U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.2151682);
    _equal(vecXYZ[0], -3.2151682, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.5157148);
      _equal(vecXYZ[1], 0.5157148, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //check sector integration points
    //sector 0
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(0U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.0635499);
    _equal(vecXYZ[0],  -3.0635499, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.3782376);
      _equal(vecXYZ[1], 0.3782376, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }


    //sector 1
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(1U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.9319000);
    _equal(vecXYZ[0], -2.9319000, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.8695617);
      _equal(vecXYZ[1], 0.8695617, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //sector 2
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(2U, 0U);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.3313780);
    _equal(vecXYZ[0], -3.3313780, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.6974228);
      _equal(vecXYZ[1], 0.6974228, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

   }

   delete feptr;
   delete fvptr;

} //end test




/**  Method:


void FiniteVolumeTraits_Test<dim>::IsoparametricLinearQuadrilateral_Test(double64 fTolerance, double64 fToleranceInternal)



Description:
Specifies the procedure of the test.

The tested element is supposed to be a uniform transformation of the original object by:
Translate(-3,0,0)
Rotate(N0, 40degrees) -2D
Shear(N0, N1, 35degrees)

This test does the following:
-construct the element.
-test the facet area, to the one measured in rhino
-test the facet normals
-test the parametric facet area
-test the volume of the sector, and check if it corresponds to the measured one
-test the sum of the volumes
-test the projection on the facet normal of a given property
-convert the facet integration points to physical space and check if they correspond to the measured values
-convert the sector integration points to physical space and check if they correspond to the measured values

@section arguments Input Arguments
double64 fTolerance - specifies the tolerance of the test
@section application Application
Runs the test for the Isoparametric Linear Quadrilateral element, for all dimensions != 1.

tested: is a test function*/
template<size_t dim>
void FiniteVolumeTraits_Test::IsoparametricLinearQuadrilateral_Test(double64 fTolerance, double64 fToleranceInternal)
{

   FiniteElement* feptr			    = new IsoparametricLinearQuadrilateral(dim);
   FiniteVolumeStencil<dim>* fvptr	= new FiniteVolumeStencil<dim>("ISOPARAMETRIC_LINEAR_QUADRILATERAL");

   csmp::Element<dim>  elmt_( feptr );
   elmt_.Assign( fvptr );

   Node<dim>  node1, node2, node3, node4;
   {
	  node1.Idx( 1 );
	  node1.x( -3. );
	  node1.y( 0. );

	  node2.Idx( 2 );
	  node2.x( -2.3680806 );
	  node2.y( 2.3583554 );

	  node3.Idx( 3 );
	  node3.x( -3.6536558 );
	  node3.y( 3.8904443 );

	  node4.Idx( 4 );
	  node4.x( -4.2855752 );
	  node4.y( 1.5320889 );
    }

  	if ( dim == 3U )
  	{
	  node1.z( 0. );
	  node2.z( 0. );
	  node3.z( 0. );
	  node4.z( 0. );
	}

	elmt_.Idx( 1 );
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );
	elmt_.Assign( 3, &node4 );

   //test facet area
   //testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
    if ( m_bTestFacetAreas )
	  {
      double64 fArea = ( elmt_ ).FacetArea( 0U );
      _info("Area is: " << fArea);
	    _equal( fArea, 1., fTolerance );

      fArea = ( elmt_ ).FacetArea( 1U );
			_info("Area is: " << fArea);
	    _equal( fArea, 1.2207746, fTolerance );

      fArea = ( elmt_ ).FacetArea( 2U );
      _info("Area is: " << fArea);
	   _equal( fArea, 1., fTolerance );

      fArea = ( elmt_ ).FacetArea( 3U );
      _info("Area is: " << fArea);
	    _equal( fArea, 1.2207746, fTolerance );

	    //mapped area
      double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
			_info("Mapped Area is: " << fAreaMapped);
	    _equal( fAreaMapped, 1., fTolerance );

      fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
			_info("Mapped Area is: " << fAreaMapped);
	    _equal( fAreaMapped, 1.2207746, fTolerance );

      fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
      _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 1., fTolerance );

      fAreaMapped = ( elmt_ ).FacetAreaMapped( 3U );
      _info("Mapped Area is: " << fAreaMapped);
	    _equal( fAreaMapped, 1.2207746, fTolerance );

	   //test mapped vs computed
		 _info("Test mapped area vs computed.");
	   for(size_t i = 0; i < 4; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearQuadrilateral_Test<" << dim<< ">: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 4; i++)
        j += ( elmt_ ).FacetArea( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 4; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		 file.close();
     }

   //test facet area
   //fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
   if ( m_bTestFacetNormals )
	 {
	   Point<dim> vecNormal, vecNormalMapped;
	   /*ignore return parameter, jacobian*/
     vecNormal = ( elmt_ ).FacetNormal(0U);
     vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);

        dumpVector<dim>("Normal is: ", vecNormal);

	   _equal( vecNormal[0], 0.7660444, fTolerance );
	   _equal( vecNormal[1], 0.6427876, fTolerance );
	   if ( dim==3 )
	     _equal( vecNormal[2], 0., fTolerance );
	   _equal( vecNormalMapped[0], 0.7660444, fTolerance );
	   _equal( vecNormalMapped[1], 0.6427876, fTolerance );
	   if ( dim==3 )
	     _equal( vecNormalMapped[2], 0., fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(1U);
     vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);

        dumpVector<dim>("Normal is: ", vecNormal);

	   _equal( vecNormal[0], -0.9659258, fTolerance );
	   _equal( vecNormal[1], 0.2588190, fTolerance );
	   if ( dim==3 )
	     _equal( vecNormal[2], 0., fTolerance );
	   _equal( vecNormalMapped[0], -0.9659258, fTolerance );
	   _equal( vecNormalMapped[1], 0.2588190, fTolerance );
	   if ( dim==3 )
	     _equal( vecNormalMapped[2], 0., fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(2U);
     vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);

        dumpVector<dim>("Normal is: ", vecNormal);

	   _equal( vecNormal[0], -0.7660444, fTolerance );
	   _equal( vecNormal[1], -0.6427876, fTolerance );
	   if ( dim==3 )
	     _equal( vecNormal[2], 0., fTolerance );
	   _equal( vecNormalMapped[0], -0.7660444, fTolerance );
	   _equal( vecNormalMapped[1], -0.6427876, fTolerance );
	   if ( dim==3 )
	     _equal( vecNormalMapped[2], 0., fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(3U);
     vecNormalMapped = ( elmt_ ).FacetNormalMapped(3U);

        dumpVector<dim>("Normal is: ", vecNormal);

	  _equal( vecNormal[0], 0.9659258, fTolerance );
	  _equal( vecNormal[1], -0.2588190, fTolerance );
	  if ( dim==3 )
	    _equal( vecNormal[2], 0., fTolerance );
	  _equal( vecNormalMapped[0], 0.9659258, fTolerance );
	  _equal( vecNormalMapped[1], -0.2588190, fTolerance );
	  if ( dim==3 )
	    _equal( vecNormalMapped[2], 0., fTolerance );

	    //test mapped vs computed
		_info("Test mapped normal vs computed.");
	  for(size_t i = 0; i < 4; i++)
	    for(size_t j = 0; j < 3; j++)
      _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearQuadrilateral_Test<3>: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 4; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 4; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		 file.close();
	 }

	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	  {
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
        double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
				_info("Area is: " << fArea);
		    _equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
          }
	  }

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	  {
		Point<dim> vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
		    //ignore return parameter, jacobian
        vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);

        dumpVector<dim>("Normal is: ", vecNormal);

		    _equal( vecNormal[0],  elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
		    _equal( vecNormal[1],  elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
		    if ( dim==3 )
		      _equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		  }
	  }

   //test sector volumes
   //testing: fT   SectorVolume( size_t sector ) const;
   if ( m_bTestSectorVolumes )
	 {
	   double64 fVolSum(0.);

     double64 fSectorVolume = ( elmt_ ).SectorVolume( 0U );
		 _info("Volume is: " << fSectorVolume);
	   _equal( fSectorVolume, 1., fTolerance );
	   fVolSum+=fSectorVolume;

     fSectorVolume = ( elmt_ ).SectorVolume( 1U );
     _info("Volume is: " << fSectorVolume);
	   _equal( fSectorVolume, 1., fTolerance );
     fVolSum+=fSectorVolume;

     fSectorVolume = ( elmt_ ).SectorVolume( 2U );
		 _info("Volume is: " << fSectorVolume);
	   _equal( fSectorVolume, 1., fTolerance );
	   fVolSum+=fSectorVolume;

     fSectorVolume = ( elmt_ ).SectorVolume( 3U );
		 _info("Volume is: " << fSectorVolume);
	   _equal( fSectorVolume, 1., fTolerance );
     fVolSum+=fSectorVolume;

		 _info("Volume of the Element: " << elmt_.Volume() << " vs Volume Sum: " << fVolSum);
	   _equal(fVolSum, elmt_.Volume(), fTolerance);
	 }

   if(m_bProjectionOnFacetNormal)
   {
	   const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
	   for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
	   {
	     VectorVariable<dim> vVariable;
	     vVariable=3.;

			 _info("LENGTH::::::::" << vVariable.Length());

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct(( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
		  _equal(fProjectionVal, fProjectionValP, fTolerance);
	   }
   }

     if(m_bRSTToXYZ)
   {
   	elmt_.CoordinateMatrix();

   	Point<dim> vecRST, vecXYZ;

	  //check facet integration points
   	//facet 0
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(0U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.0054341);
    _equal(vecXYZ[0], -3.0054341, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 1.5621999);
      _equal(vecXYZ[1], 1.5621999, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //facet 1
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(1U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.1688480);
    _equal(vecXYZ[0],  -3.1688480, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 2.5348110);
      _equal(vecXYZ[1], 2.5348110, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //facet 2
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(2U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.6482217);
    _equal(vecXYZ[0], -3.6482217, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 2.3282444);
      _equal(vecXYZ[1], 2.3282444, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //facet 3
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(3U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.4848078);
    _equal(vecXYZ[0], -3.4848078, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 1.3556333);
      _equal(vecXYZ[1], 1.3556333, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //check sector integration points
    //sector 0
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(0U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.16341395);
    _equal(vecXYZ[0], -3.16341395, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.972611073);
      _equal(vecXYZ[1], 0.972611073, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //sector 1
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(1U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.8474542);
    _equal(vecXYZ[0],  -2.8474542, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 2.1517888);
      _equal(vecXYZ[1], 2.1517888, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //sector 2
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(2U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.4902418);
    _equal(vecXYZ[0], -3.4902418, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 2.9178332);
      _equal(vecXYZ[1], 2.9178332, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }

    //sector 3
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(3U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -3.80620156);
    _equal(vecXYZ[0], -3.80620156, fTolerance);
    if ( dim > 1 )
    {
      _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 1.73865552);
      _equal(vecXYZ[1], 1.73865552, fTolerance);
    }
    if ( dim > 2 )
    {
      _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.);
      _equal(vecXYZ[2], 0., fTolerance);
    }
   }

   delete feptr;
   delete fvptr;

} //end test




/**  Method:


void FiniteVolumeTraits_Test<3>::IsoparametricLinearTetrahedron_Test(double64 fTolerance, double64 fToleranceInternal)



Description:
Specifies the procedure of the test.

The tested element is supposed to be a uniform transformation of the original object by:
Translate(-3,0,0)
Rotate(N0, 40degrees) -3D
Shear(N0, N1, 35degrees)

This test does the following:
-construct the element.
-test the facet area, to the one measured in rhino
-test the facet normals
-test the parametric facet area
-test the volume of the sector, and check if it corresponds to the measured one
-test the sum of the volumes
-test the projection on the facet normal of a given property
-convert the facet integration points to physical space and check if they correspond to the measured values
-convert the sector integration points to physical space and check if they correspond to the measured values

@section arguments Input Arguments
double64 fTolerance - specifies the tolerance of the test
@section application Application
Runs the test for the Isoparametric Linear Tetrahedron element, for 3D.

tested: is a test function*/
void FiniteVolumeTraits_Test::IsoparametricLinearTetrahedron_Test( double64 fTolerance, double64 fToleranceInternal )
{
   FiniteElement* feptr			    = new IsoparametricLinearTetrahedron();
   FiniteVolumeStencil<3U>* fvptr	= new FiniteVolumeStencil<3U>("ISOPARAMETRIC_LINEAR_TETRAHEDRON");

   Element<3U>  elmt_( feptr );
   elmt_.Assign( fvptr );

   Node<3U>  node1, node2, node3, node4;
   {
	  node1.Idx( 1 );
	  node1.x( -3. );
	  node1.y( 0. );
	  node1.z( 0. );

	  node2.Idx( 2 );
	  node2.x( -2.1691958 );
	  node2.y( 0.8889046 );
	  node2.z( -0.0995164 );

	  node3.Idx( 3 );
	  node3.x( -3.1084262 );
	  node3.y( 0.8216875 );
	  node3.z( 0.6190235 );

	  node4.Idx( 4 );
	  node4.x( -2.8077164 );
	  node4.y( -0.5393590 );
	  node4.z( 0.7790420 );
	}

	elmt_.Idx( 1 );
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );
	elmt_.Assign( 3, &node4 );

   //test facet area
   //testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
   if ( m_bTestFacetAreas )
	 {
     double64 fArea = ( elmt_ ).FacetArea( 0U );
     _info("Area is: " << setprecision(15) << fArea);
	   _equal( fArea, 0.0979896825, fTolerance );

     fArea = ( elmt_ ).FacetArea( 1U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.0760544993, fTolerance );

     fArea = ( elmt_ ).FacetArea( 2U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.0917751512, fTolerance );

     fArea = ( elmt_ ).FacetArea( 3U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.109391786, fTolerance );

       fArea = ( elmt_ ).FacetArea( 4U );
     _info("Area is: " << fArea);
	   _equal( fArea, 0.0491005164, fTolerance );

       fArea = ( elmt_ ).FacetArea( 5U );
     _info("Area is: " << fArea);
	   _equal( fArea, 0.0722797155, fTolerance );

	   //mapped area
     double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		 _info("Mapped Area is: " << setprecision(15) << fAreaMapped);
	   _equal( fAreaMapped, 0.0979896825, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.0760544993, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.0917751512, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 3U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.109391786, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 4U );
     _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.0491005164, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 5U );
     _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.0722797155, fTolerance );

	  //test mapped vs computed
		_info("Test mapped area vs computed.");
	  for(size_t i = 0; i < 6; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearTetrahedron_Test<3>: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 6; i++)
        j += ( elmt_ ).FacetArea( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 6; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		 file.close();
     }

   //test facet area
   //fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
   if ( m_bTestFacetNormals )
	 {
     Point<3U> vecNormal;

	   /*ignore return parameter, jacobian*/
     vecNormal = ( elmt_ ).FacetNormal(0U);
        dumpVector<3>("Normal is: ", vecNormal);
	   _equal( vecNormal[0], 0.8254335, fTolerance );
	   _equal( vecNormal[1], 0.2423141, fTolerance );
	   _equal( vecNormal[2], 0.5098465, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(1U);
        dumpVector<3>("Normal is: ", vecNormal);
	   _equal( vecNormal[0], -0.8835900, fTolerance );
	   _equal( vecNormal[1], 0.2535847, fTolerance );
	   _equal( vecNormal[2], 0.3936540, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(2U);
        dumpVector<3>("Normal is: ", vecNormal);
	   _equal( vecNormal[0], -0.1490923, fTolerance );
	    _equal( vecNormal[1], -0.4688692, fTolerance );
	   _equal( vecNormal[2], -0.8705936, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(3U);
        dumpVector<3>("Normal is: ", vecNormal);
	   _equal( vecNormal[0], 0.6091387, fTolerance );
	   _equal( vecNormal[1], -0.0522324, fTolerance );
	   _equal( vecNormal[2], 0.7913418, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(4U);
        dumpVector<3>("Normal is: ", vecNormal);
	   _equal( vecNormal[0], -0.2902046, fTolerance );
	   _equal( vecNormal[1], -0.5999545, fTolerance );
	   _equal( vecNormal[2], 0.7455440, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(5U);
        dumpVector<3>("Normal is: ", vecNormal);
	   _equal( vecNormal[0], 0.7325956, fTolerance );
	   _equal( vecNormal[1], -0.6743846, fTolerance );
	   _equal( vecNormal[2], 0.0922449, fTolerance );

	   //mapped normals
	   Point<3U> vecNormalMapped;
     vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);
        dumpVector<3>("Mapped Normal is: ", vecNormalMapped);
	   _equal( vecNormalMapped[0], 0.8254335, fTolerance );
	   _equal( vecNormalMapped[1], 0.2423141, fTolerance );
	   _equal( vecNormalMapped[2], 0.5098465, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);
        dumpVector<3>("Mapped Normal is: ", vecNormalMapped);
	   _equal( vecNormalMapped[0], -0.8835900, fTolerance );
	   _equal( vecNormalMapped[1], 0.2535847, fTolerance );
	   _equal( vecNormalMapped[2], 0.3936540, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);
        dumpVector<3>("Mapped Normal is: ", vecNormalMapped);
	   _equal( vecNormalMapped[0], -0.1490923, fTolerance );
	    _equal( vecNormalMapped[1], -0.4688692, fTolerance );
	   _equal( vecNormalMapped[2], -0.8705936, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(3U);
        dumpVector<3>("Mapped Normal is: ", vecNormalMapped);
	   _equal( vecNormalMapped[0], 0.6091387, fTolerance );
	   _equal( vecNormalMapped[1], -0.0522324, fTolerance );
	   _equal( vecNormalMapped[2], 0.7913418, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(4U);
        dumpVector<3>("Mapped Normal is: ", vecNormalMapped);
	   _equal( vecNormalMapped[0], -0.2902046, fTolerance );
	   _equal( vecNormalMapped[1], -0.5999545, fTolerance );
	   _equal( vecNormalMapped[2], 0.7455440, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(5U);
        dumpVector<3>("Mapped Normal is: ", vecNormalMapped);
	   _equal( vecNormalMapped[0], 0.7325956, fTolerance );
	   _equal( vecNormalMapped[1], -0.6743846, fTolerance );
	   _equal( vecNormalMapped[2], 0.0922449, fTolerance );

	   //test mapped vs computed
		_info("Test mapped normal vs computed.");
	  for(size_t i = 0; i < 6; i++)
	    for(size_t j = 0; j < 3; j++)
      _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearTetrahedron_Test<3>: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 6; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < 6; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		 file.close();
     }

	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	  {
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
        double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
				_info("Area is: " << fArea);
		    _equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
          }
	  }

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	  {
		Point<3U> vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
		    //ignore return parameter, jacobian
        vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);
            dumpVector<3>("Normal is: ", vecNormal);

		    _equal( vecNormal[0], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
		    _equal( vecNormal[1], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
		    _equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		  }
	  }

   //test sector volumes
   //testing: fT   SectorVolume( size_t sector ) const;
    if ( m_bTestSectorVolumes )
	 {
       double64 fVolSum(0.);
	   double64 fSectorVolume(0.);
	   for ( size_t iSector = 0U; iSector < elmt_.FV_Stencil()->Sectors(); iSector++ )
	   {
       fSectorVolume = ( elmt_ ).SectorVolume( iSector );
				 //_info("Volume is: " << fSectorVolume);
	       _equal( fSectorVolume, 0.041666666666667, fTolerance );
	       fVolSum+=fSectorVolume;
	   }

		 _info("Volume of the Element: " << elmt_.Volume() << "vs Volume Sum: " << fVolSum);
	   _equal(fVolSum, elmt_.Volume(), fTolerance);
     }

      if(m_bProjectionOnFacetNormal)
     {
	   const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
	   for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
	   {
	     //const CSPINDEX& prop_key
	     VectorVariable<3> vVariable;
	     vVariable=3.;

			 _info("LENGTH::::::::" << vVariable.Length());

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct(( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
  		_equal(fProjectionVal, fProjectionValP, fTolerance);
	   }
   }

   if(m_bRSTToXYZ)
   {
   	elmt_.CoordinateMatrix();

	  Point<3U> vecRST, vecXYZ;

	//check facet integration points
   	//facet 0
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(0U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.6884298 );
    _equal(vecXYZ[0],-2.68834051, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.3606430);
    _equal(vecXYZ[1], 0.360205629, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.1581760);
    _equal(vecXYZ[2], 0.158239294, fTolerance);

    //facet 1
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(1U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.71243522);
    _equal(vecXYZ[0], -2.71243522, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " <<0.542802848);
    _equal(vecXYZ[1],0.542802848, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.29580007);
    _equal(vecXYZ[2], 0.29580007, fTolerance);

    //facet 2
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(2U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.89705838);
    _equal(vecXYZ[0], -2.89705838, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.345268485);
    _equal(vecXYZ[1],0.345268485, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.317914821);
    _equal(vecXYZ[2], 0.317914821, fTolerance);

    //facet 3
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(3U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.83023398);
    _equal(vecXYZ[0], -2.83023398, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.0428137098);
    _equal(vecXYZ[1],0.0428137098, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.353474481);
    _equal(vecXYZ[2], 0.353474481, fTolerance);

    //facet 4
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(4U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.64561082);
    _equal(vecXYZ[0], -2.64561082, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " <<0.240348073);
    _equal(vecXYZ[1],0.240348073, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.331359733);
    _equal(vecXYZ[2], 0.331359733, fTolerance);

    //facet 5
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(5U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.85432869);
    _equal(vecXYZ[0], -2.85432869, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " <<0.225410929);
    _equal(vecXYZ[1],0.225410929, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.491035259);
    _equal(vecXYZ[2],0.491035259, fTolerance);

    //check sector integration points
    //sector 0
    _info("Sector 0");
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(0U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.85390822);
    _equal(vecXYZ[0], -2.85390822, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " << 0.187071956);
    _equal(vecXYZ[1],0.187071956, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.207407149);
    _equal(vecXYZ[2], 0.207407149, fTolerance);

    //sector 1
    _info("Sector 1");
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(1U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.55389558);
    _equal(vecXYZ[0], -2.55389558, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " <<0.508065295);
    _equal(vecXYZ[1],0.508065295, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.171470682);
    _equal(vecXYZ[2], 0.171470682, fTolerance);

    //sector 2
    _info("Sector 2");

    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(2U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.89306212);
    _equal(vecXYZ[0], -2.89306212, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " <<0.483792436);
    _equal(vecXYZ[1],0.483792436, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.430943413);
    _equal(vecXYZ[2], 0.430943413, fTolerance);

	  //sector 3
		_info("Sector 3");

    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(3U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0]" << vecXYZ[0] << " vs. " << -2.78447248);
    _equal(vecXYZ[0],-2.78447248, fTolerance);
    _info("vecXYZ[1]" << vecXYZ[1] << " vs. " <<-0.00769657293);
    _equal(vecXYZ[1],-0.00769657293, fTolerance);
    _info("vecXYZ[2]" << vecXYZ[2] << " vs. " << 0.488727862);
    _equal(vecXYZ[2], 0.488727862, fTolerance);

   }

   delete feptr;
   delete fvptr;

} //end test




/**  Method:


void FiniteVolumeTraits_Test<3>::IsoparametricLinearPyramid_Test(double64 fTolerance, double64 fToleranceInternal)



Description:
Specifies the procedure of the test.

The tested element is supposed to be a uniform transformation of the original object by:
Translate(-3,0,0)
Rotate(N0, 40degrees) -3D
Shear(N0, N1, 35degrees)

This test does the following:
-construct the element.
-test the facet area, to the one measured in rhino
-test the facet normals
-test the parametric facet area
-test the volume of the sector, and check if it corresponds to the measured one
-test the sum of the volumes
-test the projection on the facet normal of a given property
-convert the facet integration points to physical space and check if they correspond to the measured values
-convert the sector integration points to physical space and check if they correspond to the measured values

@section arguments Input Arguments
double64 fTolerance - specifies the tolerance of the test
@section application Application
Runs the test for the Isoparametric Linear Pyramid element, for 3D.

tested: is a test function*/
void FiniteVolumeTraits_Test::IsoparametricLinearPyramid_Test(double64 fTolerance, double64 fToleranceInternal)
{

   FiniteElement* feptr			    = new IsoparametricLinearPyramid();
   FiniteVolumeStencil<3U>* fvptr	= new FiniteVolumeStencil<3U>("ISOPARAMETRIC_LINEAR_PYRAMID");

   Element<3U>  elmt_( feptr );
   elmt_.Assign( fvptr );

   Node<3U>  node1, node2, node3, node4, node5;
   {
	  node1.Idx( 1 );
	  node1.x( -3. );
	  node1.y( 0. );
	  node1.z( 0. );

	  node2.Idx( 2 );
	  node2.x( -1. );
	  node2.y( 1.4004151 );
	  node2.z( 0. );

	  node3.Idx( 3 );
	  node3.x( -1. );
	  node3.y( 2.9325040 );
	  node3.z( 1.2855752 );

    node4.Idx( 4 );
	  node4.x( -3. );
	  node4.y( 1.5320889 );
	  node4.z( 1.2855752 );

	  node5.Idx( 5 );
	  node5.x( -2. );
	  node5.y( 0.8234644 );
	  node5.z( 1.4088321 );
	}

	elmt_.Idx( 1 );
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );
	elmt_.Assign( 3, &node4 );
	elmt_.Assign( 4, &node5 );

  const size_t nr_of_facets = elmt_.FV_Stencil()->Facets();

   //test facet area
   //testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
   if ( m_bTestFacetAreas )
	 {
       double64 fArea = ( elmt_ ).FacetArea( 0U );
     _info("Area is: " << fArea);
	   _equal( fArea, 0.25, fTolerance );

     fArea = ( elmt_ ).FacetArea( 1U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.283693754, fTolerance );

     fArea = ( elmt_ ).FacetArea( 2U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.25, fTolerance );

     fArea = ( elmt_ ).FacetArea( 3U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.283693754, fTolerance );

	   //---
#ifndef PYRAMID_TRIANGULAR_FACETS
     fArea = ( elmt_ ).FacetArea( 4U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.3919809523, fTolerance );

     fArea = ( elmt_ ).FacetArea( 5U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.3442955048, fTolerance );

     fArea = ( elmt_ ).FacetArea( 6U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.3611982773, fTolerance );

     fArea = ( elmt_ ).FacetArea( 7U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.4420105594, fTolerance );
#else
     fArea = ( elmt_ ).FacetArea( 4U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.188595316, fTolerance );

     fArea = ( elmt_ ).FacetArea( 5U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.210174764, fTolerance );

     fArea = ( elmt_ ).FacetArea( 6U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.167969043, fTolerance );

     fArea = ( elmt_ ).FacetArea( 7U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.179194629, fTolerance );

     fArea = ( elmt_ ).FacetArea( 8U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.198434434, fTolerance );

     fArea = ( elmt_ ).FacetArea( 9U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.169625972, fTolerance );

     fArea = ( elmt_ ).FacetArea( 10U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.224321642, fTolerance );

     fArea = ( elmt_ ).FacetArea( 11U );
		 _info("Area is: " << fArea);
	   _equal( fArea, 0.220047949, fTolerance );
#endif

	   //Mapped Area
     double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.25, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.283693754, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.25, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 3U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.283693754, fTolerance );

//---
#ifndef PYRAMID_TRIANGULAR_FACETS
     fAreaMapped = ( elmt_ ).FacetAreaMapped( 4U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.3919809523, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 5U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.3442955048, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 6U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.3611982773, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 7U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.4420105594, fTolerance );
#else
     fAreaMapped = ( elmt_ ).FacetAreaMapped( 4U );
     _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.188595316, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 5U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.210174764, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 6U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.167969043, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 7U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.179194629, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 8U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.198434434, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 9U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.169625972, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 10U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.224321642, fTolerance );

     fAreaMapped = ( elmt_ ).FacetAreaMapped( 11U );
		 _info("Mapped Area is: " << fAreaMapped);
	   _equal( fAreaMapped, 0.220047949, fTolerance );
#endif

     //test mapped vs computed
    _info("Test mapped area vs computed.");
	  for(size_t i = 0; i < nr_of_facets; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearPyramid_Test<3>: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < nr_of_facets; i++)
        j += ( elmt_ ).FacetArea( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < nr_of_facets; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		 file.close();
   }


   //test facet area
   //fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
   if ( m_bTestFacetNormals )
	 {
	   Point<3U> vecNormal;

     vecNormal = ( elmt_ ).FacetNormal(0U);
		 _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 1., fTolerance );
	   _equal( vecNormal[1], 0., fTolerance );
	   _equal( vecNormal[2], 0., fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(1U);
		 _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0],  -0.4726841, fTolerance );
	   _equal( vecNormal[1], 0.6750628, fTolerance );
	   _equal( vecNormal[2], 0.5664450, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(2U);
		 _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], -1., fTolerance );
	   _equal( vecNormal[1], 0., fTolerance );
	   _equal( vecNormal[2], 0., fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(3U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.4726841, fTolerance );
	   _equal( vecNormal[1], -0.6750628, fTolerance );
	   _equal( vecNormal[2], -0.5664450, fTolerance );

//---
#ifndef PYRAMID_TRIANGULAR_FACETS
     vecNormal = ( elmt_ ).FacetNormal(4U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.4813050332, fTolerance );
	   _equal( vecNormal[1], -0.3837571776, fTolerance );
	   _equal( vecNormal[2], 0.7880836845, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(5U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.06388619905, fTolerance );
	   _equal( vecNormal[1], -0.4369081457, fTolerance );
	   _equal( vecNormal[2], 0.8972345432, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(6U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.3084014915, fTolerance );
	   _equal( vecNormal[1], -0.7699360282, fTolerance );
	   _equal( vecNormal[2], 0.5586475025, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(7U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.6290817090, fTolerance );
	   _equal( vecNormal[1], -0.6291695223, fTolerance );
	   _equal( vecNormal[2], 0.4565105864, fTolerance );
#else
     //overload tolerance for this bit
     {
     const double64 fTolerance = 1.e-6;

     vecNormal = ( elmt_ ).FacetNormal(4U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.330459, fTolerance );
	   _equal( vecNormal[1], -0.314183, fTolerance );
	   _equal( vecNormal[2], 0.889992, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(5U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.601116, fTolerance );
	   _equal( vecNormal[1], -0.433791, fTolerance );
	   _equal( vecNormal[2], 0.671182, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(6U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.00797366, fTolerance );
	   _equal( vecNormal[1], -0.542791, fTolerance );
	   _equal( vecNormal[2], 0.83983, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(7U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.115274, fTolerance );
	   _equal( vecNormal[1], -0.330666, fTolerance );
	   _equal( vecNormal[2], 0.936682, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(8U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.441986, fTolerance );
	   _equal( vecNormal[1], -0.78116, fTolerance );
	   _equal( vecNormal[2], 0.440951, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(9U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.139654, fTolerance );
	   _equal( vecNormal[1], -0.725658, fTolerance );
	   _equal( vecNormal[2], 0.673733, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(10U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.662838, fTolerance );
	   _equal( vecNormal[1], -0.548723, fTolerance );
	   _equal( vecNormal[2], 0.509459, fTolerance );

     vecNormal = ( elmt_ ).FacetNormal(11U);
     _info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
	   _equal( vecNormal[0], 0.587926, fTolerance );
	   _equal( vecNormal[1], -0.704433, fTolerance );
	   _equal( vecNormal[2], 0.39764, fTolerance );

	   } //end of overloading tolerance
#endif

	   //mapped normals
	   Point<3U> vecNormalMapped;

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);
		 _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 1., fTolerance );
	   _equal( vecNormalMapped[1], 0., fTolerance );
	   _equal( vecNormalMapped[2], 0., fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);
		 _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0],  -0.4726841, fTolerance );
	   _equal( vecNormalMapped[1], 0.6750628, fTolerance );
	   _equal( vecNormalMapped[2], 0.5664450, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);
		 _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], -1., fTolerance );
	   _equal( vecNormalMapped[1], 0., fTolerance );
	   _equal( vecNormalMapped[2], 0., fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(3U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.4726841, fTolerance );
	   _equal( vecNormalMapped[1], -0.6750628, fTolerance );
	   _equal( vecNormalMapped[2], -0.5664450, fTolerance );

//---
#ifndef PYRAMID_TRIANGULAR_FACETS
     vecNormalMapped = ( elmt_ ).FacetNormalMapped(4U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.4813050332, fTolerance );
	   _equal( vecNormalMapped[1], -0.3837571776, fTolerance );
	   _equal( vecNormalMapped[2], 0.7880836845, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(5U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.06388619905, fTolerance );
	   _equal( vecNormalMapped[1], -0.4369081457, fTolerance );
	   _equal( vecNormalMapped[2], 0.8972345432, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(6U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.3084014915, fTolerance );
	   _equal( vecNormalMapped[1], -0.7699360282, fTolerance );
	   _equal( vecNormalMapped[2], 0.5586475025, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(7U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.6290817090, fTolerance );
	   _equal( vecNormalMapped[1], -0.6291695223, fTolerance );
	   _equal( vecNormalMapped[2], 0.4565105864, fTolerance );
#else
     //overload tolerance for this bit
     {
     const double64 fTolerance = 1.e-6;

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(4U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.330459, fTolerance );
	   _equal( vecNormalMapped[1], -0.314183, fTolerance );
	   _equal( vecNormalMapped[2], 0.889992, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(5U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.601116, fTolerance );
	   _equal( vecNormalMapped[1], -0.433791, fTolerance );
	   _equal( vecNormalMapped[2], 0.671182, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(6U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.00797366, fTolerance );
	   _equal( vecNormalMapped[1], -0.542791, fTolerance );
	   _equal( vecNormalMapped[2], 0.83983, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(7U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.115274, fTolerance );
	   _equal( vecNormalMapped[1], -0.330666, fTolerance );
	   _equal( vecNormalMapped[2], 0.936682, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(8U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.441986, fTolerance );
	   _equal( vecNormalMapped[1], -0.78116, fTolerance );
	   _equal( vecNormalMapped[2], 0.440951, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(9U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.139654, fTolerance );
	   _equal( vecNormalMapped[1], -0.725658, fTolerance );
	   _equal( vecNormalMapped[2], 0.673733, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(10U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.662838, fTolerance );
	   _equal( vecNormalMapped[1], -0.548723, fTolerance );
	   _equal( vecNormalMapped[2], 0.509459, fTolerance );

     vecNormalMapped = ( elmt_ ).FacetNormalMapped(11U);
     _info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
	   _equal( vecNormalMapped[0], 0.587926, fTolerance );
	   _equal( vecNormalMapped[1], -0.704433, fTolerance );
	   _equal( vecNormalMapped[2], 0.39764, fTolerance );
	   }
#endif

	   //test mapped vs computed
		 _info("Test mapped normal vs computed., tolerance: " << fToleranceInternal);

	   for(size_t i = 0; i < nr_of_facets; i++)
	   {
			_info("Normal " << i << ":");
	    for(size_t j = 0; j < 3; j++)
       _info(( elmt_ ).FacetNormal( i )[j] << " ");
			_info(" vs ");
	    for(size_t j = 0; j < 3; j++)
       _info(( elmt_ ).FacetNormalMapped( i )[j] << " ");
	   }

	   for(size_t i = 0; i < nr_of_facets; i++)
	    for(size_t j = 0; j < 3; j++)
      _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		 //compare time computed vs mapped
		 ofstream file ("speed_compare.txt", ios::out|ios::app);
		 size_t total_times(TIMES);
		 file << "\nIsoparametricLinearPyramid_Test<3>: ";

		 clock_t ticks = clock();  double64 j(0);
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < nr_of_facets; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

	   ticks = clock();
		  for(size_t t = 0; t < total_times; t++)
		   for(size_t i = 0; i < nr_of_facets; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		 ticks = clock() - ticks;
	   file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		 file.close();
   }

   	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	  {
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
        double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
				_info("Area is: " << fArea);
		    _equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
          }
	  }

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	  {
		Point<3U> vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		  {
		    //ignore return parameter, jacobian
        vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);
            dumpVector<3>("Normal is: ", vecNormal);
		    _equal( vecNormal[0],  elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
		    _equal( vecNormal[1],  elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
		    _equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		  }
	  }

   //test sector volumes
   //testing: fT   SectorVolume( size_t sector ) const;
   if ( m_bTestSectorVolumes )
	 {
	   double64 fVolSum(0.);

	   double64 fSectorVolume(0.);
	   for ( size_t iSector = 0U; iSector < elmt_.FV_Stencil()->Sectors()-1; iSector++ )
	   {
       fSectorVolume = ( elmt_ ).SectorVolume( iSector );
				 //_info("Volume is: " << fSectorVolume);
	       _equal( fSectorVolume, 1./4., fTolerance );
	       fVolSum+=fSectorVolume;
	   }

     fSectorVolume = ( elmt_ ).SectorVolume( 4U );
		 _info("Volume is: " << fSectorVolume);
	   _equal( fSectorVolume, 1./3., fTolerance );
	   fVolSum+=fSectorVolume;

		 _info("Volume of the Element: " << elmt_.Volume() << " vs Volume Sum: " << fVolSum);
	   _equal(fVolSum, elmt_.Volume(), fTolerance);
     }

   if(m_bProjectionOnFacetNormal)
   {
	   const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
	   for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
	   {
	     //const CSPINDEX& prop_key
	     VectorVariable<3> vVariable;
	     for(size_t iD= 0U; iD < 3; iD++) {
	     	vVariable.Flag(iD)=PLAIN;
	     	vVariable(iD)=3.;
	     }

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct(( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
		_equal(fProjectionVal, fProjectionValP, fTolerance);
	   }
   }

   if(m_bRSTToXYZ)
   {
   	elmt_.CoordinateMatrix();

   	Point<3U> vecRST, vecXYZ;

	  //check facet integration points
   	//facet 0
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(0U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2);
    _equal(vecXYZ[0],-2, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.03651173);
    _equal(vecXYZ[1], 1.03651173, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.463499289);
    _equal(vecXYZ[2], 0.463499289, fTolerance);

    //facet 1
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(1U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -1.55555556);
    _equal(vecXYZ[0], -1.55555556, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " <<1.68817927);
    _equal(vecXYZ[1],1.68817927, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.749182671);
    _equal(vecXYZ[2], 0.749182671, fTolerance);

    //facet 2
   	vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(2U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2);
    _equal(vecXYZ[0], -2, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.71744012);
    _equal(vecXYZ[1],1.71744012, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 1.03486605);
    _equal(vecXYZ[2], 1.03486605, fTolerance);

    //facet 3
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(3U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.44444444 );
    _equal(vecXYZ[0], -2.44444444, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.06577257);
    _equal(vecXYZ[1],1.06577257, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.749182671);
    _equal(vecXYZ[2], 0.749182671, fTolerance);

    //facet 4 - equivalent to the facet integration point for the non-planar facets.
#ifndef PYRAMID_TRIANGULAR_FACETS
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(4U, 0U);
#else
    vecRST = Point< 3>(-7./24.,-7./24.,17./48.);
#endif

    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.291666675);
    _equal(vecXYZ[0], -2.291666675, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 0.8109412000);
    _equal(vecXYZ[1], 0.8109412000, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.7266153000);
    _equal(vecXYZ[2], 0.7266153000, fTolerance);

    //facet 5
#ifndef PYRAMID_TRIANGULAR_FACETS
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(5U, 0U);
#else
    vecRST = Point< 3>(7./24.,-7./24.,17./48.);
#endif

    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -1.708333325);
    _equal(vecXYZ[0], -1.708333325, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.219395600);
    _equal(vecXYZ[1], 1.219395600, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.7266153000);
    _equal(vecXYZ[2], 0.7266153000, fTolerance);

    //facet 6
#ifndef PYRAMID_TRIANGULAR_FACETS
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(6U, 0U);
#else
    vecRST = Point< 3>(7./24.,7./24.,17./48.);
#endif

    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -1.708333325);
    _equal(vecXYZ[0], -1.708333325, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.666254875);
    _equal(vecXYZ[1], 1.666254875, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 1.101574725);
    _equal(vecXYZ[2], 1.101574725, fTolerance);

    //facet 7
#ifndef PYRAMID_TRIANGULAR_FACETS
    vecRST = elmt_.FV_Stencil()->FacetIntegrationPoint(7U, 0U);
#else
    vecRST = Point< 3>(-7./24.,7./24.,17./48.);
#endif

    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.291666675);
    _equal(vecXYZ[0], -2.291666675, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.257800450);
    _equal(vecXYZ[1], 1.257800450, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 1.101574725);
    _equal(vecXYZ[2], 1.101574725, fTolerance);

    //check sector integration points
    //sector 0
    _info("Sector 0");
	  vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(0U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.43055556);
    _equal(vecXYZ[0], -2.43055556, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 0.732281579);
    _equal(vecXYZ[1],0.732281579, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.488386154);
    _equal(vecXYZ[2], 0.488386154, fTolerance);

    //sector 1
    _info("Sector 1");
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(1U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -1.56944444);
    _equal(vecXYZ[0], -1.56944444, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.33523807);
    _equal(vecXYZ[1], 1.33523807, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 0.488386154);
    _equal(vecXYZ[2], 0.488386154, fTolerance);

    //sector 2
    _info("Sector 2");
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(2U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " <<  -1.56944444);
    _equal(vecXYZ[0], -1.56944444, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.99488745);
    _equal(vecXYZ[1], 1.99488745, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 1.04189771);
    _equal(vecXYZ[2], 1.04189771, fTolerance);

	  //sector 3
		_info("Sector 3");
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(3U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);
    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2.43055556);
    _equal(vecXYZ[0], -2.43055556, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.39193096);
    _equal(vecXYZ[1], 1.39193096, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 1.04189771);
    _equal(vecXYZ[2], 1.04189771, fTolerance);

  	//sector 4
    _info("Sector 4");
    vecRST = elmt_.FV_Stencil()->SectorIntegrationPoint(4U, 0U);
    ( elmt_ ).N_At(vecRST);
    vecXYZ = ( elmt_ ).RstToXYZ(vecRST);

    _info("vecXYZ[0] " << vecXYZ[0] << " vs. " << -2);
    _equal(vecXYZ[0],-2, fTolerance);
    _info("vecXYZ[1] " << vecXYZ[1] << " vs. " << 1.13146677);
    _equal(vecXYZ[1], 1.13146677, fTolerance);
    _info("vecXYZ[2] " << vecXYZ[2] << " vs. " << 1.04176909);
    _equal(vecXYZ[2], 1.04176909, fTolerance);
   }

   delete feptr;
   delete fvptr;

} //end test




/**  Method:


void FiniteVolumeTraits_Test<3>::IsoparametricLinearHexahedron_Test(double64 fTolerance, double64 fToleranceInternal)



Description:
Specifies the procedure of the test.

The tested element is supposed to be a uniform transformation of the original object by:
Translate(-3,0,0)
Rotate(N0, 40degrees) -3D
Shear(N0, N1, 35degrees)

This test does the following:
-construct the element.
-test the facet area, to the one measured in rhino
-test the facet normals
-test the parametric facet area
-test the volume of the sector, and check if it corresponds to the measured one
-test the sum of the volumes
-test the projection on the facet normal of a given property
-convert the facet integration points to physical space and check if they correspond to the measured values
-convert the sector integration points to physical space and check if they correspond to the measured values

@section arguments Input Arguments
double64 fTolerance - specifies the tolerance of the test
@section application Application
Runs the test for the Isoparametric Linear Hexahedron element, for 3D.

@note This set of tests uses a unit iso parametric hexahedron.  That means that it is a Hexahedron with sides 2a*2b*2c
where a,b, and c are the sides in each coordinate direction.  The volumetric center of this hexahedron is the origin
0,0,0.

 */
void FiniteVolumeTraits_Test::Test_UnitaryIsoparametricLinearHexahedron(double64 fTolerance, double64 fToleranceInternal)
{
	 FiniteElement* feptr = new IsoparametricLinearHexahedron();
   FiniteVolumeStencil<3U>* fvptr	= new FiniteVolumeStencil<3>("ISOPARAMETRIC_LINEAR_HEXAHEDRON");

   Element<3U>  elmt_( feptr );
   elmt_.Assign( fvptr );

	Node<3U>  node1, node2, node3, node4, node5, node6, node7, node8;
    {
      node1.Idx( 1 );
      node1.x( -1. );
      node1.y( -1. );
      node1.z( -1. );

      node2.Idx( 2 );
      node2.x(1. );
      node2.y(-1. );
      node2.z(-1. );

      node3.Idx( 3 );
      node3.x(1. );
      node3.y(1. );
      node3.z(-1. );

      node4.Idx( 4 );
      node4.x(-1. );
      node4.y(1. );
      node4.z(-1. );

      node5.Idx( 5 );
      node5.x( -1. );
      node5.y( -1. );
      node5.z( 1. );

      node6.Idx( 6 );
      node6.x(1. );
      node6.y(-1. );
      node6.z(1. );

      node7.Idx( 7 );
      node7.x(1. );
      node7.y(1. );
      node7.z(1. );

      node8.Idx( 8 );
      node8.x(-1. );
      node8.y(1. );
      node8.z(1. );
    }

	elmt_.Idx( 1 );
	//elmt_.FE()->CurrentID(0);
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );
	elmt_.Assign( 3, &node4 );
	elmt_.Assign( 4, &node5 );
	elmt_.Assign( 5, &node6 );
	elmt_.Assign( 6, &node7 );
	elmt_.Assign( 7, &node8 );

	//test facet area
	//testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestFacetAreas )
	{
		_info("Facet area test.  Method 1 (facetarea1 in FiniteVolumeTraits) ");
    double64 fArea = ( elmt_ ).FacetArea( 0U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 1U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 2U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 3U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 4U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 5U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 6U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 7U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 8U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 9U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 10U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 11U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

		//area mapped
		_info("Facet area test.  Method 2 (FaceAreaMapped in FiniteVolumeTraits) ");
    double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 3U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 4U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 5U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 6U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 7U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 8U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 9U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 10U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 11U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );


		//test mapped vs computed
		_info("Test mapped area vs computed.");
		for(size_t i = 0; i < 12; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearHexahedron_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetArea( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		file.close();
	}	//test facet area

	//fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
	if ( m_bTestFacetNormals )
	{
		Point<3U> vecNormal;
		/*ignore return parameter, jacobian*/
		_info(" Unit Hexahedron.  Facet Normal calculations (using FacetNormal...");
		_info("Facet 1...");
    vecNormal = ( elmt_ ).FacetNormal(0U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 2...");
    vecNormal = ( elmt_ ).FacetNormal(1U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 3...");
    vecNormal = ( elmt_ ).FacetNormal(2U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], -1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 4...");
    vecNormal = ( elmt_ ).FacetNormal(3U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], -1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 5...");
    vecNormal = ( elmt_ ).FacetNormal(4U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 6...");
    vecNormal = ( elmt_ ).FacetNormal(5U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 7...");
    vecNormal = ( elmt_ ).FacetNormal(6U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 8...");
    vecNormal = ( elmt_ ).FacetNormal(7U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 9...");
    vecNormal = ( elmt_ ).FacetNormal(8U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 10...");
    vecNormal = ( elmt_ ).FacetNormal(9U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 11...");
    vecNormal = ( elmt_ ).FacetNormal(10U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], -1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 12...");
    vecNormal = ( elmt_ ).FacetNormal(11U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], -1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );

		//mapped normal
		Point<3U> vecNormalMapped;
		_info("Unit Hexahedron.  Facet Normal calculations (using FacetNormalMapped...");
		_info("Facet 1...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 2...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );

		_info("Facet 3...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], -1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 4...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(3U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], -1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 5...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(4U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 6...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(5U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 7...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(6U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 8...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(7U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 9...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(8U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 10...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(9U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 11...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(10U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], -1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 12...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(11U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], -1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );


		//test mapped vs computed
		_info("Test mapped normal vs computed., tolerance: " << fToleranceInternal);
		for(size_t i = 0; i < 8; i++)
		{
			_info("Normal " << i << ":");
			for(size_t j = 0; j < 3; j++)
        _info(( elmt_ ).FacetNormal( i )[j] << " ");
			_info(" vs ");
			for(size_t j = 0; j < 3; j++)
        _info(( elmt_ ).FacetNormalMapped( i )[j] << " ");
		}

		for(size_t i = 0; i < 12; i++)
			for(size_t j = 0; j < 3; j++)
        _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearHexahedron_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		file.close();
	}

	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	{
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
      double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
			_info("Area is: " << fArea);
			_equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
		}
	}

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	{
		Point<3U> vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
			//ignore return parameter, jacobian
      vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);
            dumpVector<3>("Normal is: ", vecNormal);
			_equal( vecNormal[0], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
			_equal( vecNormal[1], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
			_equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		}
	}

	//test sector volumes
	//testing: fT   SectorVolume( size_t sector ) const;
	if ( m_bTestSectorVolumes )
	{
		double64 fVolSum(0.);
		double64 fSectorVolume(0.);
		for ( size_t iSector = 0U; iSector < elmt_.FV_Stencil()->Sectors(); iSector++ )
		{
      fSectorVolume = ( elmt_ ).SectorVolume( iSector );
			//_info("Volume is: " << fSectorVolume);
			_equal( fSectorVolume, 1., fTolerance );
			fVolSum+=fSectorVolume;
		}

		_info("Volume of the Element: " << elmt_.Volume() << " vs Volume Sum: " << fVolSum);
		_equal(fVolSum, elmt_.Volume(), fTolerance);
	}


	if(m_bProjectionOnFacetNormal)
	{
		const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
		for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
		{
			//const CSPINDEX& prop_key
			VectorVariable<3> vVariable;
			for(size_t iD= 0U; iD < 3; iD++) {
				vVariable.Flag(iD)=PLAIN;
				vVariable(iD)=3.;
			}

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct( ( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
			_equal(fProjectionVal, fProjectionValP, fTolerance);
		}

	}

	delete feptr;
	delete fvptr;

} //end test


// This set of tests uses a scaled-up unit iso parametric hexahedron.  That means that it is a Hexahedron with sides 2a*2b*2c
// where a,b, and c are the sides in each coordinate direction.  The volumetric center of this hexahedron is the origin
// 0,0,0.
void FiniteVolumeTraits_Test::Test_IsoparametricLinearHexahedron1(double64 fTolerance, double64 fToleranceInternal)
{
	FiniteElement* feptr			= new IsoparametricLinearHexahedron();
    FiniteVolumeStencil<3U>* fvptr	= new FiniteVolumeStencil<3>("ISOPARAMETRIC_LINEAR_HEXAHEDRON");

    Element<3U>  elmt_( feptr );
    elmt_.Assign( fvptr );

	Node<3U>  node1, node2, node3, node4, node5, node6, node7, node8;
	{
		node1.Idx( 1 );
		node1.x( -2. );
		node1.y( -2. );
		node1.z( -2. );

		node2.Idx( 2 );
		node2.x(2. );
		node2.y(-2. );
		node2.z(-2. );

		node3.Idx( 3 );
		node3.x(2. );
		node3.y(2. );
		node3.z(-2. );

		node4.Idx( 4 );
		node4.x(-2. );
		node4.y(2. );
		node4.z(-2. );

		node5.Idx( 5 );
		node5.x( -2. );
		node5.y( -2. );
		node5.z( 2. );

		node6.Idx( 6 );
		node6.x(2. );
		node6.y(-2. );
		node6.z(2. );

		node7.Idx( 7 );
		node7.x(2. );
		node7.y(2. );
		node7.z(2. );

		node8.Idx( 8 );
		node8.x(-2. );
		node8.y(2. );
		node8.z(2. );
	}

	elmt_.Idx( 1 );
	//elmt_.FE()->CurrentID(0);
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );
	elmt_.Assign( 3, &node4 );
	elmt_.Assign( 4, &node5 );
	elmt_.Assign( 5, &node6 );
	elmt_.Assign( 6, &node7 );
	elmt_.Assign( 7, &node8 );

	//test facet area
	//testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestFacetAreas )
	{
		_info("Facet area test.  Method 1 (facetarea1 in FiniteVolumeTraits) ");
    double64 fArea = ( elmt_ ).FacetArea( 0U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 1U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 2U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 3U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 4U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 5U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 6U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 7U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 8U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 9U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 10U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

    fArea = ( elmt_ ).FacetArea( 11U );
		_info("Area is: " << fArea);
		_equal( fArea, 4., fTolerance );

		//area mapped
		_info("Facet area test.  Method 2 (FaceAreaMapped in FiniteVolumeTraits) ");
    double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 3U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 4U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 5U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 6U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 7U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 8U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 9U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 10U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 11U );
		_info("Mapped Facet Area is: " << fAreaMapped);
		_equal( fAreaMapped, 4., fTolerance );


		//test mapped vs computed
		_info("Test mapped area vs computed.");
		for(size_t i = 0; i < 12; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearHexahedron_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetArea( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		file.close();
	}	//test facet area

	//fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
	if ( m_bTestFacetNormals )
	{
		Point<3U> vecNormal;
		/*ignore return parameter, jacobian*/
		_info("Unit Hexahedron.  Facet Normal calculations (using FacetNormal...");
		_info("Facet 1...");
    vecNormal = ( elmt_ ).FacetNormal(0U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 2...");
    vecNormal = ( elmt_ ).FacetNormal(1U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 3...");
    vecNormal = ( elmt_ ).FacetNormal(2U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], -1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 4...");
    vecNormal = ( elmt_ ).FacetNormal(3U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], -1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 5...");
    vecNormal = ( elmt_ ).FacetNormal(4U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 6...");
    vecNormal = ( elmt_ ).FacetNormal(5U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 7...");
    vecNormal = ( elmt_ ).FacetNormal(6U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 8...");
    vecNormal = ( elmt_ ).FacetNormal(7U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 1., fTolerance );
		_info("Facet 9...");
    vecNormal = ( elmt_ ).FacetNormal(8U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 10...");
    vecNormal = ( elmt_ ).FacetNormal(9U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], 1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 11...");
    vecNormal = ( elmt_ ).FacetNormal(10U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], -1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );
		_info("Facet 12...");
    vecNormal = ( elmt_ ).FacetNormal(11U);
        dumpVector<3>("Normal is: ", vecNormal);
		_equal( vecNormal[0], 0., fTolerance );
		_equal( vecNormal[1], -1., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );

		//mapped normal
		Point<3U> vecNormalMapped;
		_info("Unit Hexahedron.  Facet Normal calculations (using FacetNormalMapped...");
		_info("Facet 1...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 2...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 3...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], -1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 4...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(3U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], -1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 5...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(4U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 6...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(5U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 7...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(6U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 8...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(7U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 1., fTolerance );
		_info("Facet 9...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(8U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 10...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(9U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], 1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 11...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(10U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], -1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );
		_info("Facet 12...");
    vecNormalMapped = ( elmt_ ).FacetNormalMapped(11U);
        dumpVector<3>("Normal is: ", vecNormalMapped);
		_equal( vecNormalMapped[0], 0., fTolerance );
		_equal( vecNormalMapped[1], -1., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );


		//test mapped vs computed
		_info("Test mapped normal vs computed., tolerance: " << fToleranceInternal);
		for(size_t i = 0; i < 8; i++)
		{
			_info("Normal " << i << ":");
            dumpVector<3>(" ", elmt_.FacetNormal( i ));
            dumpVector<3>(" ", elmt_.FacetNormalMapped( i ));
		}

		for(size_t i = 0; i < 12; i++)
			for(size_t j = 0; j < 3; j++)
        _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearHexahedron_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		file.close();
	}

	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	{
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
      double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
			_info("Area is: " << fArea);
			_equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
		}
	}

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	{
		Point<3U> vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
			//ignore return parameter, jacobian
      vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);
            dumpVector<3>("Normal is: ", vecNormal);
			_equal( vecNormal[0], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
			_equal( vecNormal[1], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
			_equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		}
	}

	//test sector volumes
	//testing: fT   SectorVolume( size_t sector ) const;
	if ( m_bTestSectorVolumes )
	{
		double64 fVolSum(0.);
		double64 fSectorVolume(0.);
		for ( size_t iSector = 0U; iSector < elmt_.FV_Stencil()->Sectors(); iSector++ )
		{
      fSectorVolume = ( elmt_ ).SectorVolume( iSector );
			_info("Volume is: " << fSectorVolume);
			_equal( fSectorVolume, 8., fTolerance );
			fVolSum+=fSectorVolume;
		}

		_info("Volume of the Element: " << elmt_.Volume());
		_info("vs Volume Sum: " << fVolSum);
		_equal(fVolSum, elmt_.Volume(), fTolerance);
	}


	if(m_bProjectionOnFacetNormal)
	{
		//In previous versions the following part of the test used to be incorrect (or was not a real test!):
		//
		//which actually projects a vector variable (vVariable) native of physical space onto a
		//facet normal in parametric space.
		//Result from this had little sense (i.e. comparing the projections), unless the hexahedron in
		//physical space has a shape equal to the onein parametric space (can be scaled,but maintains the same normals)
		//
		//A valid test would be the integration of the vector field, using physical coordinates,
		//and comparing them to the integration in parametric space (they should yield the same value.
		//This could be done using a hexahedron without a Det(J)=1, like in this case.
		//Given that the element is linear, the Jacobian is practically constant across the whole element,
		//and hence this is simply an exchange of an area calculation for a jacobian determinant calculation.
		//in essence, both of these are the same, but there is not any immediate functionality that needs to be tested
		//here. (or compared).
		const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
		for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
		{
			//const CSPINDEX& prop_key
			VectorVariable<3> vVariable;
			for(size_t iD= 0U; iD < 3; iD++) {
				vVariable.Flag(iD)=PLAIN;
				vVariable(iD)=3.;
			}

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct( ( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
			_equal(fProjectionVal, fProjectionValP, fTolerance);
		}

	}

	delete feptr;
	delete fvptr;

} //end test


// This set of tests uses a random iso parametric hexahedron.  That means that it is a Hexahedron with geometry given by
// the corner points defined below.  The volumetric center of the element is NOT the origin.
void FiniteVolumeTraits_Test::Test_IsoparametricLinearHexahedron2(double64 fTolerance, double64 fToleranceInternal)
{
	FiniteElement* feptr			= new IsoparametricLinearHexahedron();
    FiniteVolumeStencil<3U>* fvptr	= new FiniteVolumeStencil<3>("ISOPARAMETRIC_LINEAR_HEXAHEDRON");

    Element<3U>  elmt_( feptr );
    elmt_.Assign( fvptr );

	Node<3U>  node1, node2, node3, node4, node5, node6, node7, node8;
	{
		node1.Idx( 1 );
		node1.x( -3. );
		node1.y( 0. );
		node1.z( 0. );

		node2.Idx( 2 );
		node2.x(-1.0000000 );
		node2.y( 1.4004151 );
		node2.z( 0. );

		node3.Idx( 3 );
		node3.x( -1.0000000 );
		node3.y( 2.9325040 );
		node3.z( 1.2855752 );

		node4.Idx( 4 );
		node4.x( -3.0000000 );
		node4.y( 1.5320889 );
		node4.z( 1.2855752 );

		node5.Idx( 5 );
		node5.x( -3.0000000 );
		node5.y( -1.2855752 );
		node5.z( 1.5320889 );

		node6.Idx( 6 );
		node6.x( -1.0000000 );
		node6.y( 0.1148399 );
		node6.z( 1.5320889 );

		node7.Idx( 7 );
		node7.x( -1.0000000 );
		node7.y( 1.6469287 );
		node7.z( 2.8176641 );

		node8.Idx( 8 );
		node8.x( -3.0000000 );
		node8.y( 0.2465137 );
		node8.z( 2.8176641 );
	}

	elmt_.Idx( 1 );
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );
	elmt_.Assign( 3, &node4 );
	elmt_.Assign( 4, &node5 );
	elmt_.Assign( 5, &node6 );
	elmt_.Assign( 6, &node7 );
	elmt_.Assign( 7, &node8 );

	//test facet area
	//testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestFacetAreas )
	{
    double64 fArea = ( elmt_ ).FacetArea( 0U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 1U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.13477501, fTolerance );

    fArea = ( elmt_ ).FacetArea( 2U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 3U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.13477501, fTolerance );

    fArea = ( elmt_ ).FacetArea( 4U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.09662038, fTolerance );

    fArea = ( elmt_ ).FacetArea( 5U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.09662038, fTolerance );

    fArea = ( elmt_ ).FacetArea( 6U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.09662038, fTolerance );

    fArea = ( elmt_ ).FacetArea( 7U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.09662038, fTolerance );

    fArea = ( elmt_ ).FacetArea( 8U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 9U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.13477501, fTolerance );

    fArea = ( elmt_ ).FacetArea( 10U );
		_info("Area is: " << fArea);
		_equal( fArea, 1., fTolerance );

    fArea = ( elmt_ ).FacetArea( 11U );
		_info("Area is: " << fArea);
		_equal( fArea, 1.13477501, fTolerance );

		//area mapped
    double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.13477501, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 3U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.13477501, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 4U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.09662038, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 5U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.09662038, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 6U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.09662038, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 7U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.09662038, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 8U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 9U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.13477501, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 10U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1., fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 11U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 1.13477501, fTolerance );

		//test mapped vs computed
		_info("Test mapped area vs computed.");
		for(size_t i = 0; i < 12; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearHexahedron_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetArea( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		file.close();
	}	//test facet area

	//fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
	if ( m_bTestFacetNormals )
	{
		Point<3U> vecNormal;
		/*ignore return parameter, jacobian*/

    vecNormal = ( elmt_ ).FacetNormal(0U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], 1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(1U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], -0.4726841, fTolerance );
		_equal( vecNormal[1], 0.6750628, fTolerance );
		_equal( vecNormal[2], 0.5664450, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(2U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], -1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(3U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4726841, fTolerance );
		_equal( vecNormal[1], -0.6750628, fTolerance );
		_equal( vecNormal[2], -0.5664450, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(4U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4104289, fTolerance );
		_equal( vecNormal[1], -0.5861533, fTolerance );
		_equal( vecNormal[2],  0.6985503, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(5U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4104289, fTolerance );
		_equal( vecNormal[1], -0.5861533, fTolerance );
		_equal( vecNormal[2],  0.6985503, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(6U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],   0.4104289, fTolerance );
		_equal( vecNormal[1],  -0.5861533, fTolerance );
		_equal( vecNormal[2],   0.6985503, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(7U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4104289, fTolerance );
		_equal( vecNormal[1], -0.5861533, fTolerance );
		_equal( vecNormal[2],  0.6985503, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(8U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  1., fTolerance );
		_equal( vecNormal[1],  0., fTolerance );
		_equal( vecNormal[2],  0., fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(9U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], -0.4726841, fTolerance );
		_equal( vecNormal[1], 0.6750628, fTolerance );
		_equal( vecNormal[2], 0.5664450, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(10U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], -1., fTolerance );
		_equal( vecNormal[1], 0., fTolerance );
		_equal( vecNormal[2], 0., fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(11U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4726841, fTolerance );
		_equal( vecNormal[1], -0.6750628, fTolerance );
		_equal( vecNormal[2], -0.5664450, fTolerance );

		//mapped normal
		Point<3U> vecNormalMapped;

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], 1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], -0.4726841, fTolerance );
		_equal( vecNormalMapped[1], 0.6750628, fTolerance );
		_equal( vecNormalMapped[2], 0.5664450, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], -1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(3U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4726841, fTolerance );
		_equal( vecNormalMapped[1], -0.6750628, fTolerance );
		_equal( vecNormalMapped[2], -0.5664450, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(4U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4104289, fTolerance );
		_equal( vecNormalMapped[1], -0.5861533, fTolerance );
		_equal( vecNormalMapped[2],  0.6985503, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(5U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4104289, fTolerance );
		_equal( vecNormalMapped[1], -0.5861533, fTolerance );
		_equal( vecNormalMapped[2],  0.6985503, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(6U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],   0.4104289, fTolerance );
		_equal( vecNormalMapped[1],  -0.5861533, fTolerance );
		_equal( vecNormalMapped[2],   0.6985503, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(7U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4104289, fTolerance );
		_equal( vecNormalMapped[1], -0.5861533, fTolerance );
		_equal( vecNormalMapped[2],  0.6985503, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(8U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  1., fTolerance );
		_equal( vecNormalMapped[1],  0., fTolerance );
		_equal( vecNormalMapped[2],  0., fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(9U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], -0.4726841, fTolerance );
		_equal( vecNormalMapped[1], 0.6750628, fTolerance );
		_equal( vecNormalMapped[2], 0.5664450, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(10U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], -1., fTolerance );
		_equal( vecNormalMapped[1], 0., fTolerance );
		_equal( vecNormalMapped[2], 0., fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(11U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4726841, fTolerance );
		_equal( vecNormalMapped[1], -0.6750628, fTolerance );
		_equal( vecNormalMapped[2], -0.5664450, fTolerance );

		//test mapped vs computed
		_info("Test mapped normal vs computed., tolerance: " << fToleranceInternal);

		for(size_t i = 0; i < 8; i++)
		{
			_info("Normal " << i << ":");
			for(size_t j = 0; j < 3; j++)
        _info(( elmt_ ).FacetNormal( i )[j] << " ");
			_info(" vs ");
			for(size_t j = 0; j < 3; j++)
        _info(( elmt_ ).FacetNormalMapped( i )[j] << " ");
		}

		for(size_t i = 0; i < 12; i++)
			for(size_t j = 0; j < 3; j++)
        _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearHexahedron_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 12; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		file.close();
	}

	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	{
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
      double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
			_info("Area is: " << fArea);
			_equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
		}
	}

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	{
		Point<3U> vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
			//ignore return parameter, jacobian
      vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);
			_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
			_equal( vecNormal[0], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
			_equal( vecNormal[1], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
			_equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		}
	}

	//test sector volumes
	//testing: fT   SectorVolume( size_t sector ) const;
	if ( m_bTestSectorVolumes )
	{
		double64 fVolSum(0.);
		double64 fSectorVolume(0.);
		for ( size_t iSector = 0U; iSector < elmt_.FV_Stencil()->Sectors(); iSector++ )
		{
      fSectorVolume = ( elmt_ ).SectorVolume( iSector );
			//_info("Volume is: " << fSectorVolume);
			_equal( fSectorVolume, 1., fTolerance );
			fVolSum+=fSectorVolume;
		}

		_info("Volume of the Element: " << elmt_.Volume() << " vs Volume Sum: " << fVolSum);
		_equal(fVolSum, elmt_.Volume(), fTolerance);
	}


	if(m_bProjectionOnFacetNormal)
	{
		const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
		//double64 int_physicalspace(0.);
		//double64 int_parametricspace(0.);
		//double64 j_factor(0.);

		for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
		{
			//const CSPINDEX& prop_key
			VectorVariable<3U> vVariable;
			VectorVariable<3U> vVariableP;
			VectorVariable<3U> pFacetNormal;
			for(size_t iD= 0U; iD < 3; iD++) {
				vVariable.Flag(iD)=PLAIN;
				vVariableP.Flag(iD)=PLAIN;
				vVariable(iD)=3.;
				pFacetNormal.Flag(iD)=PLAIN;
        pFacetNormal(iD)=(( elmt_ ).ParametricFacetNormal( iFacet ))[iD];
			}

			//This calculates the projection (dot product) of vVariable (3,3,3) onto the facet normal (also in physical space)
      //double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      //double64 fProjectionValP = vVariable.DotProduct(( elmt_ ).ParametricFacetNormal( iFacet ));

			//In previous versions the following part of the test used to be incorrect (or was not a real test!):
			//
			//which actually projects a vector variable (vVariable) native of physical space onto a
			//facet normal in parametric space.
			//Result from this had little sense (i.e. comparing the projections), unless the hexahedron in
			//physical space has a shape equal to the onein parametric space (can be scaled,but maintains the same normals)
			//
			//A valid test would be the integration of the vector field, using physical coordinates,
			//and comparing them to the integration in parametric space (they should yield the same value.
			//This could be done using a hexahedron without a Det(J)=1, like in this case.
			//Given that the element is linear, the Jacobian is practically constant across the whole element,
			//and hence this is simply an exchange of an area calculation for a jacobian determinant calculation.
			//in essence, both of these are the same, but there is not any immediate functionality that needs to be tested
			//here. (or compared).
			/*
			j_factor=5.741612292;
      int_physicalspace=fProjectionVal*( elmt_ ).FacetArea( iFacet );
      int_parametricspace=fProjectionValP*( elmt_ ).ParametricFacetArea( iFacet );
			cout<<"Facet N. : "<<iFacet<<" Int. PhysicSpace= "<<int_physicalspace<<endl;
			cout<<"Facet N. : "<<iFacet<<" Int. ParameSpace= "<<int_parametricspace<<endl;
			//

			//create a point to use to pass to RstToXYZ
			Point<3U> p(vVariable[0U],vVariable[1U],vVariable[2U]);
      //double64	fProjectionValP =pFacetNormal.DotProduct(( elmt_ ).RstToXYZ(p));

      Point<3U> pPoint((( elmt_ ).RstToXYZ(p))[0U],(( elmt_ ).RstToXYZ(p))[1U],(( elmt_ ).RstToXYZ(p))[2U]);

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
			_info("Projected vector Length in physical space: "<<fProjectionVal<<endl);
			vVariable(0U)=pPoint[0U];
			vVariable(1U)=pPoint[1U];
			vVariable(2U)=pPoint[2U];
			_info("Projected vector Length in parametric space: "<<vVariable.Length()<<endl);
			_equal(fProjectionVal, fProjectionValP, fTolerance);
			*/
		}

	}

	delete feptr;
	delete fvptr;

} //end test

/**  Method:










void FiniteVolumeTraits_Test<3>::IsoparametricLinearPrism_Test(double64 fTolerance)

Description:
Specifies the procedure of the test.

The tested element is supposed to be a uniform transformation of the original object by:
Translate(-3,0,0)
Rotate(N0, 40degrees) -3D
Shear(N0, N1, 35degrees)

This test does the following:
-construct the element.
-test the facet area, to the one measured in rhino
-test the facet normals
-test the parametric facet area
-test the volume of the sector, and check if it corresponds to the measured one
-test the sum of the volumes
-test the projection on the facet normal of a given property
-convert the facet integration points to physical space and check if they correspond to the measured values
-convert the sector integration points to physical space and check if they correspond to the measured values

@section arguments Input Arguments
double64 fTolerance - specifies the tolerance of the test
@section application Application
Runs the test for the Isoparametric Linear Prism element, for 3D.

tested: is a test function*/
void FiniteVolumeTraits_Test::IsoparametricLinearPrism_Test(double64 fTolerance, double64 fToleranceInternal)
{
	FiniteElement* feptr			= new IsoparametricLinearPrism();
    FiniteVolumeStencil<3U>* fvptr	= new FiniteVolumeStencil<3>("ISOPARAMETRIC_LINEAR_PRISM");

    Element<3U>  elmt_( feptr );
    elmt_.Assign( fvptr );

	Node< 3>  node1, node2, node3, node4, node5, node6;
	{
	  node1.Idx( 1 );
	  node1.x( -3. );
	  node1.y( 0. );
	  node1.z( 0. );

	  node2.Idx( 2 );
	  node2.x( -2.0000000 );
	  node2.y( 0.7002075 );
	  node2.z( 0. );

	  node3.Idx( 3 );
	  node3.x(  -3. );
	  node3.y( 0.7660444 );
	  node3.z( 0.6427876 );

    node4.Idx( 4 );
	  node4.x( -3. );
	  node4.y( -1.2855752 );
	  node4.z( 1.5320889 );

	  node5.Idx( 5 );
	  node5.x(-2. );
	  node5.y( -0.5853677 );
	  node5.z( 1.5320889 );

	  node6.Idx( 6 );
	  node6.x( -3. );
	  node6.y( -0.5195308 );
	  node6.z( 2.1748765 );

	}

	elmt_.Idx( 1 );
	elmt_.Assign( 0, &node1 );
	elmt_.Assign( 1, &node2 );
	elmt_.Assign( 2, &node3 );
	elmt_.Assign( 3, &node4 );
	elmt_.Assign( 4, &node5 );
	elmt_.Assign( 5, &node6 );

	//test facet area
	//testing: fT   FacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestFacetAreas )
	{
    double64 fArea = ( elmt_ ).FacetArea( 0U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.295435364, fTolerance );

    fArea = ( elmt_ ).FacetArea( 1U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.305527528, fTolerance );

    fArea = ( elmt_ ).FacetArea( 2U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.333553967, fTolerance );

    fArea = ( elmt_ ).FacetArea( 3U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.182770064, fTolerance );

    fArea = ( elmt_ ).FacetArea( 4U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.182770064, fTolerance );

    fArea = ( elmt_ ).FacetArea( 5U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.182770064, fTolerance );

    fArea = ( elmt_ ).FacetArea( 6U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.295435364, fTolerance );

    fArea = ( elmt_ ).FacetArea( 7U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.305527528, fTolerance );

    fArea = ( elmt_ ).FacetArea( 8U );
		_info("Area is: " << fArea);
		_equal( fArea, 0.333553967, fTolerance );

		//area mapped
    double64 fAreaMapped = ( elmt_ ).FacetAreaMapped( 0U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.295435364, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 1U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.305527528, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 2U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.333553967, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 3U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.182770064, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 4U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.182770064, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 5U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.182770064, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 6U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.295435364, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 7U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.305527528, fTolerance );

    fAreaMapped = ( elmt_ ).FacetAreaMapped( 8U );
		_info("Mapped Area is: " << fAreaMapped);
		_equal( fAreaMapped, 0.333553967, fTolerance );

	  //test mapped vs computed
		_info("Test mapped area vs computed.");
	  for(size_t i = 0; i < 9; i++)
      _equal( ( elmt_ ).FacetArea( i ), ( elmt_ ).FacetAreaMapped( i ), fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearPrism_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 9; i++)
        j += ( elmt_ ).FacetArea( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetArea: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 9; i++)
        j += ( elmt_ ).FacetAreaMapped( i );
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetAreaMapped: "<< ticks << " " << j <<endl;

		file.close();
	}

	//test facet area
	//fT   FacetNormal( size_t surface, size_t ip, fT* nrml ) const;
	if ( m_bTestFacetNormals )
	{
		Point<3U> vecNormal;
		/*ignore return parameter, jacobian*/

    vecNormal = ( elmt_ ).FacetNormal(0U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], 0.8256797, fTolerance );
		_equal( vecNormal[1], 0.4321557, fTolerance );
		_equal( vecNormal[2], 0.3626217, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(1U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  -0.8381078, fTolerance );
		_equal( vecNormal[1],  0.4178808, fTolerance );
		_equal( vecNormal[2],  0.3506436, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(2U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.0363660, fTolerance );
		_equal( vecNormal[1], -0.7655377, fTolerance );
		_equal( vecNormal[2], -0.6423624, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(3U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4104289, fTolerance );
		_equal( vecNormal[1], -0.5861533, fTolerance );
		_equal( vecNormal[2],  0.6985503, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(4U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4104289, fTolerance );
		_equal( vecNormal[1], -0.5861533, fTolerance );
		_equal( vecNormal[2],  0.6985503, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(5U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.4104289, fTolerance );
		_equal( vecNormal[1], -0.5861533, fTolerance );
		_equal( vecNormal[2],  0.6985503, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(6U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], 0.8256797, fTolerance );
		_equal( vecNormal[1], 0.4321557, fTolerance );
		_equal( vecNormal[2], 0.3626217, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(7U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0], -0.8381078, fTolerance );
		_equal( vecNormal[1],  0.4178808, fTolerance );
		_equal( vecNormal[2],  0.3506436, fTolerance );

    vecNormal = ( elmt_ ).FacetNormal(8U);
		_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
		_equal( vecNormal[0],  0.0363660, fTolerance );
		_equal( vecNormal[1], -0.7655377, fTolerance );
		_equal( vecNormal[2], -0.6423624, fTolerance );

		//normal mapped
		Point<3U> vecNormalMapped;

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(0U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], 0.8256797, fTolerance );
		_equal( vecNormalMapped[1], 0.4321557, fTolerance );
		_equal( vecNormalMapped[2], 0.3626217, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(1U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  -0.8381078, fTolerance );
		_equal( vecNormalMapped[1],  0.4178808, fTolerance );
		_equal( vecNormalMapped[2],  0.3506436, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(2U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.0363660, fTolerance );
		_equal( vecNormalMapped[1], -0.7655377, fTolerance );
		_equal( vecNormalMapped[2], -0.6423624, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(3U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4104289, fTolerance );
		_equal( vecNormalMapped[1], -0.5861533, fTolerance );
		_equal( vecNormalMapped[2],  0.6985503, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(4U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4104289, fTolerance );
		_equal( vecNormalMapped[1], -0.5861533, fTolerance );
		_equal( vecNormalMapped[2],  0.6985503, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(5U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.4104289, fTolerance );
		_equal( vecNormalMapped[1], -0.5861533, fTolerance );
		_equal( vecNormalMapped[2],  0.6985503, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(6U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], 0.8256797, fTolerance );
		_equal( vecNormalMapped[1], 0.4321557, fTolerance );
		_equal( vecNormalMapped[2], 0.3626217, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(7U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0], -0.8381078, fTolerance );
		_equal( vecNormalMapped[1],  0.4178808, fTolerance );
		_equal( vecNormalMapped[2],  0.3506436, fTolerance );

    vecNormalMapped = ( elmt_ ).FacetNormalMapped(8U);
		_info("Mapped Normal is: [" << vecNormalMapped[0] << ", " << vecNormalMapped[1] << ", " << vecNormalMapped[2] << "]");
		_equal( vecNormalMapped[0],  0.0363660, fTolerance );
		_equal( vecNormalMapped[1], -0.7655377, fTolerance );
		_equal( vecNormalMapped[2], -0.6423624, fTolerance );

		//test mapped vs computed
		_info("Test mapped normal vs computed.");
	  for(size_t i = 0; i < 9; i++)
	    for(size_t j = 0; j < 3; j++)
        _equal( ( elmt_ ).FacetNormal( i )[j], ( elmt_ ).FacetNormalMapped( i )[j], fToleranceInternal );

		//compare time computed vs mapped
		ofstream file ("speed_compare.txt", ios::out|ios::app);
		size_t total_times(TIMES);
		file << "\nIsoparametricLinearPrism_Test<3>: ";

		clock_t ticks = clock();  double64 j(0);
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 9; i++)
        j += ( elmt_ ).FacetNormal( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormal: "<< ticks << " " << j <<endl;

		ticks = clock();
		for(size_t t = 0; t < total_times; t++)
			for(size_t i = 0; i < 9; i++)
        j += ( elmt_ ).FacetNormalMapped( i )[0];
		ticks = clock() - ticks;
		file <<"\n\tCPU clock ticks used for " << total_times << " FacetNormalMapped: "<< ticks << " " << j <<endl;

		file.close();
	}

	//test parametric facet area
	//testing: fT   ParametricFacetArea( size_t surface, fT surfacePhysicalDetJ ) const;
	if ( m_bTestParametricFacetArea )
	{
		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
      double64 fArea = ( elmt_ ).ParametricFacetArea( iFacet );
			_info("Area is: " << fArea);
			_equal( fArea,  elmt_.FV_Stencil()->FacetIntegrationWeight(iFacet,0U), fTolerance );
		}
	}

	//test parametric facet normals
	//fT   ParametricFacetNormal( size_t surface, size_t ip,double64* NRML ) const
	if ( m_bTestParametricFacetNormals )
	{
		Point<3U> vecNormal;

		for ( size_t iFacet = 0U; iFacet < elmt_.FV_Stencil()->Facets(); iFacet++ )
		{
			//ignore return parameter, jacobian
      vecNormal = ( elmt_ ).ParametricFacetNormal(iFacet);
			_info("Normal is: [" << vecNormal[0] << ", " << vecNormal[1] << ", " << vecNormal[2] << "]");
			_equal( vecNormal[0], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 0), fTolerance );
			_equal( vecNormal[1], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 1), fTolerance );
			_equal( vecNormal[2], elmt_.FV_Stencil()->UnitParametricNormalComponent(iFacet, 2), fTolerance );
		}
	}


	//test sector volumes
	//testing: fT   SectorVolume( size_t sector ) const;
	if ( m_bTestSectorVolumes )
	{
		double64 fVolSum(0.);
		double64 fSectorVolume(0.);
		for ( size_t iSector = 0U; iSector < elmt_.FV_Stencil()->Sectors(); iSector++ )
		{
      fSectorVolume = ( elmt_ ).SectorVolume( iSector );
			//_info("Volume is: " << fSectorVolume);
			_equal( fSectorVolume, 0.166666666666667, fTolerance );
			fVolSum+=fSectorVolume;
		}

		_info("Volume of the Element: " << elmt_.Volume() << " vs Volume Sum: " << fVolSum);
		_equal(fVolSum, elmt_.Volume(), fTolerance);
	}

	if(m_bProjectionOnFacetNormal)
	{
		//In previous versions the following part of the test used to be incorrect (or was not a real test!):
		//
		//which actually projects a vector variable (vVariable) native of physical space onto a
		//facet normal in parametric space.
		//Result from this had little sense (i.e. comparing the projections), unless the hexahedron in
		//physical space has a shape equal to the onein parametric space (can be scaled,but maintains the same normals)
		//
		//A valid test would be the integration of the vector field, using physical coordinates,
		//and comparing them to the integration in parametric space (they should yield the same value.
		//This could be done using a hexahedron without a Det(J)=1, like in this case.
		//Given that the element is linear, the Jacobian is practically constant across the whole element,
		//and hence this is simply an exchange of an area calculation for a jacobian determinant calculation.
		//in essence, both of these are the same, but there is not any immediate functionality that needs to be tested
		//here. (or compared).
		const size_t iNrOfFacets(elmt_.FV_Stencil()->Facets());
		for ( size_t iFacet = 0U; iFacet < iNrOfFacets; iFacet++ )
		{
			//const CSPINDEX& prop_key
			VectorVariable<3> vVariable;
			for(size_t iD= 0U; iD < 3; iD++) {
	     	vVariable.Flag(iD)=PLAIN;
	     	vVariable(iD)=3.;
			}

      double64 fProjectionVal = ( elmt_ ).ProjectionOnFacetNormal( iFacet, vVariable);
      double64 fProjectionValP = vVariable.DotProduct( ( elmt_ ).ParametricFacetNormal( iFacet ));

			_info("fProjectionVal:" << fProjectionVal << " vs. fProjectionValP:" << fProjectionValP);
		  _equal(fProjectionVal, fProjectionValP, fTolerance);
		}
	}

	delete feptr;
	delete fvptr;

} //end test


} // end csmp
