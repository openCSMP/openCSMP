#include "vset_makers.h"
#include "CSMP_definitions.h"
#include "VSet.h"
#include "Model.h"

#include "IsoparametricLinearHexahedron.h"
#include "IsoparametricLinearPyramid.h"
#include "IsoparametricLinearPrism.h"
#include "IsoparametricLinearTriangle.h"
#include "IsoparametricLinearQuadrilateral.h"




///this variable controls the perturbation of the node when it is skewed, the smaller it is, the smaller the perturbation will be
#define PERTURBATION 0.0001

using namespace std;

namespace csmp {

/**
    This function tests the creation of a VSet (not related to the EFVT class).
*/
VSet<2U> test_CreateVSet(std::ostream& os)
{
    const size_t iNodes(9);
    const size_t iNrOfElements(4);
    
  	IsoparametricLinearQuadrilateral iso_quadrilateral;
  	VSet< 2>  vset( iso_quadrilateral.Nodes(),
                    iso_quadrilateral.Neighbors(),
                    iso_quadrilateral.ElementType(), iNodes, iNrOfElements );
  	
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
    //element 0
    vecNodes[0]=0;
    vecNodes[1]=1;
    vecNodes[2]=4;
    vecNodes[3]=3;
    deqElements[0]=vecNodes;
    //element 1
    vecNodes[0]=1;
    vecNodes[1]=2;
    vecNodes[2]=5;
    vecNodes[3]=4;
    deqElements[1]=vecNodes;
    //element 2
    vecNodes[0]=4;
    vecNodes[1]=5;
    vecNodes[2]=8;
    vecNodes[3]=7;
    deqElements[2]=vecNodes;
    //element 3
    vecNodes[0]=3;
    vecNodes[1]=4;
    vecNodes[2]=7;
    vecNodes[3]=6;
    deqElements[3]=vecNodes;
    
    vset.AddPlist( deqElements.begin(),deqElements.end());

    //define neighbors
    std::deque<std::vector<long64> > deqElementNeighbors(iNrOfElements);
    std::vector<long64> vecNeighbors(4);
    //element 0
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=2;
    vecNeighbors[2]=4;
    vecNeighbors[3]=LEFT_OUTSIDE;
    deqElementNeighbors[0]=vecNeighbors;
    //element 1
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=3;
    vecNeighbors[3]=1;
    deqElementNeighbors[1]=vecNeighbors;
    //element 2
    vecNeighbors[0]=2;
    vecNeighbors[1]=RIGHT_OUTSIDE;
    vecNeighbors[2]=TOP_OUTSIDE;
    vecNeighbors[3]=4;
    deqElementNeighbors[2]=vecNeighbors;
    //element 3
    vecNeighbors[0]=1;
    vecNeighbors[1]=3;
    vecNeighbors[2]=TOP_OUTSIDE;
    vecNeighbors[3]=LEFT_OUTSIDE;
    deqElementNeighbors[3]=vecNeighbors;
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//define boundaries
  	vset.AddBFlag( 0, CNR1 );
    vset.AddBFlag( 1, BOTTOM_OUTSIDE );
    vset.AddBFlag( 2, CNR2 );
    vset.AddBFlag( 3, LEFT_OUTSIDE );
    vset.AddBFlag( 5, RIGHT_OUTSIDE );
    vset.AddBFlag( 6, CNR4 );
	  vset.AddBFlag( 7, TOP_OUTSIDE );
    vset.AddBFlag( 8, CNR3 );
    
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
    
    return vset;
}








void test_Create_One_Square_VSet(std::ostream& os, VSet<2U> & vset, double64 length_of_sides, bool bSkewed )
{    
  	IsoparametricLinearQuadrilateral iso_quad;
  	
  	//elements hexahedrons
    size_t         nodes(4); // number of nodes
  	deque<size_t>  npes(1);  // number of nodes per element
    deque<size_t>  epes(1);  // elements per element
    deque<int32>   etypes(1,ISOPARAMETRIC_LINEAR_QUADRILATERAL);

    npes[0]=iso_quad.Nodes();
  	epes[0]=iso_quad.Neighbors();
  	
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int32> vecElementTypes(1);
    vecElementTypes[0]= ISOPARAMETRIC_LINEAR_QUADRILATERAL;
  	
  	vset.Resize( etypes, npes, epes, nodes, 0, 0 );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  	
  	px[0]=0;                py[0]=0;                  pz[0]=0;
  	px[1]=length_of_sides;  py[1]=0;                  pz[1]=0;
  	px[2]=length_of_sides;  py[2]=length_of_sides;    pz[2]=0;
  	px[3]=0;                py[3]=length_of_sides;    pz[3]=0;
  	
  	if( bSkewed )
  	 for (size_t i = 0; i < 8; i++)
  	  {
  	    px[i]+= (rand()%2000)*PERTURBATION;
  	    py[i]+= (rand()%2000)*PERTURBATION;
  	    pz[i]+= (rand()%2000)*PERTURBATION;
  	  }
  	  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------ELEMENTS
    //define hexahedron elements (elements 0->26), assign nodes per element
    std::deque< std::vector<size_t> > deqElements(1);
    deqElements[0].resize(nodes);
  	 
    deqElements[0][0]= 1;
    deqElements[0][1]= 2;
    deqElements[0][2]= 3;
    deqElements[0][3]= 4;
    vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    std::deque<std::vector<long64> > deqElementNeighbors(1);
    deqElementNeighbors[0].resize(4);
  	deqElementNeighbors[0][0]= BACK_OUTSIDE;
  	deqElementNeighbors[0][1]= BOTTOM_OUTSIDE;
  	deqElementNeighbors[0][2]= RIGHT_OUTSIDE;
  	deqElementNeighbors[0][3]= TOP_OUTSIDE;
  	
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.AddBFlag( 1, CNR1);
  	vset.AddBFlag( 2, CNR2);
  	vset.AddBFlag( 3, CNR3);
  	vset.AddBFlag( 4, CNR4);
  	
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
}






void test_Create_TrianglePatch_VSet( std::ostream& os, VSet<2U> & vset )
{    
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int32> vecElementTypes(1);
    vecElementTypes[0]= ISOPARAMETRIC_LINEAR_TRIANGLE;
  	
    const size_t   nodes(9); // number of nodes
  	deque<size_t>  npes(10);  // number of nodes per element
    deque<size_t>  epes(10);  // elements per element
  	IsoparametricLinearTriangle iso_tria;
    npes[0]=iso_tria.Nodes();
  	epes[0]=iso_tria.Neighbors();
    deque<int32>   etypes(1,ISOPARAMETRIC_LINEAR_TRIANGLE);
    
  	vset.Resize( etypes, npes, epes, nodes, 0, 0 );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  	
  	px[0]=0.;                 py[0]=100.;               pz[0]=0.;
  	px[1]=52.73842909594504;  py[1]=100.;               pz[1]=0.;
  	px[2]=100.;               py[2]=100.;               pz[2]=0.;
  	px[3]=53.51413499621665;  py[3]=76.29281878271485;  pz[3]=0.;
  	px[4]=100.;               py[4]=48.94918579813998;  pz[4]=0.;
  	px[5]=26.17050201164178;  py[5]=23.35089108917626;  pz[5]=0.;
  	px[6]=79.5002826553162;   py[6]=28.97475886614556;  pz[6]=0.;
  	px[7]=0.;                 py[7]=0.;                 pz[7]=0.;
  	px[8]=100.;               py[8]=0.;                 pz[8]=0.;
  	  	  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------ELEMENTS
  	// define nodes per element
  	vector<size_t> node_dummy(3);
    std::deque< std::vector<size_t> > deqElements(10,node_dummy);
    deqElements[0][0]= 0;
    deqElements[0][1]= 3;
    deqElements[0][2]= 1;
    deqElements[1][0]= 1;
    deqElements[1][1]= 3;
    deqElements[1][2]= 2;
    deqElements[2][0]= 0;
    deqElements[2][1]= 5;
    deqElements[2][2]= 3;
    deqElements[3][0]= 3;
    deqElements[3][1]= 5;
    deqElements[3][2]= 6;
    deqElements[4][0]= 3;
    deqElements[4][1]= 6;
    deqElements[4][2]= 4;
    deqElements[5][0]= 2;
    deqElements[5][1]= 3;
    deqElements[5][2]= 4;
    deqElements[6][0]= 0;
    deqElements[6][1]= 7;
    deqElements[6][2]= 5;
    deqElements[7][0]= 5;
    deqElements[7][1]= 7;
    deqElements[7][2]= 8;
    deqElements[8][0]= 5;
    deqElements[8][1]= 8;
    deqElements[8][2]= 6;
    deqElements[9][0]= 4;
    deqElements[9][1]= 6;
    deqElements[9][2]= 8;

    vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors per element
  	vector<long64> nbor_dummy(3);
    std::deque<std::vector<long64> > deqElementNeighbors(10,nbor_dummy);
    deqElementNeighbors[0][0]= 1;
    deqElementNeighbors[0][1]= TOP_OUTSIDE;
    deqElementNeighbors[0][2]= 2;
    deqElementNeighbors[1][0]= 5;
    deqElementNeighbors[1][1]= TOP_OUTSIDE;
    deqElementNeighbors[1][2]= 0;
    deqElementNeighbors[2][0]= 3;
    deqElementNeighbors[2][1]= 0;
    deqElementNeighbors[2][2]= 6;
    deqElementNeighbors[3][0]= 8;
    deqElementNeighbors[3][1]= 4;
    deqElementNeighbors[3][2]= 2;
    deqElementNeighbors[4][0]= 9;
    deqElementNeighbors[4][1]= 5;
    deqElementNeighbors[4][2]= 3;
    deqElementNeighbors[5][0]= 4;
    deqElementNeighbors[5][1]= RIGHT_OUTSIDE;
    deqElementNeighbors[5][2]= 1;
    deqElementNeighbors[6][0]= 7;
    deqElementNeighbors[6][1]= 2;
    deqElementNeighbors[6][2]= LEFT_OUTSIDE;
    deqElementNeighbors[7][0]= BOTTOM_OUTSIDE;
    deqElementNeighbors[7][1]= 8;
    deqElementNeighbors[7][2]= 6;
    deqElementNeighbors[8][0]= 9;
    deqElementNeighbors[8][1]= 3;
    deqElementNeighbors[8][2]= 7;
    deqElementNeighbors[9][0]= 8;
    deqElementNeighbors[9][1]= RIGHT_OUTSIDE;
    deqElementNeighbors[9][2]= 4;
  	
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.AddBFlag( 7, CNR1);
  	vset.AddBFlag( 8, CNR2);
  	vset.AddBFlag( 2, CNR3);
  	vset.AddBFlag( 0, CNR4);
  	vset.AddBFlag( 1, TOP_OUTSIDE );
  	vset.AddBFlag( 4, RIGHT_OUTSIDE);
    
    vset.Out(os);
    
} // end test_Create_TrianglePatch_VSet








void test_Create_One_Hexahedra_VSet(std::ostream& os, VSet<3U> & vset, bool bSkewed )
{    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	
  	//elements hexahedrons
    size_t nodes(8);         //number of nodes
  	deque<size_t>  npes(1);  //number of nodes per element
    deque<size_t>  epes(1);  //elements per element
    deque<int32>   etypes(1,ISOPARAMETRIC_LINEAR_HEXAHEDRON);
    
    npes[0]=iso_hexahedron.Nodes();
  	epes[0]=iso_hexahedron.Neighbors();
  	
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int32> vecElementTypes(1);
    vecElementTypes[0]= ISOPARAMETRIC_LINEAR_HEXAHEDRON;
  	
  	vset.Resize( etypes, npes, epes, nodes, 0, 0 );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  	
  	px[0]=0;py[0]=0;pz[0]=0;
  	px[1]=1;py[1]=0;pz[1]=0;
  	px[2]=1;py[2]=1;pz[2]=0;
  	px[3]=0;py[3]=1;pz[3]=0;
  	px[4]=0;py[4]=0;pz[4]=1;
  	px[5]=1;py[5]=0;pz[5]=1;
  	px[6]=1;py[6]=1;pz[6]=1;
  	px[7]=0;py[7]=1;pz[7]=1;
  	
  	if( bSkewed )
  	 for (size_t i = 0; i < 8; i++)
  	  {
  	    px[i]+= (rand()%2000)*PERTURBATION;
  	    py[i]+= (rand()%2000)*PERTURBATION;
  	    pz[i]+= (rand()%2000)*PERTURBATION;
  	  }
  	  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------ELEMENTS
    //define hexahedron elements (elements 0->26), assign nodes per element
    std::deque<std::vector<size_t> > deqElements(1);
    deqElements[0].resize(8);
  	 
    deqElements[0][0]= 1;
    deqElements[0][1]= 2;
    deqElements[0][2]= 3;
    deqElements[0][3]= 4;
    deqElements[0][4]= 5;
    deqElements[0][5]= 6;
    deqElements[0][6]= 7;
    deqElements[0][7]= 8;
  	vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    std::deque<std::vector<long64> > deqElementNeighbors(1);
    deqElementNeighbors[0].resize(6);
  	deqElementNeighbors[0][0]= BACK_OUTSIDE;
  	deqElementNeighbors[0][1]= BOTTOM_OUTSIDE;
  	deqElementNeighbors[0][2]= RIGHT_OUTSIDE;
  	deqElementNeighbors[0][3]= TOP_OUTSIDE;
  	deqElementNeighbors[0][4]= LEFT_OUTSIDE;
  	deqElementNeighbors[0][5]= FRONT_OUTSIDE;

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.AddBFlag( 1, CNR1);
  	vset.AddBFlag( 2, CNR2);
  	vset.AddBFlag( 3, CNR3);
  	vset.AddBFlag( 4, CNR4);
  	vset.AddBFlag( 5, CNR5);
  	vset.AddBFlag( 6, CNR6);
  	vset.AddBFlag( 7, CNR7);
  	vset.AddBFlag( 8, CNR8);
    
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
}






void test_Create_Hexahedra_VSet(std::ostream& os, VSet<3U> & vset, bool bSkewed )
{
    const size_t iNrOfElements(27);
    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
  	
  	//elements hexahedrons
    size_t nodes(iDim_i*iDim_j*iDim_k);  //number of nodes: 64 on a 4x4x4 grid
    
    //------------------------CREATE VSET
    //this is a 3D model, it is a cube of hexahedron with six pyramid elements in the middle 	
  	vset.Resize( iso_hexahedron.Nodes(),
                 iso_hexahedron.Neighbors(),
                 iso_hexahedron.ElementType(),
                 nodes, iNrOfElements );

  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
 // 	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
    
    //--------------------------ELEMENTS
    //define hexahedron elements (elements 0->26), assign nodes per element
    std::deque<std::vector<size_t> > deqElements(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 const size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
      
  	 deqElements[iElement].resize(8);
  	 
  	 deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][6]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][7]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
    }
  	
  	vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    //define neighbors
    const long64  iDim_i_(iDim_i), iDim_j_(iDim_j), iDim_k_(iDim_k), iDim_km1_2_(iDim_km1_2);
    std::deque<std::vector<long64> > deqElementNeighbors(iNrOfElements);
    for( long64 k = 0; k < iDim_k_-1; k++) //z
  	for( long64 j = 0; j < iDim_j_-1; j++) //y
  	for( long64 i = 0; i < iDim_i_-1; i++) //x
  	{
  	 const size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     deqElementNeighbors[iElement].resize(6);
  	 
     //face 0
     deqElementNeighbors[iElement][0]= k==0?BACK_OUTSIDE:(1+ iDim_km1_2_*(k-1)+(iDim_j_-1)*j+i);
     //face 1
     deqElementNeighbors[iElement][1]= j==0?BOTTOM_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*(j-1)+i);
     //face 2
     deqElementNeighbors[iElement][2]= i==iDim_i-2?RIGHT_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*j+i+1);
     //face iDim_k-1
     deqElementNeighbors[iElement][3]= j==iDim_j-2?TOP_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*(j+1)+i);
     //face 4
     deqElementNeighbors[iElement][4]= i==0?LEFT_OUTSIDE:(1+ iDim_km1_2_*k+(iDim_j_-1)*j+(i-1));
     //face 5
     deqElementNeighbors[iElement][5]= k==iDim_k-2?FRONT_OUTSIDE:(1+ iDim_km1_2_*(k+1)+(iDim_j_-1)*j+i);
      	 
    }
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  int32 bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(j==0) bBoundary=CNR4; //shouldn't this be i?
          else if(j==(iDim_j-1)) bBoundary=CNR3; //shouldn't this be i?
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i)+1);
        vset.AddBFlag( iNode, bBoundary);
      }
  	}
  	
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
}











void test_Create_Square_VSet( std::ostream& os, VSet<2U>& vset, size_t size_sides, double64 dimension, bool bSkewed )
{
  if(size_sides==1)
    test_Create_One_Square_VSet( os, vset, dimension, bSkewed );
  else
    test_Create_SlitRectangle_VSet( os, vset, size_sides, size_sides,dimension,dimension, 0, bSkewed );
}








void test_Create_SlitRectangle_VSet( std::ostream& os, VSet<2U> & vset, size_t x_dimension, size_t y_dimension, 
                                     double64 x_length, double64 y_length, size_t depth_of_slit, bool bSkewed )
{
  assert(x_dimension>0);
  assert(y_dimension>0);
  assert(depth_of_slit<x_dimension);
  assert(x_length>1.e-7);
  assert(y_length>1.e-7);
  
  IsoparametricLinearQuadrilateral iso_quadrilateral;
  	
  //node dimensions of box
  size_t iDim_i(x_dimension+1U);//j - width
  size_t iDim_j(y_dimension+1U);//i - length
  	
  size_t iNrOfElements( (iDim_i-1)*(iDim_j-1) );
    
  //elements hexahedrons
  size_t nodes(iDim_i*iDim_j);  //number of nodes: 64 on a 4x4x4 grid
  
  os <<"\ntest_Create_SlitRectangle_VSet:\n";
  os << "\n\tDim i: " << iDim_i << " Dim j: " << iDim_j << " nr of elements: " << iNrOfElements << " nodes: " << nodes;
     
  //------------------------CREATE VSET
  //this is a 3D model, it is a cube of hexahedron with six pyramid elements in the middle 	
  vset.Resize( iso_quadrilateral.Nodes(),
               iso_quadrilateral.Neighbors(),
               iso_quadrilateral.ElementType(), 
               nodes, iNrOfElements );

  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  
  	double64 delta_x = x_length/static_cast<double64>(x_dimension);
  	double64 delta_y = y_length/static_cast<double64>(y_dimension);
  	
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  if(bSkewed)
  	  {
  	    px[(iDim_i)*j+i]=i*delta_x + (rand()%2000)*PERTURBATION;
  	    py[(iDim_i)*j+i]=j*delta_y + (rand()%2000)*PERTURBATION;
  	    pz[(iDim_i)*j+i]=0.;
  	  }
  	  else
  	  {
  	    px[(iDim_i)*j+i]=i*delta_x;
  	    py[(iDim_i)*j+i]=j*delta_y;
  	    pz[(iDim_i)*j+i]=0.;
  	  }
  	}
  	
    //--------------------------ELEMENTS
    //define quad elements (elements 0->26), assign nodes per element
    std::deque<std::vector<size_t> > deqElements(iNrOfElements);
    for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 const size_t iElement((iDim_i-1)*j+i);
      
  	 deqElements[iElement].resize(4);
  	 
  	 deqElements[iElement][0]= 1+ (iDim_i)*j+i;
  	 deqElements[iElement][1]= 1+ (iDim_i)*j+i+1;
  	 deqElements[iElement][2]= 1+ (iDim_i)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ (iDim_i)*(j+1)+i;
  	}
  	
    //---------------------------------NEIGHBORS
    //define neighbors
    const long64  iDim_i_(iDim_i), iDim_j_(iDim_j);
    std::deque<std::vector<long64> > deqElementNeighbors(iNrOfElements);
    for( long j = 0; j < iDim_j_-1; j++ ) //y
  	for( long i = 0; i < iDim_i_-1; i++ ) //x
  	{
  	 const size_t iElement((iDim_i-1)*j+i);
     
     deqElementNeighbors[iElement].resize(4);
  	 
	 	 //face 1
     deqElementNeighbors[iElement][1]= j==0?BOTTOM_OUTSIDE:(1+ (iDim_i_-1)*(j-1)+i);
     //face 2
     deqElementNeighbors[iElement][2]= i==iDim_i_-2?RIGHT_OUTSIDE:(1+ (iDim_i_-1)*j+i+1);
     //face iDim_k-1
     deqElementNeighbors[iElement][3]= j==iDim_j_-2?TOP_OUTSIDE:(1+ (iDim_i_-1)*(j+1)+i);
     //face 4
     deqElementNeighbors[iElement][4]= i==0?LEFT_OUTSIDE:(1+ (iDim_i_-1)*j+(i-1));
     
     const bool over_slit  = (j == y_dimension/2       && i+1 >= x_dimension-depth_of_slit);
     const bool under_slit = (j == (y_dimension/2 - 1) && i+1 >= x_dimension-depth_of_slit);
     
     //os << "\n\nTo be over the slit: j("<<j<<") == " <<  y_dimension/2 << " and i("<<i<<") > " << x_dimension-depth_of_slit;
     //os << "\n\nTo be under the slit: j("<<j<<") == " <<  y_dimension/2-1 << " and i("<<i<<") > " << x_dimension-depth_of_slit;
     //os << "\nElement " << iElement << "-> i,j:" << i << "," << j << " under slit? " << (under_slit?"yes":"no") << " over slit? " << (over_slit?"yes":"no");
     if (over_slit)
  	  deqElementNeighbors[iElement][1]= static_cast<int32>(IRREGULAR_OUTSIDE);
     else if (under_slit)
  	  deqElementNeighbors[iElement][3]= static_cast<int32>(IRREGULAR_OUTSIDE);
  	   
    }
    
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  int32 bBoundary = NOT;
  	  
      if(j==0)
      {
        if(i==0) bBoundary=LEFT_OUTSIDE;
        else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
        else bBoundary=BOTTOM_OUTSIDE;
      }
      else if(j==(iDim_j-1))
      {
        if(i==0) bBoundary=LEFT_OUTSIDE;
        else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
        else bBoundary=TOP_OUTSIDE;
      }
      else //j in the middle
      {
        if(i==0) bBoundary=LEFT_OUTSIDE;
        else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
        else ;//do nothing: no boundary
      }
  	 
  	 if(bBoundary!=NOT)
     {
        const size_t iNode(((iDim_i)*j+i)+1);
        vset.AddBFlag( iNode, bBoundary);
     }
     }
    
    os << "\n\tIntroduce Slit Nodes...";
    
    //go over elements, introduce slit
    for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{ 
     const bool over_slit  = (j == y_dimension/2       && i >= x_dimension-depth_of_slit);
     const bool under_slit = (j == (y_dimension/2 - 1) && i >= x_dimension-depth_of_slit);
     
     const size_t iElement((iDim_i-1)*j+i);
    
     //os << "\nElement " << iElement << "-> i,j:" << i << "," << j << " under slit? " << (under_slit?"yes":"no") << " over slit? " << (over_slit?"yes":"no");
     
     if (over_slit)
  	 {
  	  if(i > x_dimension-depth_of_slit) // node is NOT at the end of the slit
  	  {
    	  //create new duplicate node (0)
    	  double64 new_px(px[(iDim_i)*j+i]), 
    	            new_py(py[(iDim_i)*j+i]), 
    	            new_pz(pz[(iDim_i)*j+i]);
    	  size_t node_number = px.size();
    	  px.push_back(new_px);
    	  py.push_back(new_py);
    	  pz.push_back(new_pz);
    	  //set new duplicate node
        deqElements[iElement][0]= 1+ node_number;
        vset.AddBFlag( 1+ node_number, IRREGULAR_OUTSIDE);
      }
      
      //create new duplicate node (1)
  	  double64 new_px = px[(iDim_i)*j+i+1]; 
  	  double64 new_py = py[(iDim_i)*j+i+1]; 
  	  double64  new_pz = pz[(iDim_i)*j+i+1];
  	  size_t node_number = px.size();
  	  px.push_back(new_px);
  	  py.push_back(new_py);
  	  pz.push_back(new_pz);
  	  //set new duplicate node
      deqElements[iElement][1]= 1+ node_number;
      
      //set boundary
  	  vset.AddBFlag( 1+ node_number, IRREGULAR_OUTSIDE);
     }
  	 
  	 if (under_slit)
  	 {
  	  //create new duplicate node (2)
  	  double64 new_px(px[(iDim_i)*(j+1)+i+1]), 
  	            new_py(py[(iDim_i)*(j+1)+i+1]), 
  	            new_pz(pz[(iDim_i)*(j+1)+i+1]);
  	  size_t node_number = px.size();
  	  px.push_back(new_px);
  	  py.push_back(new_py);
  	  pz.push_back(new_pz);
  	  //set new duplicate node
      deqElements[iElement][2]= 1+ node_number;
      //set boundary
      vset.AddBFlag( 1+ node_number, IRREGULAR_OUTSIDE);

  	  if(i > x_dimension-depth_of_slit) // node is NOT at the end of the slit
  	  {
        //create new duplicate node (3)
    	  double64 new_px = px[(iDim_i)*(j+1)+i]; 
    	  double64 new_py = py[(iDim_i)*(j+1)+i]; 
    	  double64 new_pz = pz[(iDim_i)*(j+1)+i];
    	  node_number = px.size();
    	  px.push_back(new_px);
    	  py.push_back(new_py);
    	  pz.push_back(new_pz);
    	  //set new duplicate node
        deqElements[iElement][3]= 1+ node_number;
        
        //set boundary
        vset.AddBFlag( 1+ node_number, IRREGULAR_OUTSIDE);
      }
      
  	 }	 
  	}
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
    vset.AddPlist( deqElements.begin(),deqElements.end());
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	    
  vset.EstablishZeroBasedNumbering();
  vset.Out(os);
}








//adapted ??
void test_Create_Pyramid_Hexa_VSet(std::ostream& os, VSet<3U> & vset, bool bSkewed )
{
    const size_t iNrOfElements(32/*26 hexahedrons + 6 pyramids*/);
    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	IsoparametricLinearPyramid    iso_pyramid;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
  	const size_t iPyramidsPlacement(13);
  	
  	//this is a 3D model, it is a cube of hexahedron with six pyramid elements in the middle
  	
  	//elements 0->7 are hexahedrons
  	//elements 8->13 are pyramids
    size_t nodes((iDim_i*iDim_j*iDim_k)+1);  //number of nodes: 64 on a 4x4x4 grid + 1 barycenter
  	deque<size_t>  npes(iNrOfElements);  //number of nodes per element
    deque<size_t>  epes(iNrOfElements);  //element type per element
    deque<int32>   etypes(iNrOfElements,ISOPARAMETRIC_LINEAR_HEXAHEDRON);

// NB: the pyramid elements still need to be dealt with

    for(size_t iElement = 0; iElement < 26U; iElement++)
  	{
  	  npes[iElement]=iso_hexahedron.Nodes();
  	  epes[iElement]=iso_hexahedron.Neighbors();
  	}
  	for(size_t iElement = 26U; iElement < 32U; iElement++)
  	{
  	  npes[iElement]=iso_pyramid.Nodes();
  	  epes[iElement]=iso_pyramid.Neighbors();
  	}
    
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int32> vecElementTypes(iNrOfElements);
    for(size_t iElement = 0; iElement < 26U; iElement++)
  	{
  	  vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_HEXAHEDRON;
  	}
  	for(size_t iElement = 26U; iElement < 32U; iElement++)
  	{
      vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_PYRAMID;
      etypes[iElement] = ISOPARAMETRIC_LINEAR_PRISM;
  	}

    vset.Resize( etypes, npes, epes, nodes, 0, 0 );
    vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

    
  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	  
  	}
  	
  	//add barycenter 
  	px[64]=3./2. + (rand()%2000)*PERTURBATION;
  	py[64]=3./2. + (rand()%2000)*PERTURBATION;
  	pz[64]=3./2. + (rand()%2000)*PERTURBATION;
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
    
    //--------------------------ELEMENTS
    //define hexahedron elements (elements 0->31), assign nodes per element
    std::deque<std::vector<size_t> > deqElements(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     if(i==1 && j==1 && k==1) //its the center element (6 pyramids)
  	  continue;

     if(iElement>iPyramidsPlacement)
      iElement--;
     
  	 deqElements[iElement].resize(8);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][6]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][7]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
    }
  	
  	size_t node_center(64);
    //element 26, assign nodes per element
    std::vector<size_t> vecNodes(5);
    vecNodes[0]= 1+ 21;
    vecNodes[1]= 1+ 22;
    vecNodes[2]= 1+ 26;
    vecNodes[3]= 1+ 25;
    vecNodes[4]= 1+ node_center;
    deqElements[26]=vecNodes;

    //element 27
    vecNodes[0]= 1+ 37;
    vecNodes[1]= 1+ 38;
    vecNodes[2]= 1+ 22;
    vecNodes[3]= 1+ 21;
    vecNodes[4]= 1+ node_center;
    deqElements[27]=vecNodes;

    //element 28
    vecNodes[0]= 1+ 22;
    vecNodes[1]= 1+ 38;
    vecNodes[2]= 1+ 42;
    vecNodes[3]= 1+ 26;
    vecNodes[4]= 1+ node_center;
    deqElements[28]=vecNodes;

    //element 29
    vecNodes[0]= 1+ 25;
    vecNodes[1]= 1+ 26;
    vecNodes[2]= 1+ 42;
    vecNodes[3]= 1+ 41;
    vecNodes[4]= 1+ node_center;
    deqElements[29]=vecNodes;
    
    //element 30
    vecNodes[0]= 1+ 37;
    vecNodes[1]= 1+ 21;
    vecNodes[2]= 1+ 25;
    vecNodes[3]= 1+ 41;
    vecNodes[4]= 1+ node_center;
    deqElements[30]=vecNodes;

    //element 31
    vecNodes[0]= 1+ 38;
    vecNodes[1]= 1+ 37;
    vecNodes[2]= 1+ 41;
    vecNodes[3]= 1+ 42;
    vecNodes[4]= 1+ node_center;
    deqElements[31]=vecNodes;
    
    vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    const long64  iDim_i_(iDim_i), iDim_j_(iDim_j), iDim_k_(iDim_k);
    //define neighbors
    std::deque<std::vector<long64> >  deqElementNeighbors(iNrOfElements);
    for( long64 k = 0U; k < iDim_k_-1; k++ ) //z
  	for( long64 j = 0U; j < iDim_j_-1; j++ ) //y
  	for( long64 i = 0U; i < iDim_i_-1; i++ ) //x
  	{
  	 if(i==1 && j==1 && k==1) //its the center element (6 pyramids)
  	  continue;
  	 
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     if(iElement>iPyramidsPlacement)
      iElement--;
      
     deqElementNeighbors[iElement].resize(6);
  
     //face 0
     size_t iNeighbor(1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i);
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;
     deqElementNeighbors[iElement][0U]= (k==0) ? BACK_OUTSIDE : iNeighbor;
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j-1)+i;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;
     deqElementNeighbors[iElement][1]= (j==0) ?BOTTOM_OUTSIDE:iNeighbor;
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+1;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][2]= (i==iDim_i_-2) ?RIGHT_OUTSIDE:iNeighbor;
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j+1)+i;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][3]= (j==iDim_j_-2) ?TOP_OUTSIDE:iNeighbor;
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+(i-1);
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][4]= (i==0) ?LEFT_OUTSIDE:iNeighbor;
     //face 5
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i;
     if(iNeighbor>iPyramidsPlacement)
      iNeighbor--;     
     deqElementNeighbors[iElement][5]= (k==iDim_k_-2) ?FRONT_OUTSIDE:iNeighbor;
  
    }
    
    //hexa-elements with pyramid neighbors are (6): 
    
    //element 5  - face 5 -> neighbor: 26 (under)
    deqElementNeighbors[4][5]=1+26;
    
    //element 22 - face 0 -> neighbor: 31 (over)
    deqElementNeighbors[21][0]=1+31;
    
    //element 13 - face 2 -> neighbor: 30 (left)
    deqElementNeighbors[12][2]=1+30;

    //element 14 - face 4 -> neighbor: 28 (right)
    deqElementNeighbors[13][4]=1+28;
    
    //element 11 - face 3 -> neighbor: 27 (front)
    deqElementNeighbors[10][3]=1+27;
    
    //element 16 - face 1 -> neighbor: 29 (front)
    deqElementNeighbors[15][1]=1+29;
    
    //pyramid neighbors are:
    std::vector<long64> vecNeighbors(5);
    
    //element 26
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+30;
    vecNeighbors[4]=1+4;
    deqElementNeighbors[26]=vecNeighbors;

    //element 27
    vecNeighbors[0]=1+31;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+26;
    vecNeighbors[3]=1+30;
    vecNeighbors[4]=1+10;
    deqElementNeighbors[27]=vecNeighbors;

    //element 28
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+31;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+26;
    vecNeighbors[4]=1+13;
    deqElementNeighbors[28]=vecNeighbors;

    ///element 29
    vecNeighbors[0]=1+26;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+31;
    vecNeighbors[3]=1+30;
    vecNeighbors[4]=1+15;
    deqElementNeighbors[29]=vecNeighbors;

    //element 30
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+26;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+31;
    vecNeighbors[4]=1+12;
    deqElementNeighbors[30]=vecNeighbors;

    //element 31
    vecNeighbors[0]=1+27;
    vecNeighbors[1]=1+30;
    vecNeighbors[2]=1+29;
    vecNeighbors[3]=1+28;
    vecNeighbors[4]=1+21;
    deqElementNeighbors[31]=vecNeighbors;

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  int32 bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(j==0) bBoundary=CNR4;
          else if(j==(iDim_j-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i)+1);
        vset.AddBFlag( iNode, bBoundary);
      }
  	}
  	
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
}














void test_Create_Prism_Hexa_VSet(std::ostream& os, VSet<3U> & vset, bool bSkewed )
{
    const size_t iNrOfElements(30/*24 hexahedrons + 6 prisms*/);
    
  	IsoparametricLinearHexahedron iso_hexahedron;
  	IsoparametricLinearPrism iso_prism;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	
  	//this is a 3D model, it is a cube of hexahedra with prisms in the middle
  	
  	//elements 0->23 are hexahedrons
  	//elements 8->29 are prisms
    const size_t   nodes(iDim_i*iDim_j*iDim_k);  //number of nodes: 64 on a 4x4x4 grid
  	deque<size_t>  npes(iNrOfElements);  //number of nodes per element
    deque<size_t>  epes(iNrOfElements);  //element type per element
    deque<int32>   etypes(iNrOfElements,ISOPARAMETRIC_LINEAR_HEXAHEDRON);
    
    size_t iElement = 0;
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	  if(i==1 && j==1) //its a prism
  	  { 
  	    npes[iElement]=iso_prism.Nodes();
  	    epes[iElement]=iso_prism.Neighbors();
  	    etypes[iElement] = ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    
  	    npes[iElement]=iso_prism.Nodes();
  	    epes[iElement]=iso_prism.Neighbors();
  	    etypes[iElement] = ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    continue; 
  	  }
  	
  	  npes[iElement]=iso_hexahedron.Nodes();
  	  epes[iElement]=iso_hexahedron.Neighbors();
  	  
  	  iElement++;
  	}
  	os << "e:" << iElement;
  	vset.Resize( etypes, npes, epes, nodes, 0, 0 );
    
    //--------------------------ELEMENT TYPES
  	//add element types
  	iElement = 0;
  
    vector<int32> vecElementTypes(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	  if(i==1 && j==1) //its a prism
  	  {
  	    vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    
  	    vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_PRISM;
  	    iElement++; 
  	    continue; 
  	  }
  	  
  	  vecElementTypes[iElement]= ISOPARAMETRIC_LINEAR_HEXAHEDRON;
  	  iElement++;
  	}
    vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

    
  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
//  	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
    
    //--------------------------ELEMENTS
    //define hexahedron elements (elements 0->26), assign nodes per element
    std::deque<std::vector<size_t> > deqElements(iNrOfElements);
    
    iElement = 0;
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 
     //size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     if(i==1 && j==1) //its the center element (6 pyramids)
  	 {
  	  const size_t node1 = 1+ iDim_k2*k+(iDim_j)*j+i;
  	  const size_t node2 = 1+ iDim_k2*k+(iDim_j)*j+i+1;;
  	  const size_t node3 = 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	  const size_t node4 = 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	  const size_t node5 = 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	  const size_t node6 = 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	  const size_t node7 = 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	  const size_t node8 = 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	  
  	  deqElements[iElement].resize(6);
	    deqElements[iElement][0]= node1;
  	  deqElements[iElement][1]= node2;
  	  deqElements[iElement][2]= node4;
  	  deqElements[iElement][3]= node5;
  	  deqElements[iElement][4]= node6;
  	  deqElements[iElement][5]= node8;
  	  
  	  iElement++;
      
      deqElements[iElement].resize(6);
	    deqElements[iElement][0]= node3;
  	  deqElements[iElement][1]= node4;
  	  deqElements[iElement][2]= node2;
  	  deqElements[iElement][3]= node7;
  	  deqElements[iElement][4]= node8;
  	  deqElements[iElement][5]= node6;
  	  
  	  iElement++;
      
  	  continue;
  	 }
 
  	 deqElements[iElement].resize(8);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][6]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][7]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
    
     iElement++;
    
    }
    
    vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    //define neighbors
    iElement = 0;
    std::deque<std::vector<long64> > deqElementNeighbors(iNrOfElements);
    /*for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 
  	 if(i==1 && j==1) //its a prism
  	 { iElement++; iElement++;
      continue; }
  	 
  	 deqElementNeighbors[iElement].resize(6);
  
     //face 0
     size_t iNeighbor(1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i);
     deqElementNeighbors[iElement][0]= (k==0?BACK_OUTSIDE:iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j-1)+i;
     deqElementNeighbors[iElement][1]= (j==0?BOTTOM_OUTSIDE:iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+1;
     deqElementNeighbors[iElement][2]= (i==iDim_i-2?RIGHT_OUTSIDE:iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j+1)+i;
     deqElementNeighbors[iElement][3]= (j==iDim_j-2?TOP_OUTSIDE:iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+(i-1);
     deqElementNeighbors[iElement][4]= (i==0?LEFT_OUTSIDE:iNeighbor);
     //face 5
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i;
     deqElementNeighbors[iElement][5]= (k==iDim_k-2?FRONT_OUTSIDE:iNeighbor);
     
     iElement++;
  	 
    }*/
    deqElementNeighbors[0].resize(6);
    deqElementNeighbors[0][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[0][1] = FRONT_OUTSIDE;
    deqElementNeighbors[0][2] = 1+1;
    deqElementNeighbors[0][3] = 1+3;
    deqElementNeighbors[0][4] = LEFT_OUTSIDE;
    deqElementNeighbors[0][5] = 1+10;
    
    deqElementNeighbors[1].resize(6);
    deqElementNeighbors[1][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[1][1] = FRONT_OUTSIDE;
    deqElementNeighbors[1][2] = 1+2;
    deqElementNeighbors[1][3] = 1+4;
    deqElementNeighbors[1][4] = 1+0;
    deqElementNeighbors[1][5] = 1+11;
    
    deqElementNeighbors[2].resize(6);
    deqElementNeighbors[2][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[2][1] = FRONT_OUTSIDE;
    deqElementNeighbors[2][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[2][3] = 1+6;
    deqElementNeighbors[2][4] = 1+1;
    deqElementNeighbors[2][5] = 1+12;
    
    deqElementNeighbors[3].resize(6);
    deqElementNeighbors[3][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[3][1] = 1+0;
    deqElementNeighbors[3][2] = 1+4;
    deqElementNeighbors[3][3] = 1+7;
    deqElementNeighbors[3][4] = LEFT_OUTSIDE;
    deqElementNeighbors[3][5] = 1+13;
    
    deqElementNeighbors[6].resize(6);
    deqElementNeighbors[6][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[6][1] = 1+2;
    deqElementNeighbors[6][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[6][3] = 1+9;
    deqElementNeighbors[6][4] = 1+5;
    deqElementNeighbors[6][5] = 1+16;
    
    deqElementNeighbors[7].resize(6);
    deqElementNeighbors[7][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[7][1] = 1+3;
    deqElementNeighbors[7][2] = 1+8;
    deqElementNeighbors[7][3] = BACK_OUTSIDE;
    deqElementNeighbors[7][4] = LEFT_OUTSIDE;
    deqElementNeighbors[7][5] = 1+17;
    
    deqElementNeighbors[8].resize(6);
    deqElementNeighbors[8][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[8][1] = 1+5;
    deqElementNeighbors[8][2] = 1+9;
    deqElementNeighbors[8][3] = BACK_OUTSIDE;
    deqElementNeighbors[8][4] = 1+7;
    deqElementNeighbors[8][5] = 1+18;
    
    deqElementNeighbors[9].resize(6);
    deqElementNeighbors[9][0] = BOTTOM_OUTSIDE;
    deqElementNeighbors[9][1] = 1+6;
    deqElementNeighbors[9][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[9][3] = BACK_OUTSIDE;
    deqElementNeighbors[9][4] = 1+8;
    deqElementNeighbors[9][5] = 1+19;
    
    deqElementNeighbors[10].resize(6);
    deqElementNeighbors[10][0] = 1+0;
    deqElementNeighbors[10][1] = FRONT_OUTSIDE;
    deqElementNeighbors[10][2] = 1+11;
    deqElementNeighbors[10][3] = 1+13;
    deqElementNeighbors[10][4] = LEFT_OUTSIDE;
    deqElementNeighbors[10][5] = 1+20;
    
    deqElementNeighbors[11].resize(6);
    deqElementNeighbors[11][0] = 1+1;
    deqElementNeighbors[11][1] = FRONT_OUTSIDE;
    deqElementNeighbors[11][2] = 1+12;
    deqElementNeighbors[11][3] = 1+14;
    deqElementNeighbors[11][4] = 1+10;
    deqElementNeighbors[11][5] = 1+21;
    
    deqElementNeighbors[12].resize(6);
    deqElementNeighbors[12][0] = 1+2;
    deqElementNeighbors[12][1] = FRONT_OUTSIDE;
    deqElementNeighbors[12][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[12][3] = 1+16;
    deqElementNeighbors[12][4] = 1+11;
    deqElementNeighbors[12][5] = 1+22;
    
    deqElementNeighbors[13].resize(6);
    deqElementNeighbors[13][0] = 1+3;
    deqElementNeighbors[13][1] = 1+10;
    deqElementNeighbors[13][2] = 1+14;
    deqElementNeighbors[13][3] = 1+17;
    deqElementNeighbors[13][4] = LEFT_OUTSIDE;
    deqElementNeighbors[13][5] = 1+23;
    
    deqElementNeighbors[16].resize(6);
    deqElementNeighbors[16][0] = 1+6;
    deqElementNeighbors[16][1] = 1+12;
    deqElementNeighbors[16][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[16][3] = 1+19;
    deqElementNeighbors[16][4] = 1+15;
    deqElementNeighbors[16][5] = 1+26;
    
    deqElementNeighbors[17].resize(6);
    deqElementNeighbors[17][0] = 1+7;
    deqElementNeighbors[17][1] = 1+13;
    deqElementNeighbors[17][2] = 1+18;
    deqElementNeighbors[17][3] = BACK_OUTSIDE;
    deqElementNeighbors[17][4] = LEFT_OUTSIDE;
    deqElementNeighbors[17][5] = 1+27;
    
    deqElementNeighbors[18].resize(6);
    deqElementNeighbors[18][0] = 1+8;
    deqElementNeighbors[18][1] = 1+15;
    deqElementNeighbors[18][2] = 1+19;
    deqElementNeighbors[18][3] = BACK_OUTSIDE;
    deqElementNeighbors[18][4] = 1+17;
    deqElementNeighbors[18][5] = 1+28;
    
    deqElementNeighbors[19].resize(6);
    deqElementNeighbors[19][0] = 1+9;
    deqElementNeighbors[19][1] = 1+16;
    deqElementNeighbors[19][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[19][3] = BACK_OUTSIDE;
    deqElementNeighbors[19][4] = 1+18;
    deqElementNeighbors[19][5] = 1+29;
    
    deqElementNeighbors[20].resize(6);
    deqElementNeighbors[20][0] = 1+10;
    deqElementNeighbors[20][1] = FRONT_OUTSIDE;
    deqElementNeighbors[20][2] = 1+21;
    deqElementNeighbors[20][3] = 1+23;
    deqElementNeighbors[20][4] = LEFT_OUTSIDE;
    deqElementNeighbors[20][5] = TOP_OUTSIDE;
      
    deqElementNeighbors[21].resize(6);
    deqElementNeighbors[21][0] = 1+11;
    deqElementNeighbors[21][1] = FRONT_OUTSIDE;
    deqElementNeighbors[21][2] = 1+22;
    deqElementNeighbors[21][3] = 1+24;
    deqElementNeighbors[21][4] = 1+20;
    deqElementNeighbors[21][5] = TOP_OUTSIDE;
      
    deqElementNeighbors[22].resize(6);
    deqElementNeighbors[22][0] = 1+12;
    deqElementNeighbors[22][1] = FRONT_OUTSIDE;
    deqElementNeighbors[22][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[22][3] = 1+26;
    deqElementNeighbors[22][4] = 1+21;
    deqElementNeighbors[22][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[23].resize(6);
    deqElementNeighbors[23][0] = 1+13;
    deqElementNeighbors[23][1] = 1+20;
    deqElementNeighbors[23][2] = 1+24;
    deqElementNeighbors[23][3] = 1+27;
    deqElementNeighbors[23][4] = LEFT_OUTSIDE;
    deqElementNeighbors[23][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[26].resize(6);
    deqElementNeighbors[26][0] = 1+16;
    deqElementNeighbors[26][1] = 1+22;
    deqElementNeighbors[26][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[26][3] = 1+29;
    deqElementNeighbors[26][4] = 1+25;
    deqElementNeighbors[26][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[27].resize(6);
    deqElementNeighbors[27][0] = 1+17;
    deqElementNeighbors[27][1] = 1+23;
    deqElementNeighbors[27][2] = 1+28;
    deqElementNeighbors[27][3] = BACK_OUTSIDE;
    deqElementNeighbors[27][4] = LEFT_OUTSIDE;
    deqElementNeighbors[27][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[28].resize(6);
    deqElementNeighbors[28][0] = 1+18;
    deqElementNeighbors[28][1] = 1+25;
    deqElementNeighbors[28][2] = 1+29;
    deqElementNeighbors[28][3] = BACK_OUTSIDE;
    deqElementNeighbors[28][4] = 1+27;
    deqElementNeighbors[28][5] = TOP_OUTSIDE;
    
    deqElementNeighbors[29].resize(6);
    deqElementNeighbors[29][0] = 1+19;
    deqElementNeighbors[29][1] = 1+26;
    deqElementNeighbors[29][2] = RIGHT_OUTSIDE;
    deqElementNeighbors[29][3] = BACK_OUTSIDE;
    deqElementNeighbors[29][4] = 1+28;
    deqElementNeighbors[29][5] = TOP_OUTSIDE;
    
    //prism neighbors are:
    std::vector<long64> vecNeighbors(5);
    
    //element 4
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=1+1;
    vecNeighbors[2]=1+5;
    vecNeighbors[3]=1+3;
    vecNeighbors[4]=1+14;
    deqElementNeighbors[4]=vecNeighbors;
  
    //element 5
    vecNeighbors[0]=BOTTOM_OUTSIDE;
    vecNeighbors[1]=1+8;
    vecNeighbors[2]=1+4;
    vecNeighbors[3]=1+6;
    vecNeighbors[4]=1+15;
    deqElementNeighbors[5]=vecNeighbors;

    //element 14
    vecNeighbors[0]=1+4;
    vecNeighbors[1]=1+11;
    vecNeighbors[2]=1+15;
    vecNeighbors[3]=1+13;
    vecNeighbors[4]=1+24;
    deqElementNeighbors[14]=vecNeighbors;

    //element 15
    vecNeighbors[0]=1+5;
    vecNeighbors[1]=1+18;
    vecNeighbors[2]=1+14;
    vecNeighbors[3]=1+16;
    vecNeighbors[4]=1+25;
    deqElementNeighbors[15]=vecNeighbors;

    //element 24
    vecNeighbors[0]=1+14;
    vecNeighbors[1]=1+21;
    vecNeighbors[2]=1+25;
    vecNeighbors[3]=1+23;
    vecNeighbors[4]=TOP_OUTSIDE;
    deqElementNeighbors[24]=vecNeighbors;

    //element 25
    vecNeighbors[0]=1+15;
    vecNeighbors[1]=1+28;
    vecNeighbors[2]=1+24;
    vecNeighbors[3]=1+26;
    vecNeighbors[4]=TOP_OUTSIDE;
    deqElementNeighbors[25]=vecNeighbors;

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  int32 bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(j==0) bBoundary=CNR4;
          else if(j==(iDim_j-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i)+1);
        vset.AddBFlag( iNode, bBoundary);
      }
  	}
  	
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
}












void test_Create_One_Prism_VSet(std::ostream& os, VSet<3U> & vset, bool bSkewed )
{    
  	IsoparametricLinearPrism iso_prism;
  	
  	//elements hexahedrons
    size_t nodes(iso_prism.Nodes());         //number of nodes
  	deque<size_t>  npes(1);  //number of nodes per element
    deque<size_t>  epes(1);  //elements per element
    deque<int32>      etypes(1,ISOPARAMETRIC_LINEAR_PRISM);

    npes[0]=iso_prism.Nodes();
  	epes[0]=iso_prism.Neighbors();
  	
    //--------------------------ELEMENT TYPES
  	//add element types
    vector<int32> vecElementTypes(1);
    vecElementTypes[0]= ISOPARAMETRIC_LINEAR_PRISM;
  	
  	vset.Resize( etypes, npes, epes, nodes, 0, 0 );
  	vset.AddElementTypes( vecElementTypes.begin(), vecElementTypes.end() );

  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  	
  	px[0]=0;py[0]=0;pz[0]=0;
  	px[1]=1;py[1]=0;pz[1]=0;
  	px[2]=1;py[2]=1;pz[2]=0;
  	px[3]=0;py[3]=0;pz[3]=1;
  	px[4]=1;py[4]=0;pz[4]=1;
  	px[5]=1;py[5]=1;pz[5]=1;
  	
  	if( bSkewed )
  	 for (size_t i = 0; i < 6; i++)
  	  {
  	    px[i]+= (rand()%2000)*PERTURBATION;
  	    py[i]+= (rand()%2000)*PERTURBATION;
  	    pz[i]+= (rand()%2000)*PERTURBATION;
  	  }
  	  
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
  	//--------------------------ELEMENTS
    //define prism element, assign nodes per element
    std::deque<std::vector<size_t> > deqElements(1);
    deqElements[0].resize(6);
  	 
    deqElements[0][0]= 1;
    deqElements[0][1]= 2;
    deqElements[0][2]= 3;
    deqElements[0][3]= 4;
    deqElements[0][4]= 5;
    deqElements[0][5]= 6;
  	vset.AddPlist( deqElements.begin(),deqElements.end());

     //---------------------------------NEIGHBORS
    //define neighbors
    std::deque<std::vector<long64> > deqElementNeighbors(1);
    deqElementNeighbors[0].resize(6);
  	deqElementNeighbors[0][0]= BACK_OUTSIDE;
  	deqElementNeighbors[0][1]= BOTTOM_OUTSIDE;
  	deqElementNeighbors[0][2]= RIGHT_OUTSIDE;
  	deqElementNeighbors[0][3]= TOP_OUTSIDE/*or LEFT_OUTSIDE*/;
  	deqElementNeighbors[0][4]= FRONT_OUTSIDE;

    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//----------------------------NODE BOUNDARIES
  	vset.AddBFlag( 1, CNR1);
  	vset.AddBFlag( 2, CNR2);
  	vset.AddBFlag( 3, CNR3);
  	vset.AddBFlag( 4, CNR5);
  	vset.AddBFlag( 5, CNR6);
  	vset.AddBFlag( 6, CNR7);
    
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
}






void test_Create_Prism_VSet(std::ostream& os, VSet<3U> & vset, bool bSkewed )
{
    const size_t iNrOfElements(54/*prisms*/);
    
  	IsoparametricLinearPrism iso_prism;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
  	
  	//this is a 3D model, it is a cube of hexahedron with six pyramid elements in the middle
  	
  	//elements 0->53 are prisms
    size_t nodes(iDim_i*iDim_j*iDim_k);  //number of nodes: 64 on a 4x4x4 grid

    vset.Resize( iso_prism.Nodes(),
                 iso_prism.Neighbors(),
                 iso_prism.ElementType(), 
                 nodes, iNrOfElements );
    
  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
//  	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
  	  if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
    
    //--------------------------ELEMENTS
    //define prism elements (elements 0->54), assign nodes per element
    std::deque<std::vector<size_t> > deqElements(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 //lower element
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
      
  	 deqElements[iElement].resize(6);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;

     //upper element
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i+27;

  	 deqElements[iElement].resize(6);
  	 deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][4]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][5]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 
    }
  	
  	vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    //define neighbors
    std::deque<std::vector<long64> > deqElementNeighbors(iNrOfElements);
    for(long64 k = 0U; k < iDim_k-1; k++) //z
  	for(long64 j = 0U; j < iDim_j-1; j++) //y
  	for(long64 i = 0U; i < iDim_i-1; i++) //x
  	{
  	 //lower element
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     deqElementNeighbors[iElement].resize(5);
  
     //face 0
     size_t iNeighbor(1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i);
     deqElementNeighbors[iElement][0]= (k==0?BACK_OUTSIDE:iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j-1)+i +27;
     deqElementNeighbors[iElement][1]= (j==0?BOTTOM_OUTSIDE:iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+1 +27;
     deqElementNeighbors[iElement][2]= (i==iDim_i-2?RIGHT_OUTSIDE:iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+27;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i;
     deqElementNeighbors[iElement][4]= (k==iDim_k-2?FRONT_OUTSIDE:iNeighbor);
     
     //upper element
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i+27;
     
     deqElementNeighbors[iElement].resize(5);

     //face 0
     iNeighbor = 1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i +27;
     deqElementNeighbors[iElement][0]= (k==0?BACK_OUTSIDE:iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j+1)+i;
     deqElementNeighbors[iElement][1]= (j==iDim_j-2?TOP_OUTSIDE:iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+(i-1);
     deqElementNeighbors[iElement][2]= (i==0?LEFT_OUTSIDE:iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i +27;
     deqElementNeighbors[iElement][4]= (k==iDim_k-2?FRONT_OUTSIDE:iNeighbor);
     
     
    }
    
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  int32 bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(j==0) bBoundary=CNR4;
          else if(j==(iDim_j-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i)+1);
        vset.AddBFlag( iNode, bBoundary);
      }
  	}
  	
    vset.EstablishZeroBasedNumbering();
    vset.Out(os);
}










void test_Create_Pyramid_VSet(std::ostream& os, VSet<3U> & vset, bool bSkewed )
{
  	IsoparametricLinearPyramid iso_pyramid;
  	
  	//node dimensions of box
  	const size_t iDim_k(4);//k - height
  	const size_t iDim_j(4);//j - width
  	const size_t iDim_i(4);//i - length
  	const size_t iDim_k2(iDim_k*iDim_k);//k - height
  	const size_t iDim_km1_2((iDim_k-1)*(iDim_k-1));//i - length
    
    const size_t iNrOfCells((iDim_i-1)*(iDim_j-1)*(iDim_k-1));
    const size_t iNrOfElements(iNrOfCells*6U/*pyramids*/);
     	
  	//this is a 3D model, it is a cube of 27 sets of six pyramids glued by the apex
  	
  	//elements 0->161 are pyramids
    size_t nodes(iDim_i*iDim_j*iDim_k+iNrOfCells);  //number of nodes: 64 on a 4x4x4 grid + one barycenter for each cell

    vset.Resize( iso_pyramid.Nodes(),
                 iso_pyramid.Neighbors(),
                 iso_pyramid.ElementType(),
                 nodes, iNrOfElements );
    
  	//-----------------------NODES
  	//define nodes
  	std::deque<double64> px(nodes);
  	std::deque<double64> py(nodes);
  	std::deque<double64> pz(nodes);
  
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
//  	  const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
      
      if(bSkewed)
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k + (rand()%2000)*PERTURBATION;
  	  }
  	  else
  	  {
  	    px[k*iDim_k2+(iDim_j)*j+i]=i;
  	    py[k*iDim_k2+(iDim_j)*j+i]=j;
  	    pz[k*iDim_k2+(iDim_j)*j+i]=k;
  	  }
  	}
  	
  	//add barycenters, one barycenter per cell
  	for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	  //this unused inside node boolean stays, just in case, in the future, we only want to skew internal nodes
  	  //const bool inside_node (!(i == 0 || j == 0 || k == 0 || i == iDim_i-1 || j == iDim_j-1 || k == iDim_k-1));
  	  
  	  if(bSkewed)
  	  {
    	  px[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=i + 1./2. + (rand()%2000)*PERTURBATION;
  	    py[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=j + 1./2. + (rand()%2000)*PERTURBATION;
  	    pz[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=k + 1./2. + (rand()%2000)*PERTURBATION;
      }
      else
      {
        px[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=i + 1./2.;
  	    py[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=j + 1./2.;
  	    pz[k*iDim_km1_2+(iDim_j-1)*j+i +iDim_i*iDim_j*iDim_k]=k + 1./2.;
      }
    }  	
  	
  	//load nodes
  	vset.AddXYZ( px, py, pz );
  	
    //--------------------------ELEMENTS
    //define pyramid elements (elements 0->54), assign nodes per element
    std::deque<std::vector<size_t> > deqElements(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 //we have 6 pyramids per cell
  	 
  	 //pyramid 0
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

  	 //pyramid 1
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*1;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;
  	 
  	 //pyramid 2
  	 os << "Element:(i="<<i<<",j="<<j<<",k="<<k<<")" << iDim_km1_2*k+(iDim_j-1)*j+i + 27*2 << endl;
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*2;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*j+i+1;
	   os << "Coord [" << iElement << "][0]: " << deqElements[iElement][0] << endl;
  	 deqElements[iElement][1]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
	   os << "Coord [" << iElement << "][1]: " << deqElements[iElement][1] << endl;
  	 deqElements[iElement][2]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
	   os << "Coord [" << iElement << "][2]: " << deqElements[iElement][2] << endl;
	   deqElements[iElement][3]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
	   os << "Coord [" << iElement << "][3]: " << deqElements[iElement][3] << endl;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;
	   os << "Coord [" << iElement << "][4]: " << deqElements[iElement][4] << endl;
  	 
  	 //pyramid 3
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*3;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][2]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

  	 //pyramid 4
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*4;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][1]= 1+ iDim_k2*k+(iDim_j)*j+i;
  	 deqElements[iElement][2]= 1+ iDim_k2*k+(iDim_j)*(j+1)+i;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

  	 //pyramid 5
  	 iElement = iDim_km1_2*k+(iDim_j-1)*j+i + iNrOfCells*5;
     
  	 deqElements[iElement].resize(5);
	   deqElements[iElement][0]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i+1;
  	 deqElements[iElement][1]= 1+ iDim_k2*(k+1)+(iDim_j)*j+i;
  	 deqElements[iElement][2]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i;
	   deqElements[iElement][3]= 1+ iDim_k2*(k+1)+(iDim_j)*(j+1)+i+1;
  	 deqElements[iElement][4]= 1+ iDim_km1_2*k+(iDim_j-1)*j+i+ iDim_i*iDim_j*iDim_k;

    }
  	
  	vset.AddPlist( deqElements.begin(),deqElements.end());

    //---------------------------------NEIGHBORS
    //define neighbors
    std::deque<std::vector<long64> > deqElementNeighbors(iNrOfElements);
    for(size_t k = 0U; k < iDim_k-1; k++) //z
  	for(size_t j = 0U; j < iDim_j-1; j++) //y
  	for(size_t i = 0U; i < iDim_i-1; i++) //x
  	{
  	 //pyramid 0
  	 size_t iElement(iDim_km1_2*k+(iDim_j-1)*j+i);
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     size_t iNeighbor(1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1);
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k-1)+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][4]= (k==0?BACK_OUTSIDE:iNeighbor);
     
     //pyramid 1
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j-1)+i +iNrOfCells*3;
     deqElementNeighbors[iElement][4]= (j==0?BOTTOM_OUTSIDE:iNeighbor);
     
      //pyramid 2
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i+1 +iNrOfCells*4;
     deqElementNeighbors[iElement][4]= (i==iDim_i-2?RIGHT_OUTSIDE:iNeighbor);
      
      //pyramid 3
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*(j+1)+i +iNrOfCells*1;
     deqElementNeighbors[iElement][4]= (j==iDim_j-2?TOP_OUTSIDE:iNeighbor);
     
     //pyramid 4
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i-1 +iNrOfCells*2;
     deqElementNeighbors[iElement][4]= (i==0?LEFT_OUTSIDE:iNeighbor);
     
     //pyramid 5
     iElement = iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*5;
     
     deqElementNeighbors[iElement].resize(5);
     //face 0
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*1;
     deqElementNeighbors[iElement][0]= (iNeighbor);
     //face 1
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*2;
     deqElementNeighbors[iElement][1]= (iNeighbor);
     //face 2
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*3;
     deqElementNeighbors[iElement][2]= (iNeighbor);
     //face 3
     iNeighbor = 1+ iDim_km1_2*k+(iDim_j-1)*j+i +iNrOfCells*4;
     deqElementNeighbors[iElement][3]= (iNeighbor);
     //face 4
     iNeighbor = 1+ iDim_km1_2*(k+1)+(iDim_j-1)*j+i +iNrOfCells*0;
     deqElementNeighbors[iElement][4]= (k==iDim_k-2?FRONT_OUTSIDE:iNeighbor);
     
    }
    
    vset.AddPfverts( deqElementNeighbors.begin(), deqElementNeighbors.end());
  	
  	//-----------------------------------------------------NODE BOUNDARIES
  	//define node boundaries
  	//nodes at corners:
    //nodes at edges:
  	//nodes at faces:
  	for(size_t k = 0U; k < iDim_k; k++) //z
  	for(size_t j = 0U; j < iDim_j; j++) //y
  	for(size_t i = 0U; i < iDim_i; i++) //x
  	{
  	  int32 bBoundary = NOT;
  	  
  	  if(k==0)
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR1;
          else if(i==(iDim_i-1)) bBoundary=CNR2;
          else bBoundary=EDGE1;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(j==0) bBoundary=CNR4;
          else if(j==(iDim_j-1)) bBoundary=CNR3;
          else bBoundary=EDGE3;
  	    }
  	    else //j is in the middle
  	    {
  	      if(i==0) bBoundary=EDGE4;
          else if(i==(iDim_i-1)) bBoundary=EDGE2;
          else bBoundary=BACK_OUTSIDE;
  	    }
  	  }
      else if(k==(iDim_k-1))
  	  {
  	    if(j==0)
  	    {
          if(i==0) bBoundary=CNR5;
          else if(i==(iDim_i-1)) bBoundary=CNR6;
          else bBoundary=EDGE9;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=CNR8;
          else if(i==(iDim_i-1)) bBoundary=CNR7;
          else bBoundary=EDGE11;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=EDGE12;
          else if(i==(iDim_i-1)) bBoundary=EDGE10;
          else bBoundary=FRONT_OUTSIDE;
  	    }
  	  }
  	  else //k is in the middle
  	  {
  	   if(j==0)
  	    {
          if(i==0) bBoundary=EDGE5;
          else if(i==(iDim_i-1)) bBoundary=EDGE6;
          else bBoundary=BOTTOM_OUTSIDE;
  	    }
  	    else if(j==(iDim_j-1))
  	    {
          if(i==0) bBoundary=EDGE8;
          else if(i==(iDim_i-1)) bBoundary=EDGE7;
          else bBoundary=TOP_OUTSIDE;
  	    }
  	    else //j in the middle
  	    {
  	      if(i==0) bBoundary=LEFT_OUTSIDE;
          else if(i==(iDim_i-1)) bBoundary=RIGHT_OUTSIDE;
          else ;//do nothing: no boundary
  	    }
  	  }
  	
      if(bBoundary!=NOT)
      {
        const size_t iNode((iDim_k2*k+(iDim_j)*j+i)+1);
        vset.AddBFlag( iNode, bBoundary);
      }
  	}
 
    vset.EstablishZeroBasedNumbering();
}





} // end csmp

/**
    USAGE EXAMPLE:

//Initialize Model:

VSet<3U>    mesh_container;
test_Create_Pyramid_VSet(mesh_container, true);
Model<3U>  model3D( mesh_container, true ); 
model3D.AddFaceTopology();
mesh_container.Erase();

//Initialize properties (depends on what you need...):
//model configuration
model3D.InputPropertyValue ( "permeability", 1.0e-12 );
model3D.InputPropertyValue ( "porosity", 0.25 );
model3D.InputPropertyValue ( "storativity", 1.0e-9 );
model3D.InputPropertyValue ( "fluid volume source", 0. );
model3D.InputPropertyValue ( "nodal fluid volume source", 0. );
model3D.InputPropertyValue ( "concentration", 0. );
model3D.InputPropertyValue ( "diffusivity", 1.0e-50 );

//remember to input boundary conditions as required

*/
