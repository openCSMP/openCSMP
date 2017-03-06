#include "PropertyAtPointVisitor_Test.h"

#include "Model1D.h"
#include "ANSYS_Model2D.h"
#include "ANSYS_Model3D.h"

using namespace std;

namespace csmp{

PropertyAtPointVisitor_Test::PropertyAtPointVisitor_Test( bool verbose )
 : verbose_(verbose)
{
}

PropertyAtPointVisitor_Test::~PropertyAtPointVisitor_Test()
{
}


void PropertyAtPointVisitor_Test::run()
{
    if ( verbose_ ) {
        cout<<" #######################################################################################"<<endl;
        cout << "\nRun PropertyAtPointVisitor_Test..."<<endl;
        cout<<" #######################################################################################"<<endl;
      }
  
    run1D_MeshTests();

    run2D_MeshTests();

    run3D_MeshTests();
}



void PropertyAtPointVisitor_Test::run1D_MeshTests()
{
    IsoparametricLinear1DMesh_Test( "Auto Mesh", "PropertyAtPointVisitor_Test.txt",
                                    "PropertyAtPointVisitor_Test.txt",
                                    "IsoparametricLinearLineElement's",   false, 1.0e-8 );
}



void PropertyAtPointVisitor_Test::run2D_MeshTests()
{
    IsoparametricLinear2DMesh_Test( "square1x1_tri_struct", "square1x1_tri_struct",
                                    "PropertyAtPointVisitor_Test.txt","IsoparametricLinearTriangle's", false, 1.0e-8 );
  
    IsoparametricLinear2DMesh_Test( "square1x1_quad_struct","square1x1_quad_struct",
                                    "PropertyAtPointVisitor_Test.txt","IsoparametricLinearQuadrilateral's", false, 1.0e-8 );
}



void PropertyAtPointVisitor_Test::run3D_MeshTests()
{
    IsoparametricLinear3DMesh_Test( "box1x1x1_tetra_struct", "box1x1x1_tetra_struct",
                                    "PropertyAtPointVisitor_Test.txt","IsoparametricLinearTetrahedron's",   false, 1.0e-8 );
  
    IsoparametricLinear3DMesh_Test( "box1x1x1_hexa_struct",  "box1x1x1_hexa_struct",
                                    "PropertyAtPointVisitor_Test.txt","IsoparametricLinearHexahedron's",    false, 1.0e-8 );
  
    IsoparametricLinear3DMesh_Test( "mixed_mesh",            "mixed_mesh",
                                    "PropertyAtPointVisitor_Test.txt","Mixed Mesh",                         false, 1.0e-8 );
}




void PropertyAtPointVisitor_Test::IsoparametricLinear1DMesh_Test( const char* mesh_name,
                                                                  const char* regionfile_name,
                                                                  const char* varfile_name,
                                                                  const char* mesh_type,
                                                                  bool brute_force_search,
                                                                  double64 tolerance )
{
    ///--------------
    /// 1D Mesh Test
    ///--------------

    //------------------------------------
    // Geometry
    if ( verbose_ ) cout <<"Building Model..."<<endl;
    double64    dx_min(0.1);           // min space step
    double64    dx_max(0.5);           // max space step
    double64    length(1.0);           // length
    double64    width_of_transition_zone(10);

    ErfDensity  erfc_density;
    Model1D<1U> model1DRegion( mesh_name, varfile_name, length, dx_min, dx_max, width_of_transition_zone, &erfc_density );
    Element<1>* eptr;
    vector<Element<1>*>::iterator elements_begin    = model1DRegion.Region("Model").ElementsBegin();
    vector<Element<1>*>::iterator elements_end      = model1DRegion.Region("Model").ElementsEnd();
    vector<Node<1>*>::iterator nodes_begin          = model1DRegion.Region("Model").NodesBegin();
    vector<Node<1>*>::iterator nodes_end            = model1DRegion.Region("Model").NodesEnd();
    sort(elements_begin,elements_end);
    sort(nodes_begin,nodes_end);
    if ( verbose_ ) cout <<"Finished reading mesh..."<<endl;
    // end Geometry
    // -----------------------------------

    // -----------------------------------
    // Set points
    size_t index(0);
    std::map<size_t,std::vector<double64> > pXYZ;
    pXYZ.clear();
    // Nodes
    for(vector<Node<1>*>::iterator nit=nodes_begin;nit!=nodes_end;nit++)
    {
        Point<1> node=(*nit)->Coordinate();
        pXYZ[index].push_back(node[0]);
        index++;
    }
    // BaryCenter's of the Element's
    for(vector<Element<1>*>::iterator eit=elements_begin;eit!=elements_end;eit++)
    {
        Point<1> bc=(*eit)->BaryCenter();
        pXYZ[index].push_back(bc[0]);
        index++;
    }
    // -----------------------------------

    // -----------------------------------
    // Find points
    clock_t start = clock();
    PropertyAtPointVisitor<1> pAt1DRegion(model1DRegion,pXYZ,"node variable");
    if(!brute_force_search)
        pAt1DRegion.SetBruteForceOff();
    else
        pAt1DRegion.SetBruteForceOn();
    if ( verbose_ ) cout<<"Search for "<<index<<" points ( "<<mesh_type<<" ):"<<endl;
    model1DRegion.Accept(pAt1DRegion);
    pAt1DRegion.CheckResults();
    if ( verbose_ ) cout<<"End of search:"<<endl;
    clock_t end = clock();
    OutputElapsedTime(start,end);
    // -----------------------------------

    // -----------------------------------
    // Comparison
    if ( verbose_ ) cout << "\nPropertyAtPointVisitor_Test:: Testing search algorihm in 1D ( "<<mesh_type<<" )..."<<endl;
    cout.setf(ios_base::scientific);
    std::vector<double64> NI;
    for(size_t i=0;i<pXYZ.size();i++)
    {
        eptr = pAt1DRegion.ElementThatContains(i);
        eptr->CoordinateMatrix();
        eptr->FE()->N(NI,pXYZ[i]);

        vector<double64> xyz(1,0.0);
        for ( size_t k=0; k<eptr->Nodes(); k++ )
          for ( size_t j=0; j<1; j++ )
            xyz[j] +=  NI[k]* eptr->FE()->XYZ(k,j);

        if( std::abs(pXYZ[i][0]-xyz[0])>=tolerance )
        {
            if ( verbose_ ) {
                 cout<<" Real point:("<<pXYZ[i][0]<<"), Found point:("<<xyz[0]<<")"<<endl;
                 cout<<" Found in Element["<<eptr->Idx()<<"] of type: "<< parseFiniteElementType( eptr->FE_Type() ) << endl;
              }
        }
        _equal(pXYZ[i][0],xyz[0],tolerance);
    }
    // -----------------------------------

}




void PropertyAtPointVisitor_Test::IsoparametricLinear2DMesh_Test( const char* mesh_name, const char* regionfile_name,
                                                                  const char* varfile_name, const char* mesh_type,
                                                                  bool brute_force_search, double64 tolerance )
{

    /// ---------------
    /// 2D Mesh Test
    ///----------------

    //------------------------------------
    // Geometry
    if ( verbose_ ) cout <<"Building Model..."<<endl;
    ANSYS_Model2D model2DRegion( mesh_name,regionfile_name,varfile_name);
    Element<2>* eptr;
    vector<Element<2>*>::iterator elements_begin    = model2DRegion.Region("Model").ElementsBegin();
    vector<Element<2>*>::iterator elements_end      = model2DRegion.Region("Model").ElementsEnd();
    vector<Node<2>*>::iterator nodes_begin          = model2DRegion.Region("Model").NodesBegin();
    vector<Node<2>*>::iterator nodes_end            = model2DRegion.Region("Model").NodesEnd();
    sort(elements_begin,elements_end);
    sort(nodes_begin,nodes_end);
    if ( verbose_ ) cout <<"Finished reading mesh..."<<endl;
    //------------------------------------

    //------------------------------------
    // Set points
    size_t index(0);
    std::map<size_t,std::vector<double64> > pXYZ;
    pXYZ.clear();
    // Nodes
    for(vector<Node<2>*>::iterator nit=nodes_begin;nit!=nodes_end;nit++)
    {
        Point<2> node=(*nit)->Coordinate();
        pXYZ[index].push_back(node[0]);
        pXYZ[index].push_back(node[1]);
        index++;
    }
    // BaryCenter's of the Element's
    for(vector<Element<2>*>::iterator eit=elements_begin;eit!=elements_end;eit++)
    {
        Point<2> bc=(*eit)->BaryCenter();
        pXYZ[index].push_back(bc[0]);
        pXYZ[index].push_back(bc[1]);
        index++;
    }
    //------------------------------------

    //------------------------------------
    // Find Points
    clock_t start = clock();
    PropertyAtPointVisitor<2> pAt2DRegion(model2DRegion,pXYZ,"node variable");
    if( !brute_force_search )
        pAt2DRegion.SetBruteForceOff();
    else
        pAt2DRegion.SetBruteForceOn();
    if ( verbose_ ) cout<<"Search for "<<index<<" points ( "<<mesh_type<<" ):"<<endl;
    model2DRegion.Accept(pAt2DRegion);
    pAt2DRegion.CheckResults();
    if ( verbose_ ) cout<<"End of search:"<<endl;
    clock_t end = clock();
    OutputElapsedTime(start,end);
    //------------------------------------

    //------------------------------------
    // Comparison
    if ( verbose_ ) cout << "\nPropertyAtPointVisitor_Test:: Testing search algorihm in 2D ( "<<mesh_type<<" )..."<<endl;
    cout.setf(ios_base::scientific);
    std::vector<double64> NI;
    for(size_t i=0;i<pXYZ.size();i++)
    {
        eptr = pAt2DRegion.ElementThatContains(i);
        eptr->CoordinateMatrix();
        eptr->FE()->N(NI,pXYZ[i]);

        vector<double64> xyz(3,0.0);
        for ( size_t k=0; k<eptr->Nodes(); k++ )
          for ( size_t j=0; j<2; j++ )
            xyz[j] +=  NI[k]* eptr->FE()->XYZ(k,j);

        if( std::abs(pXYZ[i][0]-xyz[0])>=tolerance || std::abs(pXYZ[i][1]-xyz[1])>=tolerance )
        {
            if ( verbose_ ) {
                 cout<<" Real point:("<<pXYZ[i][0]<<","<<pXYZ[i][1]<<"), Found point:("<<xyz[0]<<","<<xyz[1]<<")"<<endl;
                 cout<<" Found in Element["<<eptr->Idx()<<"] of type: "<< parseFiniteElementType( eptr->FE_Type() ) << endl;
              }
        }
        _equal(pXYZ[i][0],xyz[0],tolerance);
        _equal(pXYZ[i][1],xyz[1],tolerance);
    }
    //------------------------------------

}




void PropertyAtPointVisitor_Test::IsoparametricLinear3DMesh_Test( const char* mesh_name, const char* regionfile_name,
                                                                  const char* varfile_name, const char* mesh_type,
                                                                  bool brute_force_search, double64 tolerance )
{

    ///----------------
    /// 3D Mesh Test
    ///----------------

    //------------------------------------
    // Geometry
    if ( verbose_ ) cout <<"Building Model..."<<endl;
    ANSYS_Model3D model3DRegion( mesh_name,
                                 regionfile_name,
                                 varfile_name,
                                 true,true,true,false);
    Element<3>* eptr;
    vector<Element<3>*>::iterator elements_begin    = model3DRegion.Region("Model").ElementsBegin();
    vector<Element<3>*>::iterator elements_end      = model3DRegion.Region("Model").ElementsEnd();
    vector<Node<3>*>::iterator nodes_begin          = model3DRegion.Region("Model").NodesBegin();
    vector<Node<3>*>::iterator nodes_end            = model3DRegion.Region("Model").NodesEnd();
    sort(elements_begin,elements_end);
    sort(nodes_begin,nodes_end);
    if ( verbose_ ) cout <<"Finished reading mesh..."<<endl;
    //------------------------------------

    //------------------------------------
    // Set points
    size_t index(0);
    std::map<size_t,std::vector<double64> > pXYZ;
    pXYZ.clear();
    // Nodes
    for(vector<Node<3>*>::iterator nit=nodes_begin;nit!=nodes_end;nit++)
    {
        Point<3> node=(*nit)->Coordinate();
        pXYZ[index].push_back(node[0]);
        pXYZ[index].push_back(node[1]);
        pXYZ[index].push_back(node[2]);
        index++;
    }
    // BaryCenter's of the Element's
    for(vector<Element<3>*>::iterator eit=elements_begin;eit!=elements_end;eit++)
    {
        Point<3> bc=(*eit)->BaryCenter();
        pXYZ[index].push_back(bc[0]);
        pXYZ[index].push_back(bc[1]);
        pXYZ[index].push_back(bc[2]);
        index++;
    }
    //------------------------------------

    //------------------------------------
    // Find Points
    clock_t start = clock();
    PropertyAtPointVisitor<3> pAt3DRegion(model3DRegion,pXYZ,"node variable");
    if( !brute_force_search )
        pAt3DRegion.SetBruteForceOff();
    else
        pAt3DRegion.SetBruteForceOn();
    if ( verbose_ ) cout<<"Search for "<<index<<" points ("<<mesh_type<<"):"<<endl;
    model3DRegion.Accept(pAt3DRegion);
    pAt3DRegion.CheckResults();
    if ( verbose_ ) cout<<"End of search:"<<endl;
    clock_t end = clock();
    OutputElapsedTime(start,end);
    //------------------------------------

    //------------------------------------
    // Comparison
    if ( verbose_ ) cout << "\nPropertyAtPointVisitor_Test:: Testing search algorihm in 3D at Reigon ("<<mesh_type<<") ..."<<endl;
    cout.setf(ios_base::scientific);
    std::vector<double64> NI;
    for( size_t i=0;i<pXYZ.size();i++ )
    {
        eptr = pAt3DRegion.ElementThatContains(i);
        eptr->CoordinateMatrix();
        eptr->FE()->N(NI,pXYZ[i]);
        vector<double64> xyz(3,0.0);
        for ( size_t k=0; k<eptr->Nodes(); k++ )
          for ( size_t j=0; j<3; j++ )
              xyz[j] +=  NI[k]*eptr->FE()->XYZ(k,j);

        if( std::abs(pXYZ[i][0]-xyz[0])>=tolerance || std::abs(pXYZ[i][1]-xyz[1])>=tolerance || std::abs(pXYZ[i][2]-xyz[2])>=tolerance)
        {
            if ( verbose_ ) {
                 cout<<" Real point:("<<pXYZ[i][0]<<","<<pXYZ[i][1]<<","<<pXYZ[i][2]<<"), Found point:("<<xyz[0]<<","<<xyz[1]<<","<<xyz[2]<<")"<<endl;
                 cout<<" Found in Element["<<eptr->Idx()<<"] of type: "<< parseFiniteElementType( eptr->FE_Type() ) << endl;
              }
        }
        _equal(pXYZ[i][0],xyz[0],tolerance);
        _equal(pXYZ[i][1],xyz[1],tolerance);
        _equal(pXYZ[i][2],xyz[2],tolerance);
    }
    //------------------------------------

}




void PropertyAtPointVisitor_Test::OutputElapsedTime( clock_t start, clock_t end )
{

    unsigned long millisec ( (end - start) * 1000 / CLOCKS_PER_SEC);

    if ( verbose_ ) cout<<"\nElapsed Time = "<<millisec<<" ms ("<<(double64)(millisec)/1000.<<" sec; "
                        <<(double64)(millisec)/60000.<<" min; "<<(double64)(millisec)/3600000.<<" hours)"<<endl;
    cout.flush();

}

}// end namespace csmp
