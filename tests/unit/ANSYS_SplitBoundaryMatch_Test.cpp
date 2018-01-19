#include "ANSYS_SplitBoundaryMatch_Test.h"

#include "ANSYS_Model3D.h"
#include "Boundary.h"
#include "PL_Utilities.h"
#include "VTK_Interface.h"
#include "variableOperations.h"

#include "Model.h"
#include "VTU_Interface.h"

//#include "catch.hpp"

using namespace std;

namespace csmp {

/**
    tests findSplitInterfaceElements() and the creation of a SplitBoundary from ANSYS
    split model.
 
    SKM 15/01/18
*/
void ANSYS_SplitBoundaryMatch_Test::run()
{
      const bool verbose(false);
      // BINARY IO
      // ---------
      if ( verbose ) {
           cout <<"\nStart simulation of - "<<this->getName()<<endl<<endl;
           cout <<"Building ModelOutput..."<<endl;
        }
      const string variablesFile("CSMP-variables.txt");
      //                          fileset       regions-file
//      ANSYS_Model3D model( "Model_Split_wall", "Model_Split_wall", variablesFile.c_str(), true );
      ANSYS_Model3D model( "NotSplitWall", "NotSplitWall", variablesFile.c_str(), true );

    // 0. visualising the model and its regions
    // ----------------------------------------
     double64 evar(1.);
     cerr <<"\nrun: assigning values to unique regions:\n";
     for ( auto it=model.UniqueRegionsBegin(); it!=model.UniqueRegionsEnd(); ++it ) {
          cerr <<"\t'"<< (*it).first <<"'";
          (*it).second.InputPropertyValue( "element variable", makeScalar(PLAIN,evar) );
          evar += 1.;
       }
     cerr << endl;
  
     // nodes at model perimeter
     const csmp::Index nvar_key(model.Database().StorageKey("nodal variable"));
     multiset<Point<3U> > split_nodes;
     Region<3U> model_domain(model.Region("Model"));
     model_domain.InputPropertyValue( "nodal variable", makeScalar(PLAIN,0.) );
     for ( auto nit=model_domain.PerimeterNodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
          if ( (*nit)->AtBoundary() != NOT and (*nit)->AtBoundary() != REGION_BOUNDARY )
            (*nit)->Store( nvar_key, makeScalar(PLAIN,2.) );
          else {
               (*nit)->Store( nvar_key, makeScalar(PLAIN,5.) );
               split_nodes.insert( (*nit)->Coordinate() );
            }
       }
     cerr <<"\nrun: points on ANSYS split boundary:\n";
     size_t counter(0);
     for (auto pt=split_nodes.begin(); pt!=split_nodes.end(); ++pt )
       cerr << fixed << setprecision(3.) <<"\n\t"<< counter++ <<": "<< (*pt)[0] <<" "<< (*pt)[1] <<" "<< (*pt)[2];
  
     VTK_Interface<3U>  vtk_output;
     vtk_output.OutputDataToVTK(model, "test-evar", "element variable", 0U );
     vtk_output.OutputDataToVTK(model, "test-nvar", "nodal variable", 0U );
  
  
    // 1. trying to manipulate coordinates so that they are identical within a given tolerance
    // ---------------------------------------------------------------------------------------


    // 1. create element correspondance set
    // ---------------------------------------------------------------------------------------
      Catch::Timer timer;
      timer.start();
  
      // reports local ids of shared face on either side of the split boundary
      // std::pair<std::pair<Element<3U>*,size_t>,std::pair<Element<3U>*,size_t> >
      set<OppositeElements>  interface_elmt_pairs;
      _test( findSplitInterfaceElements( model.Region("Model"), interface_elmt_pairs ) );
  
      const double correspondanceIdentificionTime( timer.getElapsedMilliseconds() );
      cerr <<"\n::run: time taken to identify boundary elements: "<< correspondanceIdentificionTime <<"\n";
  

   // 2. build CSMP SplitBoundary
   // ---------------------------------------------------------------------------------------

   // 3. change some property along split boundary to verify that assignments are made correctly
   // ---------------------------------------------------------------------------------------


   // 4. write altered element properties on either side to VTK
   // ---------------------------------------------------------------------------------------
  
} // end run



} // csmp
