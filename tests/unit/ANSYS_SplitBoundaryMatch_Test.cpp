#include "ANSYS_SplitBoundaryMatch_Test.h"

#include "ANSYS_Model3D.h"
#include "Boundary.h"
#include "PL_Utilities.h"
#include "VTK_Interface.h"
#include "variableOperations.h"
#include "CSMP_mathUtilities.h"

#include "Model.h"
#include "VTU_Interface.h"

#include "catch.hpp"

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
      const string variablesFile("ANSYS_SplitBoundaryMatch_Test-variables.txt");
      //                          fileset       regions-file
//      ANSYS_Model3D model( "Model_Split_wall", "Model_Split_wall", variablesFile.c_str(), true );
//      ANSYS_Model3D model( "NotSplitWall", "NotSplitWall", variablesFile.c_str(), true );
      ANSYS_Model3D model( "Split_Edges", "Split_Edges", variablesFile.c_str(), true );
      Region<3U>    model_domain(model.Region("Model"));

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
     // painting internal boundary elements (they are flagged neither REGION_BOUNDARY nor IRREGULAR
     const double64 color(0.);
     const csmp::Index evar_key(model.Database().StorageKey("element variable"));
     for ( auto eit=model_domain.PerimeterElementsBegin(); eit!=model_domain.ElementsEnd(); ++eit )
       (*eit)->Store( evar_key, makeScalar(PLAIN,color) );

     // regularising node positions
    cerr <<"\nrun: focusing on inner box region that has "<< model.Region("INNER_BOX_VOL").Nodes() <<" perimeter nodes.\n";
    // TESTING WHETHER PRECISION HAS AN IMPACT
    if ( verbose ) {
        const double n_dec_places(1.0e-3);
        for ( auto nit=model_domain.NodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
             (*nit)->x( quantiseToScale( (*nit)->x(), n_dec_places ) );
             (*nit)->y( quantiseToScale( (*nit)->y(), n_dec_places ) );
             (*nit)->z( quantiseToScale( (*nit)->z(), n_dec_places ) );
          }
      }
  
     // nodes at model perimeter
     const csmp::Index    nvar_key(model.Database().StorageKey("nodal variable"));
     multiset<Point<3U> > split_nodes;
     model_domain.InputPropertyValue( "nodal variable", makeScalar(PLAIN,0.) );
     for ( auto nit=model_domain.PerimeterNodesBegin(); nit!=model_domain.NodesEnd(); ++nit ) {
          if ( (*nit)->AtBoundary() != NOT and (*nit)->AtBoundary() != REGION_BOUNDARY )
            (*nit)->Store( nvar_key, makeScalar(PLAIN,2.) );
          else {
               (*nit)->Store( nvar_key, makeScalar(PLAIN,5.) );
               split_nodes.insert( (*nit)->Coordinate() );
            }
       }
  
     if ( verbose ) {
          cerr <<"\nrun: points on ANSYS split boundary:\n";
          size_t counter(0);
          for (auto pt=split_nodes.begin(); pt!=split_nodes.end(); ++pt )
            cerr << fixed << setprecision(3) <<"\n\t"<< counter++ <<": "<< (*pt)[0] <<" "<< (*pt)[1] <<" "<< (*pt)[2];
       }
  
     VTK_Interface<3U>  vtk_output;
     vtk_output.OutputDataToVTK( model, "test-evar", "element variable", 0U );
     vtk_output.OutputDataToVTK( model, "test-nvar", "nodal variable", 0U );
     // output of perimeter elements to VTK
     outputRegionBoundaryToVTK( model, "Model", "model-boundary" );
  

    // 1. create element correspondance set
    // ---------------------------------------------------------------------------------------
      // reports local ids of shared face on either side of the split boundary
      // std::pair<std::pair<Element<3U>*,size_t>,std::pair<Element<3U>*,size_t> >
      set<OppositeElements>  interface_elmt_pairs;
      _test( findSplitInterfaceElements( model.Region("Model"), interface_elmt_pairs ) );
  
      // changing element property values in the discovered elements
      for ( auto eit=interface_elmt_pairs.begin(); eit!=interface_elmt_pairs.end(); ++eit ) {
           (*eit).first.first->Store( evar_key, makeScalar(PLAIN,6.) );
           (*eit).second.first->Store( evar_key, makeScalar(PLAIN,7.) );
        }
      vtk_output.OutputDataToVTK(model, "test-evar", "element variable", 1U );


   // 2. build CSMP SplitBoundary
   // ---------------------------------------------------------------------------------------
   model.DetectAndCreateSplitBoundaries();
   // checking which boundaries were created
   model.SplitBoundariesOut();
   

   // 3. change some property along split boundary to verify that assignments are made correctly
   // ------------------------------------------------------------------------------------------
   SplitBoundary<3U>& splitdomain = (*model.SplitBoundariesBegin()).second;
   // node variable "nodal variable"
   model.InputPropertyValue( "nodal variable", makeScalar(ANY,0.) );
   // inside of boundary
   splitdomain.InputNodePropertyValue( "nodal variable", makeScalar(ANY,1.), INTERIOR, INSIDE );
   splitdomain.InputNodePropertyValue( "nodal variable", makeScalar(ANY,2.), PERIMETER, INSIDE );
   // outside of boundary
   splitdomain.InputNodePropertyValue( "nodal variable", makeScalar(ANY,-1.), INTERIOR, OUTSIDE );
   splitdomain.InputNodePropertyValue( "nodal variable", makeScalar(ANY,-2.), PERIMETER, OUTSIDE );

   // split boundary variable "interface flux"
   const csmp::Index key = model.Database().StorageKey("interface flux");
   splitdomain.Store( key, makeScalar(ANY,1.0e-5) );
   _equal( splitdomain.Read(key), 1.0e-5, numeric_limits<double64>::epsilon() );

   // 4. write altered element properties on either side to VTK
   // ---------------------------------------------------------------------------------------
   vtk_output.OutputDataToVTK( model, "test-nvar", "nodal variable", 2U );

  
   // 5. see whether the split boundary survives being writting to and recovered from file
   // ------------------------------------------------------------------------------------------

} // end run



} // csmp
