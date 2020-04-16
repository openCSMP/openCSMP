#include "StatisticalAnalyzer_Example.h"

#include "CSMP_definitions.h"
#include "ANSYS_Interface.h"
#include "ModelTopology.h"
#include "VSet.h"
#include "StatisticalAnalyzer.h"
#include "Model.h"
#include "Region.h"

using namespace std;

namespace csmp{

void StatisticalAnalyzer_Example::Specifications()
{
  SetTitle( "StatisticalAnalyzer: creation of flow velocity histograms" );
  SetDifficulty( 2 );
  SetCategory( "Software Functionality" );
  AddAuthor( "Shaho" );
  AddDescription( "source in: StatisticalAnalyzer_Example.cpp" );
  AddDescription( "generate histograms with csmp::StatisticalAnalyzer" );
  AddDescription( "source file in: StatisticalAnalyzer_Example.cpp" );
  AddRequirement( "LeftRight .asc,.dat,-regions.txt" );
  AddRequirement( "StatisticalAnalyzer_Example_var.txt" );
}



/**
    Reads ANSYS 2D model and analyzes statistical properties.
*/
void StatisticalAnalyzer_Example::Run()
{
  const bool          isoparametric(true);
  ANSYS_Interface     mesh_interface(isoparametric);
  VSet<2U>            mesh_container;
  ModelTopology       mesh_topology(isoparametric);


  string model_name("undefined");
  cout << "Enter name of model: ";
  cin >> model_name;

  const bool binary_file( true );
  const bool irregular_mesh( false );
  mesh_interface.Read_ANSYS_Mesh( model_name.c_str(), mesh_container, mesh_topology, binary_file, irregular_mesh );

  Model<2U> reservoir_model( mesh_topology, mesh_container, "StatisticalAnalyzer_Example_var.txt" );

  Region<2>& regionref = reservoir_model.Region("Model");

  Index porosity_key = reservoir_model.Database().StorageKey("porosity");

  for ( vector<Element<2U>*>::iterator
    it=regionref.ElementsBegin(); it!=regionref.ElementsEnd(); it++ )
  {
      double64 porosity = NormalDistributionGenerator(0.5, 0.1, 0, 1);
      (*it)->Store( porosity_key, makeScalar(PLAIN, porosity) );
  }


  StatisticalAnalyzer<2U>                  porosity_histogram(reservoir_model);
  HistogramBins                            bins;
  map<string,pair<HistogramBins,size_t> >  porosity_results;

  /*
  option 1
  this method reads a histogram bin file with the following format

  'myfile.bins'
  8
  1.0e-5
  2.0e-5
  4.0e-5
  6.0e-5
  8.0e-5
  1.0e-4
  2.0e-4
  6.0e-4

  flux_histogram.DefineBins( "porosity_histogram.bins", bins );
  */

  // option 2
  porosity_histogram.DefineBins( 0., 1., 30, bins );

  //                                               input property
  porosity_histogram.RegionPropertyHistogramsElement( "porosity", bins, porosity_results );
  //                                                             output file name                   log10_of_bin_values
  porosity_histogram.OutputRegionPropertyAbundancePolygonsMaple( "porosity", bins, porosity_results, false );

} // Run()






/**
     Generates a normal distribution of variable values.
*/
double64 StatisticalAnalyzer_Example::NormalDistributionGenerator( double64 mean, double64 sd,
                                                                  double64 minimum, double64 maximum )
{
  static unsigned int counter( 1U );
  double64 result( 0. );

  do
  {
      double64 x1 ( (double64)rand() / (double64)RAND_MAX );
      double64 x2 ( (double64)rand() / (double64)RAND_MAX );

      if ( ( counter % 2 ) == 0 )
          result = sqrt(-2. * log(x1)) * sin(2 * 3.14159265 * x2);
      else
          result = sqrt(-2. * log(x1)) * cos(2 * 3.14159265 * x2);

      result = mean + result * sd;
      counter++;

  } while ( ( result < minimum ) || ( result > maximum ) );

  return result;
}

} // csmp
