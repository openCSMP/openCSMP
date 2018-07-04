#include "ReadingBinaries_Example.h"

#include "Model.h"
#include "Region.h"
#include "PDE_Integrator.h"
#include "CSMP_highLevelUtilities.h"
#include "VSet.h"
#include "VTK_Interface.h"
#include "BinaryFileInterface.h"
#include "CSMP_mathUtilities.h"
#include "ModelTime.h"
#include "FEM_Data.h"

#define DIM 2U

using namespace std;

namespace csmp{

void ReadingBinaries_Example::Specifications()
{
  SetTitle( "BinaryFileInterface: CSMP's native file format." );
  SetDifficulty( 2 );
  SetCategory( "Software Functionality" );
  AddAuthor( "SKM" );
  AddDescription( "source in: ReadingBinaries_Example.cpp" );
  AddDescription( "reading CSMP binary files" );
  AddRequirement( "Binaries from Example SteadyStatePressureToVset");
}


/** *****************************************************************************************

  BinaryFileInterface: Reading of CSMP binary files:
 
   - run example compiled from 'SteadyStatePressureToVset_Example.cpp' 
     to generate binary files, then read these using this program. 
     
     The input files output from this example are:

   - 'model_example2.vset'. It contains the mesh and its connectivity

   - data in CSMP native format to be visualized:
       'fluid-pressure00000001.bin',
       'velocity00000001.bin',
       'volume-flux00000001.bin' (seven zeroes before the 1).
       
   @note Do use PropertyData to store data in VSet and write it to disk in a memory
   efficient fashion because the approach presented in this example leads to padding
   of the integral types, increasing file size by up to 1/3.    
   
   @note Storing variables as FEM_Data inside a VSet has been superseded by the more efficient 
   storage method of PropertyData. There is a unit test PropertyData_Test  that shows 
   how this can be done.

   NB: YOU MUST KNOW THE DIMENSIONS OF THE MODEL YOU WANT TO READ BEFOREHAND, SINCE
      the dimensions of VSET and Model must be known at compile time!
 
***************************************************************************************** */
void ReadingBinaries_Example::Run()
{
  // ESTABLISHING OUTPUTSTREAM FROM BASECLASS
  //ostream &cout = *GetStream();

  cout <<"\nmain: Converter: CSMP binary files to VTK."<< endl;
  cout <<"\nmain: Enter name of connectivity and property database files (filenames without extensions): ";
  string  connectivity_file, physvars_file;
  cin >> connectivity_file >> physvars_file;
  physvars_file += ".txt";

// 1. Reading model geometry from binary VSet file
// -----------------------------------------------------
  VSet<DIM>  saved_model;

  cout << "\nReading model geometry from VSet... " << endl;

  double64& model_time( ModelTime::Instance().modelTime );
  saved_model.InputFrom( connectivity_file.c_str(), model_time );
  //                                                                 isoparametric
  Model<DIM>  binary_file_model( saved_model, physvars_file.c_str(), true );

// 2. Reading data from binary file sets and outputting
//    these as VTK files
// -----------------------------------------------------
  BinaryFileInterface<DIM>  binary_interface;
  VTK_Interface<DIM>        vtk_output;
  string                    data_file;

  // inifinite reading loop: enter 'done' when finished!
  for ( ; ; ) {
       cout <<"\nmain: Enter name of binary (*.bin) data file (without extension); enter 'done' when finished: ";
       cin >> data_file;
       data_file += ".bin";

       // sniffing the variable name in the input file
       string  var = binary_interface.ReadVariableName( data_file.c_str() );
       csmp::Index idx = binary_file_model.Database().StorageKey(var.c_str());

       cout <<"\nmain: reading the '"<< var <<"' data from file '"<< data_file << endl;
       cout.flush();

       // distinguishing scalar, vector and tensor data
       FEM_Data<ScalarVariable>*        scalar_ptr(0);
       FEM_Data<VectorVariable<DIM> >*  vector_ptr(0);
       FEM_Data<TensorVariable<DIM> >*  tensor_ptr(0);
           FEM_Data<ArrayVariable>*         array_ptr(0);
           FEM_Data<FlaggedArrayVariable>*  flagged_array_ptr(0);


           switch( idx.type ) {
                case SCALAR:
                     // create data container
                     scalar_ptr = new FEM_Data<ScalarVariable>();
                     // read the data from the binary file
                     binary_interface.ReadDataFrom( data_file.c_str(), binary_file_model.Database(), *scalar_ptr );
                     // input data to model
                     binary_file_model.InputVariableFrom( var.c_str(), *scalar_ptr );
                     // get rid of the temporary container
                     delete scalar_ptr;
                   break;

            case VECTOR:
                 vector_ptr = new FEM_Data<VectorVariable<DIM> >();
                 binary_interface.ReadDataFrom( data_file.c_str(), binary_file_model.Database(), *vector_ptr );
                 binary_file_model.InputVariableFrom( var.c_str(), *vector_ptr );
                 delete vector_ptr;
               break;

                case TENSOR:
                     tensor_ptr = new FEM_Data<TensorVariable<DIM> >();
                     binary_interface.ReadDataFrom( data_file.c_str(), binary_file_model.Database(), *tensor_ptr );
                     binary_file_model.InputVariableFrom( var.c_str(), *tensor_ptr );
                     delete tensor_ptr;
                   break;

               case ARRAY:
                    // create data container
                    array_ptr = new FEM_Data<ArrayVariable>();
                    // read the data from the binary file
                    binary_interface.ReadDataFrom( data_file.c_str(), binary_file_model.Database(), *array_ptr );
                    // input data to model
                    binary_file_model.InputVariableFrom( var.c_str(), *array_ptr );
                    // get rid of the temporary container
                    delete array_ptr;
                  break;

               case FLAGGEDARRAY:
                    // create data container
                    flagged_array_ptr = new FEM_Data<FlaggedArrayVariable>();
                    // read the data from the binary file
                    binary_interface.ReadDataFrom( data_file.c_str(), binary_file_model.Database(), *flagged_array_ptr );
                    // input data to model
                    binary_file_model.InputVariableFrom( var.c_str(), *flagged_array_ptr );
                    // get rid of the temporary container
                    delete flagged_array_ptr;
                  break;
           }

           // write data to VTK file
           string  output_file(var);
           for ( string::iterator i=output_file.begin(); i!=output_file.end(); i++ ) if ( *i == ' ' ) *i = '_';
           vtk_output.OutputDataToVTK( binary_file_model, output_file.c_str(), var.c_str(), static_cast<size_t>(rint(model_time)) );

       // exit for loop if file name equates to "done"
       if ( data_file == "done" ) break;
    }

  cout <<"\nmain: That's it..."<< endl;

} // Run()

} // csmp
