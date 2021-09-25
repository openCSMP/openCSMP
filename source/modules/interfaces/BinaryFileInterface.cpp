#include "BinaryFileInterface.h"
#include "binaryReadWrite.h"
#include "Exception.h"
#include "CSMP_highLevelUtilities.h"
#include "ModelTime.h"
#include "FEM_Data.h"

using namespace std;

namespace csmp {

template<size_t dim>
BinaryFileInterface<dim>::BinaryFileInterface()
 {
 }


template<size_t dim>
BinaryFileInterface<dim>::~BinaryFileInterface()
 {
 }


/** 
   writes the mesh without property data to binary file.
*/
template<size_t dim>
bool BinaryFileInterface<dim>::WriteConnectivityFile( const Model<dim>& sg, 
                                                      const char* file_name ) const
 {
    double64& model_time( ModelTime::Instance().modelTime );
    VSet<dim>  vset;
    
    sg.Mesh().OutputMeshTo( vset );
    vset.OutputTo( file_name, model_time );
    
    cout <<"\nBinaryFileInterface<"<< dim <<">::WriteConnectivityFile: '";
    cout << file_name <<"' written successfully."<< endl;

    return true;

 } // end WriteConnectivityFile
 
 


template<size_t dim>
bool BinaryFileInterface<dim>::ReadConnectivityFile( const char* file_name, 
                                                     VSet<dim>& vset,
                                                     double64& time ) const
 {
    return vset.InputFrom( file_name, time );

 } // end WriteConnectivityFile

  
  
  
  
  
/**
     Writes property values to disk using FEM_Data containers.
     The timestep is appended only if it has a positive value.
*/
template<size_t dim>
bool BinaryFileInterface<dim>::WriteDataTo( const Model<dim>& sg, 
                                            const char* file_name,
                                            const char* var_name,
                                            long timestep ) const
 {
    VARIABLE_TYPE  vtype = sg.Database().Type( var_name );
    bool           write_error(false);
    string         name( file_name );

    // file name + extension
    if ( timestep > 0 ) {
        char    num[30];
        sprintf( num, "%8ld", timestep );
        string  padded_string( num );
        replaceWhiteSpaceBy( padded_string, '0' );
        name += padded_string;
      }
    name += ".bin";
    
    // 1. opening the file
    fstream fp( name.c_str(), ios::out | ios::binary);
     
     if ( !fp.is_open() ) {
          cerr <<"\nBinaryFileInterface<" << dim;
          cerr <<">::WriteDataTo: File: "<< name;
          cerr <<" could not be opened"<< endl;
          return false;
       }

     if ( vtype == SCALAR ) {
         FEM_Data<ScalarVariable >  scalar_data;
         sg.OutputVariableTo( var_name, scalar_data );
	     if ( scalar_data.Size() != 0 ) {
	          binaryFileWrite( fp, var_name );
	          scalar_data.OutBinary( fp );
	       }
         // writing that there is no such data
         else write_error = true;
       }
     else if ( vtype == VECTOR ) {
         FEM_Data<VectorVariable<dim> > vector_data;
         sg.OutputVariableTo( var_name, vector_data );
	     if ( vector_data.Size() != 0 ) {
	          binaryFileWrite( fp, var_name );
	          vector_data.OutBinary( fp );   
	       }
         else write_error = true;
       }
     else if ( vtype == TENSOR ) {
         FEM_Data<TensorVariable<dim> > tensor_data;
         sg.OutputVariableTo( var_name, tensor_data );
	     if ( tensor_data.Size() != 0 ) {
	          binaryFileWrite( fp, var_name );
	          tensor_data.OutBinary( fp );   
	       }
         else write_error = true;
       }
     else if ( vtype == ARRAY ) {
         FEM_Data<ArrayVariable > array_data;
         sg.OutputVariableTo( var_name, array_data );
         if ( array_data.Size() != 0 ) {
              binaryFileWrite( fp, var_name );
              array_data.OutBinary( fp );
           }
         else write_error = true;
       }
     else if ( vtype == FLAGGEDARRAY ) {
         FEM_Data<FlaggedArrayVariable > flaggedarray_data;
         sg.OutputVariableTo( var_name, flaggedarray_data );
         if ( flaggedarray_data.Size() != 0 ) {
              binaryFileWrite( fp, var_name );
              flaggedarray_data.OutBinary( fp );
           }
         else write_error = true;
       }

	 fp.close();

    if ( write_error ) {
         size_t  n0(0);
         fp.write( (char*) &n0, sizeof(size_t) );

         throw csmp::Exception( ERROR, "BinaryFileInterface<dim>::WriteDataTo",
                         "Variable record was empty:", var_name );
	     return false;
      }

    cout <<"\nBinaryFileInterface<"<< dim;
    cout <<">::WriteDataTo: '"<< name <<"' done."<< endl;

    return true;
    
 } // end WriteDataTo
 




/**
    reads and returns the name of the variable which is stored in the binary file specified as
    method argument.
*/
template<size_t dim>
string  BinaryFileInterface<dim>::ReadVariableName( const char* file_name ) const
 {
     // 1. opening the file     
     string name(file_name);
     name +=".bin";

	 fstream fp(name.c_str(), ios::in | ios::binary);
     if ( !fp.is_open() )
          throw csmp::Exception( FATAL_ERROR, "BinaryFileInterface<dim>::ReadDataFrom: file: ",
                                 file_name, " could not be opened" );
   
     char  variable[NAME_STRING] = "undefined variable";

     // 2. reading the file header
     binaryFileRead( fp, variable );

     // 3. cleaning up
	   fp.close();

     return string(variable);

 } // end ReadVariableName                  





template<size_t dim>
template<typename Var>
string  BinaryFileInterface<dim>::ReadDataFrom( const char* file_name,
                                                const PropertyDatabase<dim>& pref, 
                                                FEM_Data<Var>& var_data ) const
{
     // 1. opening the file
     string name(file_name);
     name +=".bin";
  
	 fstream fp(name.c_str(), ios::in | ios::binary);
     if ( !fp.is_open() ) {
          cerr <<"\nBinaryFileInterface<" << dim;
          cerr <<">::ReadDataFrom<Var>: File: "<< name;
          cerr <<" could not be opened"<< endl;
          return string("undefined variable");
       }
       
     // verifying the existance of the variable in the database
     char  variable[200] = "undefined variable"; 
     binaryFileRead( fp, variable );
     if ( !pref.IsDefined(variable) ) {
          throw csmp::Exception( ERROR, "BinaryFileInterface<dim>::ReadDataFrom<Var>",
                         "Variable is not defined in model database:", variable );
		  fp.close();
          return string(variable);
       }
       
     /// @todo (2-D) Use template here!
     // verifying that the variable has the correct type and placement
     VARIABLE_TYPE  input_var_type(SCALAR);
     if      ( typeid(Var) == typeid(ScalarVariable) )       input_var_type = SCALAR;
     else if ( typeid(Var) == typeid(VectorVariable<dim>) )  input_var_type = VECTOR;
     else if ( typeid(Var) == typeid(TensorVariable<dim>) )  input_var_type = TENSOR;
     else if ( typeid(Var) == typeid(ArrayVariable) )        input_var_type = ARRAY;
     else if ( typeid(Var) == typeid(FlaggedArrayVariable) ) input_var_type = FLAGGEDARRAY;

     if ( input_var_type != pref.Type(variable) ) {
          throw csmp::Exception( ERROR, "BinaryFileInterface<dim>::ReadDataFrom<Var>",
                         "Supplied data object has different type than data" );
		  fp.close();
          return string(variable);
       }
       
     // 4. reading the Vdata to the file   
     // ------------------------------------
     var_data.InBinary( fp );

     // 5. cleaning up
	 fp.close();

   return string(variable);

} // end ReadDataFrom


template class BinaryFileInterface<1U>;
template class BinaryFileInterface<2U>;
template class BinaryFileInterface<3U>;

template string BinaryFileInterface<1U>::ReadDataFrom<ScalarVariable >( const char*,
                                                                        const PropertyDatabase<1>&, 
                                                                        FEM_Data<ScalarVariable>& ) const;


template string BinaryFileInterface<2U>::ReadDataFrom<ScalarVariable >( const char*,
                                                                        const PropertyDatabase<2>&, 
                                                                        FEM_Data<ScalarVariable>& ) const;


template string BinaryFileInterface<3U>::ReadDataFrom<ScalarVariable >( const char*,
                                                                        const PropertyDatabase<3>&, 
                                                                        FEM_Data<ScalarVariable>& ) const;

template string BinaryFileInterface<1U>::ReadDataFrom<VectorVariable<1U> >( const char*,
                                                                            const PropertyDatabase<1>&, 
                                                                            FEM_Data<VectorVariable<1U> >& ) const;
                                            
template string BinaryFileInterface<2U>::ReadDataFrom<VectorVariable<2U> >( const char*,
                                                                            const PropertyDatabase<2>&, 
                                                                            FEM_Data<VectorVariable<2U> >& ) const;
                                            
template string BinaryFileInterface<1U>::ReadDataFrom<TensorVariable<1U> >( const char*,
                                                                            const PropertyDatabase<1>&, 
                                                                            FEM_Data<TensorVariable<1U> >& ) const;

template string BinaryFileInterface<2U>::ReadDataFrom<TensorVariable<2U> >( const char*,
                                                                            const PropertyDatabase<2>&, 
                                                                            FEM_Data<TensorVariable<2U> >& ) const;

template string BinaryFileInterface<3U>::ReadDataFrom<VectorVariable<3U> >( const char*,
                                                                            const PropertyDatabase<3>&, 
                                                                            FEM_Data<VectorVariable<3U> >& ) const;
                                            
template string BinaryFileInterface<3U>::ReadDataFrom<TensorVariable<3U> >( const char*,
                                                                            const PropertyDatabase<3>&, 
                                                                            FEM_Data<TensorVariable<3U> >& ) const;

template string BinaryFileInterface<1U>::ReadDataFrom<ArrayVariable >( const char*,
                                                                        const PropertyDatabase<1>&,
                                                                        FEM_Data<ArrayVariable>& ) const;


template string BinaryFileInterface<2U>::ReadDataFrom<ArrayVariable >( const char*,
                                                                        const PropertyDatabase<2>&,
                                                                        FEM_Data<ArrayVariable>& ) const;


template string BinaryFileInterface<3U>::ReadDataFrom<ArrayVariable >( const char*,
                                                                        const PropertyDatabase<3>&,
                                                                        FEM_Data<ArrayVariable>& ) const;

template string BinaryFileInterface<1U>::ReadDataFrom<FlaggedArrayVariable >( const char*,
                                                                        const PropertyDatabase<1>&,
                                                                        FEM_Data<FlaggedArrayVariable>& ) const;


template string BinaryFileInterface<2U>::ReadDataFrom<FlaggedArrayVariable >( const char*,
                                                                        const PropertyDatabase<2>&,
                                                                        FEM_Data<FlaggedArrayVariable>& ) const;


template string BinaryFileInterface<3U>::ReadDataFrom<FlaggedArrayVariable >( const char*,
                                                                        const PropertyDatabase<3>&,
                                                                        FEM_Data<FlaggedArrayVariable>& ) const;

} // end namespace csmp






















