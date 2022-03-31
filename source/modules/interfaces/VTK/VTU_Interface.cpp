#include "VTU_Interface.h"
#include "VTK_Interface.h"

#include "Boundary.h"
#include "SplitBoundary.h"
#include "Region.h"
#include "Model.h"
#include "OS_Utilities.h"

#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/**
    constructor taking a reference to the model to be handled(and an optional title to the output)
     
    @todo introduce option to choose XML type:
    
    currently:  "<VTKFile type="UnstructuredGrid" version="0.9" byte_order="LittleEndian">"
    
    but one should be able to switch to 1.0 etc., however this is currently not possible because the XML
    version is not a state variable of the VTU interface.
*/
template<uint32_t dim>
VTU_Interface<dim>::VTU_Interface( const Model<dim>& model,
                                   const string& problemTitle,
                                   bool use_problem_title_as_output_folder_name )
    : model_        ( model ),
      problemTitle_ ( problemTitle ),
      toFolder_     ( use_problem_title_as_output_folder_name ),
      toSubFolder_  ( false ),
      variableNameAliases_(),
      omitZeroInFileName_ ( true ),
      elementVecAndTensDataAtCellCenters_(false),
      regionVecAndTensDataAtCellCenters_ (false),
      suffix_text_("")
{
}

/// constructor taking a reference to the model to be handled(and an optional title to the output)
template<uint32_t dim>
VTU_Interface<dim>::VTU_Interface( const Model<dim>& model,
                                   const string& problemTitle,
                                   const string& subFolderName,
                                   bool use_problem_title_as_output_folder_name )
 : model_           ( model ),
   problemTitle_    ( problemTitle ),
   subFolderName_   ( subFolderName ),
   toFolder_        ( use_problem_title_as_output_folder_name ),
   toSubFolder_     ( true ),
   omitZeroInFileName_( true ),
   elementVecAndTensDataAtCellCenters_(false),
   regionVecAndTensDataAtCellCenters_ (false),
   variableNameAliases_(),
   suffix_text_("")
{
}

/// destructor removing connectivity files from the heap
template<uint32_t dim>
VTU_Interface<dim>::~VTU_Interface()
{
  // loops over xml document entries and calls delete
  DeleteConnectivity();
}

// CONNECTIVITY FILES
// -------------------

/// allows to force recreation of connectivity map (FOR THE UPCOMING OUTPUT ONLY!!!)
/**
    @todo SKM verify this deletion process; looks like that there may be a memory leak.
*/
template<uint32_t dim>
void VTU_Interface<dim>::DeleteConnectivity()
  {
    /// delete connectivity files for Field, Node and Element data
    for( typename map<const ModelSubDomain<dim,Element>*,XML_Document*>::iterator it = regionConnectivityFiles_.begin(); it != regionConnectivityFiles_.end(); ++it )
      if( it->second )
        delete it->second;
    regionConnectivityFiles_.clear();
    for( typename map<const ModelSubDomain<dim,Face>*,XML_Document*>::iterator it = boundaryConnectivityFiles_.begin(); it != boundaryConnectivityFiles_.end(); ++it )
      if( it->second )
        delete it->second;
    boundaryConnectivityFiles_.clear();
    for( typename map<const ModelSubDomain<dim,InterFace>*,XML_Document*>::iterator it = splitBoundaryConnectivityFiles_.begin(); it != splitBoundaryConnectivityFiles_.end(); ++it )
      if( it->second )
        delete it->second;
    splitBoundaryConnectivityFiles_.clear();

    /// delete connectivity files for Element BaryCenters data
    for( typename map<const ModelSubDomain<dim,Element>*,XML_Document*>::iterator it = regionConnectivityFiles_bcd_.begin(); it != regionConnectivityFiles_bcd_.end(); ++it )
      if( it->second )
        delete it->second;
    regionConnectivityFiles_bcd_.clear();
    for( typename map<const ModelSubDomain<dim,Face>*,XML_Document*>::iterator it = boundaryConnectivityFiles_bcd_.begin(); it != boundaryConnectivityFiles_bcd_.end(); ++it )
      if( it->second )
        delete it->second;
    boundaryConnectivityFiles_bcd_.clear();
    for( typename map<const ModelSubDomain<dim,InterFace>*,XML_Document*>::iterator it = splitBoundaryConnectivityFiles_bcd_.begin(); it != splitBoundaryConnectivityFiles_bcd_.end(); ++it )
      if( it->second )
        delete it->second;
    splitBoundaryConnectivityFiles_bcd_.clear();

    /// delete connectivity files for Regions data
    for( typename map<const ModelSubDomain<dim,Element>*,XML_Document*>::iterator it = regionConnectivityFiles_rcd_.begin(); it != regionConnectivityFiles_rcd_.end(); ++it )
      if( it->second )
        delete it->second;
    regionConnectivityFiles_rcd_.clear();
    for( typename map<const ModelSubDomain<dim,Face>*,XML_Document*>::iterator it = boundaryConnectivityFiles_rcd_.begin(); it != boundaryConnectivityFiles_rcd_.end(); ++it )
      if( it->second )
        delete it->second;
    boundaryConnectivityFiles_rcd_.clear();
    for( typename map<const ModelSubDomain<dim,InterFace>*,XML_Document*>::iterator it = splitBoundaryConnectivityFiles_rcd_.begin(); it != splitBoundaryConnectivityFiles_rcd_.end(); ++it )
      if( it->second )
        delete it->second;
    splitBoundaryConnectivityFiles_rcd_.clear();

    /// delete connectivity files for Finite Element Integration Points data
    for( typename map<const ModelSubDomain<dim,Element>*,XML_Document*>::iterator it = regionConnectivityFiles_feipsd_.begin(); it != regionConnectivityFiles_feipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    regionConnectivityFiles_feipsd_.clear();
    for( typename map<const ModelSubDomain<dim,Face>*,XML_Document*>::iterator it = boundaryConnectivityFiles_feipsd_.begin(); it != boundaryConnectivityFiles_feipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    boundaryConnectivityFiles_feipsd_.clear();
    for( typename map<const ModelSubDomain<dim,InterFace>*,XML_Document*>::iterator it = splitBoundaryConnectivityFiles_feipsd_.begin(); it != splitBoundaryConnectivityFiles_feipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    splitBoundaryConnectivityFiles_feipsd_.clear();

    /// delete connectivity files for Finite Volume Sector Integration Points data
    for( typename map<const ModelSubDomain<dim,Element>*,XML_Document*>::iterator it = regionConnectivityFiles_fvsipsd_.begin(); it != regionConnectivityFiles_fvsipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    regionConnectivityFiles_fvsipsd_.clear();
    for( typename map<const ModelSubDomain<dim,Face>*,XML_Document*>::iterator it = boundaryConnectivityFiles_fvsipsd_.begin(); it != boundaryConnectivityFiles_fvsipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    boundaryConnectivityFiles_fvsipsd_.clear();
    for( typename map<const ModelSubDomain<dim,InterFace>*,XML_Document*>::iterator it = splitBoundaryConnectivityFiles_fvsipsd_.begin(); it != splitBoundaryConnectivityFiles_fvsipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    splitBoundaryConnectivityFiles_fvsipsd_.clear();

    /// delete connectivity files for Finite Volume Facet Integration Points data
    for( typename map<const ModelSubDomain<dim,Element>*,XML_Document*>::iterator it = regionConnectivityFiles_fvfipsd_.begin(); it != regionConnectivityFiles_fvfipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    regionConnectivityFiles_fvfipsd_.clear();
    for( typename map<const ModelSubDomain<dim,Face>*,XML_Document*>::iterator it = boundaryConnectivityFiles_fvfipsd_.begin(); it != boundaryConnectivityFiles_fvfipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    boundaryConnectivityFiles_fvfipsd_.clear();
    for( typename map<const ModelSubDomain<dim,InterFace>*,XML_Document*>::iterator it = splitBoundaryConnectivityFiles_fvfipsd_.begin(); it != splitBoundaryConnectivityFiles_fvfipsd_.end(); ++it )
      if( it->second )
        delete it->second;
    splitBoundaryConnectivityFiles_fvfipsd_.clear();

    /// delete connectivity files of multiblock dataset
    for( typename map<const ModelSubDomain<dim,Element>*,XML_Document*>::iterator it = regionConnectivityFiles_multiblock_.begin(); it != regionConnectivityFiles_multiblock_.end(); ++it )
      if( it->second )
        delete it->second;
    regionConnectivityFiles_multiblock_.clear();
    for( typename map<const ModelSubDomain<dim,Face>*,XML_Document*>::iterator it = boundaryConnectivityFiles_multiblock_.begin(); it != boundaryConnectivityFiles_multiblock_.end(); ++it )
      if( it->second )
        delete it->second;
    boundaryConnectivityFiles_multiblock_.clear();
    for( typename map<const ModelSubDomain<dim,InterFace>*,XML_Document*>::iterator it = splitBoundaryConnectivityFiles_multiblock_.begin(); it != splitBoundaryConnectivityFiles_multiblock_.end(); ++it )
      if( it->second )
        delete it->second;
    splitBoundaryConnectivityFiles_multiblock_.clear();
  }


template<uint32_t dim>
map<const ModelSubDomain<dim,Element>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMap( const ModelSubDomain<dim,Element>&   )
{
    return regionConnectivityFiles_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Face>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMap( const ModelSubDomain<dim,Face>&   )
{
    return boundaryConnectivityFiles_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMap( const ModelSubDomain<dim,InterFace>&   )
{
    return splitBoundaryConnectivityFiles_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Element>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapBCS( const ModelSubDomain<dim,Element>&   )
{
    return regionConnectivityFiles_bcd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Face>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapBCS( const ModelSubDomain<dim,Face>&  )
{
    return boundaryConnectivityFiles_bcd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapBCS( const ModelSubDomain<dim,InterFace>&  )
{
    return splitBoundaryConnectivityFiles_bcd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Element>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapRCS( const ModelSubDomain<dim,Element>&  )
{
    return regionConnectivityFiles_rcd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Face>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapRCS( const ModelSubDomain<dim,Face>&  )
{
    return boundaryConnectivityFiles_rcd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapRCS( const ModelSubDomain<dim,InterFace>&  )
{
    return splitBoundaryConnectivityFiles_rcd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Element>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFEIPS( const ModelSubDomain<dim,Element>&  )
{
    return regionConnectivityFiles_feipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Face>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFEIPS( const ModelSubDomain<dim,Face>&  )
{
    return boundaryConnectivityFiles_feipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFEIPS( const ModelSubDomain<dim,InterFace>&  )
{
    return splitBoundaryConnectivityFiles_feipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Element>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFVSIPS( const ModelSubDomain<dim,Element>&  )
{
    return regionConnectivityFiles_fvsipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Face>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFVSIPS( const ModelSubDomain<dim,Face>&  )
{
    return boundaryConnectivityFiles_fvsipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFVSIPS( const ModelSubDomain<dim,InterFace>&  )
{
    return splitBoundaryConnectivityFiles_fvsipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Element>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFVFIPS( const ModelSubDomain<dim,Element>&  )
{
    return regionConnectivityFiles_fvfipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Face>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFVFIPS( const ModelSubDomain<dim,Face>&  )
{
    return boundaryConnectivityFiles_fvfipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapFVFIPS( const ModelSubDomain<dim,InterFace>&  )
{
    return splitBoundaryConnectivityFiles_fvfipsd_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Element>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapMultiBlock( const ModelSubDomain<dim,Element>&  )
{
    return regionConnectivityFiles_multiblock_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,Face>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapMultiBlock( const ModelSubDomain<dim,Face>&  )
{
    return boundaryConnectivityFiles_multiblock_;
}
template<uint32_t dim>
map<const ModelSubDomain<dim,InterFace>*,XML_Document*>& VTU_Interface<dim>::GetConnectivityMapMultiBlock( const ModelSubDomain<dim,InterFace>&  )
{
    return splitBoundaryConnectivityFiles_multiblock_;
}

// INTERFACES
// ----------------

/// provide aliases for csmp variables that will be shown in vtu instead
/// returns if existing was empty
template<uint32_t dim>
bool VTU_Interface<dim>::VariableNameAliases( const map<string,string>& variableNameAliases )
{
  bool variableNameAliasesWasEmpty( variableNameAliases_.empty() );
  variableNameAliases_ = variableNameAliases;
  return variableNameAliasesWasEmpty;
}

/// for a given csmp variable, returns output alias if existing. returns variable name otherwise
template<uint32_t dim>
string VTU_Interface<dim>::FindVariableOutputAlias( const string& csmpVariableName ) const
{
  // default is equal to variable name
  string outputAliasName( csmpVariableName );
  // if alias exists in map we overwrite
  if( !variableNameAliases_.empty() )
    mapFind( csmpVariableName, outputAliasName, variableNameAliases_ );
  // returned in any case
  return outputAliasName;
}

/// allows the creation of an index-to-name correspondance.  If the resulting container is
/// empty or the array variable name does not exist upon the output call, the behaviour is the detault one (output of [#] suffix).
/// If the component name exists,it gets prefixed to the variable name for the vtu file. Works for flagged arrays, too.
/// Julian, July 2014
template<uint32_t dim>
void VTU_Interface<dim>::CreateArrayComponentPrefixNames(string array_var_name, vector<string>& index_to_name )
{
    this->array_index_to_name_.insert(make_pair(array_var_name,index_to_name));
}

/// Allows adding text after the timestep number to the vtu filename. By default, this string is empty.
template<uint32_t dim>
void VTU_Interface<dim>::SetSuffixText(string text)
{
    this->suffix_text_=text;
}

template<uint32_t dim>
const string& VTU_Interface<dim>::GetProblemTitle( ) const
{
    return this->problemTitle_;
}

/// control whether '0' is appended
template<uint32_t dim>
void VTU_Interface<dim>::OmitZeroInFileName( bool omitZeroInFileName )
{
    omitZeroInFileName_ = omitZeroInFileName;
}
template<uint32_t dim>
bool VTU_Interface<dim>::OmitZeroInFileName() const
{
    return omitZeroInFileName_;
}
template<uint32_t dim>
void VTU_Interface<dim>::OutputElementVectorAndTensorDataAtCellCenters( bool flag )
{
   elementVecAndTensDataAtCellCenters_ = flag;
}

template<uint32_t dim>
void VTU_Interface<dim>::OutputRegionVectorAndTensorDataAtRegionCenters( bool flag )
{
   regionVecAndTensDataAtCellCenters_ = flag;
}

template<uint32_t dim>
string VTU_Interface<dim>::OutputFileNamePrefix( const string& fileName )
{
    if( toSubFolder_ || toFolder_ ){
        string output_file_name;
        output_file_name  = currentDirectorySymbol();
        if( toFolder_ ){
            output_file_name += problemTitle_;
            createDirectoryIfDoesntExist( output_file_name );
            output_file_name += directorySymbol();
            if( toSubFolder_ ){
                output_file_name += subFolderName_;
                createDirectoryIfDoesntExist( output_file_name );
                output_file_name += directorySymbol();
            }
        }else{
            output_file_name += subFolderName_;
            createDirectoryIfDoesntExist( output_file_name );
            output_file_name += directorySymbol();
        }
        if( fileName.find( output_file_name ) == string::npos ){
           output_file_name += fileName;
           return output_file_name;
        }
    }
    return fileName;
}

template<uint32_t dim>
string VTU_Interface<dim>::DomainSpecificOutputFileNamePrefix( const string& fileName, const string& domainName )
{
    // write to file
    string outputName( fileName );
    // SKM_FIX: if Model is not contained in the file name
    if ( domainName.find("Model") == string::npos )
    //  if( regionNameString != "Model" )
    {
        outputName += "_";
        outputName += domainName;
    }

    if( this->suffix_text_ != string("") )
    {
        outputName += "_";
        outputName += this->suffix_text_;
    }

    replaceWhiteSpaceBy( outputName, '_' );

    return outputName;
}

template<uint32_t dim>
template<class T>
string VTU_Interface<dim>::FullOutputFileName( const string& fileName, T timestep )
{
    // write to file
    string outputName(fileName);
    string numberString;
    if( !omitZeroInFileName_ || timestep != 0 )
    {
        outputName += "_";
        numberString = number_to_string( timestep );
        outputName += numberString;
    }
    return outputName;
}

template string VTU_Interface<1U>::FullOutputFileName(const string&,int);
template string VTU_Interface<2U>::FullOutputFileName(const string&,int);
template string VTU_Interface<3U>::FullOutputFileName(const string&,int);

template string VTU_Interface<1U>::FullOutputFileName(const string&,long);
template string VTU_Interface<2U>::FullOutputFileName(const string&,long);
template string VTU_Interface<3U>::FullOutputFileName(const string&,long);

template string VTU_Interface<1U>::FullOutputFileName(const string&,size_t);
template string VTU_Interface<2U>::FullOutputFileName(const string&,size_t);
template string VTU_Interface<3U>::FullOutputFileName(const string&,size_t);

template string VTU_Interface<1U>::FullOutputFileName(const string&,double);
template string VTU_Interface<2U>::FullOutputFileName(const string&,double);
template string VTU_Interface<3U>::FullOutputFileName(const string&,double);




// CSMP DATA OUTPUT
// ----------------

namespace vtuInterfaceDispatch
{
    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<3U,int>( const csmp::Model<3U>& model, const string& propertyName, const string& subDomainName,  int timeStep )
    {
        VTK_Interface<3>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<2U,int>( const csmp::Model<2U>& model, const string& subDomainName, const string& propertyName,  int timeStep )
    {
        VTK_Interface<2>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<1U,int>( const csmp::Model<1U>&, const string&, const string&,  int )
    {
        throw csmp::Exception( ERROR, "VTU_Interface<dim>::DispatchIntegrationPointPropertyOutputToVTK_Interface", "1D not supported" );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<3U,long>( const csmp::Model<3U>& model, const string& propertyName, const string& subDomainName,  long timeStep )
    {
        VTK_Interface<3>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<2U,long>( const csmp::Model<2U>& model, const string& subDomainName, const string& propertyName,  long timeStep )
    {
        VTK_Interface<2>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<1U,long>( const csmp::Model<1U>&, const string&, const string&,  long )
    {
        throw csmp::Exception( ERROR, "VTU_Interface<dim>::DispatchIntegrationPointPropertyOutputToVTK_Interface", "1D not supported" );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<3U,size_t>( const csmp::Model<3U>& model, const string& propertyName, const string& subDomainName,  size_t timeStep )
    {
        VTK_Interface<3>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<2U,size_t>( const csmp::Model<2U>& model, const string& subDomainName, const string& propertyName,  size_t timeStep )
    {
        VTK_Interface<2>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<1U,size_t>( const csmp::Model<1U>&, const string&, const string&,  size_t )
    {
        throw csmp::Exception( ERROR, "VTU_Interface<dim>::DispatchIntegrationPointPropertyOutputToVTK_Interface", "1D not supported" );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<3U,double>( const csmp::Model<3U>& model, const string& propertyName, const string& subDomainName,  double timeStep )
    {
        VTK_Interface<3>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<2U,double>( const csmp::Model<2U>& model, const string& subDomainName, const string& propertyName,  double timeStep )
    {
        VTK_Interface<2>().OutputDataToVTK( model, subDomainName, propertyName, propertyName, timeStep );
    }

    template<>
    void dispatchCompileTimeToRuntimeVTK_Output<1U,double>( const csmp::Model<1U>&, const string&, const string&,  double )
    {
        throw csmp::Exception( ERROR, "VTU_Interface<dim>::DispatchIntegrationPointPropertyOutputToVTK_Interface", "1D not supported" );
    }

} // vtuInterfaceDispatch

/// user interface to output single property to vtu for a csmp::Region
template<uint32_t dim>
template<class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& fileName,
                                          const string& propertyName,
                                          const string& regionName,
                                          T timestep )
{
  list<string> singlePropertyList;
  singlePropertyList.push_back( propertyName );
  return OutputDataToVTU( fileName, singlePropertyList, regionName, timestep );
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const string&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const string&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const string&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const string&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const string&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const string&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const string&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const string&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const string&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const string&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const string&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const string&,double);



/// user interface to output list of properties to vtu for a csmp::Region
template<uint32_t dim>
template<class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& fileName,
                                          const list<string>& propertyNames,
                                          const string& regionName,
                                          T timestep )
{
  const Region<dim>& rref( model_.Region( regionName ) );
  return OutputDataToVTU( fileName, propertyNames, rref, timestep );
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const string&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const string&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const string&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const string&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const string&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const string&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const string&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const string&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const string&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const string&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const string&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const string&,double);

/// user interface to output vector of properties to vtu for a csmp::Region
template<uint32_t dim>
template<class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& fileName,
                                          const vector<string>& propertyNames,
                                          const string& regionName,
                                          T timestep )
{
  const Region<dim>& rref( model_.Region( regionName ) );
  return OutputDataToVTU( fileName, propertyNames, rref, timestep );
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const string&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const string&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const string&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const string&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const string&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const string&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const string&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const string&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const string&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const string&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const string&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const string&,double);

/// user interface to output set of properties to vtu for a csmp::Region
template<uint32_t dim>
template<class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& fileName,
                                          const set<string>& propertyNames,
                                          const string& regionName,
                                          T timestep )
{
  const Region<dim>& rref( model_.Region( regionName ) );
  return OutputDataToVTU( fileName, propertyNames, rref, timestep );
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const string&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const string&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const string&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const string&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const string&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const string&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const string&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const string&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const string&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const string&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const string&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const string&,double);

/// user interface to output property to vtu for a csmp::ModelSubDomain
template<uint32_t dim>
template<template <uint32_t> class CELL,class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& fileName,
                                          const string& propertyName,
                                          const ModelSubDomain<dim,CELL>& subDomain,
                                          T timestep )
{
  list<string> singlePropertyList;
  singlePropertyList.push_back( propertyName );
  return OutputDataToVTU( fileName, singlePropertyList, subDomain, timestep );
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Element>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Element>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Element>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Face>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Face>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Face>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,InterFace>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,InterFace>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,InterFace>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Element>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Element>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Element>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Face>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Face>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Face>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,InterFace>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,InterFace>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,InterFace>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Element>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Element>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Element>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Face>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Face>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Face>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,InterFace>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,InterFace>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,InterFace>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Element>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Element>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Element>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,Face>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,Face>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,Face>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<1U,InterFace>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<2U,InterFace>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const string&,const ModelSubDomain<3U,InterFace>&,double);

/// user interface to output list of properties to vtu for a csmp::ModelSubDomain
template<uint32_t dim>
template<template <uint32_t> class CELL,class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& initial_file_name,
                                          const list<string>& propertyNamesList,
                                          const ModelSubDomain<dim,CELL>& subDomain,
                                          T timestep )
{
    // avoid using duplicated properties
    set<string> propertyNames;
    for( typename list<string>::const_iterator it = propertyNamesList.begin(); it != propertyNamesList.end(); ++it)
        propertyNames.insert(*it);
    return OutputDataToVTU( initial_file_name, propertyNames, subDomain, timestep );
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Element>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Element>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Element>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Face>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Face>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Face>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,InterFace>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,InterFace>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,InterFace>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Element>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Element>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Element>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Face>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Face>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Face>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,InterFace>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,InterFace>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,InterFace>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Element>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Element>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Element>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Face>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Face>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Face>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,InterFace>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,InterFace>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,InterFace>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Element>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Element>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Element>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,Face>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,Face>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,Face>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<1U,InterFace>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<2U,InterFace>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const list<string>&,const ModelSubDomain<3U,InterFace>&,double);


/// user interface to output vector of properties to vtu for a csmp::ModelSubDomain
template<uint32_t dim>
template<template <uint32_t> class CELL,class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& initial_file_name,
                                          const vector<string>& propertyNamesList,
                                          const ModelSubDomain<dim,CELL>& subDomain,
                                          T timestep )
{
    // avoid using duplicated properties
    set<string> propertyNames;
    for( typename vector<string>::const_iterator it = propertyNamesList.begin(); it != propertyNamesList.end(); ++it)
        propertyNames.insert(*it);
    return OutputDataToVTU( initial_file_name, propertyNames, subDomain, timestep );
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Element>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Element>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Element>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Face>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Face>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Face>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,InterFace>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,InterFace>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,InterFace>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Element>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Element>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Element>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Face>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Face>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Face>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,InterFace>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,InterFace>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,InterFace>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Element>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Element>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Element>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Face>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Face>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Face>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,InterFace>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,InterFace>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,InterFace>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Element>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Element>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Element>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,Face>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,Face>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,Face>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<1U,InterFace>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<2U,InterFace>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const vector<string>&,const ModelSubDomain<3U,InterFace>&,double);


/// user interface to output set of properties to vtu for a csmp::ModelSubDomain
template<uint32_t dim>
template<template <uint32_t> class CELL,class T>
bool VTU_Interface<dim>::OutputDataToVTU( const string& initial_file_name,
                                          const set<string>& propertyNames,
                                          const ModelSubDomain<dim,CELL>& subDomain,
                                          T timestep )
{
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    string fileNamePrefix( this->OutputFileNamePrefix( initial_file_name ) );
    fileNamePrefix = DomainSpecificOutputFileNamePrefix( fileNamePrefix, DomainName(subDomain) );
    string relativeFileNamePrefix( DomainSpecificOutputFileNamePrefix( initial_file_name, DomainName(subDomain) ) );

    // split list into separate ones for node and element placed indices
    list<Index> nodeIndices;
    list<Index> elementIndices;
    list<Index> fieldDataIndices;
    list<Index> elementMatrixIndices;
    list<Index> regionIndices;
    list<Index> feipIndices;
    list<Index> fvsipIndices;
    list<Index> fvfipIndices;
    
    for( set<string>::const_iterator it = propertyNames.begin(); it != propertyNames.end(); ++it )
    {
        const PLACEMENT variablePlacement( model_.Database().Placement( it->c_str() ) );
        // the right part of this if statement was added to be able to output variables placed on the MODEL
        if( !subDomain.ValidVariable( it->c_str() ) && variablePlacement != MODEL)
        {
            if ( csmp_error.Verbose() ) {
                cerr << "\n VTU_Interface<dim>::OutputDataToVTU: Output property '" << *it <<"' ("<< parsePlacement(variablePlacement) <<")";
                cerr << " not accessible in domain "<< DomainName( subDomain ) << endl;
              }
            continue;
        }
        const VARIABLE_TYPE variableType( model_.Database().Type( it->c_str() ) );
        //if( variableType != SCALAR && variableType != VECTOR && variableType != TENSOR && variableType != ARRAY && variableType != FLAGGEDARRAY)
        //    throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", "Output for scalar, vector, tensor, array and flagged array variables only" );
        if ( variablePlacement == MODEL )
            fieldDataIndices.push_back( model_.Database().StorageKey( it->c_str() ));
        else if ( variablePlacement == REGION || variablePlacement == BOUNDARY || variablePlacement == SPLIT_BOUNDARY )
        {
            if( !regionVecAndTensDataAtCellCenters_ )
                fieldDataIndices.push_back( model_.Database().StorageKey( it->c_str() ) );
            else
            {
                switch( variableType )
                {
                case ARRAY:
                case FLAGGEDARRAY:
                case SCALAR:
                    fieldDataIndices.push_back( model_.Database().StorageKey( it->c_str() ) );
                    break;
                case VECTOR:
                case TENSOR:
                    regionIndices.push_back( model_.Database().StorageKey( it->c_str() ) ) ;
                    break;
                default: throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", "Region output for scalar, vector, tensor and array variables only" );
                }
            }
        }
        else if( variablePlacement == NODE )
            nodeIndices.push_back( model_.Database().StorageKey( it->c_str() ) );
        else if( variablePlacement == ELEMENT || variablePlacement == FACE || variablePlacement == INTER_FACE )
        {
            if( !elementVecAndTensDataAtCellCenters_ )
                elementIndices.push_back( model_.Database().StorageKey( it->c_str() ) );
            else
            {
                switch( variableType )
                {
                case ARRAY:
                case FLAGGEDARRAY:
                case SCALAR:
                    elementIndices.push_back( model_.Database().StorageKey( it->c_str() ) );
                    break;
                case VECTOR:
                case TENSOR: elementMatrixIndices.push_back( model_.Database().StorageKey( it->c_str() ) );
                    break;
                default: throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", "Element output for scalar, vector, tensor and array variables only" );
                }
            }
        }
        else if( variablePlacement == ELEMENT_INTEGRATION_POINT || variablePlacement == FACE_INTEGRATION_POINT        || variablePlacement == INTER_FACE_INTEGRATION_POINT )
            feipIndices.push_back(model_.Database().StorageKey( it->c_str() ));
        else if( variablePlacement == SECTOR_INTEGRATION_POINT  || variablePlacement == FACE_SECTOR_INTEGRATION_POINT || variablePlacement == INTER_FACE_SECTOR_INTEGRATION_POINT )
            fvsipIndices.push_back(model_.Database().StorageKey( it->c_str() ));
        else if( variablePlacement == FACET_INTEGRATION_POINT   || variablePlacement == FACE_FACET_INTEGRATION_POINT  || variablePlacement == INTER_FACE_FACET_INTEGRATION_POINT )
            fvfipIndices.push_back(model_.Database().StorageKey( it->c_str() ));
        else{
            //vtuInterfaceDispatch::dispatchCompileTimeToRuntimeVTK_Output( model_, domainName( subDomain ).c_str(), it->c_str(), timestep );
            string errormsg="Placement not supported: "+*it;
            throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", errormsg );
        }
    }

    vector<string> files;
    string outputFileName;
    bool output;

    /// scalars,arrays,vectors,tensors on model, elements and nodes
    output = ( !fieldDataIndices.empty() || !elementIndices.empty() || !nodeIndices.empty() );
    if( output )
    {
        outputFileName  = relativeFileNamePrefix;
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        files.push_back( outputFileName );
        outputFileName  = fileNamePrefix;
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        OutputFieldNodesAndElementDataToVTU( outputFileName, fieldDataIndices, nodeIndices, elementIndices, subDomain );
    }

    /// scalars,arrays,vectors,tensors on element barycenters
    output = !elementMatrixIndices.empty();
    if( output )
    {
        outputFileName  = relativeFileNamePrefix;
        outputFileName += "_bcs";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        files.push_back( outputFileName );
        outputFileName  = fileNamePrefix;
        outputFileName += "_bcs";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        OutputElementBarycentricDataToVTU( outputFileName, elementMatrixIndices, subDomain );
    }

    /// scalars,arrays,vectors,tensors on regions
    output = !regionIndices.empty();
    if( output )
    {
        outputFileName  = relativeFileNamePrefix;
        outputFileName += "_rcs";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        files.push_back( outputFileName );
        outputFileName  = fileNamePrefix;
        outputFileName += "_rcs";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        OutputRegionDataToVTU( outputFileName, regionIndices, subDomain );
    }

    /// scalars,arrays,vectors,tensors on finite element integration points
    output = !feipIndices.empty();
    if( output )
    {
        outputFileName  = relativeFileNamePrefix;
        outputFileName += "_feips";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        files.push_back( outputFileName );
        outputFileName  = fileNamePrefix;
        outputFileName += "_feips";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        OutputFiniteElementIntegrationPointsDataToVTU( outputFileName, feipIndices, subDomain );
    }

    /// scalars,arrays,vectors,tensors on finite volume sector integration points
    output = !fvsipIndices.empty();
    if( output )
    {
        outputFileName  = relativeFileNamePrefix;
        outputFileName += "_fvsips";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        files.push_back( outputFileName );
        outputFileName  = fileNamePrefix;
        outputFileName += "_fvsips";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        OutputFiniteVolumeSectorIntegrationPointsDataToVTU( outputFileName, fvsipIndices, subDomain );
    }

    /// scalars,arrays,vectors,tensors on  finite volume facet integration points
    output = !fvfipIndices.empty();
    if( output )
    {
        outputFileName  = relativeFileNamePrefix;
        outputFileName += "_fvfips";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        files.push_back( outputFileName );
        outputFileName  = fileNamePrefix;
        outputFileName += "_fvfips";
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        OutputFiniteVolumeFacetIntegrationPointsDataToVTU( outputFileName, fvfipIndices, subDomain );
    }

    /// write multiblock file *.vtm
    output = ( files.size() > 1 );
    if ( output )
    {
        outputFileName  = fileNamePrefix;
        outputFileName  = FullOutputFileName( outputFileName, timestep );
        OutputMultiBlockVTU( outputFileName, files, subDomain );
    }

    return true;
}

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Element>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Element>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Element>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Face>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Face>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Face>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,InterFace>&,int);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,InterFace>&,int);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,InterFace>&,int);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Element>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Element>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Element>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Face>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Face>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Face>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,InterFace>&,long);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,InterFace>&,long);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,InterFace>&,long);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Element>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Element>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Element>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Face>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Face>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Face>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,InterFace>&,size_t);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,InterFace>&,size_t);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,InterFace>&,size_t);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Element>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Element>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Element>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,Face>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,Face>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,Face>&,double);

template bool VTU_Interface<1U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<1U,InterFace>&,double);
template bool VTU_Interface<2U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<2U,InterFace>&,double);
template bool VTU_Interface<3U>::OutputDataToVTU(const string&,const set<string>&,const ModelSubDomain<3U,InterFace>&,double);


/// write multiblock file
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>
::OutputMultiBlockVTU( const string& fileName,
                       const vector<string>& fileNames,
                       const ModelSubDomain<dim,CELL>& subDomain )
{
    /// initialize file
    XML_Document outputFile;
    map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap
               = GetConnectivityMapMultiBlock( subDomain );
    outputFile = *ConnectivityFile<CELL>( connectivityMap, subDomain );

    outputFile.OpenNode( "VTKFile type=\"vtkMultiBlockDataSet\" version=\"0.9\" byte_order=\"LittleEndian\"" );
    /// Open multiblock section
    outputFile.OpenNode( "vtkMultiBlockDataSet");
    outputFile.OpenNode( "Block index=\"0\" name=\"Blocks\"" );

    size_t files( fileNames.size() );
    string stringNumber;
    string stringCache;
    string vtuFileName;
    for( size_t index=0; index<files; ++index )
    {
        stringNumber = number_to_string( index );
        vtuFileName  = fileNames[index];
        vtuFileName += ".vtu";
        // write DataSet section
        stringCache = "DataSet index=\""; stringCache += stringNumber;
        stringCache += "\" file=\""; stringCache += vtuFileName; stringCache += "\"";
        outputFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();
        outputFile.CloseNode( "DataSet");
    }

    /// Close multiblock section
    outputFile.CloseNode( "Block");
    outputFile.CloseNode( "vtkMultiBlockDataSet");
    outputFile.CloseNode( "VTKFile" );

    /// Close file
    CloseFile( fileName, ".vtm", outputFile );
}

template void VTU_Interface<1U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<1U,Element>&);
template void VTU_Interface<2U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<2U,Element>&);
template void VTU_Interface<3U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<3U,Element>&);

template void VTU_Interface<1U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<1U,Face>&);
template void VTU_Interface<2U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<2U,Face>&);
template void VTU_Interface<3U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<3U,Face>&);

template void VTU_Interface<1U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<1U,InterFace>&);
template void VTU_Interface<2U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<2U,InterFace>&);
template void VTU_Interface<3U>::OutputMultiBlockVTU( const string&,const vector<string>&,const ModelSubDomain<3U,InterFace>&);


/// we dispatch the output to node and element scalar/vector/tensor/array/flagged array
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>
::OutputFieldNodesAndElementDataToVTU( const string& fileName,
                   const list<Index>& fieldDataIndices,
                   const list<Index>& nodeIndices,
                   const list<Index>& elementIndices,
                   const ModelSubDomain<dim,CELL>& subDomain )
{
  /// initialize file
  XML_Document outputFile;
  map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap
             = GetConnectivityMap( subDomain );
  outputFile = *ConnectivityFile<CELL>( connectivityMap, subDomain );

  /// Open VTKFile section
  outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\"" );
  /// Open UnstructuredGrid section
  outputFile.OpenNode( "UnstructuredGrid" );

  /// 1. Regions data
  if( !fieldDataIndices.empty() )
      OutputFieldDataToVTU<CELL>(outputFile, subDomain, fieldDataIndices );

  /// 2. Nodes and Element's data

  if( !nodeIndices.empty() || !elementIndices.empty() )
  {
      /// 2.1 create connectivity, open Piece section
      EstablishConnectivityFile<CELL>( outputFile, subDomain );

      /// 2.2 point data arrays
      if( !nodeIndices.empty() )
          OutputPointDataToVTU(outputFile,subDomain,nodeIndices);

      /// 2.3 cell data arrays
      if( !elementIndices.empty() )
          OutputCellDataToVTU(outputFile,subDomain,elementIndices);

      /// 2.4 close Piece section
      outputFile.CloseNode( "Piece" );
  }

  /// Close UnstructuredGrid section
  outputFile.CloseNode( "UnstructuredGrid" );
  outputFile.CloseNode( "VTKFile" );

  /// Close file
  CloseFile( fileName, ".vtu", outputFile );

  return true;
}

template bool VTU_Interface<1U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::OutputFieldNodesAndElementDataToVTU(const string&,const list<Index>&,const list<Index>&,const list<Index>&,const ModelSubDomain<3U,InterFace>&);

template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>
::OutputElementBarycentricDataToVTU( const string& fileName,
                   const list<Index>& elementMatrixIndices,
                   const ModelSubDomain<dim,CELL>& subDomain )
{

  /// initialize file
  XML_Document outputFile;
  map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap
             = GetConnectivityMapBCS( subDomain );
  outputFile = *ConnectivityFile<CELL>( connectivityMap, subDomain );

  /// Open VTKFile section
  outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\"" );

  /// Open UnstructuredGrid section
  outputFile.OpenNode( "UnstructuredGrid" );

  /// 1. Cell centers data
  if( !elementMatrixIndices.empty() )
  {
      /// 1.1 create connectivity, open Piece section
      EstablishConnectivityFileBCPC<CELL>( outputFile, subDomain );

      /// 1.2 point data arrays
      OutputPointDataToVTU(outputFile,subDomain,elementMatrixIndices);

      /// 1.3 close Piece section
      outputFile.CloseNode( "Piece" );
  }

  /// Close UnstructuredGrid section
  outputFile.CloseNode( "UnstructuredGrid" );
  outputFile.CloseNode( "VTKFile" );

  /// Close file
  CloseFile( fileName, ".vtu", outputFile );

  return true;
}

template bool VTU_Interface<1U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::OutputElementBarycentricDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,InterFace>&);


template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>
::OutputRegionDataToVTU( const string& fileName,
                         const list<Index>& regionIndices,
                         const ModelSubDomain<dim,CELL>& subDomain )
{
  /// initialize file
  XML_Document outputFile;
  map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap
             = GetConnectivityMapRCS( subDomain );
  outputFile = *ConnectivityFile<CELL>( connectivityMap, subDomain );

  /// Open VTKFile section
  outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\"" );

  /// Open UnstructuredGrid section
  outputFile.OpenNode( "UnstructuredGrid" );

  /// 1. Region's data
  if( !regionIndices.empty() )
    {
      /// 1.1 create connectivity, open Piece section
      EstablishConnectivityFileRPC<CELL>( outputFile, subDomain );

      /// 1.2 point data arrays
      if( !regionIndices.empty() )
          OutputPointDataToVTU(outputFile,subDomain,regionIndices);

      /// 1.1. close Piece section
      outputFile.CloseNode( "Piece" );
    }

  /// Close UnstructuredGrid section
  outputFile.CloseNode( "UnstructuredGrid" );
  outputFile.CloseNode( "VTKFile" );

  /// Close file
  CloseFile( fileName, ".vtu", outputFile );

  return true;
}

template bool VTU_Interface<1U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::OutputRegionDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,InterFace>&);


template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>
::OutputFiniteElementIntegrationPointsDataToVTU( const string& fileName,
                   const list<Index>& feipIndices,
                   const ModelSubDomain<dim,CELL>& subDomain )
{
  /// initialize file
  XML_Document outputFile;
  map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap
             = GetConnectivityMapFEIPS( subDomain );
  outputFile = *ConnectivityFile<CELL>( connectivityMap, subDomain );

  /// Open VTKFile section
  outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\"" );

  /// Open UnstructuredGrid section
  outputFile.OpenNode( "UnstructuredGrid" );

  /// 1. Finite Element Integration Points
  if( !feipIndices.empty() )
  {
      /// 1.1 create connectivity, open Piece section
      EstablishConnectivityFileFEIP<CELL>( outputFile, subDomain );

      /// 1.2 point data arrays
      OutputPointDataToVTU(outputFile,subDomain,feipIndices);

      /// 1.3 close Piece section
      outputFile.CloseNode( "Piece" );
  }

  /// Close UnstructuredGrid section
  outputFile.CloseNode( "UnstructuredGrid" );
  outputFile.CloseNode( "VTKFile" );

  /// Close file
  CloseFile( fileName, ".vtu", outputFile );

  return true;
}

template bool VTU_Interface<1U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::OutputFiniteElementIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,InterFace>&);


template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>
::OutputFiniteVolumeSectorIntegrationPointsDataToVTU( const string& fileName,
                   const list<Index>& fvsipIndices,
                   const ModelSubDomain<dim,CELL>& subDomain )
{
  /// initialize file
  XML_Document outputFile;
  map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap
             = GetConnectivityMapFVSIPS( subDomain );
  outputFile = *ConnectivityFile<CELL>( connectivityMap, subDomain );

  /// Open VTKFile section
  outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\"" );

  /// Open UnstructuredGrid section
  outputFile.OpenNode( "UnstructuredGrid" );

  /// 1. Finite Volume Sector Integration Points
  if( !fvsipIndices.empty() )
  {
      /// 1.1 create connectivity, open Piece section
      EstablishConnectivityFileFVSIP<CELL>( outputFile, subDomain );

      /// 1.2 point data arrays
      OutputPointDataToVTU(outputFile,subDomain,fvsipIndices);

      /// 1.3 close Piece section
      outputFile.CloseNode( "Piece" );
  }

  /// Close UnstructuredGrid section
  outputFile.CloseNode( "UnstructuredGrid" );
  outputFile.CloseNode( "VTKFile" );

  /// Close file
  CloseFile( fileName, ".vtu", outputFile );

  return true;
}

template bool VTU_Interface<1U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,InterFace>&);


template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>
::OutputFiniteVolumeFacetIntegrationPointsDataToVTU( const string& fileName,
                   const list<Index>& fvfipIndices,
                   const ModelSubDomain<dim,CELL>& subDomain )
{
  /// initialize file
  XML_Document outputFile;
  map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap
             = GetConnectivityMapFVFIPS( subDomain );
  outputFile = *ConnectivityFile<CELL>( connectivityMap, subDomain );

  /// Open VTKFile section
  outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\"" );

  /// Open UnstructuredGrid section
  outputFile.OpenNode( "UnstructuredGrid" );

  /// 1. Finite Volume Facet Integration Points
  if( !fvfipIndices.empty() )
  {
      /// 1.1 create connectivity, open Piece section
      EstablishConnectivityFileFVFIP<CELL>( outputFile, subDomain );

      /// 1.2 point data arrays
      OutputPointDataToVTU(outputFile,subDomain,fvfipIndices);

      /// 1.3 close Piece section
      outputFile.CloseNode( "Piece" );
  }

  /// Close UnstructuredGrid section
  outputFile.CloseNode( "UnstructuredGrid" );
  outputFile.CloseNode( "VTKFile" );

  /// Close file
  CloseFile( fileName, ".vtu", outputFile );

  return true;
}

template bool VTU_Interface<1U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(const string&,const list<Index>&,const ModelSubDomain<3U,InterFace>&);


template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::OutputFieldDataToVTU( XML_Document& outputFile,
                                               const ModelSubDomain<dim,CELL>& subDomain,
                                               const list<Index>& indices )
{
    outputFile.OpenNode( "FieldData");
    // looping over properties to output
    for( list<Index>::const_iterator
         it = indices.begin(); it != indices.end(); ++it )
    {
        if( it->type == SCALAR )
          WriteFieldDataArray<ScalarVariable>( *it, outputFile );
        else if( it->type == ARRAY )
          WriteFieldDataArray<ArrayVariable>( *it, outputFile );
        else if( it->type == FLAGGEDARRAY )
          WriteFieldDataArray<FlaggedArrayVariable>( *it, outputFile );
        else if( it->type == VECTOR )
          WriteFieldDataArray<VectorVariable<dim> >( *it, outputFile );
        else if( it->type == TENSOR )
          WriteFieldDataArray<TensorVariable<dim> >( *it, outputFile );
    }
    outputFile.CloseNode( "FieldData" );
}

template void VTU_Interface<1U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<1U,Element>&,const list<Index>&);
template void VTU_Interface<2U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<2U,Element>&,const list<Index>&);
template void VTU_Interface<3U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<3U,Element>&,const list<Index>&);

template void VTU_Interface<1U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<1U,Face>&,const list<Index>&);
template void VTU_Interface<2U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<2U,Face>&,const list<Index>&);
template void VTU_Interface<3U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<3U,Face>&,const list<Index>&);

template void VTU_Interface<1U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<1U,InterFace>&,const list<Index>&);
template void VTU_Interface<2U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<2U,InterFace>&,const list<Index>&);
template void VTU_Interface<3U>::OutputFieldDataToVTU(XML_Document&,const ModelSubDomain<3U,InterFace>&,const list<Index>&);

template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::OutputPointDataToVTU( XML_Document& outputFile,
                                               const ModelSubDomain<dim,CELL>& subDomain,
                                               const list<Index>& indices )
{
  // active property
  string activeProperty;
  string variableName;

  // setting the first property as active in vtu visualizer
  activeProperty = "PointData ";
  if( indices.front().type == SCALAR  || indices.front().type == ARRAY || indices.front().type == FLAGGEDARRAY)
      activeProperty += "Scalars=";
  else if( indices.front().type == VECTOR )
      activeProperty += "Vectors=";
  else if( indices.front().type == TENSOR )
      activeProperty += "Tensors=";
  activeProperty += "\"";
  variableName = model_.Database().Name( indices.front() );
  variableName = FindVariableOutputAlias( variableName );
  // if array or flagged array, set first component as active
  if( indices.front().type == ARRAY || indices.front().type == FLAGGEDARRAY)
      variableName += "[ 0 ]";
  activeProperty += variableName; activeProperty += "\"";
  outputFile.OpenNode( activeProperty.c_str() ); activeProperty.clear();
  // looping over properties to output
  for( list<Index>::const_iterator
       it = indices.begin(); it != indices.end(); ++it )
  {
      if( it->type == SCALAR )
          WritePointDataArrayScalar<CELL>( *it, outputFile, subDomain );
      else if( it->type == ARRAY )
          WritePointDataArrayScalarArray<CELL>( *it, outputFile, subDomain );
      else if( it->type == FLAGGEDARRAY )
          WritePointDataArrayScalarFlaggedArray<CELL>( *it, outputFile, subDomain );
      else if( it->type == VECTOR )
          WritePointDataArrayVector<CELL>( *it, outputFile, subDomain );
      else if( it->type == TENSOR )
          WritePointDataArrayTensor<CELL>( *it, outputFile, subDomain );
  }
  outputFile.CloseNode( "PointData" );
  return;
}

template void VTU_Interface<1U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<1U,Element>&,const list<Index>&);
template void VTU_Interface<2U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<2U,Element>&,const list<Index>&);
template void VTU_Interface<3U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<3U,Element>&,const list<Index>&);

template void VTU_Interface<1U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<1U,Face>&,const list<Index>&);
template void VTU_Interface<2U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<2U,Face>&,const list<Index>&);
template void VTU_Interface<3U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<3U,Face>&,const list<Index>&);

template void VTU_Interface<1U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<1U,InterFace>&,const list<Index>&);
template void VTU_Interface<2U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<2U,InterFace>&,const list<Index>&);
template void VTU_Interface<3U>::OutputPointDataToVTU(XML_Document&,const ModelSubDomain<3U,InterFace>&,const list<Index>&);


template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::OutputCellDataToVTU( XML_Document& outputFile,
                                              const ModelSubDomain<dim,CELL>& subDomain,
                                              const list<Index>& indices )
{
    // active property
    string activeProperty;
    string variableName;

    // setting the first property as active in vtu visualizer
    activeProperty = "CellData ";
    if( indices.front().type == SCALAR || indices.front().type == ARRAY || indices.front().type == FLAGGEDARRAY )
        activeProperty += "Scalars=";
    else if( indices.front().type == VECTOR )
        activeProperty += "Vectors=";
    else if( indices.front().type == TENSOR )
        activeProperty += "Tensors=";
    activeProperty += "\"";
    variableName = model_.Database().Name( indices.front() );
    variableName = FindVariableOutputAlias( variableName );
    // if array or flagged array, set first component as active
    if( indices.front().type == ARRAY || indices.front().type == FLAGGEDARRAY)
        variableName += "[ 0 ]";
    activeProperty += variableName;  activeProperty += "\"";
    outputFile.OpenNode( activeProperty.c_str() ); activeProperty.clear();
    // looping over properties to output
    for( list<Index>::const_iterator it = indices.begin(); it != indices.end(); ++it )
      {
      if( it->type == SCALAR )
        WriteElementDataArrayScalar<CELL>( *it, outputFile, subDomain );
      else if( it->type == ARRAY )
        WriteElementDataArrayScalarArray<CELL>( *it, outputFile, subDomain );
      else if( it->type == FLAGGEDARRAY )
        WriteElementDataArrayScalarFlaggedArray<CELL>( *it, outputFile, subDomain );
      else if( it->type == VECTOR )
        WriteElementDataArrayVector<CELL>( *it, outputFile, subDomain );
      else if( it->type == TENSOR )
        WriteElementDataArrayTensor<CELL>( *it, outputFile, subDomain );
    }
    outputFile.CloseNode( "CellData" );
    return;
}

template void VTU_Interface<1U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<1U,Element>&,const list<Index>&);
template void VTU_Interface<2U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<2U,Element>&,const list<Index>&);
template void VTU_Interface<3U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<3U,Element>&,const list<Index>&);

template void VTU_Interface<1U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<1U,Face>&,const list<Index>&);
template void VTU_Interface<2U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<2U,Face>&,const list<Index>&);
template void VTU_Interface<3U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<3U,Face>&,const list<Index>&);

template void VTU_Interface<1U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<1U,InterFace>&,const list<Index>&);
template void VTU_Interface<2U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<2U,InterFace>&,const list<Index>&);
template void VTU_Interface<3U>::OutputCellDataToVTU(XML_Document&,const ModelSubDomain<3U,InterFace>&,const list<Index>&);


// SINGLE VERTEX OUTPUT
// --------------------

/// user interface to output vectors
template<uint32_t dim>
bool VTU_Interface<dim>::OutputVectorsToVTU( const string& fileName,
                                             const string& propertyCaption,
                                             const vector<vector<double> >& vectors )
{
  for( auto i = 0; i < vectors.size(); ++i )
    if( vectors.at( i ).size() != 3U )
      throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", "Vector required to have 3 components." );

  string fullFileName( fileName );
  fullFileName.append( ".vtu" );
  ofstream file( fullFileName.data(), ios::out );

  // header
  file << "<?xml version=\"0.9\"?>" << endl << endl;
  file << "<!--\n" << fileName << "\n-->" << endl << endl;

  // body - connectivity
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\">" << endl;
  file << "\t<UnstructuredGrid>\n\t\t<Piece NumberOfPoints=\"1\" NumberOfCells=\"1\">" << endl;
  file << "\t\t\t<Points>" << endl << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t0.0 0.0 0.0" << endl << "\t\t\t\t</DataArray>" << endl << "\t\t\t</Points>" << endl;
  file << "\t\t\t<Cells>" << endl << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t0" << endl << "\t\t\t\t</DataArray>" << endl;
  file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t1" << endl << "\t\t\t\t</DataArray>" << endl;
  file << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t" << VTK_VERTEX << endl << "\t\t\t\t</DataArray>" << endl << "\t\t\t</Cells>" << endl;

  // SKM FIX
  file.setf(ios::scientific);

  // body - values
  file << "\t\t\t<PointData Vectors=\"" << propertyCaption << " 0\">" << endl;
  for( auto i = 0; i < vectors.size(); ++i )
  {
    file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << propertyCaption << " " << i << "\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
    file << "\t\t\t\t\t";
    for( size_t ii = 0; ii < 3; ++ii )
      file << vectors.at( i ).at( ii ) << " ";
    file << endl;
    file << "\t\t\t\t</DataArray>" << endl;
  } // vectors

  file << "\t\t\t</PointData>" << endl;
  // body - close
  file << "\t\t</Piece>" << endl << "\t</UnstructuredGrid>" << endl << "</VTKFile>" << endl;
  file.close();
  return true;
}


template<uint32_t dim>
bool VTU_Interface<dim>::OutputPrincipalVectorsToVTU( const string& fileName, const string& propertyCaption,
                                                      const vector<double>& xyzLengths )
{
  if( xyzLengths.size() != 3U )
    throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", "Vector required to have 3 components." );
  vector<double> xVector( 3, 0. ); xVector.at( 0 ) = xyzLengths.at( 0 );
  vector<double> yVector( 3, 0. ); yVector.at( 1 ) = xyzLengths.at( 1 );
  vector<double> zVector( 3, 0. ); zVector.at( 2 ) = xyzLengths.at( 2 );
  vector<vector<double> > xyzVectors;
  xyzVectors.push_back( xVector );
  xyzVectors.push_back( yVector );
  xyzVectors.push_back( zVector );
  return OutputVectorsToVTU( fileName, propertyCaption, xyzVectors );
}


/// user interface to output a tensor @todo (3) Check inner vector sizes
template<uint32_t dim>
bool VTU_Interface<dim>::OutputTensorToVTU( const string& fileName, const string& propertyCaption,
                                            const vector<vector<double> >& tensor )
{
  if( tensor.size() != 3U )
    throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", "Tensor required to have 3 components." );

  string fullFileName( fileName );
  fullFileName.append( ".vtu" );
  ofstream file( fullFileName.data(), ios::out );

  // header
  file << "<?xml version=\"0.9\"?>" << endl << endl;
  file << "<!--\n" << fileName << "\n-->" << endl << endl;

  // body - connectivity
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"0.9\" byte_order=\"LittleEndian\">" << endl;
  file << "\t<UnstructuredGrid>\n\t\t<Piece NumberOfPoints=\"1\" NumberOfCells=\"1\">" << endl;
  file << "\t\t\t<Points>" << endl << "\t\t\t\t<DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t0.0 0.0 0.0" << endl << "\t\t\t\t</DataArray>" << endl << "\t\t\t</Points>" << endl;
  file << "\t\t\t<Cells>" << endl << "\t\t\t\t<DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t0" << endl << "\t\t\t\t</DataArray>" << endl;
  file << "\t\t\t\t<DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t1" << endl << "\t\t\t\t</DataArray>" << endl;
  file << "\t\t\t\t<DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t" << VTK_VERTEX << endl << "\t\t\t\t</DataArray>" << endl << "\t\t\t</Cells>" << endl;

  // SKM FIX
  file.setf(ios::scientific);

  // body - values
  file << "\t\t\t<PointData Tensors=\"" << propertyCaption << "\">" << endl;
  file << "\t\t\t\t<DataArray type=\"Float64\" Name=\"" << propertyCaption << "\" NumberOfComponents=\"9\" format=\"ascii\">" << endl;
  file << "\t\t\t\t\t";
  for( auto i = 0; i < 3; ++i )
    for( size_t ii = 0; ii < 3; ++ii )
      file << tensor.at( i ).at( ii ) << " ";

  file << endl;
  file << "\t\t\t\t</DataArray>" << endl;
  file << "\t\t\t</PointData>" << endl;
  // body - close
  file << "\t\t</Piece>" << endl << "\t</UnstructuredGrid>" << endl << "</VTKFile>" << endl;
  file.close();
  return true;
}


// WRITE VARIABLES
// --------------------

template<uint32_t dim>
void VTU_Interface<dim>::WriteScalar( XML_Document& vtu, const size_t& MAX_ENTRIES_PER_LINE,
                                      double scalarVariable, size_t& entriesOfLine, bool& newLine ) const
{
    string stringNumber;
    stringNumber = number_to_string( scalarVariable );
    vtu.InsertData( stringNumber.c_str() );
    // break or tab
    if( ++entriesOfLine == MAX_ENTRIES_PER_LINE )
    {
        entriesOfLine = 1;
        vtu.LineBreak();
        vtu.BringToLevel();
        newLine = true;
    }
    else
    {
        vtu.Tab();
        newLine = false;
    }
}

template<uint32_t dim>
void VTU_Interface<dim>::WriteVector( XML_Document& vtu, const size_t& MAX_ENTRIES_PER_LINE,
                                      const VectorVariable<dim>& vectorVariable, size_t& entriesOfLine, bool& newLine ) const
{
    string stringNumber;
    // inserting vector data
    for( auto i = 0; i < dim; ++i )
    {
      stringNumber = number_to_string( vectorVariable.Component( i ) );
      vtu.InsertData( stringNumber.c_str() );
      vtu.Tab();
    }
    if( dim != 3 )
    {
        for( auto i = dim; i < 3U; ++i )
        {
          stringNumber = number_to_string( 0. );
          vtu.InsertData( stringNumber.c_str() );
          vtu.Tab();
        }
    }
    // break or tab
    if( ++entriesOfLine == MAX_ENTRIES_PER_LINE )
    {
      entriesOfLine = 1;
      vtu.LineBreak();
      vtu.BringToLevel();
      newLine = true;
    }
    else
    {
      vtu.Tab();
      newLine = false;
    }
}

template<uint32_t dim>
void VTU_Interface<dim>::WriteTensor( XML_Document& vtu, const size_t& MAX_ENTRIES_PER_LINE,
                                      const TensorVariable<dim>& tensorVariable, size_t& entriesOfLine, bool& newLine ) const
{
    string stringNumber;
    double val;
    uint32_t rows( 3U );
    uint32_t cols( 3U );
    // inserting vector data
    for( auto column = 0; column < cols; ++column )
    {
      for( auto row = 0; row < rows; ++row )
      {
        val = ( ( row < dim ) && ( column < dim ) ? tensorVariable( row, column ) : 0.0 );
        stringNumber = number_to_string( val );
        vtu.InsertData( stringNumber.c_str() );
        vtu.Tab();
      }
    }
    // break or tab
    if( ++entriesOfLine == MAX_ENTRIES_PER_LINE )
    {
      entriesOfLine = 1;
      vtu.LineBreak();
      vtu.BringToLevel();
      newLine = true;
    }
    else
    {
      vtu.Tab();
      newLine = false;
    }
}

// WRITE FIELD DATA
// --------------------


/// writes scalar field data to xml document
template<uint32_t dim>
template<class Var>
void VTU_Interface<dim>::WriteFieldDataArray( const Index& key, XML_Document& vtu ) const
{
    size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 10 );
    string stringNumber;
    // variable info
    Var var;
    model_.Read( key, var );
    auto varSize{ key.dataDepth };
    string variableName = model_.Database().Name( key );
    variableName = FindVariableOutputAlias( variableName );
    stringNumber = number_to_string( varSize );

    string arrayTitle( "DataArray type=\"Float64\" Name=\"" );
    arrayTitle += variableName;
    arrayTitle += "\" NumberOfTuples=\"";
    arrayTitle += stringNumber;
    arrayTitle += "\" format=\"ascii\"";
    vtu.OpenNode( arrayTitle.c_str() );
    vtu.BringToLevel();
    bool newLine( false );

    for( auto i=0 ; i<varSize; ++i )
    {
        // inserting data
        stringNumber = number_to_string( var.Component( i ) );
        vtu.InsertData( stringNumber.c_str() );
        // break or tab
        if( ++entriesOfLine == MAX_ENTRIES_PER_LINE )
        {
            entriesOfLine = 1;
            vtu.LineBreak();
            vtu.BringToLevel();
            newLine = true;
        }
        else
        {
            vtu.Tab();
            newLine = false;
        }
    }

    if( !newLine )
        vtu.LineBreak();
    vtu.CloseNode( "DataArray" );
}

template void VTU_Interface<1U>::WriteFieldDataArray<ScalarVariable>(const Index&,XML_Document&) const;
template void VTU_Interface<2U>::WriteFieldDataArray<ScalarVariable>(const Index&,XML_Document&) const;
template void VTU_Interface<3U>::WriteFieldDataArray<ScalarVariable>(const Index&,XML_Document&) const;

template void VTU_Interface<1U>::WriteFieldDataArray<VectorVariable<1U> >(const Index&,XML_Document&) const;
template void VTU_Interface<2U>::WriteFieldDataArray<VectorVariable<2U> >(const Index&,XML_Document&) const;
template void VTU_Interface<3U>::WriteFieldDataArray<VectorVariable<3U> >(const Index&,XML_Document&) const;

template void VTU_Interface<1U>::WriteFieldDataArray<TensorVariable<1U> >(const Index&,XML_Document&) const;
template void VTU_Interface<2U>::WriteFieldDataArray<TensorVariable<2U> >(const Index&,XML_Document&) const;
template void VTU_Interface<3U>::WriteFieldDataArray<TensorVariable<3U> >(const Index&,XML_Document&) const;

template void VTU_Interface<1U>::WriteFieldDataArray<ArrayVariable>(const Index&,XML_Document&) const;
template void VTU_Interface<2U>::WriteFieldDataArray<ArrayVariable>(const Index&,XML_Document&) const;
template void VTU_Interface<3U>::WriteFieldDataArray<ArrayVariable>(const Index&,XML_Document&) const;

template void VTU_Interface<1U>::WriteFieldDataArray<FlaggedArrayVariable>(const Index&,XML_Document&) const;
template void VTU_Interface<2U>::WriteFieldDataArray<FlaggedArrayVariable>(const Index&,XML_Document&) const;
template void VTU_Interface<3U>::WriteFieldDataArray<FlaggedArrayVariable>(const Index&,XML_Document&) const;


// WRITE POINT DATA
// --------------------


/// writes scalar point data array to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WritePointDataArrayScalar( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 10 );
  string arrayTitle( "DataArray type=\"Float64\" Name=\"" );
  string variableName = model_.Database().Name( key );
  variableName = FindVariableOutputAlias( variableName );
  arrayTitle += variableName;
  arrayTitle += "\" NumberOfComponents=\"1\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  ScalarVariable scalarVariable;
  bool newLine( false );
  if( key.place == NODE )
  {
      // looping over regions nodes
      const auto domainNodesEnd( subDomain.NodesEnd() );
      for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
      {
          // inserting scalar data
          (*it)->Read( key, scalarVariable );
          WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
      }
  }
  else if( key.place == MODEL )
  {
      model_.Read( key, scalarVariable );
      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
  }
  else if( key.place == REGION || key.place == BOUNDARY || key.place == SPLIT_BOUNDARY )
  {
      if ( key.place == REGION ) dynamic_cast<const Region<dim>&>(subDomain).Read( key, scalarVariable );
      else if ( key.place == BOUNDARY ) dynamic_cast<const Boundary<dim>&>(subDomain).Read( key, scalarVariable );
      else if ( key.place == SPLIT_BOUNDARY ) dynamic_cast<const SplitBoundary<dim>&>(subDomain).Read( key, scalarVariable );

      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
  }
  else if( key.place == ELEMENT_INTEGRATION_POINT  || key.place == FACE_INTEGRATION_POINT  || key.place == INTER_FACE_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
      {
          auto ips = (*it)->IntegrationPoints();
          for( auto ip = 0; ip < ips; ++ip )
          {
              // inserting scalar data
              (*it)->Read( ip, key, scalarVariable );
              WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
          }
      }
  }
  else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
      {
          auto sectors = (*it)->Sectors();
          for( uint32_t sid = 0; sid < sectors; ++sid )
          {
              auto ips = (*it)->FV()->IntegrationPointsPerSector( sid );
              for( auto ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( sid, ip, key, scalarVariable );
                  WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
              }
          }
      }
  }
  else if(  key.place == FACET_INTEGRATION_POINT  || key.place == FACE_FACET_INTEGRATION_POINT  || key.place == INTER_FACE_FACET_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
      {
          const auto facets = (*it)->Facets();
          for( auto fid = 0; fid < facets; ++fid )
          {
              const auto ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
              for( auto  ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( fid, ip, key, scalarVariable );
                  WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
              }
          }
      }
  }

  if( !newLine )
      vtu.LineBreak();
  vtu.CloseNode( "DataArray" );
}

template void VTU_Interface<1U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// writes vector point data array to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WritePointDataArrayVector( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 4 );
  string arrayTitle( "DataArray type=\"Float64\" Name=\"" );
  string variableName = model_.Database().Name( key );
  variableName = FindVariableOutputAlias( variableName );
  arrayTitle += variableName; arrayTitle += "\" NumberOfComponents=\"3\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  VectorVariable<dim> vectorVariable;
  bool newLine( false );
  if( key.place == NODE )
  {
      // looping over regions nodes
      const auto domainNodesEnd( subDomain.NodesEnd() );
      for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
      {
        // acquiring vector data
        (*it)->Read( key, vectorVariable );
        WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
      }
  }
  else if( key.place == MODEL )
  {
      model_.Read( key, vectorVariable );
      WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
  }
  else if( key.place == REGION || key.place == BOUNDARY || key.place == SPLIT_BOUNDARY )
  {
      if ( key.place == REGION ) dynamic_cast<const Region<dim>&>(subDomain).Read( key, vectorVariable );
      else if ( key.place == BOUNDARY ) dynamic_cast<const Boundary<dim>&>(subDomain).Read( key, vectorVariable );
      else if ( key.place == SPLIT_BOUNDARY ) dynamic_cast<const SplitBoundary<dim>&>(subDomain).Read( key, vectorVariable );

      WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
  }
  else if( key.place == ELEMENT_INTEGRATION_POINT  || key.place == FACE_INTEGRATION_POINT  || key.place == INTER_FACE_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin();
           it != domainElementsEnd; ++it )
      {
          const auto  ips = (*it)->IntegrationPoints();
          for( auto  ip = 0; ip < ips; ++ip )
          {
              // inserting scalar data
              (*it)->Read( ip, key, vectorVariable );
              WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
          }
      }
  }
  else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin();
           it != domainElementsEnd; ++it )
      {
          const auto sectors = (*it)->Sectors();
          for( auto sid = 0; sid < sectors; ++sid )
          {
              const auto ips = (*it)->FV()->IntegrationPointsPerSector( sid );
              for( auto ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( sid, ip, key, vectorVariable );
                  WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
              }
          }
      }
  }
  else if(  key.place == FACET_INTEGRATION_POINT  || key.place == FACE_FACET_INTEGRATION_POINT  || key.place == INTER_FACE_FACET_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin();
           it != domainElementsEnd; ++it )
      {
          const auto facets = (*it)->Facets();
          for( auto fid = 0; fid < facets; ++fid )
          {
              const auto ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
              for( auto ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( fid, ip, key, vectorVariable );
                  WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
              }
          }
      }
  }

  if( !newLine )
    vtu.LineBreak();
  vtu.CloseNode( "DataArray" );
}

template void VTU_Interface<1U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WritePointDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// writes tensor point data array to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WritePointDataArrayTensor( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 3 );
  string stringNumber, arrayTitle( "DataArray type=\"Float64\" Name=\"" );
  string variableName = model_.Database().Name( key );
  variableName = FindVariableOutputAlias( variableName );
  arrayTitle += variableName; arrayTitle += "\" NumberOfComponents=\"";
  const size_t tensor_components( 9U ); // only 3D representation
  stringNumber = number_to_string( tensor_components );
  arrayTitle += stringNumber; arrayTitle += "\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  bool newLine( false );
  TensorVariable<dim> tensorVariable;
  if( key.place == NODE )
  {
      // looping over regions nodes
      const auto domainNodesEnd( subDomain.NodesEnd() );
      for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
      {
        // acquiring vector data
        (*it)->Read( key, tensorVariable );
        WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
    }
  }
  else if( key.place == MODEL )
  {
      model_.Read( key, tensorVariable );
      WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
  }
  else if( key.place == REGION || key.place == BOUNDARY || key.place == SPLIT_BOUNDARY )
  {
      if ( key.place == REGION ) dynamic_cast<const Region<dim>&>(subDomain).Read( key, tensorVariable );
      else if ( key.place == BOUNDARY ) dynamic_cast<const Boundary<dim>&>(subDomain).Read( key, tensorVariable );
      else if ( key.place == SPLIT_BOUNDARY ) dynamic_cast<const SplitBoundary<dim>&>(subDomain).Read( key, tensorVariable );

      WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
  }
  else if( key.place == ELEMENT_INTEGRATION_POINT  || key.place == FACE_INTEGRATION_POINT  || key.place == INTER_FACE_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin();
           it != domainElementsEnd; ++it )
      {
          const auto ips = (*it)->IntegrationPoints();
          for( auto ip = 0; ip < ips; ++ip )
          {
              // inserting scalar data
              (*it)->Read( ip, key, tensorVariable );
              WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
          }
      }
  }
  else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin();
           it != domainElementsEnd; ++it )
      {
          const auto sectors = (*it)->Sectors();
          for( auto sid = 0; sid < sectors; ++sid )
          {
              const auto ips = (*it)->FV()->IntegrationPointsPerSector( sid );
              for( auto ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( sid, ip, key, tensorVariable );
                  WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
              }
          }
      }
  }
  else if(  key.place == FACET_INTEGRATION_POINT  || key.place == FACE_FACET_INTEGRATION_POINT  || key.place == INTER_FACE_FACET_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin();
           it != domainElementsEnd; ++it )
      {
          const auto facets = (*it)->Facets();
          for( auto fid = 0; fid < facets; ++fid )
          {
              const auto ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
              for( auto ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( fid, ip, key, tensorVariable );
                  WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
              }
          }
      }
  }

  if( !newLine )
    vtu.LineBreak();
  vtu.CloseNode( "DataArray" );
}

template void VTU_Interface<1U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WritePointDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// writes point data array ( array of scalars ) to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WritePointDataArrayScalarArray( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
    size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 10 );
    string stringNumber;
    string arrayTitle;
    string variableName;
    string prefix;
    ArrayVariable arrayVariable;
    for( size_t component=0; component<key.dataDepth; component++)
    {
        arrayTitle = "DataArray type=\"Float64\" Name=\"";
        variableName = model_.Database().Name( key );
        variableName = FindVariableOutputAlias( variableName );
        stringNumber = number_to_string( component );
        prefix = "";
      if (this->array_index_to_name_.find(variableName) != this->array_index_to_name_.end()){
          vector<string> component_names (this->array_index_to_name_.at(variableName));
          if (component < component_names.size())
              prefix = component_names[component];
      }
      variableName = prefix + variableName + "[ " + stringNumber + " ]";
      arrayTitle += variableName;
      arrayTitle += "\" NumberOfComponents=\"1\" format=\"ascii\"";
      vtu.OpenNode( arrayTitle.c_str() );
      vtu.BringToLevel();
      bool newLine( false );
      if( key.place == NODE )
      {
          // looping over regions nodes
          const auto domainNodesEnd( subDomain.NodesEnd() );
          for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
          {
              // inserting scalar array data
              (*it)->Read( key, arrayVariable );
              WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
          }
      }
      else if( key.place == MODEL )
      {
          model_.Read( key, arrayVariable );
          WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
      }
      else if( key.place == REGION || key.place == BOUNDARY || key.place == SPLIT_BOUNDARY )
      {
          if ( key.place == REGION ) dynamic_cast<const Region<dim>&>(subDomain).Read( key, arrayVariable );
          else if ( key.place == BOUNDARY ) dynamic_cast<const Boundary<dim>&>(subDomain).Read( key, arrayVariable );
          else if ( key.place == SPLIT_BOUNDARY ) dynamic_cast<const SplitBoundary<dim>&>(subDomain).Read( key, arrayVariable );

          WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
      }
      else if( key.place == ELEMENT_INTEGRATION_POINT  || key.place == FACE_INTEGRATION_POINT  || key.place == INTER_FACE_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.ElementsEnd() );
          for( auto it = subDomain.ElementsBegin();
               it != domainElementsEnd; ++it )
          {
              const auto ips = (*it)->IntegrationPoints();
              for( auto ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( ip, key, arrayVariable );
                  WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
              }
          }
      }
      else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.ElementsEnd() );
          for( auto it = subDomain.ElementsBegin();
               it != domainElementsEnd; ++it )
          {
              const auto sectors = (*it)->Sectors();
              for( auto sid = 0; sid < sectors; ++sid )
              {
                  const auto ips = (*it)->FV()->IntegrationPointsPerSector( sid );
                  for( auto ip = 0; ip < ips; ++ip )
                  {
                      // inserting scalar data
                      (*it)->Read( sid, ip, key, arrayVariable );
                      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
                  }
              }
          }
      }
      else if(  key.place == FACET_INTEGRATION_POINT  || key.place == FACE_FACET_INTEGRATION_POINT  || key.place == INTER_FACE_FACET_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.ElementsEnd() );
          for( auto it = subDomain.ElementsBegin();
               it != domainElementsEnd; ++it )
          {
              const auto facets = (*it)->Facets();
              for( auto fid = 0; fid < facets; ++fid )
              {
                  const auto ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
                  for( auto ip = 0; ip < ips; ++ip )
                  {
                      // inserting scalar data
                      (*it)->Read( fid, ip, key, arrayVariable );
                      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
                  }
              }
          }
      }

      if( !newLine )
        vtu.LineBreak();
      vtu.CloseNode( "DataArray" );
    }
}

template void VTU_Interface<1U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;

/// writes point data array ( flagged array of scalars ) to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WritePointDataArrayScalarFlaggedArray( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
    size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 10 );
    string stringNumber;
    string arrayTitle;
    string variableName;
    string prefix;
    FlaggedArrayVariable flaggedArrayVariable;
    for( size_t component=0; component<key.dataDepth; component++)
    {
        arrayTitle = "DataArray type=\"Float64\" Name=\"";
        variableName = model_.Database().Name( key );
        variableName = FindVariableOutputAlias( variableName );
        stringNumber = number_to_string( component );
        prefix = "";
      if (this->array_index_to_name_.find(variableName) != this->array_index_to_name_.end()){
          vector<string> component_names (this->array_index_to_name_.at(variableName));
          if (component < component_names.size())
              prefix = component_names[component];
      }
      variableName = prefix + variableName + "[ " + stringNumber + " ]";
      arrayTitle += variableName;
      arrayTitle += "\" NumberOfComponents=\"1\" format=\"ascii\"";
      vtu.OpenNode( arrayTitle.c_str() );
      vtu.BringToLevel();
      bool newLine( false );
      if( key.place == NODE )
      {
          // looping over regions nodes
          const auto domainNodesEnd( subDomain.NodesEnd() );
          for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
          {
            // inserting scalar flagged array data
            (*it)->Read( key, flaggedArrayVariable );
            WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
          }
      }
      else if( key.place == MODEL )
      {
          model_.Read( key, flaggedArrayVariable );
          WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
      }
      else if( key.place == REGION || key.place == BOUNDARY || key.place == SPLIT_BOUNDARY )
      {
          if ( key.place == REGION ) dynamic_cast<const Region<dim>&>(subDomain).Read( key, flaggedArrayVariable );
          else if ( key.place == BOUNDARY ) dynamic_cast<const Boundary<dim>&>(subDomain).Read( key, flaggedArrayVariable );
          else if ( key.place == SPLIT_BOUNDARY ) dynamic_cast<const SplitBoundary<dim>&>(subDomain).Read( key, flaggedArrayVariable );

          WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
      }
      else if( key.place == ELEMENT_INTEGRATION_POINT  || key.place == FACE_INTEGRATION_POINT  || key.place == INTER_FACE_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.ElementsEnd() );
          for( auto it = subDomain.ElementsBegin();
               it != domainElementsEnd; ++it )
          {
              const auto ips = (*it)->IntegrationPoints();
              for( auto ip = 0; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( ip, key, flaggedArrayVariable );
                  WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
              }
          }
      }
      else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.ElementsEnd() );
          for( auto it = subDomain.ElementsBegin();
               it != domainElementsEnd; ++it )
          {
              const auto sectors = (*it)->Sectors();
              for( auto sid = 0; sid < sectors; ++sid )
              {
                  const auto ips = (*it)->FV()->IntegrationPointsPerSector( sid );
                  for( auto ip = 0; ip < ips; ++ip )
                  {
                      // inserting scalar data
                      (*it)->Read( sid, ip, key, flaggedArrayVariable );
                      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
                  }
              }
          }
      }
      else if(  key.place == FACET_INTEGRATION_POINT  || key.place == FACE_FACET_INTEGRATION_POINT  || key.place == INTER_FACE_FACET_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.ElementsEnd() );
          for( auto it = subDomain.ElementsBegin();
               it != domainElementsEnd; ++it )
          {
              const auto facets = (*it)->Facets();
              for( auto fid = 0; fid < facets; ++fid )
              {
                  const auto ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
                  for( auto ip = 0; ip < ips; ++ip )
                  {
                      // inserting scalar data
                      (*it)->Read( fid, ip, key, flaggedArrayVariable );
                      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
                  }
              }
          }
      }

      if( !newLine )
        vtu.LineBreak();
      vtu.CloseNode( "DataArray" );
    }
}

template void VTU_Interface<1U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WritePointDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


// WRITE ELEMENT DATA
// --------------------


/// writes cell data array to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WriteElementDataArrayScalar( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 10 );
  string arrayTitle( "DataArray type=\"Float64\" Name=\"" );
  string variableName = model_.Database().Name( key );
  variableName = FindVariableOutputAlias( variableName );
  arrayTitle += variableName; arrayTitle += "\" NumberOfComponents=\"1\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  ScalarVariable scalarVariable;
  bool newLine( false );
  // looping over domains elements
  const auto domainElementsEnd( subDomain.ElementsEnd() );
  for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
    {
      // inserting scalar data
      (*it)->Read( key, scalarVariable );
      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
    }
  if( !newLine )
    vtu.LineBreak();
  vtu.CloseNode( "DataArray" );
}

template void VTU_Interface<1U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalar(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;

/// writes vector cell data array to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WriteElementDataArrayVector( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 4 );
  string arrayTitle( "DataArray type=\"Float64\" Name=\"" );
  string variableName = model_.Database().Name( key );
  variableName = FindVariableOutputAlias( variableName );
  arrayTitle += variableName; arrayTitle += "\" NumberOfComponents=\"3\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  bool newLine( false );
  VectorVariable<dim> vectorVariable;
  // looping over regions elements
  const auto domainElementsEnd( subDomain.ElementsEnd() );
  for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
    {
      // acquiring vector data
      (*it)->Read( key, vectorVariable );
      WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
    }
  if( !newLine )
    vtu.LineBreak();
  vtu.CloseNode( "DataArray" );
}

template void VTU_Interface<1U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayVector(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// writes tensor cell data array to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WriteElementDataArrayTensor( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 3 );
  string stringNumber, arrayTitle( "DataArray type=\"Float64\" Name=\"" );
  string variableName = model_.Database().Name( key );
  variableName = FindVariableOutputAlias( variableName );
  arrayTitle += variableName; arrayTitle += "\" NumberOfComponents=\"";
  const size_t tensor_components( 9U ); // only 3D representation
  stringNumber = number_to_string( tensor_components );
  arrayTitle += stringNumber; arrayTitle += "\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  bool newLine( false );
  TensorVariable<dim> tensorVariable;
  // looping over regions nodes
  const auto domainElementsEnd( subDomain.ElementsEnd() );
  for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
    {
      // acquiring vector data
      (*it)->Read( key, tensorVariable );
      WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
    }
  if( !newLine )
    vtu.LineBreak();
  vtu.CloseNode( "DataArray" );
}

template void VTU_Interface<1U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayTensor(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// writes cell data array ( array of scalars ) to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WriteElementDataArrayScalarArray( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 10 );
  string stringNumber;
  string arrayTitle;
  string variableName;
  string prefix;
  ArrayVariable arrayVariable;
  for( size_t component=0; component<key.dataDepth; component++)
  {
      arrayTitle = "DataArray type=\"Float64\" Name=\"";
      variableName = model_.Database().Name( key );
      variableName = FindVariableOutputAlias( variableName );
      stringNumber = number_to_string( component );
      prefix = "";
      if (this->array_index_to_name_.find(variableName) != this->array_index_to_name_.end()){
          vector<string> component_names (this->array_index_to_name_.at(variableName));
          if (component < component_names.size())
              prefix = component_names[component];
      }
      variableName = prefix + variableName + "[ " + stringNumber + " ]";
      arrayTitle += variableName; arrayTitle += "\" NumberOfComponents=\"1\" format=\"ascii\"";
      vtu.OpenNode( arrayTitle.c_str() );
      vtu.BringToLevel();
      bool newLine( false );
      // looping over domains elements
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
        {
            // inserting scalar array data
            (*it)->Read( key, arrayVariable );
            WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
        }
      if( !newLine )
          vtu.LineBreak();
      vtu.CloseNode( "DataArray" );
  }
}

template void VTU_Interface<1U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalarArray(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;

/// writes cell data array ( array of scalars ) to xml document
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::WriteElementDataArrayScalarFlaggedArray( const Index& key, XML_Document& vtu, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  size_t entriesOfLine( 2 ); const size_t MAX_ENTRIES_PER_LINE( 10 );
  string stringNumber;
  string arrayTitle;
  string variableName;
  string prefix;
  FlaggedArrayVariable flaggedArrayVariable;
  for( size_t component=0; component<key.dataDepth; component++)
  {
      arrayTitle = "DataArray type=\"Float64\" Name=\"";
      variableName = model_.Database().Name( key );
      variableName = FindVariableOutputAlias( variableName );
      stringNumber = number_to_string( component );
      prefix = "";
      if (this->array_index_to_name_.find(variableName) != this->array_index_to_name_.end()){
          vector<string> component_names (this->array_index_to_name_.at(variableName));
          if (component < component_names.size())
              prefix = component_names[component];
      }
      variableName = prefix + variableName + "[ " + stringNumber + " ]";
      arrayTitle += variableName; arrayTitle += "\" NumberOfComponents=\"1\" format=\"ascii\"";
      vtu.OpenNode( arrayTitle.c_str() );
      vtu.BringToLevel();
      bool newLine( false );
      // looping over domains elements
      const auto domainElementsEnd( subDomain.ElementsEnd() );
      for( auto it = subDomain.ElementsBegin(); it != domainElementsEnd; ++it )
      {
        // inserting scalar flagged array data
        (*it)->Read( key, flaggedArrayVariable );
        WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
      }
      if( !newLine )
        vtu.LineBreak();
      vtu.CloseNode( "DataArray" );
  }
}

template void VTU_Interface<1U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::WriteElementDataArrayScalarFlaggedArray(const Index&,XML_Document&,const ModelSubDomain<3U,InterFace>&) const;









// CONNECTIVITY
// ----------------

template<uint32_t dim>
template<template <uint32_t> class CELL>
XML_Document* VTU_Interface<dim>::findConnectivityFile( map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap, const ModelSubDomain<dim,CELL>& subDomain )
{
  typename map<const ModelSubDomain<dim,CELL>*,XML_Document*>::const_iterator connectivityEntry( connectivityMap.find( &subDomain ) );
  if( connectivityEntry != connectivityMap.end() )
    return connectivityEntry->second;
  return NULL;
}
template XML_Document* VTU_Interface<1U>::findConnectivityFile<Element>(map<const ModelSubDomain<1U,Element>*,XML_Document*>&,const ModelSubDomain<1U,Element>&);
template XML_Document* VTU_Interface<2U>::findConnectivityFile<Element>(map<const ModelSubDomain<2U,Element>*,XML_Document*>&,const ModelSubDomain<2U,Element>&);
template XML_Document* VTU_Interface<3U>::findConnectivityFile<Element>(map<const ModelSubDomain<3U,Element>*,XML_Document*>&,const ModelSubDomain<3U,Element>&);

template XML_Document* VTU_Interface<1U>::findConnectivityFile<Face>(map<const ModelSubDomain<1U,Face>*,XML_Document*>&,const ModelSubDomain<1U,Face>&);
template XML_Document* VTU_Interface<2U>::findConnectivityFile<Face>(map<const ModelSubDomain<2U,Face>*,XML_Document*>&,const ModelSubDomain<2U,Face>&);
template XML_Document* VTU_Interface<3U>::findConnectivityFile<Face>(map<const ModelSubDomain<3U,Face>*,XML_Document*>&,const ModelSubDomain<3U,Face>&);

template XML_Document* VTU_Interface<1U>::findConnectivityFile<InterFace>(map<const ModelSubDomain<1U,InterFace>*,XML_Document*>&,const ModelSubDomain<1U,InterFace>&);
template XML_Document* VTU_Interface<2U>::findConnectivityFile<InterFace>(map<const ModelSubDomain<2U,InterFace>*,XML_Document*>&,const ModelSubDomain<2U,InterFace>&);
template XML_Document* VTU_Interface<3U>::findConnectivityFile<InterFace>(map<const ModelSubDomain<3U,InterFace>*,XML_Document*>&,const ModelSubDomain<3U,InterFace>&);

template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::insertConnectivityFile( map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap, XML_Document* newConnectivityFile, const ModelSubDomain<dim,CELL>& subDomain )
{
  connectivityMap.insert( make_pair( &subDomain, newConnectivityFile ) );
  return true;
}

template bool VTU_Interface<1U>::insertConnectivityFile(map<const ModelSubDomain<1U,Element>*,XML_Document*>&,XML_Document*,const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::insertConnectivityFile(map<const ModelSubDomain<2U,Element>*,XML_Document*>&,XML_Document*,const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::insertConnectivityFile(map<const ModelSubDomain<3U,Element>*,XML_Document*>&,XML_Document*,const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::insertConnectivityFile(map<const ModelSubDomain<1U,Face>*,XML_Document*>&,XML_Document*,const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::insertConnectivityFile(map<const ModelSubDomain<2U,Face>*,XML_Document*>&,XML_Document*,const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::insertConnectivityFile(map<const ModelSubDomain<3U,Face>*,XML_Document*>&,XML_Document*,const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::insertConnectivityFile(map<const ModelSubDomain<1U,InterFace>*,XML_Document*>&,XML_Document*,const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::insertConnectivityFile(map<const ModelSubDomain<2U,InterFace>*,XML_Document*>&,XML_Document*,const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::insertConnectivityFile(map<const ModelSubDomain<3U,InterFace>*,XML_Document*>&,XML_Document*,const ModelSubDomain<3U,InterFace>&);


/// returns a pointer to the regions connectivity file. creates one if not existing yet
template<uint32_t dim>
template<template <uint32_t> class CELL>
XML_Document* VTU_Interface<dim>::ConnectivityFile( map<const ModelSubDomain<dim,CELL>*,XML_Document*>& connectivityMap, const ModelSubDomain<dim,CELL>& subDomain )
{
  // looks for corresponding entry for region parameter
  XML_Document* connectivityFile( findConnectivityFile( connectivityMap, subDomain ) );

  // if there's already an entry, return the connectivity file ptr
  if( connectivityFile )
    return connectivityFile;

  // if not, we go on and create one
  XML_Document* newConnectivityFile( new XML_Document );

  EstablishConnectivityFileHeader( *newConnectivityFile );

  // create an entry for region with connectivity file ptr
  insertConnectivityFile( connectivityMap, newConnectivityFile, subDomain );

  // and return ptr
  return newConnectivityFile;
}

template XML_Document* VTU_Interface<1U>::ConnectivityFile(map<const ModelSubDomain<1U,Element>*,XML_Document*>&,const ModelSubDomain<1U,Element>&);
template XML_Document* VTU_Interface<2U>::ConnectivityFile(map<const ModelSubDomain<2U,Element>*,XML_Document*>&,const ModelSubDomain<2U,Element>&);
template XML_Document* VTU_Interface<3U>::ConnectivityFile(map<const ModelSubDomain<3U,Element>*,XML_Document*>&,const ModelSubDomain<3U,Element>&);

template XML_Document* VTU_Interface<1U>::ConnectivityFile(map<const ModelSubDomain<1U,Face>*,XML_Document*>&,const ModelSubDomain<1U,Face>&);
template XML_Document* VTU_Interface<2U>::ConnectivityFile(map<const ModelSubDomain<2U,Face>*,XML_Document*>&,const ModelSubDomain<2U,Face>&);
template XML_Document* VTU_Interface<3U>::ConnectivityFile(map<const ModelSubDomain<3U,Face>*,XML_Document*>&,const ModelSubDomain<3U,Face>&);

template XML_Document* VTU_Interface<1U>::ConnectivityFile(map<const ModelSubDomain<1U,InterFace>*,XML_Document*>&,const ModelSubDomain<1U,InterFace>&);
template XML_Document* VTU_Interface<2U>::ConnectivityFile(map<const ModelSubDomain<2U,InterFace>*,XML_Document*>&,const ModelSubDomain<2U,InterFace>&);
template XML_Document* VTU_Interface<3U>::ConnectivityFile(map<const ModelSubDomain<3U,InterFace>*,XML_Document*>&,const ModelSubDomain<3U,InterFace>&);


/// @todo (2-F) This should return filename up to user call
template<uint32_t dim>
bool VTU_Interface<dim>::CloseFile( const string& fileName, const string& extension, XML_Document& outputFile )
{
  ErrorHandler& csmp_error ( ErrorHandler::Instance() );

  // write to file
  string outputName = fileName + extension;
  outputFile.WriteToFile( outputName.c_str() );

  // io feedback
  if ( csmp_error.Verbose() )
      cout << "\nVTU_Interface<dim>::OutputDataToVTU: 'Done writing " << outputName << "'\n";
  return true;
}

/// writes common header unstructured grids
template<uint32_t dim>
void VTU_Interface<dim>::EstablishConnectivityFileHeader( XML_Document& connectivityFile ) const
  {
    // writing problem header
    connectivityFile.AddInfo( "xml version=\"0.9\"" );
    connectivityFile.AddComment( problemTitle_.c_str() );
  }

/// establishes the connectivity file for given region
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::EstablishConnectivityFile( XML_Document& connectivityFile, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  ErrorHandler& csmp_error ( ErrorHandler::Instance() );

  // io feedback
  if ( csmp_error.Verbose() )
      cout << "\nVTU_Interface<dim>::EstablishConnectivityFile: 'Writing connectivity file for sub domain " << DomainName( subDomain ) << " ...";

  // reference to the domain of concern and renumbering its node indices
  subDomain.RenumberNodes();

  // getting node and element count, creating working strings
  string stringNumber, stringCache;

  // opening piece node
  const size_t DOMAIN_NODES( subDomain.Nodes() ), DOMAIN_ELEMENTS( subDomain.Elements() );
  stringNumber = number_to_string( DOMAIN_NODES );
  stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;
  stringNumber = number_to_string( DOMAIN_ELEMENTS );
  stringCache += "\" NumberOfCells=\""; stringCache += stringNumber; stringCache += "\"";
  connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

  // points node
  connectivityFile.OpenNode( "Points" );
  connectivityFile.OpenNode( "DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
  // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
  size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
  bool newLine( false );
  connectivityFile.BringToLevel();
  // looping over all the region's nodes
  const auto domainVerticesEnd( subDomain.NodesEnd() );
  for( auto it = subDomain.NodesBegin(); it != domainVerticesEnd; ++it, ++entriesOfLine )
  {
    // writing x,y and z coordinates(tab seperated)
    stringNumber = number_to_string( (*it)->x() );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = number_to_string( yCoordinate( (*it)->Coordinate() ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = number_to_string( zCoordinate( (*it)->Coordinate() ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    // if maximum of entris per line is reached, start new one, else tab
    if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE )
    {
      entriesOfLine = 1;
      connectivityFile.LineBreak();
      connectivityFile.BringToLevel();
      newLine = true;
    }
    else
    {
      connectivityFile.Tab();
      newLine = false;
    }
  } // looping regions nodes
  // closing points node
  if( !newLine )
    connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  connectivityFile.CloseNode( "Points" );

  // cells node
  connectivityFile.OpenNode( "Cells" );
  // cell connectivity
  connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  // establish a vector with vtk element types
  vector<VTK_TYPE> elementTypesVTK;
  VTK_TYPE elementTypeVTK;
  vector<long> vtkNodeNumbering;
  // looping region's elements
  const size_t MAX_CONNECTIVITY_ENTRIES_PER_LINE( 20 ); entriesOfLine = 2;
  const auto domainSimplicesEnd( subDomain.ElementsEnd() );
  for( auto it = subDomain.ElementsBegin(); it != domainSimplicesEnd; ++it )
  {
    // storing elements VTK type
    elementTypeVTK = ElementType( *it );
    elementTypesVTK.push_back( elementTypeVTK );
    // Quadratic Hexahedron and Quadratic Wedge differ in node numbering(CSMP vs VTK)
    if( elementTypeVTK == VTK_QUADRATIC_HEXAHEDRON )
    {
      vtkNodeNumbering.clear();
      QuadraticHexahedronConnectivity( *it, vtkNodeNumbering );
      assert( vtkNodeNumbering.size() == (*it)->Nodes() );
      for( auto i = 0; i < vtkNodeNumbering.size(); ++i, ++entriesOfLine )
      {
         // writing node id to vtu document
        stringNumber = number_to_string( vtkNodeNumbering[i] );
        connectivityFile.InsertData( stringNumber.c_str() );
        // line break if limeit entries reached, else tab
        if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
        {
          entriesOfLine = 1;
          connectivityFile.LineBreak();
          connectivityFile.BringToLevel();
          newLine = true;
        }
        else
        {
          connectivityFile.Tab();
          newLine = false;
        }
      }
      // proceed to next elementactiveProperty += "\"";
      continue;
    }
    if( elementTypeVTK == VTK_QUADRATIC_WEDGE )
    {
      vtkNodeNumbering.clear();
      QuadraticWedgeConnectivity( *it, vtkNodeNumbering );
      assert( vtkNodeNumbering.size() == (*it)->Nodes() );
      for( auto i = 0; i < vtkNodeNumbering.size(); ++i, ++entriesOfLine )
      {
         // writing node id to vtu document
        stringNumber = number_to_string( vtkNodeNumbering[i] );
        connectivityFile.InsertData( stringNumber.c_str() );
        // line break if limeit entries reached, else tab
        if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
        {
          entriesOfLine = 1;
          connectivityFile.LineBreak();
          connectivityFile.BringToLevel();
          newLine = true;
        }
        else
        {
          connectivityFile.Tab();
          newLine = false;
        }
      }
      // proceed to next element
      continue;
    }
    // for all another element types: looping element's nodes
    for( auto iit = 0; iit < (*it)->Nodes(); ++iit, ++entriesOfLine )
    {
      // writing node id to vtu document
      stringNumber = number_to_string( (*it)->N(iit)->Idx() );
      connectivityFile.InsertData( stringNumber.c_str() );
      // line break if limeit entries reached, else tab
      if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
      {
        entriesOfLine = 1;
        connectivityFile.LineBreak();
        connectivityFile.BringToLevel();
        newLine = true;
      }
      else
      {
        connectivityFile.Tab();
        newLine = false;
      }
    } // element's nodes
  } // region's elements
  // close cell connectivity
  if( !newLine )
    connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  // cell offsets
  connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  // looping region's elements
  size_t offset( 0 ); entriesOfLine = 2;
  for( auto it = subDomain.ElementsBegin(); it != domainSimplicesEnd; ++it, ++entriesOfLine )
  {
    // imcrementing offset by node count and writing to data
    offset += (*it)->Nodes();
    stringNumber = number_to_string( offset );
    connectivityFile.InsertData( stringNumber.c_str() );
    // line break if limeit entries reached, else tab
    if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
    {
      entriesOfLine = 1;
      connectivityFile.LineBreak();
      connectivityFile.BringToLevel();
      newLine = true;
    }
    else
    {
      connectivityFile.Tab();
      newLine = false;
    }
  } // region's elements
  // close offsets
  if( !newLine )
    connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  // cell types
  assert( elementTypesVTK.size() == subDomain.Elements() );
  connectivityFile.OpenNode( "DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  // looping region's elements
  entriesOfLine = 2; size_t elementCount( 0 );
  for( auto it = subDomain.ElementsBegin(); it != domainSimplicesEnd; ++it, ++entriesOfLine, ++elementCount )
  {
    // imcrementing offset by node count and writing to data
    stringNumber = number_to_string( elementTypesVTK[elementCount] );
    connectivityFile.InsertData( stringNumber.c_str() );
    // line break if limeit entries reached, else tab
    if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
    {
      entriesOfLine = 1;
      connectivityFile.LineBreak();
      connectivityFile.BringToLevel();
      newLine = true;
    }
    else
    {
      connectivityFile.Tab();
      newLine = false;
    }
  } // region's elements
  // close offsets
  if( !newLine )
    connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  // closing cells node
  connectivityFile.CloseNode( "Cells" );

  // io feedback
  if ( csmp_error.Verbose() )
      cout << " done!\n";
}

template void VTU_Interface<1U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFile(XML_Document&,const ModelSubDomain<3U,InterFace>&) const;





/// establishes the connectivity file for given region using a finite element barycenter point cloud
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::EstablishConnectivityFileBCPC( XML_Document& connectivityFile, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  ErrorHandler& csmp_error ( ErrorHandler::Instance() );

  // io feedback
  if( csmp_error.Verbose() )
      cout << "\nVTU_Interface<dim>::EstablishConnectivityFile: 'Writing bary center point cloud connectivity file for sub domain " << DomainName( subDomain ) << " ...";

  // getting node and element count, creating working strings
  string stringNumber, stringCache;

  // opening piece node, number of nodes == number of elements (since berycenters)
  const size_t DOMAIN_ELEMENTS( subDomain.Elements() );
  stringNumber = number_to_string( DOMAIN_ELEMENTS );
  stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

  // no elements
  stringCache += "\" NumberOfCells=\"0\"";
  connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

  // points node
  connectivityFile.OpenNode( "Points" );
  connectivityFile.OpenNode( "DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
  // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
  size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
  bool newLine( false );
  connectivityFile.BringToLevel();
  // looping over all the region's nodes
  const auto domainSimpicesEnd( subDomain.ElementsEnd() );
  for( auto it = subDomain.ElementsBegin(); it != domainSimpicesEnd; ++it, ++entriesOfLine )
    {
    // writing x,y and z coordinates(tab seperated)
    stringNumber = number_to_string( (*it)->BaryCenter()[0] );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = number_to_string( yCoordinate( (*it)->BaryCenter() ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = number_to_string( zCoordinate( (*it)->BaryCenter() ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    // if maximum of entris per line is reached, start new one, else tab
    if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE )
      {
      entriesOfLine = 1;
      connectivityFile.LineBreak();
      connectivityFile.BringToLevel();
      newLine = true;
      }
    else
      {
      connectivityFile.Tab();
      newLine = false;
      }
    } // looping regions nodes
  // closing points node
  if( !newLine )
    connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  connectivityFile.CloseNode( "Points" );

  // cells node
  connectivityFile.OpenNode( "Cells" );
  // cell connectivity
  connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  // looping region's elements
  const size_t MAX_CONNECTIVITY_ENTRIES_PER_LINE( 20 ); entriesOfLine = 2;
  for( size_t pid = 0; pid < DOMAIN_ELEMENTS; ++pid  )
  {
    // writing node id to vtu document
    stringNumber = number_to_string( pid );
    connectivityFile.InsertData( stringNumber.c_str() );
    // line break if limeit entries reached, else tab
    if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
    {
        entriesOfLine = 1;
        connectivityFile.LineBreak();
        connectivityFile.BringToLevel();
        newLine = true;
    }
    else
    {
        connectivityFile.Tab();
        newLine = false;
    }
  } // integration points
  // close cell connectivity
  if( !newLine )
    connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  // cell offsets
  stringNumber = number_to_string( DOMAIN_ELEMENTS );
  connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  connectivityFile.InsertData( stringNumber.c_str() );
  connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  // cell types
  size_t cell_type( VTK_POLY_VERTEX );
  stringNumber = number_to_string( cell_type );
  connectivityFile.OpenNode( "DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  connectivityFile.InsertData( stringNumber.c_str() );
  connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  // closing cells node
  connectivityFile.CloseNode( "Cells" );

  // io feedback
  if( csmp_error.Verbose() )
    cout << " done!\n";
}

template void VTU_Interface<1U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileBCPC(XML_Document&,const ModelSubDomain<3U,InterFace>&) const;




/// establishes the connectivity file for given region using a region barycenter point cloud
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::EstablishConnectivityFileRPC( XML_Document& connectivityFile, const ModelSubDomain<dim,CELL>& subDomain ) const
{
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    // io feedback
    if( csmp_error.Verbose() )
        cout << "\nVTU_Interface<dim>::EstablishConnectivityFile: 'Writing bary center point cloud connectivity file for sub domain " << DomainName( subDomain ) << " ...";

    // getting node and element count, creating working strings
    string stringNumber, stringCache;

    // opening piece node, number of nodes == 1U (single point per Region)
    const size_t REGION_POINTS( 1U );
    stringNumber = number_to_string( REGION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
    // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
    size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
    bool newLine( false );
    connectivityFile.BringToLevel();
    Point<dim> pbc,pt;
    Point<dim> pt_max,pt_min;
    Point<dim> dist_to_min,dist_to_max;
    double dist = numeric_limits<double>::max();
    double dist_bc;
    subDomain.MinMaxCoordinates(pt_min,pt_max);
    // looping over all the region's nodes
    const auto domainSimplicesEnd( subDomain.ElementsEnd() );
    for( auto it  = subDomain.ElementsBegin(); it != domainSimplicesEnd; ++it )
    {
        pbc = (*it)->BaryCenter();
        dist_to_min = pbc - pt_min;
        dist_to_max = pbc - pt_max;
        dist_bc = abs( dist_to_max.Length() - dist_to_min.Length() );
        if( dist_bc < dist )
        {
            pt = pbc;
            dist = dist_bc;
        }
    }
    ++entriesOfLine;
    // writing x,y and z coordinates(tab seperated)
    stringNumber = number_to_string( pt[0] );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = number_to_string( yCoordinate( pt ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = number_to_string( zCoordinate( pt ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    // if maximum of entries per line is reached, start new one, else tab
    if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE )
    {
        entriesOfLine = 1;
        connectivityFile.LineBreak();
        connectivityFile.BringToLevel();
        newLine = true;
    }
    else
    {
        connectivityFile.Tab();
        newLine = false;
    }
    // closing points node
    if( !newLine )
    connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    connectivityFile.CloseNode( "Points" );

    // cells node
    connectivityFile.OpenNode( "Cells" );
    // cell connectivity
    size_t pid( 0 );
    stringNumber = number_to_string( pid );
    connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
    connectivityFile.BringToLevel();
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    // cell offsets
    size_t offset(1);
    stringNumber = number_to_string( offset );
    connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
    connectivityFile.BringToLevel();
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    // cell types
    size_t cell_type( VTK_POLY_VERTEX ); // VTK_VERTEX would be also fine, but this way it's more general
    stringNumber = number_to_string( cell_type );
    connectivityFile.OpenNode( "DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
    connectivityFile.BringToLevel();
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    // closing cells node
    connectivityFile.CloseNode( "Cells" );

    // io feedback
    if( csmp_error.Verbose() )
    cout << " done!\n";
}

template void VTU_Interface<1U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileRPC(XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// establishes the connectivity file for given region using a finite element integration points
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::EstablishConnectivityFileFEIP( XML_Document& connectivityFile, const ModelSubDomain<dim,CELL>& subDomain ) const
{
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    // io feedback
    if( csmp_error.Verbose() )
        cout << "\nVTU_Interface<dim>::EstablishConnectivityFile: 'Writing bary center point cloud connectivity file for sub domain " << DomainName( subDomain ) << " ...";

    // getting node and element count, creating working strings
    string stringNumber, stringCache;

    // opening piece node, number of nodes == number of element integration points
    const size_t INTEGRATION_POINTS( subDomain.IntegrationPoints() );
    stringNumber = number_to_string( INTEGRATION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
    // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
    size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
    bool newLine( false );
    connectivityFile.BringToLevel();

    // looping over all element integration points
    const auto domainSimplicesEnd( subDomain.ElementsEnd() );
    for( auto it  = subDomain.ElementsBegin(); it != domainSimplicesEnd; ++it )
    {
        const auto ips = (*it)->FE()->IntegrationPoints();
        for( auto ip = 0; ip < ips; ++ip, ++entriesOfLine )
        {
          auto pt = (*it)->IntegrationPoint( ip );
          // writing x,y and z coordinates(tab seperated)
          stringNumber = number_to_string( pt[0] );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = number_to_string( yCoordinate( pt ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = number_to_string( zCoordinate( pt ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          // if maximum of entris per line is reached, start new one, else tab
          if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE )
          {
            entriesOfLine = 1;
            connectivityFile.LineBreak();
            connectivityFile.BringToLevel();
            newLine = true;
          }
          else
          {
            connectivityFile.Tab();
            newLine = false;
          }
        } // looping regions nodes
    }
    // closing points node
    if( !newLine )
      connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    connectivityFile.CloseNode( "Points" );

      // cells node
      connectivityFile.OpenNode( "Cells" );
      // cell connectivity
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      // looping region's elements
      const size_t MAX_CONNECTIVITY_ENTRIES_PER_LINE( 20 ); entriesOfLine = 2;
      for( size_t pid = 0; pid < INTEGRATION_POINTS; ++pid  )
      {
        // writing node id to vtu document
        stringNumber = number_to_string( pid );
        connectivityFile.InsertData( stringNumber.c_str() );
        // line break if limeit entries reached, else tab
        if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
        {
            entriesOfLine = 1;
            connectivityFile.LineBreak();
            connectivityFile.BringToLevel();
            newLine = true;
        }
        else
        {
            connectivityFile.Tab();
            newLine = false;
        }
      } // integration points
      // close cell connectivity
      if( !newLine )
        connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell offsets
      stringNumber = number_to_string( INTEGRATION_POINTS );
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell types
      size_t cell_type( VTK_POLY_VERTEX );
      stringNumber = number_to_string( cell_type );
      connectivityFile.OpenNode( "DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // closing cells node
      connectivityFile.CloseNode( "Cells" );

    // io feedback
    if( csmp_error.Verbose() )
      cout << " done!\n";
}
template void VTU_Interface<1U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFEIP(XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// establishes the connectivity file for given region using a finite volume integration points
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::EstablishConnectivityFileFVSIP( XML_Document& connectivityFile, const ModelSubDomain<dim,CELL>& subDomain ) const
{
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    // io feedback
    if( csmp_error.Verbose() )
        cout << "\nVTU_Interface<dim>::EstablishConnectivityFile: 'Writing bary center point cloud connectivity file for sub domain " << DomainName( subDomain ) << " ...";

    // getting node and element count, creating working strings
    string stringNumber, stringCache;

    // opening piece node, number of nodes == number of sector integration points
    const auto INTEGRATION_POINTS{ subDomain.SectorIntegrationPoints() };
    stringNumber = number_to_string( INTEGRATION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
    // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
    size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
    bool newLine( false );
    connectivityFile.BringToLevel();

    // looping over all the region's nodes
    const auto domainSimplicesEnd( subDomain.ElementsEnd() );
    for( auto it = subDomain.ElementsBegin(); it != domainSimplicesEnd; ++it )
    {
        auto sectors = (*it)->Sectors();
        for( auto sid = 0; sid < sectors; ++sid )
        {
            const auto ips = (*it)->FV()->IntegrationPointsPerSector( sid );
            for( auto ip = 0; ip < ips; ++ip, ++entriesOfLine )
            {
              auto pt = (*it)->FV()->SectorIntegrationPoint( sid, ip );
              pt = (*it)->RstToXYZ( pt );
              // writing x,y and z coordinates(tab seperated)
              stringNumber = number_to_string( pt[0] );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = number_to_string( yCoordinate( pt ) );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = number_to_string( zCoordinate( pt ) );
              connectivityFile.InsertData( stringNumber.c_str() );
              // if maximum of entris per line is reached, start new one, else tab
              if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE )
              {
                entriesOfLine = 1;
                connectivityFile.LineBreak();
                connectivityFile.BringToLevel();
                newLine = true;
              }
              else
              {
                connectivityFile.Tab();
                newLine = false;
              }
            }
        } // looping regions nodes
    }
    // closing points node
    if( !newLine )
      connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    connectivityFile.CloseNode( "Points" );

      // cells node
      connectivityFile.OpenNode( "Cells" );
      // cell connectivity
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      // looping region's elements
      const size_t MAX_CONNECTIVITY_ENTRIES_PER_LINE( 20 ); entriesOfLine = 2;
      for( size_t pid = 0; pid < INTEGRATION_POINTS; ++pid  )
      {
        // writing node id to vtu document
        stringNumber = number_to_string( pid );
        connectivityFile.InsertData( stringNumber.c_str() );
        // line break if limeit entries reached, else tab
        if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
        {
            entriesOfLine = 1;
            connectivityFile.LineBreak();
            connectivityFile.BringToLevel();
            newLine = true;
        }
        else
        {
            connectivityFile.Tab();
            newLine = false;
        }
      } // integration points
      // close cell connectivity
      if( !newLine )
        connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell offsets
      stringNumber = number_to_string( INTEGRATION_POINTS );
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell types
      size_t cell_type( VTK_POLY_VERTEX );
      stringNumber = number_to_string( cell_type );
      connectivityFile.OpenNode( "DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // closing cells node
      connectivityFile.CloseNode( "Cells" );

    // io feedback
    if( csmp_error.Verbose() )
      cout << " done!\n";
}

template void VTU_Interface<1U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFVSIP(XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// establishes the connectivity file for given region using a finite volume integration points
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::EstablishConnectivityFileFVFIP( XML_Document& connectivityFile, const ModelSubDomain<dim,CELL>& subDomain ) const
{
    ErrorHandler& csmp_error ( ErrorHandler::Instance() );

    // io feedback
    if( csmp_error.Verbose() )
        cout << "\nVTU_Interface<dim>::EstablishConnectivityFile: 'Writing bary center point cloud connectivity file for sub domain " << DomainName( subDomain ) << " ...";

    // getting node and element count, creating working strings
    string stringNumber, stringCache;

    // opening piece node, number of nodes == number of facet integration points
    const size_t INTEGRATION_POINTS( subDomain.FacetIntegrationPoints() );
    stringNumber = number_to_string( INTEGRATION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
    // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
    size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
    bool newLine( false );
    connectivityFile.BringToLevel();

    // looping over all the region's nodes
    const auto domainSimplicesEnd( subDomain.ElementsEnd() );
    for( auto it = subDomain.ElementsBegin(); it != domainSimplicesEnd; ++it )
    {
        const auto facets = (*it)->Facets();
        for( auto fid = 0; fid < facets; ++fid )
        {
            const auto ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
            for( auto ip = 0; ip < ips; ++ip, ++entriesOfLine )
            {
              auto pt = (*it)->FV()->FacetIntegrationPoint( fid, ip );
              pt = (*it)->RstToXYZ( pt );
              // writing x,y and z coordinates(tab seperated)
              stringNumber = number_to_string( pt[0] );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = number_to_string( yCoordinate( pt ) );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = number_to_string( zCoordinate( pt ) );
              connectivityFile.InsertData( stringNumber.c_str() );
              // if maximum of entris per line is reached, start new one, else tab
              if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE )
              {
                entriesOfLine = 1;
                connectivityFile.LineBreak();
                connectivityFile.BringToLevel();
                newLine = true;
              }
              else
              {
                connectivityFile.Tab();
                newLine = false;
              }
            }
        } // looping regions nodes
    }
    // closing points node
    if( !newLine )
      connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    connectivityFile.CloseNode( "Points" );

      // cells node
      connectivityFile.OpenNode( "Cells" );
      // cell connectivity
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      // looping region's elements
      const size_t MAX_CONNECTIVITY_ENTRIES_PER_LINE( 20 ); entriesOfLine = 2;
      for( size_t pid = 0; pid < INTEGRATION_POINTS; ++pid  )
      {
        // writing node id to vtu document
        stringNumber = number_to_string( pid );
        connectivityFile.InsertData( stringNumber.c_str() );
        // line break if limeit entries reached, else tab
        if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE )
        {
            entriesOfLine = 1;
            connectivityFile.LineBreak();
            connectivityFile.BringToLevel();
            newLine = true;
        }
        else
        {
            connectivityFile.Tab();
            newLine = false;
        }
      } // integration points
      // close cell connectivity
      if( !newLine )
        connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell offsets
      stringNumber = number_to_string( INTEGRATION_POINTS );
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell types
      size_t cell_type( VTK_POLY_VERTEX );
      stringNumber = number_to_string( cell_type );
      connectivityFile.OpenNode( "DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // closing cells node
      connectivityFile.CloseNode( "Cells" );

    // io feedback
    if( csmp_error.Verbose() )
      cout << " done!\n";
}

template void VTU_Interface<1U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<1U,Element>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<2U,Element>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<3U,Element>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<1U,Face>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<2U,Face>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<3U,Face>&) const;

template void VTU_Interface<1U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<1U,InterFace>&) const;
template void VTU_Interface<2U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<2U,InterFace>&) const;
template void VTU_Interface<3U>::EstablishConnectivityFileFVFIP(XML_Document&,const ModelSubDomain<3U,InterFace>&) const;


/// converts csmp to vtk node numbering
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::QuadraticHexahedronConnectivity( const CELL<dim>* const element, vector<long>& data ) const
{
  // first 12 nodes are the same
  for( uint32_t i = 0; i < 12; ++i )
    data.push_back( element->N( i )->Idx() );
  // now the node convention differs
  data.push_back( element->N( static_cast<uint32_t>(16) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(17) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(18) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(19) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(12) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(13) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(14) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(15) )->Idx() );
}

template void VTU_Interface<1U>::QuadraticHexahedronConnectivity(const Element<1U>* const,vector<long>&) const;
template void VTU_Interface<2U>::QuadraticHexahedronConnectivity(const Element<2U>* const,vector<long>&) const;
template void VTU_Interface<3U>::QuadraticHexahedronConnectivity(const Element<3U>* const,vector<long>&) const;

template void VTU_Interface<1U>::QuadraticHexahedronConnectivity(const Face<1U>* const,vector<long>&) const;
template void VTU_Interface<2U>::QuadraticHexahedronConnectivity(const Face<2U>* const,vector<long>&) const;
template void VTU_Interface<3U>::QuadraticHexahedronConnectivity(const Face<3U>* const,vector<long>&) const;

template void VTU_Interface<1U>::QuadraticHexahedronConnectivity(const InterFace<1U>* const,vector<long>&) const;
template void VTU_Interface<2U>::QuadraticHexahedronConnectivity(const InterFace<2U>* const,vector<long>&) const;
template void VTU_Interface<3U>::QuadraticHexahedronConnectivity(const InterFace<3U>* const,vector<long>&) const;

/// converts csmp to vtk node numbering
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::QuadraticWedgeConnectivity( const CELL<dim>* const element, vector<long>& data ) const
{
  // first 8 nodes are the same
  for( auto i = 0; i < 9; ++i )
    data.push_back( element->N( i )->Idx() );
  // now the node convention differs
  data.push_back( element->N( static_cast<uint32_t>(12) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(13) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(14) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(9) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(10) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(11) )->Idx() );
}

template void VTU_Interface<1U>::QuadraticWedgeConnectivity(const Element<1U>* const,vector<long>&) const;
template void VTU_Interface<2U>::QuadraticWedgeConnectivity(const Element<2U>* const,vector<long>&) const;
template void VTU_Interface<3U>::QuadraticWedgeConnectivity(const Element<3U>* const,vector<long>&) const;

template void VTU_Interface<1U>::QuadraticWedgeConnectivity(const Face<1U>* const,vector<long>&) const;
template void VTU_Interface<2U>::QuadraticWedgeConnectivity(const Face<2U>* const,vector<long>&) const;
template void VTU_Interface<3U>::QuadraticWedgeConnectivity(const Face<3U>* const,vector<long>&) const;

template void VTU_Interface<1U>::QuadraticWedgeConnectivity(const InterFace<1U>* const,vector<long>&) const;
template void VTU_Interface<2U>::QuadraticWedgeConnectivity(const InterFace<2U>* const,vector<long>&) const;
template void VTU_Interface<3U>::QuadraticWedgeConnectivity(const InterFace<3U>* const,vector<long>&) const;


/// finds the corresponding VTK Element type for given csmp::Element
template<uint32_t dim>
template<template <uint32_t> class CELL>
VTK_TYPE VTU_Interface<dim>::ElementType( const CELL<dim>* const elmt ) const
{
  CSMP_FEM_TYPE csmpType( elmt->FE_Type() );

  switch( csmpType )
  {
  case ISOPARAMETRIC_LINEAR_BAR:
  case LINEAR_BAR:                      return VTK_LINE;
                                        break;
  case LINEAR_TRIANGLE3D:
  case ISOPARAMETRIC_LINEAR_TRIANGLE:
  case LINEAR_TRIANGLE:                 return VTK_TRIANGLE;
                                        break;
  case ISOPARAMETRIC_LINEAR_TETRAHEDRON:
  case LINEAR_TETRAHEDRON:              return VTK_TETRA;
                                        break;
  case ISOPARAMETRIC_LINEAR_PYRAMID:    return VTK_PYRAMID;
                                        break;
  case ISOPARAMETRIC_LINEAR_PRISM:      return VTK_WEDGE;
                                        break;
  case ISOPARAMETRIC_LINEAR_QUADRILATERAL:
  case LINEAR_RECTANGLE:
                                        return VTK_QUAD;
                                        break;
  case ISOPARAMETRIC_LINEAR_HEXAHEDRON:
  case LINEAR_CUBOID:
                                        return VTK_HEXAHEDRON;
                                        break;
  case ISOPARAMETRIC_QUADRATIC_BAR:
  case QUADRATIC_BAR:                   return VTK_QUADRATIC_EDGE;
                                        break;
  case ISOPARAMETRIC_QUADRATIC_TRIANGLE:
  case QUADRATIC_TRIANGLE:              return VTK_QUADRATIC_TRIANGLE;
                                        break;
  case ISOPARAMETRIC_QUADRATIC_QUADRILATERAL:
                                        return VTK_QUADRATIC_QUAD;
                                        break;
  case ISOPARAMETRIC_QUADRATIC_TETRAHEDRON:
  case QUADRATIC_TETRAHEDRON:           return VTK_QUADRATIC_TETRA;
                                        break;
  case ISOPARAMETRIC_QUADRATIC_HEXAHEDRON20:
                                        return VTK_QUADRATIC_HEXAHEDRON;
                                        break;
  case ISOPARAMETRIC_QUADRATIC_PRISM15:
                                        return VTK_QUADRATIC_WEDGE;
                                        break;
  case ISOPARAMETRIC_QUADRATIC_PYRAMID13:
                                        return VTK_QUADRATIC_PYRAMID;
                                        break;


  default:                              throw csmp::Exception( ERROR,
                                                               "VTU_Interface<dim>::ElementType",
                                                               "Unable to identify element type",
                                                               "unknown" );
  }
}

template VTK_TYPE VTU_Interface<1U>::ElementType(const Element<1U>* const) const;
template VTK_TYPE VTU_Interface<2U>::ElementType(const Element<2U>* const) const;
template VTK_TYPE VTU_Interface<3U>::ElementType(const Element<3U>* const) const;

template VTK_TYPE VTU_Interface<1U>::ElementType(const Face<1U>* const) const;
template VTK_TYPE VTU_Interface<2U>::ElementType(const Face<2U>* const) const;
template VTK_TYPE VTU_Interface<3U>::ElementType(const Face<3U>* const) const;

template VTK_TYPE VTU_Interface<1U>::ElementType(const InterFace<1U>* const) const;
template VTK_TYPE VTU_Interface<2U>::ElementType(const InterFace<2U>* const) const;
template VTK_TYPE VTU_Interface<3U>::ElementType(const InterFace<3U>* const) const;

// domain name
template<uint32_t dim>
string VTU_Interface<dim>::DomainName( const ModelSubDomain<dim,Element>& subDomain ) const
{
  for( typename map<string,csmp::Region<dim> >::const_iterator it = model_.UniqueRegionsBegin(); it != model_.UniqueRegionsEnd(); ++it )
    if( &it->second == &subDomain )
      return string( it->first + " (unique REGION)" );
  for( typename map<string,csmp::Region<dim> >::const_iterator it = model_.RegionsBegin(); it != model_.RegionsEnd(); ++it )
    if( &it->second == &subDomain )
      return string( it->first + " (REGION)" );

   return "UndefinedRegion";
 }

template<uint32_t dim>
string VTU_Interface<dim>::DomainName( const ModelSubDomain<dim,Face>& subDomain ) const
{
  for( typename map<string,Boundary<dim> >::const_iterator it = model_.BoundariesBegin(); it != model_.BoundariesEnd(); ++it )
    if( &it->second == &subDomain )
      return string( it->first + " (BOUNDARY)" );

   return "UndefinedBoundary";
 }

template<uint32_t dim>
string VTU_Interface<dim>::DomainName( const ModelSubDomain<dim,InterFace>& subDomain ) const
{
  for( typename map<string,SplitBoundary<dim> >::const_iterator it = model_.SplitBoundariesBegin(); it != model_.SplitBoundariesEnd(); ++it )
    if( &it->second == &subDomain )
      return string( it->first + " (SPLITBOUNDARY)" );

   return "UndefinedSplitBoundary";
 }

// VTU always uses 3 dimensional output
double zCoordinate( Point<3> const& p )
{ return p[2]; }
double zCoordinate( Point<2> const& )
{ return 0.; }
double zCoordinate( Point<1> const& )
{ return 0.; }
double yCoordinate( Point<3> const& p )
{ return p[1]; }
double yCoordinate( Point<2> const& p )
{ return p[1]; }
double yCoordinate( Point<1> const& )
{ return 0.; }



/**
    Reads binary model file and writes target variable to VTU file.
 
    @author SKM
    @date 2008
*/
void outputBoundariesToVTU( const char* modelBinFIleName )
{
  const string bin_file_set(modelBinFIleName);
  Model<3> model( bin_file_set );
  VTU_Interface<3> vtu( model );

  cin.ignore( numeric_limits<streamsize>::max(), '\n' );
  cout << "\n\noutputBoundariesToVTU::Enter name of boundary property to output: " << flush;
  string propertyName( "" );
  getline( cin, propertyName );

  for ( auto it = model.BoundariesBegin(); it != model.BoundariesEnd(); it++ )
    //for ( map<boundaryName,Boundary<3U> >::const_iterator
    //it=model.BoundariesBegin(); it!=model.BoundariesEnd(); it++ )
    vtu.OutputDataToVTU( modelBinFIleName,     // file name
                         propertyName.c_str(), // propertyName
                         (*it).second, 0. );   // ModelSubDomain<dim,CELL>& subDomain
}





/**
    Reads binary model file and writes desired variable to VTU file.
 
    Interactive: user is prompted for the name of the property.
 
    @author SKM
    @date 2008
 
*/
void csmpBinaryToVTU( const char* modelBinFIleName )
{
  cin.ignore( numeric_limits<streamsize>::max(), '\n' );
  cout << "\n\ncsmpBinaryToVTU::Enter name of property to output: " << flush;
  string propertyName( "" );
  getline( cin, propertyName );
  
  const string bin_file_set(modelBinFIleName);
  Model<3> model( bin_file_set );
  VTU_Interface<3> vtu( model );

  vtu.OutputDataToVTU( propertyName, propertyName, string( "Model" ), 0 );
}




template class VTU_Interface<1U>;
template class VTU_Interface<2U>;
template class VTU_Interface<3U>;

} // csmp
