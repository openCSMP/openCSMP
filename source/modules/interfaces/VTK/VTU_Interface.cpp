// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "VTU_Interface.h"
#include "VTK_Interface.h"

#include "Boundary.h"
#include "SplitBoundary.h"
#include "Region.h"
#include "Model.h"
#include "Node.h"
#include "Element.h"
#include "Face.h"
#include "InterFace.h"
#include "OS_Utilities.h"
#include "CSMP_highLevelUtilities.h"
#include "PL_Utilities.h"
#include "ErrorHandler.h"
#include "Exception.h"

using namespace std;

namespace csmp {

/**
    constructor taking a reference to the model to be handled(and an optional title to the output)
     
    @todo introduce option to choose XML type:
    
    currently:  "<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian">"
    
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

template<uint32_t dim>
void VTU_Interface<dim>::DeleteConnectivity() noexcept
{
    // Field, node and element data
    clearConnectivityMap( regionConnectivityFiles_ );
    clearConnectivityMap( boundaryConnectivityFiles_ );
    clearConnectivityMap( splitBoundaryConnectivityFiles_ );

    // Element barycentre data
    clearConnectivityMap( regionConnectivityFiles_bcd_ );
    clearConnectivityMap( boundaryConnectivityFiles_bcd_ );
    clearConnectivityMap( splitBoundaryConnectivityFiles_bcd_ );

    // Region data
    clearConnectivityMap( regionConnectivityFiles_rcd_ );
    clearConnectivityMap( boundaryConnectivityFiles_rcd_ );
    clearConnectivityMap( splitBoundaryConnectivityFiles_rcd_ );

    // Finite Element integration point data
    clearConnectivityMap( regionConnectivityFiles_feipsd_ );
    clearConnectivityMap( boundaryConnectivityFiles_feipsd_ );
    clearConnectivityMap( splitBoundaryConnectivityFiles_feipsd_ );

    // Finite Volume sector integration point data
    clearConnectivityMap( regionConnectivityFiles_fvsipsd_ );
    clearConnectivityMap( boundaryConnectivityFiles_fvsipsd_ );
    clearConnectivityMap( splitBoundaryConnectivityFiles_fvsipsd_ );

    // Finite Volume facet integration point data
    clearConnectivityMap( regionConnectivityFiles_fvfipsd_ );
    clearConnectivityMap( boundaryConnectivityFiles_fvfipsd_ );
    clearConnectivityMap( splitBoundaryConnectivityFiles_fvfipsd_ );

    // Multiblock dataset index files
    clearConnectivityMap( regionConnectivityFiles_multiblock_ );
    clearConnectivityMap( boundaryConnectivityFiles_multiblock_ );
    clearConnectivityMap( splitBoundaryConnectivityFiles_multiblock_ );
}


// ============================================================
// GetConnectivityMap — field, node and element data
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
typename VTU_Interface<dim>::template ConnectivityMap<CELL>&
VTU_Interface<dim>::GetConnectivityMap( const ModelSubDomain<dim,CELL>& )
{
    if constexpr ( std::is_same_v<CELL<dim>, Element<dim>> )
        return regionConnectivityFiles_;
    else if constexpr ( std::is_same_v<CELL<dim>, Face<dim>> )
        return boundaryConnectivityFiles_;
    else
        return splitBoundaryConnectivityFiles_;
}

// ============================================================
// GetConnectivityMapBCS — element barycentre data
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
typename VTU_Interface<dim>::template ConnectivityMap<CELL>&
VTU_Interface<dim>::GetConnectivityMapBCS( const ModelSubDomain<dim,CELL>& )
{
    if constexpr ( std::is_same_v<CELL<dim>, Element<dim>> )
        return regionConnectivityFiles_bcd_;
    else if constexpr ( std::is_same_v<CELL<dim>, Face<dim>> )
        return boundaryConnectivityFiles_bcd_;
    else
        return splitBoundaryConnectivityFiles_bcd_;
}

// ============================================================
// GetConnectivityMapRCS — region data
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
typename VTU_Interface<dim>::template ConnectivityMap<CELL>&
VTU_Interface<dim>::GetConnectivityMapRCS( const ModelSubDomain<dim,CELL>& )
{
    if constexpr ( std::is_same_v<CELL<dim>, Element<dim>> )
        return regionConnectivityFiles_rcd_;
    else if constexpr ( std::is_same_v<CELL<dim>, Face<dim>> )
        return boundaryConnectivityFiles_rcd_;
    else
        return splitBoundaryConnectivityFiles_rcd_;
}

// ============================================================
// GetConnectivityMapFEIPS — finite element integration point data
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
typename VTU_Interface<dim>::template ConnectivityMap<CELL>&
VTU_Interface<dim>::GetConnectivityMapFEIPS( const ModelSubDomain<dim,CELL>& )
{
    if constexpr ( std::is_same_v<CELL<dim>, Element<dim>> )
        return regionConnectivityFiles_feipsd_;
    else if constexpr ( std::is_same_v<CELL<dim>, Face<dim>> )
        return boundaryConnectivityFiles_feipsd_;
    else
        return splitBoundaryConnectivityFiles_feipsd_;
}

// ============================================================
// GetConnectivityMapFVSIPS — finite volume sector integration point data
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
typename VTU_Interface<dim>::template ConnectivityMap<CELL>&
VTU_Interface<dim>::GetConnectivityMapFVSIPS( const ModelSubDomain<dim,CELL>& )
{
    if constexpr ( std::is_same_v<CELL<dim>, Element<dim>> )
        return regionConnectivityFiles_fvsipsd_;
    else if constexpr ( std::is_same_v<CELL<dim>, Face<dim>> )
        return boundaryConnectivityFiles_fvsipsd_;
    else
        return splitBoundaryConnectivityFiles_fvsipsd_;
}

// ============================================================
// GetConnectivityMapFVFIPS — finite volume facet integration point data
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
typename VTU_Interface<dim>::template ConnectivityMap<CELL>&
VTU_Interface<dim>::GetConnectivityMapFVFIPS( const ModelSubDomain<dim,CELL>& )
{
    if constexpr ( std::is_same_v<CELL<dim>, Element<dim>> )
        return regionConnectivityFiles_fvfipsd_;
    else if constexpr ( std::is_same_v<CELL<dim>, Face<dim>> )
        return boundaryConnectivityFiles_fvfipsd_;
    else
        return splitBoundaryConnectivityFiles_fvfipsd_;
}

// ============================================================
// GetConnectivityMapMultiBlock — multiblock dataset index files
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
typename VTU_Interface<dim>::template ConnectivityMap<CELL>&
VTU_Interface<dim>::GetConnectivityMapMultiBlock( const ModelSubDomain<dim,CELL>& )
{
    if constexpr ( std::is_same_v<CELL<dim>, Element<dim>> )
        return regionConnectivityFiles_multiblock_;
    else if constexpr ( std::is_same_v<CELL<dim>, Face<dim>> )
        return boundaryConnectivityFiles_multiblock_;
    else
        return splitBoundaryConnectivityFiles_multiblock_;
}




// INTERFACES
// ----------------

/** provide aliases for csmp variables that will be shown in vtu instead
    @return returns if existing was empty
*/
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


/** allows the creation of an index-to-name correspondance.  If the resulting container is
    empty or the array variable name does not exist upon the output call, the behaviour is the detault one (output of [#] suffix).
    If the component name exists,it gets prefixed to the variable name for the vtu file. Works for flagged arrays, too.
    Julian, July 2014
*/
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
    if( domainName != "Model" )
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
        numberString = numberToString( timestep );
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
  if( !model_.ContainsRegion( regionName ) ) return true;
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
  if( !model_.ContainsRegion( regionName ) ) return true;
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
  if( !model_.ContainsRegion( regionName ) ) return true;
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
    set<string> propertyNames( propertyNamesList.begin(), propertyNamesList.end() );
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
    
    for( auto it = propertyNames.begin(); it != propertyNames.end(); ++it )
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
void VTU_Interface<dim>::OutputMultiBlockVTU(
        const string& fileName,
        const vector<string>& fileNames,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Initialise output document from the cached header.
    // The multiblock document is not reused between calls, so the cache
    // serves only to provide the standard XML header consistently.
    XML_Document outputFile;
    outputFile = *ConnectivityFile<CELL>(
                     GetConnectivityMapMultiBlock( subDomain ), subDomain );

    // Root VTK element
    outputFile.OpenNode(
        "VTKFile type=\"vtkMultiBlockDataSet\" "
        "version=\"1.0\" byte_order=\"LittleEndian\"" );
    outputFile.OpenNode( "vtkMultiBlockDataSet" );
    outputFile.OpenNode( "Block index=\"0\" name=\"Blocks\"" );

    // One DataSet entry per constituent VTU file
    for ( size_t index = 0u; index < fileNames.size(); ++index )
    {
        const string dataSetAttr =
            "DataSet index=\""  + to_string( index )          + "\""
            " file=\""          + fileNames[index] + ".vtu"   + "\"";

        outputFile.OpenNode(  dataSetAttr.c_str() );
        outputFile.CloseNode( "DataSet" );
    }

    outputFile.CloseNode( "Block" );
    outputFile.CloseNode( "vtkMultiBlockDataSet" );
    outputFile.CloseNode( "VTKFile" );

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


/// Dispatches output of field, node and element scalar/vector/tensor/array
/// and flagged-array data to a single VTU file.
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::OutputFieldNodesAndElementDataToVTU(
        const string&              fileName,
        const list<Index>&         fieldDataIndices,
        const list<Index>&         nodeIndices,
        const list<Index>&         elementIndices,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Initialise output document from the cached header
    XML_Document outputFile;
    outputFile = *ConnectivityFile<CELL>(
                     GetConnectivityMap( subDomain ), subDomain );

    outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" "
                         "version=\"1.0\" byte_order=\"LittleEndian\"" );
    outputFile.OpenNode( "UnstructuredGrid" );

    // 1. Field data (model / region / boundary scalars, arrays, vectors, tensors)
    if ( !fieldDataIndices.empty() )
        OutputFieldDataToVTU<CELL>( outputFile, subDomain, fieldDataIndices );

    // 2. Node and element data
    const bool hasNodeData    = !nodeIndices.empty();
    const bool hasElementData = !elementIndices.empty();

    if ( hasNodeData || hasElementData )
    {
        // 2.1 Write Points, Cells connectivity — opens the Piece node
        EstablishConnectivityFile<CELL>( outputFile, subDomain );

        // 2.2 Point data arrays
        if ( hasNodeData )
            OutputPointDataToVTU( outputFile, subDomain, nodeIndices );

        // 2.3 Cell data arrays
        if ( hasElementData )
            OutputCellDataToVTU( outputFile, subDomain, elementIndices );

        // 2.4 Close Piece — must be called exactly once per EstablishConnectivityFile
        outputFile.CloseNode( "Piece" );
    }

    outputFile.CloseNode( "UnstructuredGrid" );
    outputFile.CloseNode( "VTKFile" );

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



/// Outputs element barycentre data as a VTK point-cloud VTU file.
/// Each element contributes one point at its barycentre; no cell connectivity
/// is written (NumberOfCells="0").
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::OutputElementBarycentricDataToVTU(
        const string&              fileName,
        const list<Index>&         elementMatrixIndices,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Initialise output document from the cached header
    XML_Document outputFile;
    outputFile = *ConnectivityFile<CELL>(
                     GetConnectivityMapBCS( subDomain ), subDomain );

    outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" "
                         "version=\"1.0\" byte_order=\"LittleEndian\"" );
    outputFile.OpenNode( "UnstructuredGrid" );

    // Barycentre point-cloud data
    // Guard is required: EstablishConnectivityFileBCPC opens the Piece node,
    // so CloseNode("Piece") must only be called when it was opened.
    if ( !elementMatrixIndices.empty() )
    {
        // Opens Piece node — must be closed below
        EstablishConnectivityFileBCPC<CELL>( outputFile, subDomain );
        OutputPointDataToVTU( outputFile, subDomain, elementMatrixIndices );
        outputFile.CloseNode( "Piece" ); // closes node opened by EstablishConnectivityFileBCPC
    }

    outputFile.CloseNode( "UnstructuredGrid" );
    outputFile.CloseNode( "VTKFile" );

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




/// Outputs region-placed data as a single-point VTU file.
/// Each region contributes one point at its geometric centre; no mesh
/// connectivity is written beyond a single VTK_POLY_VERTEX cell.
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::OutputRegionDataToVTU(
        const string&              fileName,
        const list<Index>&         regionIndices,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Initialise output document from the cached header
    XML_Document outputFile;
    outputFile = *ConnectivityFile<CELL>(
                     GetConnectivityMapRCS( subDomain ), subDomain );

    outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" "
                         "version=\"1.0\" byte_order=\"LittleEndian\"" );
    outputFile.OpenNode( "UnstructuredGrid" );

    // Region centre point-cloud data.
    // Guard is required: EstablishConnectivityFileRPC opens the Piece node,
    // so CloseNode("Piece") must only be called when it was opened.
    if ( !regionIndices.empty() )
    {
        // Opens Piece node — must be closed below
        EstablishConnectivityFileRPC<CELL>( outputFile, subDomain );
        OutputPointDataToVTU( outputFile, subDomain, regionIndices );
        outputFile.CloseNode( "Piece" ); // closes node opened by EstablishConnectivityFileRPC
    }

    outputFile.CloseNode( "UnstructuredGrid" );
    outputFile.CloseNode( "VTKFile" );

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




/// Outputs finite element integration point data as a VTK point-cloud VTU file.
/// Each integration point of each element in the subdomain contributes one point;
/// no cell connectivity is written beyond a single VTK_POLY_VERTEX cell.
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::OutputFiniteElementIntegrationPointsDataToVTU(
        const string&              fileName,
        const list<Index>&         feipIndices,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Initialise output document from the cached header
    XML_Document outputFile;
    outputFile = *ConnectivityFile<CELL>(
                     GetConnectivityMapFEIPS( subDomain ), subDomain );

    outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" "
                         "version=\"1.0\" byte_order=\"LittleEndian\"" );
    outputFile.OpenNode( "UnstructuredGrid" );

    // Finite element integration point data.
    // Guard is required: EstablishConnectivityFileFEIP opens the Piece node,
    // so CloseNode("Piece") must only be called when it was opened.
    if ( !feipIndices.empty() )
    {
        // Opens Piece node — must be closed below
        EstablishConnectivityFileFEIP<CELL>( outputFile, subDomain );
        OutputPointDataToVTU( outputFile, subDomain, feipIndices );
        outputFile.CloseNode( "Piece" ); // closes node opened by EstablishConnectivityFileFEIP
    }

    outputFile.CloseNode( "UnstructuredGrid" );
    outputFile.CloseNode( "VTKFile" );

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




/// Outputs finite volume sector integration point data as a VTK point-cloud VTU file.
/// Each sector integration point of each element in the subdomain contributes one point;
/// no cell connectivity is written beyond a single VTK_POLY_VERTEX cell.
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::OutputFiniteVolumeSectorIntegrationPointsDataToVTU(
        const string&              fileName,
        const list<Index>&         fvsipIndices,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Initialise output document from the cached header
    XML_Document outputFile;
    outputFile = *ConnectivityFile<CELL>(
                     GetConnectivityMapFVSIPS( subDomain ), subDomain );

    outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" "
                         "version=\"1.0\" byte_order=\"LittleEndian\"" );
    outputFile.OpenNode( "UnstructuredGrid" );

    // Finite volume sector integration point data.
    // Guard is required: EstablishConnectivityFileFVSIP opens the Piece node,
    // so CloseNode("Piece") must only be called when it was opened.
    if ( !fvsipIndices.empty() )
    {
        // Opens Piece node — must be closed below
        EstablishConnectivityFileFVSIP<CELL>( outputFile, subDomain );
        OutputPointDataToVTU( outputFile, subDomain, fvsipIndices );
        outputFile.CloseNode( "Piece" ); // closes node opened by EstablishConnectivityFileFVSIP
    }

    outputFile.CloseNode( "UnstructuredGrid" );
    outputFile.CloseNode( "VTKFile" );

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




/// Outputs finite volume facet integration point data as a VTK point-cloud VTU file.
/// Each facet integration point of each element in the subdomain contributes one point;
/// no cell connectivity is written beyond a single VTK_POLY_VERTEX cell.
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::OutputFiniteVolumeFacetIntegrationPointsDataToVTU(
        const string&              fileName,
        const list<Index>&         fvfipIndices,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Initialise output document from the cached header
    XML_Document outputFile;
    outputFile = *ConnectivityFile<CELL>(
                     GetConnectivityMapFVFIPS( subDomain ), subDomain );

    outputFile.OpenNode( "VTKFile type=\"UnstructuredGrid\" "
                         "version=\"1.0\" byte_order=\"LittleEndian\"" );
    outputFile.OpenNode( "UnstructuredGrid" );

    // Finite volume facet integration point data.
    // Guard is required: EstablishConnectivityFileFVFIP opens the Piece node,
    // so CloseNode("Piece") must only be called when it was opened.
    if ( !fvfipIndices.empty() )
    {
        // Opens Piece node — must be closed below
        EstablishConnectivityFileFVFIP<CELL>( outputFile, subDomain );
        OutputPointDataToVTU( outputFile, subDomain, fvfipIndices );
        outputFile.CloseNode( "Piece" ); // closes node opened by EstablishConnectivityFileFVFIP
    }

    outputFile.CloseNode( "UnstructuredGrid" );
    outputFile.CloseNode( "VTKFile" );

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
                                               const ModelSubDomain<dim,CELL>&,
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
  for( size_t i{0u}; i < vectors.size(); ++i )
    if( vectors.at( i ).size() != 3U )
      throw csmp::Exception( ERROR, "VTU_Interface<dim>::OutputDataToVTU", "Vector required to have 3 components." );

  string fullFileName( fileName );
  fullFileName.append( ".vtu" );
  ofstream file( fullFileName.data(), ios::out );

  // header
  file << "<?xml version=\"1.0\"?>" << endl << endl;
  file << "<!--\n" << fileName << "\n-->" << endl << endl;

  // body - connectivity
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
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
  for( size_t i{0u}; i < vectors.size(); ++i )
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
  file << "<?xml version=\"1.0\"?>" << endl << endl;
  file << "<!--\n" << fileName << "\n-->" << endl << endl;

  // body - connectivity
  file << "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">" << endl;
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
  for( size_t i = 0; i < 3; ++i )
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
    stringNumber = numberToString( scalarVariable );
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
    for( uint32_t i = 0u; i < dim; ++i )
    {
      stringNumber = numberToString( vectorVariable.Component( i ) );
      vtu.InsertData( stringNumber.c_str() );
      vtu.Tab();
    }
    if constexpr ( dim != 3U )
      {
          for( uint32_t i = dim; i < 3U; ++i )
          {
            stringNumber = numberToString( 0. );
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
    for( uint32_t column = 0; column < cols; ++column )
    {
      for( uint32_t row = 0; row < rows; ++row )
      {
        val = ( ( row < dim ) && ( column < dim ) ? tensorVariable( row, column ) : 0.0 );
        stringNumber = numberToString( val );
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
    stringNumber = to_string( varSize );

    string arrayTitle( "DataArray type=\"Float64\" Name=\"" );
    arrayTitle += variableName;
    arrayTitle += "\" NumberOfTuples=\"";
    arrayTitle += stringNumber;
    arrayTitle += "\" format=\"ascii\"";
    vtu.OpenNode( arrayTitle.c_str() );
    vtu.BringToLevel();
    bool newLine( false );

    for( uint32_t i=0u ; i<varSize; ++i )
    {
        // inserting data
        stringNumber = numberToString( var.Component( i ) );
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
  bool newLine( false );
  if( key.place == NODE )
  {
      if constexpr ( !std::is_same_v<CELL<dim>,InterFace<dim>> ) {
          const auto domainNodesEnd( subDomain.NodesEnd() );
          // looping over nodes of Region or Boundary
          for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
            // inserting scalar data
            WriteScalar( vtu, MAX_ENTRIES_PER_LINE, (*it)->Read(key), entriesOfLine, newLine );
     }
     else { // SPLITBOUNDARY
         // inside nodes
         for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).InsideNodes().first )
           // inserting scalar data
           WriteScalar(vtu,MAX_ENTRIES_PER_LINE,nit->Read(key),entriesOfLine,newLine);
         // outside nodes
         for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).OutsideNodes().first )
           // inserting scalar data
           WriteScalar(vtu,MAX_ENTRIES_PER_LINE,nit->Read(key),entriesOfLine,newLine);

       }
  }
  else if( key.place == MODEL )
  {
      WriteScalar( vtu, MAX_ENTRIES_PER_LINE, model_.Read(key), entriesOfLine, newLine );
  }
  else if( key.place == REGION || key.place == BOUNDARY || key.place == SPLIT_BOUNDARY )
  {
      ScalarVariable scalarVariable;
      if ( key.place == REGION ) dynamic_cast<const Region<dim>&>(subDomain).Read( key, scalarVariable );
      else if ( key.place == BOUNDARY ) dynamic_cast<const Boundary<dim>&>(subDomain).Read( key, scalarVariable );
      else if ( key.place == SPLIT_BOUNDARY ) dynamic_cast<const SplitBoundary<dim>&>(subDomain).Read( key, scalarVariable );

      WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
  }
  else if( key.place == ELEMENT_INTEGRATION_POINT  || key.place == FACE_INTEGRATION_POINT  || key.place == INTER_FACE_INTEGRATION_POINT )
  {
      ScalarVariable scalarVariable;
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
      {
          const uint32_t ips = (*it)->IntegrationPoints();
          for( uint32_t ip = 0u; ip < ips; ++ip )
          {
              // inserting scalar data
              (*it)->Read( ip, key, scalarVariable );
              WriteScalar(vtu,MAX_ENTRIES_PER_LINE,scalarVariable(),entriesOfLine,newLine);
          }
      }
  }
  else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
  {
      ScalarVariable scalarVariable;
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
      {
          const uint32_t sectors = (*it)->Sectors();
          for( uint32_t sid = 0u; sid < sectors; ++sid )
          {
              const uint32_t ips = (*it)->FV()->IntegrationPointsPerSector( sid );
              for( uint32_t ip = 0u; ip < ips; ++ip )
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
      ScalarVariable scalarVariable;
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
      {
          const uint32_t facets = (*it)->Facets();
          for( uint32_t fid = 0u; fid < facets; ++fid )
          {
              const uint32_t ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
              for( uint32_t  ip = 0u; ip < ips; ++ip )
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
      if constexpr ( !std::is_same_v<CELL<dim>,InterFace<dim>> ) {
          // looping over nodes of Region or Boundary
          const auto domainNodesEnd( subDomain.NodesEnd() );
          for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
          {
            // acquiring vector data
            (*it)->Read( key, vectorVariable );
            WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
          }
     }
     else { // SPLITBOUNDARY
         // inside nodes
         for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).InsideNodes().first ) {
              nit->Read( key, vectorVariable );
              WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
           }
         // outside nodes
         for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).OutsideNodes().first ) {
              nit->Read( key, vectorVariable );
              WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
           }
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
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin();
           it != domainElementsEnd; ++it )
      {
          const uint32_t  ips = (*it)->IntegrationPoints();
          for( uint32_t  ip = 0u; ip < ips; ++ip )
          {
              // inserting scalar data
              (*it)->Read( ip, key, vectorVariable );
              WriteVector(vtu,MAX_ENTRIES_PER_LINE,vectorVariable,entriesOfLine,newLine);
          }
      }
  }
  else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin();
           it != domainElementsEnd; ++it )
      {
          const uint32_t sectors = (*it)->Sectors();
          for( uint32_t sid = 0u; sid < sectors; ++sid )
          {
              const uint32_t ips = (*it)->FV()->IntegrationPointsPerSector( sid );
              for( uint32_t ip = 0u; ip < ips; ++ip )
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
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin();
           it != domainElementsEnd; ++it )
      {
          const uint32_t facets = (*it)->Facets();
          for( uint32_t fid = 0u; fid < facets; ++fid )
          {
              const uint32_t ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
              for( uint32_t ip = 0u; ip < ips; ++ip )
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
  stringNumber = to_string( tensor_components );
  arrayTitle += stringNumber; arrayTitle += "\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  bool newLine( false );
  TensorVariable<dim> tensorVariable;
  if( key.place == NODE )
  {
      if constexpr ( !std::is_same_v<CELL<dim>,InterFace<dim>> ) {
          // looping over nodes of Region or Boundary
          const auto domainNodesEnd( subDomain.NodesEnd() );
          for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
            {
              // acquiring vector data
              (*it)->Read( key, tensorVariable );
              WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
            }
         }
      else { // SPLITBOUNDARY
         // inside nodes
         for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).InsideNodes().first ) {
              nit->Read( key, tensorVariable );
              WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
           }
         // outside nodes
         for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).OutsideNodes().first ) {
              nit->Read( key, tensorVariable );
              WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
           }
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
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin();
           it != domainElementsEnd; ++it )
      {
          const uint32_t ips = (*it)->IntegrationPoints();
          for( uint32_t ip = 0u; ip < ips; ++ip )
          {
              // inserting scalar data
              (*it)->Read( ip, key, tensorVariable );
              WriteTensor(vtu,MAX_ENTRIES_PER_LINE,tensorVariable,entriesOfLine,newLine);
          }
      }
  }
  else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
  {
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin();
           it != domainElementsEnd; ++it )
      {
          const uint32_t sectors = (*it)->Sectors();
          for( uint32_t sid = 0u; sid < sectors; ++sid )
          {
              const uint32_t ips = (*it)->FV()->IntegrationPointsPerSector( sid );
              for( uint32_t ip = 0u; ip < ips; ++ip )
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
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin();
           it != domainElementsEnd; ++it )
      {
          const uint32_t facets = (*it)->Facets();
          for( uint32_t fid = 0u; fid < facets; ++fid )
          {
              const uint32_t ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
              for( uint32_t ip = 0u; ip < ips; ++ip )
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
    for( uint32_t component=0u; component<key.dataDepth; component++)
    {
        arrayTitle = "DataArray type=\"Float64\" Name=\"";
        variableName = model_.Database().Name( key );
        variableName = FindVariableOutputAlias( variableName );
        stringNumber = to_string( component );
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
        if constexpr ( !std::is_same_v<CELL<dim>,InterFace<dim>> ) {
            // looping over nodes of Region or Boundary
            const auto domainNodesEnd( subDomain.NodesEnd() );
            for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
              {
                  // inserting scalar array data
                  (*it)->Read( key, arrayVariable );
                  WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
              }
         }
        else { // SPLITBOUNDARY
           // inside nodes
           for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).InsideNodes().first ) {
                nit->Read( key, arrayVariable );
                WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
             }
           // outside nodes
           for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).OutsideNodes().first ) {
                nit->Read( key, arrayVariable );
                WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
             }
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
          const auto domainElementsEnd( subDomain.CellsEnd() );
          for( auto it = subDomain.CellsBegin();
               it != domainElementsEnd; ++it )
          {
              const uint32_t ips = (*it)->IntegrationPoints();
              for( uint32_t ip = 0u; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( ip, key, arrayVariable );
                  WriteScalar(vtu,MAX_ENTRIES_PER_LINE,arrayVariable[component],entriesOfLine,newLine);
              }
          }
      }
      else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.CellsEnd() );
          for( auto it = subDomain.CellsBegin();
               it != domainElementsEnd; ++it )
          {
              const uint32_t sectors = (*it)->Sectors();
              for( uint32_t sid = 0u; sid < sectors; ++sid )
              {
                  const uint32_t ips = (*it)->FV()->IntegrationPointsPerSector( sid );
                  for( uint32_t ip = 0u; ip < ips; ++ip )
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
          const auto domainElementsEnd( subDomain.CellsEnd() );
          for( auto it = subDomain.CellsBegin();
               it != domainElementsEnd; ++it )
          {
              const uint32_t facets = (*it)->Facets();
              for( uint32_t fid = 0u; fid < facets; ++fid )
              {
                  const uint32_t ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
                  for( uint32_t ip = 0u; ip < ips; ++ip )
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
    for( uint32_t component=0u; component<key.dataDepth; component++)
    {
        arrayTitle = "DataArray type=\"Float64\" Name=\"";
        variableName = model_.Database().Name( key );
        variableName = FindVariableOutputAlias( variableName );
        stringNumber = to_string( component );
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
        if constexpr ( !std::is_same_v<CELL<dim>,InterFace<dim>> ) {
            // looping over nodes of Region or Boundary
            const auto domainNodesEnd( subDomain.NodesEnd() );
            for( auto it = subDomain.NodesBegin(); it != domainNodesEnd; ++it )
              {
                // inserting scalar flagged array data
                (*it)->Read( key, flaggedArrayVariable );
                WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
              }
         }
        else { // SPLITBOUNDARY
           // inside nodes
           for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).InsideNodes().first ) {
                nit->Read( key, flaggedArrayVariable );
                WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
             }
           // outside nodes
           for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).OutsideNodes().first ) {
                nit->Read( key, flaggedArrayVariable );
                WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
             }
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
          const auto domainElementsEnd( subDomain.CellsEnd() );
          for( auto it = subDomain.CellsBegin();
               it != domainElementsEnd; ++it )
          {
              const uint32_t ips = (*it)->IntegrationPoints();
              for( uint32_t ip = 0u; ip < ips; ++ip )
              {
                  // inserting scalar data
                  (*it)->Read( ip, key, flaggedArrayVariable );
                  WriteScalar(vtu,MAX_ENTRIES_PER_LINE,flaggedArrayVariable[component],entriesOfLine,newLine);
              }
          }
      }
      else if( key.place == SECTOR_INTEGRATION_POINT  || key.place == FACE_SECTOR_INTEGRATION_POINT  || key.place == INTER_FACE_SECTOR_INTEGRATION_POINT )
      {
          const auto domainElementsEnd( subDomain.CellsEnd() );
          for( auto it = subDomain.CellsBegin();
               it != domainElementsEnd; ++it )
          {
              const uint32_t sectors = (*it)->Sectors();
              for( uint32_t sid = 0u; sid < sectors; ++sid )
              {
                  const uint32_t ips = (*it)->FV()->IntegrationPointsPerSector( sid );
                  for( uint32_t ip = 0u; ip < ips; ++ip )
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
          const auto domainElementsEnd( subDomain.CellsEnd() );
          for( auto it = subDomain.CellsBegin();
               it != domainElementsEnd; ++it )
          {
              const uint32_t facets = (*it)->Facets();
              for( uint32_t fid = 0u; fid < facets; ++fid )
              {
                  const uint32_t ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
                  for( uint32_t ip = 0u; ip < ips; ++ip )
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
// ------------------

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
  const auto domainElementsEnd( subDomain.CellsEnd() );
  for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
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
  const auto domainElementsEnd( subDomain.CellsEnd() );
  for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
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
  stringNumber = to_string( tensor_components );
  arrayTitle += stringNumber; arrayTitle += "\" format=\"ascii\"";
  vtu.OpenNode( arrayTitle.c_str() );
  vtu.BringToLevel();
  bool newLine( false );
  TensorVariable<dim> tensorVariable;
  // looping over regions nodes
  const auto domainElementsEnd( subDomain.CellsEnd() );
  for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
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
  for( uint32_t component=0u; component<key.dataDepth; component++ )
  {
      arrayTitle = "DataArray type=\"Float64\" Name=\"";
      variableName = model_.Database().Name( key );
      variableName = FindVariableOutputAlias( variableName );
      stringNumber = to_string( component );
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
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
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
      stringNumber = to_string( component );
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
      const auto domainElementsEnd( subDomain.CellsEnd() );
      for( auto it = subDomain.CellsBegin(); it != domainElementsEnd; ++it )
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

// ============================================================
// findConnectivityFile
// Returns a non-owning raw pointer to the cached XML_Document
// for the given subdomain, or nullptr if not yet cached.
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
XML_Document* VTU_Interface<dim>::findConnectivityFile(
        ConnectivityMap<CELL>&              connectivityMap,
        const ModelSubDomain<dim,CELL>&     subDomain )
{
    const auto it = connectivityMap.find( &subDomain );
    if ( it != connectivityMap.end() )
        return it->second.get();    // non-owning observer; map retains ownership
    return nullptr;
}

template XML_Document* VTU_Interface<1U>::findConnectivityFile<Element>  (VTU_Interface<1U>::ConnectivityMap<Element>&,   const ModelSubDomain<1U,Element>&);
template XML_Document* VTU_Interface<2U>::findConnectivityFile<Element>  (VTU_Interface<2U>::ConnectivityMap<Element>&,   const ModelSubDomain<2U,Element>&);
template XML_Document* VTU_Interface<3U>::findConnectivityFile<Element>  (VTU_Interface<3U>::ConnectivityMap<Element>&,   const ModelSubDomain<3U,Element>&);

template XML_Document* VTU_Interface<1U>::findConnectivityFile<Face>     (VTU_Interface<1U>::ConnectivityMap<Face>&,      const ModelSubDomain<1U,Face>&);
template XML_Document* VTU_Interface<2U>::findConnectivityFile<Face>     (VTU_Interface<2U>::ConnectivityMap<Face>&,      const ModelSubDomain<2U,Face>&);
template XML_Document* VTU_Interface<3U>::findConnectivityFile<Face>     (VTU_Interface<3U>::ConnectivityMap<Face>&,      const ModelSubDomain<3U,Face>&);

template XML_Document* VTU_Interface<1U>::findConnectivityFile<InterFace>(VTU_Interface<1U>::ConnectivityMap<InterFace>&, const ModelSubDomain<1U,InterFace>&);
template XML_Document* VTU_Interface<2U>::findConnectivityFile<InterFace>(VTU_Interface<2U>::ConnectivityMap<InterFace>&, const ModelSubDomain<2U,InterFace>&);
template XML_Document* VTU_Interface<3U>::findConnectivityFile<InterFace>(VTU_Interface<3U>::ConnectivityMap<InterFace>&, const ModelSubDomain<3U,InterFace>&);


// ============================================================
// insertConnectivityFile
// Takes ownership of newConnectivityFile and stores it in the map.
// Returns true if insertion succeeded (map::emplace guarantee).
// ============================================================
template<uint32_t dim>
template<template <uint32_t> class CELL>
bool VTU_Interface<dim>::insertConnectivityFile(
        ConnectivityMap<CELL>&              connectivityMap,
        std::unique_ptr<XML_Document>       newConnectivityFile,
        const ModelSubDomain<dim,CELL>&     subDomain )
{
    const auto result = connectivityMap.emplace(
                            &subDomain, std::move( newConnectivityFile ) );
    return result.second;   // true if inserted, false if key already existed
}

template bool VTU_Interface<1U>::insertConnectivityFile<Element>  (VTU_Interface<1U>::ConnectivityMap<Element>&,   std::unique_ptr<XML_Document>, const ModelSubDomain<1U,Element>&);
template bool VTU_Interface<2U>::insertConnectivityFile<Element>  (VTU_Interface<2U>::ConnectivityMap<Element>&,   std::unique_ptr<XML_Document>, const ModelSubDomain<2U,Element>&);
template bool VTU_Interface<3U>::insertConnectivityFile<Element>  (VTU_Interface<3U>::ConnectivityMap<Element>&,   std::unique_ptr<XML_Document>, const ModelSubDomain<3U,Element>&);

template bool VTU_Interface<1U>::insertConnectivityFile<Face>     (VTU_Interface<1U>::ConnectivityMap<Face>&,      std::unique_ptr<XML_Document>, const ModelSubDomain<1U,Face>&);
template bool VTU_Interface<2U>::insertConnectivityFile<Face>     (VTU_Interface<2U>::ConnectivityMap<Face>&,      std::unique_ptr<XML_Document>, const ModelSubDomain<2U,Face>&);
template bool VTU_Interface<3U>::insertConnectivityFile<Face>     (VTU_Interface<3U>::ConnectivityMap<Face>&,      std::unique_ptr<XML_Document>, const ModelSubDomain<3U,Face>&);

template bool VTU_Interface<1U>::insertConnectivityFile<InterFace>(VTU_Interface<1U>::ConnectivityMap<InterFace>&, std::unique_ptr<XML_Document>, const ModelSubDomain<1U,InterFace>&);
template bool VTU_Interface<2U>::insertConnectivityFile<InterFace>(VTU_Interface<2U>::ConnectivityMap<InterFace>&, std::unique_ptr<XML_Document>, const ModelSubDomain<2U,InterFace>&);
template bool VTU_Interface<3U>::insertConnectivityFile<InterFace>(VTU_Interface<3U>::ConnectivityMap<InterFace>&, std::unique_ptr<XML_Document>, const ModelSubDomain<3U,InterFace>&);




/// Returns a non-owning raw pointer to the cached connectivity document for
// the given subdomain, creating and caching it if it does not yet exist.
// Ownership of the document remains with the connectivity map throughout.
template<uint32_t dim>
template<template <uint32_t> class CELL>
XML_Document* VTU_Interface<dim>::ConnectivityFile(
        ConnectivityMap<CELL>&          connectivityMap,
        const ModelSubDomain<dim,CELL>& subDomain )
{
    // Return cached document if one already exists for this subdomain
    XML_Document* existing = findConnectivityFile( connectivityMap, subDomain );
    if ( existing )
        return existing;

    // No cached document yet — allocate one and write the standard header
    auto newDoc = std::make_unique<XML_Document>();
    EstablishConnectivityFileHeader( *newDoc );

    // Save a non-owning observer before transferring ownership to the map
    XML_Document* observer = newDoc.get();
    insertConnectivityFile( connectivityMap, std::move( newDoc ), subDomain );

    return observer;
}

template XML_Document* VTU_Interface<1U>::ConnectivityFile<Element>  (VTU_Interface<1U>::ConnectivityMap<Element>&,   const ModelSubDomain<1U,Element>&);
template XML_Document* VTU_Interface<2U>::ConnectivityFile<Element>  (VTU_Interface<2U>::ConnectivityMap<Element>&,   const ModelSubDomain<2U,Element>&);
template XML_Document* VTU_Interface<3U>::ConnectivityFile<Element>  (VTU_Interface<3U>::ConnectivityMap<Element>&,   const ModelSubDomain<3U,Element>&);

template XML_Document* VTU_Interface<1U>::ConnectivityFile<Face>     (VTU_Interface<1U>::ConnectivityMap<Face>&,      const ModelSubDomain<1U,Face>&);
template XML_Document* VTU_Interface<2U>::ConnectivityFile<Face>     (VTU_Interface<2U>::ConnectivityMap<Face>&,      const ModelSubDomain<2U,Face>&);
template XML_Document* VTU_Interface<3U>::ConnectivityFile<Face>     (VTU_Interface<3U>::ConnectivityMap<Face>&,      const ModelSubDomain<3U,Face>&);

template XML_Document* VTU_Interface<1U>::ConnectivityFile<InterFace>(VTU_Interface<1U>::ConnectivityMap<InterFace>&, const ModelSubDomain<1U,InterFace>&);
template XML_Document* VTU_Interface<2U>::ConnectivityFile<InterFace>(VTU_Interface<2U>::ConnectivityMap<InterFace>&, const ModelSubDomain<2U,InterFace>&);
template XML_Document* VTU_Interface<3U>::ConnectivityFile<InterFace>(VTU_Interface<3U>::ConnectivityMap<InterFace>&, const ModelSubDomain<3U,InterFace>&);




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
    connectivityFile.AddInfo( "xml version=\"1.0\"" );
    connectivityFile.AddComment( problemTitle_.c_str() );
  }



// TODO: broken for splitboundary output
/// establishes the connectivity file for given region
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::EstablishConnectivityFile( XML_Document& connectivityFile, const ModelSubDomain<dim,CELL>& subDomain ) const
{
  ErrorHandler& csmp_error ( ErrorHandler::Instance() );

  // io feedback
  if ( csmp_error.Verbose() )
      cout << "\nVTU_Interface<dim>::EstablishConnectivityFile: 'Writing connectivity file for sub domain " << DomainName( subDomain ) << " ...";

  // getting node and element count, creating working strings
  string stringNumber, stringCache;

  // opening piece node
  const size_t DOMAIN_NODES( subDomain.RenumberNodes() ), DOMAIN_ELEMENTS( subDomain.Cells() );
  assert ( DOMAIN_NODES > 0 );
  stringNumber = to_string( DOMAIN_NODES );
  stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;
  stringNumber = to_string( DOMAIN_ELEMENTS );
  stringCache += "\" NumberOfCells=\""; stringCache += stringNumber; stringCache += "\"";
  connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

  // points node
  connectivityFile.OpenNode( "Points" );
  connectivityFile.OpenNode( "DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
  // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
  size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
  bool newLine( false );
  connectivityFile.BringToLevel();


  // =====================================================================
  //  POINTS: Writing Node coordinates triples
  // =====================================================================
  // REGION and BOUNDARY domains
  if constexpr ( !std::is_same_v<CELL<dim>,InterFace<dim>> ) {
      const auto domainVerticesEnd( subDomain.NodesEnd() );
      for( auto it = subDomain.NodesBegin(); it != domainVerticesEnd; ++it, ++entriesOfLine )
        {
          // writing x, y and z coordinates(tab seperated)
          stringNumber = numberToString( (*it)->x() );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( yCoordinate( (*it)->Coordinate() ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( zCoordinate( (*it)->Coordinate() ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          // if maximum of entris per line is reached, start new one, else tab
          if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE ) {
              entriesOfLine = 1;
              connectivityFile.LineBreak();
              connectivityFile.BringToLevel();
              newLine = true;
            }
          else {
              connectivityFile.Tab();
              newLine = false;
            }
        } // looping nodes
    }
  // SPLITBOUNDARY domains
  else {
     // finding unit normal to lower-dimensional subdomain
     auto unrml = subDomain.AverageUnitNormal();
     // approximating the dimensions of the SplitBoundary
     Point<dim> pmin, pmax; pmin = 1.0e30; pmax = -1.0e30;
     for ( auto it=subDomain.PerimeterCellsBegin(); it!=subDomain.CellsEnd(); ++it ) {
          pmin = min( pmin, (*it)->BaryCenter() );
          pmax = max( pmax, (*it)->BaryCenter() );
       }
     const double node_offset = distance( pmin, pmax ) * 1.0e-3; // 0.1%
     // inside nodes
     for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).InsideNodes().first ) {
          // writing x, y and z coordinates(tab separated)
          stringNumber = numberToString( nit->x() + unrml[0] * node_offset );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( yCoordinate( nit->Coordinate() + unrml * node_offset ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( zCoordinate( nit->Coordinate() + unrml * node_offset ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          // if maximum of entris per line is reached, start new one, else tab
          if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE ) {
              entriesOfLine = 1;
              connectivityFile.LineBreak();
              connectivityFile.BringToLevel();
              newLine = true;
            }
          else {
              connectivityFile.Tab();
              newLine = false;
            }
          ++entriesOfLine;
       }
     // outside nodes
     for ( const auto& nit : dynamic_cast<const SplitBoundary<dim>&>(subDomain).OutsideNodes().first ) {
          // writing x, y and z coordinates(tab separated)
          stringNumber = numberToString( nit->x() );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( yCoordinate( nit->Coordinate() ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( zCoordinate( nit->Coordinate() ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          // if maximum of entris per line is reached, start new one, else tab
          if( entriesOfLine == MAX_COORDINATE_ENTRIES_PER_LINE ) {
              entriesOfLine = 1;
              connectivityFile.LineBreak();
              connectivityFile.BringToLevel();
              newLine = true;
            }
          else {
              connectivityFile.Tab();
              newLine = false;
            }
          ++entriesOfLine;
       }
    } // end SplitBoundary subdomains
    
  // closing points node
  if( !newLine ) connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  connectivityFile.CloseNode( "Points" );


  // =====================================================================
  //    CELLS: CONNECTIVITY (nodes per element like in plist)
  // =====================================================================
  connectivityFile.OpenNode( "Cells" );
  // cell connectivity
  connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  // establish a vector with vtk element types
  vector<VTK_TYPE> elementTypesVTK;
  vector<size_t>   vtkNodeNumbering;
  const size_t MAX_CONNECTIVITY_ENTRIES_PER_LINE( 20 ); entriesOfLine = 2;
  const auto domainSimplicesEnd( subDomain.CellsEnd() );

  // ------------------------------------
  // Nodes of REGION and BOUNDARY objects
  // ------------------------------------
  if constexpr ( !std::is_same_v<CELL<dim>,InterFace<dim>> )
    {
      for( auto it = subDomain.CellsBegin(); it != domainSimplicesEnd; ++it )
        {
          // CSMP_FINITE_ELEMENT_TYPE
          // storing elements VTK type
          VTK_TYPE elementTypeVTK = ElementType( *it );
          elementTypesVTK.push_back( elementTypeVTK );
          // Quadratic Hexahedron and Quadratic Wedge differ in node numbering(CSMP vs VTK)
          if ( elementTypeVTK == VTK_QUADRATIC_HEXAHEDRON ) {
              vtkNodeNumbering.clear();
              QuadraticHexahedronConnectivity( *it, vtkNodeNumbering );
              assert( vtkNodeNumbering.size() == (*it)->Nodes() );
              for( size_t i{0u}; i < vtkNodeNumbering.size(); ++i, ++entriesOfLine ) {
                   // writing node id to vtu document
                  stringNumber = to_string( vtkNodeNumbering[i] );
                  connectivityFile.InsertData( stringNumber.c_str() );
                  // line break if limeit entries reached, else tab
                  if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
                      entriesOfLine = 1;
                      connectivityFile.LineBreak();
                      connectivityFile.BringToLevel();
                      newLine = true;
                    }
                  else {
                      connectivityFile.Tab();
                      newLine = false;
                    }
                }
              // proceed to next elementactiveProperty += "\"";
              continue;
            }
          if( elementTypeVTK == VTK_QUADRATIC_WEDGE ) {
              vtkNodeNumbering.clear();
              QuadraticWedgeConnectivity( *it, vtkNodeNumbering );
              assert( vtkNodeNumbering.size() == (*it)->Nodes() );
              for( size_t i{0u}; i < vtkNodeNumbering.size(); ++i, ++entriesOfLine ) {
                   // writing node id to vtu document
                  stringNumber = to_string( vtkNodeNumbering[i] );
                  connectivityFile.InsertData( stringNumber.c_str() );
                  // line break if limeit entries reached, else tab
                  if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
                      entriesOfLine = 1;
                      connectivityFile.LineBreak();
                      connectivityFile.BringToLevel();
                      newLine = true;
                    }
                  else {
                      connectivityFile.Tab();
                      newLine = false;
                    }
                }
              // proceed to next element
              continue;
            }
          // for all other element types: looping element's nodes
          for( uint32_t iit = 0u; iit < (*it)->FE()->Nodes(); ++iit, ++entriesOfLine ) {
              // writing node id to vtu document
              stringNumber = to_string( (*it)->N(iit)->Idx() );
              connectivityFile.InsertData( stringNumber.c_str() );
              // line break if limeit entries reached, else tab
              if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
                  entriesOfLine = 1;
                  connectivityFile.LineBreak();
                  connectivityFile.BringToLevel();
                  newLine = true;
                }
              else {
                  connectivityFile.Tab();
                  newLine = false;
                }
            } // end for nodes
        } // end for elements
    } // end REGION and BOUNDARY subdomains
    
  // ---------------------
  // SPLITBOUNDARY objects
  // ---------------------
  else {
      for( const auto& iface : subDomain.CellVector() )
        {
          // CSMP_FINITE_ELEMENT_TYPE
          VTK_TYPE elementTypeVTK = ElementType( iface );
          // translating InterFace objects into element types: line->quad, tria->prism, quad->hex
          switch( elementTypeVTK ) {
               // linear FEM
               case VTK_LINE: elementTypeVTK = VTK_QUAD;
                 break;
               case VTK_TRIANGLE: elementTypeVTK = VTK_WEDGE;
                 break;
               case VTK_QUAD: elementTypeVTK = VTK_HEXAHEDRON;
                 break;
               // quadratic FEM
               case VTK_QUADRATIC_EDGE: elementTypeVTK = VTK_QUADRATIC_QUAD;
                 break;
               case VTK_QUADRATIC_TRIANGLE: elementTypeVTK = VTK_QUADRATIC_WEDGE;
                 break;
               case VTK_QUADRATIC_QUAD: elementTypeVTK = VTK_QUADRATIC_HEXAHEDRON;
                 break;
             default:
               throw csmp::Exception( ERROR, "EstablishConnectivityFile", "InterFace FE-type not recognised");
          }
          elementTypesVTK.push_back( elementTypeVTK );
          // Quadratic Hexahedron and Quadratic Wedge differ in node numbering(CSMP vs VTK)
          // (for InterFace, the inside nodes become the back plane and the outside nodes the front
          if ( elementTypeVTK == VTK_QUADRATIC_HEXAHEDRON ) {
              vtkNodeNumbering.clear();
              QuadraticHexahedronConnectivity( iface, vtkNodeNumbering );
              assert( vtkNodeNumbering.size() == iface->Nodes() );
              for( size_t i{0u}; i < vtkNodeNumbering.size(); ++i, ++entriesOfLine ) {
                   // writing node id to vtu document
                  stringNumber = to_string( vtkNodeNumbering[i] );
                  connectivityFile.InsertData( stringNumber.c_str() );
                  // line break if limeit entries reached, else tab
                  if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
                      entriesOfLine = 1;
                      connectivityFile.LineBreak();
                      connectivityFile.BringToLevel();
                      newLine = true;
                    }
                  else {
                      connectivityFile.Tab();
                      newLine = false;
                    }
                }
              // proceed to next elementactiveProperty += "\"";
              continue;
            }
          if( elementTypeVTK == VTK_QUADRATIC_WEDGE ) {
              vtkNodeNumbering.clear();
              QuadraticWedgeConnectivity( iface, vtkNodeNumbering );
              assert( vtkNodeNumbering.size() == iface->Nodes() );
              for( size_t i{0u}; i < vtkNodeNumbering.size(); ++i, ++entriesOfLine ) {
                   // writing node id to vtu document
                  stringNumber = to_string( vtkNodeNumbering[i] );
                  connectivityFile.InsertData( stringNumber.c_str() );
                  // line break if limeit entries reached, else tab
                  if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
                      entriesOfLine = 1;
                      connectivityFile.LineBreak();
                      connectivityFile.BringToLevel();
                      newLine = true;
                    }
                  else {
                      connectivityFile.Tab();
                      newLine = false;
                    }
                }
              // proceed to next element
              continue;
            }
          // for all other element types: looping element's nodes (InterFace case considered)
          uint32_t node{0};
          for( uint32_t iit = 0u; iit < iface->Nodes(); ++iit, ++entriesOfLine ) {
              // reading the nodes on the inside of the interface first
              if ( iit < iface->FE()->Nodes() )  {
                   iface->CurrentSide( INSIDE );
                   node = iit;
                }
              else {
                   iface->CurrentSide( OUTSIDE );
                   node = iit - iface->FE()->Nodes();
                }
              // writing node id to vtu document
              stringNumber = to_string( iface->N(node)->Idx() );
              connectivityFile.InsertData( stringNumber.c_str() );
              // line break if limeit entries reached, else tab
              if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
                  entriesOfLine = 1;
                  connectivityFile.LineBreak();
                  connectivityFile.BringToLevel();
                  newLine = true;
                }
              else {
                  connectivityFile.Tab();
                  newLine = false;
                }
            } // end nodes
        } // end InterFace objects
  
   } // end SPLITBOUNDARY
    
  // close cell connectivity
  if( !newLine ) connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  
  
  // =====================================================================
  //    CELLS: OFFSETS (nodes per element matching type specs below)
  // =====================================================================
  // (works for all types of subdomains)
  connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  // looping subdomains's elements
  size_t offset( 0 ); entriesOfLine = 2;
  for( auto it = subDomain.CellsBegin(); it != domainSimplicesEnd; ++it, ++entriesOfLine ) {
      // imcrementing offset by node count and writing to data
      offset += (*it)->Nodes(); // NB: this gives correct number for INTERFACE
      stringNumber = to_string( offset );
      connectivityFile.InsertData( stringNumber.c_str() );
      // line break if limeit entries reached, else tab
      if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
          entriesOfLine = 1;
          connectivityFile.LineBreak();
          connectivityFile.BringToLevel();
          newLine = true;
        }
      else {
          connectivityFile.Tab();
          newLine = false;
        }
    } // region's elements
  // close offsets
  if( !newLine ) connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  

  // =====================================================================
  //    CELL_TYPES: VTK_ FE types
  // =====================================================================
  // (InterFace objects are treated as Quads, Prisms or Hexahedra, see above)
  assert( elementTypesVTK.size() == subDomain.Cells() );
  connectivityFile.OpenNode( "DataArray type=\"UInt8\" Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  // looping region's elements
  entriesOfLine = 2;
  size_t elementCount{0ul};
  for( auto it = subDomain.CellsBegin(); it != domainSimplicesEnd; ++it, ++entriesOfLine, ++elementCount )
    {
      // imcrementing offset by node count and writing to data
      stringNumber = to_string( elementTypesVTK[elementCount] );
      connectivityFile.InsertData( stringNumber.c_str() );
      // line break if limeit entries reached, else tab
      if( entriesOfLine == MAX_CONNECTIVITY_ENTRIES_PER_LINE ) {
          entriesOfLine = 1;
          connectivityFile.LineBreak();
          connectivityFile.BringToLevel();
          newLine = true;
        }
      else {
          connectivityFile.Tab();
          newLine = false;
        }
    } // subdomain elements
    
  // close offsets
  if( !newLine ) connectivityFile.LineBreak();
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
  const size_t DOMAIN_ELEMENTS( subDomain.Cells() );
  stringNumber = to_string( DOMAIN_ELEMENTS );
  stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

  // no elements
  stringCache += "\" NumberOfCells=\"0\"";
  connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

  // points node
  connectivityFile.OpenNode( "Points" );
  connectivityFile.OpenNode( "DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
  // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
  size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
  bool newLine( false );
  connectivityFile.BringToLevel();
  // looping over all the region's nodes
  const auto domainSimpicesEnd( subDomain.CellsEnd() );
  for( auto it = subDomain.CellsBegin(); it != domainSimpicesEnd; ++it, ++entriesOfLine )
    {
    // writing x,y and z coordinates(tab seperated)
    stringNumber = numberToString( (*it)->BaryCenter()[0] );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = numberToString( yCoordinate( (*it)->BaryCenter() ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = numberToString( zCoordinate( (*it)->BaryCenter() ) );
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
    stringNumber = to_string( pid );
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
  stringNumber = to_string( DOMAIN_ELEMENTS );
  connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
  connectivityFile.BringToLevel();
  connectivityFile.InsertData( stringNumber.c_str() );
  connectivityFile.LineBreak();
  connectivityFile.CloseNode( "DataArray" );
  // cell types
  size_t cell_type( VTK_POLY_VERTEX );
  stringNumber = to_string( cell_type );
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
    stringNumber = to_string( REGION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
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
    const auto domainSimplicesEnd( subDomain.CellsEnd() );
    for( auto it  = subDomain.CellsBegin(); it != domainSimplicesEnd; ++it )
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
    stringNumber = numberToString( pt[0] );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = numberToString( yCoordinate( pt ) );
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.Tab();
    stringNumber = numberToString( zCoordinate( pt ) );
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
    stringNumber = to_string( pid );
    connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\"");
    connectivityFile.BringToLevel();
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    // cell offsets
    size_t offset(1);
    stringNumber = to_string( offset );
    connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
    connectivityFile.BringToLevel();
    connectivityFile.InsertData( stringNumber.c_str() );
    connectivityFile.LineBreak();
    connectivityFile.CloseNode( "DataArray" );
    // cell types
    size_t cell_type( VTK_POLY_VERTEX ); // VTK_VERTEX would be also fine, but this way it's more general
    stringNumber = to_string( cell_type );
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
    stringNumber = to_string( INTEGRATION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
    // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
    size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
    bool newLine( false );
    connectivityFile.BringToLevel();

    // looping over all element integration points
    const auto domainSimplicesEnd( subDomain.CellsEnd() );
    for( auto it  = subDomain.CellsBegin(); it != domainSimplicesEnd; ++it )
    {
        const uint32_t ips = (*it)->FE()->IntegrationPoints();
        for( uint32_t ip{0u}; ip < ips; ++ip, ++entriesOfLine )
        {
          auto pt = (*it)->IntegrationPoint( ip );
          // writing x,y and z coordinates(tab seperated)
          stringNumber = numberToString( pt[0] );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( yCoordinate( pt ) );
          connectivityFile.InsertData( stringNumber.c_str() );
          connectivityFile.Tab();
          stringNumber = numberToString( zCoordinate( pt ) );
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
        stringNumber = to_string( pid );
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
      stringNumber = to_string( INTEGRATION_POINTS );
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell types
      size_t cell_type( VTK_POLY_VERTEX );
      stringNumber = to_string( cell_type );
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
    const size_t INTEGRATION_POINTS{ subDomain.SectorIntegrationPoints() };
    stringNumber = to_string( INTEGRATION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
    // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
    size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
    bool newLine( false );
    connectivityFile.BringToLevel();

    // looping over all the region's nodes
    const auto domainSimplicesEnd( subDomain.CellsEnd() );
    for( auto it = subDomain.CellsBegin(); it != domainSimplicesEnd; ++it )
    {
        const uint32_t sectors = (*it)->Sectors();
        for( uint32_t sid = 0u; sid < sectors; ++sid )
        {
            const uint32_t ips = (*it)->FV()->IntegrationPointsPerSector( sid );
            for( uint32_t ip = 0u; ip < ips; ++ip, ++entriesOfLine )
            {
              auto pt = (*it)->FV()->SectorIntegrationPoint( sid, ip );
              pt = (*it)->RstToXYZ( pt );
              // writing x,y and z coordinates(tab seperated)
              stringNumber = numberToString( pt[0] );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = numberToString( yCoordinate( pt ) );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = numberToString( zCoordinate( pt ) );
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
      for( size_t pid = 0u; pid < INTEGRATION_POINTS; ++pid  )
      {
        // writing node id to vtu document
        stringNumber = to_string( pid );
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
      stringNumber = to_string( INTEGRATION_POINTS );
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell types
      const size_t cell_type( VTK_POLY_VERTEX );
      stringNumber = to_string( cell_type );
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
    stringNumber = to_string( INTEGRATION_POINTS );
    stringCache = "Piece NumberOfPoints=\""; stringCache += stringNumber;

    // no elements
    stringCache += "\" NumberOfCells=\"1\"";
    connectivityFile.OpenNode( stringCache.c_str() ); stringCache.clear(); stringNumber.clear();

    // points node
    connectivityFile.OpenNode( "Points" );
    connectivityFile.OpenNode( "DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"3\" format=\"ascii\"");
    // running count of data entries per line and given maximum(new line beyond that) starting at current indentation
    size_t entriesOfLine( 2 ); const size_t MAX_COORDINATE_ENTRIES_PER_LINE( 5 );
    bool newLine( false );
    connectivityFile.BringToLevel();

    // looping over all the region's nodes
    const auto domainSimplicesEnd( subDomain.CellsEnd() );
    for( auto it = subDomain.CellsBegin(); it != domainSimplicesEnd; ++it )
    {
        const uint32_t facets = (*it)->Facets();
        for( uint32_t fid = 0u; fid < facets; ++fid )
        {
            const uint32_t ips = (*it)->FV()->IntegrationPointsPerFacet( fid );
            for( uint32_t ip = 0u; ip < ips; ++ip, ++entriesOfLine )
            {
              auto pt = (*it)->FV()->FacetIntegrationPoint( fid, ip );
              pt = (*it)->RstToXYZ( pt );
              // writing x,y and z coordinates(tab seperated)
              stringNumber = numberToString( pt[0] );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = numberToString( yCoordinate( pt ) );
              connectivityFile.InsertData( stringNumber.c_str() );
              connectivityFile.Tab();
              stringNumber = numberToString( zCoordinate( pt ) );
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
      for( size_t pid = 0u; pid < INTEGRATION_POINTS; ++pid  )
      {
        // writing node id to vtu document
        stringNumber = to_string( pid );
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
      stringNumber = to_string( INTEGRATION_POINTS );
      connectivityFile.OpenNode( "DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\"");
      connectivityFile.BringToLevel();
      connectivityFile.InsertData( stringNumber.c_str() );
      connectivityFile.LineBreak();
      connectivityFile.CloseNode( "DataArray" );
      // cell types
      const int cell_type( VTK_POLY_VERTEX );
      stringNumber = to_string( cell_type );
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
void VTU_Interface<dim>::QuadraticHexahedronConnectivity( const CELL<dim>* const element, vector<size_t>& data ) const
{
  if ( !data.empty() ) data.clear();
  data.reserve(20);
  // first 12 nodes are the same
  for( uint32_t i{0u}; i < 12; ++i )
    data.push_back( element->N( i )->Idx() );
  // now the node convention differs
  data.push_back( element->N( 16 )->Idx() );
  data.push_back( element->N( 17 )->Idx() );
  data.push_back( element->N( 18 )->Idx() );
  data.push_back( element->N( 19 )->Idx() );
  data.push_back( element->N( 12 )->Idx() );
  data.push_back( element->N( 13 )->Idx() );
  data.push_back( element->N( 14 )->Idx() );
  data.push_back( element->N( 15 )->Idx() );
}

template void VTU_Interface<1U>::QuadraticHexahedronConnectivity(const Element<1U>* const,vector<size_t>&) const;
template void VTU_Interface<2U>::QuadraticHexahedronConnectivity(const Element<2U>* const,vector<size_t>&) const;
template void VTU_Interface<3U>::QuadraticHexahedronConnectivity(const Element<3U>* const,vector<size_t>&) const;

template void VTU_Interface<1U>::QuadraticHexahedronConnectivity(const Face<1U>* const,vector<size_t>&) const;
template void VTU_Interface<2U>::QuadraticHexahedronConnectivity(const Face<2U>* const,vector<size_t>&) const;
template void VTU_Interface<3U>::QuadraticHexahedronConnectivity(const Face<3U>* const,vector<size_t>&) const;

template void VTU_Interface<1U>::QuadraticHexahedronConnectivity(const InterFace<1U>* const,vector<size_t>&) const;
template void VTU_Interface<2U>::QuadraticHexahedronConnectivity(const InterFace<2U>* const,vector<size_t>&) const;
template void VTU_Interface<3U>::QuadraticHexahedronConnectivity(const InterFace<3U>* const,vector<size_t>&) const;




/// converts csmp to vtk node numbering
template<uint32_t dim>
template<template <uint32_t> class CELL>
void VTU_Interface<dim>::QuadraticWedgeConnectivity( const CELL<dim>* const element, vector<size_t>& data ) const
{
  data.reserve(15);
  // first 8 nodes are the same
  for( uint32_t i{0}; i < 9; ++i )
    data.push_back( element->N( i )->Idx() );
  // now the node convention differs
  data.push_back( element->N( static_cast<uint32_t>(12) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(13) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(14) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(9) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(10) )->Idx() );
  data.push_back( element->N( static_cast<uint32_t>(11) )->Idx() );
}

template void VTU_Interface<1U>::QuadraticWedgeConnectivity(const Element<1U>* const,vector<size_t>&) const;
template void VTU_Interface<2U>::QuadraticWedgeConnectivity(const Element<2U>* const,vector<size_t>&) const;
template void VTU_Interface<3U>::QuadraticWedgeConnectivity(const Element<3U>* const,vector<size_t>&) const;

template void VTU_Interface<1U>::QuadraticWedgeConnectivity(const Face<1U>* const,vector<size_t>&) const;
template void VTU_Interface<2U>::QuadraticWedgeConnectivity(const Face<2U>* const,vector<size_t>&) const;
template void VTU_Interface<3U>::QuadraticWedgeConnectivity(const Face<3U>* const,vector<size_t>&) const;

template void VTU_Interface<1U>::QuadraticWedgeConnectivity(const InterFace<1U>* const,vector<size_t>&) const;
template void VTU_Interface<2U>::QuadraticWedgeConnectivity(const InterFace<2U>* const,vector<size_t>&) const;
template void VTU_Interface<3U>::QuadraticWedgeConnectivity(const InterFace<3U>* const,vector<size_t>&) const;





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
   return subDomain.Name() + "_(REGION)";
 }

template<uint32_t dim>
string VTU_Interface<dim>::DomainName( const ModelSubDomain<dim,Face>& subDomain ) const
{
   return subDomain.Name() + "_(BOUNDARY)";
 }

template<uint32_t dim>
string VTU_Interface<dim>::DomainName( const ModelSubDomain<dim,InterFace>& subDomain ) const
{
   return subDomain.Name() + "_(SPLITBOUNDARY)";
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
