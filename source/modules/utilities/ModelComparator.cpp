// SPDX-FileCopyrightText: © 2026 The openCSMP project
//
// SPDX-License-Identifier: LGPL-3.0-only

#include "ModelComparator.h"

#include "Model.h"
#include "VSet.h"
#include "Exception.h"

#include <cmath>

using namespace std;

namespace csmp
{

template<uint32_t dim>
double ModelComparator<dim>::CompareModels( const Model<dim>& model1, const Model<dim>& model2,
                                            const char* property1, const char* property2,
                                            const char* region ) const
{
    double error( 0. );

    const Index model1Key( model1.Database().StorageKey( property1 ) );
    const Index model2Key( model2.Database().StorageKey( property2 ) );

    if( model1Key.type != model2Key.type || model1Key.place != model2Key.place )
        throw Exception( ERROR,
                         "ModelComparator::CompareModels",
                         "Variables seem to be of different nature" );

    if( model1Key.type == SCALAR )
        error = CompareRegionScalarVariable( model1, model2, model1Key, model2Key, region );
    else
        throw Exception( ERROR,
                         "ModelComparator::CompareModels",
                         "Comparison for Variable type not implemented yet." );

    return error;
}

template<uint32_t dim>
double ModelComparator<dim>::CompareVSets( VSet<dim>& vset1, VSet<dim>& vset2,
                                           const char* property1, const char* property2,
                                           const char* variablesFile1, const char* variablesFile2,
                                           const char* region ) const
{
    double error( 0. );

    Model<dim> model1( vset1, variablesFile1 );
    Model<dim> model2( vset2, variablesFile2 );

    error = CompareModels( model1, model2, property1, property2, region );

    return error;
}


template<uint32_t dim>
double ModelComparator<dim>::CompareVSets( const char* vset1File, const char* vset2File,
                                           const char* property1, const char* property2,
                                           const char* variablesFile1, const char* variablesFile2,
                                           const char* region ) const
{
    double error( 0. );

    VSet<dim> vset1, vset2;
    double modelTime( 0. );
    vset1.InputFrom( vset1File, modelTime );
    vset2.InputFrom( vset2File, modelTime );

    error = CompareVSets( vset1, vset2, property1, property2, variablesFile1, variablesFile2, region );

    return error;
}




template<uint32_t dim>
double ModelComparator<dim>::CompareRegionScalarVariable( const Model<dim>& model1, const Model<dim>& model2,
                                                          const Index model1Key, const Index model2Key,
                                                          const char* region ) const
{
    deque<ScalarVariable> model1Values;
    deque<ScalarVariable> model2Values;

    const Region<dim>& rref1( model1.Region( region ) );
    const Region<dim>& rref2( model2.Region( region ) );

    if( rref1.Nodes() != rref2.Nodes() )
        throw Exception( ERROR,
                         "ModelComparator::CompareScalarNodalVariable",
                         "Regions seem to be of different size." );

    if( model1Key.place == NODE )
    {
        ReadRegionNodalScalarVariable( rref1, model1Key, model1Values );
        ReadRegionNodalScalarVariable( rref2, model2Key, model2Values );
    }
    else if( model1Key.place == ELEMENT )
    {
        ReadRegionElementScalarVariable( rref1, model1Key, model1Values );
        ReadRegionElementScalarVariable( rref2, model2Key, model2Values );
    }
    else
    {
        throw Exception( ERROR,
                         "ModelComparator::CompareScalarNodalVariable",
                         "Comparison for this variable placement not implemented." );
    }

    double error( CompareScalarDequesL2( model1Values, model2Values ) );

    return error;
}

/// L2 norm comparison of two deques containing scalar variables
template<uint32_t dim>
double ModelComparator<dim>::CompareScalarDequesL2( const std::deque<ScalarVariable>& deque1,
                                                    const std::deque<ScalarVariable>& deque2 ) const
{
    double error( 0. );

    if( deque1.size() != deque2.size() )
        throw Exception( ERROR,
                         "ModelComparator::CompareScalarDeques",
                         "Deques unequal in size." );


    for( size_t i = 0; i < deque1.size(); ++i ) error += std::pow( std::fabs( deque1[i]() - deque2[i]() ), 2 );
 
    return std::sqrt( error );
}


template<uint32_t dim>
void ModelComparator<dim>::ReadRegionElementScalarVariable( const Region<dim>& region,
                                                            const Index propKey,
                                                            std::deque<ScalarVariable>& scalarDeque ) const
{
    ScalarVariable scalarValue( PLAIN, 0. );
    const auto elementsEnd( region.CellsEnd() );
    for( auto it = region.CellsBegin(); it != elementsEnd; ++it )
    {
        (*it)->Read( propKey, scalarValue );
        scalarDeque.push_back( scalarValue );
    } // elements of region



}

template<uint32_t dim>
void ModelComparator<dim>::ReadRegionNodalScalarVariable( const Region<dim>& region,
                                                          const Index propKey,
                                                          std::deque<ScalarVariable>& scalarDeque ) const
{
    ScalarVariable scalarValue( PLAIN, 0. );
    const auto nodesEnd( region.NodesEnd() );
    for( auto it = region.NodesBegin(); it != nodesEnd; ++it )
    {
        (*it)->Read( propKey, scalarValue );
        scalarDeque.push_back( scalarValue );
    } // nodes of region



}





///////////////////////////////////////////////////////////////////////////////////////////////
/// Comparison is done for identically equal meshes but might be different numbering of nodes//
///////////////////////////////////////////////////////////////////////////////////////////////



template<uint32_t dim>
double ModelComparator<dim>::CompareModelsRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                                           const char* property1, const char* property2,
                                                           const char* region )
{
    double error( 0. );

    const Index model1Key( model1.Database().StorageKey( property1 ) );
    const Index model2Key( model2.Database().StorageKey( property2 ) );

    if( model1Key.type != model2Key.type or model1Key.place != model2Key.place)
        throw Exception( ERROR,
                         "ModelComparator::CompareModelsAtPossiblyDifferentsMeshes",
                         "Variables seem to be of different nature" );

    if( model1Key.type == SCALAR )
        error = CompareRegionScalarVariableRenumberedNodes( model1, model2, property1, property2, region );
    else
        throw Exception( ERROR,
                         "ModelComparator::CompareModelsRenumberedNodes",
                         "Comparison for Variable type not implemented yet." );

    return error;
}

template<uint32_t dim>
double ModelComparator<dim>::CompareVSetsRenumberedNodes( VSet<dim>& vset1, VSet<dim>& vset2,
                                                          const char* property1, const char* property2,
                                                          const char* variablesFile1, const char* variablesFile2,
                                                          const char* region )
{
    double error( 0. );

    Model<dim> model1( vset1, variablesFile1 );
    Model<dim> model2( vset2, variablesFile2 );

    error = CompareModelsRenumberedNodes( model1, model2, property1, property2, region );

    return error;
}


template<uint32_t dim>
double ModelComparator<dim>::CompareVSetsRenumberedNodes( const char* vset1File, const char* vset2File,
                                                          const char* property1, const char* property2,
                                                          const char* variablesFile1, const char* variablesFile2,
                                                          const char* region )
{
    double error( 0. );

    VSet<dim> vset1, vset2;
    double modelTime( 0. );
    vset1.InputFrom( vset1File, modelTime );
    vset2.InputFrom( vset2File, modelTime );

    error = CompareVSetsRenumberedNodes( vset1, vset2, property1, property2, variablesFile1, variablesFile2, region );

    return error;
}





template<uint32_t dim>
double ModelComparator<dim>::CompareRegionScalarVariableRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                                                         const char* property1, const char* property2,
                                                                         const char* region )
{
    map<Point<dim>, ScalarVariable > points_and_values1, points_and_values2;

    const Region<dim>& rref1( model1.Region( region ) );
    const Region<dim>& rref2( model2.Region( region ) );

    if( rref1.Nodes() != rref2.Nodes() )
        throw Exception( ERROR,
                         "ModelComparator::CompareRegionScalarVariableRenumberedNodes",
                         "Regions seem to be of different size." );

    const Index model1Key( model1.Database().StorageKey( property1 ) );
    const Index model2Key( model1.Database().StorageKey( property2 ) );

    if( model1Key.place == NODE )
    {
        ReadRegionNodalScalarVariableAndNodeCoordinates( rref1, model1Key, points_and_values1 );
        ReadRegionNodalScalarVariableAndNodeCoordinates( rref2, model2Key, points_and_values2 );

    }
    else
    {
        throw Exception( ERROR,
                         "ModelComparator::CompareRegionScalarVariableRenumberedNodes",
                         "Comparison for this variable placement not implemented." );
    }

    double error( 0. );
    ScalarVariable scalarValue( PLAIN, 0. );
    // L2 norm comparison

    typename std::map<Point<dim>, ScalarVariable >::iterator pvit1,pvit2,pvit_find;

    for ( pvit1=points_and_values1.begin(), pvit2 =points_and_values2.begin();pvit1!=points_and_values1.end(); pvit1++,pvit2++ )
    {

        if( (*pvit1).first!=(*pvit2).first){

            pvit_find = points_and_values2.find((*pvit1).first);

            if(pvit_find==points_and_values2.end()){
                throw Exception( ERROR,
                             "ModelComparator::CompareRegionScalarVariableRenumberedNodes",
                             "The node cannot be found. Probably the meshes are different" );
            }else{
                error += ( (*pvit1).second() - (*pvit_find).second() )*( (*pvit1).second() - (*pvit_find).second() );
            }
        }else
            error += ( (*pvit1).second() - (*pvit2).second() )*( (*pvit1).second() - (*pvit2).second() );

    }


    return std::sqrt( error );

}

template<uint32_t dim>
void ModelComparator<dim>::ReadRegionNodalScalarVariableAndNodeCoordinates(const Region<dim>& region,
                                                          const Index propKey,
                                                          map<Point<dim>, ScalarVariable >& points_and_values)
{
    ScalarVariable scalarValue( PLAIN, 0. );
    const auto nodesEnd( region.NodesEnd() );
    for( auto it = region.NodesBegin(); it != nodesEnd; ++it )
      {
          (*it)->Read( propKey, scalarValue );

          points_and_values[(*it)->Coordinate()] = scalarValue ;

      } // nodes of region

}

/////////////////////////////////////////////////////////////////////////////////////////
/// Comparison is done for meshes that can be in general different from each other     //
/////////////////////////////////////////////////////////////////////////////////////////

template<uint32_t dim>
double ModelComparator<dim>::CompareModelsAtPoints( Model<dim>& model1, Model<dim>& model2,
                                              const char* property1, const char* property2, const char* region )
{
    double error( 0. );

    Index model1Key( model1.Database().StorageKey( property1 ) );
    Index model2Key( model2.Database().StorageKey( property2 ) );

    if( model1Key.type != model2Key.type || model1Key.place != model2Key.place)
        throw Exception( ERROR,
                         "ModelComparator::CompareModelsAtPoints",
                         "Variables seem to be of different nature" );

    if( model1Key.type == SCALAR )
        error = CompareRegionScalarVariableAtPoints( model1, model2, property1, property2, region );
    else
    {
        cout<<"\nERROR:Variable has not valid for comparison type:"
            <<"\nFirst type: " <<parseType( model1Key.type )
            <<"\nSecond type: "<<parseType( model2Key.type )<<"\n";
        throw Exception( ERROR,
                         "ModelComparator::CompareModelsAtPossiblyDifferentsMeshes",
                         "Comparison for Variable type not implemented yet." );
    }

    return error;
}

template<uint32_t dim>
double ModelComparator<dim>::CompareVSetsAtPoints( VSet<dim>& vset1, VSet<dim>& vset2,
                                                   const char* property1, const char* property2,
                                                   const char* variablesFile1, const char* variablesFile2, const char* region )
{
    double error( 0. );

    Model<dim> model1( vset1, variablesFile1 );
    Model<dim> model2( vset2, variablesFile2 );

    error = CompareModelsAtPoints( model1, model2, property1, property2, region );

    return error;
}


template<uint32_t dim>
double ModelComparator<dim>::CompareVSetsAtPoints( const char* vset1File, const char* vset2File,
                                                   const char* property1, const char* property2,
                                                   const char* variablesFile1, const char* variablesFile2, const char* region )
{
    double error( 0. );

    VSet<dim> vset1, vset2;
    double modelTime( 0. );
    vset1.InputFrom( vset1File, modelTime );
    vset2.InputFrom( vset2File, modelTime );

    error = CompareVSetsAtPoints( vset1, vset2, property1, property2, variablesFile1, variablesFile2, region );

    return error;
}





template<uint32_t dim>
double ModelComparator<dim>::CompareRegionScalarVariableAtPoints( Model<dim>& model1, Model<dim>& model2,
                                                                  const char* property1, const char* property2, const char* region )
{
    std::vector<ScalarVariable> values;
    std::map<size_t, std::vector<double> > points;

    const Region<dim>& rref1( model1.Region( region ) );
    const Region<dim>& rref2( model2.Region( region ) );

    if( rref1.Nodes() != rref2.Nodes() )
        throw Exception( ERROR,
                         "ModelComparator::CompareRegionScalarVariableAtPoints",
                         "Regions seem to be of different size." );

    Index model1Key( model1.Database().StorageKey( property1 ) );
    Index model2Key( model1.Database().StorageKey( property2 ) );

    if( model1Key.place == NODE )
    {
        ReadRegionNodalScalarVariableAndNodeCoordinates( rref1, model1Key, values, points );
    }
    else
    {
        cout<<"\nERROR:Variable has not valid for comparison placement:"
            <<"\nFirst place: " <<parsePlacement( model1Key.place )
            <<"\nSecond place: "<<parsePlacement( model2Key.place )<<"\n";
        throw Exception( ERROR,
                         "ModelComparator::CompareRegionScalarVariableAtPoints",
                         "Comparison for this variable placement not implemented." );
    }

    double error( 0. );

    // finds values of certain property at given points
    PropertyAtPointVisitor<dim> pAt(model2, points, property2);
    pAt.SetBruteForceOff();
    model2.Accept(pAt);
    pAt.CheckResults();

    // L2 norm comparison
    size_t numberOfPoints( points.size() );
    ScalarVariable scalarValue( PLAIN, 0. );
    for (size_t i{0U}; i< numberOfPoints; i++ )
    {
        pAt.PropertyValueAt(i, scalarValue);
        error += (scalarValue() - values[i]())*(scalarValue() - values[i]());
    }


    return std::sqrt( error );

}


template<uint32_t dim>
void ModelComparator<dim>::ReadRegionNodalScalarVariableAndNodeCoordinates( const Region<dim>& region,
                                                                            const Index propKey,
                                                                            std::vector<ScalarVariable> & values,
                                                                            map<size_t, std::vector<double> >& points)
{
    size_t index(0U);
    ScalarVariable scalarValue( PLAIN, 0. );
    const auto nodesEnd( region.NodesEnd() );

    values.clear();
    values.resize(region.Nodes());
    for( auto it = region.NodesBegin(); it != nodesEnd; ++it )
    {
        (*it)->Read( propKey, scalarValue );
        for(auto i{0U};i<dim;i++)
            points[index].push_back((*it)->operator[](i));

        values[index] = scalarValue ;
        index++;

    } // nodes of region

}


// explicits
template class ModelComparator<3U>;
template class ModelComparator<2U>;

} // csmp
