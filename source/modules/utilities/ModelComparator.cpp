#include "ModelComparator.h"

#include "Model.h"
#include "VSet.h"
#include "Exception.h"

#include <cmath>

using namespace std;

namespace csmp
{

/////////////////////////////////////////////////////////////////////////////////////////
/// Comparison is done for identically equal meshes with identical numberings of nodes//
/////////////////////////////////////////////////////////////////////////////////////////



template<size_t dim>
double64 ModelComparator<dim>::CompareModels( const Model<dim>& model1, const Model<dim>& model2,
                                              const char* property1, const char* property2, const char* region ) const
{
    double64 error( 0. );

    const Index model1Key( model1.Database().StorageKey( property1 ) );
    const Index model2Key( model2.Database().StorageKey( property2 ) );

    if( model1Key.type != model2Key.type || model1Key.place != model2Key.place)
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareModels",
                         "Variables seem to be of different nature" );

    if( model1Key.type == SCALAR )
        error = CompareRegionScalarVariable( model1, model2, model1Key, model2Key, region );
    else
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareModels",
                         "Comparison for Variable type not implemented yet." );

    return error;
}

template<size_t dim>
double64 ModelComparator<dim>::CompareVSets( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                                             const char* variablesFile1, const char* variablesFile2,
                                             bool isoparametric, const char* region ) const
{
    double64 error( 0. );

    Model<dim> model1( vset1, variablesFile1, isoparametric );
    Model<dim> model2( vset2, variablesFile2, isoparametric );

    error = CompareModels( model1, model2, property1, property2, region );

    return error;
}


template<size_t dim>
double64 ModelComparator<dim>::CompareVSets( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                                             const char* variablesFile1, const char* variablesFile2,
                                             bool isoparametric, const char* region ) const
{
    double64 error( 0. );

    VSet<dim> vset1, vset2;
    double64 modelTime( 0. );
    vset1.InputFrom( vset1File, modelTime );
    vset2.InputFrom( vset2File, modelTime );

    error = CompareVSets( vset1, vset2, property1, property2, variablesFile1, variablesFile2, isoparametric, region );

    return error;
}


template<size_t dim>
double64 ModelComparator<dim>::CompareRegionScalarVariable( const Model<dim>& model1, const Model<dim>& model2,
                                                            const Index model1Key, const Index model2Key, const char* region ) const
{
    deque<ScalarVariable> model1Values;
    deque<ScalarVariable> model2Values;

    const Region<dim>& rref1( model1.Region( region ) );
    const Region<dim>& rref2( model2.Region( region ) );

    if( rref1.Nodes() != rref2.Nodes() )
        throw Exception( CSMP_ERROR,
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
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareScalarNodalVariable",
                         "Comparison for this variable placement not implemented." );
    }

    double64 error( CompareScalarDequesL2( model1Values, model2Values ) );

    return error;
}

/// L2 norm comparison of two deques containing scalar variables
template<size_t dim>
double64 ModelComparator<dim>::CompareScalarDequesL2( const std::deque<ScalarVariable>& deque1,
                                                      const std::deque<ScalarVariable>& deque2 ) const
{
    double64 error( 0. );

    if( deque1.size() != deque2.size() )
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareScalarDeques",
                         "Deques unequal in size." );


    for( size_t i = 0; i < deque1.size(); ++i )
    {
        error += std::pow( std::fabs( deque1[i].Value() - deque2[i].Value() ), 2 );
    }

    return std::sqrt( error );
}


template<size_t dim>
void ModelComparator<dim>::ReadRegionElementScalarVariable( const Region<dim>& region,
                                                            const Index propKey,
                                                            std::deque<ScalarVariable>& scalarDeque ) const
{
    ScalarVariable scalarValue( PLAIN, 0. );
    const typename std::vector<Element<dim>*>::const_iterator elementsEnd( region.ElementsEnd() );
    for( typename std::vector<Element<dim>*>::const_iterator it = region.ElementsBegin(); it != elementsEnd; ++it )
    {
        (*it)->Read( propKey, scalarValue );
        scalarDeque.push_back( scalarValue );
    } // elements of region



}

template<size_t dim>
void ModelComparator<dim>::ReadRegionNodalScalarVariable( const Region<dim>& region,
                                                          const Index propKey,
                                                          std::deque<ScalarVariable>& scalarDeque ) const
{
    ScalarVariable scalarValue( PLAIN, 0. );
    const typename std::vector<Node<dim>* >::const_iterator nodesEnd( region.NodesEnd() );
    for( typename std::vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != nodesEnd; ++it )
    {
        (*it)->Read( propKey, scalarValue );
        scalarDeque.push_back( scalarValue );
    } // nodes of region



}





///////////////////////////////////////////////////////////////////////////////////////////////
/// Comparison is done for identically equal meshes but might be different numbering of nodes//
///////////////////////////////////////////////////////////////////////////////////////////////



template<size_t dim>
double64 ModelComparator<dim>::CompareModelsRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                              const char* property1, const char* property2, const char* region )
{
    double64 error( 0. );

    const Index model1Key( model1.Database().StorageKey( property1 ) );
    const Index model2Key( model2.Database().StorageKey( property2 ) );

    if( model1Key.type != model2Key.type or model1Key.place != model2Key.place)
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareModelsAtPossiblyDifferentsMeshes",
                         "Variables seem to be of different nature" );

    if( model1Key.type == SCALAR )
        error = CompareRegionScalarVariableRenumberedNodes( model1, model2, property1, property2, region );
    else
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareModelsRenumberedNodes",
                         "Comparison for Variable type not implemented yet." );

    return error;
}

template<size_t dim>
double64 ModelComparator<dim>::CompareVSetsRenumberedNodes( VSet<dim>& vset1, VSet<dim>& vset2, const char* property1, const char* property2,
                                             const char* variablesFile1, const char* variablesFile2,
                                             bool isoparametric, const char* region )
{
    double64 error( 0. );

    Model<dim> model1( vset1, variablesFile1, isoparametric );
    Model<dim> model2( vset2, variablesFile2, isoparametric );

    error = CompareModelsRenumberedNodes( model1, model2, property1, property2, region );

    return error;
}


template<size_t dim>
double64 ModelComparator<dim>::CompareVSetsRenumberedNodes( const char* vset1File, const char* vset2File, const char* property1, const char* property2,
                                             const char* variablesFile1, const char* variablesFile2,
                                             bool isoparametric, const char* region )
{
    double64 error( 0. );

    VSet<dim> vset1, vset2;
    double64 modelTime( 0. );
    vset1.InputFrom( vset1File, modelTime );
    vset2.InputFrom( vset2File, modelTime );

    error = CompareVSetsRenumberedNodes( vset1, vset2, property1, property2, variablesFile1, variablesFile2, isoparametric, region );

    return error;
}





template<size_t dim>
double64 ModelComparator<dim>::CompareRegionScalarVariableRenumberedNodes( Model<dim>& model1, Model<dim>& model2,
                                                                    const char* property1, const char* property2, const char* region )
{
    map<Point<dim>, ScalarVariable > points_and_values1, points_and_values2;

    const Region<dim>& rref1( model1.Region( region ) );
    const Region<dim>& rref2( model2.Region( region ) );

    if( rref1.Nodes() != rref2.Nodes() )
        throw Exception( CSMP_ERROR,
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
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareRegionScalarVariableRenumberedNodes",
                         "Comparison for this variable placement not implemented." );
    }

    double64 error( 0. );
    ScalarVariable scalarValue( PLAIN, 0. );
    // L2 norm comparison

    typename std::map<Point<dim>, ScalarVariable >::iterator pvit1,pvit2,pvit_find;

    for ( pvit1=points_and_values1.begin(), pvit2 =points_and_values2.begin();pvit1!=points_and_values1.end(); pvit1++,pvit2++ )
    {

        if( (*pvit1).first!=(*pvit2).first){

            pvit_find = points_and_values2.find((*pvit1).first);

            if(pvit_find==points_and_values2.end()){
                throw Exception( CSMP_ERROR,
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

template<size_t dim>
void ModelComparator<dim>::ReadRegionNodalScalarVariableAndNodeCoordinates(const Region<dim>& region,
                                                          const Index propKey,
                                                          map<Point<dim>, ScalarVariable >& points_and_values)
{
    ScalarVariable scalarValue( PLAIN, 0. );
    const typename std::vector<Node<dim>* >::const_iterator nodesEnd( region.NodesEnd() );

    for( typename std::vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != nodesEnd; ++it )
    {
        (*it)->Read( propKey, scalarValue );

        points_and_values[(*it)->Coordinate()] = scalarValue ;

    } // nodes of region

}

/////////////////////////////////////////////////////////////////////////////////////////
/// Comparison is done for meshes that can be in general different from each other     //
/////////////////////////////////////////////////////////////////////////////////////////

template<size_t dim>
double64 ModelComparator<dim>::CompareModelsAtPoints( Model<dim>& model1, Model<dim>& model2,
                                              const char* property1, const char* property2, const char* region )
{
    double64 error( 0. );

    Index model1Key( model1.Database().StorageKey( property1 ) );
    Index model2Key( model2.Database().StorageKey( property2 ) );

    if( model1Key.type != model2Key.type || model1Key.place != model2Key.place)
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareModelsAtPoints",
                         "Variables seem to be of different nature" );

    if( model1Key.type == SCALAR )
        error = CompareRegionScalarVariableAtPoints( model1, model2, property1, property2, region );
    else
    {
        cout<<"\nERROR:Variable has not valid for comparison type:"
            <<"\nFirst type: " <<parseType( model1Key.type )
            <<"\nSecond type: "<<parseType( model2Key.type )<<"\n";
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareModelsAtPossiblyDifferentsMeshes",
                         "Comparison for Variable type not implemented yet." );
    }

    return error;
}

template<size_t dim>
double64 ModelComparator<dim>::CompareVSetsAtPoints( VSet<dim>& vset1, VSet<dim>& vset2,
                                                     const char* property1, const char* property2,
                                                     const char* variablesFile1, const char* variablesFile2,
                                                     bool isoparametric, const char* region )
{
    double64 error( 0. );

    Model<dim> model1( vset1, variablesFile1, isoparametric );
    Model<dim> model2( vset2, variablesFile2, isoparametric );

    error = CompareModelsAtPoints( model1, model2, property1, property2, region );

    return error;
}


template<size_t dim>
double64 ModelComparator<dim>::CompareVSetsAtPoints( const char* vset1File, const char* vset2File,
                                                     const char* property1, const char* property2,
                                                     const char* variablesFile1, const char* variablesFile2,
                                                     bool isoparametric, const char* region )
{
    double64 error( 0. );

    VSet<dim> vset1, vset2;
    double64 modelTime( 0. );
    vset1.InputFrom( vset1File, modelTime );
    vset2.InputFrom( vset2File, modelTime );

    error = CompareVSetsAtPoints( vset1, vset2, property1, property2, variablesFile1, variablesFile2, isoparametric, region );

    return error;
}





template<size_t dim>
double64 ModelComparator<dim>::CompareRegionScalarVariableAtPoints( Model<dim>& model1, Model<dim>& model2,
                                                                    const char* property1, const char* property2, const char* region )
{
    std::vector<ScalarVariable> values;
    std::map<size_t, std::vector<double64> > points;

    const Region<dim>& rref1( model1.Region( region ) );
    const Region<dim>& rref2( model2.Region( region ) );

    if( rref1.Nodes() != rref2.Nodes() )
        throw Exception( CSMP_ERROR,
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
        throw Exception( CSMP_ERROR,
                         "ModelComparator::CompareRegionScalarVariableAtPoints",
                         "Comparison for this variable placement not implemented." );
    }

    double64 error( 0. );

    // finds values of certain property at given points
    PropertyAtPointVisitor<dim> pAt(model2, points, property2);
    pAt.SetBruteForceOff();
    model2.Accept(pAt);
    pAt.CheckResults();

    // L2 norm comparison
    size_t numberOfPoints( points.size() );
    ScalarVariable scalarValue( PLAIN, 0. );
    for (size_t i=0; i< numberOfPoints; i++ )
    {
        pAt.PropertyValueAt(i, scalarValue);
        error += (scalarValue() - values[i]())*(scalarValue() - values[i]());
    }


    return std::sqrt( error );

}


template<size_t dim>
void ModelComparator<dim>::ReadRegionNodalScalarVariableAndNodeCoordinates(const Region<dim>& region,
                                                          const Index propKey,
                                                          std::vector<ScalarVariable> & values,
                                                          map<size_t, std::vector<double64> >& points)
{
    size_t index(0U);
    ScalarVariable scalarValue( PLAIN, 0. );
    const typename std::vector<Node<dim>* >::const_iterator nodesEnd( region.NodesEnd() );

    values.clear();
    values.resize(region.Nodes());
    for( typename std::vector<Node<dim>*>::const_iterator it = region.NodesBegin(); it != nodesEnd; ++it )
    {
        (*it)->Read( propKey, scalarValue );
        for(size_t i=0U;i<dim;i++)
            points[index].push_back((*it)->operator[](i));

        values[index] = scalarValue ;
        index++;

    } // nodes of region

}


// explicits
template class ModelComparator<3U>;
template class ModelComparator<2U>;

} // csmp
