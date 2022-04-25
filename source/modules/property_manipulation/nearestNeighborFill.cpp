//
//  nearestNeighborInterpolation.cpp
//  Open CSMP++
//
//  Created by Stephan Matthai on 15/4/2022.
//

#include "nearestNeighborFill.h"

#include "Model.h"
#include "Region.h"
#include "ErrorHandler.h"

using namespace std;

namespace csmp {

/**
    replaces no-data values of target variable with nearest-neighbor values until there are none left, by default NAN's are no-data values
*/
template<uint32_t dim>
void nearestNeighborFill( Model<dim>& model, const char* target_region, const char* variable, double no_data_value )
 {
    ErrorHandler&  csmp_error( ErrorHandler::Instance() );
    if ( !model.Database().IsDefined(variable) ) {
        csmp_error.Note( WARNING, "nearestNeighborFill:", variable, "is not defined, nothing could be done.");
        return;
      }
    csmp::Index var_key = model.Database().StorageKey(variable);
    if ( var_key.type != SCALAR ) {
        csmp_error.Note( WARNING, "nearestNeighborFill:", "method currently only handles scalars, nothing could be done.");
        return;
      }
    if ( var_key.place != ELEMENT && var_key.place != NODE ) {
        csmp_error.Note( WARNING, "nearestNeighborFill:", "method currently only handles variables placed on element or nodes, nothing could be done.");
        return;
      }
    if ( !model.ContainsRegion(target_region) ) {
        csmp_error.Note( WARNING, "nearestNeighborFill:", target_region, "is not defined, nothing could be done.");
        return;
      }
    Region<dim>& gref(model.Region(target_region));
  
    // -------------------------------------------------------------------------------------------------------------
    // 1. processing element variables
    // -------------------------------------------------------------------------------------------------------------
    if ( var_key.place == ELEMENT ) {
         // counting no-data values, and memorizing pointers to elements with such values
         set<Element<dim>*> elementsMissingDataValues, filledValues;
         for ( typename vector<Element<dim>*>::const_iterator it=gref.CellsBegin(); it!=gref.CellsEnd(); it++ )
           if ( isnan((*it)->Read(var_key)) || fabs(no_data_value-(*it)->Read(var_key)) < numeric_limits<double>::epsilon() )
             elementsMissingDataValues.insert( (*it) );
          
         if ( elementsMissingDataValues.empty() ) {
              csmp_error.Note( INFO, "nearestNeighborFill:", "all elements have valid data values, nothing was done.");
              return;
           }
         // starting nearest neighbor-fill loop
         else {
              cout <<"\nnearestNeighborFill: detected "<< elementsMissingDataValues.size() <<" elements with no-data values; filling these now.\n";
              do {
                 // looping over the element neighbors, collecting data values for later weighted averaging,
                 // using the inverse of the barycenter to barycenter distance as weighting factor
                 for ( typename set<Element<dim>*>::iterator
                       it=elementsMissingDataValues.begin(); it!=elementsMissingDataValues.end(); it++ )
                   {
                      set<pair<double,double> > valuesAndWeights;
                      for ( auto i{0U}; i<(*it)->Neighbors(); i++ )
                        if ( (*it)->Neighbor(i) != NULL )
                          {
                             Point<dim> bctr = (*it)->BaryCenter();
                             double   nval = (*it)->Neighbor(i)->Read(var_key);
                             // if the neighbor element exists and has a valid variable value, the barycentric distance is determined and stored
                             if ( !isnan(nval) && fabs(no_data_value-nval) > numeric_limits<double>::epsilon() ) {
                                  double   distance  = bctr.DistanceTo( (*it)->Neighbor(i)->BaryCenter() );
                                  valuesAndWeights.insert( make_pair(nval,1./distance) );
                               }
                          }
                      // assigned a weighted average to the element if possible
                      if ( !valuesAndWeights.empty() ) {
                           double  sumOfWeights(0.), sumOfWeightedVals(0.);
                           for ( set<pair<double,double> >::const_iterator
                                 vit=valuesAndWeights.begin(); vit!=valuesAndWeights.end(); vit++ ) {
                                sumOfWeightedVals += (*vit).first * (*vit).second;
                                sumOfWeights      += (*vit).second;
                             }
                           (*it)->Store( var_key, makeScalar(PLAIN,sumOfWeightedVals/sumOfWeights) );
                           filledValues.insert( (*it) );
                        }
                   }
                  // removing the successfully filled elements from elementsMissingDataValues
                  for ( typename set<Element<dim>*>::const_iterator it=filledValues.begin(); it!=filledValues.end(); it++ )
                    elementsMissingDataValues.erase( (*it) );
                  filledValues.clear();
                }
              while ( !elementsMissingDataValues.empty() );
          }
        return;
      } // end element variable processing


    // -------------------------------------------------------------------------------------------------------------
    // 2. processing node variables
    // -------------------------------------------------------------------------------------------------------------
    if ( var_key.place == NODE ) {
         // counting no-data values, and memorizing pointers to elements with such values
         set<Node<dim>*> nodesMissingDataValues, filledValues;
         for ( typename vector<Node<dim>*>::const_iterator it=gref.NodesBegin(); it!=gref.NodesEnd(); it++ )
           if ( isnan((*it)->Read(var_key)) || fabs(no_data_value-(*it)->Read(var_key)) < numeric_limits<double>::epsilon() )
             nodesMissingDataValues.insert( (*it) );
          
         if ( nodesMissingDataValues.empty() ) {
              csmp_error.Note( INFO, "nearestNeighborFill:", "all nodes have valid data values, nothing was done.");
              return;
           }
         // starting nearest neighbor-fill loop
         else {
              cout <<"\nnearestNeighborFill: detected "<< nodesMissingDataValues.size() <<" nodes with no-data values; filling these now.\n";
              do {
                 // looping over the element neighbors, collecting data values for later weighted averaging,
                 // using the inverse of the barycenter to barycenter distance as weighting factor
                 for ( typename set<Node<dim>*>::iterator
                       it=nodesMissingDataValues.begin(); it!=nodesMissingDataValues.end(); it++ )
                   {
                      set<pair<double,double> > valuesAndWeights;
                      for ( auto i{0U}; i<(*it)->Neighbors(); i++ )
                        if ( (*it)->Neighbor(i) != NULL )
                          {
                             Point<dim> nxyz = (*it)->Coordinate();
                             double   nval = (*it)->Neighbor(i)->Read(var_key);
                             // if the neighbor element exists and has a valid variable value, the barycentric distance is determined and stored
                             if ( !isnan(nval) && fabs(no_data_value-nval) > numeric_limits<double>::epsilon() ) {
                                  double   distance  = nxyz.DistanceTo( (*it)->Neighbor(i)->Coordinate() );
                                  valuesAndWeights.insert( make_pair(nval,1./distance) );
                               }
                          }
                      // assigned a weighted average to the element if possible
                      if ( !valuesAndWeights.empty() ) {
                           double  sumOfWeights(0.), sumOfWeightedVals(0.);
                           for ( set<pair<double,double> >::const_iterator
                                 vit=valuesAndWeights.begin(); vit!=valuesAndWeights.end(); vit++ ) {
                                sumOfWeightedVals += (*vit).first * (*vit).second;
                                sumOfWeights      += (*vit).second;
                             }
                           (*it)->Store( var_key, makeScalar(PLAIN,sumOfWeightedVals/sumOfWeights) );
                           filledValues.insert( (*it) );
                        }
                   }
                  // removing the successfully filled elements from elementsMissingDataValues
                  for ( typename set<Node<dim>*>::const_iterator it=filledValues.begin(); it!=filledValues.end(); it++ )
                    nodesMissingDataValues.erase( (*it) );
                  filledValues.clear();
                }
              while ( !nodesMissingDataValues.empty() );
          }
        return;
      } // end node variable processing

 } // end nearestNeighborFill
 
template void nearestNeighborFill( Model<1U>&, const char*, const char*, double );
template void nearestNeighborFill( Model<2U>&, const char*, const char*, double );
template void nearestNeighborFill( Model<3U>&, const char*, const char*, double );



} // end csmp
