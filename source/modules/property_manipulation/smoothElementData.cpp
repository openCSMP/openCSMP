//
//  smoothElementData.cpp
//
//  Created by Stephan Matthai on 4/2/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "smoothElementData.h"
#include "Region.h"
#include "Model.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"

using namespace std;

namespace csmp {


/**
    Checks the dip angle of the face normal; if the normal is close to horizontal
    the function returns true because the face is subvertical.
*/
inline bool isSubvertical( const Element<2U>* const eptr, size_t face )
 {
    assert( eptr != nullptr );
    assert( face < eptr->Faces() );

    vector<double64>  face_nrml, vertical({0.,1.});
    eptr->UnitNormalToFace( face, face_nrml );
   
    // limit on the dot-product value from the permitted deviation angle of the normal
    const double64 n_degrees(60.);
    const double64 max_length( cos( degreesToRadians( n_degrees ) ) ); // positive
  
    // finding the dip of the face normal vector ignoring its azimuth
    // (if the dot-product is zero, the face is exactly vertical)
    const double64 dot_product = face_nrml[0] * vertical[0] + face_nrml[1] * vertical[1];
    if ( fabs(dot_product) >= max_length ) return false;

    return true;
   
 } // end isSubvertical (2D)
 
 
 inline bool isSubvertical( const Element<3U>* const eptr, size_t face )
 {
    assert( eptr != nullptr );
    assert( face < eptr->Faces() );

    vector<double64>  face_nrml, vertical({0.,1.,0.});
    eptr->UnitNormalToFace( face, face_nrml );
   
    // limit on the dot-product value from the permitted deviation angle of the normal
    const double64 n_degrees(60.);
    const double64 max_length( cos( degreesToRadians( n_degrees ) ) ); // positive
  
    // finding the dip of the face normal vector ignoring its azimuth
    // (if the dot-product is zero, the face is exactly vertical)
    const double64 dot_product = face_nrml[0] * vertical[0] + face_nrml[1] * vertical[1] + face_nrml[2] * vertical[2];
    if ( fabs(dot_product) >= max_length ) return false;

    return true;
   
 } // end isSubvertical



/**
       Smoothes the target (element property) using a variaty of approaches and writes the model to disk afterwards.
         
       @attention variables file is assumed to be part of the binary file set.
       
       @todo perhaps use sqrt() of values for a milder form of value compression
*/
template<size_t dim>
void smoothElementData( Model<dim>& model, 
                        const std::string& region_to_be_smoothed, 
                        const std::string& variable_name,
                        int number_of_smoothing_cycles, 
                        bool log10_smoothing, bool in_plane_smoothing )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
     // 1. reading CSMP binary_native model from file
     // ---------------------------------------------
     Region<dim>&  domain(model.Region(region_to_be_smoothed));
     csmp::Index   key = model.Database().StorageKey(variable_name.c_str());
     double64      min_val_database, max_val_database;
     model.Database().RangeOf( variable_name.c_str(), min_val_database, max_val_database );
     bool interpolation_problem(false);

     if ( key.type != SCALAR ) {
          csmp_error.notice( ERROR, "smoothElementData", "method works only for scalar variables; not done");
          return;
       }
      if ( key.place != ELEMENT ) {
          csmp_error.notice( ERROR, "smoothElementData", "for the smoothing, the property variables must be placed on the element");
          return;
       }
   
      // 2. logarithmitizing the values, if possible (no zeroes)
      // -------------------------------------------------------
      bool can_be_logarithmically_smoothed(true), created_log_prop(false);
      if ( log10_smoothing ) {
           // can this be done
           if ( min_val_database <= 0. ) {
                cerr <<"\n\n value range of '"<< variable_name <<"'': ";
                cerr << min_val_database <<" - "<< max_val_database <<"\n";
                csmp_error.notice( ERROR, "smoothElementData:", "log10 of 0 or negative number is undefined; smoothing data as is...");
                can_be_logarithmically_smoothed = false;
             }
           // logarithmic values 
           VectorVariable<dim>  vc;
           TensorVariable<dim>  ts;
           string log_prop_name( string("log10 of ") + variable_name );
           created_log_prop = !model.Database().IsDefined( log_prop_name.c_str() );
           min_val_database = log10(min_val_database);
           max_val_database = log10(max_val_database);
           const csmp::Index log_key = model.CreateProperty( log_prop_name.c_str(), "SI", key.type, key.place, key.dataDepth,
                                                             log10(min_val_database), log10(max_val_database) );
           // taking the decadic logarithm of values
        for ( typename vector<Element<dim>*>::iterator it=domain.ElementsBegin(); it!=domain.PerimeterElementsBegin(); ++it ) {
                 if ( key.type == SCALAR ) {
                      (*it)->Store( log_key, makeScalar( (*it)->Status(key), log10( (*it)->Read(key)) ) );
                   }
                 else if ( key.type == VECTOR ) { 
                       (*it)->Read( key, vc );
                       for ( size_t i=0U; i<dim; ++i )
                         vc(i) = log10( vc[i] );
                       (*it)->Store( log_key, vc );
                   }
                 else if ( key.type == TENSOR ) {
                      (*it)->Read( key, ts );
                      for ( size_t i=0U; i<dim; ++i )
                        for ( size_t j=0U; j<dim; ++j )
                          ts(i,j) = log10( ts(i,j) );
                      (*it)->Store( log_key, ts );
                   }
                 else 
                 csmp_error.notice( ERROR, "smoothElementData", variable_name, "property type cannot be logarithmitised.");
             }
           // using the log10 of the value as opposed to original values                                                 
           key = log_key;
        }
   
   
      // 3. standard smoothing by extrapolation of properties to nodes and back
      // ----------------------------------------------------------------------
      if ( in_plane_smoothing == false )
        {
           const string node_prop_name( string("nodal ") + variable_name );
           const bool   node_prop_created( !model.Database().IsDefined( node_prop_name.c_str() ) );
           const csmp::Index log_key = model.CreateProperty( node_prop_name.c_str(), "SI", key.type, NODE, key.dataDepth,
                                                             min_val_database, max_val_database );
           // smoothing by extrapolation and interpolation
           for ( int cycle=1U; cycle <= number_of_smoothing_cycles; cycle++ ) {
                model.ExtrapolateElementToNodeProperty( variable_name.c_str(), node_prop_name.c_str() );
                model.InterpolateNodeToElementProperty( node_prop_name.c_str(), variable_name.c_str() );
             }                                                   
           if ( node_prop_created )
             model.DeleteProperty( node_prop_name.c_str() );
        }
   
      // 3. in-plane smoothing  of the data in multiple iterations
      // ---------------------------------------------------------
      else {
          vector<double64> smoothed_vals(domain.InteriorElements());
          for ( int cycle=1U; cycle <= number_of_smoothing_cycles; cycle++ )
            {
              size_t elmt(0U);
            for ( typename vector<Element<dim>*>::iterator it=domain.ElementsBegin(); it!=domain.PerimeterElementsBegin(); ++it,  ++elmt )
                {
                  // the new value taken is the volume-weighted mean average of the element neighbors and its own value
                  // current element
                  double64 average(0.), sum_of_weights(0.), min_val(0.), max_val(0.);
                  double64 elmt_volume = (*it)->Volume();
                  double64 val((*it)->Read(key));
                  min_val = min( min_val, val );
                  max_val = max( max_val, val );
                  average        += val * elmt_volume;
                  sum_of_weights += elmt_volume;
                  // sampling the element neighbors as long as they are inside of the target region
                  // and they sit across a subvertical face
                  for ( size_t nbor=0U; nbor<(*it)->Neighbors(); ++nbor )
                    if ( isSubvertical( (*it), nbor ) && domain.Contains((*it)->Neighbor(nbor)) )
                      {
                         elmt_volume = (*it)->Neighbor(nbor)->Volume();
                         val = (*it)->Read(key);
                         min_val = min( min_val, val );
                         max_val = max( max_val, val );
                         average      += (*it)->Neighbor(nbor)->Read(key) * elmt_volume;
                         sum_of_weights += elmt_volume;
                      }
                  // averaging
                  average /= sum_of_weights;
                  // storing interim values
                  smoothed_vals[elmt] = average;
                }
              
              // storing the smoothed values at the end of smoothing cycle
              for ( size_t i=0U; i<domain.InteriorElements(); ++i ) {
                  if ( smoothed_vals[i] >= min_val_database && smoothed_vals[i] <= max_val_database )
                     domain.E(i)->Store( key, makeScalar( domain.E(i)->Status(key), smoothed_vals[i] ) );
                  else {
                       cerr <<"\n\tElement: "<< domain.E(i)->Idx() <<": smoothing cycle "<< cycle;
                       cerr <<": computed out-of-range average value for '"<< variable_name <<"': "<< smoothed_vals[i];
                       interpolation_problem = true;
                    }
                }
              
            } // end smoothing cycles
       
         if ( interpolation_problem )
          throw csmp::Exception( ERROR, "smoothElementData",
                                 variable_name, "values computed during smoothing are out of range" );
                               
      } // end in_plane_smoothing

     // 4. returning the smoothed values into the non-logarithmic form and into the original variable
     // ---------------------------------------------------------------------------------------------
     if ( log10_smoothing and can_be_logarithmically_smoothed ) {
          const csmp::Index    original_key = model.Database().StorageKey(variable_name.c_str());
          VectorVariable<dim>  vc;
          TensorVariable<dim>  ts;

          // taking the decadic logarithm of values
       for ( typename vector<Element<dim>*>::iterator it=domain.ElementsBegin(); it!=domain.PerimeterElementsBegin(); ++it ) {
                if ( key.type == SCALAR ) {
                     (*it)->Store( original_key, makeScalar( (*it)->Status(key), pow( 10., (*it)->Read(key)) ) );
                  }
                else if ( key.type == VECTOR ) { 
                      (*it)->Read( key, vc );
                      for ( size_t i=0U; i<dim; ++i )
                        vc(i) = pow( 10., vc[i] );
                      (*it)->Store( original_key, vc );
                  }
                else if ( key.type == TENSOR ) {
                     (*it)->Read( key, ts );
                     for ( size_t i=0U; i<dim; ++i )
                       for ( size_t j=0U; j<dim; ++j )
                         ts(i,j) = pow( 10., ts(i,j) );
                     (*it)->Store( original_key, ts );
                  }
            }

          // deleting the auxiliary variable before the model is written to file  
          const string log_prop_name( string("log10 of ") + variable_name );
          if ( created_log_prop ) model.DeleteProperty( log_prop_name.c_str() );
       }

    // 5. saving the model to disk  (using the extension -smoothed#, where # is the number of iterations
    // -------------------------------------------------------------------------------------------------
    string  output_fileset( model.Name() );
    if ( log10_smoothing ) output_fileset +="log10-";
    output_fileset +="-smoothed";
    output_fileset += to_string(number_of_smoothing_cycles);
    
    // file output
    model.OutputToBinaryFile( output_fileset.c_str() );

 } // end smoothElementData

template void smoothElementData( Model<2U>&, const string&, const string&, int, bool, bool );

} // end csmp
