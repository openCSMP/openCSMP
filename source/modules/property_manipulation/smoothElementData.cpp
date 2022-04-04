//
//  smoothElementData.cpp
//
//  Created by Stephan Matthai on 4/2/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "smoothElementData.h"
#include "Region.h"
#include "Model.h"
#include "Element.h"
#include "Node.h"
#include "ErrorHandler.h"
#include "CSMP_mathUtilities.h"
#include "VTU_Interface.h"

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

    vector<double>  face_nrml, vertical({0.,1.});
    eptr->UnitNormalToFace( face, face_nrml );
   
    // limit on the dot-product value from the permitted deviation angle of the normal
    const double n_degrees(60.);
    const double max_length( cos( degreesToRadians( n_degrees ) ) ); // positive
  
    // finding the dip of the face normal vector ignoring its azimuth
    // (if the dot-product is zero, the face is exactly vertical)
    const double dot_product = face_nrml[0] * vertical[0] + face_nrml[1] * vertical[1];
    if ( fabs(dot_product) >= max_length ) return false;

    return true;
   
 } // end isSubvertical (2D)
 
 
 inline bool isSubvertical( const Element<3U>* const eptr, size_t face )
 {
    assert( eptr != nullptr );
    assert( face < eptr->Faces() );

    vector<double>  face_nrml, vertical({0.,1.,0.});
    eptr->UnitNormalToFace( face, face_nrml );
   
    // limit on the dot-product value from the permitted deviation angle of the normal
    const double n_degrees(60.);
    const double max_length( cos( degreesToRadians( n_degrees ) ) ); // positive
  
    // finding the dip of the face normal vector ignoring its azimuth
    // (if the dot-product is zero, the face is exactly vertical)
    const double dot_product = face_nrml[0] * vertical[0] + face_nrml[1] * vertical[1] + face_nrml[2] * vertical[2];
    if ( fabs(dot_product) >= max_length ) return false;

    return true;
   
 } // end isSubvertical



/**
       Smoothes the target (element property) using a variaty of approaches and writes the model to disk afterwards.
         
       @attention variables file is assumed to be part of the binary file set.
       
       @todo perhaps use sqrt() of values for a milder form of value compression
*/
template<uint32_t dim>
void smoothElementData( Model<dim>& model, 
                        const std::string& region_to_be_smoothed, 
                        const std::string& variable_name,
                        int number_of_smoothing_cycles, 
                        bool log10_smoothing, bool in_plane_smoothing,
                        bool output_model_to_binary )
 {
     ErrorHandler&  csmp_error( ErrorHandler::Instance() );
 
     // 1. reading CSMP binary_native model from file
     // ---------------------------------------------
     Region<dim>&  domain(model.Region(region_to_be_smoothed));
     csmp::Index   key = model.Database().StorageKey(variable_name.c_str());
     double      min_val_database, max_val_database;
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
        for ( auto it=domain.CellsBegin(); it!=domain.PerimeterCellsBegin(); ++it ) {
                 if ( key.type == SCALAR ) {
                      (*it)->Store( log_key, makeScalar( (*it)->Status(key), log10( (*it)->Read(key)) ) );
                   }
                 else if ( key.type == VECTOR ) { 
                       (*it)->Read( key, vc );
                       for ( auto i{0U}; i<dim; ++i )
                         vc(i) = log10( vc[i] );
                       (*it)->Store( log_key, vc );
                   }
                 else if ( key.type == TENSOR ) {
                      (*it)->Read( key, ts );
                      for ( auto i{0U}; i<dim; ++i )
                        for ( auto j{0U}; j<dim; ++j )
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
                model.ExtrapolateCellToNodeProperty( variable_name.c_str(), node_prop_name.c_str() );
                model.InterpolateNodeToCellProperty( node_prop_name.c_str(), variable_name.c_str() );
             }                                                   
           if ( node_prop_created )
             model.DeleteProperty( node_prop_name.c_str() );
        }
   
      // 3. in-plane smoothing  of the data in multiple iterations
      // ---------------------------------------------------------
      else {
          vector<double> smoothed_vals(domain.InteriorCells());
          for ( int cycle=1U; cycle <= number_of_smoothing_cycles; cycle++ )
            {
              size_t elmt(0U);
              for ( auto it=domain.CellsBegin(); it!=domain.PerimeterCellsBegin(); ++it,  ++elmt )
                  {
                    // the new value taken is the volume-weighted mean average of the element neighbors and its own value
                    // current element
                    double average(0.), sum_of_weights(0.), min_val(0.), max_val(0.);
                    double elmt_volume = (*it)->Volume();
                    double val((*it)->Read(key));
                    min_val = min( min_val, val );
                    max_val = max( max_val, val );
                    average        += val * elmt_volume;
                    sum_of_weights += elmt_volume;
                    // sampling the element neighbors as long as they are inside of the target region
                    // and they sit across a subvertical face
                    for ( auto nbor=0U; nbor<(*it)->Neighbors(); ++nbor )
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
                for ( size_t i{0U}; i<domain.InteriorCells(); ++i ) {
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
       for ( auto it=domain.CellsBegin(); it!=domain.PerimeterCellsBegin(); ++it ) {
                if ( key.type == SCALAR ) {
                     (*it)->Store( original_key, makeScalar( (*it)->Status(key), pow( 10., (*it)->Read(key)) ) );
                  }
                else if ( key.type == VECTOR ) { 
                      (*it)->Read( key, vc );
                      for ( auto i{0U}; i<dim; ++i )
                        vc(i) = pow( 10., vc[i] );
                      (*it)->Store( original_key, vc );
                  }
                else if ( key.type == TENSOR ) {
                     (*it)->Read( key, ts );
                     for ( auto i{0U}; i<dim; ++i )
                       for ( auto j{0U}; j<dim; ++j )
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
    if ( output_model_to_binary ) {
         string  output_fileset( model.Name() );
         if ( log10_smoothing ) output_fileset +="log10-";
         output_fileset +="-smoothed";
         output_fileset += to_string(number_of_smoothing_cycles);
         // CSMP-native file output
         model.OutputToBinaryFile( output_fileset.c_str() );
      }

 } // end smoothElementData

template void smoothElementData( Model<2U>&, const string&, const string&, int, bool, bool, bool );
template void smoothElementData( Model<3U>&, const string&, const string&, int, bool, bool, bool );



/**
     Smoothing of porosity and permeability values by extrapolation to the node and interpolation back to the elements.
           - a single smoothing pass is applied.
           
     (other methods were tested with less success, see function body)

    developed for SMOOTHING OF CRC3-CRC2 2D CROSS-SECTION (SKM5/2/2020)
*/
template<uint32_t dim>
void smoothPorosityAndPermeabilityDistribution( Model<dim>& model )
  {
     const int  smoothing_passes(1);
     const bool log10_smoothing(false);  // creates more patchy pattern
     const bool in_plane_smoothing(false);

     smoothElementData( model, "Model", "permeability",  smoothing_passes, log10_smoothing, in_plane_smoothing );
     smoothElementData( model, "Model", "porosity",  smoothing_passes, log10_smoothing, in_plane_smoothing );
     
     printRangeOfVariable( model, "permeability" );
     printRangeOfVariable( model, "porosity" );
 
 } // end smooth data

template void smoothPorosityAndPermeabilityDistribution( Model<2U>& );
template void smoothPorosityAndPermeabilityDistribution( Model<3U>& );





/// Edoardo Pezulli's neighbor extrapolation based method, that avoids Element and Node objects located at the region boundary
template<uint32_t dim>
void spreadPropertiesOfInitialisedCellsAcross( typename vector<Element<dim>*>::iterator first,
                                               typename vector<Element<dim>*>::iterator last )
 {
    assert( first != last );
    sort( first, last );
    
    set<Node<dim>*>  new_nodes;
     
    const typename vector<csmp::Element<dim>*>::iterator elmts_end(last);
    for ( typename vector<Element<dim>*>::iterator it=first; it!=elmts_end; ++it )
      for ( auto i{0U}; i<(*it)->Nodes(); ++ i ) {
           assert( (*it)->N(i) != nullptr );
           new_nodes.insert( (*it)->N(i) );
        }
    
    //for (size_t e = 0; e < new_elmts.size(); ++e )
    typename vector<csmp::Element<dim>*>::iterator it(first);
    while( it != elmts_end ) 
      {
         //Set elements to search for starting with original un-initialized element
         // std::set<Element<dim>*> elmts_to_search{new_elmts[e]};
         set<Element<dim>*> elmts_to_search;
         elmts_to_search.insert( (*it) );
         bool found_initialized_el = false;   //condition when closest neighour uninitialized is found
         //Element traversal for uninitialized node
         while (found_initialized_el == false){
           //Set of neighbors to search next (if not found)
           std::set<Element<dim>*> elmts_searching_next;
           //Going over current element set
           for (typename set<Element<dim>*>::iterator eit = elmts_to_search.begin();
                eit != elmts_to_search.end(); ++eit){
             //Searching all neighbors for unitialized element
             const auto n_nbors{(*eit)->Neighbors()};
             for ( uint32_t nbor = 0; nbor < n_nbors; ++nbor) {
               if ( (*eit)->Neighbor(nbor) != nullptr){                 //Only consider existing neighbors
                 //if neighbor element is not new it must be initilized
                 // TODO: search elements in sorted vector using binary search (vector specific method rather than algorithm)
                 if ( find( first, elmts_end, (*eit)->Neighbor(nbor) ) == elmts_end ){
                   //Found an initialized neighbor
                   Element<dim>* found_el = (*eit)->Neighbor(nbor);
                   found_initialized_el = true;
                   //Performing copy of data of found element to original element
                   (*it)->LVS( found_el->LVS() );    //THIS WAS THE AIM: TODO: discuss intention here

                   //getting nodes of closest uni-initialized element (assumed to be most relevant to original el)
                   std::vector<Node<dim>*> closest_nodes = (*eit)->NodeVector();
                   //Now initializing all new nodes in the new element
                   for ( uint32_t n_new = 0 ; n_new < (*it)->Nodes(); ++n_new){
                     if ( find(new_nodes.begin(), new_nodes.end(), (*it)->N(n_new)) != new_nodes.end() ){
                       //Looking for node with same boundary flag if possible
                       BOX_BOUNDARY new_node_bound = (*it)->N(n_new)->AtBoundary();
                       //Getting all potential nodes and their box boundaries
                       std::multimap<BOX_BOUNDARY,Node<dim>*> potential_nodes;
                       for (auto n_found = 0; n_found < found_el->Nodes(); ++n_found){
                           Node<dim>* potential_node = found_el->N(n_found);             //getting potential node from initialized element
                           //if node is also within an unitialized elm (then its closest)
                           if ( std::find(closest_nodes.begin(), closest_nodes.end(), potential_node) != closest_nodes.end() ){
                             //then insert
                             potential_nodes.insert( make_pair(potential_node->AtBoundary(), potential_node));
                           }
                         }
                       //Now match node which matches boundary of new node
                       if (potential_nodes.find(new_node_bound) != potential_nodes.end() ){
                           assert( (*it)->N(n_new) != potential_nodes.find(new_node_bound)->second );         //check nodes arnt the same
                           (*it)->N(n_new)->LVS( potential_nodes.find(new_node_bound)->second->LVS() );      //THIS WAS THE AIM (and/or similar copying below)
                         }
                       //If nodes dont match, non boundary nodes have priority
                       else if ( potential_nodes.find( NOT ) != potential_nodes.end() ){
                         //else if we have an interior node - that gets used instead
                         (*it)->N(n_new)->LVS( potential_nodes.find(NOT)->second->LVS() );
                       } else if ( potential_nodes.find( INTERNAL ) != potential_nodes.end() ){
                         //Internal boundaries are also internal
                         (*it)->N(n_new)->LVS( potential_nodes.find(INTERNAL)->second->LVS() );
                       } else
                         throw (csmp::Exception( ERROR, "spreadPropertiesOfInitialisedCellsAcross",
                                                "Only Box Boundary nodes exist on element!"));
                       } //end of if node is new
                     } //end of node initialization

                     break; //end of neighbor search
                 } else {
                   //then neighbor is also uninitialized, we add to neighbor search
                   elmts_searching_next.insert((*eit)->Neighbor(nbor) );
                 }
               }//end of if nullptr
             } // end of neigbor search
             if (found_initialized_el == true )
               break; //stop elm serach if we found initialized element
           }//end of current element traversal

           //set neighbors as next to search
           elmts_to_search = elmts_searching_next;

         }//end of while loop

        it++;

       }//end of new elm initialization
       
       
//   }//end of if initialize

 } // end spreadPropertiesOfInitialisedCellsAcross

template void spreadPropertiesOfInitialisedCellsAcross<1U>( std::vector<Element<1U>*>::iterator, std::vector<Element<1U>*>::iterator );
template void spreadPropertiesOfInitialisedCellsAcross<2U>( std::vector<Element<2U>*>::iterator, std::vector<Element<2U>*>::iterator );
template void spreadPropertiesOfInitialisedCellsAcross<3U>( std::vector<Element<3U>*>::iterator, std::vector<Element<3U>*>::iterator );



} // end csmp
