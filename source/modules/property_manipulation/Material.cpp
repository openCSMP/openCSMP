//
//  Material.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 13/6/20.
//  Copyright © 2020 Stephan Matthai. All rights reserved.
//

#include "Material.h"
#include "Model.h"

using namespace std;

namespace csmp {

/// transforms rational (-1,0..5..n) property values into material IDs stored on the elements
template<size_t dim>
void material_IDs_FromPropertyValues( Model<dim>& model, const std::string& elmt_prop_name )
 {
    const csmp::Index key(model.Database().StorageKey(elmt_prop_name.c_str()));
    if ( key.type!= SCALAR || key.place != ELEMENT )
      csmp::Exception( ERROR, "material_IDs_FromPropertyValues", "material identifier property must be a scalar placed on the element");

    Region<dim>& model_domain(model.Region("Model"));
    int32 mtrl_min(INT_MAX), mtrl_max(INT_MIN);
    
    for ( typename vector<Element<dim>*>::iterator 
          it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
        const double64 value = (*it)->Read(key);
        // checking that the property value can indeed be converted into an integer in a meaningful range
        // using the modulus operator % to determine whether the number has a decimal fraction
        if ( fmod( value, 1. ) != 0 ) {
             cerr <<"\n\tproperty value: "<< value << endl;
             csmp::Exception( ERROR, "material_IDs_FromPropertyValues", "material ID contains decimal places and can therefore not be converted to integer.");
          }
        const int32 material_identifier = static_cast<int32>(value);
        mtrl_min = min( mtrl_min, material_identifier );
        mtrl_max = max( mtrl_max, material_identifier ); 
        (*it)->Material_ID( material_identifier );
      }
  
    cout <<"\nmaterial_IDs_FromPropertyValues: successfully initialised material ID values in the range "<< mtrl_min <<" to "<< mtrl_max << endl;
    
 } // end material_IDs_FromPropertyValues

template void material_IDs_FromPropertyValues( Model<1U>&, const string& );
template void material_IDs_FromPropertyValues( Model<2U>&, const string& );
template void material_IDs_FromPropertyValues( Model<3U>&, const string& );






/// transforms rational (-1,0..5..n) property values into material IDs stored on the elements
template<size_t dim>
void propertyValuesPFromMaterial_IDs( Model<dim>& model, const std::string& elmt_prop_name )
 {
    const csmp::Index key(model.Database().StorageKey(elmt_prop_name.c_str()));
    if ( key.type!= SCALAR || key.place != ELEMENT )
      csmp::Exception( ERROR, "propertyValuesPFromMaterial_IDs", elmt_prop_name,
                              "must be a scalar placed on the element");

    Region<dim>& model_domain(model.Region("Model"));
    
    for ( typename vector<Element<dim>*>::iterator 
          it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         const double64 material_identifier = static_cast<double64>( (*it)->Material_ID() );
         (*it)->Store( key, makeScalar( (*it)->Status(key), material_identifier ) );
      }
  
    cout <<"\npropertyValuesPFromMaterial_IDs: successfully initialised '"<< elmt_prop_name <<"' from material ID values."<< endl;
    printRangeOfVariable( model, elmt_prop_name.c_str() );
    
 } // end propertyValuesPFromMaterial_IDs
 
template void propertyValuesPFromMaterial_IDs( Model<1U>&, const string& );
template void propertyValuesPFromMaterial_IDs( Model<2U>&, const string& );
template void propertyValuesPFromMaterial_IDs( Model<3U>&, const string& );
 
 
 
 
 

/**
    Removes stair-steps in material boundaries where possible; stair-steps are identified by two or more element faces on the outside of a material domain.
    
    Elements with completely isolated material IDs wiill also be removed.
    
    @attention Method uses BOX_BOUNDARY flags, i.e. AtBoundary() to determine whether an element is placed on a model boundary or interior
  
*/
template<size_t dim>
void smoothMaterialInterfaces( Model<dim>& model )
 {
    Region<dim>& model_domain(model.Region("Model"));
    size_t       modified_interior_elmts(0U),
                 modified_boundary_elmts(0U);
    
    for ( typename vector<Element<dim>*>::iterator 
          it=model_domain.ElementsBegin(); it!=model_domain.ElementsEnd(); ++it ) {
         if ( !isTriangular( (*it)->FE_Type() ) )
           csmp::Exception( ERROR, "smoothMaterialInterfaces", "thus far, this method has only been implemented for triangles.");
         if ( !(*it)->IsSurfaceElement() )
           csmp::Exception( ERROR, "smoothMaterialInterfaces", "thus far, this method only works for surface type elements.");
           
         // removing completely isolated elements
         const int32 mtrl_ID = (*it)->Material_ID();
         bool  nbor_with_same_ID(false);
         for ( size_t i=0U; i<(*it)->Neighbors(); ++i )
           if ( (*it)->Neighbor(i) != nullptr && (*it)->Neighbor(i)->Material_ID() == mtrl_ID ) {
                nbor_with_same_ID = true;
                break;
             }
         if ( !nbor_with_same_ID ) 
           cerr <<"\nsmoothMaterialInterfaces: element "<< (*it)->Idx() <<" has spatially isolated material ID. What shall be done?\n";
           
         // application to elements inside of Model domain
         if ( (*it)->AtBoundary() == NOT ) {
              // checking neighbors of adjacent pair of element faces: if 2 are the same their material key is assigned to current element
              if ( (*it)->Neighbor(0)->Material_ID() == (*it)->Neighbor(1)->Material_ID() &&
                   (*it)->Neighbor(0)->Material_ID() != (*it)->Material_ID() ) {
                    (*it)->Material_ID( (*it)->Neighbor(0)->Material_ID() );
                    modified_interior_elmts++;
                }

              else if ( (*it)->Neighbor(1)->Material_ID() == (*it)->Neighbor(2)->Material_ID() &&
                   (*it)->Neighbor(1)->Material_ID() != (*it)->Material_ID() ) {
                    (*it)->Material_ID( (*it)->Neighbor(1)->Material_ID() );
                    modified_interior_elmts++;
                }

              else if ( (*it)->Neighbor(2)->Material_ID() == (*it)->Neighbor(0)->Material_ID() &&
                   (*it)->Neighbor(2)->Material_ID() != (*it)->Material_ID() ) {
                    (*it)->Material_ID( (*it)->Neighbor(2)->Material_ID() );
                    modified_interior_elmts++;
                }
           }
         // boundary cases - the element must be isolated and is subsumed by its neighbors if they have the same material
         else {
              if ( (*it)->Neighbor(0) == nullptr ) {
                  if ( (*it)->Neighbor(1)->Material_ID() == (*it)->Neighbor(2)->Material_ID() &&
                       (*it)->Neighbor(1)->Material_ID() != (*it)->Material_ID() ) {
                        (*it)->Material_ID( (*it)->Neighbor(1)->Material_ID() );
                        modified_boundary_elmts++;
                    }
                }
              else if ( (*it)->Neighbor(1) == nullptr ) {
                  if ( (*it)->Neighbor(2)->Material_ID() == (*it)->Neighbor(0)->Material_ID() &&
                       (*it)->Neighbor(2)->Material_ID() != (*it)->Material_ID() ) {
                        (*it)->Material_ID( (*it)->Neighbor(2)->Material_ID() );
                        modified_boundary_elmts++;
                    }
                }
              else if ( (*it)->Neighbor(2) == nullptr ) {
                  if ( (*it)->Neighbor(1)->Material_ID() == (*it)->Neighbor(0)->Material_ID() &&
                       (*it)->Neighbor(1)->Material_ID() != (*it)->Material_ID() ) {
                        (*it)->Material_ID( (*it)->Neighbor(0)->Material_ID() );
                        modified_boundary_elmts++;
                    }
                }
           }
      }
  
   cout <<"\n\nsmoothMaterialInterfaces: modified "<< modified_interior_elmts <<" interior and ";
   cout << modified_boundary_elmts <<" boundary element material IDs.\n";

 } // end smoothMaterialInterfaces

template void smoothMaterialInterfaces( Model<1U>& );
template void smoothMaterialInterfaces( Model<2U>& );
template void smoothMaterialInterfaces( Model<3U>& );


} // end csmp
