//
//  MaximumDifference.cpp
//  CSMP_GitHub_UnitTests
//
//  Created by Stephan Matthai on 26/9/21.
//  Copyright © 2021 Stephan Matthai. All rights reserved.
//

#include "MaximumDifference.h"
#include "Model.h"
#include "Region.h"

namespace csmp {

// INLINE METHODS

template<size_t dim>
void MaximumDifference<dim>::Calculate()
 {
    Var1.AssignTo( var1 );
    Var2.AssignTo( var2 );
    
    if ( std::fabs(var1() - var2()) > max_difference )
      max_difference = std::fabs(var1() - var2());
     
 } // end Calculate
 
 
template class MaximumDifference<1U>;
template class MaximumDifference<2U>;
template class MaximumDifference<3U>;



template<size_t dim>
double64 maximumDifference( const Model<dim>& sg,
                            const char* new_property, const char* old_property, bool normalise )
 {
   const bool      absolute(!normalise);
   double64        residual(0.), temp;
   ScalarVariable  new_prop, old_prop;
   
   const csmp::Region<dim>&  sgref(sg.Region("Model"));
   csmp::Index  new_key = sg.Database().StorageKey(new_property);
   csmp::Index  old_key = sg.Database().StorageKey(old_property);
   
   if ( new_key.place != old_key.place ) {
       throw csmp::Exception( ERROR, "maximumResidual",
                             "Properties to be compared do not have the same placement, residual cannot be computed, returning 0.0...!");
     }

   if ( new_key.place == NODE ) {
       for ( typename std::vector<Node<dim>*>::const_iterator
             nit = sgref.NodesBegin(); nit != sgref.NodesEnd(); nit++ ) {
           (*nit)->Read( new_key, new_prop );
           (*nit)->Read( old_key, old_prop );
           if ( !absolute ) temp = std::fabs( ( new_prop() - old_prop() ) / new_prop() );
           else             temp = std::fabs( new_prop() - old_prop() );
           if ( temp > residual ) residual = temp;
         }
     }

   else if ( new_key.place == ELEMENT ) {
       for ( typename std::vector<Element<dim>*>::const_iterator
             eit = sgref.ElementsBegin(); eit != sgref.ElementsEnd(); eit++ ) {
           (*eit)->Read( new_key, new_prop );
           (*eit)->Read( old_key, old_prop );
           if ( !absolute ) temp = std::fabs( ( new_prop() - old_prop() ) / new_prop() );
           else             temp = std::fabs( new_prop() - old_prop() );
           if ( temp > residual ) residual = temp;
         }
     }

   else if ( new_key.place == ELEMENT_INTEGRATION_POINT ) {
       for ( typename std::vector<Element<dim>*>::const_iterator
             eit = sgref.ElementsBegin(); eit != sgref.ElementsEnd(); eit++ )
         for ( size_t i=0U; i<(*eit)->IntegrationPoints(); i++ )
           {
             (*eit)->Read( i, new_key, new_prop );
             (*eit)->Read( i, old_key, old_prop );
             if ( !absolute ) temp = std::fabs( ( new_prop() - old_prop() ) / new_prop() );
             else             temp = std::fabs( new_prop() - old_prop() );
             if ( temp > residual ) residual = temp;
           }
     }
   
   else {
       throw csmp::Exception( ERROR, "maximumResidual",
                             "Cannot identify placement of variable, residual cannot be computed, returning 0.0...!");
     }

   return residual;
}

template double64 maximumDifference( const Model<1U>&, const char*, const char*, bool );
template double64 maximumDifference( const Model<2U>&, const char*, const char*, bool );
template double64 maximumDifference( const Model<3U>&, const char*, const char*, bool );




} // end csmp
